"""
The OpenCell design (1,236 transcripts) run against a cell-line variant set, both ways.

The point of this test is the *comparison*, not a fixture: --variant_genome reads the
pre-substituted genome pickles, --variant_set patches each slice at fetch time from the sidecar,
and the two must produce byte-identical output.  6_validate_runtime_patch.py already asserts that at
the get_seq() level on random windows; this asserts it at the design level, on real homology arms,
real guide choices, real recoding and real recut CFD, across a whole production-sized job.

It also re-runs the same input with no variant flag and diffs against the reference ground truth,
because "the feature is inert unless asked for" is the constraint the whole design rests on.

    python protoSpaceJAM/tests/run_full_test_variant.py

Slow: three full OpenCell runs, launched concurrently when the `protoSpaceJAM` console script is
installed.  Skips itself when the variant set has not been built here.
"""

import csv
import os
import shutil
import subprocess
import unittest

from protoSpaceJAM.protoSpaceJAM import main as pJAM

VARIANT_SET = "GRCh38_KOLF2.1J_hap1"

_manifest = os.path.join(
    "protoSpaceJAM", "genome_files", "variant_sets", VARIANT_SET, "manifest.json"
)

# the OpenCell design as run_full_test.py runs it; only outdir and the variant flag change
_OPENCELL_ARGS = {
    "path2csv": os.path.join("input", "OpenCell_protospaceX_design.csv"),
    "ssODN_max_size": 200,
    "Npayload": "ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT",
    "Cpayload": "GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
    "Strand_choice": "NonTargetStrand",
    "recode_order": "PAM_first",
    "test_mode": True,
}


def _args(outdir, **extra):
    args = dict(_OPENCELL_ARGS)
    args["outdir"] = outdir
    args.update(extra)
    return args


def _cli_args(outdir, **extra):
    """The same job as a `protoSpaceJAM` command line."""
    argv = ["protoSpaceJAM"]
    for key, value in _args(outdir, **extra).items():
        argv += [f"--{key}", str(value)]
    return argv


def _run_jobs(jobs):
    """
    Run several full OpenCell designs, concurrently when the console script is installed.

    Three sequential OpenCell runs are about three CPU-hours; they share nothing but read-only
    precomputed assets, so they parallelize perfectly.  `pip install -e .` puts a `protoSpaceJAM`
    entry point on PATH -- when it is missing (a bare source checkout) fall back to running them
    one after another in-process, which is slow but has no extra requirement.

    `jobs` is a list of (outdir, extra_kwargs).
    """
    if shutil.which("protoSpaceJAM") is None:
        for outdir, extra in jobs:
            pJAM(_args(outdir, **extra))
        return

    procs = []
    for outdir, extra in jobs:
        log = open(outdir + ".log", "w")
        procs.append(
            (outdir, log, subprocess.Popen(_cli_args(outdir, **extra), stdout=log, stderr=log))
        )
    failures = []
    for outdir, log, proc in procs:
        rc = proc.wait()
        log.close()
        if rc != 0:
            failures.append(f"{outdir} exited {rc}; see {outdir}.log")
    if failures:
        raise RuntimeError("; ".join(failures))


@unittest.skipUnless(
    os.path.isfile(_manifest),
    f"variant set {VARIANT_SET} not built (see precompute/variant_sets/4_make_variant_genome.py)",
)
class test_OpenCell_variant_aware(unittest.TestCase):

    sidecar_out = os.path.join("tests", "OpenCell_design_variant_set")
    genome_out = os.path.join("tests", "OpenCell_design_variant_genome")
    reference_out = os.path.join("tests", "OpenCell_design_reference")

    @classmethod
    def setUpClass(cls):
        print("running the OpenCell design against " + VARIANT_SET + " (full variant test)...")
        try:
            os.chdir(os.path.join("protoSpaceJAM"))
        except FileNotFoundError:
            raise FileNotFoundError(
                "protoSpaceJAM directory not found, are you in the repo's root directory?"
            )
        if not os.path.isfile("protoSpaceJAM.py"):
            raise FileNotFoundError(
                "protoSpaceJAM.py not found in current directory, are you in the "
                "protoSpaceJAM directory?"
            )

        for d in (cls.sidecar_out, cls.genome_out, cls.reference_out):
            if os.path.exists(d):
                shutil.rmtree(d)

        _run_jobs([
            (cls.sidecar_out, {"variant_set": VARIANT_SET}),
            (cls.genome_out, {"variant_genome": VARIANT_SET}),
            (cls.reference_out, {}),
        ])

    def _rows(self, outdir):
        with open(os.path.join(outdir, "result.csv")) as fh:
            return list(csv.DictReader(fh))

    def test_if_generated_results(self):
        for d in (self.sidecar_out, self.genome_out, self.reference_out):
            self.assertTrue(os.path.isfile(os.path.join(d, "result.csv")), f"{d}/result.csv")
        for d in (self.sidecar_out, self.genome_out):
            self.assertTrue(
                os.path.isfile(os.path.join(d, "variants_report.csv")),
                f"{d}/variants_report.csv",
            )

    def test_sidecar_matches_substituted_genome(self):
        """The whole point: two routes to the same sequence, one job, zero differences."""
        for filename in ("result.csv", "variants_report.csv"):
            with open(os.path.join(self.sidecar_out, filename)) as f1, \
                 open(os.path.join(self.genome_out, filename)) as f2:
                self.assertEqual(
                    f1.read(), f2.read(), f"{filename} differs between the two variant paths"
                )

    def test_reference_run_is_unchanged(self):
        """No variant flag must still reproduce the committed reference ground truth."""
        expected = os.path.join(
            "tests", "GroundTruths", "OpenCell_design_full_recoding_result.csv"
        )
        with open(expected) as f1, open(os.path.join(self.reference_out, "result.csv")) as f2:
            self.assertEqual(
                f1.read().strip().replace("\r\n", "\n"),
                f2.read().strip().replace("\r\n", "\n"),
            )

    def test_variant_columns_populated(self):
        """A full-callset variant set must actually fire, on a job this size."""
        rows = self._rows(self.sidecar_out)
        self.assertTrue(rows)
        with_donor_variants = [r for r in rows if r.get("variants_in_donor")]
        with_guide_variants = [r for r in rows if r.get("variants_in_protospacer_PAM")]
        penalized = [r for r in rows if float(r.get("variant_weight") or 1) < 1.0]
        print(
            f"  {len(rows)} designs: {len(with_donor_variants)} with a variant in the donor, "
            f"{len(with_guide_variants)} in the protospacer/PAM, {len(penalized)} penalized"
        )
        self.assertTrue(with_donor_variants)
        self.assertTrue(with_guide_variants)
        self.assertTrue(penalized)

    def test_variant_run_differs_from_reference(self):
        """Personalizing the sequence has to change something, or nothing is being applied."""
        ref = self._rows(self.reference_out)
        var = self._rows(self.sidecar_out)
        self.assertEqual(len(ref), len(var), "the two runs designed different numbers of rows")
        differing = sum(1 for a, b in zip(ref, var) if a["DNA donor"] != b["DNA donor"])
        print(f"  {differing}/{len(ref)} donors differ from the reference run")
        self.assertTrue(differing, "no donor changed; the variant set is not being applied")


if __name__ == "__main__":
    unittest.main()
