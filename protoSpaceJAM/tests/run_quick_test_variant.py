"""
protoSpaceJAM against a cell-line variant set, on six loci chosen to cover the weight table.

The input rows, and what each one is here to exercise (genotypes are czML1175_KOLF2.1J-Parental
hap1, read off the filtered joint callset):

    ENST00000463591 C   PRDM16-204   1:3386696 A>G het at the stop codon
                                     -> distal_het 0.75 and seed_het 0.25 on ranks 2/3
    ENST00000514428 N   RERE-218     1:8624343 C>T het inside a PAM
                                     -> weight 0, PAM_differs_from_reference
    1:946247                         a lone hom-alt substitution: exactly one base of the donor
                                     differs from the reference, distal_hom 0.9 on rank 3
    1:978812                         hom-alt in the seed of all three guides -> seed_hom 0.5
    1:1168310                        a 2 bp insertion at the edit position and a 35 bp deletion
                                     301 bp away: one indel_selected_hap survives into the final
                                     donor, the other is centered out
    1:1157315                        a 40 bp region where this sample has no genotype call at all
                                     -> no_call_region, the reason that exists because the
                                     single-sample VCF was extracted with --exclude-uncalled

Loci 3-6 were chosen against czML1175 *and* against the filtered callset.  Do not pick them from
an upstream VCF: the previous fixture used 1:977899, which is a GQ=5 call and vanishes the moment
a quality filter is applied.
"""

import os
import shutil
import unittest

from protoSpaceJAM.protoSpaceJAM import main as pJAM

VARIANT_SET = "GRCh38_KOLF2.1J_hap1"

# The variant assets are built by precompute/variant_sets/4_make_variant_genome.py and live under
# genome_files/, which is gitignored.  Skip rather than fail when they have not been built here.
_variant_set_dir = os.path.join(
    "protoSpaceJAM", "genome_files", "variant_sets", VARIANT_SET, "manifest.json"
)
_variant_genome_dir = os.path.join(
    "protoSpaceJAM", "genome_files", "fa_pickle", VARIANT_SET
)


def _base_args(outdir):
    return {
        "path2csv": os.path.join("input", "test_input_variant.csv"),
        "outdir": outdir,
        "num_gRNA_per_design": 3,
        "ssODN_max_size": 200,
        "test_mode": True,
    }


@unittest.skipUnless(
    os.path.isfile(_variant_set_dir),
    f"variant set {VARIANT_SET} not built (see precompute/variant_sets/4_make_variant_genome.py)",
)
class test_variant_aware_design(unittest.TestCase):

    runtime_out = os.path.join("tests", "quick_variant_test_result")
    oracle_out = os.path.join("tests", "quick_variant_test_result_materialized")

    @classmethod
    def setUpClass(cls):
        print("running protoSpaceJAM against a variant set (quick variant test)...")
        try:
            os.chdir(os.path.join("protoSpaceJAM"))
        except FileNotFoundError:
            raise FileNotFoundError("protoSpaceJAM directory not found, are you in the repo's root directory?")

        if not os.path.isfile("protoSpaceJAM.py"):
            raise FileNotFoundError("protoSpaceJAM.py not found in current directory, are you in the protoSpaceJAM directory?")

        for d in (cls.runtime_out, cls.oracle_out):
            if os.path.exists(d):
                shutil.rmtree(d)

        # the runtime path: substitutions applied to each slice as get_seq() fetches it
        args = _base_args(cls.runtime_out)
        args["variant_set"] = VARIANT_SET
        pJAM(args)

        # the same designs off the pre-substituted genome pickles.  This is the oracle the
        # runtime path is defined against; if the two ever disagree, one of them is wrong.
        cls.oracle_available = os.path.isdir(
            os.path.join("genome_files", "fa_pickle", VARIANT_SET)
        )
        if cls.oracle_available:
            args = _base_args(cls.oracle_out)
            args["variant_genome"] = VARIANT_SET
            pJAM(args)

    @classmethod
    def tearDownClass(cls):
        keep_results = True
        if not keep_results:
            for d in (cls.runtime_out, cls.oracle_out):
                if os.path.exists(d):
                    shutil.rmtree(d)

    def _compare_to_ground_truth(self, filename, ground_truth):
        expected_path = os.path.join("tests", "GroundTruths", ground_truth)
        new_path = os.path.join(self.runtime_out, filename)
        self.assertTrue(os.path.isfile(new_path), f"{new_path} was not written")
        with open(expected_path) as f1, open(new_path) as f2:
            content1 = f1.read().strip().replace("\r\n", "\n")
            content2 = f2.read().strip().replace("\r\n", "\n")
        self.assertEqual(content1, content2)

    def test_if_generated_results(self):
        self.assertTrue(os.path.isfile(os.path.join(self.runtime_out, "result.csv")))
        self.assertTrue(os.path.isfile(os.path.join(self.runtime_out, "variants_report.csv")))

    def test_compare_results(self):
        print("comparing results...")
        self._compare_to_ground_truth("result.csv", "quick_variant_result.csv")

    def test_compare_variants_report(self):
        print("comparing variants report...")
        self._compare_to_ground_truth("variants_report.csv", "quick_variant_variants_report.csv")

    def test_runtime_patch_matches_materialized_genome(self):
        """--variant_set and --variant_genome are two routes to the same sequence."""
        if not self.oracle_available:
            self.skipTest(f"materialized genome {VARIANT_SET} not built")
        for filename in ("result.csv", "variants_report.csv"):
            with open(os.path.join(self.runtime_out, filename)) as f1, \
                 open(os.path.join(self.oracle_out, filename)) as f2:
                self.assertEqual(f1.read(), f2.read(), f"{filename} differs between the two paths")

    def _result_rows(self):
        import csv

        with open(os.path.join(self.runtime_out, "result.csv")) as fh:
            return [r for r in csv.DictReader(fh) if r.get("variant_source")]

    def test_variant_columns_present_and_populated(self):
        """The columns exist, and at least one design actually reports a variant."""
        rows = self._result_rows()
        self.assertTrue(rows, "no rows carried the variant columns")
        for r in rows:
            self.assertEqual(r["variant_source"], VARIANT_SET)
            self.assertEqual(r["variant_haplotype"], "1")
        self.assertTrue(
            any(r["variants_in_donor"] for r in rows),
            "the test input is supposed to hit variants in at least one donor",
        )
        self.assertTrue(
            any(r["variants_in_protospacer_PAM"] for r in rows),
            "the test input is supposed to hit at least one protospacer/PAM",
        )
        self.assertTrue(
            any(float(r["variant_weight"]) < 1.0 for r in rows),
            "no guide was penalized, so the penalty path is untested",
        )

    def test_weight_table_is_covered(self):
        """
        Every region x zygosity combination the guide penalty can produce is reached.

        Without this, a change that silently stopped classifying (say) seed variants would still
        pass the fixture comparison -- it would just quietly regenerate a new "correct" answer the
        next time someone refreshed the ground truths.
        """
        flags = set()
        for r in self._result_rows():
            flags.update(w for w in r["variant_warnings"].split(";") if w)
        for expected in ("distal_het", "seed_het", "pam_het", "distal_hom", "seed_hom"):
            self.assertIn(expected, flags, f"{expected} is not exercised by the test input")

    def test_unrepresented_variant_reasons_are_covered(self):
        """The two ways a window can be untrustworthy while the SNVs are all correct."""
        import csv

        with open(os.path.join(self.runtime_out, "variants_report.csv")) as fh:
            reasons = {r["skip_reason"] for r in csv.DictReader(fh)}
        self.assertIn("indel_selected_hap", reasons)
        self.assertIn(
            "no_call_region",
            reasons,
            "no-call regions are not being reported; is the set built without --nocall_bed?",
        )


if __name__ == "__main__":

    unittest.main()
