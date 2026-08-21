"""
protoSpaceJAM against a cell-line variant set, on the nine designs whose guides were verified by
hand.

These are *every* OpenCell design whose gRNA recognition site (protospacer + PAM) carries a
non-reference base in czML1175_KOLF2.1J-Parental hap1 -- nine designs, nine bases, nine VCF
records, out of 1,201 designable rows (0.7%).  They are the cases where cell line mode changes
the oligo you would order, so they are the cases worth pinning.

    Entry  gene           chr:pos        region  off  ref>var  GT   GQ  DP  weight
      63   RABAC1-201     19:41959295    distal    4  C>T      1|0  36  32  0.75
     114   NEK9-201       14:75084551    distal    3  C>T      1/1  58  61  0.90
     146   RABGAP1L-201   1:174219147    distal    5  G>A      1/1  56  81  0.90
     191   ACTR2-201      2:65268753     distal    7  G>A      1/1  54  73  0.90
     273   BLVRA-201      7:43771165     seed     11  G>A      1/1  52  67  0.50
     293   MAP1LC3B-201   16:87392422    distal    7  T>C      1/1  53  48  0.90
     866   ARHGEF7-204    13:111217690   seed     18  T>C      1/1  55  39  0.50
     966   CDK1-203       10:60794005    seed     11  A>G      1/1  57  76  0.50
    1059   ULK3-201       15:74843106    distal    4  T>C      1|0  32  45  0.75

`off` is the 0-based offset along the protospacer in guide orientation (0 is PAM-distal, 19 abuts
the PAM); `ref>var` are the genomic bases, which appear complemented in the six minus-strand
guides.  Seven of nine are hom-alt with zero reference reads, the two hets are phased onto the
selected haplotype, all nine carry GQ >= 32, and none sits in a homopolymer.  Entry 1059
(ULK3-201) is the only design in the whole OpenCell run where cell line mode picked a *different*
guide than the reference run.

How the fixtures were made, and what makes them trustworthy
-----------------------------------------------------------
`tests/GroundTruths/quick_variant_*.csv` are this pipeline's own output -- but before they were
promoted to fixtures, all nine guides were checked against `variant/manual_eval_set/`, which is
independent evidence: `extract_grna_diffs.py` slices the guide window out of both genome pickles
(`fa_pickle/GRCh38/` and `fa_pickle/GRCh38_KOLF2.1J_hap1/`), reverse-complements on the minus
strand and diffs base by base against the source VCF, never reading protoSpaceJAM's own
`gRNA_seq` columns.  Protospacer, PAM, reference protospacer, guide coordinates, variant weight
and the reported variant string all matched for all nine.

**If you ever regenerate these fixtures, re-run that check first.**  Regenerating them blind turns
whatever the code now emits into the new "correct" answer, which is exactly the failure mode the
manual verification exists to prevent.

What this input does not cover
------------------------------
No variant landed in a PAM in this run, and none of the nine sits in a no-call region, so
`pam_het`, `pam_hom`, `seed_het` and `no_call_region` are not exercised here.  `indel_selected_hap`
is -- two of the nine have an unrepresentable indel in their homology arms.
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
    """
    The OpenCell parameters, verbatim from run_full_test_variant.py.

    Guide choice has to be reproduced exactly or the verified sequences do not apply, and ranking
    is sensitive to these.  num_gRNA_per_design stays at its default of 1, as the OpenCell run had
    it: a penalized guide only reaches rank 1 when there was no better alternative, which is what
    makes these nine the interesting ones.
    """
    return {
        "path2csv": os.path.join("input", "test_input_variant.csv"),
        "outdir": outdir,
        "ssODN_max_size": 200,
        "Npayload": "ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT",
        "Cpayload": "GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
        "Strand_choice": "NonTargetStrand",
        "recode_order": "PAM_first",
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
        Every region x zygosity combination *these nine reach* is still being classified.

        Without this, a change that silently stopped classifying (say) seed variants would still
        pass the fixture comparison -- it would just quietly regenerate a new "correct" answer the
        next time someone refreshed the ground truths.

        The nine carry four distal_hom, three seed_hom and two distal_het, and no PAM variant at
        all, so `pam_het`, `pam_hom` and `seed_het` cannot be asserted from this input.  They are
        reachable through the weight table in util/variant_annot.py and through the full OpenCell
        run; if you need them pinned in a quick test, they need loci chosen for it.
        """
        flags = set()
        for r in self._result_rows():
            flags.update(w for w in r["variant_warnings"].split(";") if w)
        for expected in ("distal_het", "distal_hom", "seed_hom"):
            self.assertIn(expected, flags, f"{expected} is not exercised by the test input")
        self.assertIn("gRNA_seq_differs_from_reference", flags,
                      "all nine guides are supposed to differ from the reference sequence")

    def test_every_design_has_a_variant_in_its_guide(self):
        """
        The defining property of this input: all nine, no exceptions.

        A design dropping out, or the ranking quietly preferring an unpenalized guide, would
        otherwise show up only as a fixture diff someone might refresh away.
        """
        rows = self._result_rows()
        self.assertEqual(len(rows), 9, "expected exactly nine designs")
        for r in rows:
            with self.subTest(design=r["ID"]):
                self.assertTrue(r["variants_in_protospacer_PAM"],
                                f"{r['ID']} was chosen for having a variant in its guide")
                self.assertLess(float(r["variant_weight"]), 1.0)
                self.assertEqual(r["spec_score_stale"], "True")
                self.assertTrue(r["gRNA_seq_ref"], "the reference guide sequence should be reported")
                self.assertNotEqual(r["gRNA_seq"], r["gRNA_seq_ref"])

    def test_unrepresented_variant_reasons_are_covered(self):
        """
        A window can be untrustworthy while every SNV in it is correct.

        Only `indel_selected_hap` is reachable from these nine; `no_call_region` is not, so it is
        not asserted here.  If both need pinning in a quick test, that needs a locus chosen for it.
        """
        import csv

        with open(os.path.join(self.runtime_out, "variants_report.csv")) as fh:
            reasons = {r["skip_reason"] for r in csv.DictReader(fh)}
        self.assertIn("indel_selected_hap", reasons,
                      "two of the nine have an unrepresentable indel in their homology arms")


if __name__ == "__main__":

    unittest.main()
