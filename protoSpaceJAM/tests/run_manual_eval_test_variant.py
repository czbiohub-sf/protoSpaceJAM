"""
Cell-line (variant) mode checked against the manual eval set.

The nine designs here are *every* OpenCell design whose gRNA recognition site carries a
non-reference base in czML1175_KOLF2.1J-Parental hap1 -- nine designs, nine bases, nine VCF
records, out of 1,201 designable rows (0.7%).  They are the cases where cell line mode changes
the oligo you would order, so they are the cases worth pinning.

Where the expected values come from
-----------------------------------
`variant/manual_eval_set/extract_grna_diffs.py`, which does **not** read protoSpaceJAM's
`gRNA_seq` / `gRNA_seq_ref` columns.  It slices the guide window out of both genome pickles
(`fa_pickle/GRCh38/` and `fa_pickle/GRCh38_KOLF2.1J_hap1/`), reverse-complements on the minus
strand, and diffs base by base against the source VCF.  So CASES below is independent evidence,
not a snapshot of what this pipeline happened to emit -- which is the only reason it is worth
asserting against.  `variant/` is gitignored and local; the table is transcribed here so the test
does not depend on it surviving.

What this covers that run_quick_test_variant.py does not
--------------------------------------------------------
Real designs, and the guide *sequence* rather than the weight bookkeeping: that the personalized
protospacer differs from the reference at exactly the base the VCF says, at the right offset, in
the right orientation.  Entry 1059 (ULK3-201) additionally pins the one design in the whole
OpenCell run where cell line mode picked a *different* guide.

What it does not cover, and why the six-locus test stays
--------------------------------------------------------
No variant landed in a PAM in this run, and none of the nine involves an indel or a no-call
region.  So `pam_het`, `seed_het`, `indel_selected_hap` and `no_call_region` are reachable only
through `run_quick_test_variant.py`.  The two tests are complements; deleting either loses
coverage.

Genotype notes for the set: seven of nine are hom-alt with zero reference reads (AD 0,N at
DP 39-81), two are phased hets on the selected haplotype with balanced support, all nine carry
GQ >= 32, and none sits in a homopolymer.
"""

import csv
import os
import shutil
import unittest

from protoSpaceJAM.protoSpaceJAM import main as pJAM

VARIANT_SET = "GRCh38_KOLF2.1J_hap1"

_variant_set_manifest = os.path.join(
    "protoSpaceJAM", "genome_files", "variant_sets", VARIANT_SET, "manifest.json"
)

CASES = [
    dict(
        entry=63, gene='RABAC1-201', ensembl_id='ENST00000222008', terminus='N',
        chrom='19', gRNA_start=41959299, gRNA_end=41959277, strand='-',
        reference_protospacer='CGCAGACATGGCAGCGCAGA',
        variant_protospacer='CGCAAACATGGCAGCGCAGA',
        pam='AGG',
        pos=41959295, offset=4, region='distal',
        ref_base='C', variant_base='T',
        vcf_ref='C', vcf_alt='T', zygosity='het',
        gt='1|0', gq=36, dp=32,
        weight=0.75, same_guide_in_reference_run=True,
    ),
    dict(
        entry=114, gene='NEK9-201', ensembl_id='ENST00000238616', terminus='C',
        chrom='14', gRNA_start=75084548, gRNA_end=75084570, strand='+',
        reference_protospacer='CTACAGGCTCAGGAGACTAG',
        variant_protospacer='CTATAGGCTCAGGAGACTAG',
        pam='AGG',
        pos=75084551, offset=3, region='distal',
        ref_base='C', variant_base='T',
        vcf_ref='C', vcf_alt='T', zygosity='hom',
        gt='1/1', gq=58, dp=61,
        weight=0.9, same_guide_in_reference_run=True,
    ),
    dict(
        entry=146, gene='RABGAP1L-201', ensembl_id='ENST00000251507', terminus='N',
        chrom='1', gRNA_start=174219142, gRNA_end=174219164, strand='+',
        reference_protospacer='AGTTTGCAGAACTGAAATGG',
        variant_protospacer='AGTTTACAGAACTGAAATGG',
        pam='AGG',
        pos=174219147, offset=5, region='distal',
        ref_base='G', variant_base='A',
        vcf_ref='G', vcf_alt='A', zygosity='hom',
        gt='1/1', gq=56, dp=81,
        weight=0.9, same_guide_in_reference_run=True,
    ),
    dict(
        entry=191, gene='ACTR2-201', ensembl_id='ENST00000260641', terminus='C',
        chrom='2', gRNA_start=65268760, gRNA_end=65268738, strand='-',
        reference_protospacer='GGTATGACGGGAACAAGCTT',
        variant_protospacer='GGTATGATGGGAACAAGCTT',
        pam='TGG',
        pos=65268753, offset=7, region='distal',
        ref_base='G', variant_base='A',
        vcf_ref='G', vcf_alt='A', zygosity='hom',
        gt='1/1', gq=54, dp=73,
        weight=0.9, same_guide_in_reference_run=True,
    ),
    dict(
        entry=273, gene='BLVRA-201', ensembl_id='ENST00000265523', terminus='N',
        chrom='7', gRNA_start=43771176, gRNA_end=43771154, strand='-',
        reference_protospacer='ACTCACCTCTGCATTCATCT',
        variant_protospacer='ACTCACCTCTGTATTCATCT',
        pam='TGG',
        pos=43771165, offset=11, region='seed',
        ref_base='G', variant_base='A',
        vcf_ref='G', vcf_alt='A', zygosity='hom',
        gt='1/1', gq=52, dp=67,
        weight=0.5, same_guide_in_reference_run=True,
    ),
    dict(
        entry=293, gene='MAP1LC3B-201', ensembl_id='ENST00000268607', terminus='N',
        chrom='16', gRNA_start=87392415, gRNA_end=87392437, strand='+',
        reference_protospacer='AGATCCCTGCACCATGCCGT',
        variant_protospacer='AGATCCCCGCACCATGCCGT',
        pam='CGG',
        pos=87392422, offset=7, region='distal',
        ref_base='T', variant_base='C',
        vcf_ref='T', vcf_alt='C', zygosity='hom',
        gt='1/1', gq=53, dp=48,
        weight=0.9, same_guide_in_reference_run=True,
    ),
    dict(
        entry=866, gene='ARHGEF7-204', ensembl_id='ENST00000375736', terminus='N',
        chrom='13', gRNA_start=111217708, gRNA_end=111217686, strand='-',
        reference_protospacer='CAGTTGATTGTTGCTATTAT',
        variant_protospacer='CAGTTGATTGTTGCTATTGT',
        pam='CGG',
        pos=111217690, offset=18, region='seed',
        ref_base='T', variant_base='C',
        vcf_ref='T', vcf_alt='C', zygosity='hom',
        gt='1/1', gq=55, dp=39,
        weight=0.5, same_guide_in_reference_run=True,
    ),
    dict(
        entry=966, gene='CDK1-203', ensembl_id='ENST00000395284', terminus='C',
        chrom='10', gRNA_start=60794016, gRNA_end=60793994, strand='-',
        reference_protospacer='CTATCTGTTGATATAACATA',
        variant_protospacer='CTATCTGTTGACATAACATA',
        pam='TGG',
        pos=60794005, offset=11, region='seed',
        ref_base='A', variant_base='G',
        vcf_ref='A', vcf_alt='G', zygosity='hom',
        gt='1/1', gq=57, dp=76,
        weight=0.5, same_guide_in_reference_run=True,
    ),
    dict(
        entry=1059, gene='ULK3-201', ensembl_id='ENST00000440863', terminus='N',
        chrom='15', gRNA_start=74843110, gRNA_end=74843088, strand='-',
        reference_protospacer='CCGGAATGGCGGGGCCCGGC',
        variant_protospacer='CCGGGATGGCGGGGCCCGGC',
        pam='TGG',
        pos=74843106, offset=4, region='distal',
        ref_base='T', variant_base='C',
        vcf_ref='T', vcf_alt='C', zygosity='het',
        gt='1|0', gq=32, dp=45,
        weight=0.75, same_guide_in_reference_run=False,
    ),
]


# The OpenCell parameters, verbatim from run_full_test_variant.py.  Guide choice has to be
# reproduced exactly or the expected sequences do not apply, and guide ranking is sensitive to
# these.  num_gRNA_per_design is left at its default of 1, as the OpenCell run had it: a penalized
# guide only surfaces at rank 1 when there was no better alternative, which is what makes these
# nine the interesting ones.
_OPENCELL_ARGS = {
    "path2csv": os.path.join("input", "test_input_manual_eval.csv"),
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


def _rows_by_design(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    return {(r["ID"], r["terminus"]): r for r in rows}


@unittest.skipUnless(
    os.path.isfile(_variant_set_manifest),
    f"variant set {VARIANT_SET} not built (see precompute/variant_sets/4_make_variant_genome.py)",
)
class test_manual_eval_set(unittest.TestCase):

    # Two runs, not three.  There is deliberately no --variant_genome run here: runtime-patch
    # vs materialized-genome equivalence is asserted over these very designs by
    # run_full_test_variant.py (they are a subset of OpenCell, diffed byte-identically), over six
    # other loci by run_quick_test_variant.py, and over 1e6 windows by
    # precompute/variant_sets/6_validate_runtime_patch.py.  A fourth assertion of the same
    # property on nine designs buys nothing and costs a third of this test's runtime.
    variant_out = os.path.join("tests", "manual_eval_variant")
    reference_out = os.path.join("tests", "manual_eval_reference")

    @classmethod
    def setUpClass(cls):
        print("running the manual eval set (9 designs) against " + VARIANT_SET + "...")
        try:
            os.chdir("protoSpaceJAM")
        except FileNotFoundError:
            raise FileNotFoundError("protoSpaceJAM directory not found, are you in the repo's root directory?")
        if not os.path.isfile("protoSpaceJAM.py"):
            raise FileNotFoundError("protoSpaceJAM.py not found, are you in the protoSpaceJAM directory?")

        for d in (cls.variant_out, cls.reference_out):
            if os.path.exists(d):
                shutil.rmtree(d)

        pJAM(_args(cls.variant_out, variant_set=VARIANT_SET))
        pJAM(_args(cls.reference_out))          # the same nine designs with no variant flag

        cls.variant = _rows_by_design(os.path.join(cls.variant_out, "result.csv"))
        cls.reference = _rows_by_design(os.path.join(cls.reference_out, "result.csv"))

    def _case_row(self, case):
        key = (case["ensembl_id"], case["terminus"])
        self.assertIn(key, self.variant, f"{case['gene']} produced no design")
        return self.variant[key]

    def test_all_nine_designs_are_designable(self):
        self.assertEqual(len(CASES), 9)
        for case in CASES:
            row = self._case_row(case)
            self.assertTrue(row["gRNA_seq"], f"{case['gene']} returned no guide")

    def test_personalized_guide_sequence(self):
        """The emitted protospacer+PAM is the cell line's, not the reference's."""
        for case in CASES:
            row = self._case_row(case)
            with self.subTest(gene=case["gene"]):
                self.assertEqual(row["gRNA_seq"], case["variant_protospacer"])
                self.assertEqual(row["PAM"], case["pam"])
                self.assertEqual(row["gRNA_seq_ref"], case["reference_protospacer"])

    def test_exactly_one_base_differs_at_the_expected_offset(self):
        """
        Re-derive the diff here instead of trusting the reported pair.

        Two different frames meet in this assertion, which is the whole point of making it:

        - `offset` is 0-based along the protospacer in *guide* orientation, so 0 is PAM-distal
          and 19 abuts the PAM.  Getting the orientation wrong on a minus-strand guide would
          still yield a one-base diff, just at 19 - offset, so the index pins the orientation.
        - `ref_base` / `variant_base` are the *genomic* (plus-strand) bases, straight off the VCF.
          On a minus-strand guide the protospacer carries their complements.  Six of the nine
          guides are minus-strand, so an implementation that forgot to complement would pass the
          index check and fail here.
        """
        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        for case in CASES:
            row = self._case_row(case)
            ref, var = row["gRNA_seq_ref"], row["gRNA_seq"]
            expect_ref, expect_var = case["ref_base"], case["variant_base"]
            if case["strand"] == "-":
                expect_ref, expect_var = complement[expect_ref], complement[expect_var]
            with self.subTest(gene=case["gene"]):
                self.assertEqual(len(ref), len(var))
                diffs = [i for i, (a, b) in enumerate(zip(ref, var)) if a != b]
                self.assertEqual(diffs, [case["offset"]],
                                 f"{case['gene']}: expected one diff at {case['offset']}")
                self.assertEqual(ref[case["offset"]], expect_ref)
                self.assertEqual(var[case["offset"]], expect_var)

    def test_guide_coordinates_and_strand(self):
        for case in CASES:
            row = self._case_row(case)
            with self.subTest(gene=case["gene"]):
                self.assertEqual(int(float(row["gRNA_start"])), case["gRNA_start"])
                self.assertEqual(int(float(row["gRNA_end"])), case["gRNA_end"])
                self.assertEqual(row["chr"], case["chrom"])

    def test_variant_weight_and_region(self):
        """Weight follows region x zygosity: distal_hom 0.9, distal_het 0.75, seed_hom 0.5."""
        for case in CASES:
            row = self._case_row(case)
            with self.subTest(gene=case["gene"]):
                self.assertAlmostEqual(float(row["variant_weight"]), case["weight"], places=6)
                flags = set(f for f in row["variant_warnings"].split(";") if f)
                self.assertIn(f"{case['region']}_{case['zygosity']}", flags)
                self.assertIn("gRNA_seq_differs_from_reference", flags)

    def test_the_variant_is_reported_against_the_guide(self):
        for case in CASES:
            row = self._case_row(case)
            expected = (f"{case['chrom']}:{case['pos']}"
                        f"{case['vcf_ref']}>{case['vcf_alt']}({case['zygosity']})")
            with self.subTest(gene=case["gene"]):
                self.assertEqual(row["variants_in_protospacer_PAM"], expected)
                self.assertEqual(row["variant_source"], VARIANT_SET)
                self.assertEqual(row["variant_haplotype"], "1")
                self.assertEqual(row["spec_score_stale"], "True",
                                 "the off-target score is still the reference one")

    def test_no_variant_landed_in_a_pam_in_this_set(self):
        """Documented property of the set; if it ever breaks, the README is wrong too."""
        for case in CASES:
            row = self._case_row(case)
            with self.subTest(gene=case["gene"]):
                self.assertNotIn("pam_hom", row["variant_warnings"])
                self.assertNotIn("pam_het", row["variant_warnings"])

    def test_reference_mode_is_unaffected(self):
        """
        Negative control: the same nine designs with no variant flag.

        Eight keep the same guide, so the reference run must emit the *reference* protospacer --
        proving the difference above comes from cell line mode and not from something else that
        changed.  Entry 1059 (ULK3-201) is the one design where cell line mode picked a different
        guide, so there the two runs must disagree.
        """
        for case in CASES:
            key = (case["ensembl_id"], case["terminus"])
            self.assertIn(key, self.reference, f"{case['gene']} missing from the reference run")
            ref_row = self.reference[key]
            with self.subTest(gene=case["gene"]):
                self.assertFalse(ref_row.get("variant_source"),
                                 "reference mode must not emit variant columns")
                if case["same_guide_in_reference_run"]:
                    self.assertEqual(ref_row["gRNA_seq"], case["reference_protospacer"])
                    self.assertEqual(ref_row["PAM"], case["pam"])
                else:
                    self.assertNotEqual(
                        ref_row["gRNA_seq"], case["variant_protospacer"],
                        f"{case['gene']} is the design where cell line mode should have "
                        "changed the guide choice, but both runs picked the same guide")

    def test_variants_report_carries_every_guide_base(self):
        """
        Every guide base reaches variants_report.csv, in the right region, one row each.

        The report spells guide regions `guide_seed` / `guide_distal` (`region` there also takes
        `donor_arm`, since the guide sits inside the homology arms and each of these variants is
        reported once in each frame).
        """
        path = os.path.join(self.variant_out, "variants_report.csv")
        self.assertTrue(os.path.isfile(path), "variants_report.csv was not written")
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
        guide_rows = {(r["chr"], int(r["pos"])): r for r in rows
                      if r["region"].startswith("guide_")}
        self.assertEqual(len(guide_rows), len(CASES),
                         "one guide-region row per design, no more and no fewer")
        for case in CASES:
            key = (case["chrom"], case["pos"])
            with self.subTest(gene=case["gene"]):
                self.assertIn(key, guide_rows)
                r = guide_rows[key]
                self.assertEqual(r["region"], f"guide_{case['region']}")
                self.assertEqual(r["ref"], case["vcf_ref"])
                self.assertEqual(r["alt"], case["vcf_alt"])
                self.assertEqual(r["zygosity"], case["zygosity"])
                self.assertEqual(r["applied"], "True")
                self.assertEqual(r["sequence_reliable"], "True")


if __name__ == "__main__":

    unittest.main()
