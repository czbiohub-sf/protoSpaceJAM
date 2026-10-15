"""SNP recoding regressions using synthetic donors, without genome downloads."""

import itertools
import random
import unittest
from types import SimpleNamespace

from protoSpaceJAM.util.hdr import HDR_flank


def make_donor(phase=1, payload="A", strand=1, guide_same=True, mode="full", donor_type="ssODN", order="PAM_first", max_size=None, flip=False, noncoding=False, personalized=False, position_type=None, split_codon=False, tag_type="start"):
    """Build independent, reproducible gene-oriented arms and their real genomic coordinates."""
    arm_length = 81 + phase - 1
    rng = random.Random(13)
    sequence = "".join(rng.choice(("GCT", "GAA", "ACC", "TTC", "CAG", "GGT")) for _ in range(2 * arm_length // 3 + 1))[:2 * arm_length]
    codon_start = arm_length - phase + 1
    sequence = sequence[:codon_start] + "GCTGCTGCT" + sequence[codon_start + 9:]
    if personalized:
        sequence = sequence[:codon_start + 2] + "A" + sequence[codon_start + 3:]
    guide_left = arm_length - 10
    pam_start = guide_left + 20 if guide_same else guide_left
    sequence = sequence[:pam_start] + ("AGG" if guide_same else "CCT") + sequence[pam_start + 3:]
    coordinates = [1000 + strand * i for i in range(len(sequence))]
    guide_strand = strand if guide_same else -strand
    guide_start = coordinates[guide_left if guide_same else guide_left + 22]
    cut = guide_start + 16 * guide_strand
    phases = [0 if noncoding or position_type in ("5UTR", "3UTR") else strand * (i % 3 + 1) for i in range(len(sequence))]
    if split_codon:
        phases[arm_length:arm_length + 6] = [strand, 0, 0, 0, 2 * strand, 3 * strand]
    position_types = {"1": {}}
    if position_type:
        interval = tuple(sorted((coordinates[arm_length - 2], coordinates[arm_length + 6])))
        position_types["1"]["synthetic"] = {interval: position_type}
        if "junction" in position_type:
            position_types["1"]["synthetic"][(min(coordinates), max(coordinates))] = "cds"
    variant_ctx = None
    if personalized:
        variant = SimpleNamespace(pos=coordinates[codon_start + 2], end=coordinates[codon_start + 2], ref="T", alt="A", applied=True, makes_sequence_unreliable=False, zygosity="hom", short=lambda: "synthetic variant")
        variant_ctx = SimpleNamespace(variants_in=lambda *args: [variant])
    donor = HDR_flank(
        left_flk_seq=sequence[:arm_length], right_flk_seq=sequence[arm_length:],
        left_flk_coord_lst=[coordinates[0], coordinates[arm_length - 1]],
        right_flk_coord_lst=[coordinates[arm_length], coordinates[-1]],
        left_flk_phases=phases[:arm_length], right_flk_phases=phases[arm_length:],
        type=tag_type, ENST_ID="synthetic", ENST_strand=strand, ENST_chr="1",
        gStart=guide_start, gStrand=guide_strand, InsPos=coordinates[arm_length if tag_type == "stop" else arm_length - 1],
        CutPos=cut, Cut2Ins_dist=cut - coordinates[arm_length - 1],
        payload_type="SNP", tag="", SNP_payload=payload, name="synthetic",
        ssODN_max_size=max_size, loc2posType=position_types,
        recoding_args={"cfdThres": 0.03, "recode_order": order, "recoding_off": mode == "off", "recoding_stop_recut_only": mode == "stop", "recoding_coding_region_only": False},
        Donor_type=donor_type, Strand_choice="NonCodingStrand" if flip else "CodingStrand",
        syn_check_args={"check_enzymes": "", "CustomSeq2Avoid": "", "MinArmLenPostTrim": 1000},
        coordinate_without_ENST=noncoding, variant_ctx=variant_ctx,
    )
    return donor, sequence, arm_length


class SNPProtectedCodonTests(unittest.TestCase):
    def test_phase_and_payload_boundaries(self):
        cases = [
            ("123123", "123123", 1, {6, 7, 8}),
            ("1231", "23123", 1, {3, 4, 5}),
            ("12312", "3123", 1, {3, 4, 5}),
            ("1231", "23123123", 4, set(range(3, 9))),
            ("123", "123123", 3, {3, 4, 5}),
            ("123", "123123", 6, set(range(3, 9))),
            ("2", "3123", 1, {0, 1}),
            ("1231", "2", 1, {3, 4}),
            ("000", "00000", 2, {3, 4}),
            ("123", "000123", 2, {3, 4}),
            ("123", "123", 0, set()),
        ]
        for left, right, length, expected in cases:
            with self.subTest(left=left, right=right, payload_length=length):
                self.assertEqual(HDR_flank._snp_protected_arm_indices(left, right, length), expected)

    def test_split_codons_protect_coding_siblings_across_introns(self):
        cases = [
            ("123", "100023123", 1, {3, 7, 8}),
            ("1231000", "23123", 1, {3, 7, 8}),
            ("12312", "0003123", 4, {3, 4, 5, 6, 7, 8}),
            ("123", "1002003123", 1, {3, 6, 9}),
        ]
        for left, right, length, expected in cases:
            with self.subTest(left=left, right=right):
                self.assertEqual(HDR_flank._snp_protected_arm_indices(left, right, length), expected)

    def test_shared_silent_mutation_preserves_masked_bases_and_case(self):
        donor = HDR_flank.__new__(HDR_flank)
        for protected in ("GCT", "gCt"):
            with self.subTest(protected=protected):
                mutated = donor.get_silent_mutations("GCT" + protected + "GCT", "123777123")
                self.assertEqual(mutated[3:6], protected)
                self.assertNotEqual(mutated[:3].upper(), "GCT")
                self.assertNotEqual(mutated[6:].upper(), "GCT")


class SNPWindowScanTests(unittest.TestCase):
    def test_both_scan_strands_preserve_protected_internal_codon(self):
        sequence = "GCT" * 6 + "GCAGG"
        phases = "123" * 2 + "777" + "123" * 4 + "12"
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                donor = HDR_flank.__new__(HDR_flank)
                donor.left_flk_seq = "A" * 10
                donor.tag = "T"
                donor.cfdThres = 0
                donor.info_p6 = ""
                donor.phase6recoded = False
                donor.Donor_postMut = sequence
                donor.gRNA_seq = donor.revcom(sequence) if reverse else sequence
                scan = donor.slide_win_mutation_noncoding if reverse else donor.slide_win_mutation_coding
                scan(donor.gRNA_seq, phases[::-1] if reverse else phases)
                self.assertEqual(donor.Donor_postMut[6:9].upper(), sequence[6:9])
                self.assertNotEqual(donor.Donor_postMut[:6].upper(), sequence[:6])


class SNPFinalDonorTests(unittest.TestCase):
    def assert_final_donor(self, donor, baseline, boundary, payload, phase=1, mode="full", noncoding=False):
        """Compare delivered nucleotides with an independent reference-plus-payload expectation."""
        expected = baseline[:boundary] + payload + baseline[boundary + len(payload):]
        if noncoding:
            protected = set(range(boundary, boundary + len(payload)))
        else:
            first_codon = boundary - phase + 1
            last_codon_end = ((boundary + len(payload) + 2) // 3) * 3
            protected = set(range(first_codon, last_codon_end))
        final = donor.Donor_final.upper()
        if getattr(donor, "strand_flipped", False):
            final = donor.revcom(final)
        offset = int(donor.centering_start) if donor.ssODN_max_size is not None and donor.Donor_type == "ssODN" else 0
        if donor.Donor_type == "ssODN" and donor.ssODN_max_size is not None:
            self.assertLessEqual(len(final), donor.ssODN_max_size)
            self.assertEqual(len(final), int(donor.centering_end) - offset)
        else:
            self.assertEqual(len(final), len(baseline))
        for index in protected:
            self.assertTrue(offset <= index < offset + len(final), "Protected codon was trimmed out of this fixture")
            self.assertEqual(final[index - offset], expected[index], "Protected nucleotide changed at arm index %s" % index)
        if mode == "off":
            self.assertEqual(final, expected[offset:offset + len(final)])
        return final, protected, offset

    def test_final_donors_across_phases_strands_modes_and_types(self):
        options = itertools.product((1, 2, 3), (1, -1), (True, False), ("full", "stop", "off"), ("ssODN", "dsDNA"), ("PAM_first", "protospacer_first"))
        for phase, strand, guide_same, mode, donor_type, order in options:
            with self.subTest(phase=phase, strand=strand, guide_same=guide_same, mode=mode, donor_type=donor_type, order=order):
                payload = "ATTA" if phase == 2 else "A"
                donor, baseline, boundary = make_donor(phase=phase, payload=payload, strand=strand, guide_same=guide_same, mode=mode, donor_type=donor_type, order=order)
                self.assert_final_donor(donor, baseline, boundary, payload, phase, mode)

    def test_centering_and_strand_flip_preserve_protected_codons(self):
        for strand, phase in itertools.product((1, -1), (1, 2, 3)):
            with self.subTest(strand=strand, phase=phase):
                donor, baseline, boundary = make_donor(phase=phase, payload="ATTA", strand=strand, max_size=100, flip=True)
                self.assertTrue(donor.strand_flipped)
                self.assert_final_donor(donor, baseline, boundary, "ATTA", phase)

    def test_stop_targets_use_the_existing_right_arm_boundary(self):
        for strand, donor_type in itertools.product((1, -1), ("ssODN", "dsDNA")):
            with self.subTest(strand=strand, donor_type=donor_type):
                donor, baseline, boundary = make_donor(phase=2, payload="ATTA", tag_type="stop", strand=strand, donor_type=donor_type)
                self.assert_final_donor(donor, baseline, boundary, "ATTA", phase=2)

    def test_nontranscript_inputs_protect_payload(self):
        for mode, donor_type in itertools.product(("full", "stop", "off"), ("ssODN", "dsDNA")):
            with self.subTest(mode=mode, donor_type=donor_type):
                donor, baseline, boundary = make_donor(payload="ATTA", noncoding=True, mode=mode, donor_type=donor_type)
                self.assert_final_donor(donor, baseline, boundary, "ATTA", mode=mode, noncoding=True)

    def test_coding_annotations_include_protected_bases_on_both_arms(self):
        for donor_type in ("ssODN", "dsDNA"):
            with self.subTest(donor_type=donor_type):
                donor, baseline, boundary = make_donor(phase=2, payload="ATTA", donor_type=donor_type)
                _, protected, _ = self.assert_final_donor(donor, baseline, boundary, "ATTA", phase=2)
                for index in protected:
                    self.assertTrue(any(start <= index < end for start, end in donor.Donor_features["coding_coord"]), "Missing coding annotation at %s" % index)
                    self.assertTrue(any(start <= index < end for start, end in donor.Donor_features["ORF_coord"]), "Missing ORF annotation at %s" % index)

    def test_personalized_arm_is_the_protected_reference(self):
        for donor_type in ("ssODN", "dsDNA"):
            with self.subTest(donor_type=donor_type):
                donor, baseline, boundary = make_donor(personalized=True, donor_type=donor_type)
                self.assertEqual(baseline[boundary + 2], "A")
                self.assert_final_donor(donor, baseline, boundary, "A")
                self.assertEqual(donor.variants_in_donor, "synthetic variant")
                self.assertTrue(donor.Donor_features["variant_coord"])

    def test_annotation_masks_cannot_overwrite_snp_protection(self):
        for position_type, strand in itertools.product(("5UTR", "3UTR", "-3_to_+6bp_of_exon_intron_junction"), (1, -1)):
            with self.subTest(position_type=position_type, strand=strand):
                donor, baseline, boundary = make_donor(payload="ATTA", position_type=position_type, strand=strand)
                self.assert_final_donor(donor, baseline, boundary, "ATTA", noncoding="UTR" in position_type)
                phases = donor.left_flk_phases + donor.right_flk_phases
                for index in donor.snp_protected_indices:
                    self.assertEqual(phases[index], "7")

    def test_final_donor_preserves_split_codon_across_intron(self):
        for strand, donor_type in itertools.product((1, -1), ("ssODN", "dsDNA")):
            with self.subTest(strand=strand, donor_type=donor_type):
                donor, baseline, boundary = make_donor(split_codon=True, strand=strand, donor_type=donor_type)
                self.assertEqual(donor.Donor_final[boundary].upper(), "A")
                self.assertEqual(donor.Donor_final[boundary + 4:boundary + 6].upper(), baseline[boundary + 4:boundary + 6])
                self.assertEqual(donor.snp_protected_indices, {boundary, boundary + 4, boundary + 5})

    def test_invariant_rejects_wrong_protected_base_and_wrong_length(self):
        donor, baseline, boundary = make_donor()
        valid = baseline[:boundary] + "A" + baseline[boundary + 1:]
        donor._validate_snp_protected_sequence(valid)
        outside = "A" if valid[0] != "A" else "C"
        donor._validate_snp_protected_sequence(outside + valid[1:])
        wrong = valid[:boundary + 2] + "A" + valid[boundary + 3:]
        with self.assertRaises(ValueError):
            donor._validate_snp_protected_sequence(wrong)
        with self.assertRaises(ValueError):
            donor._validate_snp_protected_sequence(valid[:-1])


if __name__ == "__main__":
    unittest.main()
