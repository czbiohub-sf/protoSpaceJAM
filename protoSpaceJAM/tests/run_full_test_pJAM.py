import unittest
import os
import shutil

from protoSpaceJAM.protoSpaceJAM import main as pJAM


class test_with_OpenCell_design(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("running protoSpaceJAM tests using the OpenCell design (full test)...")

        #change directory to protoSpaceJAM
        try:
            os.chdir(os.path.join('protoSpaceJAM'))
        except FileNotFoundError:
            raise FileNotFoundError("protoSpaceJAM directory not found, are you in the repo's root directory?")

        #check current directiory
        if not os.path.isfile("protoSpaceJAM.py"):
            raise FileNotFoundError("protoSpaceJAM.py not found in current directory, are you in the protoSpaceJAM directory?")

        #Remove previous results
        if os.path.exists("tests/OpenCell_design"):
            shutil.rmtree("tests/OpenCell_design")

        #run protoSpaceJAM
        test_args = {
                    "path2csv": os.path.join("input", "OpenCell_protospaceX_design.csv"),
                    "outdir": os.path.join("tests", "OpenCell_design"),
                    "ssODN_max_size": 200,
                    "Npayload": "ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT",
                    "Cpayload": "GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
                    "Strand_choice": "NonTargetStrand",
                    "recode_order": "PAM_first"}
        pJAM(test_args)

    @classmethod
    def tearDownClass(cls):
        #Remove results
        if os.path.exists(os.path.join("tests", "OpenCell_design")):
            shutil.rmtree(os.path.join("tests", "OpenCell_design"))

    def test_if_generated_results(self):
        #check if results were generated
        self.assertTrue(os.path.isfile(os.path.join("tests", "OpenCell_design", "result.csv")))

    def test_compare_results(self):
        print("comparing results...")

        #define path to results
        ExpectedResPath = os.path.join("tests","GroundTruths", "OpenCell_design_full_reocoding_result.csv")
        NewResPath = os.path.join("tests","OpenCell_design","result.csv")

        #compare results
        with open(ExpectedResPath, 'r') as file1, open(NewResPath, 'r') as file2:
            content1 = file1.read().strip().replace('\r\n', '\n')
            content2 = file2.read().strip().replace('\r\n', '\n')

        self.assertEqual(content1, content2)


if __name__ == '__main__':
    unittest.main()