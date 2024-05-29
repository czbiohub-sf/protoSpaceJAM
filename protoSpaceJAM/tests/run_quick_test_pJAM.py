import unittest
import os
import shutil

from protoSpaceJAM.protoSpaceJAM import main as pJAM

class test_with_OpenCell_design(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("running protoSpaceJAM tests using only four designs (quick test)...")
        #change directory to protoSpaceJAM
        try:
            os.chdir(os.path.join('protoSpaceJAM'))
        except FileNotFoundError:
            raise FileNotFoundError("protoSpaceJAM directory not found, are you in the repo's root directory?")

        #check current directiory
        if not os.path.isfile("protoSpaceJAM.py"):
            raise FileNotFoundError("protoSpaceJAM.py not found in current directory, are you in the protoSpaceJAM directory?")

        #Remove previous results
        if os.path.exists(os.path.join("tests", "quick_test_result")):
            shutil.rmtree(os.path.join("tests", "quick_test_result"))

        #run protoSpaceJAM
        test_args = {
                    "path2csv": os.path.join("input", "test_input.csv"),
                    "outdir": os.path.join("tests","quick_test_result"),
                    "ssODN_max_size": 200,
                    "Npayload": "ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT",
                    "Cpayload": "GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
                    "Strand_choice": "NonTargetStrand",
                    "recode_order": "PAM_first",
                    "test_mode": True}
        pJAM(test_args)

    @classmethod
    def tearDownClass(cls):
        keep_results = True
        #Remove results
        if (not keep_results) and os.path.exists(os.path.join("tests", "quick_test_result")):
            shutil.rmtree(os.path.join("tests", "quick_test_result"))
        #pass    #keep results


    def test_if_generated_results(self):
        #check if results were generated
        self.assertTrue(os.path.isfile(os.path.join("tests", "quick_test_result", "result.csv")))

    def test_compare_results(self):
        print("comparing results...")

        #define path to results
        ExpectedResPath = os.path.join("tests","GroundTruths", "quick_result.csv")
        NewResPath = os.path.join("tests","quick_test_result","result.csv")

        #compare results
        with open(ExpectedResPath, 'r') as file1, open(NewResPath, 'r') as file2:
            content1 = file1.read().strip().replace('\r\n', '\n')
            content2 = file2.read().strip().replace('\r\n', '\n')
        #The alterntive: self.assertTrue(filecmp.cmp(NewResPath, ExpectedResPath)) won't work b/c of differences in newlines vs carriage returns

        self.assertEqual(content1, content2)

if __name__ == '__main__':


    unittest.main()