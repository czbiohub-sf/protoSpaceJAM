import unittest
import filecmp
import os
import subprocess
import shutil

class test_with_OpenCell_design(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("running protoSpaceJAM tests using only four designs (quick test)...")

        #check if current directiory is
        if not os.path.isfile("main.py"):
            raise FileNotFoundError("main.py not found in current directory, are you in the protoSpaceJAM directory?")

        #Remove previous results
        if os.path.exists("tests/quick_test"):
            shutil.rmtree("tests/quick_test")

        #run protoSpaceJAM
        subprocess.run(["python", "main.py",
                       "--path2csv", "input/test_protoSpaceJAM.csv",
                       "--outdir", "tests/quick_test",
                       "--ssODN_max_size", "200",  "--recoding_full",
                        "--Npayload", "ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT",
                        "--Cpayload", "GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG",
                        "--Strand_choice", "NonTargetStrand",
                        "--recode_order", "PAM_first"])

    @classmethod
    def tearDownClass(cls):
        #Remove results
        if os.path.exists("tests/quick_test"):
            shutil.rmtree("tests/quick_test")

    def test_if_generated_results(self):
        #check if results were generated
        self.assertTrue(os.path.isfile("tests/quick_test/result.csv"))

    def test_compare_results(self):
        print("comparing results...")

        #define path to results
        ExpectedResPath = os.path.join("tests","GroundTruths", "quick_result.csv")
        NewResPath = os.path.join("tests","quick_test","result.csv")

        #compare results
        with open(ExpectedResPath, 'r') as file1, open(NewResPath, 'r') as file2:
            content1 = file1.read().strip().replace('\r\n', '\n')
            content2 = file2.read().strip().replace('\r\n', '\n')
        #The alterntive: self.assertTrue(filecmp.cmp(NewResPath, ExpectedResPath)) won't work b/c of differences in newlines vs carriage returns

        self.assertEqual(content1, content2)

if __name__ == '__main__':
    unittest.main()