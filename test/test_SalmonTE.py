import unittest
import SalmonTE
import logging
import os
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


class MyTestCase(unittest.TestCase):
        #self.assertEqual(True, True)

    def test_collecting_fastq_files_rrcutler(self):
        # This is a test regarding the issue: https://github.com/LiuzLab/SalmonTE/issues/12#issuecomment-477872509
        self.input_info = SalmonTE.collect_FASTQ_files(["example/SRR7290434_R1.fq", "example/SRR7290434_R2.fq"])        
        self.assertEqual(self.input_info["paired"], True)
        self.assertEqual(self.input_info["num_fastq"], 1)


    def test_collecting_fastq_files_single_end(self):
        self.input_info = SalmonTE.collect_FASTQ_files(["example/CTRL_1_R1.fastq", "example/TARDBP_1_R1.fastq"])
        self.assertEqual(self.input_info["paired"], False)
        self.assertEqual(self.input_info["num_fastq"], 2)
        self.assertEqual(set(["CTRL_1_R1.fastq", "TARDBP_1_R1.fastq"]), set(self.input_info["file_list"]))


    def test_run_salmon(self):
        def test_collecting_fastq_files_paired_ends(self):
            self.input_info = SalmonTE.collect_FASTQ_files(["example/*.fastq"])
            self.assertEqual(self.input_info["paired"], True)
            self.assertEqual(self.input_info["num_fastq"], 4)

        test_collecting_fastq_files_paired_ends(self)
        self.input_info["--outpath"] = os.path.join(os.getcwd(), "SalmonTE_output")
        self.input_info["--reference"] = os.path.join(os.getcwd(), "reference/hs")
        self.input_info["--num_threads"] = 4
        self.input_info["--exprtype"] = "TPM"
        SalmonTE.run_salmon(self.input_info)

if __name__ == '__main__':
    unittest.main()
