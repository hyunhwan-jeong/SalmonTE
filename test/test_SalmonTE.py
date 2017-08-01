import unittest
import SalmonTE
import logging

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


class MyTestCase(unittest.TestCase):
    def test_collecting_fastq_files_paired_ends(self):

        self.input_info = SalmonTE.collect_FASTQ_files(["../example"])
        self.assertEqual(self.input_info["paired"], True)
        self.assertEqual(self.input_info["num_fastq"], 4)
        #self.assertEqual(True, True)


    def test_collecting_fastq_files_single_end(self):
        self.input_info = SalmonTE.collect_FASTQ_files(["../example/CTRL_1_R1.fastq", "../example/TARDBP_1_R1.fastq"])
        self.assertEqual(self.input_info["paired"], False)
        self.assertEqual(self.input_info["num_fastq"], 2)


if __name__ == '__main__':
    unittest.main()
