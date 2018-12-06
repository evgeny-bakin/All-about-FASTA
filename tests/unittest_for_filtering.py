import unittest
import os
import os,sys,inspect
current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)
import module_polina
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser


fasta_test =  os.path.join(os.path.dirname(__file__), 'fasta_test.fasta')
fastq_test = os.path.join(os.path.dirname(__file__), 'fastq_test.fastq')
fasta_res =  os.path.join(os.path.dirname(__file__), 'fasta_res')
fastq_res = os.path.join(os.path.dirname(__file__), 'fastq_res')

class TestFilteringFasta(unittest.TestCase):

    def test_fasta_min_length(self):
        module_polina.min_length(fasta_test, 1665, fasta_res, "fasta")
        counter = sum(1 for title, seq in SimpleFastaParser(open(fasta_res)))
        self.assertEqual(counter, 6)

    def test_fasta_N(self):
        module_polina.delete_N(fasta_test, 0, fasta_res, "fasta")
        counter = sum(1 for title, seq in SimpleFastaParser(open(fasta_res)))
        self.assertEqual(counter, 5)

    def test_fasta_motif(self):
        module_polina.delete_motif(fasta_test, "aaaag", fasta_res, "fasta")
        counter = sum(1 for title, seq in SimpleFastaParser(open(fasta_res)))
        self.assertEqual(counter, 1)

    def test_fasta_deduplicate(self):
        module_polina.deduplicate(fasta_test, 0, fasta_res, "fasta")
        counter = sum(1 for title, seq in SimpleFastaParser(open(fasta_res)))
        self.assertEqual(counter, 4)

    def test_fasta_min_qual(self):
        module_polina.min_quality(fasta_test, "32:80:phred33", fasta_res, "fasta")
        self.assertFalse(os.path.exists(fasta_res))

    def tearDown(self):
        os.system("rm -r  fasta_res")


class TestFilteringFastq(unittest.TestCase):

    def test_fastq_min_length(self):
        module_polina.min_length(fastq_test, 113, fastq_res, "fastq")
        counter = sum(1 for title, seq, qual in FastqGeneralIterator(open(fastq_res)))
        self.assertEqual(counter, 8)

    def test_fastq_N(self):
        module_polina.delete_N(fastq_test, 0, fastq_res, "fastq")
        counter = sum(1 for title, seq, qual in FastqGeneralIterator(open(fastq_res)))
        self.assertEqual(counter, 6)

    def test_fastq_motif(self):
        module_polina.delete_motif(fastq_test, 'ttttg', fastq_res, "fastq")
        counter = sum(1 for title, seq, qual in FastqGeneralIterator(open(fastq_res)))
        self.assertEqual(counter, 2)

    def test_fastq_deduplicate(self):
        module_polina.deduplicate(fastq_test, 0, fastq_res, "fastq")
        counter = sum(1 for title, seq, qual in FastqGeneralIterator(open(fastq_res)))
        self.assertEqual(counter, 8)

    def test_fastq_qual(self):
        module_polina.min_quality(fastq_test, "32:90:phred33", fastq_res, "fastq")
        counter = sum(1 for title, seq, qual in FastqGeneralIterator(open(fastq_res)))
        self.assertEqual(counter, 6)

    def tearDown(self):
        os.system("rm -r  fastq_res")


if __name__ == '__main__':
    unittest.main()