import unittest
import os
import module_polina
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator


fasta_test =  os.path.join(os.path.dirname(__file__), 'fasta_test.fasta')
fastq_test = os.path.join(os.path.dirname(__file__), 'fastq_test.fastq')
fasta_res =  os.path.join(os.path.dirname(__file__), 'fasta_res')
fastq_res = os.path.join(os.path.dirname(__file__), 'fastq_res')

class TestFiltering(unittest.TestCase):

    def test_fasta_min_length(file):
        module_polina.min_length(fasta_test, 1665, fasta_res, "fasta")
        counter = 0
        with open(fasta_res) as f:
            counter = sum(1 for _ in f)
        assert(counter/2 == 6)

    def test_fastq_min_length(file):
        module_polina.min_length(fastq_test, 113, fastq_res, "fastq")
        counter = 0
        for title, seq, qual in FastqGeneralIterator(open(fastq_res)):
            counter += 1
        assert(counter == 8)


if __name__ == '__main__':
    unittest.main()