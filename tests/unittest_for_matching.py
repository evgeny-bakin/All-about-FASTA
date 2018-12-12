import unittest
import os
import module_alena
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser

fasta_test1 = os.path.join(os.path.dirname(__file__), 'fasta_test1.fasta')
fasta_test2 = os.path.join(os.path.dirname(__file__), 'fasta_test2.fasta')
fasta_result =  os.path.join(os.path.dirname(__file__), 'fasta_result.fasta')

class TestMatchingFasta(unittest.TestCase):
    
    def test_join_fasta(self):
        module_alena.join_sequences(fasta_test1, fasta_test2, fasta_result, 'fasta')
        counter = sum(1 for title, seq in SimpleFastaParser(open(fasta_result)))
        self.assertEqual(counter, 13)
        
    def test_overlap_fasta(self):
        module_alena.overlap_sequences(fasta_test1, fasta_test2, fasta_result, 'fasta')
        counter = sum(1 for title, seq in SimpleFastaParser(open(fasta_result)))
        self.assertEqual(counter, 3)
        
    def test_subtract_fasta(self):
        module_alena.subtract_sequences(fasta_test1, fasta_test2, fasta_result, 'fasta')
        counter = sum(1 for title, seq in SimpleFastaParser(open(fasta_result)))
        self.assertEqual(counter, 10)

class TestMatchingFastq(unittest.TestCase):
    
    def test_join_fastq(self):
        pass
    def test_overlap_fastq(self):
        pass
    def test_subtract_fastq(self):
        pass
    
    
if __name__ == '__main__':
    unittest.main()