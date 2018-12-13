import unittest
import os
import basic_statistics
from basic_statistics import *
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser


fasta_test =  os.path.join(os.path.dirname(__file__), 'fasta_test.fasta')
fastq_test = os.path.join(os.path.dirname(__file__), 'test1000.fastq')
fasta_out =  os.path.join(os.path.dirname(__file__), 'fasta_report')
fastq_out = os.path.join(os.path.dirname(__file__), 'fastq_report')
fasta_stat_out =  os.path.join(os.path.dirname(__file__), 'fasta_stat_output')
fastq_stat_out = os.path.join(os.path.dirname(__file__), 'fastq_stat_output')
qual_graph = os.path.join(os.path.dirname(__file__), 'average_quality.png')

class TestBasicStatisticsFastq(unittest.TestCase):
    def test_fastq_basic_statistics(self):
        basic_statistics.basic_statistics(fastq_test, "verbose", fastq_stat_out, "fastq")
        counter = sum(1 for line in (open(fastq_stat_out)))
        self.assertEqual(counter, 6)

    def test_fastq_quality_score(self):
        basic_statistics.quality_score(fastq_test, "phred33", fastq_out, "fastq")
        counter = sum(1 for row in (open(fastq_out)))
        self.assertEqual(counter, 152)
        self.assertTrue(os.path.exists(qual_graph))

    def test_fastq_gc_content_analysis(self):
        basic_statistics.gc_content_analysis(fastq_test, "3d", fastq_out, "fastq")
        rows = [str(row) for row in (open(fastq_out))]
        counter = len(rows)
        last_row = rows[-1]
        self.assertEqual(counter, 1002)
        self.assertEqual(last_row, "Average_GC,44.19,\n")


    def test_n50(self):
        basic_statistics.n50(fastq_test, "150", fastq_out, "fastq")
        lines = [str(line) for line in (open(fastq_out))]
        last_line = lines[-1]
        self.assertEqual(last_line, "N50: 151\n")

class TestBasicStatisticsFasta(unittest.TestCase):

    def test_fasta_basic_statistics(self):
        basic_statistics.basic_statistics(fasta_test,  "verbose", fasta_stat_out, "fasta")
        counter = sum(1 for line in (open(fasta_stat_out)))
        self.assertEqual(counter, 6)

    def test_fasta_gc_content_analysis(self):
        basic_statistics.gc_content_analysis(fasta_test, "3d", fasta_out, "fasta")
        rows = [str(row) for row in (open(fasta_out))]
        counter = len(rows)
        last_row = rows[-1]
        self.assertEqual(counter, 11)
        self.assertEqual(last_row, "Average_GC,42.22,\n")


    def test_n50(self):
        basic_statistics.n50(fasta_test, "150", fasta_out, "fasta")
        lines = [str(line) for line in (open(fasta_out))]
        last_line = lines[-1]
        counter = sum(1 for title, seq in SimpleFastaParser(open(fasta_out)))
        self.assertEqual(last_line, "N50: 1665\n")

if __name__ == '__main__':
    unittest.main()