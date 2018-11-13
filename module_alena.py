#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 18:08:22 2018

@author: alena
"""

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio import SeqIO
import csv
from Bio import SeqIO
import time
import progressbar
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def gc_content_analysis (input_file, output_file):
    if file_type == "fasta" or file_type == "fastq":
        dict = {}
        sequence = str(input_file)
        analysis = output_file
        my_seq = Seq(sequence, IUPAC.unambiguous_dna)
        gc_content = GC(my_seq)
        dict['GC-content'] = gc_content
        print('GC content of analysed sequence is', round(gc_content, 2), '%')
        if gc_content < 50:
            return
        else:
            square = Seq(sequence).count('GG')
            dict['Number of "GG" motifs'] = square
            if square >= 4:
                print('G-quadruplex structures can be found in analysed sequence!')
                dict['3D structures occurence'] = 'TRUE'
            else:
                dict['3D structures occurence'] = 'FALSE'
                return
        with open(analysis, 'w') as csvFile:
                result = csv.writer(csvFile, delimiter = ',')
                for key, value in dict.items():
                    result.writerow([key, value])


def join_sequences_seqio(input_file, parameters, output_file):
    file_type = 'fasta'
    print('Reading your first file...')
    print()
    result = tuple()
    with open(input_file, 'r'):
        result = [seq_record for seq_record in SeqIO.parse(input_file, file_type) if seq_record.description not in result]
        print('Reading your second file...')
    with open(parameters, 'r'):
        result = [seq_record for seq_record in SeqIO.parse(parameters, file_type) if seq_record.description not in result]
    with open(output_file, 'w'):
        bar = progressbar.ProgressBar().start()
        for i in range(len(result)):
            time.sleep(0.0)
            bar.update(i)
            bar.finish()
            print()
        print('Writing joined files to a', output_file, 'file...')
        SeqIO.write(result, output_file, file_type)
    return
