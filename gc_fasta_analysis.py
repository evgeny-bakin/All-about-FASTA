#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 18:08:22 2018

@author: alena
"""

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import csv

def GC_analysis(args_input_file, args_output_file):
    if file_type == "fasta" or if file_type == "fastq":
        dict = {}
        sequence = str(args_input_file)
        analysis = args_output_file
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
    
