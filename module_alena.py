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


from Bio import SeqIO
import time
import progressbar
def join_sequences(input_file, parameters, output_file):
    print('Reading your first file...')
    print()
    with open(input_file, 'r'):
        record1 = {}
        print('These sequences were identified as', file_type,':')
        print()
        for record in SeqIO.parse(input_file, file_type):
            print(record.description)
            record1[record.id] = record
    print('Reading your second file...')
    print()
    with open(parameters, 'r'):
        record2 = {}
        print('These sequences were identified as', file_type,':')
        print()
        for record in SeqIO.parse(parameters, file_type):
            print(record.description)
            record2[record.id] = record
    print()
    print('Joining your files, please, wait...')
    print()
    with open(output_file, 'w'):
        result = list(record1.values())
        bar = progressbar.ProgressBar().start()
        for i in range(len(result)):
            time.sleep(0.1)
            bar.update(i)
            for k in record2:
                if k not in record1:
                    result.append(record2[k])
                else:
                    pass
        bar.finish()
        print()
        print('Writing joined files to a', output_file, 'file...')
        SeqIO.write(result, output_file, file_type)
    return
