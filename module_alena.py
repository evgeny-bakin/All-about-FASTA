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


def join_sequences(input_file, parameters, output_file, file_type):
    print('Reading your first file...')
    print()
    with open(input_file, 'r'):
        file1 = tuple(seq_record for seq_record in SeqIO.parse(input_file, file_type))
        file1_set = {seq_record.name for seq_record in file1}
    print('Reading your second file...')
    print()
    with open(parameters, 'r'):
        file2 = tuple(seq_record for seq_record in SeqIO.parse(parameters, file_type))
        file2_set = {seq_record.name for seq_record in file2}
    print('Writing joined files to a', output_file, 'file...')
    second_unique_set = file2_set.difference(file1_set)
    result = list()
    for seq_record in file1:
        if seq_record.name in file1_set:
            result.append(seq_record)
            file1_set.difference_update({seq_record.name})
    for seq_record in file2:
        if seq_record.name in second_unique_set:
            result.append(seq_record)
            second_unique_set.difference_update({seq_record.name})
    with open(output_file, 'w'):
        SeqIO.write(result, output_file, file_type)
    return

def overlap_sequences(input_file, parameters, output_file, file_type):
    print('Reading your first file...')
    print()
    with open(input_file, 'r'):
        file1 = tuple(seq_record for seq_record in SeqIO.parse(input_file, file_type))
        file1_set = {seq_record.name for seq_record in file1}
    print('Reading your second file...')
    print()
    with open(parameters, 'r'):
        file2 = tuple(seq_record for seq_record in SeqIO.parse(parameters, file_type))
        file2_set = {seq_record.name for seq_record in file2}
    print('Writing joined files to a', output_file, 'file...')
    overlapped_set = file2_set.intersection(file1_set)
    result = list()
    for seq_record in file1:
        if seq_record.name in overlapped_set:
            result.append(seq_record)
            overlapped_set.difference_update({seq_record.name})
    with open(output_file, 'w'):
        SeqIO.write(result, output_file, file_type)
    return
def subtract_sequences(input_file, parameters, output_file, file_type):
    print('Reading your first file...')
    print()
    with open(input_file, 'r'):
        file1 = tuple(seq_record for seq_record in SeqIO.parse(input_file, file_type))
        file1_set = {seq_record.name for seq_record in file1}
    print('Reading your second file...')
    print()
    with open(parameters, 'r'):
        file2 = tuple(seq_record for seq_record in SeqIO.parse(parameters, file_type))
        file2_set = {seq_record.name for seq_record in file2}
    print('Writing joined files to a', output_file, 'file...')
    subtracted_set = file2_set.symmetric_difference(file1_set)
    result = list()
    for seq_record in file1:
        if seq_record.name in subtracked_set:
            result.append(seq_record)
            subtracted_set.difference_update({seq_record.name})
    for seq_record in file2:
        if seq_record.name in subtracked_set:
            result.append(seq_record)
            subtracted_set.difference_update({seq_record.name})
    with open(output_file, 'w'):
        SeqIO.write(result, output_file, file_type)
    return
