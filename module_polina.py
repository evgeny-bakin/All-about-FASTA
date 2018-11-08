import Bio
import progressbar
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from progressbar import *
from time import sleep

def complement_reverse_sequence(input_file, output_file, file_type):
    
    with open(input_file, 'r'):
        if file_type == "fasta" or "fastq":
            for seq_record in SeqIO.parse(input_file, file_type):
                sequence = str(repr(input_file))
    
                complement_sequence = Seq(sequence, IUPAC.unambiguous_dna).complement()
                reverse_complement_sequence = Seq(sequence, IUPAC.unambiguous_dna).reverse_complement()
        
                print('Complement sequence is: \n', complement_sequence, file = open(output_file, 'a'))
                print('Reverse complement sequence is: \n', reverse_complement_sequence, file = open(output_file, 'a'))
    
    return


def delete_reads_shorter(input_file, output_file, parameters):
    print('\n Reading your file... \n')
    with open(input_file, 'r'):
        if file_type == 'fastq':
            reads_len = [len(seq_record.seq) for seq_record in SeqIO.parse(input_file, 'fastq')]
            reads_cnt = len(reads_len)
            reads ={}
            for seq_record in SeqIO.parse(input_file, 'fastq'):
                if len(seq_record.seq) in reads:
                    reads[len(seq_record.seq)].append(seq_record.seq)
                else:
                    reads.update({len(seq_record.seq): seq_record.seq})
            print('Found {} reads'.format(reads_cnt))

            print('\n Searching reads shorter then {} and writing them to the output file \n'.format(int(parameters)))
            pbar = progressbar.ProgressBar.start()

            print('Deleting reads shorter then {}. This reads will be written to {}'.format(int(parameters), output_file))
            for i in range(reads_cnt):
                time.sleep(0.1)
                pbar.update(i)
                for key in reads:
                    if key < int(parameters):
                        print(reads[key], file = open(output_file, 'a')) #в выходной файл записываем отброшенные ридыи удаляем их из словаря
                        del reads[key]
                print('Reads longer then {} are in {}'.format(int(parameters), input_file))
                SeqIO.write(*reads.values(), input_file, file_type) #оставляем во входном файле только риды, длиннее Х

            pbar.finish()

        else:
            print('Deleting reads shorter than X is possible for fastq files.')
    return
