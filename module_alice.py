'''
This module contains "fast_report" function, providing some basic statistics for our FASTA/FASTQC file
'''
import Bio
import statistics
from Bio import SeqIO

    
def fast_report(input_file, output_file, output_mode = "quiet"):
    with open(input_file, "r"):
        read_length = [len(seq_record.seq) for seq_record in SeqIO.parse(input_file, file_type)]
        count_reads = len(read_length)
        if file_type == "fasta" and count_reads == 1:
            print("Length:", read_length, file = open(output_file, "a"))
            if output_mode == "verbose":
                print("The sequence length is {}".format(read_length[0]))
        elif file_type == "fastq":
            print("Reads:", count_reads,
                  "Min length:", min(read_length),
                  "Max length:", max(read_length),
                  "Prevailing length:", statistics.mode(read_length), 
                   file = open(output_file, "a"))
            if output_mode == "verbose":
                print("The file contains {} reads. The read lengths range from {} to {}. The length of the most reads is {}.".format(count_reads, min(read_length), max(read_length), statistics.mode(read_length)))
        
    return                
