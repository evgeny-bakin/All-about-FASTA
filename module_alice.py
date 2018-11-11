'''
This module contains "fast_report" function, providing some basic statistics for our FASTA/FASTQC file
'''
import Bio
import statistics
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import progressbar as pb
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
    
def quality_score(input_file, output_file):
   
    widgets = ['Assesment of quality score: ', pb.Percentage(), ' ', 
            pb.Bar(marker=pb.RotatingMarker()), ' ', pb.ETA()]
    counter = 0
        
    with open( input_file, "r"):
        if file_type == "fastq":
            print("Data exploring...")
            qualities = [seq_record.letter_annotations["phred_quality"] for seq_record in SeqIO.parse(input_file, file_type)]
            max_read_length = max([len(q) for q in qualities])
            base_numbers = [i for i in range(1,(max_read_length+1))]
            
            per_base_qualities = np.zeros((max_read_length,), dtype=int)
            per_base_quantities = np.zeros((max_read_length,), dtype=int)
            
            timer_1 = pb.ProgressBar(widgets=widgets, maxval=(len(qualities))).start()
    
            for i in qualities:
                counter+=1
                timer_1.update(counter)
                for j in range(len(i)):
                    per_base_qualities[j]+= i[j]
                    per_base_quantities[j] += 1
            
            timer_1.finish()
            print("Counting the average quality...")
            average_per_base_quality = [int(i) for i in per_base_qualities/per_base_quantities]
            
            print("All calculations are finished. Drawing the plot...")
            plt.bar(base_numbers, average_per_base_quality, color=['green'])
            plt.xlabel("base number")
            plt.ylabel("average quality")
            plt.autoscale(enable=True, axis='both', tight=None)
            plt.savefig('average_quality.png')
            print("The plot has been saved as 'average_quality.png'")
    
        else:
            print("Quality score function is available only for fastq files")
            
    with open(output_file,'w') as output_file_1:
        writer = csv.writer(output_file_1, delimiter='\t',lineterminator='\n',)
        writer.writerow(['base number'] + base_numbers)
        writer.writerow(['average quality'] + average_per_base_quality)
    
    return
    