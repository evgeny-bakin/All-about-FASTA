'''
This module contains "fast_report" function, providing some basic statistics for our FASTA/FASTQC file
'''
import Bio
import statistics
import csv
import matplotlib.pyplot as plt
import pandas as pd
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
            max_read_length = max([len(seq_record.seq) for seq_record in SeqIO.parse(input_file, file_type)])
            base_numbers = [i for i in range(1,(max_read_length+1))]
            
            timer_1 = pb.ProgressBar(widgets=widgets, maxval=(max_read_length+1)).start()
            timer_2 = pb.ProgressBar(widgets=widgets, maxval=(len(qualities))).start()
    
            per_base_qualities = {}
            for i in range(1,(max_read_length+1)):
                timer_1.update(i)
                per_base_qualities['base'+str(i)] = []
            timer_1.finish()
            for i in qualities:
                counter+=1
                timer_2.update(counter)
                for j in range(1,(len(i)+1)):
                    per_base_qualities['base'+str(j)].append(i[j-1]) 
            timer_2.finish()
            print("Counting the average quality...")
            average_per_base_quality = [ int(statistics.mean(per_base_qualities[key])) for key in per_base_qualities]
            
            print("All calculations are finished. Drawing the plot...")
            plt.bar(base_numbers, average_per_base_quality, color=['green'])
            plt.xlabel("base number")
            plt.ylabel("average quality")
            plt.xticks(range(0, max_read_length))
            plt.yticks(range(0,65))
            plt.show() 
        else:
            print("Quality score function is available only for fastq files")
    with open(output_file,'w') as output_file_1:
        writer = csv.writer(output_file_1, delimiter='\t',lineterminator='\n',)
        writer.writerow(['base number'] + base_numbers)
        writer.writerow(['average quality'] + average_per_base_quality)
    
    return
    