'''
This module contains "basic_statistics" function, providing some basic statistics for our FASTA/FASTQ file and "quality_score" function, performing assesment of average per base quality
'''
import Bio
import statistics
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import progressbar as pb
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
import time


def basic_statistics(input_file, output_mode, output_file, file_type):

    with open(input_file, "r"):
        read_length = [len(seq_record.seq) for seq_record in SeqIO.parse(input_file, file_type)]
        count_reads = len(read_length)
        total_length = sum(read_length)
        read_mean = statistics.mean(read_length)
        read_mode = statistics.mode(read_length)
        mode_percent = (read_length.count(read_mode) / count_reads) * 100

        if count_reads == 1:
            print("Length:", read_length, file = open(output_file, "w"))
            if output_mode == "verbose":
                print("There is one sequence in your file. The sequence length is {}".format(read_length[0]))
        else:
            print("Reads:", count_reads, "\n"
                  "Min length:", min(read_length), "\n"
                  "Max length:", max(read_length), "\n"
                  "Average length:", round(read_mean, 2), "\n"
                  "Prevailing length:", read_mode, "({}%)".format(round(mode_percent, 1)), "\n"
                  "Total length:", total_length,
                   file = open(output_file, "w"))
            if output_mode == "verbose":
                print("The file contains {} reads. The read lengths range from {} to {}. The most frequent length is {} ({}% of all reads). The total length is {}.".format(count_reads,
                min(read_length), max(read_length), statistics.mode(read_length), round(mode_percent, 1), total_length))

    return


def quality_score(input_file, phred, output_file, file_type):
    start_time = time.clock()
    widgets = ['Assesment of quality score: ', pb.Percentage(), ' ',
               pb.Bar(marker=pb.RotatingMarker()), ' ', pb.ETA()]
    counter = 0

    with open(input_file, "r"):
        if file_type == "fastq":
            print("Data exploring...")

            if phred == "0" or "33":  #0 goes as an argument if phred hasn't been specified
                encoding = "phred33"
                num = 33
            elif phred == "64":
                encoding = "phred64"
                num = 64
            else:
                print("Unsupportable encoding")

            print("Encoding: {}".format(encoding))

            qualities = []
            for _, _, qual in FastqGeneralIterator(open(input_file)):
                    qualities.append([ord(sym) - num for sym in qual])
            qualities = tuple(qualities)
            max_read_length = max([len(q) for q in qualities])
            base_numbers = [i for i in range(1, (max_read_length + 1))]

            per_base_qualities = np.zeros((max_read_length,), dtype=int)
            per_base_quantities = np.zeros((max_read_length,), dtype=int)
            add_one = np.ones((max_read_length,), dtype=int)

            timer_1 = pb.ProgressBar(widgets=widgets, maxval=(len(qualities))).start()

            for i in qualities:
                counter += 1
                timer_1.update(counter)
                i = np.asarray(i)
                per_base_qualities[:len(i)] += i
                per_base_quantities[:len(i)] += add_one[:len(i)]

            timer_1.finish()
            print("Counting the average quality...")
            average_per_base_quality = per_base_qualities / per_base_quantities
            print(time.clock() - start_time, "seconds_calculations")

            print("All calculations are finished. Drawing the plot...")
            plt.bar(base_numbers, average_per_base_quality, color=['green'])
            plt.xlabel("base number")
            plt.ylabel("average quality")
            plt.autoscale(enable = True, axis='both', tight = None)
            plt.savefig('average_quality.png')
            print("The plot has been saved as 'average_quality.png'")

            df_quality = pd.DataFrame(average_per_base_quality, base_numbers, columns=["average quality"])
            df_quality.to_csv(output_file, sep='\t')

        else:
            print("Quality score function is available only for fastq files")


    print("Time:", time.clock() - start_time, "seconds")
    return

# the following function has been written by Alena Kizenko
def gc_content_analysis(input_file, parameters, output_file, file_type):
    total_GC = 0
    total_length = 0
    if file_type == 'fasta':
        if parameters == 'false':
            with open(output_file, 'w', newline='') as csvfile:
                fieldnames = ['Sequence_ID', 'GC-content']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for title, seq in SimpleFastaParser(open(input_file)):
                    gc_content = GC(seq)
                    gc_content = round(gc_content, 2)
                    total_GC += Seq(seq).count('C')
                    total_GC += Seq(seq).count('G')
                    total_length += len(seq)
                    writer.writerow({'Sequence_ID': title, 'GC-content': gc_content})
                average_GC = total_GC / total_length * 100
                average_GC = round(average_GC, 2)
                writer.writerow({'Sequence_ID': 'Average_GC', 'GC-content': average_GC})

        elif parameters == '3d':
            with open(output_file, 'w', newline='') as csvfile:
                fieldnames = ['Sequence_ID', 'GC-content', '3D structures occurence']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for title, seq in SimpleFastaParser(open(input_file)):
                    gc_content = GC(seq)
                    gc_content = round(gc_content, 2)
                    total_GC += Seq(seq).count('C')
                    total_GC += Seq(seq).count('G')
                    total_length += len(seq)
                    if gc_content > 50 and Seq(seq).count('GG') >= 4:
                        writer.writerow(
                            {'Sequence_ID': title, 'GC-content': gc_content, '3D structures occurence': 'TRUE'})
                    else:
                        writer.writerow(
                            {'Sequence_ID': title, 'GC-content': gc_content, '3D structures occurence': 'FALSE'})
                average_GC = total_GC / total_length * 100
                average_GC = round(average_GC, 2)
                writer.writerow({'Sequence_ID': 'Average_GC', 'GC-content': average_GC})


    if file_type == 'fastq':
        if parameters == 'false':
            with open(output_file, 'w', newline='') as csvfile:
                fieldnames = ['Sequence_ID', 'GC-content']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for title, seq, qual in FastqGeneralIterator(open(input_file)):
                    gc_content = GC(seq)
                    gc_content = round(gc_content, 2)
                    total_GC += Seq(seq).count('C')
                    total_GC += Seq(seq).count('G')
                    total_length += len(seq)
                    writer.writerow({'Sequence_ID': title, 'GC-content': gc_content})
                average_GC = total_GC / total_length * 100
                average_GC = round(average_GC, 2)
                writer.writerow({'Sequence_ID': 'Average_GC', 'GC-content': average_GC})
        elif parameters == '3d':
            with open(output_file, 'w', newline='') as csvfile:
                fieldnames = ['Sequence_ID', 'GC-content', '3D_structures_occurence']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for title, seq, qual in FastqGeneralIterator(open(input_file)):
                    gc_content = GC(seq)
                    gc_content = round(gc_content, 2)
                    total_GC += Seq(seq).count('C')
                    total_GC += Seq(seq).count('G')
                    total_length += len(seq)
                    if gc_content > 50 and Seq(seq).count('GG') >= 4:
                        writer.writerow(
                            {'Sequence_ID': title, 'GC-content': gc_content, '3D_structures_occurence': 'TRUE'})
                    else:
                        writer.writerow(
                            {'Sequence_ID': title, 'GC-content': gc_content, '3D_structures_occurence': 'FALSE'})
                average_GC = total_GC / total_length * 100
                average_GC = round(average_GC, 2)
                writer.writerow({'Sequence_ID': 'Average_GC', 'GC-content': average_GC})
    print('GC-content scores of input sequences have been saved to {}'.format(output_file))
    print('Average GC-content of input sequences is {}'.format(average_GC))


# parameters(optional) should include either treshold (minimal contig length used for calculation) - one figure or metric's name(N50/NA50),genome_size and (optionally) treshold  . Example - NA50:100000:150

def n50(input_file, parameters, output_file, file_type):

    metric = "N50"
    treshold = 0

    if parameters != "0":
        parameters = parameters.split(":")
        if len(parameters) == 1:
            treshold = int(parameters[0])
        else:
            metric = parameters[0]
            genome_size = int(parameters[1])
            if len(parameters) == 3:
                treshold = int(parameters[2])


    contig_length = [len(seq_record.seq) for seq_record in SeqIO.parse(input_file, file_type)]
    contig_length = sorted(contig_length)[::-1]
    count_contigs = len(contig_length)
    print("Total number of contigs: {}".format(count_contigs))

    if metric == "NA50":
        total_length = genome_size
    else:
        total_length = sum(contig_length)

    half_length = total_length/2
    summ = 0

    for i in contig_length:
        if i > treshold:
            while summ < half_length:
                summ+=i
            if summ >= half_length:
                print("{}: {}".format(metric,i))
                print("\n", "Assembly quality assessment for {}".format(input_file), "\n"
                "{}: {}".format(metric, i),
                file=open(output_file, "a"))
                break
    print("The report has been saved to {}".format(output_file))




