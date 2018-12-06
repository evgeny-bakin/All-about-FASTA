import Bio
import re
import memory_profiler
from memory_profiler import profile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser


def complement_reverse_sequence(input_file, output_file):
    
    with open(input_file, 'r'):
        if file_type == "fasta" or "fastq":
            for seq_record in SeqIO.parse(input_file, file_type):
                sequence = str(repr(input_file))
    
                complement_sequence = Seq(sequence, IUPAC.unambiguous_dna).complement()
                reverse_complement_sequence = Seq(sequence, IUPAC.unambiguous_dna).reverse_complement()
        
                print('Complement sequence is: \n', complement_sequence, file = open(output_file, 'a'))
                print('Reverse complement sequence is: \n', reverse_complement_sequence, file = open(output_file, 'a'))
    
    return


def delete_reads_shorter_tuple(input_file, parameters, output_file, file_type):

    print('\nReading your file... \n')
    sequence_hadle = SeqIO.parse(input_file, "{}".format(file_type))

    print('Searching reads shorter than {} \n'.format(int(parameters)))
    long_reads = tuple(seq_record for seq_record in sequence_hadle if len(seq_record.seq) >= int(parameters))


    SeqIO.write(long_reads, output_file, "{}".format(file_type))
    print('Reads longer than {} are written to {} \n'.format(int(parameters), output_file))

    return


def min_length(input_file, parameters, output_file, file_type):
    print('\nReading your file... \n')

    print('Searching reads shorter than {} \n'.format(int(parameters)))

    if file_type == "fastq":
        handle = open(output_file, "w")
        for title, seq, qual in FastqGeneralIterator(open(input_file)):
            if len(seq) >= int(parameters):
                handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        handle.close()

    else:
        handle = open(output_file, "w")
        for title, seq in SimpleFastaParser(open(input_file)):
            if len(seq) >= int(parameters):
                handle.write(">%s\n%s\n" % (title, seq))
        handle.close()

    print('Reads longer than {} are written to {} \n'.format(int(parameters), output_file))

    return


def delete_N(input_file, parameters, output_file, file_type):
    print('\nReading your file... \n')

    print('Searching reads containing N \n')

    if file_type == "fastq":
        handle = open(output_file, "w")
        for title, seq, qual in FastqGeneralIterator(open(input_file)):
            if 'N' not in seq.upper():
                handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        handle.close()

    else:
        handle = open(output_file, "w")
        for title, seq in SimpleFastaParser(open(input_file)):
            if 'N' not in seq.upper():
                handle.write(">%s\n%s\n" % (title, seq))
        handle.close()

    print('Reads containing N are written to {} \n'.format(output_file))

    return


def delete_motif(input_file, parameters, output_file, file_type):
    print('\nReading your file... \n')

    print('Searching reads containing motif {} \n'.format(parameters))

    if file_type == "fastq":
        handle = open(output_file, "w")
        for title, seq, qual in FastqGeneralIterator(open(input_file)):
            if not re.findall(r'{}'.format(parameters.upper()), seq.upper()):
                handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        handle.close()

    else:
        handle = open(output_file, "w")
        for title, seq in SimpleFastaParser(open(input_file)):
            if not re.findall(r'{}'.format(parameters.upper()), seq.upper()):
                handle.write(">%s\n%s\n" % (title, seq))
        handle.close()

    print('Reads without motif {} are written to {} \n'.format(parameters, output_file))

    return


def deduplicate(input_file, parameters, output_file, file_type):
    print('Reading your file... \n')

    with open(input_file, 'r'):
        reads = tuple(seq_record for seq_record in SeqIO.parse(input_file, file_type))
        reads_set = {seq_record.name for seq_record in reads}

    print("Searching duplicates \n")

    result = list()
    for seq_record in reads:
        if seq_record.name in reads_set:
            result.append(seq_record)
            reads_set.difference_update({seq_record.name})

    with open(output_file, 'w'):
        SeqIO.write(result, output_file, file_type)

    print(f"Data without duplicates are written to {output_file} \n")

    return


def min_quality(input_file, parameters, output_file, file_type):

    """
    В параметре через двоеточия последовательно передаются - минимальньное качество, доля оснований в процентах, phred,
    e.g. 20:30:phred33
    """

    if file_type == "fastq":

        print('Reading your file... \n')

        print("Searching poor quality reads \n")

        cmd = parameters.split(":")
        min_quality = float(cmd[0])
        bases_ratio = float(cmd[1])
        phred = cmd[2]

        if phred == "phred33":
            n = 33

        elif phred == "phred64":
            n = 64

        else:
            raise ValueError("Your phred is not correct. We work with phred33 or phred64")

        handle = open(output_file, "w")
        for title, seq, qual in FastqGeneralIterator(open(input_file)):
            counter = 0
            for sym in qual:
                if ord(sym) - n >= min_quality:
                    counter += 1
            if counter*100/len(qual) >= bases_ratio:
                handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        handle.close()

        print(f"Reads with quality better {min_quality} are written to {output_file}")

    else:
        print("\n Deleting poor_quality reads available only for fastq")

    return
