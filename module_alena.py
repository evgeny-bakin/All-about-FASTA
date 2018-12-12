from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser


def join_sequences(input_file, parameters, output_file, file_type):
#the result file will include joined sequences from both files without duplicates
    if file_type == 'fastq':
        print('Reading your first file...')
        print()
        with open(input_file) as handle:
            file1_set = {title for title, seq, qual in FastqGeneralIterator(handle)}
        print('Reading your second file...')
        print()
        with open(input_file) as handle2:
            file2_set = {title for title, seq, qual in FastqGeneralIterator(handle2)}
        print('Joining your files...')
        print()
        second_unique_set = file2_set.difference(file1_set)
        print('Writing joined files to a', output_file, 'file...')
        print()
        with open(output_file, 'w') as out_handle:
            with open(input_file) as in_handle:
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    if title in file1_set:
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        file1_set.difference_update({title})
            with open(parameters) as in_handle2:
                for title, seq, qual in FastqGeneralIterator(in_handle2):
                    if title in second_unique_set:
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        second_unique_set.difference_update({title})
    elif file_type == 'fasta':
        print('Reading your first file...')
        print()
        file1_set = {title for title, seq in SimpleFastaParser(open(input_file))}
        print('Reading your second file...')
        print()
        file2_set = {title for title, seq in SimpleFastaParser(open(parameters))}
        print('Joining your files...')
        print()
        second_unique_set = file2_set.difference(file1_set)
        print('Writing joined files to a', output_file, 'file...')
        print()
        with open(output_file, 'w') as out_handle:
            with open(input_file) as in_handle:
                for title, seq in SimpleFastaParser(in_handle):
                    if title in file1_set:
                        out_handle.write(">%s\n%s\n" % (title, seq))
                        file1_set.difference_update({title})
            with open(parameters) as in_handle2:
                for title, seq in SimpleFastaParser(in_handle2):
                    if title in second_unique_set:
                        out_handle.write(">%s\n%s\n" % (title, seq))
                        second_unique_set.difference_update({title})
    return

def overlap_sequences(input_file, parameters, output_file, file_type):
#the result file will include sequences, which are the same for both files
    if file_type == 'fastq':
        print('Reading your first file...')
        print()
        with open(input_file) as handle:
            file1_set = {title for title, seq, qual in FastqGeneralIterator(handle)}
        print('Reading your second file...')
        print()
        with open(input_file) as handle2:
            file2_set = {title for title, seq, qual in FastqGeneralIterator(handle2)}
        print('Overlapping your files...')
        print()
        overlapped_set = file2_set.intersection(file1_set)
        print('Writing overlapped files to a', output_file, 'file...')
        print()
        with open(output_file, 'w') as out_handle:
            with open(input_file) as in_handle:
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    if title in overlapped_set:
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        overlapped_set.difference_update({title})
            with open(parameters) as in_handle2:
                for title, seq, qual in FastqGeneralIterator(in_handle2):
                    if title in overlapped_set:
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        overlapped_set.difference_update({title})
    elif file_type == 'fasta':
        print('Reading your first file...')
        print()
        file1_set = {title for title, seq in SimpleFastaParser(open(input_file))}
        print('Reading your second file...')
        print()
        file2_set = {title for title, seq in SimpleFastaParser(open(parameters))}
        print('Overlapping your files...')
        print()
        overlapped_set = file2_set.intersection(file1_set)
        print(overlapped_set)
        print('Writing overlapped files to a', output_file, 'file...')
        print()
        with open(output_file, 'w') as out_handle:
            with open(input_file) as in_handle:
                for title, seq in SimpleFastaParser(in_handle):
                    if title in overlapped_set:
                        out_handle.write(">%s\n%s\n" % (title, seq))
                        overlapped_set.difference_update({title})
            with open(parameters) as in_handle2:
                for title, seq in SimpleFastaParser(in_handle2):
                    if title in overlapped_set:
                        out_handle.write(">%s\n%s\n" % (title, seq))
                        overlapped_set.difference_update({title})
    return
def subtract_sequences(input_file, parameters, output_file, file_type):
#the result file will include sequences, which are unique for both files
    if file_type == 'fastq':
        print('Reading your first file...')
        print()
        with open(input_file) as handle:
            file1_set = {title for title, seq, qual in FastqGeneralIterator(handle)}
        print('Reading your second file...')
        print()
        with open(input_file) as handle2:
            file2_set = {title for title, seq, qual in FastqGeneralIterator(handle2)}
        print('Subtracting your files...')
        print()
        subtracted_set = file2_set.symmetric_difference(file1_set)
        print('Writing subtracted files to a', output_file, 'file...')
        print()
        with open(output_file, 'w') as out_handle:
            with open(input_file) as in_handle:
                for title, seq, qual in FastqGeneralIterator(in_handle):
                    if title in subtracted_set:
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        subtracted_set.difference_update({title})
            with open(parameters) as in_handle2:
                for title, seq, qual in FastqGeneralIterator(in_handle2):
                    if title in subtracted_set:
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        subtracted_set.difference_update({title})
    elif file_type == 'fasta':
        print('Reading your first file...')
        print()
        file1_set = {title for title, seq in SimpleFastaParser(open(input_file))}
        print('Reading your second file...')
        print()
        file2_set = {title for title, seq in SimpleFastaParser(open(parameters))}
        print('Subtracting your files...')
        print()
        subtracted_set = file2_set.symmetric_difference(file1_set)
        print('Writing subtracted files to a', output_file, 'file...')
        print()
        with open(output_file, 'w') as out_handle:
            with open(input_file) as in_handle:
                for title, seq in SimpleFastaParser(in_handle):
                    if title in subtracted_set:
                        out_handle.write(">%s\n%s\n" % (title, seq))
                        subtracted_set.difference_update({title})
            with open(parameters) as in_handle2:
                for title, seq in SimpleFastaParser(in_handle2):
                    if title in subtracted_set:
                        out_handle.write(">%s\n%s\n" % (title, seq))
                        subtracted_set.difference_update({title})
    return
