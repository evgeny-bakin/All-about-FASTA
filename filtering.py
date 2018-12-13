import Bio
import re
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser


def min_length(input_file, parameters, output_file, file_type):

    try:
        p = int(parameters)
        print('\nReading your file... \n')

        print('Searching reads shorter than {} \n'.format(int(p)))

        if file_type == "fastq":
            handle = open(output_file, "w")
            for title, seq, qual in FastqGeneralIterator(open(input_file)):
                if len(seq) >= int(p):
                    handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            handle.close()

        else:
            handle = open(output_file, "w")
            for title, seq in SimpleFastaParser(open(input_file)):
                if len(seq) >= int(p):
                    handle.write(">%s\n%s\n" % (title, seq))
            handle.close()

        print('Reads longer than {} are written to {} \n'.format(int(p), output_file))

    except ValueError:
        print("Your value of min length is not correct. Expected number")

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
    pattern = re.compile(r'[AaTtCcGg]+')
    if not pattern.fullmatch(parameters):
        raise ValueError('Motif must include only ATGC bases and must be non 0')

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

    print("Searching duplicates \n")

    if file_type == "fastq":
        reads_set = {title for title, seq, qual in FastqGeneralIterator(open(input_file))}
        handle = open(output_file, "w")
        for title, seq, qual in FastqGeneralIterator(open(input_file)):
            if title in reads_set:
                handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                reads_set.remove(title)
        handle.close()

    else:
        reads_set = {title for title, seq in SimpleFastaParser(open(input_file))}
        handle = open(output_file, "w")
        for title, seq in SimpleFastaParser(open(input_file)):
            if title in reads_set:
                handle.write(">%s\n%s\n" % (title, seq))
                reads_set.remove(title)
        handle.close()

    print(f"Data without duplicates are written to {output_file} \n")

    return


def min_quality(input_file, parameters, output_file, file_type):
    cmd = parameters.split(":")
    phred = cmd[2]
    if phred == "phred33":
        n = 33

    elif phred == "phred64":
        n = 64

    else:
        raise ValueError("Your phred is not correct. We work with phred33 or phred64")

    if file_type == "fastq":
        cmd = parameters.split(":")
        try:
            min_quality = float(cmd[0])
            bases_ratio = float(cmd[1])

            handle = open(output_file, "w")
            for title, seq, qual in FastqGeneralIterator(open(input_file)):
                counter = 0
                for sym in qual:
                    if ord(sym) - n >= min_quality:
                        counter += 1
                if counter * 100 / len(qual) >= bases_ratio:
                    handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
            handle.close()

            print(f"Reads with quality better {min_quality} are written to {output_file}")

        except ValueError:
            print("Your min quality or bases ratio value is not correct. Numbers expected")

    else:
        print("\n Deleting poor_quality reads available only for fastq")

    return
