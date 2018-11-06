import Bio
import pandas
import progressbar
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from progressbar import *

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

    with open(input_file, 'r'):

        count_reads = len([len(seq_record.seq) for seq_record in SeqIO.parse(input_file, file_type)])

        for i in range(count_reads + 1):
            widgets = ['Calculate length for each read in file...', Percentage(), ' ',
                       Bar(marker='0', left='[', right=']'),
                       ' ', ETA(), ' ', FileTransferSpeed()]
            pbar1 = ProgressBar(widgets=widgets, maxval=count_reads)
            pbar1.start()
            reads = pd.Series(seq_record.seq for seq_record in SeqIO.parse(input_file, 'fasta'))
            length_of_reads = reads.str.len()
            pbar1.update(i)
        pbar1.finish()

        print('Search reads shorter than {}'.format(parameters), 'and creating otput file')
        for i in range(length_of_reads):
            widgets = ['Search reads shorter than {}'.format(parameters), 'and creating otput file', Percentage(), ' ',
                       Bar(marker='0', left='[', right=']'),
                       ' ', ETA(), ' ', FileTransferSpeed()]
            pbar2 = ProgressBar(widgets=widgets, maxval = count_reads)
            pbar2.start()
            if length_of_reads[i] < int(parameters):
                print(reads[i], file = open(output_file, 'w'))
            pbar2.update(i)
        pbar2.finish()
    return
