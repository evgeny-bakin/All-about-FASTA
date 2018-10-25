import Bio 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def complement_reverse_sequence(input_file, output_file, *options):
    
    with open(input_file, 'r'):
     #   if file_type == ("fasta" or "fastq"):
        #for seq_record in SeqIO.parse(input_file, check_file_format(args.input_file)):
        sequence = str(repr(input_file))
    
    complement_sequence = Seq(sequence, IUPAC.unambiguous_dna).complement()
    reverse_complement_sequence = Seq(sequence, IUPAC.unambiguous_dna).reverse_complement()
        
    print('Complement sequence is: \n', complement_sequence, file = open(output_file, 'a'))
    print('Reverse complement sequence is: \n', reverse_complement_sequence, file = open(output_file, 'a'))
    
    return