import Bio 
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def complement_reverse_sequence(args.input_file, args.output_file):
    sequence = str(args.input_file)
    complement_sequence = Seq(sequence, IUPAC.unambiguous_dna).complement()
    reverse_complement_sequence = Seq(sequence, IUPAC.unambiguous_dna).reverse_complement()
        
    print('Complement sequence is: \n', complement_sequence, file = open(args.output_file, 'a'))
    print('Reverse complement sequence is: \n', reverse_complement_sequence, file = open(args.output_file, 'a'))
    
    return

        
        
        
