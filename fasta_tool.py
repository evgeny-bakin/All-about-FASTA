import argparse
import os
import inspect
from Bio import SeqIO
import module_evg

def check_input_file(parser, file_name):
    full_file_name = os.path.abspath(file_name)
    if not os.path.exists(full_file_name):
        parser.error("The file {} does not exist!".format(file_name))
    else:
        return full_file_name
        
def check_file_format(file_name):
    with open(file_name, "r") as handle:
        try:
            for record in SeqIO.parse(handle, "fastq"):
                rec_id = record.id
                return "fastq"
        except:
            pass
        handle.seek(0)
        try:
            for record in SeqIO.parse(handle, "fasta"):
                rec_id = record.id
                return "fasta"
        except:
                pass
        
        print("Unknown file format! Only FASTA/FASTQ files are available.")
        exit()
    return

def make_arguments_parser():
    parser = argparse.ArgumentParser(
        description='Simple tool for basic manipulations with FASTA/FASTQ files.')
        
    parser.add_argument('-i', dest = 'input_file', 
                        help='Input FASTA/FASTQ file', 
                        required=True,
                        type=lambda x: check_input_file(parser, x))

    parser.add_argument('-f', dest = 'function', 
                        help='Function for FASTQ/FASTQ processing', 
                        required=True,
                        type=str)
                        
    parser.add_argument('-p', dest = 'parameters', 
                        help='Parameters of requested function', 
                        required=False,
                        type=str)
                        
    parser.add_argument('-o', dest = 'output_file', 
                        help='Output file name (without extension)', 
                        required=False,
                        type=str,
                        default='')
                                                    
    return parser.parse_args()

if __name__ == "__main__":
    
    args = make_arguments_parser()
    
    file_type = check_file_format(args.input_file)
    print("File type is {}.".format(file_type))
    
    all_functions = [item[0] for item in 
                            inspect.getmembers(module_evg, inspect.isfunction)]
                            
    if args.function not in all_functions:
        print("Function '{}' is unavailable!".format(args.function))
        
    exec('{}({},{},{})'.format(args.function, args.input_file, args.parameters, args.output_file))