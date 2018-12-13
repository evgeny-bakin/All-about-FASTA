# breakFAST
Development of a tool for FASTA/FASTQ processing.
Bioinformatics institute, Autumn-2018.

## Usage

`python3 ./breakFAST.py -i <input file> -f <function> -p <parameters> -o <output_file>`

## Optional arguments
```bash
 -i, --input        Path to fasta/fastq file
 -f, --function     Name of chosen function
 -p, --parameters   Desired parameters
 -o, --output       Path to output file
 ````

## Module Filtering
This module is developed for quick filtering of fasta/fastq files. 
* **min_length** allows delete reads with length less then X. Good reads will be written into the new file.

    example of usage:
    
    `python3 fasta_tool.py -i <input file> -f min_length -p <X> -o <output file>`
    
    `-i` path to fasta or fastq file with raw data
    
    `-p` desired value of minimum reads length (default 0)
    
    `-o` name of output file where will be written reads with length equal or longer then X (default 0)
    
* **delete_N** allows delete reads which contain unknown nucleotide N. 

    `python3 fasta_tool.py -i <input file> -f delete_N -o <output file>`
    
    `-i` path to fasta or fastq file with raw data
    
    `-o` name of output file where will be written reads without unknown nucleotides
    
* **delete_motif**  allows delete reads which contain a given motif

    `python3 fasta_tool.py -i <input file> -f delete_motif -p <motif> -o <output file>`
    
    `-i` path to fasta or fastq file with raw data
    
    `-p` motif - string with sequence of nucleotides. It can be both lowercase and uppercase characters
    
    `-o` name of output file where will be written reads without a given motif
    
* **deduplicate** allows delete the same reads

    `python3 fasta_tool.py -i <input file> -f deduplicate -o <output file>`
    
    `-i` path to fasta or fastq file with raw data
    
    `-o` name of output file where will be written reads from input files without repites
    
* **min_quality** allows delete poor quality reads

    `python3 fasta_tool.py -i <input file> -f min_quality -p <min_qual:base_ratio:phred> -o <output file>`
    
    `-i` path to fasta or fastq file with raw data
    
    `-p` three necessary parameters passed through "**:**". The first is a value of minimal quality per base. 
    The second is a value of bases ratio with quality more or equal then first number. 
    The third is a value of phred (33 or 64 avaliable).
    e.g. "32:80:phred33", "30:85:phred64" etc. 
    
    `-o` name of output file where will be written reads with good quality
   ## Module Matching
   This module has been developed for quick matching fasta/fastq files.
   
   * **join_sequences** allows to join two fasta/fastq files in one and remove duplicates, in case they occur
   
   `python3 fasta_tool.py -i <first file> -f join_sequences -p <second file> -o <output file>`
   
   `-i` path to first fasta or fastq file
   
   `-p` path to second fasta or fastq file
   
   `-o` path to output file where result set of joined sequences will be stored
   
   * **overlap_sequences** allows to overlap two fasta/fastq files and find sequences, which are the same for both files
   
   `python3 fasta_tool.py -i <first file> -f overlap_sequences -p <second file> -o <output file>`
   
   `-i` path to first fasta or fastq file
   
   `-p` path to second fasta or fastq file
   
   `-o` path to output file where result set of overlapped sequences will be stored
   
   * **subtract_sequences** allows to subtract two fasta/fastq files in one and find unique sequences for each file
   
   `python3 fasta_tool.py -i <first file> -f subtract_sequences -p <second file> -o <output file>`
   
   `-i` path to first fasta or fastq file
   
   `-p` path to second fasta or fastq file
   
   `-o` path to output file where result set of subtracted sequences will be stored
   
   ## Contributors\
   
   This tool was developed by Alena Kizenko\[1], Alice Morshneva\[1] and Polina Pavlova\[1] under supervision of Evgeny Bakin\[2].
   \[1] Saint-Petersburg State University, Saint-Petersburg, Russia
   \[2] Bioinformatics institute, Saint-Petersburg, Russia
   
