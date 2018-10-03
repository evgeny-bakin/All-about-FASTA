# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# я не понимаю, почему код не работает

import argparse


parser = argparse.ArgumentParser(description="Analyse your text")
parser.add_argument("i", "input", type=argparse.FileType(),
                    help="path to your input file")
parser.add_argument("-o", "--output",
                    help = "path to your output file")
parser.add_argument("-w", "--words",
                    help="number of words")
parser.add_argument("-l", "--lines",
                    help="number of lines")

    
args = parser.parse_args()  
text = args.input.read()
    
if args.words:
        print("Number of words:" + len(text.split())
if args.lines:
        print("Number of lines:" + text.count('\n'))
    
    
