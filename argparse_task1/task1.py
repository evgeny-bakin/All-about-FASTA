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
parser.add_argument("-he", "--hello",
                    help="how many times word Hello appears in text")
parser.add_argument("-a", "--lettera",
                    help="how many times letter a appears in text")
    
args = parser.parse_args()  
text = args.input.read()
    
if args.words:
        print("Number of words:" + len(text.split())
if args.lines:
        print("Number of lines:" + text.count('\n'))
if args.hello:
        print("Word "Hello" appears" + map((text.count("hello")), text.lower()) + "times")
if args.lettera:
        print("Letter "a" appears" + map(text.count("a"), text.lower()) + "times")
    
    
