#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 19:45:54 2018

@author: alena
"""
import argparse
import csv
import os.path
from collections import Counter
import functools
import operator

def createParser ():
    parser = argparse.ArgumentParser()
    parser.add_argument ('-i', '--input', type=argparse.FileType('r'),
                         help ='path to input file')
    parser.add_argument ('-o', '--output', type=argparse.FileType('w'),
                         help ='path to output file')
    parser.add_argument ('-w', '--words', action = 'store_true',
                         help = 'count number of words in a file')
    parser.add_argument ('-l', '--lines', action = 'store_true',
                         help = 'count number of lines in a file')
    parser.add_argument ('-f', '--findmeaning', action = 'store_true',
                         help = 'find what is text about')
    parser.add_argument ('-m', '--mostcommon', action = 'store_true',
                         help = 'the most common letter')
    return parser
 
if __name__ == '__main__':
    parser = createParser()
    namespace = parser.parse_args()
    if os.path.exists(namespace.input.name):
        file = namespace.input.read()
        file = file.lower()
        dict = {}
        wordcount = len(file.split())
        if namespace.words:
            dict['Number of words'] = wordcount
            print('Number of words:', wordcount)
        if namespace.lines:
            linecount = file.count('\n')
            dict['Number of lines'] = 1
            print('Number of lines:', linecount)
        if namespace.findmeaning:
            meaning = len(list(filter(lambda x: x == 'love', file.split())))
            if meaning > 2:
                dict['Text theme'] = 'love story'
                print("Your text is about love!")
            else:
                dict['Text theme'] = 'undefined'
                print('Too difficult to understand the meaning of your text.')       
        if namespace.mostcommon:
            if wordcount == 0:
                dict['Most common word'] = 'None'
                print('Cannot recognize most common word')
            else:
                lst = []
                [lst.append(i[0]) for i in Counter(file.split()).most_common(5)]
                dict['Most common words'] = lst
                print('These are 5 most common words used in your text:')
                print(list(map(lambda x: x.upper(), lst)))
                print ('The concatenated word:', functools.reduce(operator.add,lst))                
        if namespace.output:
            analysis = namespace.output
        else:
            analysis = 'analysis.csv'
        with open(analysis, 'w') as csvFile:
            result = csv.writer(csvFile, delimiter = ',')
            for key, value in dict.items():
                result.writerow([key, value])
