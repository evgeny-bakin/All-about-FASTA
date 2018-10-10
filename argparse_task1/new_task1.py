#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 19:45:54 2018

@author: alena
"""
import argparse

def createParser ():
    parser = argparse.ArgumentParser()
    parser.add_argument ('-i', '--input', type=argparse.FileType('r'),
                         help ='path to input file')
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
    file = namespace.input.read()
    
    if namespace.words:
        wordcount = len(file.split())
        if wordcount == 1:
            print('There is only 1 word in your text.')
        elif wordcount == 0:
            print('Sorry, there are no words in your text.')
        else:
            print('There are', wordcount, 'words in your text.')
    if namespace.lines:
        linecount = file.count('\n')
        if linecount == 1:
            print('There is only 1 line in your text.')
        elif linecount == 0:
            print('Sorry, there are no words in your text.')
        else:
            print('There are', linecount, 'lines in your text.')
    if namespace.findmeaning:
        
        file = file.lower()
        if 'love' in file.split():
            print("Your text is about love!")
        else:
            print('Too difficult to understand the meaning of your text.')
        
    if namespace.mostcommon:
        wordcount = len(file.split())
        if wordcount == 1:
            print('There is only 1 word in your text.')
        elif wordcount == 0:
            print('Sorry, there are no words in your text.')
        else:
            from collections import Counter
            file = file.split()
            Counter = Counter(file) 
            mostcom = Counter.most_common(5) 
            print('These are 5 most common words used in your text:')
            for i in mostcom:
                print(i, end = " ")