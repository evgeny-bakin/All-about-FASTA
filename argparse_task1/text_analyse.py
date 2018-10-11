import argparse
import csv
import os.path
from functools import reduce

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Analysis of the text file.')
    
    parser.add_argument('file', type = argparse.FileType('r'), help='add file for analysis ')
    
    parser.add_argument('-l', '--lines', help = 'count lines in file', action = 'store_true')
    
    parser.add_argument('-lt', '--letters', help = 'count letters in file', action = 'store_true')
    
    parser.add_argument('-w', '--words', help = 'count words in file', action = 'store_true')
    
    parser.add_argument('-c', '--commas', help = 'count commas in file', action = 'store_true')
    
    parser.add_argument('-o', '--output', type = str, help = 'output .csv file with results', 
                        default = 'result_of_analysis')
    
    args = parser.parse_args()    
    
    if os.path.exists(args.file.name):
        analyse = args.file.read()      
        
        res = []
        
        if args.lines:
            line = analyse.count('\n')
            res.append(['lines', line])
        
        if args.words:
            word = list(map(lambda analyse: len(analyse.split()), analyse.split()))
            res.append(['words', len(word)])
    
        if args.commas:
            comma = len(list(filter(lambda analyse: ',' in analyse, analyse.split())))
            res.append(['commas', comma])
            
        if args.letters:  
            letter = reduce(lambda x, y: x + y, [len(l) for l in analyse.split()])
            res.append(['letters', letter])
        
        result = args.output
                    
        with open (result, 'w') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(res)
               
        print('The result was written into the file',result)