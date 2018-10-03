import argparse
import csv

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Analysis of text file.')
    
    parser.add_argument('file', type = argparse.FileType('r'), help='add file for analysis ')
    
    parser.add_argument('-l', '--lines', help = 'count lines in file', action = 'store')
    
    parser.add_argument('-lt', '--letters', type = int, help = 'count letters in file', action = 'store')
    
    parser.add_argument('-w', '--words', type = int, help = 'count words in file', action = 'store')
    
    parser.add_argument('-c', '--commas',type = int, help = 'count commas in file', action = 'store')
    
    parser.add_argument('-o', '--output', type = argparse.FileType('w'), help = 'output file with results')
    
    args = parser.parse_args()    
    
    analyse = args.file.read()    
      
    args.lines = analyse.count('\n')
    
    args.words = len(analyse.split())
    
    args.commas = analyse.count(',')
    
    letter = []
    letters_analyse = analyse.split()
    for l in letters_analyse:
        letter.append(len(l))
    args.letters = sum(letter)
                   
    result = args.output
    
    res = [['lines', args.lines], ['words', args.words], ['letters', args.letters],
           ['commas', args.commas]]
    
    with open ('result', 'w', newline = '') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(res)
    csvfile.close()
    
    print(result)
        

        
        
    
    
    
    
    

    
    
                
        
    
    
    
    
    
    