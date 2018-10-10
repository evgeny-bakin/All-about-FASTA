import csv
import argparse
import os.path
import functools

if __name__ == "__main__":
       
    parser = argparse.ArgumentParser(description='Welcome to the Lingvist! The Lingvist is a program that can be used to analyse text files')
    
    parser.add_argument('file', type = argparse.FileType(), 
                        help = 'provide a text file for analysis')
    parser.add_argument('-o', '--output', type = str,
                        help = "specify the csv file name for writing a report (lingv_report.csv by default)")
    parser.add_argument('-wc', '--words', 
                        help = 'counts words in the choosen text file', action = 'store_true')
    parser.add_argument('-lc', '--lines', 
                        help = 'counts lines in the choosen text file', action = 'store_true')
    parser.add_argument('-la', '--length', 
                        help = 'analyses how long the choosen text is', action = 'store_true')
    parser.add_argument('-cc', '--characters', 
                        help = 'counts printed symbols in the choosen text file', action = 'store_true')
    parser.add_argument('-a', '--articles', 
                        help = 'counts articles in the choosen text file', action = 'store_true')
    parser.add_argument('-fa', '--full_analysis', 
                        help = 'use this option to count all text parameters mentioned above: amount of words, lines, printed characters in the choosen text file and also analyse its length', action = 'store_true')
    
    args = parser.parse_args()
    
    if os.path.exists(args.file.name):
        print ('The following text will be analysed:')
        text = args.file.read()
        print (text)
        
        if args.output:
            report = args.output
        else:
            report = 'lingv_report.csv'
        
        wordslist = text.split()
        words = len(wordslist)
        lines = text.count("\n") 
        lengths = list(map(lambda x: len(x), wordslist))
        articles = len(list(filter(lambda x: x == "a" or x == "an" or x == "the", wordslist)))
        max_word_length = functools.reduce(lambda a,b: a if a > b else b, lengths)
        the_longest_words = ",".join([wordslist[a] for a in range(words) if len(wordslist[a]) == max_word_length])
        
        characters = sum( [len(word) for word in wordslist] )
        if characters == 0:
            analysis = "File is empty"
        if 0 < characters <= 100:
            analysis = "The text is extremely short"
        elif 100 < characters <= 1000:
            analysis = "The text is short"
        elif 1000 < characters <= 10000:
            analysis = "The text's length is average"
        else:
            analysis = "The text is long"
        
        if args.full_analysis:
            args.words = args.lines = args.characters = args.length = args.the_longest_words = args.articles = 'TRUE'
            
        with open(report , 'w') as csvFile:
            writer = csv.writer(csvFile)
            csvWords = [["Words", words]]
            csvLines = [["Lines", lines]]
            csvCharacters = [["Characters", characters]]
            csvLength = [["Length", analysis]]
            csvLongest = [["The longest words", the_longest_words]]
            csvArticles = [["Articles", articles]]
           
            if args.words:
                writer.writerows(csvWords)
                print ("Words:", words)
            if args.lines:
                writer.writerows(csvLines)
                print("Lines:", lines)
            if args.characters:
                writer.writerows(csvCharacters)
                print("Characters:", characters)
            if args.articles:
                writer.writerows(csvArticles)
                print("Articles:", articles)
            if args.the_longest_words:
                writer.writerows(csvLongest)
                print("The longest words:", the_longest_words)
            if args.length:
                writer.writerows(csvLength)
                print(analysis)
                
         
        arguments = (args.words, args.lines, args.lines, args.characters, args.length, args.articles, args.the_longest_words, args.full_analysis)       
        if all(arguments) == False:
            print ()
            print ('You have not specified any option. Please, specify, what you would like to count.')
        else:
            print()
            print('The results were recorded in the csv file named', report)
            print('Thank you for using the Lingvist!') 