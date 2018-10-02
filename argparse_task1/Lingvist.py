import csv
import argparse

parser = argparse.ArgumentParser(description='Welcome to the Lingvist! The Lingvist is a program that can be used to analyse text files in .txt format')
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
parser.add_argument('-fa', '--full_analysis', 
                    help = 'use this option to count all text parameters mentioned above: amount of words, lines, printed characters in the choosen text file and also analyse its length', action = 'store_true')

args = parser.parse_args()
print ('The following text will be analysed:')
text = args.file.read()
print (text)

if args.output:
    report = args.output
else:
    report = 'lingv_report.csv'

wordslist=text.split()
words = len(text.split())
lines = text.count("\n") 
characters = 0
for word in wordslist:
    characters += len(word)
if characters == 0:
    analysis = "File is empty"
if 0 < characters <= 100:
    analysis = "The text is extremely short"
if 100 < characters <= 1000:
    analysis = "The text is short"
if 1000 < characters <= 10000:
    analysis = "The text's length is average"
if characters > 10000:
    analysis = "The text is long"

if args.full_analysis:
    args.words = 'TRUE'
    args.lines = 'TRUE'
    args.characters = 'TRUE'
    args.length = 'TRUE'

with open(report , 'w') as csvFile:
    writer = csv.writer(csvFile)
    if args.words:
        csvWords = [["Words", words]]
        writer.writerows(csvWords)
        print ("Words:", words)
    if args.lines:
        csvLines = [["Lines", lines]]
        writer.writerows(csvLines)
        print("Lines:", lines)
    if args.characters:
        csvCharacters = [["Characters", characters]]
        writer.writerows(csvCharacters)
        print("Characters:", characters)
    if args.length:
        csvLength = [["Length", analysis]]
        writer.writerows(csvLength)
        print(analysis)
        
csvFile.close()

if (args.words or args.lines or args.lines or args.characters or args.length or args.full_analysis) == False:
    print ()
    print ('You have not specified any option. Please, specify, what you would like to count.')
else:
    print()
    print('The results were recorded in the csv file named', report)
    print('Thank you for using the Lingvist!')
    

