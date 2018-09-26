import argparse

import csv

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='My super script.')
    
    parser.add_argument('-n', '--number', help='Number to be squared', metavar='Int', type=int, default=3)

    parser.add_argument('-s', '--string', help='String for length counting', metavar='Str', type=str, default='')

    parser.add_argument('-f', '--flag', help='Make cube or not?', action='store_true')

    args = parser.parse_args()    
    
    print(args.number)
    
    print(args.number**2)
    
    print(args.string)

    print(len(args.string))    
    
    print(args.flag)
    
    if args.flag:
        print(args.number**3)        


