#!/usr/bin/env python3

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Flatten bed with frequencies to one line per rNMP')
    parser.add_argument('bed', type=argparse.FileType('r'), help='Input BED file')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output to file')
    parser.add_argument('-c', default=4, type=int, help='Column number for frequency, default=4')
    args = parser.parse_args()

    # read data
    for l in args.bed:
        ws = l.rstrip().split('\t')
        if len(ws) >= 6:
            for i in range(int(ws[args.c-1])):
                args.o.write(l)
    print('Done!')

if __name__ == '__main__':
    main()
