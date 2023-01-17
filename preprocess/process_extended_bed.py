#!/usr/bin/env python3

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Process rNMP bed files aligned to extended hmtDNA genome')
    parser.add_argument('bed', type=argparse.FileType('r'), help='rNMP bed file')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    parser.add_argument('-l', type=int, default=16569, help='chrM length (16,569 nt)')
    args = parser.parse_args()

    # load all rNMPs
    data = []
    for l in args.bed:
        ws = l.rstrip().split('\t')
        if len(ws) < 6:
            continue
        ws[1] = int(ws[1]) % args.l
        ws[2] = ws[1] + 1
        data.append(ws)
    
    # output
    data.sort(key=lambda x:(x[1], x[5]))
    for d in data:
        d[1] = str(d[1])
        d[2] = str(d[2])
        args.o.write('\t'.join(d) + '\n')


    print('Done!')


if __name__ == '__main__':
    main()

        