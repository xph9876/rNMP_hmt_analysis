#!/usr/bin/env python3

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Merge the bed files for chrM and control region alignment')
    parser.add_argument('bed_mt', type=argparse.FileType('r'), help='rNMP bed file for linear chrM')
    parser.add_argument('bed_CR', type=argparse.FileType('r'), help='rNMP bed file for control region')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    parser.add_argument('--cr_name', default='chrM_CR', help='Control region sequence name (chrM_CR)')
    parser.add_argument('--mt_name', default='chrM', help='chrM sequence name (chrM)')
    parser.add_argument('--cr_start', default=16024, type=int, help='1-based start coordinate of control region (16,024)')
    parser.add_argument('--chrM_length', default=16569, type=int, help='chrM length (16,569)')
    args = parser.parse_args()

    # load all rNMPs in CR
    data = []
    processed = set()
    for l in args.bed_CR:
        ws = l.rstrip().split('\t')
        if len(ws) < 6:
            continue
        if ws[0] != args.cr_name:
            continue
        loc = (int(ws[1]) + args.cr_start - 1) % args.chrM_length
        st = ws[5]
        rd = ws[3]
        data.append((loc, st, rd))
        processed.add(rd)

    # process rNMPs in chrM
    for l in args.bed_mt:
        ws = l.rstrip().split('\t')
        if len(ws) < 6:
            continue
        if ws[0] != args.mt_name:
            continue
        loc = int(ws[1])
        rd = ws[3]
        st = ws[5]
        if rd in processed:
            continue
        else:
            data.append((loc, st, rd))

    data.sort()
    for d in data:
        args.o.write(f'{args.mt_name}\t{d[0]}\t{d[0]+1}\t{d[2]}\t.\t{d[1]}\n')

    print('Done!')


if __name__ == '__main__':
    main()

        