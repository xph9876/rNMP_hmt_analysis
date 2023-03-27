#!/usr/bin/env python3

from collections import defaultdict
import argparse
import sys

# read chromosome sizes
def read_fai(fr):
    data = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        data[ws[0]] = int(ws[1])
    return data


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Calculate the DNA-seq coverage at each coordinate')
    parser.add_argument('sam',type=argparse.FileType('r'), help='Aligned filtered sam file')
    parser.add_argument('fai',type=argparse.FileType('r'), help='Corresponding fasta index file')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output to file')
    parser.add_argument('--covered_only', action='store_true', help='Only output covered items')
    parser.add_argument('--compact', action='store_true', help='Output in compacted format, in this way, all chromsomes are treated as linear')
    parser.add_argument('--max_insert_size', default=1000, help='Maximum insert size allowed (1000nt)')
    args = parser.parse_args()

    # get chrom sizes
    sizes = read_fai(args.fai)

    # Use the algorithom for "merge interval" question
    # Initiate each chromosome and strand
    data = {}
    for k, v in sizes.items():
        for st in '+-':
            data[(k, st)] = defaultdict(int)
    # iteration
    for l in args.sam:
        # Skip header
        if l[0] == '@':
            continue
        # fetch informations
        ws = l.rstrip('\n').split('\t')
        flag = int(ws[1])
        chrom = ws[2]
        loc = int(ws[3])-1
        tlen = int(ws[8])
        # skip second read
        if tlen <= 0:
            continue
        # filter on the read length
        tlen %= sizes[chrom]
        if tlen > args.max_insert_size:
            continue
        # determine strand
        # first read
        if flag % 128 // 64:
            if flag % 32 // 16:
                st = '+'
            else:
                st = '-'
        elif flag % 32 // 16:
            st = '-'
        else:
            st = '+'
        data[(chrom, st)][loc] += 1
        data[(chrom, st)][loc+tlen] -= 1

    # calc coverage
    coverages = defaultdict(lambda: defaultdict(int))
    for chrom, st in data.keys():
        size = sizes[chrom]
        curr = 0
        freq = 0
        if args.compact:
            for loc in sorted(data[(chrom, st)].keys()):
                if loc >= size :
                    continue
                if freq:
                    args.o.write(f'{chrom}\t{curr}\t{loc}\tTrue\t{freq}\t{st}\n')
                elif not args.covered_only:
                    args.o.write(f'{chrom}\t{curr}\t{loc}\tFalse\t{freq}\t{st}\n')
                curr = loc
                freq += data[(chrom, st)][loc]
        else:
            for loc in sorted(data[(chrom, st)].keys()):
                for i in range(curr, loc):
                    coverages[(chrom, st)][i % size] += freq
                curr = loc
                freq += data[(chrom, st)][loc]
            for i in range(size):
                freq = coverages[(chrom, st)][i]
                if freq:
                    args.o.write(f'{chrom}\t{i}\t{i+1}\tTrue\t{freq}\t{st}\n')
                elif not args.covered_only:
                    args.o.write(f'{chrom}\t{i}\t{i+1}\tFalse\t{freq}\t{st}\n')

    print('Done!')

if __name__ == '__main__':
    main()
