#!/usr/bin/env python3

import argparse
import os
import numpy.random as random

def read_fai(fr):
    res = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        res[ws[0]] = int(ws[1])
    return res

def main():
    parser = argparse.ArgumentParser(description='Generate a random GTF containing random elements with the same length of input GTF')
    parser.add_argument('gtf', type=argparse.FileType('r'), help='Input gtf file')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Input fasta index file')
    parser.add_argument('-o', default='.', help='Output to file')
    args = parser.parse_args()

    if not os.path.isdir(args.o):
        os.mkdir(args.o)

    # read chromosome length
    chrom_len = read_fai(args.fai)

    # read data
    i = 0
    for l in args.gtf:
        ws = l.rstrip().split('\t')
        if len(ws) > 8:
            i += 1
            annots = ws[-1].split(';')
            geneid = annots[0].rstrip('\"').replace('gene_id \"', '')
            start = int(ws[3])
            end = int(ws[4])
            strand = ws[6]
            rand_start = random.randint(chrom_len[ws[0]] - end + start)
            while(rand_start >= start and rand_start < end):
                rand_start = random.randint(chrom_len[ws[0]] - end + start)
            with open(f'{args.o}/rd{i}.bed', 'w') as fw:
                fw.write('\t'.join([ws[0], str(rand_start), str(rand_start + end - start),\
                        f'rd{i}', geneid, ws[6]]) + '\n')
    print('Done!')

if __name__ == '__main__':
    main()
