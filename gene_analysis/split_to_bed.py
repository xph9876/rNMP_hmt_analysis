#!/usr/bin/env python3

import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Split a gtf file to one sequence per bed')
    parser.add_argument('gtf', type=argparse.FileType('r'), help='Input gtf file')
    parser.add_argument('-o', default='.', help='Output folder')
    args = parser.parse_args()

    if not os.path.isdir(args.o):
        os.mkdir(args.o)

    # read data
    for l in args.gtf:
        ws = l.rstrip().split('\t')
        if len(ws) > 8:
            annots = ws[-1].split(';')
            geneid = annots[0].rstrip('\"').replace('gene_id \"', '')
            with open(f'{args.o}/{geneid}.bed', 'w') as fw:
                fw.write('\t'.join([ws[0], ws[3], ws[4], geneid, ws[2], ws[6]]) + '\n')
    print('Done!')

if __name__ == '__main__':
    main()
