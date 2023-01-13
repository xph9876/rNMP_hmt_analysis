#!/usr/bin/env python3

import argparse
from collections import defaultdict
import sys

# read gtf
def read_gtf(fr):
    data = {}
    # read data
    for l in fr:
        ws = l.rstrip().split('\t')
        if len(ws) > 8:
            annots = ws[-1].split('; ')
            gene_id = annots[0].rstrip('\"').replace('gene_id \"', '')
            tx_id = 'None'
            gene_name = gene_id
            for k in annots:
                if k.startswith('transcript_id'):
                    tx_id = k.rstrip('\"').replace('transcript_id \"', '')
                elif k.startswith('gene_name'):
                    gene_name = k.rstrip('\"').replace('gene_name \"', '')
            data[gene_id] = [ws[0], ws[3], ws[4], ws[6], ws[2], gene_name, gene_id, tx_id]
    return data

# read order
def read_libinfo(fr, c):
    data = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        name = ws.pop(c)
        data[name] = ws
    return data

# read data
def read_data(fr):
    data = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        data[(ws[0],ws[1])] = ws[2]
    return data

# read mito count
def read_mito_count(fr):
    data = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        if len(ws) != 2:
            continue
        data[ws[0]] = float(ws[1])
    return data

def main():
    parser = argparse.ArgumentParser(description='Add library and cds info to raw count file')
    parser.add_argument('raw', type=argparse.FileType('r'), help='Raw rNMP count for each cds')
    parser.add_argument('info', type=argparse.FileType('r'), help='List of library information')
    parser.add_argument('gtf', type=argparse.FileType('r'), help='GTF annotation')
    parser.add_argument('mito_count', type=argparse.FileType('r'), help='CSV file of mitochondrial rNMP count for each library')
    parser.add_argument('-c', type=int, default=1, help='Column number of library name in libinfo, default=1')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output to file')
    args = parser.parse_args()

    args.c -= 1

    # read
    cds_info = read_gtf(args.gtf)
    libinfo = read_libinfo(args.info, args.c)
    data = read_data(args.raw)
    mito_count = read_mito_count(args.mito_count)

    # output
    args.o.write('Library\tGenotype\tGene_name\tChromosome\tStart\tEnd\tStrand\tLength' + \
            '\tGene_type\tGene_id\tTranscript_id\trNMP_count\tRPB\tPPB\tCount_rank\tRPB_rank\n')
    for lib, v1 in libinfo.items():
        values = []
        for cds, v2 in cds_info.items():
            gene_len = int(v2[2])-int(v2[1])
            if (lib, cds) not in data:
                continue
            count = float(data[(lib, cds)])
            rpb = count/gene_len
            ppb = rpb/mito_count[lib]
            values.append([lib, v1[0], v2[5]] + v2[:4] + [str(gene_len), v2[4], cds, v2[7], str(int(count)),\
                    str(rpb), str(ppb)])
        count_ranks = [x[2] for x in sorted(values, key=lambda x: -float(x[-3]))]
        rpb_ranks = [x[2] for x in sorted(values, key=lambda x: -float(x[-2]))]
        for v in values:
            args.o.write('\t'.join(v) + '\t' + str(count_ranks.index(v[2])+1) + '\t' + str(rpb_ranks.index(v[2])+1) + '\n')

    print('Done!')

if __name__ == '__main__':
    main()
