#!/usr/bin/env python3

import argparse
from collections import defaultdict
import sys
import pandas as pd

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
            data[gene_id] = [ws[0], int(ws[3]), int(ws[4]), ws[6], ws[2], gene_name, gene_id, tx_id]
    df = pd.DataFrame.from_dict(data, orient='index', columns=['Chromosome', 'Start', 'End', 'Strand', 'Gene_type', 'Gene_name', 'Gene_id', 'Transcript_id'])
    return df

# read order
def read_libinfo(fr, c):
    data = {}
    i = 0
    for l in fr:
        i += 1
        ws = l.rstrip().split('\t')
        name = ws.pop(c)
        data[name] = ['-'.join(ws), i]
    data = pd.DataFrame.from_dict(data, orient='index', columns=['Genotype', 'Order'])
    return data

# read data
def read_data(fr):
    df = pd.read_csv(fr, sep='\t')
    df[['Chrom', 'Library', 'Gene_id', 'st']] = df.Sample.str.split('_', expand=True)
    df = df.drop(columns=['Chrom', 'st', 'Sample'])
    return df

# read mito count
def read_mito_count(fr):
    data = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        if len(ws) != 2:
            continue
        data[ws[0]] = float(ws[1])
    return data

# background count for each gene
def read_mono(fr, s):
    data = {}
    header = fr.readline().rstrip('\n').split('\t')
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        name = ws[0].split(':')[1]
        data[name] = [float(x) for x in ws[1:]] + [s]
    if s == 'nontemplate':
        return pd.DataFrame.from_dict(data, orient='index', columns=['A_bg', 'C_bg', 'G_bg', 'T_bg', 'rNMP_strand'])
    else:
        return pd.DataFrame.from_dict(data, orient='index', columns=['T_bg', 'G_bg', 'C_bg', 'A_bg', 'rNMP_strand'])

# function to generate loc columns in df
def generate_loc(x):
    name = f'{x.Start}-{x.End}({x.Strand})'
    return name


def main():
    parser = argparse.ArgumentParser(description='Add library and cds info to raw count file')
    parser.add_argument('raw', type=argparse.FileType('r'), help='Raw rNMP count for each cds')
    parser.add_argument('info', type=argparse.FileType('r'), help='List of library information')
    parser.add_argument('gtf', type=argparse.FileType('r'), help='GTF annotation')
    parser.add_argument('mito_count', type=argparse.FileType('r'), help='CSV file of mitochondrial rNMP count for each library')
    parser.add_argument('bg_count', type=argparse.FileType('r'), help='TSV file of mitochondrial rNMP count for each gene')
    parser.add_argument('-c', type=int, default=1, help='Column number of library name in libinfo, default=1')
    parser.add_argument('-s', choices={'template', 'nontemplate'}, default='nontemplate', help='Input strand: template/(nontemplate)')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output to file')
    args = parser.parse_args()
    args.c -= 1

    # read
    cds_info = read_gtf(args.gtf)
    libinfo = read_libinfo(args.info, args.c)
    df = read_data(args.raw)
    mito_count = read_mito_count(args.mito_count)
    bg_count = read_mono(args.bg_count, args.s)

    # add information
    df['Total_chrM_rNMPs'] = df.Library.map(mito_count)
    df['Total'] = df['A'] + df['C'] + df['G'] + df['T']
    df = df.merge(libinfo, left_on='Library', right_index=True)
    df = df.merge(cds_info, on='Gene_id')
    df['loc'] = df.apply(generate_loc, axis=1)
    df = df.merge(bg_count, left_on='loc', right_index=True)
    df['Total_bg'] = df['End'] - df['Start']
    df['Length'] = df['End'] - df['Start']

    # Calculate RPB and PPB
    for fea in ['A', 'C', 'G', 'T', 'Total']:
        df[f'RPB_{fea}'] = df[fea]/df[f'{fea}_bg']
        df[f'PPB_{fea}'] = df[f'RPB_{fea}']/df['Total_chrM_rNMPs']
    df = df.rename(columns={
        'RPB_Total':'RPB', 'PPB_Total':'PPB'
        })

    # Get library order
    df = df.sort_values(by=['Order', 'Length'], ascending=[True, False])
    # remove useless columns and rows
    df = df.drop(columns=['loc', 'Total_bg', 'Order'])
    df = df.dropna(axis=0)


    # output
    df.to_csv(args.o, sep='\t', index=False)


    print('Done!')

if __name__ == '__main__':
    main()
