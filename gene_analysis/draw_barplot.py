#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from matplotlib.ticker import FuncFormatter

# generate an order dict of genes
def read_genes(f):
    types = {}
    lengths = {}
    for l in f:
        ws = l.split('\t')
        if len(ws) < 3:
            continue
        types[ws[0]] = ws[1]
        lengths[ws[0]] = int(ws[2])
    return types, lengths


# generate tick labels
@FuncFormatter
def to_scientific(x, pos):
    if not x:
        return '0'
    s = f'{x:e}'
    ws = s.split('e')
    a = str(float(ws[0]))
    b = str(int(ws[1]))
    return f'{a}E{b}'

    


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw bargraph to compare PPB of each human mt gene')
    parser.add_argument('data',type=argparse.FileType('r'), help='CDS CSV data file')
    parser.add_argument('genes', type=argparse.FileType('r'), help='List of genes')
    parser.add_argument('-o', default='human_genes', help='output basename')
    parser.add_argument('--baseline', default='ND5', help='Baseline gene name, (ND5)')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.data,sep='\t')
    libs = df.Library.unique()
    gene_types, gene_lengths = read_genes(args.genes)

    # set color
    colors = {'tDNA':'#1b9e77','CDS':'#d95f02','Coding':'#d95f02','rDNA':'#7570b3', 'ncDNA':'#e41a1c'}
    gene_colors = {k:colors[v] for k, v in gene_types.items()}

    # append gene type and gene length to csv
    df['Gene_type'] = df['Gene_name'].map(gene_types)
    df['Gene_type'] = pd.Categorical(df['Gene_type'], ['CDS', 'Coding', 'rDNA','ncDNA','tDNA'])
    df = df.sort_values(by=['Gene_type', 'Length'])
    df = df[~df.Gene_type.isin({'ncDNA', 'tDNA'})].copy()

    # normalize on baseline genes
    assert args.baseline in df.Gene_name.unique(), f'Baseline gene {args.baseline} not in the input files'
    nd5 = df[df.Gene_name == args.baseline][['Library', 'PPB']]
    df = df.merge(nd5, how='left', left_on='Library', right_on='Library', suffixes=('', '_base'))
    df['Enrichment_factor'] = df.PPB/df.PPB_base
    df.to_csv('test.tsv', sep='\t')

    # barplot
    sns.set(font_scale=1.5, style='ticks')
    fig, ax = plt.subplots(figsize=(6, 6))
    plt.subplots_adjust(left=0.1, top=0.98, bottom=0.4, right=0.98)
    sns.barplot(x='Gene_name', y='Enrichment_factor', ci='sd', data=df, \
        errwidth=1.3, capsize=0.3, ax=ax, palette=gene_colors)
    sns.stripplot(x='Gene_name', y='Enrichment_factor', dodge=True, color='black', size=2, data=df, ax=ax)
    sns.despine()
    ax.set_xticklabels([f'{gene_lengths[x.get_text()]}bp-{x.get_text()}' for x in ax.get_xticklabels()], rotation=90)
    # plt.ylim([0, 2])
    # ax.yaxis.set_major_formatter(to_scientific)
    plt.ylabel('')
    plt.xlabel('')
    fig.savefig(args.o + '_bar.png')

    print(f'Bar plot saved as {args.o}_bar.png')

if __name__ == '__main__':
    main()
