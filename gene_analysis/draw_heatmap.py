#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# generate an order dict of genes
def read_order(f):
    l = [i.rstrip('\n') for i in f if len(i) >= 2]
    return {l[i]:i for i in range(len(l))}

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw graphs for CDS of human mt ')
    parser.add_argument('data',type=argparse.FileType('r'), help='CDS CSV data file')
    parser.add_argument('-l', type=argparse.FileType('r'), help='Add an list of gene order')
    parser.add_argument('--no_cbar', action='store_false', help='Do not draw color bar')
    parser.add_argument('-o', help='output basename')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.data,sep='\t')
    libs = df.Library.unique()
    # follow the order of gene list
    if args.l:
        gene_order = read_order(args.l)
        df['Order'] = df['Gene_name'].map(gene_order)
        df.sort_values(by='Order', ascending=True, inplace=True)
    else:
        df.sort_values(by='Length', ascending=False, inplace=True)

    # generate labels
    gene_names = df.Gene_name.unique()
    lengths = df.Length.unique()
    genotypes = df.Genotype.unique()
    len_dict = pd.Series(df.Length.values, index=df.Gene_name).to_dict()
    geno_dict = pd.Series(df.Genotype.values, index=df.Library).to_dict()
    ylabel_texts = [f'{len_dict[x]}nt {x}' for x in gene_names]
    xlabel_texts = [f'{geno_dict[x]}-{x}' for x in libs]

    for fea in ['EF'] + [f'EF_{r}' for r in 'ACGT']:
        data = df.pivot(index="Gene_name", columns="Library", values=fea).fillna(0)
        # sort
        data = data.loc[gene_names, libs]
        # set color settings for graph
        sns.set(font_scale=1.5)
        if len(gene_names) < 15:
            fig, ax = plt.subplots(figsize=(15,10))
            plt.subplots_adjust(left=0.18, top=0.99, bottom=0.5, right=0.95)
        else:
            fig, ax = plt.subplots(figsize=(15,12))
            plt.subplots_adjust(left=0.18, top=0.99, bottom=0.3, right=0.95)
        if fea == 'EF':
            sns.heatmap(data, vmin=0, vmax=3.4, ax=ax, cmap='rocket')
        else:
            sns.heatmap(data, vmin=0, vmax=5.4, ax=ax, cmap='rocket')
        ax.set_yticklabels(ylabel_texts, rotation='horizontal')
        ax.set_xticklabels(xlabel_texts, rotation='vertical')
        
        plt.xlabel('')
        plt.ylabel('')
        # show or save
        if not args.o:
            plt.show()
        else:
            fig.savefig(args.o + f'_{fea}_heatmap.png')
            print(f'Heatmap is saved to {args.o}_{fea}_heatmap.png')

if __name__ == '__main__':
    main()
