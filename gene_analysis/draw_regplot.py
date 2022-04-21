#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw line chart to show the relation of gene size and PPB')
    parser.add_argument('data',type=argparse.FileType('r'), help='CDS CSV data file')
    parser.add_argument('-o', default='human_genes', help='output basename')
    args = parser.parse_args()

    # read data
    df = pd.read_csv(args.data,sep='\t')
    libs = df.Library.unique()
    genotypes = df.Genotype.unique()
    genes = df.Gene_name.unique()

    # normalization on baseline
    base = df.groupby(['Library']).mean().reset_index()[['Library', 'PPB']]
    df = df.merge(base, left_on='Library', right_on='Library', how='left', suffixes=('','_baseline'))
    df['PPB_norm'] = df['PPB']/df['PPB_baseline']

    # draw
    sns.set(style='ticks', font_scale=2.5)
    for genotype in genotypes:
        data = df[df.Genotype==genotype]
        for fea in ['PPB_norm']:
            fig, ax = plt.subplots(figsize=(6,6))
            plt.subplots_adjust(left=0.1, top=0.8, bottom=0.1)
            g = sns.lmplot(data=data, x='Length', y=fea, legend=False)
            sns.despine()
            r, p = spearmanr(data['Length'], data[fea])
            plt.xlim((0, 2000))
            plt.ylim((0, 3))
            plt.xlabel('')
            plt.ylabel('')
            plt.title(f'r={r:.2f} p={p:.2g}')
            plt.savefig(f'{args.o}_{genotype}_{fea}_reg.png')
            plt.close('all')

    print('Done!')


if __name__ == '__main__':
    main()
