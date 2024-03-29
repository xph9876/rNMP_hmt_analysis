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

    # draw
    sns.set(style='ticks', font_scale=2.5)
    for genotype in genotypes:
        data = df[df.Genotype==genotype]
        for fea in ['EF'] + [f'EF_{r}' for r in 'ACGT']:
            fig, ax = plt.subplots(figsize=(6,6))
            plt.subplots_adjust(left=0.1, top=0.8, bottom=0.1)
            g = sns.lmplot(data=data, x='Length', y=fea, legend=False)
            sns.despine()
            r, p = spearmanr(data['Length'], data[fea])
            plt.xlim((0, 2000))
            if fea == 'EF':
                plt.ylim((0, 2.4))
            plt.xlabel('')
            plt.ylabel('')
            plt.title(f'r={r:.2f} p={p:.2g}')
            plt.savefig(f'{args.o}_{genotype}_{fea}_reg.png')
            plt.close('all')

    print('Done!')


if __name__ == '__main__':
    main()
