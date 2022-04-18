#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster import hierarchy as hc
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
    genotypes = list(df.Genotype.unique())
    len_dict = pd.Series(df.Length.values, index=df.Gene_name).to_dict()
    geno_dict = pd.Series(df.Genotype.values, index=df.Library).to_dict()
    ylabel_texts = [f'{len_dict[x]}nt {x}' for x in gene_names]
    xlabel_texts = [f'{geno_dict[x]}-{x}' for x in libs]
    colors = sns.hls_palette(len(genotypes), s=1, l=0.5)
    color_dict = {x:colors[genotypes.index(x.split('-')[0])] for x in xlabel_texts}

    for fea in ["Count_rank","RPB_rank",'PPB']:
        data = df.pivot("Library", "Gene_name", fea)
        # sort
        data = data.loc[libs, gene_names]
        # cluster
        z = hc.linkage(data, method='average', optimal_ordering=False)
        # set color settings for graph
        sns.set(font_scale=2, style='ticks')
        fig, ax = plt.subplots(figsize=(7,6))
        plt.subplots_adjust(left=0.02, top=0.99, bottom=0.4, right=0.98)
        hc.dendrogram(z, color_threshold=0, labels=xlabel_texts)
        sns.despine()
        for l in ax.xaxis.get_ticklabels():
            l.set_rotation('vertical')
            l.set_fontsize('15')
            text = l.get_text()
            l.set_color(color_dict[text])

        for l in ax.yaxis.get_ticklabels():
            l.set_visible(False)
        # plt.xlabel('')
        # plt.ylabel('')
    #     ax.set_title('Heatmap for dinucleotides')
        # ax.set_xticklabels(sample, rotation='vertical')
        # ax.set_yticklabels(label_texts,rotation='horizontal')

        # # change top of colorbar to '0.5-1'
        # cax = plt.gcf().axes[-1]
        # color_labels = cax.get_ymajorticklabels()
        # color_labels_texts = [i.get_text() for i in color_labels]
        # color_labels_texts[-1] += ' - 1'
        # cax.set_yticklabels(color_labels_texts)

        # show or save
        if not args.o:
            plt.show()
        else:
            fig.savefig(args.o + f'_{fea}_cluster.png')
            print(f'Cluster results is saved to {args.o}_{fea}_cluster.png')

if __name__ == '__main__':
    main()
