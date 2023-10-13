#!/usr/bin/env python3

import pandas as pd
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import sys


# draw figures
def draw(df, fea, out):
    # set color
    color = defaultdict(lambda:'#AAAAAA')
    for geno in ['DLTB-8', 'DLTB-P', 'TLTB-8', 'TLTB-P']:
        color[geno] = '#AAAAFF'
    sns.set(style='ticks', font_scale=2.5)
    fig, ax = plt.subplots(figsize=(6, 8), dpi=300)
    plt.subplots_adjust(bottom=0.45)
    sns.barplot(x='Sample', y=fea, data=df, palette=color, ax=ax)
    sns.despine()
    xmin, xmax = ax.get_xlim()
    plt.plot([xmin,xmax], [1,1], 'k--')
    labels = ax.get_xticklabels()
    ax.set_xticklabels(labels, rotation=90)
    plt.ylabel('')
    plt.xlabel('')
    fig.savefig(out)


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw sequencing depth bar plot for 7S DNA')
    parser.add_argument('tsv', type=argparse.FileType('r'), help='7S sequencing depth tsv file')
    parser.add_argument('-o', default='7s', help='Output base name, (7s)')
    args = parser.parse_args()

    # calculate coverage
    df = pd.read_csv(args.tsv, sep='\t')

    # Draw figures
    draw(df, '7S_EF', args.o + '_ef.png')
    draw(df, '7S_rev_EF', args.o + '_rev_ef.png')


    print('Done!')



if __name__ == '__main__':
    main()
