#!/usr/bin/env python3

import seaborn as sns
import argparse
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import PercentFormatter
from statannotations.Annotator import Annotator


def get_category(data):
    return '-'.join(list(data.split('-')[:-1]))

# generate dataframe from tsv file
def main():
    parser = argparse.ArgumentParser(description='Generate barplot with strand split file')
    parser.add_argument('fwd', type=argparse.FileType('r'), help='rNMP count for forward strand')
    parser.add_argument('rev', type=argparse.FileType('r'), help='rNMP count for reverse strand')
    parser.add_argument('-l', type=float, default=7, help='Width of the figure, (7)')
    parser.add_argument('--no_annot', action='store_true', help='Do not draw statistical annotations')
    parser.add_argument('-o', help='Output plot name')
    args = parser.parse_args()

    if not args.o:
        args.o = args.tsv.name.split('.')[0] + '_strand_split.png'

    # define color palette
    pal = ['#FB5156', '#5555FF']

    # generate df
    fwd = pd.read_csv(args.fwd, sep='\t')
    fwd['Strand'] = 'Forward'
    rev = pd.read_csv(args.rev, sep='\t')
    rev['Strand'] = 'Reverse'
    df = pd.concat([fwd, rev])
    dtypes = { x:'float' for x in ['A','C','G','T']}
    df = df.astype(dtypes)
    df['Total'] = df['A'] + df['C'] + df['G'] + df['T']
    both_dict = df.groupby('Sample').sum()['Total'].to_dict()
    df['Both'] = df.Sample.map(both_dict)
    df['Ratio'] = df['Total']/df['Both']
    df['Celltype'] = df.Sample.apply(get_category)

    # draw
    sns.set(font_scale=1.5, style='ticks')
    fig, ax = plt.subplots(figsize=(args.l,6), dpi=300)
    plt.subplots_adjust(left=0.16, right=1, top=0.98, bottom=0.4)
    sns.barplot(x='Celltype', y='Ratio', hue='Strand', data=df, palette=pal,\
        capsize=0.3, errwidth=1.2, ax=ax)
    sns.swarmplot(x='Celltype', y='Ratio', hue='Strand', data=df, color='k', \
        size=3.5, dodge=True,ax=ax)

    # figure settings
    ax.get_legend().remove()
    sns.despine()
    plt.ylabel('')
    plt.xlabel('')
    plt.setp(ax.patches, linewidth=1)
    plt.ylim((0,1))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))

    if not args.no_annot:
        # draw annotations
        pairs = []
        counts = df.groupby('Celltype').Ratio.count()
        for k, v in counts.to_dict().items():
            if v >= 6:
                pairs.append(((k, 'Forward'), (k, 'Reverse')))
        annot = Annotator(ax, pairs, data=df, x='Celltype', y='Ratio', hue='Strand')
        annot.configure(test='t-test_ind', verbose=2)
        annot.apply_test()
        annot.annotate()

    # tick labels
    xticklabels = []
    for l in ax.get_xticklabels():
        geno = l.get_text()
        count = len(df[df.Celltype == geno])/2
        xticklabels.append(geno + f'\n(N = {count:.0f})')
    rotation = 70 if args.l < 10 else 0
    ax.set_xticklabels(xticklabels,fontsize=16, rotation=rotation)
    plt.savefig(args.o)

    print('Done!')

if __name__ == '__main__':
    main()
