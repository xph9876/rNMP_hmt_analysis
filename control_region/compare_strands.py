#!/usr/bin/env python3

from collections import defaultdict
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator

# read celltypes
def read_celltypes(fr):
    data = {}
    celltypes = []
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        if ws[0] == 'Control':
            continue
        if len(ws) == 2:
            data[ws[1]] = ws[0]
        if not celltypes or celltypes[-1] != ws[0]:
            celltypes.append(ws[0])
    return data, celltypes

# read data
def load_data(same, oppo, celltypes, celltype_list):
    data = defaultdict(lambda: {'+':0, '-':0})
    # skip header
    same.readline()
    oppo.readline()
    for l in same:
        ws = l.rstrip('\n').split('\t')
        if len(ws) < 5:
            continue
        data[ws[0]]['+'] = sum([float(x) for x in ws[1:]])
    for l in oppo:
        ws = l.rstrip('\n').split('\t')
        if len(ws) < 5:
            continue
        data[ws[0]]['-'] = sum([float(x) for x in ws[1:]])
    # calculate percentage
    df = []
    for k, v in data.items():
        if k not in celltypes:
            continue
        sm = v['+'] + v['-']
        df.append([k, celltypes[k], 'Light', v['+']/sm])
        df.append([k, celltypes[k], 'Heavy', v['-']/sm])
    df = pd.DataFrame(df, columns=['Library', 'Celltype', 'Strand', 'Ratio'])
    df.Strand = pd.Categorical(df.Strand, ['Light', 'Heavy'])
    df.Celltype = pd.Categorical(df.Celltype, celltype_list)
    df = df.sort_values(by='Celltype')
    return df

# draw barplot
def draw_barplot(data, out, annotation):
    colors = {'Light':'#FB5156', 'Heavy':'#5555FF'}
    sns.set(style='ticks', font_scale=1.5)
    fig, ax = plt.subplots(figsize=(7,6), dpi=300)
    plt.subplots_adjust(left=0.16, top=0.98, right=1, bottom=0.4)
    sns.barplot(x='Celltype', y='Ratio', hue='Strand', data=data, errorbar='sd', errwidth=1.2, \
        capsize=0.3, palette=colors, ax=ax)
    sns.swarmplot(x='Celltype', y='Ratio', hue='Strand', data=data, dodge=True, \
        color='k', size=3.5, ax=ax)
    plt.legend().set_visible(False)
    sns.despine()
    plt.ylim([0,1.05])

    # draw annotations
    if annotation:
        pairs = []
        counts = data.groupby('Celltype').Ratio.count()
        for k, v in counts.to_dict().items():
            if v >= 6:
                pairs.append(((k, 'Light'), (k, 'Heavy')))
        annot = Annotator(ax, pairs, data=data, x='Celltype', y='Ratio', hue='Strand')
        annot.configure(test='t-test_ind', verbose=2)
        annot.apply_test()
        annot.annotate()

    # tick labels
    xticklabels = []
    for l in ax.get_xticklabels():
        geno = l.get_text()
        count = len(data[data.Celltype == geno])/2
        xticklabels.append(geno.replace(' ', '\n') + f'\n(N = {count:.0f})')
    ax.set_xticklabels(xticklabels, fontsize=12, rotation=70)

    # labels
    plt.xlabel('')
    plt.ylabel('')

    fig.savefig(f'{out}_barplot.png')


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Compare rNMP incorporation counts in the same and opposite strand')
    parser.add_argument('same', type=argparse.FileType('r'), help='Raw composition file for the same strand')
    parser.add_argument('oppo', type=argparse.FileType('r'), help='Raw composition file for the opposite strand')
    parser.add_argument('celltypes', type=argparse.FileType('r'), help='List of library informations')
    parser.add_argument('-o', default='out', help='Output figure basename')
    parser.add_argument('--no_annot', action='store_true', help='Do not draw statistical annotations')
    args = parser.parse_args()

    # annotations for each cell types
    celltypes, celltype_list = read_celltypes(args.celltypes)

    # load data
    data = load_data(args.same, args.oppo, celltypes, celltype_list)

    # output plots
    draw_barplot(data, args.o, not args.no_annot)

    print('Done!')

if __name__ == '__main__':
    main()
