#!/usr/bin/env python3

from collections import defaultdict
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

# read celltypes
def read_celltypes(fr):
    data = {}
    celltypes = []
    for l in fr:
        ws = l.rstrip('\n').split('\t')
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
def draw_barplot(data, out):
    colors = {'Light':'#FB5156', 'Heavy':'#5555FF'}
    sns.set(style='ticks', font_scale=1.5)
    fig, ax = plt.subplots(figsize=(15,6))
    plt.subplots_adjust(left=0.06, top=0.92, right=0.99, bottom=0.2)
    sns.barplot(x='Celltype', y='Ratio', hue='Strand', data=data, ci="sd", errwidth=2, \
        capsize=0.2, palette=colors, ax=ax)
    sns.swarmplot(x='Celltype', y='Ratio', hue='Strand', data=data, dodge=True, \
        color='k', ax=ax)
    plt.legend().set_visible(False)
    sns.despine()
    plt.ylim([0,1.05])
    # tick labels
    labels = [x.get_text().replace(' ','\n') for x in ax.get_xticklabels()]
    ax.set_xticklabels(labels)
    ax.set_xlabel('Cell type')
    ax.set_ylabel('Ratio')

    fig.savefig(f'{out}_barplot.png')

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Compare rNMP incorporation counts in the same and opposite strand')
    parser.add_argument('same', type=argparse.FileType('r'), help='Raw composition file for the same strand')
    parser.add_argument('oppo', type=argparse.FileType('r'), help='Raw composition file for the opposite strand')
    parser.add_argument('celltypes', type=argparse.FileType('r'), help='List of library informations')
    parser.add_argument('-o', default='out', help='Output figure basename')
    args = parser.parse_args()

    # annotations for each cell types
    celltypes, celltype_list = read_celltypes(args.celltypes)

    # load data
    data = load_data(args.same, args.oppo, celltypes, celltype_list)

    # output plots
    draw_barplot(data, args.o)

    print('Done!')

if __name__ == '__main__':
    main()
