#!/usr/bin/env python3

import seaborn as sns
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# Load library information
def load_info(fr):
    info = {}
    for i, l in enumerate(fr):
        ws = l.rstrip().split('\t')
        info[ws[1]] = (i, ws[0])
    return info 

# load background frequence
def get_bg(fr):
    data = {}
    header = fr.readline().rstrip('\n').split('\t')
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        for i in range(1, len(ws)):
            data[header[i]] = float(ws[i])
    return data



# get category from info
def get_category(x, info):
    if x in info:    
        return info[x][1]
    return ''


# get order from info
def get_order(x, info):
    if x in info:    
        return info[x][0]
    return ''


# Get difference
def get_diff(x, fea):
    return x[f'{fea}_norm_Light'] - x[f'{fea}_norm_Heavy']

# Get preferred:
def get_preferred(x):
    if x.Total_Light > x.Total_Heavy:
        return 'Light'
    return 'Heavy'

# Calculate z score
def calc_z(x, fea):
    feas = [x[f'{feas}_diff_norm'] for feas in list('ACGT')]
    return (x[f'{fea}_diff_norm'] - np.mean(feas))/np.std(feas)


# Load count
def load_data(fwd, rev, info, l_bg, h_bg):
    # load fwd
    fwd = pd.read_csv(fwd, sep='\t')
    fwd['Strand'] = 'Light'
    for k, v in l_bg.items():
        fwd[f'{k}_bg'] = v
    # load rev
    rev = pd.read_csv(rev, sep='\t')
    rev['Strand'] = 'Heavy'
    for k, v in h_bg.items():
        rev[f'{k}_bg'] = v
    # Add total
    df = pd.concat([fwd, rev])
    df['Total'] = df['A'] + df['C'] + df['G'] + df['T']
    # Both strand data
    both = df.groupby('Sample').sum(numeric_only=True).reset_index()
    both['Strand'] = 'Both'
    df = pd.concat([df, both])
    # Normalize
    for fea in ['A', 'C', 'G', 'T']:
        df[f'{fea}_norm'] = df[fea]/df[f'{fea}_bg']
    # Merge
    light = df[df.Strand == 'Light'].drop(columns='Strand')
    heavy = df[df.Strand == 'Heavy'].drop(columns='Strand')
    both = df[df.Strand == 'Both'].drop(columns='Strand')
    df = light.merge(heavy, on='Sample', suffixes=('_Light', '_Heavy'))
    df = df.merge(both, on='Sample')
    df['Preferred'] = df.apply(get_preferred, axis=1)
    # Difference of light/heavy after normalization
    for fea in ['A', 'C', 'G', 'T']:
        df[f'{fea}_diff'] = df.apply(get_diff, axis=1, fea=fea)
    # add info
    df['Celltype'] = df.Sample.apply(get_category, info=info)
    df['Order'] = df.Sample.apply(get_order, info=info)
    df = df[df.Celltype!=''].sort_values(by=['Order'], ascending=[True])
    # Calculate contribution
    data = []
    for i, x in df.iterrows():
        data.append([x.Sample])
        diff_percent = abs(x.Total_Light-x.Total_Heavy)/(x.Total_Light+x.Total_Heavy)  
        if x.Preferred == 'Light':
            total = sum(max(x[f'{fea}_diff'],0) for fea in 'ACGT')
            data[-1] += [
                max(0, x[f'{fea}_diff']/total) * diff_percent * 100 for fea in 'ACGT'
            ]
        else:
            total = sum(min(x[f'{fea}_diff'],0) for fea in 'ACGT')
            data[-1] += [
                -max(0, x[f'{fea}_diff']/total)* diff_percent * 100 for fea in 'ACGT'
            ]
    df = pd.DataFrame(data, columns=['Sample', 'A', 'C', 'G', 'U']).set_index('Sample')
    return df


# Draw heatmap
def draw(df, out):
    sns.set(font_scale=1.8, style='ticks')
    fig, ax = plt.subplots(figsize=(4,10), dpi=300)
    plt.subplots_adjust(left=0.25, right=0.8, top=0.98, bottom=0.05)
    sns.heatmap(df, 
        cmap=sns.diverging_palette(145, 300, s=60, as_cmap=True), 
        ax=ax, 
        vmin=-39, 
        vmax=39
    )
    plt.ylabel('')
    plt.xlabel('')
    fig.savefig(f'{out}_contribution.png')


# generate dataframe from tsv file
def main():
    parser = argparse.ArgumentParser(description='Generate heatmap with rNMP contribution in either strand')
    parser.add_argument('fwd_count', type=argparse.FileType('r'), help='rNMP count for light strand')
    parser.add_argument('rev_count', type=argparse.FileType('r'), help='rNMP count for heavy strand')
    parser.add_argument('fwd_bg', type=argparse.FileType('r'), help='background count for light strand')
    parser.add_argument('rev_bg', type=argparse.FileType('r'), help='background count for heavy strand')
    parser.add_argument('info', type=argparse.FileType('r'), help='Library information')
    parser.add_argument('-o', default='strand_split_contribution', help='Output plot name')
    args = parser.parse_args()

    # get library information
    info = load_info(args.info)

    # get background frequence
    l_bg = get_bg(args.fwd_bg)
    h_bg = get_bg(args.rev_bg)

    # Format data
    df = load_data(args.fwd_count, args.rev_count, info, l_bg, h_bg)

    # Draw figures
    draw(df, args.o)

    print('Done!')

if __name__ == '__main__':
    main()
