#!/usr/bin/env python3

from collections import defaultdict
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# replication origin coordinate
cr1_s = 16024
cr1_e = 16569
cr2_s = 0
cr2_e = 576
oh_s = 110
oh_e = 441


# read total reads
def read_total(fr):
    total = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        if len(ws) == 2:
            total[ws[0]] = int(ws[1])
    return total


# read celltypes
def read_celltypes(fr):
    data = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        if len(ws) == 2:
            data[ws[1]] = ws[0]
    return data


# load rNMPs in bed file and calculate moving avg
def load_bed(same, oppo, total, flank):
    global cr1_s, cr1_e, cr2_s, cr2_e, oh_s, oh_e
    df = []
    beds = [[fr, 'Light'] for fr in same]
    beds += [[fr, 'Heavy'] for fr in oppo]
    for bed, st in beds:
        # Filter out libraries
        name = bed.name.split('/')[-1].split('.')[0]
        if name not in total:
            continue
        # initialization    
        data = np.zeros(cr1_e-cr1_s + cr2_e-cr2_s)
        for l in bed:
            ws = l.rstrip().split('\t')
            if len(ws) < 6:
                continue
            s = int(ws[1])
            e = int(ws[2])
            if ws[0] == 'chrM':
                if st == 'Light' and ws[5] == '+':
                    if cr1_s <= s < cr1_e:
                        data[s-cr1_s] += 1
                    elif cr2_s <= s <cr2_e:
                        data[s+cr1_e-cr1_s] += 1
                elif st == 'Heavy' and ws[5] == '-':
                    if cr1_s <= s < cr1_e:
                        data[cr1_s - e] += 1
                    elif cr2_s <= s < cr2_e:
                        data[cr2_e - e] += 1
        data /= total[name]
        data /= data.shape[0]
        data = np.convolve(data, np.ones(flank*2+1)/(flank*2+1), mode='same')
        for i in range(data.shape[0]):
            if st == 'Heavy':
                d = [name, st, len(data)-1-i, data[i]]
            else:
                d = [name, st, i, data[i]]
            if st == 'Light':
                if i < cr1_e - cr1_s + oh_s:
                    d.append('Upstream')
                elif i < cr1_e - cr1_s + oh_e:
                    d.append('OriH')
                else:
                    d.append('Downstream')
            else:
                if i < cr2_e - oh_e:
                    d.append('Upstream')
                elif cr2_e - oh_e <= i < cr2_e - oh_s:
                    d.append('OriH')
                else:
                    d.append('Downstream')
            df.append(d)
    df = pd.DataFrame(df, columns=['Library', 'Strand', 'Position', 'Moving_avg', 'Region'])
    return df


# draw histogram
def draw_linecharts(data, out):
    global oh_s, oh_r, cr1_s, cr1_e, cr2_s, cr2_e
    colors = ['blue', 'red', 'blue']
    sns.set(style='ticks', font_scale=1.5)
    celltypes = data['Celltype'].dropna().unique()
    for celltype in celltypes:
        df = data[data.Celltype == celltype].sort_values(by='Position')
        df.to_csv(f'{out}_{celltype}.csv')
        fig, axs = plt.subplots(2,1, figsize=(10,8))
        plt.subplots_adjust(left=0.08, top=0.92, right=0.99, bottom=0.05, hspace=0.14)
        # draw heavy and light separately
        for i, st in enumerate(['Light', 'Heavy']):
            df_curr = df[df.Strand == st]
            sns.lineplot(x='Position', y='Moving_avg', hue='Region', hue_order=['Upstream', 'OriH', 'Downstream'],\
                data=df_curr, ci="sd", palette=colors, ax=axs[i], legend=False)
            axs[i].set_xlim((0, cr2_e - cr2_s + cr1_e - cr1_s + 1))
            axs[i].set_ylim((-0.95e-7, 6e-7))
            axs[i].set_ylabel(f'{st} PPB')
            if st == 'Light':
                axs[i].set_xlabel('')
                sns.despine(ax = axs[i])
            else:
                axs[i].invert_yaxis()
                axs[i].set_xlabel('Position')
                axs[i].xaxis.set_ticklabels([])
                axs[i].yaxis.get_offset_text().set_visible(False)
                sns.despine(ax = axs[i], top=False, bottom=True)
        # tick labels
        n = len(df.Library.unique())
        plt.suptitle(f'{celltype}   N={n}')
        fig.savefig(f'{out}_{celltype}.png')
        plt.close()



def main():
    # argparse
    parser = argparse.ArgumentParser(description='Calculate the distribution of rNMPs in human mt replication origin region (MT-CR)')
    parser.add_argument('total',type=argparse.FileType('r'), help='Tsv file of total rNMP counts for each libraries')
    parser.add_argument('celltypes', type=argparse.FileType('r'), help='List of cell types')
    parser.add_argument('--same', type=argparse.FileType('r'), nargs='+', help='Bed files on the same strand of MT-CR')
    parser.add_argument('--oppo', type=argparse.FileType('r'), nargs='+', help='Bed files on the opposite strand of MT-CR')
    parser.add_argument('-f', type=int, default=25, help='Flank length at both direction, default = 25 nt')
    parser.add_argument('-o', default='mt_cr', help='Output basename')
    args = parser.parse_args()

    assert len(args.same) == len(args.oppo) and len(args.same) >=1, 'There should be at least one pair of same/oppo bed files'

    # total reads
    total = read_total(args.total)

    # annotations for each cell types
    celltypes = read_celltypes(args.celltypes)

    # load data
    ppb = load_bed(args.same, args.oppo, total, args.f)

    # add annotations
    ppb['Celltype'] = ppb['Library'].map(celltypes)

    # output plots
    draw_linecharts(ppb, args.o)

    print('Done!')

if __name__ == '__main__':
    main()
