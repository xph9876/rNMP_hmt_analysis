#!/usr/bin/env python3

from collections import defaultdict
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import matplotlib.ticker as ticker
from curlyBrace import curlyBrace

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
def draw_linecharts(df, out, pal, draw_legend=False):
    colors = {x:sns.color_palette(pal)[i] for i, x in enumerate(df.Celltype.unique())}
    global oh_s, oh_r, cr1_s, cr1_e, cr2_s, cr2_e
    sns.set(style='ticks', font_scale=2)
    df = df.sort_values(by='Position')
    df.to_csv(f'{out}.csv')
    fig, axs = plt.subplots(2,1, figsize=(12,8),dpi=300)
    if draw_legend:
        plt.subplots_adjust(left=0.05, top=0.92, right=0.7, bottom=0.05, hspace=0.18)
    else:
        plt.subplots_adjust(left=0.05, top=0.92, right=0.95, bottom=0.05, hspace=0.18)
    # draw heavy and light separately
    for i, st in enumerate(['Light', 'Heavy']):
        for geno in df.Celltype.unique():
            df_curr = df[(df.Strand == st) & (df.Celltype == geno)]
            axs[i].plot(df_curr.Position, df_curr.Moving_avg, color=colors[geno], label=geno)
        axs[i].set_xlim((0, cr2_e - cr2_s + cr1_e - cr1_s + 1))
        axs[i].set_ylim((0, 4.8e-7))
        locs = [16100-cr1_s, 16300-cr1_s, cr1_e-cr1_s, cr1_e+200-cr1_s, cr1_e+400-cr1_s]
        axs[i].xaxis.set_major_locator(ticker.FixedLocator(locs))
        if st == 'Light':
            axs[i].set_xlabel('')
            sns.despine(ax = axs[i])
            axs[i].xaxis.set_ticklabels(['16100', '16300', '16569/1', '200', '400'])
            # add LSP
            axs[i].plot((cr1_e + 407 - cr1_s, cr1_e + 407 - cr1_s), (0, 0.05e-7), 'k-')
        else:
            axs[i].invert_yaxis()
            axs[i].set_xlabel('Position')
            axs[i].xaxis.set_ticklabels([])
            axs[i].yaxis.get_offset_text().set_visible(False)
            sns.despine(ax = axs[i], top=False, bottom=True)
            # add HSP
            axs[i].plot((cr1_e + 560 - cr1_s, cr1_e + 560 - cr1_s), (0, 0.05e-7), 'k-')
    if draw_legend:
        plt.legend(loc=(1.02, 0.2))
    else:
        plt.legend('',frameon=False)
    # highlight OriH region
    curlyBrace(
        fig, axs[0], 
        p1=(cr1_e+oh_s-cr1_s, 1.6e-7), 
        p2=(cr1_e+oh_e-cr1_s, 1.6e-7),
        str_text='OriH',
        color='black'
        )
    curlyBrace(
        fig, axs[0], 
        p1=(16106-cr1_s, 3.9e-7), 
        p2=(cr1_e+191-cr1_s, 3.9e-7),
        k_r=0.05,
        str_text='D-loop',
        color='black'
        )
    # tick labels
    fig.savefig(f'{out}.png')
    plt.close()



def main():
    # argparse
    parser = argparse.ArgumentParser(description='Calculate the distribution of rNMPs in human mt replication origin region (MT-CR), combined by celltypes')
    parser.add_argument('total',type=argparse.FileType('r'), help='Tsv file of total rNMP counts for each libraries')
    parser.add_argument('celltypes', type=argparse.FileType('r'), help='List of cell types')
    parser.add_argument('-f', type=int, default=25, help='Flank length at both direction, default = 25 nt')
    parser.add_argument('-o', default='mt_cr', help='Output basename')
    parser.add_argument('--same', type=argparse.FileType('r'), nargs='+', help='Bed files on the same strand of MT-CR')
    parser.add_argument('--oppo', type=argparse.FileType('r'), nargs='+', help='Bed files on the opposite strand of MT-CR')
    parser.add_argument('--selected', default=['CD4T', 'hESC-H9','DLTB', 'TLTB', 'WB-GTP control', 'WB-GTP PTSD', 'HCT116', 'HEK293T'], nargs='+', help='Selected genotypes, (All WT > 3)')
    parser.add_argument('--palette', default='Dark2', help='Color paletter used to generate plots, (Set1)')
    parser.add_argument('--draw_legend', action='store_true', help='Draw figure legend')
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
    ppb = ppb[ppb.Celltype.isin(args.selected)].dropna().copy()
    ppb = ppb.groupby(['Celltype', 'Strand', 'Position', 'Region']).mean(numeric_only=True).reset_index()
    ppb.to_csv('test.csv')

    # output plots
    draw_linecharts(ppb, args.o, args.palette, args.draw_legend)

    print('Done!')

if __name__ == '__main__':
    main()
