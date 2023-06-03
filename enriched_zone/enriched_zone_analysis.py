#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys


# get mitochondrial DNA size
def get_mt_size(fai, name):
    # get mt length
    for l in fai:
        ws = l.split('\t')
        if ws[0] == name:
            size = int(ws[1])
            return size
    return None


# load rNMP coordinates
def load_bed(bed, size, info, bin, name, selected):
    data = {}
    counts = {}
    for fr in bed:
        fs = fr.name.split('/')[-1].split('.')[0]
        if fs not in info:
            print(f'Library {fs} doesn\'t have related information, skipped')
            continue
        counts[fs] = {'+':0, '-':0}
        geno = info[fs]
        if selected and geno not in selected:
            continue
        if fs not in data:
            data[fs] = {'+':[0]*(size//bin+1), '-':[0]*(size//bin+1)}
        for l in fr:
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 6:
                continue
            if ws[0] != name:
                continue
            loc = int(ws[1])
            st = ws[5]
            counts[fs][st] += 1
            # only last bin have a different bin size
            if loc // bin == size // bin:
                b = size % bin
            else:
                b = bin
            data[fs][st][loc//bin] += 1/ (b/size)
    df = []
    for fs, v in data.items():
        for st, efs in v.items():
            for i, ef in enumerate(efs):
                df.append([fs, info[fs], st, i*bin, min(i*bin+bin, size), ef/counts[fs][st]])
    df = pd.DataFrame(df, columns=['Sample', 'Celltype', 'Strand', 'Start','End','EF'])
    return df


# read order
def read_libinfo(fr, c):
    data = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        name = ws.pop(c)
        data[name] = '_'.join(ws)
    return data


# get rNMP enriched zones
def get_rez(ef, threshold, sample_threshold):
    rezs = ef[ef.EF > threshold].copy()
    # common enriched zones
    n_sample = len(ef.Sample.unique())
    n_celltype = len(ef.Celltype.unique())
    common = rezs.groupby(by=['Strand', 'Start', 'End']).agg(
        {
            'Sample':['nunique','unique'],
            'Celltype':'nunique',
            'EF':['mean','median']
        }
    ).reset_index()
    common.columns = ['Strand', 'Start', 'End', 'Sample_count', 'Sample', 'Celltype_count', 'EF_mean', 'EF_median']
    common = common[(common.Celltype_count == n_celltype) & (common.Sample_count > n_sample * sample_threshold)].copy()
    common = common.sort_values(by=['Strand', 'Start'])
    names = []
    curr = {'+':1, '-':1}
    for _, x in common.iterrows():
        st = x.Strand
        if st == '+':
            names.append(f'REZ{curr[st]}L')
        else:
            names.append(f'REZ{curr[st]}H')
        curr[st] += 1
    common['Name'] = names
    common = common[['Name', 'Strand', 'Start', 'End', 'Sample_count', 'EF_mean', 'EF_median', 'Sample']]
    return rezs, common


# draw circular plot for REZ
def draw(
        data, common, size, bin, out, order, palette, 
        tick_interval=2000, scale_pos=np.pi, draw_legend=False, no_annot=False,
        show_loc=False
    ):
    # sort and assign colors
    for c in data.Celltype.unique():
        if c not in order:
            order.append(c)
    color_pool = list(sns.color_palette(palette)) + sns.color_palette('Set2')
    colors = {x:color_pool[i] for i, x in enumerate(order)}
    orderd = {order[i]:i for i in range(len(order))}
    data['order'] = data.Celltype.map(orderd)
    data = data.sort_values(by=['order','Sample'])
    # initialize figures
    for st in ['+', '-']:
        sns.set(style='ticks', font_scale=5)
        if draw_legend:
            fig, _ = plt.subplots(figsize=(15,10), dpi=300)
            plt.subplots_adjust(left=0, right=0.7)
        else:
            fig, _ = plt.subplots(figsize=(10,10), dpi=300)
            margin = 0.07
            plt.subplots_adjust(left=margin, right=1-margin, top=1-margin, bottom=margin)
        ax = plt.subplot(111, polar=True)
        # calc width and angles
        width = np.pi * 2 / (size//bin + 1)
        angles = [x*width for x in range(size//bin + 1)] + [0]
        # each genotype
        i = 0
        grouped = set()
        for fs in data.Sample.unique():
            # prepare data for bar plot
            d = data[(data.Sample==fs) & (data.Strand==st)]
            geno = d.Celltype.unique()[0]
            height = [0] * len(angles)
            bottom = [i] * len(angles)
            for _, x in d.iterrows():
                height[x.Start//bin] += 1
            # plot forward
            if geno not in grouped:
                bars = ax.bar(angles, height, width=width, bottom=bottom, align='edge', color=colors[geno], linewidth=0.1, label=geno)
                grouped.add(geno)
            else:
                bars = ax.bar(angles, height, width=width, bottom=bottom, align='edge', linewidth=0.1, color=colors[geno])
            i += 1
        # annotation
        if not no_annot:
            common_st = common[common.Strand == st]
            for _, x in common_st.iterrows():
                loc = angles[x.Start//bin] + width / 2
                r_marker = loc * 360 / (2 * np.pi)
                ax.scatter(
                    [loc], [i+4], 
                    marker=(3, 0, -r_marker+180), color='k', s=500
                )
        # theta ticks
        ax.xaxis.set_tick_params(size=2, width=2)
        ax.set_xticks([x*tick_interval*2*np.pi/size for x in range(size//tick_interval+1)])
        ax.set_xticklabels([int(x*tick_interval/1000) for x in range(size//tick_interval+1)])
        # limits
        plt.ylim((-0.5*i, 1.25*i))
        ax.set_rorigin(-0.5*i)
        # plot a line for zero
        ax.plot([0]*50, np.linspace(0, i*1.1, 50), 'k', linewidth=2)
        # plot other small ticks
        for x in range(1, size//tick_interval+1):
            loc = x * tick_interval * np.pi * 2 / size
            ax.plot((loc, loc), (0, i*1.1), color='k', linewidth=1)
        # plot position
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        plt.setp(ax.spines.values(), linewidth=0)
        # legend
        if draw_legend:
            plt.legend(loc=(1.1,0.1))
        # other settings
        ax.grid(False)
        ax.yaxis.set_visible(False)
        if not show_loc:
            ax.xaxis.set_visible(False)
        st_name = 'light' if st == '+' else 'heavy'
        fig.savefig(f'{out}_REZ_{st_name}.png')
        plt.close('all')



def main():
    # argparse
    parser = argparse.ArgumentParser(description='Analyze enriched zones for human mtDNA')
    parser.add_argument('bed', type=argparse.FileType('r'), nargs='+', help='rNMP incorporation bed files')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Fasta index file')
    parser.add_argument('libinfo', type=argparse.FileType('r'), help='Library information')
    parser.add_argument('-o', default='enriched_zone', help='Output base name, (enriched_zones)')
    parser.add_argument('-c', type=int, default=2, help='Col num for FS number, default=2')
    parser.add_argument('-b', type=int, default=200, help='Bin size, default=200nt')
    parser.add_argument('-t', type=int, default=2000, help='Tick interval, default=2,000nt')
    parser.add_argument('--ef_threshold', type=float, default=1, help='Enrichment factor threshold for enriched regions, default=1')
    parser.add_argument('--sample_threshold', type=float, default=0.8, help='The minimum library ratio threshold of common enriched regions, default=0.8')
    parser.add_argument('--mt_name', default='chrM', help='Mitochondria name in reference genome, default=chrM')
    parser.add_argument('--selected', default=None, nargs='+', help='Selected paricular cell type(s)')
    parser.add_argument('--order', default=['CD4T', 'hESC-H9','DLTB', 'TLTB', 'WB-GTP control', 'WB-GTP PTSD', 'HCT116', 'HEK293T'], nargs='+', help='Order in hue')
    parser.add_argument('--palette', default='Dark2', help='Color paletter used to generate plots, (Dark2)')
    parser.add_argument('--legend', action='store_true', help='Draw legend')
    parser.add_argument('--no_annot', action='store_true', help='Do not annotate common REZs')
    parser.add_argument('--show_loc', action='store_true', help='Draw location indicators')
    args = parser.parse_args()
    args.c -= 1

    # mt size
    size = get_mt_size(args.fai, args.mt_name)

    # load information
    info = read_libinfo(args.libinfo, args.c)

    # load data
    data = load_bed(args.bed, size, info, args.b, args.mt_name, args.selected)

    # enriched zones
    rezs, common = get_rez(data, args.ef_threshold, args.sample_threshold)
    rezs.to_csv(f'{args.o}_rezs.tsv', sep='\t', index=False)
    common.to_csv(f'{args.o}_common_rezs.tsv', sep='\t', index=False)

    # draw
    draw(
        rezs, common, size, args.b, args.o, args.order, args.palette, 
        tick_interval=args.t, draw_legend=args.legend, no_annot=args.no_annot,
        show_loc=args.show_loc
    )

    print('Done!')

if __name__ == '__main__':
    main()
