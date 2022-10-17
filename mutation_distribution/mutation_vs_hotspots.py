#!/usr/bin/env python3

from collections import defaultdict
import numpy as np
import argparse
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


# load mutations
def load_vcf(frs):
    data = {}
    for fr in frs:
        for l in fr:
            if l[0] == '#':
                continue
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 6:
                continue
            loc = int(ws[1])
            if loc not in data:
                data[loc] = 0
            data[loc] += 1
    return [(x, data[x]) for x in sorted(data.keys())]


# load hotspots
def load_hotspots(fr):
    data = []
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        if len(ws) < 3:
            continue
        data.append((int(ws[0]), float(ws[4])))
    return sorted(data)


# draw circular barplot
def draw(muts, hotspots, size, out, scale_pos=1.1*np.pi, tick_interval=2000):
    # colors = {'+':'#FB5156', '-':'#5555FF'}
    sns.set(style='ticks', font_scale=2)
    fig, _ = plt.subplots(figsize=(6,6), dpi=300)
    ax = plt.subplot(111, polar=True)
    width = 30 / size * 2 * np.pi
    # limits
    mx = max([x[1] for x in hotspots])
    mx_muts = max([x[1] for x in muts])
    # plot mutations
    for_bars = ax.bar(
                    [x[0]/size*2*np.pi for x in muts], 
                    [x[1] / mx_muts * mx * 0.3 for x in muts],
                    width=width, 
                    # color=colors['+'], 
                    # edgecolor=colors['+'], 
                    color=['red'],
                    linewidth=0.2
                    )
    # plot hotspots
    for_bars = ax.bar(
                    [x[0]/size*2*np.pi for x in hotspots], 
                    [-x[1] for x in hotspots],
                    width=width, 
                    # color=colors['+'], 
                    # edgecolor=colors['+'], 
                    color=['green'],
                    linewidth=0.2
                    )
    # theta ticks
    ax.xaxis.set_tick_params(size=2, width=2)
    ax.set_xticks([x*tick_interval*2*np.pi/size for x in range(size//tick_interval+1)])
    ax.set_xticklabels([int(x*tick_interval/1000) for x in range(size//tick_interval+1)])
    # limits
    plt.ylim((-mx, mx*0.5))
    ax.set_rorigin(-mx*0.6)
    # plot a line for zero
    plt.plot((0,0), (0, mx*0.5), color='k', linewidth=4)
    # plot other small ticks
    for x in range(1, size//tick_interval+1):
        loc = x * tick_interval * np.pi * 2 / size
        plt.plot((loc, loc), (0, mx*0.5), color='k', linewidth=2)
    # draw an r scale
    ax.yaxis.set_visible(False)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    plt.setp(ax.spines.values(), linewidth=0)
    # other settings
    ax.grid(False)
    fig.savefig(out)
    plt.close('all')


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Compare the mutation and hotspots in mtDNA')
    parser.add_argument('vcf', type=argparse.FileType('r'), nargs='+', help='VCF file(s) for mtDNA mutations')
    parser.add_argument('hotspots', type=argparse.FileType('r'), help='TSV file for mtDNA rNMP hotspots')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Fasta index file')
    parser.add_argument('-o', default='mt.png', help='Output figure name, (mt.png)')
    parser.add_argument('--mt_name', default='chrM', help='Mitochondria name in reference genome (chrM)')
    args = parser.parse_args()

    # mt size
    size = get_mt_size(args.fai, args.mt_name)

    # load data
    muts = load_vcf(args.vcf)
    hotspots = load_hotspots(args.hotspots)

    # draw
    draw(muts, hotspots, size, args.o)

    print('Done!')

if __name__ == '__main__':
    main()
