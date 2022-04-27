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


# load rNMP coordinates
def load_bed(bed, size, bin, name):
    data = {'+':[0]*(size//bin+1), '-':[0]*(size//bin+1)}
    for l in bed:
        ws = l.rstrip('\n').split('\t')
        if len(ws) < 6:
            continue
        if ws[0] != name:
            continue
        loc = int(ws[1])
        st = ws[5]
        data[st][loc//bin] += 1
    return data


# draw circular barplot
def draw(data, size, bin, out, tick_interval=2000, scale_pos=1.1*np.pi):
    colors = {'+':'#FB5156', '-':'#5555FF'}
    sns.set(style='ticks', font_scale=3)
    fig, _ = plt.subplots(figsize=(20,20))
    ax = plt.subplot(111, polar=True)
    # calc width and angles
    width = np.pi * 2 / len(data['+'])
    angles = [x*width for x in range(len(data['+']))]
    # plot forward
    for_bars = ax.bar(angles, data['+'], width=width, color=colors['+'], edgecolor=colors['+'], linewidth=0.2)
    # plot reverse
    rev = [-x for x in data['-']]
    rev_bars = ax.bar(angles, rev, width=width, color=colors['-'], edgecolor=colors['-'], linewidth=0.2)
    # theta ticks
    ax.xaxis.set_tick_params(size=2, width=2)
    ax.set_xticks([x*tick_interval*2*np.pi/size for x in range(size//tick_interval+1)])
    ax.set_xticklabels([int(x*tick_interval/1000) for x in range(size//tick_interval+1)])
    # limits
    mx = max(max(data['+']), max(data['-']))
    plt.ylim((-mx/2, mx/2))
    ax.set_rorigin(-mx*0.8)
    # plot a line for zero
    plt.plot((0,0), (0, mx*1.2), color='k', linewidth=2)
    # plot other small ticks
    for x in range(1, size//tick_interval+1):
        loc = x * tick_interval * np.pi * 2 / size
        plt.plot((loc, loc), (0, mx*1.2), color='k', linewidth=1)
    # draw an r scale
    ax.yaxis.set_visible(False)
    # plt.plot((scale_pos - 3/180*np.pi, scale_pos + 3/180*np.pi), (max(data['+']), max(data['+'])), color='k', linewidth=1.5)
    # plt.plot((scale_pos - 5/180*np.pi, scale_pos + 5/180*np.pi), (-max(data['-']), -max(data['-'])), color='k', linewidth=1.5)
    # plt.plot((scale_pos, scale_pos), (max(data['+']), -max(data['-'])), color='k', linewidth=1.5)
    # plt.text(scale_pos, max(data['+'])*1.05, str(max(data['+'])), ha='center', va='top', rotation=scale_pos/np.pi*180+180, rotation_mode='anchor')
    # plt.text(scale_pos, -max(data['-'])*0.95, str(max(data['-'])), ha='center', va='bottom', rotation=scale_pos/np.pi*180+180, rotation_mode='anchor')
    # plot position
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    plt.setp(ax.spines.values(), linewidth=0)
    # other settings
    ax.grid(False)
    fig.savefig(out)
    plt.close('all')



def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw circular barplot for rNMP incorporation in mitonchondrial dNA')
    parser.add_argument('bed', type=argparse.FileType('r'), help='rNMP incorporation bed file')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Fasta index file')
    parser.add_argument('-o', default='mt.png', help='Output figure name, (mt.png)')
    parser.add_argument('-b', type=int, default=100, help='Bin size, default=100nt')
    parser.add_argument('-t', type=int, default=2000, help='Tick interval, default=2,000nt')
    parser.add_argument('--mt_name', default='chrM', help='Mitochondria name in reference genome, default=chrM')
    args = parser.parse_args()

    # mt size
    size = get_mt_size(args.fai, args.mt_name)

    # load data
    data = load_bed(args.bed, size, args.b, args.mt_name)

    # draw
    draw(data, size, args.b, args.o, tick_interval=args.t)

    print('Done!')

if __name__ == '__main__':
    main()
