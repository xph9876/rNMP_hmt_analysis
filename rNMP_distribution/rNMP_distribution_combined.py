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
def load_bed(bed, size, info, counts, bin, name, selected):
    data = {}
    geno_counts = {k:0 for k in info.values()}
    for fr in bed:
        fs = fr.name.split('/')[-1].split('.')[0]
        if fs not in info or fs not in counts:
            print(f'Library {fs} doesn\'t have related information or counts, skipped')
            continue
        geno = info[fs]
        if geno not in selected:
            continue
        geno_counts[geno] += 1
        if geno not in data:
            data[geno] = {'+':[0]*(size//bin+1), '-':[0]*(size//bin+1)}
        for l in fr:
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 6:
                continue
            if ws[0] != name:
                continue
            loc = int(ws[1])
            st = ws[5]
            data[geno][st][loc//bin] += 1/counts[fs]*size/bin
    for geno in data.keys():
        for st in data[geno].keys():
            for i in range(len(data[geno][st])):
                data[geno][st][i] /= geno_counts[geno]
    return data


# read order
def read_libinfo(fr, c):
    data = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        name = ws.pop(c)
        data[name] = '_'.join(ws)
    return data


# read mito count
def read_mito_count(fr):
    data = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        if len(ws) != 2:
            continue
        data[ws[0]] = float(ws[1])
    return data


# draw circular barplot
def draw(data, size, bin, out, tick_interval=2000, scale_pos=np.pi):
    colors = {x:sns.color_palette('Set1')[i] for i, x in enumerate(data)}
    linestyle = {'+':'-', '-':'-'}
    sns.set(style='ticks', font_scale=2.25)
    fig, _ = plt.subplots(figsize=(15,10))
    plt.subplots_adjust(left=0, right=0.7)
    ax = plt.subplot(111, polar=True)
    # calc width and angles
    width = np.pi * 2 / (size//bin + 1)
    angles = [x*width for x in range(size//bin + 1)] + [0]
    # each genotype
    for geno, d in data.items():
        # plot forward
        for_bars = ax.plot(angles, d['+']+[d['+'][0]], linestyle=linestyle['+'], color=colors[geno], label=geno)
        # plot reverse
        rev = [-x for x in d['-']+[d['-'][0]]]
        rev_bars = ax.plot(angles, rev, linestyle=linestyle['-'], color=colors[geno])
    # theta ticks
    ax.xaxis.set_tick_params(size=2, width=2)
    ax.set_xticks([x*tick_interval*2*np.pi/size for x in range(size//tick_interval+1)])
    ax.set_xticklabels([int(x*tick_interval/1000) for x in range(size//tick_interval+1)])
    # limits
    mx = max(max(max(d['+']) for d in data.values()), max(max(d['-']) for d in data.values()))
    plt.ylim((-mx, mx))
    # ax.set_rorigin(-2*mx)
    # plot a line for zero
    ax.plot([0]*50, np.linspace(0, mx*1.2, 50), 'k', linewidth=2)
    # plot other small ticks
    for x in range(1, size//tick_interval+1):
        loc = x * tick_interval * np.pi * 2 / size
        ax.plot((loc, loc), (0, mx*1.2), color='k', linewidth=1)
    # draw an r scale
    ax.yaxis.set_visible(False)
    ax.plot(np.linspace(0, np.pi*2, num=50), [0]*50, color='k', linewidth=3)
    ax.plot(np.linspace(0, np.pi*2, num=50), [1]*50, color='k', linestyle='dotted', linewidth=1.5)
    ax.plot(np.linspace(0, np.pi*2, num=50), [-1]*50, color='k', linestyle='dotted', linewidth=1.5)
    plt.plot((scale_pos - 3/180*np.pi, scale_pos + 3/180*np.pi), (-1, -1), color='k', linewidth=1.5)
    plt.plot((scale_pos - 5/180*np.pi, scale_pos + 5/180*np.pi), (1, 1), color='k', linewidth=1.5)
    plt.plot((scale_pos, scale_pos), (-1, 1), color='k', linewidth=1.5)
    plt.text(scale_pos, 1*1.05, '1', ha='center', va='top', rotation=scale_pos/np.pi*180+180, rotation_mode='anchor')
    plt.text(scale_pos, -1*0.95, '1', ha='center', va='bottom', rotation=scale_pos/np.pi*180+180, rotation_mode='anchor')
    # plot position
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    plt.setp(ax.spines.values(), linewidth=0)
    plt.legend(loc=(1.1,0.1))
    # other settings
    ax.grid(False)
    fig.savefig(out)
    plt.close('all')



def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw combined circular barplot for rNMP incorporation in mitonchondrial DNA')
    parser.add_argument('bed', type=argparse.FileType('r'), nargs='+', help='rNMP incorporation bed files')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Fasta index file')
    parser.add_argument('libinfo', type=argparse.FileType('r'), help='Library information')
    parser.add_argument('count', type=argparse.FileType('r'), help='rNMP count for each library')
    parser.add_argument('-o', default='mt.png', help='Output figure name, (mt.png)')
    parser.add_argument('-c', type=int, default=2, help='Col num for FS number, default=2')
    parser.add_argument('-b', type=int, default=100, help='Bin size, default=100nt')
    parser.add_argument('-t', type=int, default=2000, help='Tick interval, default=2,000nt')
    parser.add_argument('--mt_name', default='chrM', help='Mitochondria name in reference genome, default=chrM')
    parser.add_argument('--selected', default=['CD4T', 'HEK293T', 'hESC-H9', 'DLTB', 'TLTB'], nargs='+', help='Selected genotypes, (All WT > 3)')
    args = parser.parse_args()
    args.c -= 1

    # mt size
    size = get_mt_size(args.fai, args.mt_name)

    # load information
    info = read_libinfo(args.libinfo, args.c)

    # load mito count
    counts = read_mito_count(args.count)

    # load data
    data = load_bed(args.bed, size, info, counts, args.b, args.mt_name, args.selected)

    # draw
    draw(data, size, args.b, args.o, tick_interval=args.t)

    print('Done!')

if __name__ == '__main__':
    main()
