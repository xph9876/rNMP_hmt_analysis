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


# calculate mean coverage of each bin
def calc_coverage(frs, bin, size, mt_name):
    freqs = {'+':{}, '-':{}}
    for fr in frs:
        name = fr.name.split('.')[0].split('/')[-1]
        total = {}
        for st in '+-':
            freqs[st][name] = [0 for i in range(size//bin+1)]
            total[st] = 0
        for l in fr:
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 6:
                continue
            if ws[0] != mt_name:
                continue
            s = int(ws[1])
            e = int(ws[2])
            freq = float(ws[4])
            st = ws[5]
            if freq == 0:
                continue
            freqs[st][name][s//bin] += freq
            total[st] += freq
        for st in '+-':
            freqs[st][name] = np.asarray(freqs[st][name])
            if total[st] != 0:
                freqs[st][name] = freqs[st][name] / total[st] * len(freqs[st][name])
    data = []
    for st in freqs.keys():
        if st == '+':
            strand = 'Light'
        else:
            strand = 'Heavy'
        for k, v in freqs[st].items():
            for i in range(len(v)):
                data.append([k, strand, i*bin, v[i]])
            data.append([k, strand, len(v) * bin, v[0]])
    freqs = pd.DataFrame(data, columns=['Sample', 'Strand', 'Location', 'Frequency'])
    freqs.to_csv('test.csv')
    return freqs


# draw circular plot for REZ
def draw(
        df, size, bin, out, tick_interval=2000,
        scale_pos=np.pi, draw_legend=False
    ):
    sns.set(style='ticks', font_scale=3)
    for st in ['Light', 'Heavy']:
        data = df[df.Strand == st].copy()
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
        data['Angles'] = data['Location'] // bin * width
        # draw
        sns.lineplot(
            x='Angles', y='Frequency', hue='Sample', 
            data=data, ax=ax, linewidth=2
        )
        # theta ticks
        ax.xaxis.set_tick_params(size=2, width=2)
        ax.set_xticks([x*tick_interval*2*np.pi/size for x in range(size//tick_interval+1)])
        ax.set_xticklabels([int(x*tick_interval/1000) for x in range(size//tick_interval+1)])
        # limits
        ax.set_ylim([0, 2.3])
        ax.set_rorigin(-6)
        # plot a line for 0 and 1
        nspikes = 100
        ax.plot(np.linspace(0, 2*np.pi, nspikes), [0]*nspikes, 'k', linewidth=2)
        ax.plot(np.linspace(0, 2*np.pi, nspikes), [1]*nspikes, 'k--', linewidth=1)
        # plot ticks
        ax.plot([0]*5, np.linspace(0, 2, 5), 'k', linewidth=2)
        for x in range(1, size//tick_interval+1):
            loc = x * tick_interval * np.pi * 2 / size
            ax.plot((loc, loc), (0, 2), color='k', linewidth=1)
        # plot position
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        plt.setp(ax.spines.values(), linewidth=0)
        # legend
        if draw_legend:
            plt.legend(loc=(1.1,0.1))
        else:
            ax.get_legend().remove()
        # other settings
        ax.grid(False)
        ax.yaxis.set_visible(False)
        ax.set_xlabel('')
        ax.set_ylabel('')
        fig.savefig(f'{out}_{st}.png', transparent=True)
        plt.close('all')



def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw background coverage for DNA seq')
    parser.add_argument('bed', type=argparse.FileType('r'), nargs='+', help='Coverage bed file')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Fasta index file')
    parser.add_argument('-o', default='coverage', help='Output base name, (coverage)')
    parser.add_argument('-b', type=int, default=200, help='Bin size, default=200nt')
    parser.add_argument('-t', type=int, default=2000, help='Tick interval, default=2,000nt')
    parser.add_argument('--mt_name', default='chrM', help='Mitochondria name in reference genome, default=chrM')
    parser.add_argument('--draw_legend', action='store_true', help='Draw figure legend')
    args = parser.parse_args()

    # mt size
    size = get_mt_size(args.fai, args.mt_name)

    # calculate coverage
    coverages = calc_coverage(args.bed, args.b, size, args.mt_name)

    # draw figures
    draw(
        coverages, size, args.b, args.o, 
        tick_interval=args.t, draw_legend=args.draw_legend
    )

    print('Done!')

if __name__ == '__main__':
    main()
