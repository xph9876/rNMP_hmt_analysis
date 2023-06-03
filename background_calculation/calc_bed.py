#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys


# calculate mean coverage of each bin
def calc_coverage(frs):
    freqs = {}
    for fr in frs:
        name = fr.name.split('/')[-1].split('.')[0].replace('-chrm', '').replace('DNA-seq-', '')
        freqs[name] = {'7s':0, 'Total':0, 'Other':0}
        for l in fr:
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 3:
                continue
            freqs[name]['Total'] += int(ws[4])
            loc = int(ws[2])
            if ws[5] == '-' and (loc <= 191 or loc > 16106):
                freqs[name]['7s'] += int(ws[4])
            else:
                freqs[name]['Other'] += int(ws[4])
    data = []
    for name, v in freqs.items():
        data.append([name, v['Total'], v['7s'], v['Other']])
    df = pd.DataFrame(data, columns=['Sample', 'Total', '7S','Other'])
    df['7S_portion'] = df['7S'] / df['Total']
    df['7S_EF'] = df['7S'] / df['Total']/((191+16569-16106)/(16569*2))
    df['Other_portion'] = df['Other'] / df['Total']
    return df


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
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    args = parser.parse_args()

    # calculate coverage
    coverages = calc_coverage(args.bed)
    coverages.to_csv(args.o, sep='\t', index=False)


    print('Done!')

if __name__ == '__main__':
    main()
