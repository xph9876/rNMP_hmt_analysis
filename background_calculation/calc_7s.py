#!/usr/bin/env python3

import pandas as pd
import argparse
import sys


# calculate mean coverage of each bin
def calc_coverage(frs):
    freqs = {}
    for fr in frs:
        name = fr.name.split('/')[-1].split('.')[0].replace('-chrm', '').replace('DNA-seq-', '')
        freqs[name] = {
            '7s':0, 'Total':0, 'Other':0
        }
        for l in fr:
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 3:
                continue
            freqs[name]['Total'] += int(ws[4])
            loc = int(ws[2])
            # 7s DNA
            if (loc <= 191 or loc > 16106):
                freqs[name]['7s'] += int(ws[4])
            else:
                freqs[name]['Other'] += int(ws[4])
    data = []
    for name, v in freqs.items():
        data.append([
            name, v['Total'], v['7s'], v['Other']
        ])
    df = pd.DataFrame(data, columns=['Sample', 'Total', '7S','Other'])
    df['7S_portion'] = df['7S'] / df['Total']
    df['7S_EF'] = df['7S'] / df['Total']/((191+16569-16106)/16569)
    df['Other_portion'] = df['Other'] / df['Total']
    return df


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
