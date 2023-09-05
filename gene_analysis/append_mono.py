#!/usr/bin/env python3

import argparse
import sys
import pandas as pd

def read_mono(fr):
    data = {}
    header = fr.readline().rstrip('\n').split('\t')
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        name = ws[0].split(':')[1]
        data[name] = [int(x) for x in ws[1:]]
    return data


def add_mono(x, freq, id):
    name = f'{x.Start}-{x.End}({x.Strand})'
    return freq[name][id]


def main():
    parser = argparse.ArgumentParser(description='Append mononucleotide frequency to rNMP in hmtDNA gene file')
    parser.add_argument('tsv', type=argparse.FileType('r'), help='Input rNMP count in hmtDNA gene tsv file')
    parser.add_argument('mono', type=argparse.FileType('r'), help='Input mono count file')
    parser.add_argument('-r', action='store_true', help='Append the reverse strand of input')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output to file')
    args = parser.parse_args()

    # read tsv data
    data = pd.read_csv(args.tsv, sep='\t')

    # read freq
    freq = read_mono(args.mono)

    # add columns
    if not args.r:
        data['A'] = data.apply(add_mono, axis=1, freq=freq, id=0)
        data['C'] = data.apply(add_mono, axis=1, freq=freq, id=1)
        data['G'] = data.apply(add_mono, axis=1, freq=freq, id=2)
        data['T'] = data.apply(add_mono, axis=1, freq=freq, id=3)
    else:
        data['A'] = data.apply(add_mono, axis=1, freq=freq, id=3)
        data['C'] = data.apply(add_mono, axis=1, freq=freq, id=2)
        data['G'] = data.apply(add_mono, axis=1, freq=freq, id=1)
        data['T'] = data.apply(add_mono, axis=1, freq=freq, id=0)

    data.to_csv(args.o, index=False, sep='\t')

    print('Done!')

if __name__ == '__main__':
    main()
