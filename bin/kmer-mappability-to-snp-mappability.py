#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument('--k', type=int, required=True)
parser.add_argument('kmer_mappability')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

KMER_MAPPABILITY = args.kmer_mappability # '/net/topmed10/working/porchard/rnaseq/work/test-crossmap/work/60/ea9a1e7f64ba860c7616ddf035ddc4/out/mappability_75mer_2mismatch.bed'
K = args.k
SPAN = (K * 2) - 1


class ChunkedFileReader:

    def __init__(self, fh, field):
        assert(isinstance(field, int))
        self.fh = fh
        self.field = field
        self.current_value = None
        self.next_set = []
        self.eof = False


    def __next__(self):
        # read until eof or until the field changes
        lines = self.next_set.copy()
        self.next_set = []
        if self.eof:
            raise StopIteration

        while True:
            line = self.fh.readline()
            if line == '':
                self.eof = True
                return lines
            line = line.rstrip()
            new_value = line.split()[self.field]
            if self.current_value is None or new_value == self.current_value:
                self.current_value = new_value
                lines.append(line)
            else:
                self.next_set.append(line)
                self.current_value = new_value
                return lines

    def __iter__(self):
        return self



with open(KMER_MAPPABILITY, 'r') as fh:
    for chunk in ChunkedFileReader(fh, 0):
        pos = []
        value = []
        chrom = chunk[0].split('\t')[0]
        logging.info(f'Processing chrom {chrom}')
        for i in chunk:
            chrom, start, end, score = i.split('\t')
            start, end, score = int(start), int(end), float(score)
            for p in range(start, end):
                pos.append(p)
                value.append(score)
        chr_series = pd.Series(value, index=pos)
        pos_scores = chr_series.rolling(SPAN, min_periods=0, center=True, closed='both').mean()
        for i, v in zip(pos_scores.index, pos_scores.values):
            print(f'{chrom}\t{i}\t{i+1}\t{v}')

