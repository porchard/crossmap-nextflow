#!/usr/bin/env python3.8

import os
import sys
from Bio import SeqIO

fasta = sys.argv[1]
chroms = sys.argv[2:]

seqs = dict()

with open(fasta, 'r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        sys.stderr.write('Loaded sequence {}\n'.format(record.id))
        seqs[record.id] = record

sys.stderr.write('Finished loading sequences\n')

for record_id, record in sorted(seqs.items()):
    if record_id not in chroms:
        sys.stderr.write('Skipping sequence {}\n'.format(record_id))
        continue
    sys.stderr.write('Writing sequence {}\n'.format(record_id))
    SeqIO.write(record, sys.stdout, 'fasta')
