#!/usr/bin/env python3.8

import os
import sys
from Bio import SeqIO

fasta = sys.argv[1]

seqs = dict()

with open(fasta, 'r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        sys.stderr.write('Loaded sequence {}\n'.format(record.id))
        seqs[record.id] = record

sys.stderr.write('Finished loading sequences\n')

for record_id, record in sorted(seqs.items()):
    sys.stderr.write('Writing sequence {}\n'.format(record_id))
    CHRFASTA = f'{record_id}.fa'
    with open(CHRFASTA, 'w') as f:
        SeqIO.write(record, f, 'fasta')
