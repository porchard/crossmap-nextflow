# Crossmap

This is a NextFlow implementation of the pipeline described at: https://github.com/battle-lab/crossmap

# Running

Requires NextFlow >= 20.10 and Singularity v3.

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz && zcat gencode.v26.annotation.gtf.gz > gencode.v26.annotation.gtf

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.primary_assembly.genome.fa.gz && zcat GRCh38.primary_assembly.genome.fa.gz | perl -pe 's/^(>.*?)\s+.*$/$1/' > hg38.fa

nextflow run -resume --results /path/to/results --mismatch 2 --exon_k 75 --utr_k 36 --fasta hg38.fa --gtf gencode.v26.annotation.gtf /path/to/main.nf
```


# Notes
Using GENCODE v26 and the hg38 fasta file from GENCODE (as in above commands):
* Per-gene mappability scores are generally similar to hg38 scores generated by the Battle lab (see [here](notebooks/compare-gene-mappability.ipynb)). Large differences are especially (but not solely) evident in the chrX/chrY PAR regions, suggesting that PAR masking is at least partially responsible for the differences. Other differences might be due to e.g. differences in software versions used.
* Cross-mapping scores are generally similar to hg38 scores generated by the Battle lab (see [here](notebooks/compare-crossmapping.ipynb))
* Per-SNP scores are generally similar to hg38 scores generated by the Battle lab (see [here](notebooks/compare-snp-mappability.ipynb))