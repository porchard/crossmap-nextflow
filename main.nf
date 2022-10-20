#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// KMER MAPPABILITY PARAMS
FULL_FASTA = params.fasta
GENE_ANNOT_FN = params.gtf // must be unzipped.
MISMATCH = params.mismatch // 2
EXON_K = params.exon_k // 76
UTR_K = params.utr_k // 36

// Per-gene score params
KS = (EXON_K == UTR_K) ? [EXON_K] : [EXON_K, UTR_K]

def count_genes_in_gtf (f) {
    gene_count = 0
    no_gene_count = 0
    file(f).eachLine {  str -> (!str.contains('^#') && str.contains('\tgene\t')) ? gene_count++ : no_gene_count++}
    return gene_count
}

// probably no need to change if have reasonable resources available
MAX_CHR=7 // maximum number of chromosomes to load in memory at a time
MAX_GENE_ALIGNMENT=200 // maximum number of genes to align before cleaning alignments
DIR_NAME_LEN=12 // length of the sub-directory names. First dir_len letters from gene
N_GENE_PER_CROSSMAP_BATCH = 200

process split_fasta {

    container 'library://porchard/default/crossmap:20220311'
    memory '10 GB'

    input:
    path(fasta)

    output:
    path("*.fa")

    """
    mkdir split-fasta
    split-fasta.py $fasta
    """

}

process make_k_mappability {

    cpus 16
    time '48h'
    container 'library://porchard/default/crossmap:20220311'

    input:
    path(full_fasta)
    each k

    output:
    tuple val(k), path("out")

    script:
    gem_index_pref = 'ref_gem_index'
    gem_mappability_pref = "mappability_" + k + "mer_" + MISMATCH + "mismatch"

    """
    # index the reference genome (execute once for one fasta file)
    mkdir -p out
    gem-indexer -T ${task.cpus} -c dna -i $full_fasta -o out/${gem_index_pref}
    # compute mappability -- approx time: ~6 hr
    gem-mappability -m $MISMATCH -T ${task.cpus} -I out/${gem_index_pref}.gem -l $k -o out/${gem_mappability_pref}
    """

}

process make_bed {

    time '48h'
    publishDir "${params.results}/kmer-mappability"
    container 'library://porchard/default/crossmap:20220311'
    memory '50 GB'

    input:
    tuple val(k), path(outdir)

    output:
    tuple val(k), path("${outdir}/mappability_${k}mer_${MISMATCH}mismatch.bed")

    script:
    gem_index_pref = 'ref_gem_index'
    gem_mappability_pref = "mappability_" + k + "mer_" + MISMATCH + "mismatch"

    """
    # convert mappability file to bed, step by step
    gem-2-wig -I ${outdir}/${gem_index_pref}.gem -i ${outdir}/${gem_mappability_pref}.mappability -o ${outdir}/${gem_mappability_pref}
    wigToBigWig ${outdir}/${gem_mappability_pref}.wig ${outdir}/${gem_mappability_pref}.sizes ${outdir}/${gem_mappability_pref}.bigWig
    bigWigToBedGraph ${outdir}/${gem_mappability_pref}.bigWig ${outdir}/${gem_mappability_pref}.bed
    """

}


process make_bowtie_index {

    memory '50 GB'
    container 'library://porchard/default/crossmap:20220311'

    input:
    path(fasta)

    output:
    path("genome*.ebwt")

    """
    /opt/bowtie-1.3.1-linux-x86_64/bowtie-build $fasta genome
    """

}


process process_annotation_data {

    memory '15 GB'
    container 'library://porchard/default/crossmap:20220311'

    input:
    path(gene_annot_fn)

    output:
    path('annot/annot.exon_utr.txt')


    """
    mkdir -p annot
    mkdir -p logs
    gtf_to_txt.R --gtf $GENE_ANNOT_FN -f 'exon,UTR' -o annot/annot.exon_utr.txt 2>&1 | tee logs/gtf_to_txt.log
    """

}

process kmer_mappability_to_gene_mappability {

    memory '30 GB'
    container 'library://porchard/default/crossmap:20220311'
    publishDir "${params.results}/gene-mappability"

    input:
    path(exon_utr_annot_fn)
    path(exon_kmer_mappability_fn)
    path(utr_kmer_mappability_fn)

    output:
    path('gene_mappability/gene_mappability.txt')

    """
    mkdir -p gene_mappability
    mkdir -p logs
    compute_mappability.R --annot $exon_utr_annot_fn --k_exon $EXON_K --k_utr $UTR_K --kmap_exon $exon_kmer_mappability_fn --kmap_utr $utr_kmer_mappability_fn --verbose 1 -o gene_mappability/gene_mappability.txt 2>&1 | tee logs/compute_mappability.log
    """

}

process generate_ambiguous_kmers {

    memory '30 GB'
    container 'library://porchard/default/crossmap:20220311'

    input:
    path(mappability_fn)
    path(exon_utr_annot_fn)
    path(exon_kmer_mappability_fn)
    path(utr_kmer_mappability_fn)
    path("split_fasta/*")

    output:
    path('ambiguous_kmer_dir'), emit: kmer_dir
    path("ambiguous_kmer_dir/*/*.kmer.txt"), emit: kmer_txt

    script:
    mappability_th1 = 0
    mappability_th2 = 1
    
    """
    mkdir ambiguous_kmer_dir
    mkdir -p logs
    generate_ambiguous_kmers.R --mappability $mappability_fn --genome split_fasta --annot $exon_utr_annot_fn --k_exon $EXON_K --k_utr $UTR_K --kmap_exon $exon_kmer_mappability_fn --kmap_utr $utr_kmer_mappability_fn --th1 $mappability_th1 --th2 $mappability_th2 --dir_name_len $DIR_NAME_LEN --verbose 1 -o ambiguous_kmer_dir 2>&1 | tee logs/ambiguous-kmers.log
    """

}

process make_ambiguous_kmer_fastas {

    container 'library://porchard/default/crossmap:20220311'
    memory '10 GB'

    input:
    path(x)

    output:
    path(x)

    """
    for fn in ${x}/*/*.kmer.txt
    do
        fasta_fn=\$(echo \$fn | sed 's/.txt\$/.fa/g')
        awk -v i=-1 '{i += 1 ; print ">"i ; print}' < \$fn > \$fasta_fn
    done
    """

}

process preprocess_cross_mappability {

    container 'library://porchard/default/crossmap:20220311'
    memory '30 GB'

    input:
    path(mappability_fn)
    path(exon_utr_annot_fn)
    path(ambiguous_kmer_dir)
    path(bowtie_index)

    output:
    tuple path(mappability_fn), path(exon_utr_annot_fn), path(ambiguous_kmer_dir), path(bowtie_index), path('ambiguous_kmers_alignment'), path('cross_mappability')

    """
    mkdir -p ambiguous_kmers_alignment
    mkdir -p cross_mappability
    compute_cross_mappability.R --annot $exon_utr_annot_fn --mappability $mappability_fn --kmer $ambiguous_kmer_dir --align ambiguous_kmers_alignment --index genome --n1 1 --n2 $N_GENE_PER_CROSSMAP_BATCH --mismatch $MISMATCH --max_chr $MAX_CHR --max_gene $MAX_GENE_ALIGNMENT --initonly TRUE --dir_name_len $DIR_NAME_LEN --verbose 1 -o cross_mappability 2>&1 | tee compute_cross_mappability_1_init.log
    """

}

process process_cross_mappability {

    container 'library://porchard/default/crossmap:20220311'
    memory '30 GB'
    tag "$n1"
    errorStrategy 'ignore'

    input:
    tuple path(mappability_fn), path(exon_utr_annot_fn), path(ambiguous_kmer_dir), path(bowtie_index), path(ambiguous_kmers_alignment), path(cross_mappability)
    each n1

    output:
    path("crossmap.txt")

    script:
    n2 = n1 + N_GENE_PER_CROSSMAP_BATCH - 1

    """
    compute_cross_mappability.R --annot $exon_utr_annot_fn --mappability $mappability_fn --kmer $ambiguous_kmer_dir --align $ambiguous_kmers_alignment --index genome --n1 $n1 --n2 $n2 --mismatch $MISMATCH --max_chr $MAX_CHR --max_gene $MAX_GENE_ALIGNMENT --initonly FALSE --dir_name_len $DIR_NAME_LEN --verbose 1 -o $cross_mappability
    cat *.crossmap.txt > crossmap.txt
    """

}

process concat_cross_mappability {

    publishDir "${params.results}/cross_mappability"

    input:
    path("*.crossmap.txt")

    output:
    path('crossmap.txt')

    """
    cat *.crossmap.txt > crossmap.txt
    """

}


workflow {

    full_fasta = Channel.fromPath(FULL_FASTA)
    gene_annot_fn = Channel.fromPath(GENE_ANNOT_FN)
    
    split_fastas = split_fasta(full_fasta).flatten().toSortedList()
    bowtie_index = make_bowtie_index(full_fasta)

    kmer_mappability = make_k_mappability(full_fasta, KS) | make_bed // k, bed
    exon_kmer_mappability_fn = kmer_mappability.filter({it -> it[0] == EXON_K}).map({it -> it[1]})
    utr_kmer_mappability_fn = kmer_mappability.filter({it -> it[0] == UTR_K}).map({it -> it[1]})

    annot_txt = process_annotation_data(gene_annot_fn)
    gene_mappability = kmer_mappability_to_gene_mappability(annot_txt, exon_kmer_mappability_fn, utr_kmer_mappability_fn)
    ambig_kmers = generate_ambiguous_kmers(gene_mappability, annot_txt, exon_kmer_mappability_fn, utr_kmer_mappability_fn, split_fastas)
    
    ambig_kmers_with_fasta = make_ambiguous_kmer_fastas(ambig_kmers.kmer_dir)

    NUMBER_GENES = count_genes_in_gtf(GENE_ANNOT_FN)
    n1 = (1..NUMBER_GENES).step(N_GENE_PER_CROSSMAP_BATCH)
    
    process_cross_mappability(preprocess_cross_mappability(gene_mappability, annot_txt, ambig_kmers_with_fasta, bowtie_index), n1).toSortedList() | concat_cross_mappability

}