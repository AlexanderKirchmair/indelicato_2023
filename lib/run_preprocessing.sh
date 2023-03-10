#!/bin/bash
# Preprocessing of RNA-seq data using the nf-core/rnaseq pipeline
# Usage: bash -i run_rnaseq.sh ['samplesheet.csv'] ['resultsdir'] ['tmpdir']

## Input (samplesheet)
in=$(readlink -f $1)

## Output (results dir)
out=$(readlink -f $2)

## Pipeline options
wd=$(readlink -f $3)
config=$(readlink -f ./preprocessing/nf.config)

## Genome
fasta=$(readlink -f ./data/genome/gencode.GRCh38.fa)
gtf=$(readlink -f ./data/genome/gencode.v38.gtf)
star_idx=$(readlink -f ./data/genome/gencode.v38.star2.7.9a)

nextflow run nf-core/rnaseq -r 3.6 \
  -profile singularity,frda \
  -w $wd \
  -c $config \
  --input $in \
  --outdir $out \
  --fasta $fasta \
  --gtf $gtf \
  --star_index $star_idx \
  --aligner 'star_salmon' \
  --gencode \
  --skip_trimming \
  --remove_ribo_rna false
