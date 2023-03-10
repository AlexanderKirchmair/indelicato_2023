#!/usr/bin/bash

#$ -S /usr/bin/bash
#$ -pe smp 8
#$ -cwd
#$ -V
#$ -N FASTQC
#$ -o logs/$JOB_NAME.log
#$ -e logs/$JOB_NAME.err

mkdir logs
mkdir data/FRDA_01_FASTQC

fastqc -o './data/FRDA_01_FASTQC' -a './tables/adapters.tsv' -t 8 ./data/FRDA_01_RAW/*fastq.gz
