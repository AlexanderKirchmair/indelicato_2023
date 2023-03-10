#!/usr/bin/bash

#$ -S /usr/bin/bash
#$ -pe smp 1
#$ -cwd
#$ -V
#$ -N SNPCHECK
#$ -o logs/$JOB_NAME.log
#$ -e logs/$JOB_NAME.err

# NGSCheckMate
# git clone https://github.com/parklab/NGSCheckMate.git
# edit ncm.conf

in=$(readlink -f $1)
out=$(readlink -f $2)
bed=$(readlink -f ~/tools/NGSCheckMate/SNP/SNP_GRCh38_hg38_wChr.bed)
snp=$(readlink -f ~/tools/NGSCheckMate/SNP/SNP.pt)

# BAM/VCF mode
NGSCheckMate=$(readlink -f ~/tools/NGSCheckMate/ncm.py)
python2.7 $NGSCheckMate -B -l $in -bed $bed -O $out

# FASTQ mode
# NGSCheckMate=$(readlink -f ~/tools/NGSCheckMate/ncm_fastq.py)
# python2.7 $NGSCheckMate -l $in -pt $snp -O $out -p 1
