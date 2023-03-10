#!/usr/bin/bash

#$ -S /usr/bin/bash
#$ -pe smp 4
#$ -cwd
#$ -V
#$ -N TRIMMING
#$ -o logs/$JOB_NAME-$JOB_ID.log
#$ -e logs/$JOB_NAME-$JOB_ID.err

in=$(readlink -f $1)
out=$(basename $in)
out=$2/$out

cutadapt -q 20 -O 10 -e 0.15 -m 10 -j 4 -a A{100} --action=trim -o $out $in
