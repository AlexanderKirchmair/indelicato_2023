#!/bin/bash

mkdir data/FRDA_02_TRIMMED

for file in data/FRDA_01_RAW/*fastq.gz
do
  qsub lib/run_trimming.sh $file data/FRDA_02_TRIMMED
done

