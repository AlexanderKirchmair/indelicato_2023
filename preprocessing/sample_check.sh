#!/usr/bin/bash

mkdir samplecheck

ls -d $PWD/data/FRDA_03_NF_RESULTS/star_salmon/*bam > samplecheck/files.txt
qsub lib/run_NGSCheckMate.sh samplecheck/files.txt samplecheck
mv r_script.r.Rout samplecheck/
