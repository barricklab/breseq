#!/bin/bash
gdtools APPLY -r ../data/lambda/lambda.5.gbk -f FASTA -o output.fasta input.gd
mason illumina -s 1 -sq -hm 1 -hs 0 -hi 0 -hnN --no-N -n 50 -N 5000 -o ../data/lambda/lambda.short_sequence_repeats.fastq output.fasta
