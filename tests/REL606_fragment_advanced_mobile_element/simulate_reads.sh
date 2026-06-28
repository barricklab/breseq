#!/bin/bash
gdtools APPLY -r ../data/REL606/REL606.fragment.gbk -f FASTA -o output.fasta input.gd
mason illumina -s 1 -sq -hm 1 -hs 0 -hi 0 -hnN --no-N -n 36 -N 40000 -o sim_reads.fastq output.fasta
gzip -f sim_reads.fastq
mv sim_reads.fastq.gz ../data/REL606/REL606.advanced_mobile_element.fastq.gz
