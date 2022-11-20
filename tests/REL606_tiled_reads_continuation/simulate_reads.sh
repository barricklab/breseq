#!/bin/bash
gdtools APPLY -f FASTA -r ../data/REL606/REL606.fragment.gbk -o mutated.fna input.gd
breseq SIMULATE-READS -l 50 -r mutated.fna -m tiled
rm mutated.fna
