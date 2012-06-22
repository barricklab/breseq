#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
GENBANK=${SELF}/../data/REL606/REL606.fragment.gbk
MUTATED=${SELF}/REL606.mutated.fna
OUTPUT=${SELF}/REL606.fragment.2.fastq
GENOMEDIFF=${SELF}/input.gd


## This script runs commands to generate the fastq data for this test
## Add new mutations to make the test more comprehensive
genomediff MUTATE -r ${GENBANK} -o ${MUTATED} ${GENOMEDIFF}
fastq_utils SIMULATE -l 36 -c 35 -o ${OUTPUT} ${MUTATED}

