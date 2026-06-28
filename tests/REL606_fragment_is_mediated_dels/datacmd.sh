#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
GENBANK=${SELF}/../data/REL606/REL606.fragment.gbk
MUTATED=${SELF}/REL606.mutated.fna
OUTPUT=${SELF}/../data/REL606/REL606.fragment.3.fastq.gz
GENOMEDIFF=${SELF}/input.gd

## Apply all mutations: two existing IS-mediated DELs plus new double IS deletion
${TESTBINPREFIX}/gdtools APPLY -o ${MUTATED} -r ${GENBANK} ${GENOMEDIFF}

## Simulate reads: ~35x coverage, 36 bp, no errors/mutations, single-end
## Combined mutated genome ≈ 23600 bp → N ≈ 35 * 23600 / 36 ≈ 22944
wgsim -e 0 -r 0 -R 0 -X 0 -N 22944 -1 36 -2 36 ${MUTATED} /dev/stdout /dev/null | gzip > ${OUTPUT}
