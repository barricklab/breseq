#!/bin/bash

#NOTE: This test must be run AFTER lambda_mixed_pop (which generates its input)

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="contingency_loci.csv"
EXPECTED_OUTPUTS[0]="expected.csv"

TESTCMD="\
    ${BRESEQ} CL-TABULATE \
    -b ${SELF}/../lambda_mixed_pop/data/reference.bam \
    -f ${SELF}/../lambda_mixed_pop/data/reference.fasta \
    -r ${SELF}/../lambda_mixed_pop/data/reference.gff3 \
    -m 5 \
    -s \
    -o ${SELF}/contingency_loci.csv \
    "

do_test $1 ${SELF}
