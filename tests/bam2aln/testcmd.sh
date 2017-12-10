#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output.txt"
EXPECTED_OUTPUTS[0]="expected.txt"
REFERENCE_ARG="-r ${DATADIR}/bull/bull_1.gbk"

TESTCMD="\
    ${BRESEQ} BAM2ALN \
    --format TXT \
    -f ${DATADIR}/bull/bull_1.fasta \
    -b ${DATADIR}/bull/bull_1.bam \
    -n 50 \
    -o ${SELF}/${CURRENT_OUTPUTS[0]} \
    -r rachael:57-167 \
    "

do_test $1 ${SELF}
