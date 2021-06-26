#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.tab"
EXPECTED_OUTPUTS[0]="${SELF}/expected.tab"
REFERENCE_ARG="-r ${DATADIR}/bull/bull_1.gbk"

TESTCMD="\
    ${BRESEQ} BAM2COV \
    -t \
    -a \
    --summary-json-file $SELF/summary.json \
    -f ${DATADIR}/bull/bull_1.fasta \
    -b ${DATADIR}/bull/bull_1.bam \
    -o ${SELF}/output \
    -r rachael:657-2167 \
    "

do_test $1 ${SELF}
