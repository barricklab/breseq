#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda-contig.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    -c ${DATADIR}/lambda/lambda-contig.gbk \
    ${DATADIR}/lambda/lambda_mixed_population.fastq \
    "

do_test $1 ${SELF}
