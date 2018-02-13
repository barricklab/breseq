#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk -r ${DATADIR}/lambda/other.gbk -r ${DATADIR}/REL606/REL606.fragment.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq \
    "

do_test $1 ${SELF}
