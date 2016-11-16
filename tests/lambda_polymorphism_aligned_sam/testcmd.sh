#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -p \
    --aligned-sam \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.sam \
    "

do_test $1 ${SELF}
