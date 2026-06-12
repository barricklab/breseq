#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda_bad_orfs.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -p \
    --no-junction-prediction \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq \
    "

do_test $1 ${SELF}
