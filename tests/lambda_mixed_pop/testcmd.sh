#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    -g ${SELF}/header.gd \
    --genbank-field-for-seq-id VERSION \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq \
    "

do_test $1 ${SELF}
