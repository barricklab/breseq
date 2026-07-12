#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.1-2.gbk -r ${DATADIR}/lambda/lambda.3.gbk -r ${DATADIR}/lambda/lambda.4.gbk -r ${DATADIR}/lambda/lambda.5.gbk"


TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    --unzipped-html \
    -l 50 \
    ${DATADIR}/lambda/empty.fastq \
    ${DATADIR}/lambda/only_bad.fastq \
    ${DATADIR}/lambda/lambda_mixed_population.A.fastq.gz \
    ${DATADIR}/lambda/lambda_mixed_population.B.fastq.gz \
    ${DATADIR}/lambda/lambda_mixed_population.3.fastq.gz \
    ${DATADIR}/lambda/lambda_mixed_population.4.fastq.gz \
    ${DATADIR}/lambda/lambda_mixed_population.5.fastq.gz \
    "

do_test $1 ${SELF}
