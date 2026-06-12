#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output output/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
#REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

TESTCMD=" \
    ln -s \"$(cd ${DATADIR}/lambda && pwd)/lambda.gbk\" \"${SELF}/lambda space.gbk\"  ; \
    ln -s \"$(cd ${DATADIR}/lambda && pwd)/lambda_mixed_population.fastq\" \"${SELF}/lambda mixed population space.fastq\" ; \
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF}/output\ output \
    -g ${SELF}/header.gd \
    --genbank-field-for-seq-id version \
    -r \"${SELF}/lambda space.gbk\" \
    ${SELF}/lambda\ mixed\ population\ space.fastq ;\
    rm \"${SELF}/lambda space.gbk\"  ; \
    rm \"${SELF}/lambda mixed population space.fastq\" ; \
    if ls -d \"${SELF}/output output\"/0[1-9]_* >/dev/null 2>&1; then \
        echo 'Intermediate stage directories were not cleaned up!'; \
        exit 1; \
    fi \
    "

do_test $1 ${SELF}
