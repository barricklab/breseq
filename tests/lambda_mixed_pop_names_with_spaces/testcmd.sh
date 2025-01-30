#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output output/output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
#REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

TESTCMD=" \
    cp ${DATADIR}/lambda/lambda.gbk \"${SELF}/lambda space.gbk\"  ; \
    cp ${DATADIR}/lambda/lambda_mixed_population.fastq \"${SELF}/lambda mixed population space.fastq\" ; \
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF}/output\ output \
    -g ${SELF}/header.gd \
    --genbank-field-for-seq-id version \
    -r \"${SELF}/lambda space.gbk\" \
    ${SELF}/lambda\ mixed\ population\ space.fastq ;\
    rm \"${SELF}/lambda space.gbk\"  ; \
    rm \"${SELF}/lambda mixed population space.fastq\" ; \
    "

do_test $1 ${SELF}
