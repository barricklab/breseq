#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../../etc/common.sh

testcmd() {
    ${BRESEQ} --no-junction \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${DATADIR}/lambda/lambda_mixed_population.1.fastq \
        ${DATADIR}/lambda/lambda_mixed_population.2.fastq \
        ${DATADIR}/lambda/lambda_mixed_population.3.fastq \
        ${DATADIR}/lambda/lambda_mixed_population.4.fastq \
        ${DATADIR}/lambda/lambda_mixed_population.5.fastq
}

do_test $1 ${SELF}
