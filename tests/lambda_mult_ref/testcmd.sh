#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../../etc/common.sh

testcmd() {
    ${BRESEQ} \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.1.gbk \
        -r ${DATADIR}/lambda/lambda.2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -r ${DATADIR}/lambda/lambda.4.gbk \
        -r ${DATADIR}/lambda/lambda.5.gbk \
        ${DATADIR}/lambda/lambda_mixed_population.fastq
}

do_test $1 ${SELF}
