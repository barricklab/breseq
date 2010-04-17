#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../../etc/common.sh

testcmd() {
    ${BRESEQ} --no-junction \
        -o ${SELF} \
        -r ${SELF}/data/lambda.gbk \
        ${SELF}/data/lambda_mixed_population.fastq
}

do_test $1 ${SELF}
