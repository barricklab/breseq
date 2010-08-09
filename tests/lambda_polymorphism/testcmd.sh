#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

testcmd() {
    ${BRESEQ} \
    	--polymorphism-prediction \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.gbk \
        --polymorphism-frequency-cutoff=0.25 \
        ${DATADIR}/lambda/lambda_mixed_population.fastq
}

do_test $1 ${SELF}
