#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

TESTCMD="\
    ${BRESEQ} \
    	--trim-read-ends \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${DATADIR}/lambda/lambda_mixed_population.fastq \
    "

do_test $1 ${SELF}
