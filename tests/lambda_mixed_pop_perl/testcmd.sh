#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

testcmd() {
    ${BRESEQ} \
    	--perl-calc-trims \
    	--perl-identify-mutations \
    	--perl-error-count \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${DATADIR}/lambda/lambda_mixed_population.fastq
}

do_test $1 ${SELF}