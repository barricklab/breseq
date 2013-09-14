#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

TESTCMD="\
    ${BRESEQ} \
        -o ${SELF} \
        -c ${DATADIR}/lambda/lambda-contig.gbk \
        ${DATADIR}/lambda/lambda_mixed_population.fastq \
    "

do_test $1 ${SELF}
