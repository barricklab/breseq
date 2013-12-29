#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="expected.gd"

TESTCMD="\
    ${BRESEQ} \
    	--polymorphism-prediction \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.1-2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -r ${DATADIR}/lambda/lambda.4.gbk \
        -r ${DATADIR}/lambda/lambda.5.gbk \
        ${DATADIR}/lambda/lambda_mixed_population.1.fastq \
    	${DATADIR}/lambda/lambda_mixed_population.2.fastq \
        ${DATADIR}/lambda/lambda_mixed_population.3.fastq \
        ${DATADIR}/lambda/lambda_mixed_population.4.fastq \
        ${DATADIR}/lambda/lambda_mixed_population.5.fastq \
    "

do_test $1 ${SELF}
