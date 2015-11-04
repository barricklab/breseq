#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

TESTCMD="\
    ${BRESEQ} \
    	--polymorphism-prediction \
    	--polymorphism-frequency-cutoff 0.2 \
    	-b 20 \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        ${DATADIR}/lambda/lambda_mixed_population.fastq \
    "

do_test $1 ${SELF}
