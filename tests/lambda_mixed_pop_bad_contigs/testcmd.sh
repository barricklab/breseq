#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="expected.gd"

TESTCMD="\
    ${BRESEQ} \
    -o ${SELF} \
    -r ${DATADIR}/lambda/lambda.gbk \
  	-r ${DATADIR}/lambda/other.gbk \
  	-r ${DATADIR}/REL606/REL606.fragment.gbk \
    ${DATADIR}/lambda/lambda_mixed_population.fastq \
    "

do_test $1 ${SELF}
