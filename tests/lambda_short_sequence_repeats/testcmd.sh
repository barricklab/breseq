#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="expected.gd"

TESTCMD=" \
		${BRESEQ} \
		-b 0 \
        -o ${SELF} \
        -r ${DATADIR}/lambda/lambda.5.gbk
        ${DATADIR}/lambda/lambda.short_sequence_repeats.fastq
	"

do_test $1 ${SELF}
