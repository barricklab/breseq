#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.5.gbk"


TESTCMD=" \
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -b 0 \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda.short_sequence_repeats.fastq
	"

do_test $1 ${SELF}
