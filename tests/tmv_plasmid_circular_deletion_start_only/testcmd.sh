#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/tmv_plasmid/tmv-plasmid-truncate-end.gbk"

TESTCMD=" \
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${DATADIR}/tmv_plasmid/D3-9_1P.fastq.gz \
    ${DATADIR}/tmv_plasmid/D3-9_2P.fastq.gz \
	"

do_test $1 ${SELF}
