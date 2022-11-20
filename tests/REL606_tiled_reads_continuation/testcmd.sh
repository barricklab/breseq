#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/REL606/REL606.fragment.gbk"


TESTCMD=" \
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -p \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    -s ${DATADIR}/REL606/REL606.is.gbk \
    ${DATADIR}/REL606/REL606.tiled_reads_continuation.fastq
	"

do_test $1 ${SELF}
