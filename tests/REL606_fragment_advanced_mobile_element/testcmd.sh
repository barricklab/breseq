#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output/evidence/annotated.gd"
EXPECTED_OUTPUTS[0]="expected.gd"
REFERENCE_ARG="-r ${DATADIR}/REL606/REL606.fragment.gbk"

TESTCMD=" \
		${BRESEQ} \
		-b 0 \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        ${DATADIR}/REL606/REL606.advanced_mobile_element.fastq
	"

do_test $1 ${SELF}
