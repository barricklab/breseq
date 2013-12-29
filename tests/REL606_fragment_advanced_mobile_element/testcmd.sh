#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

OUTPUT_CHECKS[0]="output/evidence/annotated.gd expected.gd"

TESTCMD=" \
		${BRESEQ} \
		-b 0 \
        -o ${SELF} \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${DATADIR}/REL606/REL606.advanced_mobile_element.fastq
	"

do_test $1 ${SELF}
