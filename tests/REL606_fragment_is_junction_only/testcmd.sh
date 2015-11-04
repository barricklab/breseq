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
      -s ${DATADIR}/REL606/REL606.is.gbk \
      ${DATADIR}/REL606/REL606.junction_only.fastq
	"

do_test $1 ${SELF}
