#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

OUTPUT_CHECKS[0]="output.gff3 expected.gff3"

TESTCMD="\
    ${GDTOOLS} \
        APPLY
        -f GFF3
        -o ${SELF}/output.gff3 \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
