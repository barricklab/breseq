#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="output.gff3"
EXPECTED_OUTPUTS[0]="expected.gff3"

TESTCMD="\
    ${GDTOOLS} \
        APPLY
        -f GFF3
        -o ${SELF}/output.gff3 \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        -r ${DATADIR}/REL606/REL606.fragment.2.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
