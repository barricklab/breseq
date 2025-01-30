#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gff3"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gff3"

CURRENT_OUTPUTS[1]="${SELF}/output.gd"
EXPECTED_OUTPUTS[1]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        APPLY \
        -f GFF3 \
        -o ${SELF}/output.gff3 \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        --applied-gd ${SELF}/output.gd \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
