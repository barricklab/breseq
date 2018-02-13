#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gff3"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gff3"

TESTCMD="\
    ${BRESEQ} \
        CONVERT-REFERENCE \
        -f GFF \
        -o ${SELF}/output.gff3 \
        ${DATADIR}/pDCAF3/pDCAF3.gbk \
    "

do_test $1 ${SELF}
