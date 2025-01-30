#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gbk"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gbk"

TESTCMD="\
    ${BRESEQ} \
        CONVERT-REFERENCE \
        -f GENBANK \
        -o ${SELF}/output.gbk \
        ${DATADIR}/pDCAF3/KC619530.gbk \
        && perl -i -pe 's/^(LOCUS.+)\d{2}-\w{3}-\d{4}$/\$1XX-XX-XXXX/g' ${CURRENT_OUTPUTS[0]} \
    "

do_test $1 ${SELF}


