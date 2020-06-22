#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gbk"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gbk"

TESTCMD="\
    ${GDTOOLS} \
        APPLY \
        -f GENBANK \
        -o ${SELF}/output.gbk \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/input.gd \
        && sed -i '' '1s/[[:digit:]]\{2\}-.\{3\}-[[:digit:]]\{4\}/XX-XX-XXXX/' ${CURRENT_OUTPUTS[0]} \
    "

do_test $1 ${SELF}
