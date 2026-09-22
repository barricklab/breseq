#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# A DEL that crosses the origin of a circular sequence, after other mutations have
# already shifted its position. GFF3 output has the features and the sequence.

CURRENT_OUTPUTS[0]="${SELF}/output.gff3"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gff3"

TESTCMD="\
    ${GDTOOLS} \
        APPLY \
        -f GFF3 \
        -o ${SELF}/output.gff3 \
        -r ${DATADIR}/tmv_plasmid/tmv-plasmid.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
