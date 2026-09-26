#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# The block map from the coordinates of the mutated sequence back to those of the reference,
# on the gdtools_apply_3 input (DEL, INV, MOB, SNP, INS, AMP, SUB, CON, chained INS).
CURRENT_OUTPUTS[0]="${SELF}/output.tsv"
EXPECTED_OUTPUTS[0]="${SELF}/expected.tsv"

TESTCMD="\
    ${GDTOOLS} \
        APPLY \
        -f FASTA \
        -o ${SELF}/output.fasta \
        --coordinate-map ${SELF}/output.tsv \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/../gdtools_apply_3/input.gd \
    "

do_test $1 ${SELF}
