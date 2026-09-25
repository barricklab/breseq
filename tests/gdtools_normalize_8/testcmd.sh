#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# --con-minimum-mutations raises the bar for a gene conversion, and --no-gene-conversion is the
# switch that turns the step off.
#
# The input is gdtools_normalize_7's: SNP 10028 C (explained by no copy of the repeat) and SNP
# 10210 G (the donor copy's base). Requiring two explained mutations per CON means nothing
# qualifies, so both SNPs must come through unchanged.

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        NORMALIZE \
        --con-minimum-mutations 2 \
        -o ${SELF}/output.gd \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/../gdtools_normalize_7/input.gd \
    "

do_test $1 ${SELF}
