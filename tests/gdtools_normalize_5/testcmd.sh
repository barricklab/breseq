#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# NORMALIZE replaces a cluster of mutations with the gene conversion (CON) that explains it.
#
# REL606.fragment.gbk carries two forward-strand copies of a 768-bp repeat (9767..10534 and
# 31302..32069) that differ at exactly two positions: 10028 T / 31563 A and 10210 T / 31745 G.
# input.gd has the two SNPs that make the first copy identical to the second. Together they are
# one recombination event, and after it the recipient matches the donor over the whole repeat, so
# NORMALIZE should write one CON with the maximal extent -- the whole identical stretch, which is
# exactly the annotated repeat, not just the span of the two SNPs: CON 9767 768 REL606-5:31302-32069.
# The APPLY check at the end of NORMALIZE confirms the CON reproduces the same genome as the SNPs.

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        NORMALIZE \
        -o ${SELF}/output.gd \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
