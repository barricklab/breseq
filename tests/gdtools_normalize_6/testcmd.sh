#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# A gene conversion whose donor is on the other strand of another sequence, with indels.
#
# gene_conversion.fna has two contigs. 'recip' carries copy A of a 200-bp sequence between
# unrelated flanks; 'donor' carries the REVERSE COMPLEMENT of copy B, which differs from A at six
# places: two near its ends that stay unconverted (so the tract is bounded inside the repeat), and
# four inner ones that input.gd applies to recip: a SNP, a 2-bp INS, a 1-bp DEL, and a SNP. After
# those four the recipient matches the donor without gaps from one boundary difference to the
# other, so NORMALIZE must write one CON covering everything between those two differences (179
# reference bases, recip 161..339), whose region runs backwards (start > end = reverse complement)
# and is 180 bp long, one longer than the tract because the cluster nets +1 base. The APPLY check
# confirms the bookkeeping.

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        NORMALIZE \
        -o ${SELF}/output.gd \
        -r ${DATADIR}/gene_conversion/gene_conversion.fna \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
