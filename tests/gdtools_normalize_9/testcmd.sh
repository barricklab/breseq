#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# A de novo mutation inside a converted tract splits it into two CONs around the mutation.
#
# Same repeat pair as gdtools_normalize_5. The SNPs at 10028 A and 10210 G are both the donor
# copy's bases, but the SNP at 10100 is to a base neither copy carries. It breaks the identical
# stretch, so the two donor-explained SNPs cannot share one CON without the de novo SNP lying
# inside it -- which VALIDATE would reject, since 'within' is not allowed on a CON. Expected: a
# CON for each explained SNP, each extending from its side of the repeat up to the base next to
# the de novo SNP (9767..10099 and 10101..10534), with the SNP at 10100 kept between them.

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
