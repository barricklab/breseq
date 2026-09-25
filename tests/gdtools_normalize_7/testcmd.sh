#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# A cluster only partly explained by a donor: the explained part becomes a CON, the rest stays.
#
# Same repeat pair as gdtools_normalize_5 (9767..10534 vs 31302..32069, differing at 10028 and
# 10210). Here the SNP at 10028 is to C, which neither copy carries, while 10210 G is the donor's
# base. The identical stretch around 10210 stops at 10028 on the left and at the end of the repeat
# on the right, so only 10210 is explained. With the default of one mutation per CON that lone SNP
# becomes a CON whose tract is that whole stretch (CON 10029 506 REL606-5:31564-32069), and the SNP
# at 10028 is left alone just outside it. gdtools_normalize_8 runs the same input with
# --con-minimum-mutations 2 and expects no CON.

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
