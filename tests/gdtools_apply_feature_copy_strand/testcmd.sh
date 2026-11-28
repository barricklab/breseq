#!/bin/bash
#
# Features copied by APPLY keep their own strand unless the COPY is inverted.
#
# repeat_feature_1() used to decide whether to mirror a copied feature from the feature's own
# strand rather than from the orientation of the copy. An AMP copies its unit forward, so every
# minus-strand feature of the new unit came out mirrored end over end within the unit and flipped
# to plus (seen on a 175-kb IS6-bounded duplication in Thermincola CM166461, where it also left a
# spurious IS copy at the start of the unit). The same test sat under INT and under the nested
# features of a MOB.
#
# input.gd, on REL606.fragment.gbk:
#   AMP 4001 2300 x2  unit 4001-6300 holds mltD complement(4031..5389) and gloB
#                     complement(5461..6216). Both copies of both genes must stay minus-strand,
#                     at the same offsets in each unit.
#   INT forward donor 4001-5500: mltD keeps minus strand.
#   INT reverse donor 5500-4001: the copy is inverted, so mltD becomes plus strand, mirrored.
#
# check.sh states those outcomes from the GFF3, so a golden rebuild cannot quietly accept a mirror.

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gff3"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gff3"

TESTCMD="\
    ${GDTOOLS} \
        APPLY \
        -f GFF3 \
        -o ${SELF}/output.gff3 \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/input.gd \
    && bash ${SELF}/check.sh ${SELF}/output.gff3 \
    "

do_test $1 ${SELF}
