#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Same colliding insertions as gdtools_normalize_3, but written the way a hand-authored file
# usually is: with insert_position left off.
#
# Parsing fills in insert_position=1 and flags '_dont_print_insert_position' so the field is
# dropped again on write. That flag records that the input line omitted the field, not that the
# field is unwanted forever -- so it used to discard the value the INS back-off had just
# renumbered, emitting two entries identical except for their id and losing an insertion.
#
# Expected: the insertion still holding the implicit 1 prints without the field (unchanged), the
# renumbered one prints insert_position=2, and NORMALIZE exits 0.

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
