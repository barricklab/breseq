#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Failure test: an INS insert_position below 1 must be rejected, while RA keeps 0.
#
# insert_position numbers the pileup columns at a reference position: 0 is the reference base
# itself and 1 is the first inserted base after it. An INS is the inserted bases, so it has no
# column 0 -- which is also why a blank insert_position is filled in with 1 on load. RA does use
# 0, for a base that is not in an insertion, so the floor cannot be set on the shared field type;
# diff_entry_load_field_variable_types narrows it for INS alone.
#
# The field was typed as a plain Integer, which set no floor at all, and a negative did reach
# real output: APPLY negates insert_position to carry chained insertions through an INV in
# reverse order, and VALIDATE accepted the result.
#
# The third line is the control: an RA at column 0 must still pass.
EXPECTED_EXIT_CODE=1

TESTCMD="\
    ${GDTOOLS} \
        VALIDATE \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
