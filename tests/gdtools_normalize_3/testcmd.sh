#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Covers the INS back-off's insert_position renumbering walk, which had no test.
#
# Both insertions sit in the same 7-base A homopolymer (REL606-5 5567..5573), so both normalize
# to position 5573 and collide. The second one processed backs off onto the first and the walk
# renumbers it, which is the only path that decrements the reverse_iterator inside the loop --
# and the only one that can step past rbegin().
#
# Expected: the two insertions share position 5573 with insert_position 1 and 2, and NORMALIZE
# exits 0. Without distinct insert_positions the two entries are byte-identical, which the
# post-normalization validity check rejects as a duplicate.

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
