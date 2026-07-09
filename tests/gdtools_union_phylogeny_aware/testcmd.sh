#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# Two input GD files contain a mutation that is identical in every field except
# for its phylogeny_id= tag. With --phylogeny-aware (-p), UNION must keep both
# entries and preserve their phylogeny_id tags. This is a regression test for a
# bug where merge() stripped phylogeny_id in phylogeny-aware mode, collapsing the
# two mutations into one and dropping the tag from the output.
TESTCMD="\
    ${GDTOOLS} \
        UNION \
        -p \
        -o ${SELF}/output.gd \
        ${SELF}/input_1.gd \
        ${SELF}/input_2.gd \
    "

do_test $1 ${SELF}
