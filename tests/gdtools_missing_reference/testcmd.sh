#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Failure test: a gdtools subcommand must exit non-zero when a -r reference does not
# exist, and must say so before it starts work.
#
# This pins the gdtools half of the path-role rollout. gdtools had essentially no
# existence checks of its own -- a missing file surfaced as cGenomeDiff::read's
# "Could not open file for reading" assertion (or, for a reference, deep inside
# LoadFiles) rather than as an option error -- so without this test the markers on
# the ~25 gdtools subcommands are unexercised.
EXPECTED_EXIT_CODE=1

TESTCMD="\
    ${GDTOOLS} \
        ANNOTATE \
        -o ${SELF}/output.gd \
        -r ${SELF}/does_not_exist.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
