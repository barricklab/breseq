#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# The path checks run on EVERY run, not only under --dry-run. A missing -r must fail
# at option-parsing time: before the output directory is created, and long before
# stage 01 would have been the first thing to try opening the file.
#
# There is deliberately no EXPECTED_EXIT_CODE here. The '!' negation lives inside
# TESTCMD so that BOTH halves are asserted -- "breseq failed" AND "nothing was
# created". Bash parses '! a && b' as '(! a) && b', so the command as a whole
# succeeds only if both hold, and EXPECTED_EXIT_CODE stays at its default of 0.
TESTCMD="\
    ! ${BRESEQ} \
        -o ${SELF}/output \
        -r ${SELF}/does_not_exist.gbk \
        ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    && [ ! -e ${SELF}/output ] \
    "

do_test $1 ${SELF}
