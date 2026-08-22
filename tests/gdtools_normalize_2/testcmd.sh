#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Failure test: NORMALIZE must reject an ambiguous Genome Diff, as APPLY already does.
#
# input.gd has a plain SNP sitting on bases the DEL removes, with no 'within' or 'before' to say
# which applies first. The two orderings give different genomes, so there is no single coordinate
# to normalize to. NORMALIZE used to accept this and shift the DEL anyway, which produced a
# silently different deletion; it now fails with the same error APPLY gives.
EXPECTED_EXIT_CODE=1

TESTCMD="\
    ${GDTOOLS} \
        NORMALIZE \
        -o ${SELF}/output.gd \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
