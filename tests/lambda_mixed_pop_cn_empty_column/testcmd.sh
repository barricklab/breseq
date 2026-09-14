#!/bin/bash
#
# breseq must survive an empty cell in CNery's CNV.csv, which is what pandas writes for a NaN.
#
# It did not. Every optional column was read with from_string<double>(), whose ASSERT(!s.empty())
# ends the run, so one missing number in one window of one diagnostic column threw away a whole
# pipeline with "FATAL ERROR / Attempt to convert empty string". An absent COLUMN was already
# handled; an empty CELL says exactly the same thing and now takes the same path.
#
# The expected .gd is the one the unmodified CNery produces, which is the real assertion here: the
# holes are filled from the same chain CNery itself falls back along (ori-ter-corrected coverage ->
# GC-corrected -> raw), so they must not move a call. Filling them with a neutral 0.0 would read as
# a deletion and show up here as a changed relative_coverage.

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

# cnery_stub/CNery is the real CNery plus one blanked cell per file; breseq finds it by PATH.
STUBDIR=`cd ${SELF}/cnery_stub && pwd`

TESTCMD="\
    PATH=\"${STUBDIR}:\${PATH}\" \
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    --predict-copy-number \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    "

do_test $1 ${SELF}
