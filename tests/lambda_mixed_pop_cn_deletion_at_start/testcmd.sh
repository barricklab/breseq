#!/bin/bash

# Copy number analysis on a reference whose FIRST window is called CN != 1.
#
# CNery's break_pts.csv opened the first segment of every sequence at Startpos 0 (its
# _segments_from_path() seeded start_pos = 0 before it had looked at a window). breseq read that
# through into a CN evidence entry, and a CN with start = 0 fails GenomeDiff validation -- so the
# run died at the Output stage, several stages after the file that actually caused it:
#
#   >>ERROR: Expected positive integral value for field 9: [start] instead of [0].
#
# It stayed hidden because copy number 1 is the baseline breseq drops, so while a genome opened in
# state 1 the bad coordinate was filtered out before it could reach a .gd. What makes row 1 survive
# is a first window called CN != 1 -- reproduced here with a reference that opens INSIDE a real
# deletion (lambda 21701-48502, where 21701-27700 is deleted in this read set, so the first 6000
# bases have no coverage at all). The other way in is a library too thin to call anything, where
# CNery normalises against 1.0 and returns CN 0 genome-wide.
#
# Fixed in CNery, by its commit "Open the first segment of a sequence at 1, not 0". breseq does not
# work around the old behaviour -- read_cnery_segments() rejects a start below 1 outright -- so this
# test REQUIRES a CNery at or past that commit. dev-environment.yml pins cnery-prerelease to an
# exact commit build, so that is a pin question, not a user question: if this test fails with
# "Segment starting before the first base of the sequence (0) on row 1", the pin is behind the fix.
#
# What it checks is the outcome rather than the mechanism: a CN entry starting at 1, picked up as
# evidence by the deletion that explains it.

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.deletion_at_start.fna"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    --predict-copy-number \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    "

do_test $1 ${SELF}
