#!/bin/bash
#
# Copy number (CN) evidence in POLYMORPHISM MODE, which is a different call than the consensus one
# lambda_mixed_pop_cn_evidence covers: breseq hands CNery --polymorphism-mode whenever it is itself
# run with -p, and CNery then decodes over a continuous grid and writes fractional copy numbers.
#
# The .gd alone cannot prove any of that happened. On this read set the continuous grid lands on the
# same segmentation the integer one does, so this test's expected.gd is -- correctly -- the same
# shape as lambda_polymorphism's, with copy_number=0 on the one CN entry. That is worth asserting
# (it is where a fractional 0.0 collapsing back to "0" is checked end to end), but on its own it
# would pass just as well if the -p never reached CNery at all. check_polymorphism_mode.sh is what
# closes that gap, by looking at the dtype CNery wrote into break_pts.csv.
#
# -k because that file lives in 09_copy_number_variation/, which is deleted when the pipeline
# finishes. --predict-copy-number is deprecated and CN is on by default, but it is what sets
# copy_number_explicitly_requested, which makes a MISSING CNery fatal instead of a warning-and-skip
# -- a CN test that passes quietly on a machine without CNery is worth nothing, which is why the
# three consensus CN tests all still carry the flag too.
#
# No --genbank-field-for-seq-id, deliberately: that keeps the seq id NC_001416 and the coordinates
# directly comparable with tests/lambda_polymorphism/expected.gd, which is the file this one should
# be diffed against whenever it is rebuilt.

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"

TESTCMD="\
    ${BRESEQ} \
    -k \
    ${BRESEQ_TEST_THREAD_ARG} \
    -p \
    -o ${SELF} \
    -g ${SELF}/header.gd \
    --predict-copy-number \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    && ${SELF}/check_polymorphism_mode.sh ${SELF} \
    "

do_test $1 ${SELF}
