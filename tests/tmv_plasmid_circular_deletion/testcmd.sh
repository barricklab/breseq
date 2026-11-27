#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/tmv_plasmid/tmv-plasmid.gbk"

# A deletion across the origin of a circular plasmid, on real paired reads.
#
# This is also the NEGATIVE control for missing pair (MP) evidence, which used to be a separate
# test, tmv_plasmid_missing_pairs: the identical run with --predict-missing-pairs in place of
# --predict-discordant-pairs. Both flags are deprecated and inert now that DP and MP are predicted
# by default, so the two tests had become one run with one golden. What that test asserted is
# asserted here: the singleton marking in resolve_alignments, the MP candidate-region seeding in
# identify_mutations and the MP prediction stage all run on this data, and MP does NOT fire,
# because both of this dataset's breakpoints are reference-mapped and its unmappable mates are
# just a flat background. expected.gd has no MP line; one appearing fails the test.
# tests/lambda_pair_distances is the positive control.
TESTCMD=" \
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    --predict-discordant-pairs \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${DATADIR}/tmv_plasmid/D3-9_1P.fastq.gz \
    ${DATADIR}/tmv_plasmid/D3-9_2P.fastq.gz \
	"

do_test $1 ${SELF}
