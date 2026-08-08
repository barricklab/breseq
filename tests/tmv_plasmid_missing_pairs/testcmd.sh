#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/tmv_plasmid/tmv-plasmid.gbk"

# Exercises --predict-missing-pairs end to end: the singleton marking in resolve_alignments, the MP
# candidate-region seeding in identify_mutations, and the MP prediction stage. This dataset's
# breakpoints are both reference-mapped, so it is a smoke test rather than a positive control -- a
# real MP call needs paired reads over genuinely novel inserted sequence.
TESTCMD=" \
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    --predict-missing-pairs \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    ${DATADIR}/tmv_plasmid/D3-9_1P.fastq.gz \
    ${DATADIR}/tmv_plasmid/D3-9_2P.fastq.gz \
	"

do_test $1 ${SELF}
