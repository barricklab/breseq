#!/bin/bash

# Copy number analysis on a run where one reference sequence gets no reads at all, and another is
# junction-only. Both used to break the CNery step, and neither shows up in the .gd comparison:
#
#   * A reference with no coverage fails its coverage-distribution fit, which sets the deletion
#     propagation cutoff to -1, which made identify_mutations_pileup::pileup_callback() return
#     before writing a single row. CNery was handed a coverage table holding nothing but its header.
#   * A junction-only reference is never pileup'd at all, so it has no coverage table -- but its
#     path was still passed to CNery, which raises FileNotFoundError on a missing input and took the
#     whole breseq run down with it.
#
# The checks below are what actually pin those two, so they are part of the test command rather than
# the expected .gd: the failure modes are a missing file and an empty one, not a moved call.

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
# lambda.gbk gets the reads; other.gbk (687 bp, an unrelated anemone GFP mRNA) gets none.
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk -r ${DATADIR}/lambda/other.gbk -s ${DATADIR}/REL606/REL606.fragment.gbk"

# -k so 08_mutation_identification/ survives for the checks below.
TESTCMD="\
    ${BRESEQ} \
    -k \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    --predict-copy-number \
    ${REFERENCE_ARG} \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    && ${SELF}/check_coverage_tables.sh ${SELF} \
    "

do_test $1 ${SELF}
