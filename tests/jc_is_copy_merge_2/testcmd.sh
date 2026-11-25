#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=3
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# Two candidate junctions that resolve to ONE JC entry, one PASSING and one REJECTED.
#
# The same data as jc_is_copy_merge_1 (see there and tests/data/annamia_is630/README), with
# --junction-score-cutoff lowered to 0.25 so that the weaker candidate (pos_hash_score 4,
# neg_log10_pos_hash_p_value 0.4) is rejected as COVERAGE_EVENNESS_SKEW while the stronger one
# (12, 0.1) passes. This is the arrangement that stopped a full-genome run with a fatal duplicate
# entry in jc_evidence.gd.
#
# The passing candidate is kept and the rejected one is dropped without adding its reads: they were
# written back to the reference BAM when it was rejected, so the survivor has no merged_keys.
TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    --junction-score-cutoff 0.25 \
    -o ${SELF} \
    -r ${DATADIR}/annamia_is630/AP025630_is630_windows.gbk.gz \
    ${DATADIR}/annamia_is630/DRR285493_subset_1.fastq.gz \
    ${DATADIR}/annamia_is630/DRR285493_subset_2.fastq.gz \
    "

do_test $1 ${SELF}
