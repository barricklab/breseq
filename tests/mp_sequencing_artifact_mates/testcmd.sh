#!/bin/bash
#
# Missing pair (MP) evidence does not count reads whose never-aligned mate is a sequencing artifact
# (unaligned_mate_is_sequencing_artifact, resolve_alignments.cpp; --missing-pair-maximum-mate-
# base-fraction, --missing-pair-maximum-mate-run-fraction).
#
# On two-colour Illumina a read that runs into a long homopolymer or a strong hairpin fails into
# poly-G or a low-complexity tail and aligns nowhere, and the mates it leaves behind made MP evidence
# with no insertion present: Trueperella CP097247 at a 12-bp C run, Varibaculum AP041144 at a
# 20-bp-stem hairpin, both confirmed reference-consistent by assembly and nanopore reads.
#
# tests/data/mp_artifact (make_mp_artifact_data.py): a 20-kb reference; reads from a sample genome
# carrying a REAL novel 1,500-bp insertion after 14000; plus 60 pairs ending just left of 6000 whose
# mate 2 has failed into G-rich read-through. The same reads are run twice:
#   output_baseline/  gate OFF: an MP at ~6000 (the case is live) as well as the one at ~14000
#   data/             gate at its defaults: only the MP at ~14000; the artifact mates set aside
# check.sh states both outcomes. The golden is the gated run's.

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=2
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/mp_artifact/mp_artifact.fna"
READS="${DATADIR}/mp_artifact/mp_artifact_R1.fastq.gz ${DATADIR}/mp_artifact/mp_artifact_R2.fastq.gz"

TESTCMD=" \
    ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -o ${SELF}/output_baseline \
        ${REFERENCE_ARG} \
        --missing-pair-maximum-mate-base-fraction 0 \
        --missing-pair-maximum-mate-run-fraction 0 \
        ${READS} \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        ${READS} \
    && bash ${SELF}/check.sh ${SELF}/output_baseline/output/output.gd ${SELF}/data/annotated.gd ${SELF}/data/reference.bam \
    "

do_test $1 ${SELF}
