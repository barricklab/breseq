#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
# Large real-world data, downloaded on demand and md5-verified rather than
# committed. See tests/long_ltee_clone/testcmd.sh for what TEST_DATA does and
# why it has to be declared here, before common.sh is sourced.
TEST_DATA=ena_SRR030258,ltee_REL606
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# LTEE clone REL10938: population Ara-1, generation 40,000, mutT hypermutator
# (Barrick et al. 2009). Paired-end Illumina Genome Analyzer, 2x36 bp,
# 4,205,734 read pairs, ~65x.
#
# WHY THIS DATA SET: it is the oldest era in the collection -- 2009-vintage
# Genome Analyzer reads at a length no other long test uses -- and the only
# short-insert paired library whose insert size is actually recorded (ENA
# nominal_length 135). With 2x36 reads that means the mates never overlap and
# every pair leaves a ~63 bp unsequenced gap, which is the regime the pair-based
# evidence types are least exercised in elsewhere. Its curated genome diff has
# 656 mutations including 10 MOB and 1 AMP, so it is a real test of junction and
# copy-number prediction and not just of SNP calling.
#
# Observed median pair distance from this run: 137 (MAD 11), against ENA's recorded
# nominal_length of 135 -- so the metadata is trustworthy here, and this really is
# the ~135 bp end of the insert-size range. Compare long_ltee_ara_m3_32k_mp2800 at the
# other end.
#
# EXPERIMENTAL PAIR EVIDENCE: this test deliberately turns on the four
# experimental predictors so that DP/MP/PD/SC are covered by a real data set and
# not only by the small lambda fixtures. Two consequences:
#   * --predict-soft-clipping lowers require-match-fraction from 0.9 to 0.5
#     unless it is given explicitly (settings.cpp), which changes alignment
#     acceptance globally -- so this golden is NOT comparable to
#     long_ltee_clone's, which runs stock options.
#   * this golden will move whenever the DP/MP/PD/SC code changes. That is the
#     point of the test; rebuild it and review the diff rather than treating the
#     failure as a surprise.
#
# See tests/long_ltee_clone/testcmd.sh for the notes that apply to every long
# test: rebuild with BRESEQ_TEST_DATA_DIR unset so the committed expected.gd
# keeps repo-relative paths, and note that #=MAPPED-BASES / #=MAPPED-READS are
# compared, so a bowtie2 version change fails this test even with breseq
# unchanged.
REFERENCE_ARG="-r ${DOWNLOADDIR}/ltee_REL606/REL606.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    --predict-copy-number \
    --predict-discordant-pairs \
    --predict-missing-pairs \
    --predict-pair-distance \
    --predict-soft-clipping \
    ${DOWNLOADDIR}/ena_SRR030258/SRR030258_1.fastq.gz \
    ${DOWNLOADDIR}/ena_SRR030258/SRR030258_2.fastq.gz \
    "

do_test $1 ${SELF}
