#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=7
# Large real-world data, downloaded on demand and md5-verified rather than
# committed. See tests/long_ltee_clone/testcmd.sh for what TEST_DATA does and
# why it has to be declared here, before common.sh is sourced.
TEST_DATA=ena_SRR2584534,ltee_REL606
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# LTEE clone REL11393: population Ara+1, generation 50,000, non-mutator
# (Tenaillon et al. 2016). Paired-end Illumina HiSeq 2000, 2x101 bp, 1,720,420
# read pairs, ~75x.
#
# WHY THIS DATA SET: structural variation. It is the last generation of the
# experiment and, because the population never became a mutator, 50,000
# generations of selection produced a mutation set that is structural rather
# than substitutional -- the curated genome diff has 42 MOB, 26 DEL, 3 AMP and
# the only INV of any clone considered for these tests. That combination
# exercises junction prediction, IS-element handling and copy-number calling
# together, which no other data set here does.
#
# It shares a platform family with long_ltee_clone (2x101 HiSeq 2000) but
# nothing else: different population, five times the generations, an unrelated
# mutation profile, and the experimental pair predictors turned on.
#
# Insert size is not recorded in ENA for the 2016 submissions, so it was measured:
# this run reports a median pair distance of 232 (data/summary.json,
# pair_distance_median). The 2x101 and 2x150 libraries are NOT the same prep -- 232
# against 573 for long_ltee_ara_p3_30k_pe150 -- so the two tests differ in insert
# size as well as read length, which is why both are worth keeping.
#
# EXPERIMENTAL PAIR EVIDENCE: like the other paired long tests, this one turns on
# the four experimental predictors. Note that --predict-soft-clipping lowers
# require-match-fraction from 0.9 to 0.5 unless set explicitly (settings.cpp),
# so this golden is NOT comparable to long_ltee_clone's stock-options output,
# and it will move whenever the DP/MP/PD/SC code changes -- rebuild and review
# the diff rather than treating that as a surprise.
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
    ${DOWNLOADDIR}/ena_SRR2584534/SRR2584534_1.fastq.gz \
    ${DOWNLOADDIR}/ena_SRR2584534/SRR2584534_2.fastq.gz \
    "

do_test $1 ${SELF}
