#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
# Large real-world data, downloaded on demand and md5-verified rather than
# committed. See tests/long_ltee_clone/testcmd.sh for what TEST_DATA does and
# why it has to be declared here, before common.sh is sourced.
TEST_DATA=ena_SRR2588848,ltee_REL606
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# LTEE clone REL10411: population Ara+3, generation 30,000, mismatch-repair
# deficient hypermutator (Tenaillon et al. 2016). Paired-end Illumina HiSeq
# 2500, 2x150 bp, 1,809,983 read pairs, ~82x.
#
# WHY THIS DATA SET: the longest reads in the suite, and a mutational spectrum
# that is mostly SHORT INDELS instead of substitutions. Losing mismatch repair
# means indels accumulate in homopolymer runs, and the curated genome diff shows
# it: 126 INS and 59 DEL alongside the SNPs. That is the regime where indel
# calling and the homopolymer error model are actually under load, and no other
# long test puts them there.
#
# Insert size is not recorded in ENA for the 2016 submissions (nominal_length is
# blank, and the NCBI SRA XML has no NOMINAL_LENGTH either), so it had to be
# measured: this run reports a median pair distance of 573
# (data/summary.json, pair_distance_median). That is a genuinely different regime
# from long_ltee_ara_m1_40k_pe36 at 137 and long_ltee_ara_m3_32k_mp2800 at ~2.8 kb, so
# the suite spans roughly 137 / 573 / 2800 bases of insert.
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
    ${DOWNLOADDIR}/ena_SRR2588848/SRR2588848_1.fastq.gz \
    ${DOWNLOADDIR}/ena_SRR2588848/SRR2588848_2.fastq.gz \
    "

do_test $1 ${SELF}
