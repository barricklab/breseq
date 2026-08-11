#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
# Large real-world data, downloaded on demand and md5-verified rather than
# committed. See tests/long_ltee_clone/testcmd.sh for what TEST_DATA does and
# why it has to be declared here, before common.sh is sourced.
TEST_DATA=ena_SRR1536190,ltee_REL606
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# LTEE clone REL10985: population Ara+6, generation 40,000, hypermutator
# (Raeside et al. 2014). Single-end Illumina Genome Analyzer IIx, 36 bp,
# 5,541,345 reads, ~43x.
#
# WHY THIS DATA SET: it pairs the lowest depth in the suite with the heaviest
# mutation load. 43x is the shallow end of what these tests cover -- deep enough
# to call confidently, shallow enough that anything which quietly raises an
# evidence threshold shows up here first -- and the curated genome diff has
# 2,205 SNPs, so the read-alignment evidence path is being asked to make
# thousands of real calls rather than a couple of dozen. It is also a third
# distinct sequencing submission (2014 rearrangement study, Genome Analyzer IIx)
# and, with long_ltee_ara_m3_38k_se36, the second unpaired data set.
#
# Stock breseq options: the three pair predictors require paired mapping, so
# there is nothing to turn on here.
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
    ${DOWNLOADDIR}/ena_SRR1536190/SRR1536190.fastq.gz \
    "

do_test $1 ${SELF}
