#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
# Large real-world data, downloaded on demand and md5-verified rather than
# committed. See tests/long_ltee_clone/testcmd.sh for what TEST_DATA does and
# why it has to be declared here, before common.sh is sourced.
TEST_DATA=ena_SRR098038,ltee_REL606
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# LTEE clone ZDB107: population Ara-3, generation 38,000, Cit+ hypermutator
# (Blount et al. 2012). Single-end Illumina Genome Analyzer II, 36 bp,
# 7,046,631 reads, ~55x.
#
# WHY THIS DATA SET: two things no other long test has. First, it is UNPAIRED --
# every other real-data long test is paired, so without this one the entire
# single-end path through a full-size genome is untested. Second, it is the Cit+
# lineage: its curated genome diff carries 3 AMP entries (the rnk-citT
# amplification), 29 MOB and 33 DEL, which is the densest copy-number and mobile
# element content available anywhere in the LTEE collection. If CN, AMP or
# junction prediction regresses, this is the test that should notice.
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
    ${DOWNLOADDIR}/ena_SRR098038/SRR098038.fastq.gz \
    "

do_test $1 ${SELF}
