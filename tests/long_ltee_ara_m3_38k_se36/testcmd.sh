#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=7
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
# This is also THE test for predicting a deletion between two unannotated paralogs
# from missing coverage alone: DEL 2032711 23293 between manB and cpsG. Being
# single-end there is no discordant pair, no read crosses the breakpoint in a way
# that can be split into a junction, and 36 bp reads that do cross tie between the
# two copies and are dropped as redundant -- so there is not one RA anywhere near
# the locus either. Coverage is the only evidence, and the crossover is recovered
# by reading which copy the surviving reads came from at each column where the two
# copies differ. If this DEL disappears, that whole path is dead; no other test
# exercises it on real data without a DP to fall back on.
#
# The position is exact to the limit of the data. The crossover falls in a 123 bp
# stretch (2032588-2032711) where the two copies happen to be identical, so any
# start in that interval yields the same sequence; 2032711 is the rightmost one,
# matching the right-shift convention normalize_to_sequence() applies to DEL. The
# SIZE, 23293, is exact.
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
    ${DOWNLOADDIR}/ena_SRR098038/SRR098038.fastq.gz \
    "

do_test $1 ${SELF}
