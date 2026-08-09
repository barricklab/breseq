#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
# Large real-world data, downloaded on demand and md5-verified rather than
# committed. See tests/long_ltee_clone/testcmd.sh for what TEST_DATA does and
# why it has to be declared here, before common.sh is sourced.
TEST_DATA=ena_SRR098033,ltee_REL606
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# LTEE clone ZDB30: population Ara-3, generation 32,000 (Blount et al. 2012).
# Paired-end Illumina Genome Analyzer II, 2x35 bp, 25,400,140 read pairs.
#
# WHY THIS DATA SET: insert size, and nothing else comes close. ENA records
# nominal_length 2834 for this library, so with 2x35 reads each pair spans about
# 2.76 kb of sequence that was never read. Every other paired data set in the
# LTEE collection -- and every other long test -- is a few hundred bases at
# most, so this is the only place the pair-distance and discordant-pair machinery
# sees a kilobase-scale span. Alongside long_ltee_ara_m1_40k_pe36 (135 bp, also
# recorded) it brackets the range: the two libraries whose insert sizes we
# actually know, twenty-fold apart.
#
# It is deliberately a second Ara-3 clone (long_ltee_ara_m3_38k_se36 is ZDB107 at
# 38,000 generations). Duplicating the population is the price of the insert
# size; the two clones' curated genome diffs are unalike anyway -- this one is
# MOB-heavy with no AMP, that one is AMP-heavy.
#
# COVERAGE CAP: as served this run is ~384x, far above the 40-80x band the rest
# of these tests were chosen for, and at 2x35 bp that is 50 million reads. '-l 80'
# caps what breseq PARSES at roughly 80-fold. The limit is applied inside
# normalize_fastq_paired (see read_file_base_limit in breseq_cmdline.cpp), so it
# truncates within the paired set -- it takes the head of the files in order and
# never splits a pair, which is what makes the golden reproducible. The full
# 1.39 GB is still downloaded and md5-verified: the manifest checksum has to
# cover exactly the bytes ENA serves, which is why the data is not subsampled on
# disk instead. This is also the largest download of any test in the tree.
#
# EXPERIMENTAL PAIR EVIDENCE: this test runs the three PAIR predictors, and it is
# the only place any of them see a kilobase-scale span. What they currently produce
# here is MP 98, PD 6 and DP 0 -- the zero is expected and worth pinning: the
# distance distribution is unimodal about 2837 with a discordance cutoff of ~5381,
# so essentially no pair is an outlier INDIVIDUALLY, which is precisely the regime
# PD exists for and DP cannot see. If a change ever makes DP fire in bulk here, it
# is calling library structure, not variants. The golden will move whenever the
# DP/MP/PD code changes; that is the point, so rebuild it and review the diff
# rather than treating the failure as a surprise.
#
# SOFT CLIPPING IS DELIBERATELY OFF HERE, unlike the other three paired long tests.
# On a mate-pair library SC does not measure the genome, it measures the library
# prep, and it buries everything else:
#
#   * 15,086 SC entries against 21-46 on the other paired tests, 4,351 of them
#     accepted, spread evenly at roughly one per 1,064 bp of the 4.6 Mb genome, and
#     97% with consensus_fraction 1.0000. The golden was 4.1 MB, about 90% SC.
#   * They are real, and they are circularization junctions. Matching each clipped
#     tail back to the reference puts 2,610 of them 2.4-3.3 kb from their own SC
#     position -- against 2 expected by chance, a 1305x enrichment -- and 91% land
#     within 2.4-5 kb, centred on this library's 2837 bp fragment size. A soft clip
#     here is the far end of the same fragment, joined during circularization.
#
# Trimming would not have helped: the 2011 Illumina mate-pair v2 junction is a
# direct genomic-to-genomic blunt ligation with no adapter in it. Scanning 400,000
# raw reads finds zero Nextera Mate Pair junction adapter and 0.01% TruSeq
# read-through, so there is nothing for fastp to cut, and NxTrim-style splitters key
# on exactly the junction adapter this protocol lacks. SC is behaving correctly; the
# artifact is simply indistinguishable from structural variation without pairing
# information, which is why it is the pair predictors that earn their keep here.
#
# Consequence to be aware of when comparing goldens: --predict-soft-clipping lowers
# require-match-fraction from 0.9 to 0.5 unless set explicitly (settings.cpp). This
# test does NOT pass it, so it maps at the stock 0.9 while the other three paired
# long tests map at 0.5.
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
    -l 80 \
    --predict-discordant-pairs \
    --predict-missing-pairs \
    --predict-pair-distance \
    ${DOWNLOADDIR}/ena_SRR098033/SRR098033_1.fastq.gz \
    ${DOWNLOADDIR}/ena_SRR098033/SRR098033_2.fastq.gz \
    "

do_test $1 ${SELF}
