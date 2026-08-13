#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=7
# Large real-world data, downloaded on demand and md5-verified rather than
# committed (see tests/test_data_manifest.tsv and tests/fetch_test_data.sh).
# Declared here, next to TEST_CORES and *before* sourcing common.sh, so that
# tests/Snakefile can parse it with the same anchored regex it uses for
# TEST_CORES/TEST_DEPENDS and schedule the downloads as real job-graph nodes
# (deduplicated across tests). common.sh's do_test acts on it too, so
# './tests/test.sh test long_ltee_clone' -- which bypasses Snakemake entirely --
# fetches the same data.
TEST_DATA=ena_SRR2589061,ltee_REL606
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# A full, real Escherichia coli genome: Long-Term Evolution Experiment clone
# REL4540A (population Ara-5, generation 10,000, non-mutator; Tenaillon et al.
# 2016), 2,384,764 paired-end Illumina read pairs at roughly 100x against the
# REL606 ancestor. This is the case none of the small, fast tests cover -- a
# complete 4.6 Mb genome with real error profiles, real IS-element activity and
# real coverage structure -- which is exactly why it is a 'long' test and not part
# of 'make test'. The same clone is the worked example in
# docs/tutorial-curation-ecoli-ltee.md.
#
# The reference is the Barrick lab's curated REL606.gbk (pinned to a commit of
# barricklab/LTEE), NOT the NCBI RefSeq assembly this test originally used. The
# RefSeq GenBank annotates no IS elements at all, and breseq predicts MOB from
# annotated mobile elements, so against RefSeq this clone's IS-element activity
# came out as unassigned JC evidence and zero MOB entries. See the comment on the
# ltee_REL606 row in tests/test_data_manifest.tsv.
#
# ${DOWNLOADDIR} is ./tests/data/downloads unless BRESEQ_TEST_DATA_DIR points
# elsewhere. Either is fine for the comparison -- do_check filters out
# ^(#=COMMAND|#=CREATED|#=PROGRAM|#=READSEQ|#=REFSEQ) (DIFF_IGNORE in common.sh),
# which are the only lines a .gd records input paths on. But REBUILD THIS TEST'S
# expected.gd WITH BRESEQ_TEST_DATA_DIR UNSET, so the committed file keeps
# repo-relative paths instead of someone's absolute home directory.
#
# Note that #=MAPPED-BASES / #=MAPPED-READS *are* compared, so a bowtie2 version
# change will fail this test even with breseq itself unchanged. That is
# deliberate for a test whose job is to detect drift on real data;
# dev-environment.yml pins bowtie2 for exactly this reason.
REFERENCE_ARG="-r ${DOWNLOADDIR}/ltee_REL606/REL606.gbk"

TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    ${REFERENCE_ARG} \
    --predict-copy-number \
    ${DOWNLOADDIR}/ena_SRR2589061/SRR2589061_1.fastq.gz \
    ${DOWNLOADDIR}/ena_SRR2589061/SRR2589061_2.fastq.gz \
    "

do_test $1 ${SELF}
