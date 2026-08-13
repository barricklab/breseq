#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
# The second half of this test reads data/reference.bam produced by lambda_mult_ref_read, which is
# the only test whose BAM has more than one @RG. The Snakefile parses TEST_DEPENDS to enforce that
# ordering when running in parallel; common.sh ignores it (serial runners rely on sorted discovery
# order, where lambda_mult_ref_read already comes first).
TEST_DEPENDS=lambda_mult_ref_read
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.no_read_groups.tab"
EXPECTED_OUTPUTS[0]="${SELF}/expected.no_read_groups.tab"
CURRENT_OUTPUTS[1]="${SELF}/output.multiple_read_groups.tab"
EXPECTED_OUTPUTS[1]="${SELF}/expected.multiple_read_groups.tab"

# Covers both ends of --per-read-group:
#
# 1. bull_1.bam has no @RG lines at all (it predates read groups and is not written by breseq), so
#    every read falls back to group 0 and exactly one RG-0 column set is emitted. Its values must
#    equal the unprefixed totals, since there is nothing else for the reads to belong to.
#
# 2. lambda_mult_ref_read runs breseq over 7 read files, so its BAM carries several @RG lines --
#    including two for read files that contributed no aligned reads, whose columns are therefore
#    all zero. Per position the per-group columns must sum to the unprefixed totals. (That test
#    also uses several reference files, hence the -0 suffix on the seq id.)
TESTCMD="\
    ${BRESEQ} BAM2COV \
        --format TSV \
        --per-read-group \
        -f ${DATADIR}/bull/bull_1.fasta \
        -b ${DATADIR}/bull/bull_1.bam \
        -o ${CURRENT_OUTPUTS[0]} \
        -r rachael:657-767 \
    && ${BRESEQ} BAM2COV \
        --format TSV \
        --per-read-group \
        -f ${SELF}/../lambda_mult_ref_read/data/reference.fasta \
        -b ${SELF}/../lambda_mult_ref_read/data/reference.bam \
        -o ${CURRENT_OUTPUTS[1]} \
        -r NC_001416-0:8000-8110 \
    "

do_test $1 ${SELF}
