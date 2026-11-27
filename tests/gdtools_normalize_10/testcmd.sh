#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# The gene conversion search must survive a seed that begins with 16 T's.
#
# The k-mer index behind it packs a 16-base key into the high 32 bits of each entry. The end of
# a key's range used to be found as (key + 1) << 32, which wraps to 0 for the all-T key, the only
# one that fills all 32 bits; the lookup then ran off the end of the index and breseq died with a
# segmentation fault right after "Predicting gene conversions...".
#
# input.gd has two SNPs that each complete a 20-base homopolymer. The seed centered on the first
# is T x 20 (forward-strand lookup); the seed centered on the second ends in A x 19, so its
# reverse complement begins with T x 16 (reverse-strand lookup).
#
# Expected: NORMALIZE exits 0 and both SNPs come through unchanged.

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        NORMALIZE \
        -o ${SELF}/output.gd \
        -r ${DATADIR}/gene_conversion/gene_conversion_homopolymer.fna \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
