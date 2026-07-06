#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gff3"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gff3"

# Regression test for GenBank records without a FEATURES section (issue #421)
# and for merging split annotation/sequence GenBank files. CONVERT-REFERENCE
# merges all inputs into one output, so one command exercises:
#   no_features.gbk    -- multi-record file whose 2nd record has no FEATURES
#                         (header -> ORIGIN -> sequence), which previously
#                         mis-fired the "Multiple LOCUS lines" error.
#   features_only.gbk  -- FEATURES table but no sequence (ends at "//")
#   sequence_only.gbk  -- same seq_id, sequence but no FEATURES
#                         => the two are merged, like GenBank+FASTA loading.
TESTCMD="\
    ${BRESEQ} \
        CONVERT-REFERENCE \
        -f GFF \
        -o ${SELF}/output.gff3 \
        ${DATADIR}/genbank_no_features/no_features.gbk \
        ${DATADIR}/genbank_no_features/features_only.gbk \
        ${DATADIR}/genbank_no_features/sequence_only.gbk \
    "

do_test $1 ${SELF}
