#!/bin/bash

# A GFF3 whose feature rows name the contig as the file spells it, `p(1)`, must load as one
# contig: safe_seq_id_name used to apply to ##sequence-region and the FASTA header and not to
# column 1 of a feature row, so this file registered a second contig with no sequence and
# failed with "No sequence found for reference: p(1)".

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gff3"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gff3"

TESTCMD="\
    ${BRESEQ} \
        CONVERT-REFERENCE \
        -f GFF3 \
        -o ${SELF}/output.gff3 \
        ${DATADIR}/unsafe_seq_ids/unsafe_seq_id.gff3 \
    "

do_test $1 ${SELF}
