#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# A GenBank CONTIG (assembly join) line does not mean the record has no sequence.
# NCBI's RefSeq "CON" division records for complete genomes -- which is what the
# assembly FTP site serves as <assembly>_genomic.gbff.gz -- carry a CONTIG join
# AND a full ORIGIN block. breseq used to treat CONTIG as "no sequence" and skip
# to the end of the record, silently discarding the sequence of every such
# reference and then failing with a confusing "No sequence found for reference".
#
# The fixture is synthesized at test time by inserting a CONTIG line immediately
# before ORIGIN in a small committed reference, so no new binary file is needed
# and the input stays in step with the uncompressed tests. awk (not sed) does the
# insertion because inserting a newline in a sed replacement is not portable
# between GNU and BSD/macOS.
#
# expected.fasta is therefore required to be identical to the sequence obtained
# from the unmodified file: if the CONTIG line causes the sequence to be dropped
# again, the conversion either fails outright or emits an empty sequence.

CURRENT_OUTPUTS[0]="${SELF}/output.fasta"
EXPECTED_OUTPUTS[0]="${SELF}/expected.fasta"

TESTCMD="\
    awk '/^ORIGIN/ && !seen { print \"CONTIG      join(KC619530.1:1..8920)\"; seen=1 } { print }' \
        ${DATADIR}/pDCAF3/pDCAF3.gbk > ${SELF}/output.contig.gbk \
    && ${BRESEQ} \
        CONVERT-REFERENCE \
        -f FASTA \
        -o ${SELF}/output.fasta \
        ${SELF}/output.contig.gbk \
    "

do_test $1 ${SELF}
