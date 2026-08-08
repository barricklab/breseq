#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Reference files may be gzip-compressed, in any supported format. Compression is
# detected from the file's magic bytes rather than its name, so files can be used
# exactly as downloaded -- NCBI's assembly FTP site, for example, serves
# references only as '*_genomic.gbff.gz' / '*_genomic.fna.gz'.
#
# The inputs are gzipped here at test time rather than committed, both to avoid
# adding binary files to git and so the comparison is against the very same bytes
# the uncompressed tests use: this test's expected.gff3 must stay identical to
# tests/breseq_convert_reference_1's apart from the ##original-file-name line,
# which records the input path.
#
# Covers both reference readers: GenBank (cReferenceSequences::ReadGenBank, whose
# helpers take std::istream so a gzip stream can be passed) and FASTA (cFastaFile,
# which derives from flexgzfstream).

CURRENT_OUTPUTS[0]="${SELF}/output.gff3"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gff3"
CURRENT_OUTPUTS[1]="${SELF}/output.fasta"
EXPECTED_OUTPUTS[1]="${SELF}/expected.fasta"

# The gzipped inputs are named 'output.*' so that common.sh's do_clean, which
# removes <testdir>/output*, tidies them up along with the results.
#
# The '##original-file-name' line is stripped before comparison because it records
# the input path, and that path's leading './' depends on how the test was
# invoked -- 'make test' runs testcmd.sh as 'tests/<name>/testcmd.sh' while
# tests/test.sh runs it as './tests/<name>/testcmd.sh', so SELF (and therefore
# the recorded path) differs between the two. It is not what this test is about.
TESTCMD="\
    gzip -c ${DATADIR}/pDCAF3/pDCAF3.gbk > ${SELF}/output.pDCAF3.gbk.gz \
    && gzip -c ${DATADIR}/REL606/REL606.fragment.fna > ${SELF}/output.REL606.fragment.fna.gz \
    && ${BRESEQ} \
        CONVERT-REFERENCE \
        -f GFF \
        -o ${SELF}/output.raw.gff3 \
        ${SELF}/output.pDCAF3.gbk.gz \
    && grep -v '^##original-file-name' ${SELF}/output.raw.gff3 > ${SELF}/output.gff3 \
    && ${BRESEQ} \
        CONVERT-REFERENCE \
        -f FASTA \
        -o ${SELF}/output.fasta \
        ${SELF}/output.REL606.fragment.fna.gz \
    "

do_test $1 ${SELF}
