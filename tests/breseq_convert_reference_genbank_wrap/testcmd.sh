#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gbk"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gbk"

# Writing GenBank must wrap free-text qualifiers only at spaces (and /translation
# without inserting any), so that reading the written file back gives the same
# values. Writing it a second time from the first output checks this: any value
# changed by the round trip changes the second file, and the cmp fails the test.
TESTCMD="\
    ${BRESEQ} \
        CONVERT-REFERENCE \
        -f GENBANK \
        -o ${SELF}/output.gbk \
        ${DATADIR}/genbank_wrap/genbank_wrap.gbk \
    && ${BRESEQ} \
        CONVERT-REFERENCE \
        -f GENBANK \
        -o ${SELF}/output_roundtrip.gbk \
        ${SELF}/output.gbk \
    && perl -i -pe 's/^(LOCUS.+)\d{2}-\w{3}-\d{4}$/\$1XX-XX-XXXX/g' ${CURRENT_OUTPUTS[0]} ${SELF}/output_roundtrip.gbk \
    && cmp ${SELF}/output.gbk ${SELF}/output_roundtrip.gbk \
    "

do_test $1 ${SELF}
