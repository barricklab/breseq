#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gbk"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gbk"

CURRENT_OUTPUTS[1]="${SELF}/output.gd"
EXPECTED_OUTPUTS[1]="${SELF}/expected.gd"

# input.gd also carries a header line, '#=TOPOLOGY<TAB>REL606-5<TAB>CIRCULAR', that sets one
# sequence's topology in the output. This reference has three records and all three are linear, so
# expected.gbk is where "one record of a multi-record reference changed, alongside the mutations in
# the same file" is asserted: REL606-5's LOCUS line must read circular and sfGFP's must stay linear.
# expected.gd asserts the line is carried through to the applied GD.
#
# The second command writes the same thing as GFF3, where topology is the Is_circular attribute of
# each sequence's region feature. That attribute is filled in when the reference is read, so it has
# to be rewritten from the sequence's current topology on output; it once was not, and the GFF3
# said linear for a sequence the GenBank output called circular.
#
# The third command is the error path: bad_topology.gd names a seq_id that is not in the
# reference. APPLY must fail, say which seq_id, and write no output.
TESTCMD="\
    ${GDTOOLS} \
        APPLY \
        -s sfGFP \
        -s REL606-5 \
        -f GENBANK \
        -o ${SELF}/output.gbk \
        -r ${DATADIR}/lambda/lambda.gbk \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        -r ${DATADIR}/sfGFP/sfGFP.gbk \
        --applied-gd ${SELF}/output.gd \
        ${SELF}/input.gd \
        && perl -i -pe 's/^(LOCUS.+)\d{2}-\w{3}-\d{4}$/\$1XX-XX-XXXX/g' ${CURRENT_OUTPUTS[0]} \
    && ${GDTOOLS} \
        APPLY \
        -s sfGFP \
        -s REL606-5 \
        -f GFF3 \
        -o ${SELF}/output.topology.gff3 \
        -r ${DATADIR}/lambda/lambda.gbk \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        -r ${DATADIR}/sfGFP/sfGFP.gbk \
        ${SELF}/input.gd \
    && test \$(grep -c '^REL606-5[[:space:]].*[[:space:]]region[[:space:]].*Is_circular=true' ${SELF}/output.topology.gff3) -eq 1 \
    && test \$(grep -c '^sfGFP[[:space:]].*[[:space:]]region[[:space:]].*Is_circular=false' ${SELF}/output.topology.gff3) -eq 1 \
    && ! ${GDTOOLS} \
        APPLY \
        -f GENBANK \
        -o ${SELF}/output.bad_topology.gbk \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/bad_topology.gd \
        > ${SELF}/output.bad_topology.log 2>&1 \
    && grep -q 'NOT-A-SEQUENCE' ${SELF}/output.bad_topology.log \
    && test ! -e ${SELF}/output.bad_topology.gbk \
    "

do_test $1 ${SELF}
