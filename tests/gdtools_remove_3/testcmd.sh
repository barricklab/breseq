#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        REMOVE \
        --comment \
        -c \"snp_type!=intergenic\" \
        -c \"snp_type!=synonymous\" \
        -c FAKE_FIELD==UNDEFINED \
        -o ${SELF}/output.gd \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
