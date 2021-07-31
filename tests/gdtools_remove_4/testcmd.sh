#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        REMOVE \
        -v \
        -e \
        -c ASSIGNED==1 \
        -o ${SELF}/intermediate.gd \
        ${SELF}/input.gd \
        ; \
    ${GDTOOLS} \
        REMOVE \
        -v \
        -e \
        -c type!=JC \
        -o ${SELF}/output.gd \
        ${SELF}/intermediate.gd \
        ;
    rm ${SELF}/intermediate.gd \
    "

do_test $1 ${SELF}
