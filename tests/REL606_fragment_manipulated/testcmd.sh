#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

testcmd() {
    ${BRESEQ} \
        -o ${SELF} \
        -n REL606_fragment_manipulated \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${DATADIR}/REL606/REL606.fragment.1.fastq
}

do_test $1 ${SELF}
