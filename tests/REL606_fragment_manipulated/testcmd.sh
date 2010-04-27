#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../../etc/common.sh

testcmd() {
    ${BRESEQ} \
        -o ${SELF} \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${DATADIR}/REL606/REL606.fragment.1.fastq
}

do_test $1 ${SELF}
