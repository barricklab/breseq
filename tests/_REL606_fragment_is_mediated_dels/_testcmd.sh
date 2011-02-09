#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

TESTCMD="\
    ${BRESEQ} \
        --force-quality-scores \
        -o ${SELF} \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${DATADIR}/REL606/REL606.fragment.3.fastq \
    "

do_test $1 ${SELF}
