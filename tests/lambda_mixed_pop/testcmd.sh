#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh


TESTCMD="\
    ${BRESEQ} \
    -o ${SELF} \
    -r ${DATADIR}/lambda/lambda.gbk \
    ${DATADIR}/lambda/lambda_mixed_population.fastq \
    "

##     --perl-identify-candidate-junctions \
##    --per_position_file \


do_test $1 ${SELF}
