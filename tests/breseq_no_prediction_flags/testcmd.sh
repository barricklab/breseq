#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

# All four new opt-out flags at once. This is what proves the opt-outs actually work: CN, DP,
# MP and PD are on by default now, so without them this run would produce CN evidence (the
# default-on flip is otherwise only ever exercised in the positive direction).
#
# The trailing count is the real assertion. The golden alone would not catch a regression that
# made one of these flags a no-op, because the reviewer would have to notice one line appearing
# in a 100-line .gd -- whereas "an evidence type that must be absent is present" fails loudly.
# Uses cut -f1 (tab is its default delimiter) rather than a grep pattern containing a literal
# tab, which does not survive being stored in TESTCMD and eval'd.
TESTCMD="\
    ${BRESEQ} \
    ${BRESEQ_TEST_THREAD_ARG} \
    -o ${SELF} \
    --no-copy-number-prediction \
    --no-discordant-pair-prediction \
    --no-missing-pair-prediction \
    --no-pair-distance-prediction \
    -r ${DATADIR}/lambda/lambda.gbk \
    ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
    && test \$(cut -f1 ${SELF}/data/annotated.gd | grep -cE '^(CN|DP|MP|PD)\$') -eq 0 \
    "

do_test $1 ${SELF}
