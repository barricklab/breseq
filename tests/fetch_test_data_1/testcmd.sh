#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# Offline self-test of tests/fetch_test_data.sh -- see run_checks.sh for what it
# pins down. Deliberately a NORMAL test, not a 'long' one: it uses curl's file://
# support against small files already in the repository, so it needs no network
# and runs in well under a second. That is what lets CI verify the download,
# checksum-mismatch and failure paths on every pull request, even though the real
# data sets those paths serve are only fetched by 'make test-long'.
#
# Note it sets no TEST_DATA of its own (it drives the fetcher directly, with its
# own generated fixture manifest), so it never touches the real manifest and
# tests/Snakefile never schedules a download for it.

CURRENT_OUTPUTS[0]="${SELF}/output.txt"
EXPECTED_OUTPUTS[0]="${SELF}/expected.txt"

TESTCMD="\
    ${SELF}/run_checks.sh ${SELF} \
    "

do_test $1 ${SELF}
