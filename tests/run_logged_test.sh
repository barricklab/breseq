#!/bin/bash
#
# Run a single test, redirecting its output to a log file and recording
# pass/fail + elapsed time to a result file. Used by the Snakefile-based
# parallel test runner (see tests/Snakefile and 'make test').
#
# $1 == test directory (e.g. tests/lambda_mixed_pop)
# $2 == path to write the log file
# $3 == path to write the result file
#
TESTDIR="$1"
LOG="$2"
RESULT="$3"

TESTNAME=$(basename "${TESTDIR}")

START=$(date '+%s')
"${TESTDIR}/testcmd.sh" test > "${LOG}" 2>&1
STATUS=$?
END=$(date '+%s')
ELAPSED=$((END - START))

# Always record a result, regardless of pass/fail, so the summary step can
# report on every test that ran. This file is intentionally not declared as
# a Snakemake 'output:' so it survives Snakemake's cleanup of failed jobs'
# outputs.
if [[ ${STATUS} -eq 0 ]]; then
	echo "PASS ${TESTNAME} ${ELAPSED}" > "${RESULT}"
else
	echo "FAIL ${TESTNAME} ${ELAPSED}" > "${RESULT}"
fi

exit ${STATUS}
