#!/bin/bash
#
# Run breseq consistency tests directly, without Snakemake -- e.g. to debug
# a single test interactively with output streamed live to the terminal, or
# as a serial fallback to 'make test'.
#
# Usage: ./tests/run.sh <test_name>|all

SELF=`dirname ${BASH_SOURCE}`
TESTNAME="$1"

if [[ -z "${TESTNAME}" ]]; then
	echo "Usage: $0 <test_name>|all" >&2
	exit 1
fi

if [[ "${TESTNAME}" == "all" ]]; then
	START=$(date '+%s')
	for TESTCMD in `export LC_ALL=POSIX; find ${SELF} -mindepth 2 -maxdepth 2 -name testcmd.sh | sort`; do
		TESTDIR=`dirname ${TESTCMD}`
		NAME=`basename ${TESTDIR}`
		[[ "${NAME}" == _* ]] && continue
		[[ "${NAME}" == *_disabled ]] && continue
		echo "=== run: ${NAME} ==="
		${SELF}/run_logged_test.sh ${TESTDIR} ${TESTDIR}/test.log ${TESTDIR}/test.result
	done
	END=$(date '+%s')
	${SELF}/print_test_summary.sh ${SELF} $((END - START))
else
	if [[ ! -e "${SELF}/${TESTNAME}/testcmd.sh" ]]; then
		echo "No such test: ${TESTNAME}" >&2
		exit 1
	fi
	"${SELF}/${TESTNAME}/testcmd.sh" test
fi
