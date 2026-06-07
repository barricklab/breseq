#!/bin/bash
#
# Rebuild a test's expected output (expected.gd) from scratch: cleans the
# test directory, re-runs breseq/gdtools, and copies the fresh output over
# expected.gd. Use this after intentional changes to breseq's behavior. To
# instead promote an existing (e.g. failed-run) output without re-running,
# use './tests/build.sh'. See do_test's 'rebuild' case in common.sh.
#
# Usage: ./tests/rebuild.sh <test_name>|all

SELF=`dirname ${BASH_SOURCE}`
ACTION=rebuild
TESTNAME="$1"

if [[ -z "${TESTNAME}" ]]; then
	echo "Usage: $0 <test_name>|all" >&2
	exit 1
fi

if [[ "${TESTNAME}" == "all" ]]; then
	for TESTCMD in `export LC_ALL=POSIX; find ${SELF} -mindepth 2 -maxdepth 2 -name testcmd.sh | sort`; do
		NAME=`basename $(dirname ${TESTCMD})`
		[[ "${NAME}" == _* ]] && continue
		echo "=== ${ACTION}: ${NAME} ==="
		"${TESTCMD}" ${ACTION}
	done
else
	if [[ ! -e "${SELF}/${TESTNAME}/testcmd.sh" ]]; then
		echo "No such test: ${TESTNAME}" >&2
		exit 1
	fi
	"${SELF}/${TESTNAME}/testcmd.sh" ${ACTION}
fi
