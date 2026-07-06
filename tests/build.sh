#!/bin/bash
#
# Regenerate a test's expected output (expected.gd) from whatever output is
# currently present in its test directory -- e.g. to promote a failed
# 'make test' run's actual (but correct) output to be the new expected
# result. This does NOT re-run breseq/gdtools; for that, use
# './tests/rebuild.sh' instead. See do_build in common.sh.
#
# Usage: ./tests/build.sh <test_name>|all

SELF=`dirname ${BASH_SOURCE}`
ACTION=build
TESTNAME="$1"

if [[ -z "${TESTNAME}" ]]; then
	echo "Usage: $0 <test_name>|all" >&2
	exit 1
fi

if [[ "${TESTNAME}" == "all" ]]; then
	for TESTCMD in `export LC_ALL=POSIX; find ${SELF} -mindepth 2 -maxdepth 2 -name testcmd.sh | sort`; do
		NAME=`basename $(dirname ${TESTCMD})`
		[[ "${NAME}" == _* ]] && continue
		[[ "${NAME}" == *_disabled ]] && continue
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
