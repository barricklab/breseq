#!/bin/bash
#
# Print a pass/fail + timing summary for a test run, reading the
# tests/<name>/test.result files written by tests/run_logged_test.sh, and
# exit with a status suitable for continuous integration (e.g. GitHub
# Actions): 0 if every test passed, 1 if any test failed.
#
# $1 == tests directory (e.g. "tests")
# $2 == total wall-clock elapsed seconds for the whole run (optional)

TESTSDIR="$1"
TOTAL_ELAPSED="${2:-0}"

format_time() {
	local s=$1
	printf '%02d:%02d:%02d' $((s / 3600)) $(((s % 3600) / 60)) $((s % 60))
}

PASS_COUNT=0
FAIL_COUNT=0

echo ""
echo "========================================================================="
echo "TEST SUMMARY"
echo "========================================================================="

for RESULT in $(find "${TESTSDIR}" -mindepth 2 -maxdepth 2 -name test.result | sort); do
	read STATUS NAME ELAPSED < "${RESULT}"
	if [[ "${STATUS}" == "PASS" ]]; then
		PASS_COUNT=$((PASS_COUNT + 1))
		printf "  %-6s %-50s %s\n" "${STATUS}" "${NAME}" "$(format_time ${ELAPSED})"
	else
		FAIL_COUNT=$((FAIL_COUNT + 1))
		printf " >>%-6s %-50s %s\n" "${STATUS}" "${NAME}" "$(format_time ${ELAPSED})"
	fi
done

echo "-------------------------------------------------------------------------"
printf "  %d passed, %d failed -- total time %s\n" "${PASS_COUNT}" "${FAIL_COUNT}" "$(format_time ${TOTAL_ELAPSED})"
echo "========================================================================="
echo ""

if [[ ${FAIL_COUNT} -eq 0 ]]; then
	echo "RESULT: ALL TESTS PASSED"
	exit 0
else
	echo "RESULT: ${FAIL_COUNT} TEST(S) FAILED"
	exit 1
fi
