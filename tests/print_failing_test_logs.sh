#!/bin/bash
#
# Print the captured log of every test that did NOT pass in the last run, reading
# the tests/<name>/test.result files written by tests/run_logged_test.sh and the
# tests/<name>/test.log files Snakemake pointed that script at.
#
# This exists because a failing 'make test' otherwise says nothing diagnosable in
# CI: tests/Snakefile sends each test's output to tests/<name>/test.log, and a
# GitHub Actions job log only ever showed "(command exited with non-zero exit
# code)". A real failure -- a breseq FATAL ERROR with its stack trace, or a
# golden-file diff -- was reachable only by reproducing the run locally, which
# for a platform-specific failure means it was not reachable at all.
#
# Run it by hand after a failed 'make test' too; it is the same output either way
# ('::group::' lines are a GitHub Actions fold marker and harmless in a terminal).
#
# $1 == tests directory (e.g. "tests")
# $2 == max lines to tail from each log (optional; default 300, or $BRESEQ_TEST_LOG_TAIL)
#
# Always exits 0: this is a reporter, not a gate. tests/print_test_summary.sh
# owns the pass/fail exit status for the run.

TESTSDIR="${1:-tests}"
TAIL_LINES="${2:-${BRESEQ_TEST_LOG_TAIL:-300}}"

# Total cap on the extracted-banner section per log. A run can produce many
# banners -- one 'Failed check' per compared output file, or, if breseq dies
# inside its own segfault handler, a FATAL ERROR block per repetition -- and the
# point of this section is to put the first error on screen, not to reproduce the
# log (which the CI artifact keeps in full anyway).
BANNER_MAX_LINES="${BRESEQ_TEST_LOG_BANNER_MAX:-400}"

# Banners worth pulling out of the middle of a log, where a plain tail would miss
# them. The first two are written by my_error_handler() (src/breseq/common.h);
# the last two by do_breseq/do_check (tests/common.sh).
BANNER_PATTERN='> FATAL ERROR <|> STACK TRACE <|^XXXXXXXX|^Failed check'

FAILED_ANY=0

for RESULT in $(find "${TESTSDIR}" -mindepth 2 -maxdepth 2 -name test.result | sort); do
	read STATUS NAME ELAPSED < "${RESULT}"
	[[ "${STATUS}" == "PASS" ]] && continue

	FAILED_ANY=1
	LOG="$(dirname "${RESULT}")/test.log"

	echo "::group::${STATUS} ${NAME} -- ${LOG}"

	if [[ ! -e "${LOG}" ]]; then
		echo "No log file at ${LOG} -- the test produced none, or it was cleaned."
	elif [[ ! -s "${LOG}" ]]; then
		echo "Log file ${LOG} is empty."
	else
		# The banner and its surroundings first, so the actual error is at the top
		# of the fold even in a log whose tail is unrelated cleanup chatter.
		if grep -qE "${BANNER_PATTERN}" "${LOG}"; then
			echo "----- error banners in ${LOG} (first ${BANNER_MAX_LINES} lines) -----"
			grep -nE -B 5 -A 25 "${BANNER_PATTERN}" "${LOG}" | head -n "${BANNER_MAX_LINES}"
			echo ""
		fi
		echo "----- last ${TAIL_LINES} lines of ${LOG} -----"
		tail -n "${TAIL_LINES}" "${LOG}"
	fi

	echo "::endgroup::"
done

if [[ ${FAILED_ANY} -eq 0 ]]; then
	echo "No failing tests found under ${TESTSDIR} -- nothing to report."
fi

exit 0
