#!/bin/bash
#
# This script is a front-end to the testing infrastructure for Breseq.  The main
# body of code lives in ./common.sh.
#
# $1 == test action
# $2 == testdir ('tests' is an alias for all tests)
# $3 == '' | long | all   (which extra categories to include)
#
# Batch discovery always skips '_'-prefixed directories. Two extra categories
# are opt-in via $3:
#   ''     -> skip both 'long'-named and '*_disabled' directories (normal set)
#   long   -> also include 'long'-named directories
#   all    -> include everything: 'long'-named AND '*_disabled' directories.
#             This is a clean-time escalation (used by 'make clean-tests' to
#             wipe every test dir's output); there is no equivalent 'make'
#             target that RUNS '*_disabled' tests -- run those by name instead.
#
# Example commands:   ./test.sh test tests            #Run normal tests
#                     ./test.sh test tests long       #Run normal tests and long tests
#                     ./test.sh clean tests all       #Clean every test dir (incl. disabled)
#
# load the common testing tools (relative to this script):

TESTDIR=`dirname ${BASH_SOURCE}`
. ${TESTDIR}/common.sh

#For measuring elapsed time of 'make test'
start_timer=$(date '+%s')

# if ${TESTDIR}/$2/${TESTEXEC} exists, then we're running a single test.  otherwise,
# we're running a batch of tests.

if [[ -e ${TESTDIR}/$2/${TESTEXEC} ]]; then
	${TESTDIR}/$2/${TESTEXEC} $1
else
	for i in `export LC_ALL=POSIX; find $2 -name ${TESTEXEC} | sort`; do

		NAME=`basename $(dirname $i)`

		# skip test directories that begin with an underscore!
		[[ "${NAME}" == _* ]] && continue

		# skip '*_disabled' directories unless running the 'all' category
		[[ "${NAME}" == *_disabled && "$3" != "all" ]] && continue

		# only do 'long' tests when requested ('long' or 'all')
		if [[ "${NAME}" == *long* && "$3" != "long" && "$3" != "all" ]]; then
			continue
		fi

		echo "TEST:" $i $1
		$i $1

	done
fi

#For measuring elapsed time of 'make test'
end_timer=$(date '+%s')
elapsed_time=$((end_timer - start_timer))
seconds=$((elapsed_time % 60))
minutes=$(((elapsed_time / 60) % 60))

printf "Elapsed time for 'make test': %02d:%02d\n" $minutes $seconds
