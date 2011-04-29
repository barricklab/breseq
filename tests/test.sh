#!/bin/bash
#
# This script is a front-end to the testing infrastructure for Breseq.  The main
# body of code lives in ./common.sh.
#
# $1 == test action
# $2 == testdir

# load the common testing tools (relative to this script):
TESTDIR=`dirname ${BASH_SOURCE}`
. ${TESTDIR}/common.sh

# validate input params:
if [[ ($# -ne 2) || (! -d $2) ]]; then
    do_usage
fi

# if $2/${TESTEXEC} exists, then we're running a single test.  otherwise,
# we're running a batch of tests.
if [[ -e $2/${TESTEXEC} ]]; then
	do_test $1 $2
else
	for i in `find $2 -name ${TESTEXEC}`; do

		# skip test directories that begin with an underscore!
		if [[ ${i:${#TESTDIR}+1:1} != '_' ]]; then 
			echo "TEST:" $i $1
			$i $1
		fi
	done
fi
