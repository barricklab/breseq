#!/bin/bash
#
# This script is a front-end to the testing infrastructure for Breseq.  The main
# body of code lives in ./common.sh.
#
# $1 == test action
# $2 == testdir ('tests' is an alias for all tests)
# $3 == long
#
#
# Example commands:   ./test.sh test tests            #Run normal tests
#                     ./test.sh build tests long      #Run all tests, including long tests
#                     ./test.sh test tests long       #Run normal tests and long tests
#
# load the common testing tools (relative to this script):

echo "3 " $3

TESTDIR=`dirname ${BASH_SOURCE}`
. ${TESTDIR}/common.sh

#For measuring elapsed time of 'make test'
start_timer=$(date '+%s') 

# if ${TESTDIR}/$2/${TESTEXEC} exists, then we're running a single test.  otherwise,
# we're running a batch of tests.

if [[ -e ${TESTDIR}/$2/${TESTEXEC} ]]; then
	${TESTDIR}/$2/${TESTEXEC} $1
else
	for i in `find $2 -name ${TESTEXEC}`; do

		# skip test directories that begin with an underscore!
		if [[ (! $i =~ "/_") && (! $i =~ "^_") ]]; then 
			#Only do long tests if requested
			if [[ ($3 = "long") ]] || [[ (! $i =~ "/long") && (! $i =~ "^long") ]];  then 
				echo "TEST:" $i $1
				$i $1
			fi
		fi

	done
fi

#For measuring elapsed time of 'make test'
end_timer=$(date '+%s') 
elapsed_time=$((end_timer - start_timer))
seconds=$((elapsed_time % 60))
minutes=$(((elapsed_time / 60) % 60))

printf "Elapsed time for 'make test': %02d:%02d\n" $minutes $seconds
