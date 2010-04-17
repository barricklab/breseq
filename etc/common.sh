#!/bin/bash
#
# Common testing functionality.

# define common variables; paths must be relative the location of this script!
COMMONDIR=`dirname ${BASH_SOURCE}`
BRESEQ=${COMMONDIR}/../src/perl/breseq.pl

# testing in breseq is based on matching sha1 hashes; which files should be hashed?
# this pattern must be suitable for "find".
FILE_PATTERN='-name *.tab -or -name *.html ! -name summary.html'
HASH=`which shasum`
EXPECTED=expected.sha1
TESTEXEC=testcmd.sh


# build the list of hashes
# $1 == testdir
#
do_build() {
    pushd $1 > /dev/null
    for i in `find . ${FILE_PATTERN}`; do
        ${HASH} $i
    done > ${EXPECTED}
    popd > /dev/null
}


# check current hashes against expected values.
# $1 == testdir
#
do_check() {
    if [[ ! -e $1/${EXPECTED} ]]; then
				echo "Test $1 does not have a list of expected values (${EXPECTED})."
        do_usage
    fi
    pushd $1 > /dev/null
    if ! ${HASH} --check ${EXPECTED} > /dev/null; then
        popd > /dev/null
        echo "Failed check: $1"
        exit -1
    fi
    popd > /dev/null
    echo "Passed check: $1"
    exit 0
}


# clean a test
# $1 == testdir
#
do_clean() {
	rm -R $1/0* $1/output
}


# main test-running method
# $1 == action
# $2 == testdir
# $3 == command-line to run test
#
do_test() {
	case "$1" in
    build)
        do_build $2
    ;;
    check)
        do_check $2
    ;;
    clean)
        do_clean $2
    ;;
    test)
        do_clean $2
        testcmd # defined by $2/${TESTEXEC}
        do_check $2
    ;;
    *)
        do_usage
    ;;
	esac
}

# print usage information.
#
do_usage() {
    echo "Usage: $0 build|check|clean|test <testdir>" >&2
    exit -1
}

