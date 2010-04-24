#!/bin/bash
#
# Common testing functionality for breseq.
#
# Testing in breseq is based on matching sha1 hashes of files ('cause
# the data itself is too big to want to store in hg).  But, some of the output
# includes dates, so we only hash certain files.

# common variables:
# paths must either be relative to the location of this script or absolute.
COMMONDIR=`dirname ${BASH_SOURCE}`
# path to breseq:
BRESEQ="perl ${COMMONDIR}/../src/perl/breseq.pl"
# path to test data:
DATADIR=${COMMONDIR}/../tests/data
# this is a find-compatible list of files that we'll hash:
FILE_PATTERN='-name *.tab -or -name *.html ! -name summary.html'
# executable used to hash files:
HASH=`which sha1sum`
# name of file containing expected hash values:
EXPECTED=expected.sha1
# name of testexec file
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
        echo "Building expected values for test $1..."
        do_build $1
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
# testcmd == function that must be defined in the test command file.
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
        rebuild)
            do_clean $2
            testcmd
            do_build $2
        ;;
        test)
            testcmd
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
    echo "Usage: $0 build|check|clean|rebuild|test <testdir>" >&2
    exit -1
}

