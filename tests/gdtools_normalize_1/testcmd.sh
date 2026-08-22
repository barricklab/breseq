#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

# NORMALIZE must not let a 'within' mutation steer a neighboring deletion's coordinate.
#
# input.gd puts a SNP inside the new IS1 copy created by the MOB (within=1). For 'within=<MOB id>'
# the POSITION field is not a reference coordinate -- it is the MOB's position plus an offset into
# the inserted element -- so that SNP occupies no reference bases and cannot legitimately bound the
# DEL. It nonetheless falls inside the DEL's reference interval (20030..20529), which used to send
# the DEL to 20300-500 = 19800: 230 bases left of where the author put it, deleting sequence it was
# never equivalent to, and tracking the SNP's coordinate so the same authored DEL landed elsewhere
# in every sample that happened to carry a mutation inside it.
#
# Expected: the DEL only right-shifts within its equivalence window (20030 -> 20032), and NORMALIZE
# exits 0 -- the post-normalization validity check would otherwise reject the output.

CURRENT_OUTPUTS[0]="${SELF}/output.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"

TESTCMD="\
    ${GDTOOLS} \
        NORMALIZE \
        -o ${SELF}/output.gd \
        -r ${DATADIR}/REL606/REL606.fragment.gbk \
        ${SELF}/input.gd \
    "

do_test $1 ${SELF}
