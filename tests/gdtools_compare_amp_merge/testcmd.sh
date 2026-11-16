#
# AMPs at the same endpoints but different new_copy_number are one row in a compare table, and each
# sample's cell is its MEAN copy number, 1 + frequency * (new_copy_number - 1): here a fixed x2 (2), a
# x2 in 75% of the population (1.75) and a fixed x3 (3). GD/PHYLIP/FASTA output is not merged.
#
SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.csv"
EXPECTED_OUTPUTS[0]="${SELF}/expected.csv"
CURRENT_OUTPUTS[1]="${SELF}/output.html"
EXPECTED_OUTPUTS[1]="${SELF}/expected.html"

TESTCMD="\
    ${GDTOOLS} \
        COMPARE \
        -f TABLE \
        -o ${SELF}/output.csv \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${SELF}/input_1.gd \
        ${SELF}/input_2.gd \
        ${SELF}/input_3.gd \
    && ${GDTOOLS} \
        COMPARE \
        -f HTML \
        -o ${SELF}/output.html \
        -r ${DATADIR}/lambda/lambda.gbk \
        ${SELF}/input_1.gd \
        ${SELF}/input_2.gd \
        ${SELF}/input_3.gd \
    "

do_test $1 ${SELF}
