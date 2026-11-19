#
# A sample with no mutations still gets its own (empty) column in a compare table by default, so a
# table shows which samples were checked and came up clean (#283). --collapse-samples omits it.
#
SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/output.html"
EXPECTED_OUTPUTS[0]="${SELF}/expected.html"
CURRENT_OUTPUTS[1]="${SELF}/output_collapsed.html"
EXPECTED_OUTPUTS[1]="${SELF}/expected_collapsed.html"

TESTCMD="\
    ${GDTOOLS} \
        COMPARE \
        -f HTML \
        -o ${SELF}/output.html \
        -r ${DATADIR}/lambda/lambda.1-2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -r ${DATADIR}/lambda/lambda.4.gbk \
        -r ${DATADIR}/lambda/lambda.5.gbk \
        ${SELF}/../gdtools_compare_3/input_1.gd \
        ${SELF}/empty.gd \
    && ${GDTOOLS} \
        COMPARE \
        -f HTML \
        --collapse-samples \
        -o ${SELF}/output_collapsed.html \
        -r ${DATADIR}/lambda/lambda.1-2.gbk \
        -r ${DATADIR}/lambda/lambda.3.gbk \
        -r ${DATADIR}/lambda/lambda.4.gbk \
        -r ${DATADIR}/lambda/lambda.5.gbk \
        ${SELF}/../gdtools_compare_3/input_1.gd \
        ${SELF}/empty.gd \
    "

do_test $1 ${SELF}
