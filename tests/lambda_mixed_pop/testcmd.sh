#!/bin/bash
#
# The baseline consensus-mode run on the lambda mixed population reads. One run carries three
# things that used to be three tests on the same reads with the same golden:
#
# (1) The consensus calls themselves, with a header .gd and seq ids taken from the VERSION field.
#     tests/lambda_mixed_pop_cl_tabulate reads this run's data/reference.{bam,fasta,gff3}.
#
# (2) Copy number (CN) evidence in consensus mode (was lambda_mixed_pop_cn_evidence).
#     --predict-copy-number is deprecated and CN is on by default, but it is what sets
#     copy_number_explicitly_requested, which makes a MISSING CNery fatal instead of a
#     warning-and-skip -- a CN test that passes quietly on a machine without CNery is worth
#     nothing. -k keeps 09_copy_number_variation/ for inspection when the CN line moves.
#
# (3) The four deprecated opt-in flags (was breseq_deprecated_predict_flags). They must still be
#     accepted, must warn, and must change nothing -- they are now what the run does anyway. A
#     deprecated option that errors out breaks every existing script and pipeline that passes it,
#     which is the whole point of keeping it registered (DEPRECATED_OPTION in anyoption.h: parsed
#     normally, never shown in help). The grep -c asserts all four warnings fire. Without it the
#     test would still pass if a flag were silently dropped from the option table, since AnyOption
#     exits on an UNKNOWN option but this golden would be identical either way. "Change nothing"
#     was checked when the tests were merged: adding the flags left this expected.gd untouched.
#
# (4) That the pileups split across threads really were split. Their results do not depend on the
#     chunking, so a change that quietly sent them back down the serial path would leave every
#     golden untouched. The second grep requires the coverage, error-count and evidence pileups
#     to each report more than one chunk ("chunks", plural), which this 48.5-kb reference gets
#     at -j 4 only because kMinimumPileupChunkLength is small.
#
# stderr goes to a file for those greps and is then replayed, so a failure still shows up in
# test.log.

SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"
COMPARE_ARG="--genbank-field-for-seq-id VERSION"

TESTCMD="\
    { ${BRESEQ} \
        -k \
        ${BRESEQ_TEST_THREAD_ARG} \
        -o ${SELF} \
        -g ${SELF}/header.gd \
        --genbank-field-for-seq-id VERSION \
        --predict-copy-number \
        --predict-discordant-pairs \
        --predict-missing-pairs \
        --predict-pair-distance \
        ${REFERENCE_ARG} \
        ${DATADIR}/lambda/lambda_mixed_population.fastq.gz \
        2> ${SELF}/output.stderr.txt ; \
      BRESEQ_STATUS=\$? ; \
      cat ${SELF}/output.stderr.txt >&2 ; \
      test \$BRESEQ_STATUS -eq 0 ; } \
    && test \$(grep -c 'option is DEPRECATED' ${SELF}/output.stderr.txt) -eq 4 \
    && test \$(grep -c 'Piled up [0-9]* chunks of the reference' ${SELF}/output.stderr.txt) -ge 3 \
    "

do_test $1 ${SELF}
