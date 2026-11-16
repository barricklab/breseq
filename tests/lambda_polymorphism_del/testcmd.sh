#
# Polymorphism-mode deletions that leave no missing coverage, because only part of the population
# carries them, predicted from their junction (and, when CN can resolve them, the relative coverage).
#
# A simulated population: 50% of genomes carry deleted.gd (a 3 kb deletion, which CN sees at half
# copy, and an 80 bp one, too small for CN's 100 bp windows), 50% are the reference. Reads are drawn
# at 40x from each. Expect both DELs at frequency ~0.5: the 3 kb one from jc_frequency and
# cn_frequency combined, the 80 bp one from jc_frequency alone. The reads are generated at test time,
# into output_sim/ (removed by do_clean), so no read data is committed.
#
SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
REFERENCE_ARG="-r ${DATADIR}/lambda/lambda.gbk"
SIMDIR="${SELF}/output_sim"

TESTCMD="\
    mkdir -p ${SIMDIR} \
    && ${GDTOOLS} APPLY \
        -f FASTA \
        ${REFERENCE_ARG} \
        -o ${SIMDIR}/deleted.fasta \
        ${SELF}/deleted.gd \
    && ${BRESEQ} SIMULATE-READS \
        ${REFERENCE_ARG} \
        -c 40 -l 100 --seed 3 \
        -o ${SIMDIR}/reference_reads.fastq \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SIMDIR}/deleted.fasta \
        -c 40 -l 100 --seed 4 \
        -o ${SIMDIR}/deleted_reads.fastq \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -p \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        ${SIMDIR}/reference_reads.fastq \
        ${SIMDIR}/deleted_reads.fastq \
    "

do_test $1 ${SELF}
