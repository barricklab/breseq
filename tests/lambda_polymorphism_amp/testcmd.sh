#
# Polymorphism-mode amplifications with a calculated frequency (barricklab/conda#2: the AMPs built
# from JC + CN evidence carried no frequency, and the predictor asserted).
#
# A simulated population: 75% of genomes carry amplified.gd (a 5 kb block x2 and a 3 kb block x3),
# 25% are the reference. Reads are drawn at 60x from the amplified genome and 20x from the reference.
# Expect both AMPs at frequency ~0.75, the first with new_copy_number 2 and the second with 3 -- the
# smallest whole copy number that keeps the frequency <= 1 for a mean copy number of 2.5 -- each with
# the jc_frequency and cn_frequency they were combined from. The reads are generated at test time,
# into output_sim/ (removed by do_clean), so no read data is committed.
#
# The reference-genome reads are deliberately NOT named reference_*: stage 02 removes its bowtie2
# index with "rm -f 02_reference_alignment/reference*", which also deletes such a read file's
# reference_reads.stage1.bam.
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
        -o ${SIMDIR}/amplified.fasta \
        ${SELF}/amplified.gd \
    && ${BRESEQ} SIMULATE-READS \
        ${REFERENCE_ARG} \
        -c 20 -l 100 --seed 1 \
        -o ${SIMDIR}/ancestor_reads.fastq \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SIMDIR}/amplified.fasta \
        -c 60 -l 100 --seed 2 \
        -o ${SIMDIR}/amplified_reads.fastq \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -p \
        -o ${SELF} \
        ${REFERENCE_ARG} \
        ${SIMDIR}/ancestor_reads.fastq \
        ${SIMDIR}/amplified_reads.fastq \
    "

do_test $1 ${SELF}
