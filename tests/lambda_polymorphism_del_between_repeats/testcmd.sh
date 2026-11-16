#
# Polymorphism-mode deletion between two direct copies of a repeat, which leaves neither missing
# coverage (only part of the population has it) nor a unique junction (the crossover lies inside the
# repeat), so its only evidence is a CN region below single copy whose ends sit at the two copies.
#
# The reference is lambda with two IS150 copies inserted (insert_two_is150.gd, the element taken from
# REL606.is.gbk), so both are annotated repeat_regions. A simulated population: 50% of genomes carry
# deleted.gd, the recombination of those copies (the sequence between them plus the second copy),
# 50% do not. Expect DEL 16446 5446 between=IS150 at frequency ~0.5, from cn_frequency alone.
#
# The deletion is kept to about 10% of the genome: an 11 kb one in the middle of this 51 kb linear
# sequence was partly absorbed by CNery's coverage-gradient correction (0.82 vs 1.2-1.3 around it).
#
SELF=`dirname ${BASH_SOURCE}`
TEST_CORES=4
. ${SELF}/../common.sh

CURRENT_OUTPUTS[0]="${SELF}/data/annotated.gd"
EXPECTED_OUTPUTS[0]="${SELF}/expected.gd"
SIMDIR="${SELF}/output_sim"

TESTCMD="\
    mkdir -p ${SIMDIR} \
    && ${GDTOOLS} APPLY \
        -f GENBANK \
        -r ${DATADIR}/lambda/lambda.gbk \
        -r ${DATADIR}/REL606/REL606.is.gbk \
        -s NC_001416 \
        -o ${SIMDIR}/lambda_two_is150.gbk \
        ${SELF}/insert_two_is150.gd \
    && ${GDTOOLS} APPLY \
        -f FASTA \
        -r ${SIMDIR}/lambda_two_is150.gbk \
        -o ${SIMDIR}/deleted.fasta \
        ${SELF}/deleted.gd \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SIMDIR}/lambda_two_is150.gbk \
        -c 40 -l 100 --seed 5 \
        -o ${SIMDIR}/ancestor_reads.fastq \
    && ${BRESEQ} SIMULATE-READS \
        -r ${SIMDIR}/deleted.fasta \
        -c 40 -l 100 --seed 6 \
        -o ${SIMDIR}/deleted_reads.fastq \
    && ${BRESEQ} \
        ${BRESEQ_TEST_THREAD_ARG} \
        -p \
        -o ${SELF} \
        -r ${SIMDIR}/lambda_two_is150.gbk \
        ${SIMDIR}/ancestor_reads.fastq \
        ${SIMDIR}/deleted_reads.fastq \
    "

do_test $1 ${SELF}
