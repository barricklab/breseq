#!/bin/bash

SELF=`dirname ${BASH_SOURCE}`
. ${SELF}/../common.sh

${BRESEQ} CONVERT-REFERENCE --format FASTA -o ${DATADIR}/bull/bull_1.fasta ${DATADIR}/bull/bull_1.gbk
samtools faidx ${DATADIR}/bull/bull_1.fasta

bowtie2-build ${DATADIR}/bull/bull_1.fasta ${DATADIR}/bull/bull_1_index
bowtie2 -x ${DATADIR}/bull/bull_1_index -U ${DATADIR}/bull/bull_1.fastq -S ${DATADIR}/bull/bull_1.sam

samtools view -Sb ${DATADIR}/bull/bull_1.sam > ${DATADIR}/bull/bull_1.unsorted.bam
samtools sort --threads 4 ${DATADIR}/bull/bull_1.unsorted.bam > ${DATADIR}/bull/bull_1.bam
samtools index ${DATADIR}/bull/bull_1.bam
rm ${DATADIR}/bull/bull_1_index*
rm ${DATADIR}/bull/bull_1.unsorted.bam
rm ${DATADIR}/bull/bull_1.sam
