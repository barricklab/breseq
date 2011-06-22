The files dm3_3R_4766911_4767130.* are a spliced alignments from an
RNASeq experiment (40bp Illumina single end) measuring alternative
splicing in fly.  They align with chromosome 3R (genome build: dm3)
spanning a region holding a small (3nt) exon on the gene 'skap'
(CG11963).  Some of the spliced alignments join three exons, and cross
two introns; the CIGAR lines match /\d+M124N3M91N\d+M/, where
the internal 3M is the 3bp exon. and the flanking 124N and 91N are the
introns.


NOTE: These files were copied from the examples directory in the
Samtools distribution version 0.1.3. The original README follows.

ORIGINAL 00README.txt

File ex1.fa contains two sequences cut from the human genome
build36. They were extracted with command:

  samtools faidx human_b36.fa 2:2043966-2045540 20:67967-69550

Sequence names were changed manually for simplicity. File ex1.sam.gz
contains MAQ alignments exatracted with:

  (samtools view NA18507_maq.bam 2:2044001-2045500;
   samtools view NA18507_maq.bam 20:68001-69500)

and processed with `samtools fixmate' to make it self-consistent as a
standalone alignment.

To try samtools, you may run the following commands:

  samtools faidx ex1.fa                 # index the reference FASTA
  samtools import ex1.fa.fai ex1.sam.gz ex1.bam   # SAM->BAM
  samtools index ex1.bam                # index BAM
  samtools tview ex1.bam ex1.fa         # view alignment
  samtools pileup -cf ex1.fa ex1.bam    # pileup and consensus
  samtools pileup -cf ex1.fa -t ex1.fa.fai ex1.sam.gz

