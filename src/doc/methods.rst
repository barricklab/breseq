Methods
==============

This section describes the algorithms used by :program:`breseq` in more detail for those who want to know how it predicts mutations.

Read Alignment
----------------

By default :program:`breseq` uses :program:`SSAHA2` to map reads to the reference genome sequence with the following alignment options:: 

   ssaha2 -kmer 13 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 1 ...

These parameters mean that only alignments to the reference genome with exact matches of at least 13 bp will be reported, and that attempts to extend these seed alignments that allow gaps and mismatches will use the cross_match algorithm with a word size of one base.

Currently the distance contraints available in paired-end or mate-paired read data are not used during read alignment or as a source of evidence supporting mutational events. These data sets are treated as single-end read data.

New junction (NJ) evidence
-----------------------------

First :program:`breseq` searches for mosaic read alignments that may indicate new junctions between disjoint regions of the reference sequence. 

In a pre-processing step, all read alignments with insertions or deletions of ≥2 bp are split into the separate sub-alignments because gaps larger than a single bp can be problematic for generating accurate and consistent alignments. Next, for each read that has multiple alignments to the reference, all combinations of these alignments are tested to find cases where: (1) there are at least 9 referent bases in the alignments that are not shared by both matches.

In a pre-processing step, all read alignments with insertions or deletions of ≥2 nt are split into the separate sub-alignments. Then every read with multiple alignments to the reference sequence is examined for candidate junctions that it might support. These are compiled and the consensus best number are analyzed.


Read alignment (RA) evidence
------------------------------

Read end trimming
*****************

Alignments of the ends of short reads can be with respect. :program:`breseq` uses a conservative strategy that ignores possibly ambiguous.

Bayesian SNP caller
********************

Thus, :program:`breseq` will only find indels of 1 nt length and base substitutions by this method. 


Missing coverage (MC) evidence
------------------------------

Calibrated by fitting a negative binomial distribution to censored distribution.

*********


Mutational event prediction
---------------------------

Calibrated by fitting a negative binomial distribution to censored distribution.


Limitations
-----------

:program:`breseq` cannot find some types of mutations:

`Mutations in repeat regions` 
	In genomic regions where the only mapped reads also map equally well to other locations in the genome, it is not possible to call mutations. This is an inherent limitation of short-read data. These regions are reported as 'UN' evidence, so that the user can distinguish where in the genome there was not sufficient coverage of uniquely mapped reads to call mutations.
`Chromosomal inversions and rearrangements through repeat sequences`
   These types of mutations cannot be detected when they involve sequence repeats on the order of the read length. Reads that span repeats and uniquely align in the reference sequence on each end are necessary to detect them. :program:`breseq` does not use mate-paired or paired end information to identify these kinds of mutational events.
