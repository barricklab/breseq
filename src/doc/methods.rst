Methods
==============

This section describes the algorithms used by :program:`breseq` in more detail.

Read Alignment
----------------

By default :program:`breseq` uses :program:`SSAHA2` to map reads to the reference genome sequence with the following alignment options:: 

   ssaha2 -kmer 13 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 1 ...

These parameters mean that `only alignments to the reference genome with exact matches of at least 13 bp will be reported`, and that attempts to extend these seed alignments that allow gaps and mismatches will use the cross_match algorithm with a word size of one base. 

Currently the distance contraints available in paired-end or mate-paired libraries are not used during read alignment or as a source of evidence supporting mutations. These data sets are treated as single-end reads.

:program:`breseq` keeps track of two kinds of read alignments:

`unique read matches` 
	Where a read aligns best to only one location in the reference sequence.
`degenerate read matches`
    Where a read aligns best to multiple, equivalent locations in the reference sequence.
    
For many calculations, :program:`breseq` is concerned with:

`unique-only reference positions` 
	Bases in the reference sequence that ``do not`` overlap any degenerate read matches. 

.. _new-junction-evidence:   
    
New junction (NJ) evidence
-----------------------------

First :program:`breseq` searches for mosaic read alignments that may indicate new junctions between disjoint regions of the reference sequence. 

Identifying candidate junctions
*******************************

In a pre-processing step, all read alignments with insertions or deletions of ≥2 bp are split into the separate sub-alignments. This step ensures that this type of mutation will be found by as new junction evidence, rather than read alignment evidence. This strategy tends to be more accurate because gaps larger than a single bp can be problematic for generating accurate and consistent read alignments, especially when there are indels involving simple sequence repeats. 

Next, for each read that has multiple alignments to the reference, all pairs of these alignments are tested to find cases where: 

#. Neither is completely contained within the other in read coordinates.
#. One contains at least 10 read bases that do not overlap the other. 
#. There are at most 20 bp unique to the read between matches to the reference.

If a pair of alignments passes this test, :program:`breseq` generates the putative sequence of the new junction from the reference sequence and any intervening base pairs that are unique to the read. In cases where the two read alignments overlap because some of its sequence could be assigned to match either location in the reference genome, breseq attempts to assign as much of this overlap to the side of the read that maps uniquely to the reference genome as possible. It considers a location in the reference genome non-unique if it has degenerate reference matches. If ``repeat_region`` annotation exists in the reference genome, then program:`breseq` prefers to have junctions exactly overlap their boundaries as possible, because this represents the simplest mutation possible in cases of new transposon insertions.

This candidate junction sequence includes as many flanking reference bases on each end as the longest read in the entire dataset. (So, if the dataset consists of 36-bp reads, and the two read alignments overlapped by five bases, the sequence for this candidate junction would be 36 x 2 - 5 = 67 bases.

After processing all reads in this manner, :program:`breseq` combines all candidate junctions that match the same reference sequence and calculates a "position-hash" score. This score is a count of the number of `different` start positions in the reference sequence that are observed among the reads that support a candidate junction. This scoring scheme favors candidate junctions supported by reads that are evenly distributed on each strand of the reference genome and evenly distributed at different positions relative to the junction point. (Pathological candidates tend to be supported only by reads that barely overlap the junction and are all on one strand.)

:program:`breseq` sorts all candidate junctions according to this score, breaking ties by favoring junctions with more reads, and retains the top-scoring candidate junctions until adding new candidates would cause their cumulative length to exceed 0.1x the total reference sequence length or thier number to exceed 5000.

Scoring and accepting junctions
*******************************

New junctions may be supported by reads that do not overlap both sides sufficiently to seed alignments during mapping. To include these, :program:`breseq` performs a second :program:`SSAHA2` read-mapping step where it aligns all reads to the file of candidate junction sequences. Then, for each read, it determined whether its best alignment is to a junction candidate or to the reference sequence. A position-hash score is calculated again for each candidate junction by counting the number of different start positions that are observed among the reads that map to a candidate junction. Candidate junctions are accepted as evidence for a mutational event if their position-hash score exceeds a threshold calculated from the distribution of read coverage across reference sites that do not have any degenerate read matches.

This score cutoff is calculated as follows (R-style pseudocode)::

   max_junction_score = int(2 * avg_read_length);
   pr_at_least_one = 1-pnbinom(0, size = nb_fit_size/max_junction_score, mu = nb_fit_mu/max_junction_score);
   junction_accept_cutoff = qbinom(0.01, max_junction_score, pr_at_least_one);

This calculation begins by fitting a censored negative-binomial (overdispersed Poisson) to the distribution of read coverage at unique-only reference positions as described under :ref:`read-coverage`. From this distribution :program:`breseq` calculates the chance that at least one read will start any given position, on a given strand, assuming that level of average coverage. The chance of observing a given position-hash score is then calculated according to the binomial distribution assuming 2 x read length trials, and this chance per trial of observing a read in this register. By default a cutoff of 0.01 is used to establish the cutoff for this test.

For junctions that pass this scoring cutoff, the ends of reads aligning to the junction are re-added as split sub-alignments to the BAM alignment database, resolving ambiguous alignments, such that each read base aligns to only one reference nucleotide. These reads can be recognized in read alignment output because they are renamed with suffixes -M1 and -M2 for the two portions.

.. _read-alignment-evidence:

Read alignment (RA) evidence
------------------------------

:program:`breseq` calls base substitution mutations


Read end trimming
*****************

Alignments of of short reads can be ambiguous with respect to insertion and deletion mutations :ref:`ambiguous-end-figure`. :program:`breseq` uses a conservative strategy to ignores possibly ambiguous bases at the ends of reads when calling mutations within read alignments.

It examines the reference sequence for perfect sequence repeats with lengths of 1-18 bases. Then for each position in the reference it determines how many bases must be trimmed from the end of a read beginning or ending at that position until the remaining bases are unambiguously aligned with respect to possible mutations causing changes in repeat sequences of these lengths. The minimum number of bases trimmed at each end of any read is 1, because one can never unambiguously know if another copy of that base was inserted by a mutation.

Here is an example, showing the logic of end trimming:

.. figure:: images/end_trimming_example.png
   :width: 450px
   :align: center
   
   **End trimming.**

This example shows the number of bases that will be trimmer from the left and right ends of a read if its match to the reference genome begins or ends on that base. (Note that the strand of the genome that the read matches makes no difference!)  The green, blue, and yellow highlight where the numbers come from for three case. 

For green, a read with its left end aligned to this position is not informative with respect to how many AG copies there are in the sequenced genome. Therefore, it is only unambiguously aligned at the bases starting CAT-, and the first four bases will be trimmed. Similarly, a read with its right end aligned to the green position cannot tell how many TA copies there are. It will only be unambiguously aligned through -CTT, and it's last four bases will be trimmed.

Trimming these ends enables more accurate mutation prediction because reads extending into these repeats from either side, but not completely crossing them, would otherwise indicate that there is evidence *against* an indel mutation in the repeat. 

For example, consider this mutation, which involves insertion of a new AGC at a site where there are already two AGC copies:

.. figure:: images/missed_mutation_no_trimming.png
   :width: 600px
   :align: center
   
   **Indel mutation prediction aided by end trimming.**
	
This image shows reads 1-6 aligned to the reference genome with and without end trimming (lowercase letters in reads). Two reads cross the entire AGCx2 repeat and show that a third AGC has been inserted.

Without end trimming, two reads on the top strand that do not cross the new AGC insertion, contradict that there was any change to the sequence here when they are aligned to the reference. With end trimming, these bases are ignored because they are ambiguous with respect to possible insertions, like the event that happened, or deletion of one AGC copy.


Base quality re-calibration
***************************

In the FASTQ input files, each read base has been assigned a quality score by the built-in pipeline for a given sequencing technology. Base quality re-calibration using covariates such as identity of the reference base, identity of the mismatch base, base position within the read, and neighboring base identities can significantly improve these error rate estimates[McKenna2010]_.

:program:`breseq` uses an empirical error model that is trained by assuming that nearly all of the disagreements between mapped reads and the reference genome are due to sequencing errors and not bona fide differences in the sample from the reference genomes. It simply counts the number of times that each base or a single-base gap is observed in a read opposite each base or a single-base gap. These counts are further binned by the quality score of the read base. (The quality score of the next aligned base in the read is used for single-base deletions). A pseudocount of 1 is added to all categories to prevent zero observations. These error counts are converted to error rates by dividing the count in each cell by the sum across that quality score.

This plot shows a typical empirical error model fit to Illumina Genome Analyzer data. Note that :program:`breseq` is agnostic about exactly what quality score scheme is used by the input FASTQ. The X-axis assumes that `Sanger FASTQ format <http://en.wikipedia.org/wiki/FASTQ_format>`_ was used when converting ASCII characters to numbers, but it will not affect the results if input quality score format used a different offset, e.g. `Solexa FASTQ format <http://en.wikipedia.org/wiki/FASTQ_format>`_.

Haploid Bayesian SNP caller
***************************

At each alignment position, :program:`breseq` calculates the Bayesian posterior probability of possible sample bases given the observed read bases. Specifically, it uses a haploid model with five possible base states (A, T, C, G, and a gap), assumes a uniform prior probability of each state, and uses the empirical error model derived during base quality re-calibration to update the prior with each read base observation. 

Thus, at a given alignment position, the log10 ratio of the posterior probability that the sample has a certain base state b\ :sub:`x` versus the probability that the sample has a different base is: 

:math:`L(b_x) = \sum\limits_{i=1}^{n}\{\log_{10}[E(b_x, b_i, q_i)] - log_{10}[1 - E(b_x, b_i, q_i)]\}`

Where there are *n* reads aligned to this position, b\ :sub:`i` is the base observed in the *i*\ th read, q\ :sub:`i` is the quality of this base, and *E* is the probability of observing this read base given its quality score at a reference position with base b\ :sub:`x` according to the empirical error model.

|breseq| determines the base with the highest *L*\ , and records read alignment evidence if this base is different from the reference base. This evidence is assigned *L* as a new quality score for the alignment column and this base prediction.

Recall that :program:`breseq` will only find indels of 1 nt length and base substitutions as read alignment evidence, because all read alignments with gaps of ≥2 bases were split in a pre-processing step. Longer indels are identified from :ref:`new-junction-evidence`.

.. _unknown-base-evidence:

Unknown base (UN) evidence
--------------------------

When there is insufficient evidence to call a base at a reference position, :program"`breseq` reports this base as "unknown". Contiguous stretches of unknown bases are output and shown in the results. Explicitly marking bases as unknown is useful when analyzing many similar genomes; it allows one to ascertain when a mutation found in certain data sets may have been missed in others due to low coverage and/or poor data quality in some data sets.

.. _missing-coverage-evidence:

Missing coverage (MC) evidence
------------------------------

As :program:`breseq` traverses read alignments it predicts deletions as it encounters genomic regions missing and low coverage.

.. _read-coverage:

Coverage distribution model
***************************

If read sequences were randomly distributed across the entire reference sequence, then the number of positions with a given read coverage depth would follow a Poisson distribution. In practice, the actual read coverage depth distribution deviates from this idealized expectation in two ways:

First, it is generally overdispersed relative to a Poisson distribution, e.g., there are more positions with higher and lower coverage than expected. This may represent a bias in the steps used to prepare a DNA fragment library or sequencing differences that cause more reads originating in certain regions of the genome to fail quality filtering steps. This overdispersion occurs even when re-sequencing a known genome. In fact, there is often a fingerprint of coverage bias where specific stretches consistently have higher or lower coverage than average across different instrument runs and DNA preparations.

Second, there may be real mutations in the sequenced genome that affect the observed coverage distribution, such as large deletions and duplications. Deletions will add weight to low end of the distribution because they cause reference positions to have zero or very low coverage. Non-zero coverage commonly is present in practice because there may be a small amount of contaminating DNA from a different sample that does not have this deletion or a small number of reads with errors may spuriously map to the deleted region. Duplications and amplifications will add weight to the distribution at higher coverage values. Deletions tend to be more common that amplifications during laboratory evolution.

:program:`breseq` fits a  negative-binomial distribution (an overdispersed Poisson distribution) to the read coverage depth observed at unique-only reference positions. It uses left censored data to mitigate the effects of deleted regions on the overall fit. The threshold for censoring is determined by by first finding the read-depth with the maximum representaton in the distribution after moving-average smoothing with a window size of 5. Positions with coverage less than half this maximum read-depth are ignored during fitting.

In this example, circles represent the number of positions in the reference with a given read coverage depth. Data points that were censored during fitting are shown in red. The solid line is the least-squares best negative binomial fit to the data, and the dashed line is the best Poisson fit.

Seed and extend algorithm
*************************

From the fit coverage distribution, :program:`breseq` calibrates how it will call deletions. Deletion predictions are initiated at every position with unique-only coverage of zero. They are extended in each direction and merged until unique coverage exceeds a threshold calculated from the overall coverage distribution for the reference sequence. This cutoff is the the minimum threshold coverage *t* that satisfies the following relationship:

:math:`F(t) > 0.05\times\sqrt{L}`, 

where *F* is the negative binomial cumulative distribution function with best-fit mean and size parameters and *L* is the reference sequence length. 

In some cases there is ambiguity concerning the size of missing coverage regions because they encompass or overlap regions with degenerate read matches. Even if a specific example of a repetitive region is deleted, there will still appear to be coverage because exact copies still exist elsewhere in the genome.

:program:`breseq` assumes that any regions with degenerate coverage that occur wholly within a region of low unique coverage (defined as above) have been deleted with the flanking sequences. If a region of degenerate coverage overlaps one end of the missing region prediction, then that end is assigned a range of possible reference positions. They reflect the two extreme possibilities that (1) the entire contiguous repetitive region is missing and (2) the entire contiguous repetitive region is still there. To determine the latter boundary, the same extend algorithm and threshold used for unique-only coverage are applied to the degenerate coverage depth normalized to the number of matches each degenerate read had in the original genome.

This example shows things....

Mutational event prediction
---------------------------

The previous sections describe **evidence** for mutations. :program:`breseq` next tries to predict biological **mutational events** from this evidence.

Base substitutions
******************

``RA evidence = SNP or SUB mutation``


When the quality score of RA evidence is greater than a cutoff of log10 the total reference sequence length + 6, a base substitution mutation is called. When only single base is affected, |breseq| calls a SNP mutation. When multiple SNPS occur adjacent to each other or in conjunction with indels (see below), |breseq| calls a substitution (SUB) mutation.

Short insertions and deletions
*******************************

``RA or JC evidence = INS, DEL, or SUB mutation``

For single-base insertions and deletions, RA evidence with gap characters is used to call mutations as in the case of base substitutions. For slightly longer insertions and deletions, for which missing coverage evidence may not exist, these events may be predicted solely on the basis of new junctions joining them.

Large deletions 
*************************

``MC+JC evidence = DEL mutation``

Missing coverage typically indicates a large deletion event. When a junction also exists that will join the precise endpoints, |breseq| predicts a deletion mutation.

Mobile element insertions
*******************************

``JC+JC evidence = MOB mutation``

When two junctions exist that would join positions close by in the reference sequence to the ends of an annotated ``repeat_region``, |breseq| predicts a mobile element insertion. It further tries to shift the ends of the junctions such that they align best with the ends of the mobile element. 

Duplications
*************

``JC evidence = AMP mutation``

If new junction evidence connects a region of the genome to a region on the same strand upstream, then it typically indicates that the intervening bases have been duplicated. |breseq| currently does not predict copy number changes from coverage, so coverage should be manually examined to see if these mutations are real.

Chromosomal inversions
*************************

``JC+JC evidence = INV mutation``

If new junctions that reciprocally join two distance regions of a genome occur in unique sequences (see below), |breseq| predicts an inversion mutation.

Orphan evidence
******************

Evidence that is not assigned by |breseq| to any of the methods above is shown in a separate section of the output so that it can be manually examined.

Limitations
-----------

Even given perfect data, |breseq| cannot find some types of mutations:

`Novel sequences, not existing in reference sequences`
   Because |breseq| maps reads to the reference sequences, it will not find new sequences that have been inserted into the genome or new extrachromosomal DNA fragments such as plasmids. Reads that do not map to the reference genome are dumped to an output filesuitable for *de novo* assembly, so that they can be examined with other tools.
`Mutations in repeat regions` 
   In genomic regions where the only mapped reads also map equally well to other locations in the genome, it is not possible to call mutations. This is an inherent limitation of short-read data. These regions are reported as 'UN' evidence, so that the user can distinguish where in the genome there was not sufficient coverage of uniquely mapped reads to call mutations.
`Chromosomal inversions and rearrangements through repeat sequences`
   These types of mutations cannot be detected when they involve sequence repeats on the order of the read length. Reads that span repeats and uniquely align in the reference sequence on each end are necessary to detect them. :program:`breseq` does not use mate-paired or paired end information to identify these kinds of mutational events.
   
References
-----------

.. [McKenna2010]  McKenna et al. (2010) The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research*  **20**:1297-1303.
