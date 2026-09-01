This section describes the algorithms used by _breseq_ that are shared across evidence types,
and how the evidence it gathers is turned into mutation predictions.

The detection algorithm for each individual kind of evidence has its own page under
[Evidence Types](evidence-overview.md).

# Read mapping

_breseq_ uses [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2)
to map reads to the reference genome sequence.

_breseq_ does not use the distance constraints available in paired-end or mate-paired
libraries during read *alignment* — each read is mapped on its own merits. It does use
them afterwards as a source of evidence: [DP](evidence-dp.md), [MP](evidence-mp.md) and
[PD](evidence-pd.md) are all built from how the two mates of a pair were placed relative
to one another. Those three evidence types are experimental and off by default, so
without them a paired data set is effectively treated as single-end reads.

_breseq_ keeps track of two kinds of read alignments:

**unique read matches**\
Where a read aligns best to only one location in the reference sequence.

**repeat read matches**\
Where a read aligns equivalently to multiple locations in the reference
sequence (because the repeats are longer than the read length).

For some calculations, _breseq_ is concerned with:

**unique-only reference positions**\
Positions in the reference sequence that do not overlap any *repeat read
matches*.

# New junction evidence (JC)

Moved to [JC: New junction evidence](evidence-jc.md), which describes how junction
candidates are identified, how they are scored and accepted, and why one is rejected.

# Read alignment evidence (RA)

Moved to [RA: Read alignment evidence](evidence-ra.md), which describes the consensus
and polymorphism models, the statistical filters, and why an item is rejected.
<span id="polymorphism-prediction"></span>

# Read end trimming

The ends of alignments of short reads to a reference sequence can be
ambiguous with respect to insertion and deletion mutations. _breseq_
uses a conservative strategy to ignore these bases when calling
mutations.

_breseq_ examines the reference sequence for perfect sequence
repeats with lengths of 1-18 bases. Then for each position in the
reference it determines how many bases must be trimmed from the end of a
read beginning or ending at that position until the remaining bases are
unambiguously aligned with respect to possible mutations causing changes
in sequence repeats of these lengths. The minimum number of bases
trimmed at each end of any read is 1, because one can never
unambiguously know if another copy of that base was inserted by a
mutation.

<figure>
<img src="../images/end_trimming_example.png" class="align-center" width="450" alt="Example of alignment end trimming." /><figcaption aria-hidden="true"><strong>Example of alignment end trimming.</strong></figcaption>
</figure>

This example shows the number of bases that will be trimmed from the
left and right ends of a read if its match to the reference genome
begins or ends on that base. (Note that the strand of the genome that
the read matches makes no difference!) The green, blue, and yellow
highlight the repeats where the numbers come from for three test cases.

For green, a read with its left end aligned to this position is not
informative with respect to how many AG copies there are in the
sequenced genome. Therefore, it is only unambiguously aligned at the
bases starting CAT-, and the first four bases will be trimmed.
Similarly, a read with its right end aligned to the green position
cannot tell how many TA copies there are. It will only be unambiguously
aligned through -CTT, and its last four bases will be trimmed.

Trimming ends in this way enables more accurate mutation predictions
because reads extending into these repeats from either side, but not
completely crossing them, could otherwise be misinterpreted as evidence
*against* a mutation.

For example, consider this mutation, which involves insertion of a new
AGC at a site where there are already two AGC copies:

<figure>
<img src="../images/missed_mutation_no_trimming.png" class="align-center" width="600" alt="Indel mutation prediction aided by end trimming." /><figcaption aria-hidden="true"><strong>Indel mutation prediction aided by end trimming.</strong></figcaption>
</figure>

This image shows reads 1-6 aligned to the reference genome with and
without end trimming (lowercase letters in reads). Two reads cross the
entire AGCx2 repeat and show that a third AGC has been inserted.

Without end trimming, two reads on the top strand that do not cross the
new AGC insertion, contradict that there was any change to the sequence
here when they are aligned to the reference. With end trimming, these
bases are ignored because they are ambiguous with respect to possible
insertions, like the event that happened, or deletion of one AGC copy.

# Base quality re-calibration

In the FASTQ input files, each read base has been assigned a quality
score by the normal pipeline for a given sequencing technology. Base
quality re-calibration using covariates such as identity of the
reference base, identity of the mismatch base, base position within the
read, and neighboring base identities can significantly improve these
error rate estimates (see the [Bibliography](bibliography.md)).

_breseq_ uses an empirical error model that is trained by assuming
that nearly all of the disagreements between mapped reads and the
reference genome are due to sequencing errors and not bona fide
differences between the sample and the reference: it simply counts the
number of times that each base or a single-base gap is observed in a
read opposite each base or a single-base gap. These counts are further
binned by the quality score of the read base. (The quality score of the
next aligned base in the read is used for single-base deletions). A
pseudocount of one is added to counts in all categories, and these error
counts are converted to error rates by dividing the count in each cell
by the sum across that base quality score.

<figure>
<img src="../images/error_rates.png" class="align-center" width="600" height="400" alt="Example of re-calibrated error rates." /><figcaption aria-hidden="true"><strong>Example of re-calibrated error rates.</strong></figcaption>
</figure>

This plot shows a typical empirical error model fit to Illumina Genome
Analyzer data. Notice that the rate of single-base deletions is much
lower than the rate of any base miscall. Base qualities normally do not
give information about the rates of indel mutations, and this
re-calibration step allows _breseq_ to estimate the rates of these
sequencing errors.

Recall that _breseq_ requires input in [Sanger FASTQ
format](https://en.wikipedia.org/wiki/FASTQ_format). Therefore the
expected total error rate (*E*) at a given
quality score (*Q*) before re-calibration
is:

$E=10^{-\\frac{Q}{10}}$

# Unknown base evidence (UN)

Moved to [UN: Unknown base evidence](evidence-un.md).

# Missing coverage evidence (MC)

Moved to [MC: Missing coverage evidence](evidence-mc.md), which describes the seed-and-extend
algorithm and how ambiguous boundaries in repeats are handled.

# Read coverage distribution

If read sequences were randomly distributed across the entire reference
sequence, then the number of positions with a given depth of read
coverage would follow a Poisson distribution. In practice, the actual
read coverage depth distribution deviates from this idealized
expectation in at least two ways:

First, it is generally overdispersed relative to a Poisson distribution,
e.g., there are more positions with higher and lower coverage than
expected. This may represent a bias in the steps used to prepare a DNA
fragment library or sequencing differences that cause more reads
originating in certain regions of the genome to fail quality filtering
steps. This overdispersion occurs even when re-sequencing a known
genome. In fact, there is often a fingerprint of coverage bias where
specific stretches consistently have higher or lower coverage than
average across different instrument runs and DNA sample preps.

Second, there may be real mutations in the sample that affect the
observed coverage distribution, such as large deletions and
duplications. Deletions will add weight to the low end of the
distribution because they cause reference positions to have zero or very
low coverage. Non-zero coverage in true deletions is sometimes present
in practice because there may be a small amount of contaminating DNA
from a different sample that does not have this deletion or high error
rate reads may spuriously map there. Duplications and amplifications
will add weight to the distribution at higher coverage values.

For a normal sample, _breseq_ attempts to fit a negative binomial
distribution (an overdispersed Poisson distribution) to the read
coverage depth observed at unique-only reference positions for each
reference sequence (e.g., chromosome). It uses left censored data to
mitigate the effects of deleted regions on the overall fit. The
threshold for censoring is determined by first finding the read depth
with the maximum representaton in the distribution after smoothing using
a moving average window size of 5 bases. Positions with coverage less
than half this maximal read depth are ignored during fitting.

<figure>
<img src="../images/coverage_distribution.png" class="align-center" width="500" height="428" alt="Example of coverage distributon fit." /><figcaption aria-hidden="true"><strong>Example of coverage distributon fit.</strong></figcaption>
</figure>

In this example of real data, circles represent the number of positions
in the reference with a given depth of read coverage. Data points that
were censored during fitting are shown in red. The solid line is the
least-squares best fit of a negative binomial distribution, and the
dashed line is the best Poisson fit.

If a draft genome sequence is used as a reference, it may have short
contigs for which this distribution cannot be fit. You should use the
`-c` option in place of the `-r` option for this reference file to
notify _breseq_ that this is the case so that it will fit the
coverage distribution of all reference sequences in that input file
together (e.g., as one chromosome).

It is possible that the fitting procedure will fail for certain highly
biased data or when coverage is very low for a certain reference
sequence. For example, if you have done a pull-down of only certain
regions of a chromosome (like in exon sequencing). In this case,
_breseq_ will fall back to a rougher estimate of the coverage and
cutoffs for calling deletions or it may call the entire reference
sequence as deleted (and not call mutations in it). If you are doing
targeted sequencing, you should use the `-t` option so that _breseq_
will call mutations in these sequences no matter what coverage
distribution looks like (naturally, deletion mutations will not be
called in this case).

# Mutation prediction

The [Evidence Types](evidence-overview.md) pages describe **evidence** for mutations.
_breseq_ next tries to predict biologically relevant **mutational events** from that
evidence. These rules are summarized in each section using
[GenomeDiff](genomediff-file-format.md) abbreviations for types of mutations and evidence.

## Base substitutions

*RA evidence = SNP or SUB mutation*

Base substitution mutations are called from RA evidence. When only a
single base is affected, _breseq_ calls a base substitution (SNP)
mutation. When multiple base substitutions occur adjacent to each other
or in conjunction with indels (see below), _breseq_ calls a
substitution (SUB) mutation.

## Short insertions and deletions

*RA or JC evidence = INS, DEL, or SUB mutation*

For single-base insertions and deletions, RA evidence with gap
characters is used to call mutations as in the case of base
substitutions. For longer insertions and deletions, for which missing
coverage evidence may not exist, these events may be predicted solely on
the basis of new junctions joining them.

## Large deletions

*MC+JC evidence = DEL mutation*

Missing coverage typically indicates a large deletion event. When a
junction also exists that precisely joins compatible endpoints,
_breseq_ predicts a deletion (DEL) mutation.

*MC evidence = DEL mutation, between homologous copies*

A deletion that occurred between two near-identical copies of a sequence
leaves no junction at all: a read crossing the breakpoint aligns just as
well to either copy, so there is nothing to split. When the two copies
are an annotated `repeat_region` (an IS element, say), the annotation is
enough to place the deletion. When they are not — paralogous genes, an
rRNA operon, any unannotated repeat — _breseq_ works the deletion out
from the reads instead.

What survives such a deletion is a single hybrid copy: the left copy up
to the crossover point, the right copy after it. The offset that aligns
the two copies is also the size of the deletion, and it is recovered
from the reference sequence, seeded by the boundaries of the missing
coverage. The crossover is then located at the columns where the two
copies actually differ: on each of them, coverage has moved from one
copy to the other, and the position where that flips is the breakpoint.
Reads that tie between the two copies are counted here even though base
calling discards them, since which copy a read came from is exactly the
question being asked.

This is reported like any other large deletion, with a `between=` field
naming the two homologous features. Use
`--skip-homologous-DEL-prediction` to turn it off.

## Mobile element insertions

*JC+JC evidence = MOB mutation*

When two junctions exist that would join positions close by in the
reference sequence to the ends of an annotated `repeat_region`,
_breseq_ predicts a mobile element insertion (MOB). It further tries
to shift the ends of the junctions such that they align best with the
ends of the mobile element.

## Duplications

*JC evidence = AMP mutation*

If new junction evidence connects a region of the genome to a region
upstream on the same strand, then it typically indicates that the
intervening bases have been duplicated and _breseq_ predicts a
duplication. _breseq_ currently does not use evidence from changes
in read coverage depth to predict copy number, so coverage should be
manually examined to verify this class of mutations.

## Other evidence

"Orphan" evidence that passed scoring thresholds but is not assigned by
_breseq_ to any of the mutational events above is shown in a
separate section of the output so that it can be manually examined.
_breseq_ also displays some "marginal" evidence that fails the
established cutoffs, but stil has some support, on a separate results
page.

# Limitations

Even given perfect data, _breseq_ cannot find some types of
mutations:

**Novel sequences, not existing in the reference**\
Because _breseq_ maps reads to reference sequences, it cannot reconstruct
entirely novel sequences that have been inserted into the genome, or novel
extrachromosomal DNA fragments such as plasmids. Reads that do not map to the
reference genome are dumped to an output file suitable for de novo assembly, so
that they can be examined with other software programs.

With paired data, [MP](evidence-mp.md) evidence can at least tell you *where* such an
insertion begins, since a fragment crossing into novel sequence leaves one mate
unmappable. It reports the insertion point and nothing about the insert itself, so it
locates the problem rather than solving it.

**Mutations in repeat regions**\
In genomic regions where the only mapped reads also match equally well
to other locations in the genome, it is not possible to call mutations.
This is an inherent limitation of short-read data. These regions are
reported as [UN](evidence-un.md) evidence, so that the user can distinguish where in
the genome there was not sufficient coverage of uniquely mapped reads to
call mutations.

**Chromosomal inversions and rearrangements through repeat sequences**\
These types of mutations are difficult to detect when they involve sequence
repeats on the order of the read length, because [JC](evidence-jc.md) evidence needs
reads that span the repeat and align uniquely on each end.

With paired data this is no longer an absolute limit. [DP](evidence-dp.md) evidence
detects rearrangements from pairs whose mates are individually misplaced — including
between two different reference sequences, which is how a translocation or plasmid
integration is found — and [PD](evidence-pd.md) detects events too small to make any
individual pair unusual. Neither resolves a breakpoint to the base, so the two are
complements to `JC` rather than replacements for it. All three pair-based types are
experimental and off by default.

# Annotated bibliography

More information about the methods used by _breseq_ is available in
these publications:

-   Barrick, J.E., Yu, D.S., Yoon, S.H., Jeong, H, Oh, T.K., Schneider,
    D., Lenski, R.E., and Kim, J.F. (2009) Genome evolution and
    adaptation in a long-term experiment with *Escherichia coli*.
    *Nature* **461**:1243-1247. **Methods used by an early version of
    breseq are described in the supplemental materials.** doi:
    [10.1038/nature08480](https://doi.org/10.1038/nature08480)
-   Barrick, J.E., Lenski, R.E. (2009) Genome-wide mutational diversity
    in an evolving population of *Escherichia coli*. *Cold Spring Harb.
    Symp. Quant. Biol.* **74**:119-129. **Early description of
    polymorphism mode for single-nucleotide variants and small indels.**
    doi:
    [10.1101%2Fsqb.2009.74.018](https://doi.org/10.1101%2Fsqb.2009.74.018)
-   Deatherage, D.E., Barrick, J.E. (2014) Identification of mutations
    in laboratory-evolved microbes from next-generation sequencing data
    using *breseq*. *Methods Mol. Biol.* **1151**: 165–188. **Tutorial
    and practical guide to running breseq and interpreting the output.**
    doi:
    [10.1007/978-1-4939-0554-6_12](https://doi.org/10.1007/978-1-4939-0554-6_12)
-   Barrick, J.E., Colburn, G., Deatherage D.E., Traverse, C.C., Strand,
    M.D., Borges, J.J., Knoester, D.B., Reba, A., Meyer, A.G.(2014)
    Identifying structural variation in haploid microbial genomes from
    short-read resequencing data using *breseq*. *BMC Genomics*
    **15**:1039. **Detailed description of methods used to predict
    structural variation.** doi:
    [10.1186/1471-2164-15-1039](https://doi.org/10.1186/1471-2164-15-1039)
-   Deatherage, D.E., Traverse, C.C., Wolf, L.N., Barrick, J.E. (2015)
    Detecting rare structural variation in evolving microbial
    populations from new sequence junctions using *breseq*. *Front.
    Genet.* **5**:468. **Detailed description of methods used to predict
    polymorphic structural variation.** doi:
    [doi.org/10.3389/fgene.2014.00468](https://doi.org/10.3389/fgene.2014.00468)
