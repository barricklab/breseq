# RA: Read alignment evidence

A reference position where the pileup of aligned reads disagrees with the reference base. `RA` is the
workhorse evidence type: it finds base substitutions and small indels, and it is the source of every
`SNP`, `SUB`, and short `INS`/`DEL` mutation. It is on by default.

## At a glance

| | |
|---|---|
| Section in `summary.html` | `Read Alignment Evidence` |
| Metrics / gates tables | *(none — settings are listed in that section)* |
| Banner in `index.html` | *(none — `RA` is always attached to a mutation)* |
| Sort order in the GenomeDiff file | 10 |
| Enabled by | on by default; turn off with `--no-read-alignment-prediction` |
| Requires | nothing |
| Promotes to | `SNP`, `SUB`, `INS`, `DEL` |
| Rejected items visible in HTML | yes, on `marginal.html` (top 20) |

## The signal

At each reference position, some number of reads are aligned, each contributing an observed base and
a quality score for it. If the sample really carries a different base there, most reads will say so.
If the sample matches the reference, the occasional disagreement is a sequencing error.

Distinguishing the two is not a matter of counting. A base observed in three reads at high quality is
better evidence than the same base in six reads at low quality, and the error rate itself varies by
instrument, by cycle, by base and by strand. So `RA` works from an **empirical error model**
calibrated on the run's own data, and asks a probabilistic question rather than a threshold one.

Because alignments with indels longer than 2 bases are split during pre-processing (see
[`JC`](evidence-jc.md)), `RA` typically only finds indels of at most 2 bases. Longer ones are
junction evidence.

## Two modes

`RA` operates in one of two fundamentally different modes.

**CONSENSUS mode** (the default) is appropriate for re-sequencing a clonal haploid genome. It expects
every variant allele to be present in 100% of the sample.

**POLYMORPHISM (metagenomic) mode**, enabled with `-p`/`--polymorphism-prediction`, analyses a mixed
population evolved from a common ancestor. It predicts variants at frequencies between 0% and 100%
where a mixture model is well supported.

!!! warning
    Polymorphism prediction is prone to false positives. There are many biases in NGS data, and since
    _breseq_ analyses one sample at a time it cannot account for all of them. Expect to curate the
    output, and expect the statistical filters below to matter far more than they do in consensus
    mode.

Both modes compute both a consensus score and a polymorphism score for every alignment column; what
differs is how those scores are used to decide.

## From reads to evidence

### Base quality re-calibration

Before any calling, _breseq_ builds an **empirical error model** from the data itself: for each
combination of reference base, read base, and quality score, how often does that observation occur at
positions where the sample almost certainly matches the reference? This gives an error probability
*E(b\_x, b\_i, q\_i)* that reflects the actual instrument and library rather than the nominal meaning
of the quality score. Error rates are re-calibrated separately for each input FASTQ file, which is
why independently generated data sets should be supplied in separate files.

### Consensus score (Bayesian SNP caller)

At each alignment position, _breseq_ calculates the Bayesian posterior probability of each possible
sample base given the observed read bases. It uses a haploid model with five states (A, T, C, G, and
a gap), a uniform prior, and the empirical error model to update that prior with each observation.

The log10 ratio of the posterior probability that the sample has base *b\_x* versus a different base
is:

$L(b_x) = \sum\limits_{i=1}^{n}\{\log_{10}[E(b_x, b_i, q_i)] - \log_{10}[1 - E(b_x, b_i, q_i)]\}$

where there are *n* reads aligned to this position, *b\_i* is the base observed in the *i*th read,
and *q\_i* its quality.

_breseq_ takes the base with the highest *L*, and records `RA` evidence if it differs from the
reference base. The evidence is assigned log10 *L* minus the log10 of the cumulative length of all
reference sequences, as a **consensus E-value score**. Subtracting the genome length is what makes
the score an expectation *per genome* rather than per position — the same convention the other
evidence types' scores use. See [`score`](evidence-overview.md#score).

### Polymorphism score (mixed allele model)<span id="polymorphism-prediction"></span>

Next, _breseq_ tests whether the reads support a **mixture** of a major and minor variant, against
the hypothesis that all disagreements are sequencing errors.

It computes the chance of generating the observed alignment under each hypothesis that the sample is
100% of each of the four bases or a gap, takes the two highest-probability states, and tests a
mixture model allowing them at any intermediate fraction. How the mixture is specified depends on the
mode:

1.  In **consensus** mode, only the raw frequency implied by the read counts of the major and variant
    alleles is tested.
2.  In **polymorphism** mode, the maximum-likelihood allele frequency is found to a precision of
    0.000001, taking into account the observed bases *and* their quality scores.

The one-allele and two-allele models are then compared by a likelihood-ratio test: twice the natural
logarithm of the ratio of their probabilities, against a chi-squared distribution with 1 degree of
freedom. As with the consensus score, the resulting p-value is converted to a **polymorphism E-value
score** by multiplying by the total number of reference positions.

## Gates

`RA` prints no gates table. Its thresholds are listed as an `option` / `value` table in the
`Read Alignment Evidence` section of `summary.html`, which reports the `Mode`, the `Ploidy`, and then
each cutoff below — with `OFF` where a cutoff is disabled. The consensus and polymorphism cutoffs are
listed separately because they are independent settings.

### Score cutoffs

`--consensus-score-cutoff` and `--polymorphism-score-cutoff` set the E-value each score must reach.
Failing either produces `SCORE_CUTOFF`.

### Frequency cutoffs

`--consensus-frequency-cutoff` and `--polymorphism-frequency-cutoff` are applied to a **95% confidence
bound** on the variant allele frequency, not to the point estimate. Which bound depends on what is
being asked, because the two modes assume opposite defaults:

- In **consensus** mode the **lower** bound is tested — is the variant confidently the majority
  allele?
- In **polymorphism** mode the **upper** bound is tested for consensus calls — can a fixed
  interpretation still not be ruled out?

The interval comes from the fitted allele model, so it widens for low coverage or poor base quality
rather than tracking read count alone. A variant is therefore never rejected merely for shallow
coverage, only for being confidently below the cutoff. See
[`freq` and `range`](evidence-overview.md#freq-and-range).

### Coverage requirements

Four cutoffs per mode, all defaulting to 0 (off): minimum variant coverage, minimum total coverage,
and each of those per strand. They produce `VARIANT_COVERAGE`, `TOTAL_COVERAGE`,
`VARIANT_STRAND_COVERAGE` and `TOTAL_STRAND_COVERAGE` respectively.

### Strand bias

Fisher's exact test for the hypothesis that the top/bottom strand distribution of reads supporting
the major base differs from that of reads supporting the minor base. Compared against
`--polymorphism-strand-bias-cutoff` (default 0.05); failing produces `FISHER_STRAND`.

A significant result may indicate a sequencing-error hotspot on one strand generating a false
positive. This happens frequently in real data.

In practice most problem predictions of this kind have zero or a handful of reads on one strand and
many on the other supporting the minor variant. The test can fail to reject a false positive when
coverage of the minor variant is low enough that even an entirely one-strand distribution is not
significant; `--polymorphism-minimum-variant-coverage-each-strand` deals with that case more directly
by requiring at least one supporting read on each strand.

Conversely, at high coverage there may be so many observations that a statistically significant bias
is detected simply because library prep is slightly more efficient on one strand in a given sequence
context, even with good coverage of all strand/base combinations. **Use this option with caution above
about 1000 reads.**

### Quality score bias

A one-sided Kolmogorov–Smirnov test for whether base qualities supporting the minor variant are
suspiciously lower than those supporting the major variant. Compared against
`--polymorphism-quality-bias-cutoff`; failing produces `KS_BASE_QUALITY`.

**This test is OFF by default**, because it also rejects strand-balanced, well-supported calls.

### Homopolymer filters

Applying the error model per column over-predicts indel polymorphisms in homopolymer stretches. If
there are 10 A's in a row and any one is deleted, the gap aligns to the rightmost possible position,
so every one of the ten deletions looks like the same mutation — making the observed rate ten times
what the per-column error model expects. A few reads can then reach significance by the
likelihood-ratio test. Similar logic applies to insertions.

`--polymorphism-reject-indel-homopolymer-length` (and the consensus equivalent) filters these,
producing `INDEL_HOMOPOLYMER`. A value of 5 gives reasonable results for *E. coli*. These false
predictions generally also have very low frequencies (<2%) for the minor indel variant.

`--polymorphism-reject-surrounding-homopolymer-length` rejects base substitutions that *create* a
homopolymer, producing `SURROUNDING_HOMOPOLYMER`. The homopolymer must begin and end after the
changed base: with a setting of 5, TATTT→TTTTT is rejected but ATTTT→TTTTT is not.

### Flowcharts

<figure>
<img src="../images/mutation_calling_settings.png" class="align-center" width="800" height="448" />
</figure>

<figure>
<img src="../images/consensus_mode_RA_flowchart.png" class="align-center" width="600" height="572" />
</figure>

<figure>
<img src="../images/polymorphism_mode_RA_flowchart.png" class="align-center" width="800" height="525" />
</figure>

## Why items are rejected

`RA` is unusual in carrying **two** reject fields rather than one, because each item is judged twice —
once as a candidate consensus mutation and once as a candidate polymorphism:

- `consensus_reject=` — why it was not accepted as a consensus mutation.
- `polymorphism_reject=` — why it was not accepted as a polymorphism.

An item can be rejected as one and accepted as the other; the `prediction` field records which
interpretation won. In the HTML the two render as *Rejected as consensus:* and
*Rejected as polymorphism:* respectively.

| `reject` value | Shown in the report as | What it means here |
|---|---|---|
| `SCORE_CUTOFF` | E-value score below prediction cutoff. | The consensus or polymorphism E-value did not reach its cutoff. Shares its constant with `SC`, but the statistic is entirely different — here it is a Bayesian model comparison, not a background rate. |
| `FREQUENCY_CUTOFF` | Frequency below/above cutoff threshold. | The tested confidence bound fell the wrong side of the cutoff. Note "below/above": which direction depends on the mode and on which bound is being tested. |
| `FISHER_STRAND` | Biased read strand distribution supporting prediction. | The variant reads' strand split differs from the reference reads'. Typically a one-strand sequencing-error hotspot. |
| `KS_BASE_QUALITY` | Biased base quality scores supporting prediction. | The variant reads have systematically lower base qualities than the reference reads. Off by default. |
| `VARIANT_COVERAGE` | Variant not supported by required number of total reads. | Fewer than the required reads support the variant. |
| `TOTAL_COVERAGE` | Genome position does not have required minimum number of aligned reads. | The position is too shallowly covered overall. |
| `VARIANT_STRAND_COVERAGE` | Variant not supported by required number of reads on each strand. | The variant lacks support on one strand. Often more useful than `FISHER_STRAND` at low coverage. |
| `TOTAL_STRAND_COVERAGE` | Genome position does not have required minimum number of aligned reads on each strand. | The position is too shallowly covered on one strand. |
| `INDEL_HOMOPOLYMER` | Polymorphic indel expands or contracts a homopolymer stretch. | Probable homopolymer slippage rather than a real indel. |
| `SURROUNDING_HOMOPOLYMER` | Polymorphic base substitution creates a homopolymer stretch. | The substitution would create a homopolymer run, a known error mode. |
| `POLYMORPHIC_INDEL` | Indel polymorphism suppressed by --polymorphism-no-indels. | Set when `--polymorphism-no-indels` is in force. |

Items near a contig end carry `ignore=CONTIG_END` and are dropped from the report entirely.

## GenomeDiff fields

Positional fields: `seq_id`, `position`, `insert_position`, `ref_base`, `new_base` — see
[GenomeDiff File Format](genomediff-file-format.md#ra-read-alignment-evidence). `insert_position` is
the number of bases inserted after the reference position to reach this base; 0 refers to the
reference base itself.

Notable `name=value` pairs:

*   **prediction** — `consensus`, `polymorphism`, or `unknown`. Which interpretation the item was
    finally given.
*   **score** — the E-value score for the winning interpretation.
*   **consensus_score**, **polymorphism_score** — the two scores, where both were computed.
*   **major_base**, **minor_base** — the two most probable alleles.
*   **major_cov**, **minor_cov**, **new_cov**, **ref_cov**, **total_cov** — read counts as
    `forward/reverse`. The strand split is written into every one of these, which is what makes a
    one-strand artifact visible at a glance.
*   **major_frequency**, **frequency**, **frequency_lower**, **frequency_upper** — the allele
    frequency and the confidence interval the cutoffs are applied to.
*   **fisher_strand_p_value**, **ks_quality_p_value** — the two bias tests.
*   **bias_e_value**, **bias_p_value** — the combined bias statistics.
*   **consensus_reject**, **polymorphism_reject** — the two reject fields described above.
*   **snp_type** — `synonymous`, `nonsynonymous`, `nonsense`, `intergenic`, `pseudogene`, or
    `noncoding`, for a substitution in an annotated feature.

## In the HTML report

Accepted `RA` items appear attached to the `SNP`, `SUB`, `INS` or `DEL` mutation they support;
rejected ones appear on `marginal.html` under `Marginal read alignment evidence`.

`* link`\
Links to a results page showing the alignment of reads to this position.

`seq id`\
Identifier for the reference sequence where the change is located.

`position`\
Position in the reference sequence of the substitution, insertion, or deletion. It consists of two
parts: the reference position, and the insert position within it.

`ref`, `new`\
The reference base and the new base supported by the evidence.

`freq`\
The variant allele frequency.

`range`\
Confidence limits on `freq`, from the fitted allele model, so they widen for low coverage or poor
base quality rather than tracking read count alone. **This is the interval the frequency cutoffs are
applied to** — not the point estimate in `freq`. See
[`freq` and `range`](evidence-overview.md#freq-and-range).

`score (cons/poly)`\
The consensus and polymorphism E-value scores. The base-10 logarithm ratio of the posterior
probability that this position is the called base to the probability that it is any other, minus the
log10 of the total number of positions in all reference sequences. Higher means more evidence.

`reads`\
The number of reads overlapping the mutation. Note that unaligned portions of reads (lowercase bases
on a white background), read ends trimmed because their alignments may be ambiguous (lowercase bases
on a coloured background), and positions with very low base quality (highlighted yellow) are **not**
counted in this coverage number.

`annotation, gene, product`\
Description of the change's effects. The format is the same as in
[Mutation Display](output.md#mutation-display).

<figure>
<img src="../images/ra_1.png" width="750" />
</figure>

Partial alignment of reads showing that most support a base substitution. The `>` and `<` for each
named read indicate the strand of the reference sequence that it matched (top and bottom
respectively).

## Options

`RA` is always on. Its options divide into consensus and polymorphism groups, listed in
`breseq --help` under `Consensus Read Alignment (RA) Evidence Options` and
`Polymorphism Read Alignment (RA) Evidence Options`.

`-p, --polymorphism-prediction`\
Predict polymorphic mutations. Add this when analysing mixed population (metagenomic) samples.

`--consensus-score-cutoff <float>` (default 10)\
Log10 E-value cutoff for consensus base substitutions and small indels.

`--consensus-frequency-cutoff <float>` (default: consensus mode 0.50, polymorphism mode 0.95)\
Only predict consensus mutations when a 95% confidence bound on the variant allele frequency is at or
above this value. See [Frequency cutoffs](#frequency-cutoffs) for which bound is tested.

`--polymorphism-score-cutoff <float>` (default: consensus mode 10, polymorphism mode 2)\
Log10 E-value cutoff for the test of polymorphism versus no polymorphism.

`--polymorphism-frequency-cutoff <float>` (default: consensus mode 0.10, polymorphism mode 0.05)\
Only predict polymorphisms when the **lower** 95% confidence bound on the minor variant allele
frequency is at or above this value. Variants that cannot be shown to lie below
`--consensus-frequency-cutoff` are predicted as consensus mutations instead.

`--polymorphism-strand-bias-cutoff <float>` (default 0.05)\
P-value criterion for Fisher's exact test for strand bias. `0` = OFF. Produces `FISHER_STRAND`.

`--polymorphism-quality-bias-cutoff <float>` (default OFF)\
P-value criterion for the K–S test for base-quality bias. Produces `KS_BASE_QUALITY`.

`--polymorphism-bias-cutoff <float>`\
Shorthand setting both bias cutoffs. Either specific option overrides it.

`--polymorphism-reject-indel-homopolymer-length <int>` (default OFF)\
Reject indel polymorphisms that could result from expansion or contraction of homopolymer repeats of
this length or greater. Produces `INDEL_HOMOPOLYMER`.

`--polymorphism-reject-surrounding-homopolymer-length <int>` (default OFF)\
Reject polymorphic base substitutions that create a homopolymer of this length or more. Produces
`SURROUNDING_HOMOPOLYMER`.

`--polymorphism-no-indels`\
Do not predict small insertion/deletion polymorphisms from read alignment or new junction evidence.
Produces `POLYMORPHIC_INDEL`.

Each mode also has `--{consensus,polymorphism}-minimum-{variant,total}-coverage` and their
`-each-strand` variants, all defaulting to 0 except
`--polymorphism-minimum-variant-coverage-each-strand`, which defaults to 2.

## Worked examples

### An accepted consensus substitution

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
RA	109	.	REL606	2450	0	G	T	frequency=1.000e+00	gene_name=thrA
	major_base=T	major_cov=33/26	major_frequency=1.000e+00	new_cov=33/26
	prediction=consensus	ref_cov=0/0	score=222.9	snp_type=synonymous	total_cov=33/26
```

*(Fields trimmed and wrapped for display; a real `.gd` line is one tab-separated row.)*

This is what an unambiguous clonal substitution looks like. `insert_position` is 0, so it is a
substitution at the reference base rather than an insertion. `ref_cov=0/0` — **not one read** still
carries the reference G — while `new_cov=33/26` shows 59 reads supporting T, split sensibly between
the two strands. The frequency is 1.0 and `prediction=consensus`.

The `score` of 222.9 is an E-value: given the empirical error model, seeing this pileup at a position
that really carried G would be expected once in 10²²² genomes. The default consensus cutoff is 10.

`snp_type=synonymous` tells you this is a silent change in *thrA*, which matters biologically but not
to the calling.

### An item rejected as consensus, kept as a polymorphism

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
RA	142	.	REL606	2818186	1	.	C	consensus_reject=FREQUENCY_CUTOFF
	fisher_strand_p_value=4.72841e-01	frequency=1.765e-01	frequency_lower=1.007e-01
	frequency_upper=2.747e-01	ks_quality_p_value=7.62455e-01	major_base=.
	major_cov=21/21	major_frequency=8.235e-01	minor_base=C	minor_cov=3/6
	new_cov=3/6	prediction=polymorphism	ref_cov=21/21	score=30.5	total_cov=24/27
```

This item shows the two-field rejection scheme at work. It carries `consensus_reject=FREQUENCY_CUTOFF`
and **no** `polymorphism_reject`, with `prediction=polymorphism`: it was rejected as a consensus
mutation and accepted as a polymorphism.

Why. `insert_position` is 1, so this is an inserted C after position 2818186, and `major_base=.`
means the majority of reads show no insertion at all. The frequency is about 18%, and
`frequency_upper` reaches only 0.27 — far below the 0.50 consensus frequency cutoff. So the confidence
interval rules out a fixed interpretation, and `FREQUENCY_CUTOFF` is the correct verdict for the
*consensus* question.

The polymorphism question is different, and this item answers it well. Both bias tests are
unremarkable: `fisher_strand_p_value` is around 0.47, nowhere near the 0.05 cutoff, and the strand
splits bear that out — `minor_cov=3/6` has the variant on both strands, and `major_cov=21/21` is
perfectly balanced. `ks_quality_p_value` shows no quality bias either. With `score=30.5`, this is a
genuine minority allele rather than a one-strand error hotspot.

!!! tip "Read the strand splits first"
    Every `RA` coverage field is written `forward/reverse`. A minor allele whose `minor_cov` has a
    zero on one side is the classic false positive, and it is visible before any p-value. This is
    also why `--polymorphism-minimum-variant-coverage-each-strand` defaults to 2 while the other
    coverage cutoffs default to 0.

## See also

- [Evidence overview](evidence-overview.md) — how the nine types divide the space
- [JC: New junction evidence](evidence-jc.md) — where indels longer than 2 bases go
- [Base substitutions](methods.md#base-substitutions) — how `RA` becomes a `SNP` or `SUB`
- [Short insertions and deletions](methods.md#short-insertions-and-deletions)
