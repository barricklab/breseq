# SC: Soft clipping evidence

A position where reads stop aligning part-way along their length, and the unaligned tails **agree
with each other** on what comes next. That agreement is the signal: reads clipped for uninteresting
reasons differ from one another, while reads crossing a real breakpoint all carry the same donor
sequence. It is experimental and off by default; enable it with `--predict-soft-clipping`. It is not
currently promoted to any mutation type, so accepted `SC` items appear as unassigned evidence.

## At a glance

| | |
|---|---|
| Section in `summary.html` | *(none of its own; see the tables below)* |
| Tables in `summary.html` | `Soft clipping (SC) evidence metrics` + `... gates` |
| Banner in `index.html` | `Unassigned soft clipping evidence` |
| Sort order in the GenomeDiff file | 15 |
| Enabled by | `--predict-soft-clipping` |
| Requires | nothing (works on single-end data) |
| Promotes to | *(nothing — reported as evidence only)* |
| Rejected items visible in HTML | yes, on `marginal.html` (top 20) |

!!! warning "`--predict-soft-clipping` changes read mapping"
    Enabling `SC` lowers `--require-match-fraction` from 0.9 to 0.5 unless you set it explicitly.
    A read only produces a clipped tail if the aligner is willing to place it while leaving part
    unaligned, so `SC` needs a permissive match fraction to see anything at all. The side effect is
    that **every** evidence type in the run is computed from a slightly different set of alignments.
    Two runs that differ only in this flag are not directly comparable.

## The signal

An aligner that cannot match a read end-to-end has two choices: reject the read, or place the part
that does match and mark the rest as *soft clipped* — present in the read, absent from the alignment.

Soft clipping happens for many reasons, most of them dull: adapter read-through at the end of a short
fragment, low-quality tails, a stray chimeric molecule. What distinguishes a structural variant is
that the clipped tails **agree**. If a mobile element inserted at this base, every read crossing the
insertion point stops aligning at the same base and continues into the same donor sequence. Adapter
and quality artifacts produce tails that disagree with one another.

`SC` is essentially [`JC`](evidence-jc.md) with one half missing. `JC` needs a read whose *both*
halves align somewhere in the reference; `SC` fires when only one half does. In exchange for that
weaker requirement it gives up knowing what the other side is: the clipped tail is real sequence, but
it is not itself aligned, so `SC` reports the breakpoint base and the tail's consensus and stops
there.

## What this evidence cannot see

- **It is one-sided.** A real breakpoint produces two `SC` items, one for reads clipped on each side,
  on opposite strands.
- **It cannot say where the tail came from.** The `clipped_sequence` is reported verbatim; matching
  it back to the reference is left to you. That is a genuinely useful thing to do by hand — see the
  mate-pair example below.
- **It cannot distinguish novel sequence from displaced sequence.** A clipped tail is sequence the
  aligner declined to place. Whether it exists elsewhere in the reference, or nowhere at all, `SC`
  does not know. Sequence demonstrably absent from the reference is [`MP`](evidence-mp.md)'s signal.

## From reads to evidence

1. **Tabulate clipping.** For every reference position, count reads clipped at that base in each
   direction and on each strand, along with the reads that align straight through it. A read must
   have at least `--soft-clipping-minimum-bases` clipped to count as a clip event — and the same
   number of aligned bases on both sides of a position to count *against* clipping there.
2. **Build a tail consensus.** Among the clipped reads at a position, compute the consensus of their
   clipped tails and how many agree with it. A read counts as agreeing when at least
   `--soft-clipping-consensus-base-fraction` of the compared bases match; at the default of 0.95 over
   12 bases, that allows one mismatch.
3. **Measure the background.** Across the whole reference, measure the rate at which agreeing clip
   events occur per read opportunity, and how unevenly that rate is spread. Positions that are
   themselves heavily clipped are excluded from this estimate, so a genuine breakpoint cannot define
   the background it is judged against.
4. **Score and gate.** Compute the genome-wide E-value and apply the gates below.

## Metrics

The `Soft clipping (SC) evidence metrics` table in `summary.html` reports these row for row. None of
them accepts or rejects anything — they are what the [gates](#gates) below are measured against.

### clipped bases required

`--soft-clipping-minimum-bases`, default 12. How many bases must be clipped for the event to count.
Note this is dual-purpose: it also sets how much aligned sequence a read must have on **both** sides
of a position to count against clipping there, so raising it well beyond the read length leaves no
reads to compare against.

### clipping background

The genome-wide rate of agreeing clip events per read opportunity, **measured across the whole
reference rather than assumed**. This is the null probability the score tests against.
`--soft-clipping-minimum-rate` puts a floor under it, preventing an unusually clean library from
making every clipped position significant.

### background dispersion

The over-dispersion ρ of that clipping rate, with a Pearson φ, the number of tested positions, and
how many positions were trimmed from the estimate.
`--soft-clipping-dispersion-trim-frequency` excludes positions at or above a given clipped fraction;
`--soft-clipping-maximum-dispersion` caps the fit. `none (binomial)` means a plain binomial null was
used.

### one-strand clip events

**The diagnostic row.** The percentage of the run's clip events that had every clipped read on a
single strand. Reads clipped at a real breakpoint arrive from both strands; the dominant artifacts
are always the read's 3' end, so they appear on only one strand for a given clip direction.

A value near 100% means the run's clipping is essentially all artifact, and it is the first thing to
check when a library predicts implausibly many `SC` items.

## Gates

The `Soft clipping (SC) evidence gates` table lists each decision, the rule in force, and the
`reject=` value an item picks up when it fails.

### strand gate

`--soft-clipping-fisher-strand-p-value-cutoff`, default 0.05. Rejects as `FISHER_STRAND`. Fisher's exact test comparing the
strand split of the clipped reads against that of the reads reading through the position.

**On real data this gate, not the score, does most of the work** — see the real tables below,
where it accounts for 16 of the 19 rejections while the score accounts for one.

### clipped tail complexity

Rejects as `LOW_COMPLEXITY_TAIL`. `run < X, one base < Y` — rejects a consensus tail that is a homopolymer
(`--soft-clipping-maximum-tail-homopolymer-fraction`, default 0.66) or nearly one base even when not
in a single run (`--soft-clipping-maximum-tail-base-fraction`, default 0.75).

This is the companion to the strand gate. A dark-cycle poly-G tail agrees with itself *perfectly*, so
it sails through the consensus test; the complexity gate is what catches it at low counts where the
strand test lacks power.

### tail consensus

`--soft-clipping-consensus-fraction-cutoff`, default 0.5. Rejects as `CLIPPED_TAIL_CONSENSUS`. The
fraction of clipped reads at the position that must agree on one consensus tail. This is the gate
that expresses the type's core idea: reads clipped for uninteresting reasons disagree with each
other, while reads spanning a real breakpoint carry the same donor sequence.

### local frequency

`--soft-clipping-frequency-cutoff`, defaulting to `--polymorphism-frequency-cutoff`. Rejects as
`FREQUENCY_BELOW_CUTOFF`. Applied to the **lower confidence bound** on the clipped fraction, so an
item can be rejected at a frequency that reads as above the cutoff. See
[`freq` and `range`](evidence-overview.md#freq-and-range).

Note this is distinct from `RA`'s `FREQUENCY_CUTOFF`: soft clipping is only ever rejected for being
too *infrequent*. A high clipped fraction always means something was detected.

### score cutoff

`--soft-clipping-score-cutoff`, default 3. Rejects as `SCORE_CUTOFF`. Minus log10 of the expected
number of positions anywhere in the reference showing agreement this strong by chance, given the
background rate and dispersion.
See [the shared definition of `score`](evidence-overview.md#score).

### outcome

Items reported and accepted, and how many were rejected by each gate above.

### A real pair of tables

From `long_ltee_ara_p1_50k_pe101`, a 2x101 paired library of an *E. coli* clone:

**Soft clipping (SC) evidence metrics**

| metric | value | basis |
|---|---|---|
| clipped bases required | 12 bases | also how much aligned reference a read must have on BOTH sides of a position to count as reading through it |
| clipping background | 0.0032% | agreeing clip events over read opportunities, measured across the whole reference rather than assumed |
| background dispersion | ρ = 0.00175 | Pearson φ = 1.08 over 8969504 (position, direction) pairs (mean 49 reads), 44 trimmed at a clipped fraction of 0.25 or above |
| one-strand clip events | 99.6% | 13840 of 13901 agreeing clip events sit at positions that saw only one read strand |

**Soft clipping (SC) evidence gates**

| gate | rule | rejects as | basis |
|---|---|---|---|
| strand gate | p ≥ 0.050 | `FISHER_STRAND` | Fisher's exact test of the clipped reads' strand split against that of the reads reading through |
| clipped tail complexity | run < 0.660, one base < 0.75 | `LOW_COMPLEXITY_TAIL` | as fractions of the compared clipped bases |
| tail consensus | fraction ≥ 0.50 | `CLIPPED_TAIL_CONSENSUS` | the fraction of clipped reads agreeing on one consensus tail |
| local frequency | ≥ 10.0% | `FREQUENCY_BELOW_CUTOFF` | exact lower confidence bound on the clipped fraction, not the fraction itself |
| score cutoff | ≥ 3.0 | `SCORE_CUTOFF` | −log10 of the expected number of positions reaching this by chance |
| outcome | 21 reported | | 2 accepted, 16 rejected (strand), 1 rejected (clipped tail complexity), 1 rejected (score), 1 rejected (local frequency) |

Read the `one-strand clip events` metric first: **99.6%** of this run's clip events sat at positions
that saw only one read strand. That is a library whose clipping is essentially all end-of-read
artifact, and the `outcome` row shows the consequence — 16 of the 19 rejections are the strand gate,
against one apiece for tail complexity, score and frequency.

This is the concrete form of the claim made above: on real data the strand test, not the score, is
what separates signal from artifact. Note also how small the surviving set is. Twenty-one positions
reported out of nearly fourteen thousand clip events, and two accepted.

## Why items are rejected

An `SC` item commonly carries several reasons at once — the artifacts that trip one gate usually trip
others.

| `reject` value | Shown in the report as | What it means here |
|---|---|---|
| `FISHER_STRAND` | Biased read strand distribution supporting prediction. | The clipped reads' strand split differs from that of the reads reading through. On `SC` this is the single most informative gate: it is the signature of an end-of-read artifact. Note the decoded sentence is worded for `RA`, where the test compares variant against reference reads; here it compares *clipped* against *spanning* reads. |
| `LOW_COMPLEXITY_TAIL` | Clipped read tails are a homopolymer or nearly one base; typical of dark-cycle (poly-G) or adapter read-through, not of donor sequence. | The consensus tail carries no information. Catches what the consensus test structurally cannot. |
| `CLIPPED_TAIL_CONSENSUS` | Clipped read tails do not agree on a consensus sequence. | The tails disagree, so they are not all continuing into one donor sequence. The ordinary verdict on adapter and quality clipping. |
| `SCORE_CUTOFF` | E-value score below prediction cutoff. | Agreement this strong is not surprising given the run's own background clipping rate. Shares its constant with `RA`, but the underlying statistic is entirely different. |
| `FREQUENCY_BELOW_CUTOFF` | Lower confidence bound on frequency below cutoff threshold. | Distinct from `RA`'s `FREQUENCY_CUTOFF`: soft clipping is only ever rejected for being too *infrequent*. What is compared is the lower bound of the `range` column, so an item can be rejected at a frequency that reads as above the cutoff. See [`freq` and `range`](evidence-overview.md#freq-and-range). |
| `NEARBY_BETTER_SOFT_CLIPPING` | Stronger soft-clipping evidence in the same direction within a few bases; probably the same breakpoint. | Two candidates describe one breakpoint; the weaker is dropped. |

Items near a contig end carry `ignore=CONTIG_END` and are dropped from the report entirely.

## GenomeDiff fields

The positional fields are `seq_id`, `position` and `strand` — see
[GenomeDiff File Format](genomediff-file-format.md#sc-soft-clipping-evidence).

Notable `name=value` pairs:

*   **clipped_sequence** — the consensus of the clipped tails. The most useful field on the item:
    matching it back to the reference by hand is often what identifies the event.
*   **consensus_fraction** — the fraction of clipped reads at this position agreeing with that
    consensus. Gated by `--soft-clipping-consensus-fraction-cutoff`.
*   **agree_read_count**, **agree_read_count_forward**, **agree_read_count_reverse** — clipped reads
    agreeing with the consensus, and their strand split. **A zero in one of the two strand fields is
    the fastest way to spot an artifact by eye.**
*   **read_count** — all clipped reads at the position, agreeing or not.
*   **spanning_read_count_forward**, **spanning_read_count_reverse** — reads aligning straight
    through the position, by strand. These are the comparison group for the strand test.
*   **total_count** — clipped plus spanning reads.
*   **fisher_strand_p_value** — the strand test's p-value, compared against
    `--soft-clipping-fisher-strand-p-value-cutoff`.
*   **log10_e_value** — the genome-wide score.
*   **frequency**, **frequency_lower**, **frequency_upper** — the clipped fraction and its 95%
    confidence bounds. The cutoff is applied to `frequency_lower`.

## In the HTML report

Accepted items appear on `index.html` under `Unassigned soft clipping evidence`; rejected items on
`marginal.html` under `Marginal soft clipping evidence`, sorted from high to low score.

`seq id`\
Identifier for the reference sequence.

`position`\
The last aligned base before the clip.

`direction`\
Which side of `position` the clipped tails extend towards.

`clipped`\
`read_count` — all clipped reads here.

`agree`\
`agree_read_count` — those agreeing with the tail consensus.

`agree +/-`\
The strand split of the agreeing reads. A zero on either side means every clipped read came from one
strand — the signature of an end-of-read artifact.

`total`\
Clipped plus spanning reads.

`freq`, `range`\
The clipped fraction and its confidence interval; the cutoff is applied to the interval's lower
bound.

`score`\
`log10_e_value`. See [`score`](evidence-overview.md#score).

`strand p`\
`fisher_strand_p_value`.

`clipped seq`\
The consensus clipped tail.

`annotation`, `gene`, `product`\
Annotation at the clip position.

## Options

!!! warning "Experimental"
    `SC` prediction is experimental, off by default, and marked HIGHLY EXPERIMENTAL in the command
    line help. Enabling it also changes `--require-match-fraction`, as described at the top of this
    page.

`--predict-soft-clipping`\
Predict soft clipping (SC) evidence: positions where reads are unexpectedly soft clipped at their
ends, which may indicate unannotated structural variation.

`--soft-clipping-minimum-bases <int>` (default 12)\
Minimum clipped bases at a read end to count as a clip event, and the aligned sequence required on
both sides of a position to count against clipping there. Sets the
[clipped bases required](#clipped-bases-required) gate.

`--soft-clipping-score-cutoff <float>` (default 3)\
Log10 E-value cutoff. `0` = OFF. Produces `SCORE_CUTOFF`.

`--soft-clipping-consensus-base-fraction <float>` (default 0.95)\
Fraction of clipped bases that must match for a read to count as sharing the position's consensus
tail. Multiplied by the bases compared and rounded down, so the default allows one mismatch in 12.
`0` = OFF, count every clipped read.

`--soft-clipping-consensus-fraction-cutoff <float>` (default 0.5)\
Minimum fraction of clipped reads that must agree on the consensus for the evidence to be accepted.
`0` = OFF. Produces `CLIPPED_TAIL_CONSENSUS`.

`--soft-clipping-fisher-strand-p-value-cutoff <float>` (default 0.05)\
Reject evidence whose clipped reads are distributed across the strands differently from the reads
reading through. `0` = OFF. Sets the [strand gate](#strand-gate) and produces `FISHER_STRAND`.

`--soft-clipping-maximum-tail-homopolymer-fraction <float>` (default 0.66)\
Reject evidence whose consensus tail contains a single-base run at least this fraction of its length.
At the default `--soft-clipping-minimum-bases`, this rejects a run of 8 or more of the 12 bases
compared. `0` = OFF. Produces `LOW_COMPLEXITY_TAIL`.

`--soft-clipping-maximum-tail-base-fraction <float>` (default 0.75)\
Reject evidence whose consensus tail is at least this fraction a single base, even when that base is
not in one run. `0` = OFF.

`--soft-clipping-frequency-cutoff <float>` (default: follows `--polymorphism-frequency-cutoff`)\
Minimum clipped fraction for acceptance. `0` = OFF. Produces `FREQUENCY_BELOW_CUTOFF`.

`--soft-clipping-minimum-read-count <int>` (default 3)\
Minimum consensus-supporting clipped reads for a position to be reported at all. `0` = OFF.

`--soft-clipping-dispersion-trim-frequency <float>` (default 0.25)\
Positions at or above this clipped fraction are excluded from the background dispersion estimate, so
genuine breakpoints do not inflate the null they are tested against. `0` = OFF.

`--soft-clipping-maximum-dispersion <float>` (default 0.005)\
Upper bound on the estimated dispersion of the clipping rate. `0` = OFF, use a plain binomial null.

`--soft-clipping-minimum-rate <float>` (default 1e-5)\
Lower bound on the estimated background clipping rate. `0` = OFF.

## Worked examples

### An accepted item

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
SC	1107	.	REL606	2103918	1	agree_read_count=9	agree_read_count_forward=7
	agree_read_count_reverse=2	clipped_sequence=CCAGCCAGCCAG	consensus_fraction=0.9000
	fisher_strand_p_value=1.33583e-01	frequency=0.2647	frequency_lower=0.1456
	frequency_upper=0.4165	gene_name=ECB_01992	log10_e_value=7.4	read_count=10
	spanning_read_count_forward=11	spanning_read_count_reverse=13	total_count=34
```

*(Fields trimmed to those under discussion, and wrapped for display; a real `.gd` line is one
tab-separated row.)*

Reading it: the clipped reads arrive on **both** strands, so `fisher_strand_p_value` is
unremarkable and the strand gate passes — this is what a real breakpoint looks like, and the spanning
reads are likewise split evenly. `consensus_fraction` shows nearly all the clipped reads agree on the
same tail, and that tail is a repeating `CCAG` motif rather than a homopolymer, so it carries real
sequence information and clears the complexity gate. `log10_e_value` is well above the default cutoff
of 3.

The frequency is around a quarter, with the interval comfortably above the cutoff — a minority of
molecules at this position carry the breakpoint, which is what you expect from a population sample
rather than a clone.

### A rejected item: the canonical dark-cycle artifact

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
SC	1035	.	REL606	412573	1	agree_read_count=9	agree_read_count_forward=9
	agree_read_count_reverse=0	clipped_sequence=GGGGGGGGGGGG	consensus_fraction=0.6923
	fisher_strand_p_value=1.22813e-05	frequency=0.2143	frequency_lower=0.1166
	frequency_upper=0.3442	gene_name=apbA	log10_e_value=6.5	read_count=13
	reject=FISHER_STRAND,LOW_COMPLEXITY_TAIL	spanning_read_count_forward=5
	spanning_read_count_reverse=24	total_count=42
```

This item is rejected, and it is instructive that it is rejected **despite a `log10_e_value` above
the cutoff and a frequency above the cutoff**. Neither of those gates is what catches it.

Two things give it away. First, `agree_read_count_reverse` is zero: every single clipped read came
from the forward strand, while the reads spanning the position are mostly reverse. That asymmetry is
what `fisher_strand_p_value` measures, and it is several orders of magnitude past the 0.05 cutoff.
A real breakpoint is crossed by molecules sequenced in both directions.

Second, `clipped_sequence` is a pure poly-G run. On this instrument chemistry, a dead cycle reads as
G, so the 3' tail of a failing read becomes a string of Gs. Such tails agree with each other
*perfectly* — note `consensus_fraction` is respectable — which is exactly why the consensus test
cannot catch them and why `LOW_COMPLEXITY_TAIL` exists.

Because a poly-G tail is always the read's 3' end, it can only ever appear on one strand for a given
clip direction. The two rejection reasons are therefore two views of the same underlying artifact,
which is why they so often appear together.

### When SC measures the library instead of the genome

The `long_ltee_ara_m3_32k_mp2800` test deliberately leaves `--predict-soft-clipping` **off**, and the
reason is worth understanding before enabling `SC` on any mate-pair library.

A mate-pair protocol works by circularising a long fragment, so that two sequences originally
kilobases apart end up adjacent. That join is a genuine sequence junction, and reads crossing it are
genuinely soft clipped. With `SC` on, that run produces many thousands of items — roughly one per
kilobase of genome, most with a perfect tail consensus — and they are all *real*. Matching the
clipped tails back to the reference places them at the library's fragment size from their own
position, which is the circularisation junction, not a variant.

`SC` is behaving correctly here. The problem is that the artifact is indistinguishable from
structural variation without pairing information, which is why on such libraries the pair-based
predictors do the useful work instead. Trimming does not help: this protocol's junction is a blunt
genomic-to-genomic ligation with no adapter in it for a trimmer to find.

!!! tip "Reading a surprising SC count"
    Go straight to the `one-strand clip events` row of the metrics table. If it is near 100%, the run's
    clipping is dominated by end-of-read artifact and the strand gate is doing its job. If it is low
    but the count is still enormous, suspect library structure — a mate-pair or otherwise chimeric
    protocol — and check where the clipped tails map back to.

## See also

- [Evidence overview](evidence-overview.md) — how `SC` compares with `JC`, `DP`, `MP` and `PD`
- [JC: New junction evidence](evidence-jc.md) — the same signal with both halves aligned
- [MP: Missing pair evidence](evidence-mp.md) — sequence demonstrably absent from the reference
