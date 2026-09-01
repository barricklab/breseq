# MP: Missing pair evidence

A point where a **novel sequence** — one present in neither the reference nor any candidate junction
— has been inserted into the genome. `MP` is the only evidence type that can see sequence absent from
the reference entirely. It is experimental and off by default; enable it with
`--predict-missing-pairs`. It is not currently promoted to any mutation type, so accepted `MP` items
appear as unassigned evidence.

## At a glance

| | |
|---|---|
| Section in `summary.html` | *(none of its own; see the tables below)* |
| Tables in `summary.html` | `Missing pair (MP) evidence metrics` + `... gates` |
| Banner in `index.html` | `Unassigned missing pair evidence` |
| Sort order in the GenomeDiff file | 17 |
| Enabled by | `--predict-missing-pairs` |
| Requires | paired reads |
| Promotes to | *(nothing — reported as evidence only)* |
| Rejected items visible in HTML | yes, on `marginal.html` (top 20) |

## The signal

A sequencing library is made of fragments, each read from both ends. When a fragment straddles the
boundary between reference sequence and inserted novel sequence, one mate lands in reference
sequence, where it maps normally, and the other lands inside the insert — which is not in the
reference, so it maps nowhere at all.

The mapped mates therefore **pile up facing the insertion point**, all on the same strand, all with
missing partners. That pile-up is what `MP` detects. It marks where the insert begins; it says
nothing about what the insert is or how long it might be.

Because the signal is one-sided, a real insertion normally produces **two** `MP` items, one at each
shoulder, on opposite strands.

## What this evidence cannot see

This is precisely the case neither [`JC`](evidence-jc.md) nor [`DP`](evidence-dp.md) can reach:

- `JC` needs a read split across the breakpoint with **both** halves mapping to the reference. If one
  half is novel sequence, there is nothing for it to align to.
- `DP` needs **both** mates mapped, so that a distance and orientation can be computed. A mate in
  novel sequence has neither.
- [`SC`](evidence-sc.md) sees a partially-aligning read, which is a related but different signal — a
  clipped tail is sequence the aligner declined to place, not sequence demonstrably absent from the
  reference.

In exchange, `MP` gives up almost everything else. It does not locate the breakpoint to the base, it
does not report a size, and it cannot name the inserted sequence. It tells you *there is something
here that is not in your reference*, and that is all — which is often exactly the thing you most need
to know.

!!! note "A mate that failed to align is not the same as a mate flagged unmapped"
    `MP` counts only reads whose mate produced **no alignment at all**. It deliberately does *not*
    count every read flagged `BAM_FMUNMAP`, because that flag also covers mates the aligner placed
    and _breseq_ then rejected on `--require-match-fraction`. Those are partially-aligning reads —
    `SC`'s signal — not evidence that a mate's sequence is missing from the reference. On a short-read
    library the two populations differ by orders of magnitude, and conflating them produces a flood
    of false positives. See the worked example below.

## From reads to evidence

1. **Mark pairs.** During alignment resolution each read records whether its mate ever produced an
   alignment, in a dedicated BAM tag. This is the fact `MP` counts, distinct from `BAM_FMUNMAP`.
2. **Measure the background.** Across *every* reference column, _breseq_ measures the rate at which a
   mate fails to align, and how unevenly that rate is spread. This is the null the score is tested
   against, and it is fitted to the run rather than assumed.
3. **Seed candidate regions.** A region is opened where, within one paired-mapping-distance window,
   at least `--missing-pair-seed` reads on one strand have unmapped mates *and* they make up at
   least `--missing-pair-seed-fraction` of that strand's reads there. Both conditions are sensitivity
   filters only; neither decides the call.
4. **Place the boundary.** The insertion point is put at the last retained reference base on the
   flank the supporting reads sit on, and the strand records which side the insert lies on.
5. **Rescan and count.** Within the counting window, count the supporting reads
   (`unpaired_read_count`), every crossing-strand read on the kept flank (`window_read_count`), the
   distinct start positions (`distinct_read_count`), and the pairs whose mate maps *past* the
   insertion point (`spanning_pair_count`).
6. **Score and gate.** Compute the genome-wide score and apply the gates below.

## Metrics

The `Missing pair (MP) evidence metrics` table in `summary.html` reports these row for row. None of
them accepts or rejects anything — they are what the [gates](#gates) below are measured against.

### counting window

The width, in bases, over which reads around a candidate are counted. Derived from the library's
paired-mapping distance.

### mate-unmapped background

The genome-wide rate at which a read's mate fails to align, **measured over every reference column,
not assumed**. This is the null probability the score tests against. `--missing-pair-minimum-rate`
puts a floor under it (default 0.0001), which only binds on a library so clean that the measured rate
is essentially zero, where it would otherwise make the test infinitely sensitive. The `basis` column
says when the floor was applied.

### background dispersion

How unevenly that background rate is spread across the genome, as an over-dispersion ρ, reported
alongside a Pearson φ. A value of `none (binomial)` means the background was even enough to use a
plain binomial null. Two options shape this: `--missing-pair-dispersion-trim-frequency` excludes
windows at or above a given local mate-unmapped fraction, so that a genuine insertion cannot define
the background it is judged against; `--missing-pair-maximum-dispersion` caps the fit against
degenerating on a small or unusual reference.

The `basis` column notes that **columns overlap, so the column count is not an independent sample
size** — worth remembering before reading the dispersion as a precise quantity.

### independent chances

`2 strands × reference length / counting window`, the number of independent opportunities for a
false positive. This is the multiplier that turns a per-position probability into a genome-wide
E-value.

### candidate seed

The `--missing-pair-seed-fraction` threshold, reported as `fraction ≥ X`. A sensitivity filter that
stops a fixed count from seeding continuously at any decent coverage. **This is not what decides an
`MP` call** — the score is.

## Gates

The `Missing pair (MP) evidence gates` table lists each decision, the rule in force, and the
`reject=` value an item picks up when it fails.

### score cutoff

`--missing-pair-score-cutoff`, default 3. Rejects as `MISSING_PAIR_SCORE`. The score is minus the
log10 of the expected number of
positions anywhere in the reference where this many reads on one strand would lose their mates by
chance, given the background rate and dispersion above. 0 means one such position expected per
genome; 3 means one per thousand genomes. See
[the shared definition of `score`](evidence-overview.md#score).

**This is the test that decides an `MP` call.** Every other quantity on the item is a local sanity
check.

### supporting reads

`--missing-pair-minimum-reads`, default 3. Rejects as `MISSING_PAIR_COUNT`. Counts only reads whose
mate produced no alignment at all — see the note above.

### distinct start positions

`--missing-pair-minimum-distinct`, default 2. Rejects as `MISSING_PAIR_DUPLICATES`, so that PCR
duplicates of one molecule cannot carry a prediction.

### local frequency

`--missing-pair-frequency-cutoff`, defaulting to `--polymorphism-frequency-cutoff`. Rejects as
`MISSING_PAIR_FREQUENCY`. Applied to the lower confidence bound on unpaired / (unpaired + spanning)
reads. This says how much of the sample carries the insertion, not whether it is there.

### outcome

How many items were examined, accepted, and rejected by each gate above — including how many were
`dropped with no supporting read left in the rescan window`.

### A real pair of tables

From `long_ltee_ara_p1_50k_pe101`, a 2x101 paired library of an *E. coli* clone:

**Missing pair (MP) evidence metrics**

| metric | value | basis |
|---|---|---|
| counting window | 232 bases | one median fragment length — beyond that a read has no mate that could have reached the breakpoint |
| mate-unmapped background | 0.544% | measured over every reference column, not assumed |
| background dispersion | ρ = 0.00088 | Pearson φ = 1.08 over 9090064 columns (mean 88 reads), 8640 trimmed at a local fraction of 0.25 or above. Columns overlap, so the column count is not an independent sample size. |
| independent chances | 39912 | 2 strands x reference length / counting window |
| candidate seed | fraction ≥ 0.25 | 1 regions seeded — a sensitivity filter only |

**Missing pair (MP) evidence gates**

| gate | rule | rejects as | basis |
|---|---|---|---|
| score cutoff | ≥ 3.0 | `MISSING_PAIR_SCORE` | −log10 of the expected number of positions where this many reads would lose their mates by chance |
| supporting reads | ≥ 3 | `MISSING_PAIR_COUNT` | reads with a mate that produced no alignment at all |
| distinct start positions | ≥ 2 | `MISSING_PAIR_DUPLICATES` | so PCR duplicates of one molecule cannot carry a prediction |
| local frequency | ≥ 10.0% | `MISSING_PAIR_FREQUENCY` | exact lower confidence bound on unpaired / (unpaired + spanning) reads |
| outcome | 1 examined | | 0 accepted, 1 rejected (score) |

`mate-unmapped background` is the metric to read first: on this library about half a percent of reads
have a mate that never aligned. That is the rate every candidate is judged against, and it is
*measured*, not assumed — which is what lets the same cutoff behave sensibly on libraries whose
backgrounds differ by orders of magnitude.

`background dispersion` reports a Pearson φ very close to 1, meaning the background is nearly
binomial: mate-unmapped reads are scattered fairly evenly rather than clumped. On a library where
they clump, φ rises, the null widens, and correspondingly more supporting reads are needed to clear
the same score.

This run seeded one candidate and rejected it on score, which is the ordinary outcome for a clone
whose genome matches its reference. `MP` firing rarely is the expected behaviour, not a failure — see
the second worked example.

## Why items are rejected

| `reject` value | Shown in the report as | What it means here |
|---|---|---|
| `MISSING_PAIR_SCORE` | Missing pair score below the genome-wide false-positive cutoff. | The decisive gate. This many reads losing their mates in this window is not surprising given the run's own background rate. Usually the correct verdict on a library with a high or patchy mate-unmapped rate. |
| `MISSING_PAIR_COUNT` | Too few reads with unmapped mates support this position. | Fewer than `--missing-pair-minimum-reads` supporting reads. |
| `MISSING_PAIR_DUPLICATES` | Supporting reads with unmapped mates start at too few distinct positions. | The support comes from too few distinct molecules — probably PCR duplicates of one. |
| `MISSING_PAIR_FREQUENCY` | Missing pair local frequency below cutoff. | Of the molecules that could have contradicted the call, too small a fraction lost their mates. A statement about how much of the sample carries the insertion, not about whether it is there. |
| `NEARBY_BETTER_MISSING_PAIR` | A better-scoring missing pair prediction lies within one paired-mapping distance. | Two candidates describe the same shoulder; the weaker one is dropped. |

Items at a contig end carry `ignore=CONTIG_END` rather than a rejection, and are dropped from the
report entirely — near a sequence end the statistic is not measurable. See
[accepted, rejected, ignored](evidence-overview.md#accepted-rejected-ignored).

## GenomeDiff fields

The positional fields are specified in
[GenomeDiff File Format](genomediff-file-format.md#mp-missing-pair-evidence): `seq_id`, `position`
and `strand`.

Notable `name=value` pairs:

*   **score** *\<float>* — the genome-wide E-value score described above. **This is the test that
    decides an `MP` call.** The null it is measured against is fitted to the run itself and reported
    in the metrics table.
*   **unpaired_read_count** — supporting reads: mapped, on the crossing strand, with the flank on the
    kept side, and with a mate that aligned nowhere.
*   **window_read_count** — *every* crossing-strand read on the kept flank within the counting
    window. This is the denominator **score** uses, and the one the genome-wide null is measured
    over.
*   **spanning_pair_count** — pairs whose mate aligns *past* the insertion point, into the sequence
    the insert would occupy. Those molecules demonstrably carry reference sequence there, so they are
    evidence against the call.
*   **total_read_count** — `unpaired_read_count + spanning_pair_count`. This is the denominator
    **frequency** uses, and it is deliberately *not* `window_read_count`: a pair whose mate never
    reached the insertion point cannot speak to the call either way.
*   **frequency**, **frequency_lower**, **frequency_upper** — of the molecules that could have
    contradicted this call, the fraction that instead lost their mate, with 95% confidence bounds.
    This says how much of the sample carries the insertion; it does *not* say whether the insertion
    is there, which is **score**'s job.
*   **distinct_read_count** — distinct outer coordinates among the supporting reads, so that PCR
    duplicates of one molecule cannot carry a prediction.
*   **candidate_unpaired_count** — peak in-window count while the candidate region was open. A seed
    diagnostic; nothing is gated on it.
*   **redundant** — present when a majority of the supporting reads mapped to more than one place.

## In the HTML report

Accepted items appear on `index.html` under `Unassigned missing pair evidence`; rejected items appear
on `marginal.html` under `Marginal missing pair evidence`.

`seq id`\
Identifier for the reference sequence containing the insertion point.

`position`\
The last retained reference base on the flank the supporting reads sit on.

`direction`\
Which side of `position` the insert lies on.

`unpaired`\
`unpaired_read_count` — the supporting reads.

`distinct`\
`distinct_read_count` — how many distinct molecules those reads represent.

`spanning`\
`spanning_pair_count` — pairs arguing *against* the call.

`total`\
`unpaired + spanning`, the denominator of `freq`.

`window`\
`window_read_count` — every crossing-strand read in the window, the denominator `score` uses.

`freq`, `range`\
Local variant frequency and its 95% confidence interval. The cutoff is applied to the interval, not
to `freq`; see [`freq` and `range`](evidence-overview.md#freq-and-range).

`score`\
The genome-wide score. See [`score`](evidence-overview.md#score).

`gene`, `product`\
Annotation at the insertion point.

## Options

!!! warning "Experimental"
    `MP` prediction is experimental and off by default. Enable it with `--predict-missing-pairs`.
    It requires paired-mapping, which is the default for paired input.

`--predict-missing-pairs`\
Predict missing read-pair (MP) evidence: places where reads pile up whose mates did not map anywhere,
the signature of a novel sequence inserted into the genome.

`--missing-pair-seed <int>` (default 3)\
Minimum number of reads with unmapped mates within a paired-mapping-distance window required to seed
an MP candidate region. Sets the [candidate seed](#candidate-seed) gate.

`--missing-pair-seed-fraction <float>` (default 0.25)\
Minimum fraction of the reads on one strand within a paired-mapping-distance window whose mates did
not map, required to seed a candidate region. A sensitivity filter that keeps a fixed count from
seeding continuously at any decent coverage.

`--missing-pair-minimum-reads <int>` (default 3)\
Only accept MP evidence supported by at least this many reads with unmapped mates. Produces
`MISSING_PAIR_COUNT`.

`--missing-pair-minimum-distinct <int>` (default 2)\
Only accept MP evidence whose supporting reads start at at least this many distinct positions, so
that PCR duplicates of one molecule cannot carry a prediction. Produces `MISSING_PAIR_DUPLICATES`.

`--missing-pair-score-cutoff <float>` (default 3)\
Only accept MP evidence whose score is at or above this value. `0` = OFF. Sets the
[score cutoff](#score-cutoff) gate and produces `MISSING_PAIR_SCORE`.

`--missing-pair-minimum-rate <float>` (default 0.0001)\
Floor on the genome-wide rate at which a read's mate fails to align, used as the null the score tests
against. Only binds on a library so clean that the measured rate is essentially zero. Sets the
[mate-unmapped background](#mate-unmapped-background) gate.

`--missing-pair-dispersion-trim-frequency <float>` (default 0.25)\
Exclude windows at or above this local mate-unmapped fraction when measuring how unevenly the
background is spread, so a real insertion cannot define the background it is judged against.
`0` = OFF.

`--missing-pair-maximum-dispersion <float>` (default 0.05)\
Cap on the fitted over-dispersion of the background. Guards against a degenerate fit on a small or
unusual reference. `0` = OFF.

`--missing-pair-frequency-cutoff <float>` (default: follows `--polymorphism-frequency-cutoff`)\
Only accept MP evidence when the lower 95% confidence bound on its local variant frequency is at or
above this value. A secondary check on how much of the sample carries the insertion, not a test of
whether it is there. `0` = OFF. Produces `MISSING_PAIR_FREQUENCY`.

## Worked examples

### An accepted item

```text title="tests/long_ltee_ara_m1_40k_pe36/expected.gd"
MP	1951	.	REL606	2039143	1	candidate_unpaired_count=7	distinct_read_count=3
	frequency=1.0000	frequency_lower=0.8074	frequency_upper=1.0000	gene_name=wbbA
	gene_position=coding (552/747 nt)	score=13.0	spanning_pair_count=0
	total_read_count=14	unpaired_read_count=14	window_read_count=81
```

*(Fields trimmed to those under discussion, and wrapped for display; a real `.gd` line is one
tab-separated row.)*

Reading it: every read that could have argued against this call instead lost its mate —
`spanning_pair_count` is zero, so `total_read_count` equals `unpaired_read_count` and the frequency
pins at 1.0. The `score` is far above the default cutoff of 3, meaning a pile-up this size is
extremely unlikely to arise from this run's background rate. `distinct_read_count` confirms the
support is not PCR duplicates of a single molecule. The insertion point sits inside *wbbA*, a
glycosyltransferase in the O-antigen cluster — a region where insertions are common in the LTEE.

Note what the item does *not* say: nothing about what was inserted, or how long it is. Only that at
this base, on this strand, the reference stops describing the sample.

### A whole run that correctly predicts nothing

The clearest lesson about `MP` comes from a library where it fires zero times. The
`long_ltee_ara_m3_32k_mp2800` test is a 2×35 mate-pair library, and it predicts **no** `MP` items at
all. That is the right answer, not a loss of sensitivity, and the reason is the distinction drawn
above:

Because the reads are short, `--require-match-fraction 0.9` demands that 32 of 35 bases align. Over a
million reads in that run's BAM carry the `BAM_FMUNMAP` flag — but only about 0.6% of those have a
mate that bowtie2 truly placed nowhere. The other 99.4% are mates the aligner *did* place and
_breseq_ then rejected for aligning over too little of their length. Those are partially-aligning
reads, which is [`SC`](evidence-sc.md)'s signal, not evidence that a mate's sequence is missing from
the reference.

When `MP` counted every `BAM_FMUNMAP` read, this test produced dozens of items, and the handful that
looked most convincing — frequency 1.0, no spanning pairs — were exactly the positions where every
crossing mate happened to be rejected. A frequency of 1.0 with zero opposing pairs looks identical
whether the numerator is real or an artifact; only the genome-wide score can tell them apart, and
only if its numerator counts the right thing.

A genuinely clonal insertion costs essentially every spanning fragment its mate, which in a window
this size means thousands of reads, not a handful. Counting only never-aligned mates makes that
distinction visible, and the seed then never fires here at all.

!!! tip "Reading a surprising MP count"
    If a run predicts far more `MP` items than seems plausible, read the
    `mate-unmapped background` and `background dispersion` rows of the metrics table first. A high
    background with high dispersion means the library, not the genome, is generating the signal.

## See also

- [Evidence overview](evidence-overview.md) — how `MP` compares with `JC`, `SC`, `DP` and `PD`
- [PD: Pair distance evidence](evidence-pd.md) — the other collective-pair statistic
- [SC: Soft clipping evidence](evidence-sc.md) — partially-aligning reads, the signal `MP` excludes
