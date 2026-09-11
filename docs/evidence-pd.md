# PD: Pair distance evidence

A point where the read pairs whose unsequenced middle gap spans it map at a systematically different
distance from the one the sequencing library predicts. Pairs mapping **farther** apart than expected
mean the sample is missing reference sequence there; pairs mapping **closer together** mean sequence
was added. It is predicted by default on paired-end data; turn it off with
`--no-pair-distance-prediction`. It is not currently promoted to any mutation type, so accepted `PD`
items appear as unassigned evidence.

## At a glance

| | |
|---|---|
| Section in `summary.html` | *(none of its own; see the tables below)* |
| Tables in `summary.html` | `Pair distance (PD) evidence metrics` + `... gates` |
| Banner in `index.html` | `Unassigned pair distance evidence` |
| Sort order in the GenomeDiff file | 18 |
| Enabled by | on by default; turn off with `--no-pair-distance-prediction` |
| Requires | paired reads, within a single reference sequence |
| Promotes to | *(nothing — reported as evidence only)* |
| Rejected items visible in HTML | yes, on `marginal.html` (top 20) |

## The signal

A paired library has a characteristic fragment size: the two reads of a pair are sequenced from
opposite ends of one molecule, and the distance between where they map is a draw from that
distribution. Insert an extra kilobase into the sample between the two reads and the pair now maps
*closer* together than the molecule really was, because the inserted bases are not in the reference
to be counted. Delete a kilobase and the pair maps *farther* apart.

`PD` tests this collectively. For each position it asks whether the population of pairs whose
unsequenced middle spans that point is shifted, as a group, away from the library's distribution.

**No individual pair has to be unusual.** That is the whole point, and it is what separates `PD` from
its neighbours:

- [`JC`](evidence-jc.md) needs a read split across the breakpoint with both halves mapping.
- [`DP`](evidence-dp.md) tests each pair *on its own* — wrong orientation, or a distance above the
  discordant cutoff — and applies no cutoff at all on the short side. Events of a few hundred bases
  are therefore invisible to it in **both** directions.

`PD` covers exactly that gap. On a library with a few-hundred-base fragment size, a 300 bp deletion
shifts every spanning pair by 300 bases — far too little to make any one of them an outlier, and
easily enough to be obvious in the aggregate.

## What this evidence cannot see

`PD` is a statistic about a *distribution*, and it pays for that in resolution:

- **It cannot place a breakpoint to the base.** `position_range` is typically tens of bases wide. The
  exception is when the call snaps onto a validated split-read junction inside its interval, in which
  case the coordinates are exact and the item carries `snapped_to_junction`.
- **It cannot identify inserted sequence.** For an insertion it reports the *length* of what was
  added, never its identity. Where the length is consistent with an annotated repeat family the item
  may carry a `repeat_name`, but that is an inference from size, not an observation — and when two
  families are closer together than the estimator's scatter, the item honestly carries a
  `repeat_size_candidates` list instead of picking one.
- **It is single-sequence by construction.** The statistic *is* the insert distribution, which is
  undefined across two reference sequences, so `PD` only considers pairs whose mates map to the same
  sequence. Cross-sequence events are [`DP`](evidence-dp.md)'s business.

Because `PD` cannot resolve a breakpoint to the base, a lone `PD` item is not enough to promote an
insertion to a `MOB` mutation, which requires base-pair placement. Insertions that `PD` alone can see
therefore remain evidence, annotated with the repeat family they are consistent with.

!!! note "PD and DP describing the same breakpoint"
    Where a `PD` and a `DP` item describe the same breakpoint, the `DP` is removed. `PD` uses the
    whole pair population where `DP` uses only its tail, so it is the better-supported statement
    about the same event.

## From reads to evidence

1. **Fit the library.** Measure the paired-mapping distance distribution — its median, spread and
   orientation — from concordant pairs genome-wide. This is reported in the
   `Paired-End Mapping Distance Information` section of `summary.html`.
2. **Exclude multi-mapped pairs.** Pairs whose distance was *selected* rather than measured — where a
   mate mapped equally well to several repeat copies and one was picked by a per-locus vote — are
   excluded. Their distance error is coherent across every pair at that locus, which looks exactly
   like a real collective shift. The `multi-mapped pairs excluded` metric reports what this cost.
3. **Seed candidate regions.** At each position compute a rank-sum statistic over all pairs whose gap
   covers it: under the null each covering pair's distance quantile is uniform, so the mean quantile
   is a sensitive detector of a shift too small to make any single pair an outlier. Regions where
   |z| exceeds the seed threshold are opened.
4. **Estimate the shift.** Within a region, estimate `size_shift` by profile likelihood, together
   with an interval. Positive means reference sequence is missing from the sample; negative means
   sequence was added.
5. **Place the sides.** For a deletion the two sides bracket exactly the reference bases removed. For
   an insertion the sides are adjacent, because `PD` cannot place inserted sequence in the reference.
   If a validated junction lies inside the interval, snap to it.
6. **Score and gate.** Compute the genome-wide E-value and apply the gates below.

## Metrics

The `Pair distance (PD) evidence metrics` table in `summary.html` reports these row for row. None of
them accepts or rejects anything — they are what the [gates](#gates) below are measured against.

### paired-mapping distance

The library's fitted distance distribution. If the `basis` column notes that *fragments are shorter
than two reads*, only the minority of pairs that have a gap at all can carry `PD` evidence — the rest
overlap, leaving no unsequenced middle for a breakpoint to fall in.

### independent chances

`reference length / mean covering gap`, where the mean covering gap is the distance over which the
set of pairs covering a position turns over. This is the number of genuinely independent tests, and
it is smaller than the reference length because neighbouring positions are covered by largely the
same pairs.

### multi-mapped pairs excluded

The percentage of pairs excluded at step 2 above. A high value means a repeat-rich reference, a long
insert, or both. This gate exists because ambiguously-placed pairs were a documented source of false
`PD` calls, particularly in rRNA operons where a per-locus cluster vote makes the distance error
coherent across every pair at the locus.

### candidate seed

`|z| ≥ X` for the rank-sum statistic, together with how many regions were seeded. Set by
`--pair-distance-seed-z`; at the default of 0 it is *derived from the reference length* so that the
genome-wide false-seed rate stays near 1%. A sensitivity filter only — what decides a `PD` call is
the score.

### fitted null model

The fitted σ of the seed statistic's null tail, the range it was fitted over, and an IQR inflation
factor. This is the null the score is measured against, fitted to this run's own candidate regions
rather than assumed.

## Gates

The `Pair distance (PD) evidence gates` table lists each decision, the rule in force, and the
`reject=` value an item picks up when it fails.

### score cutoff

`--pair-distance-score-cutoff`, default 3. Rejects as `PAIR_DISTANCE_SCORE`. Minus log10 of the
expected number of `PD` regions anywhere in this reference that would reach this seed |z| by chance.
See
[the shared definition of `score`](evidence-overview.md#score). A score below 0 means the region is
expected to occur by chance somewhere in the genome and it is discarded outright; a score below the
cutoff but above 0 is kept as marginal evidence.

### size shift

Rejects as `PAIR_DISTANCE_SIZE`. The profile-likelihood interval on `size_shift` must not span zero,
or the direction of the event — sequence added or removed — is undetermined.

### supporting pairs

`--pair-distance-minimum-pairs`, default 3. Rejects as `PAIR_DISTANCE_COUNT`.

### distinct fragment ends

`--pair-distance-minimum-distinct`, default 2. Rejects as `PAIR_DISTANCE_DUPLICATES`, so that PCR
duplicates of one molecule cannot carry a prediction.

### local frequency

`--pair-distance-frequency-cutoff`, defaulting to `--polymorphism-frequency-cutoff`. Rejects as
`PAIR_DISTANCE_FREQUENCY`. Applied to the lower confidence bound on shifted / (shifted + normal)
pairs covering the point.

### outcome

Items accepted, rejected by each gate above, and `dropped as expected by chance somewhere in this
genome`.

### A real pair of tables

From `long_ltee_ara_p1_50k_pe101`, a 2x101 paired library of an *E. coli* clone:

**Pair distance (PD) evidence metrics**

| metric | value | basis |
|---|---|---|
| paired-mapping distance | 232 median, FR | read length 101 |
| independent chances | 22548 | reference length / mean covering gap (205 bases), the distance over which the set of pairs covering a position turns over |
| multi-mapped pairs excluded | 0.78% | 11078 of 1425305 pairs had at least one mate with more than one placement, so their distance was selected rather than measured |
| candidate seed | \|z\| ≥ 3.00 | 1403 regions seeded — a sensitivity filter only, set low so these regions are almost all noise and can serve as the null |
| fitted null model | σ = 1.15 | fitted to the number of independent regions over \|z\| = 2.25 to 4.50; the interquartile spread of z reads 1.08x, which describes the bulk and not the tail |

**Pair distance (PD) evidence gates**

| gate | rule | rejects as | basis |
|---|---|---|---|
| score cutoff | ≥ 3.0 | `PAIR_DISTANCE_SCORE` | −log10 of the expected number of PD regions reaching this seed \|z\| by chance |
| size shift | interval excludes 0 | `PAIR_DISTANCE_SIZE` | the profile-likelihood interval must not span zero |
| supporting pairs | ≥ 3 | `PAIR_DISTANCE_COUNT` | shifted read pairs required |
| distinct fragment ends | ≥ 2 | `PAIR_DISTANCE_DUPLICATES` | so PCR duplicates of one molecule cannot carry a prediction |
| local frequency | ≥ 10.0% | `PAIR_DISTANCE_FREQUENCY` | exact lower confidence bound on shifted / (shifted + normal) pairs |
| outcome | 195 examined | | 2 accepted, 193 dropped as expected by chance somewhere in this genome |

The `candidate seed` metric and the `outcome` row together show how `PD` is meant to work, and why
the numbers look alarming until you read them properly. **1403 regions were seeded and 2 were
accepted.** That is not a 99.9% failure rate; the seed is deliberately set low precisely *so that*
the regions it opens are almost all noise. They are the sample from which the null is fitted — the
`fitted null model` metric is that fit. A seed tuned to admit only real events would leave nothing to
calibrate against.

This is exactly why the seed sits among the metrics rather than the gates. It shapes the null; it
does not reject anything, and no item ever carries a rejection naming it.

`multi-mapped pairs excluded` is the metric to watch on a repeat-rich reference. Here it is under 1%.
On a long-insert library against the same genome it runs several times higher, and every excluded
pair is one whose distance error would have been coherent with its neighbours — indistinguishable, to
this statistic, from a real collective shift.

## Why items are rejected

| `reject` value | Shown in the report as | What it means here |
|---|---|---|
| `PAIR_DISTANCE_SCORE` | Pair distance score below the genome-wide false-positive cutoff. | The decisive gate. A collective shift this size is not surprising given how the seed statistic behaves across this run's own candidate regions. |
| `PAIR_DISTANCE_SIZE` | Size of the mapping distance shift is not distinguishable from zero. | The profile-likelihood interval on `size_shift` includes zero, so the direction of the shift — sequence added or removed — cannot be determined. |
| `PAIR_DISTANCE_COUNT` | Too few read pairs with a shifted mapping distance support this position. | Fewer than `--pair-distance-minimum-pairs` shifted pairs. |
| `PAIR_DISTANCE_DUPLICATES` | Supporting read pairs span too few distinct fragment ends. | The support comes from too few distinct molecules. |
| `PAIR_DISTANCE_FREQUENCY` | Pair distance local frequency below cutoff. | Of the pairs that decisively support or oppose the call, too small a fraction support it. |
| `PAIR_DISTANCE_INCONSISTENT` | Supporting read pairs cannot all be spanning an event of the estimated size at one position. | The geometry does not close: no single event of the estimated size at a single position could have produced this set of pairs. |
| `NEARBY_BETTER_PAIR_DISTANCE` | A better-supported pair distance prediction lies within one paired-mapping distance. | Two candidates describe the same event; the weaker one is dropped. |

Items at a contig end carry `ignore=CONTIG_END` and are dropped from the report entirely.

## GenomeDiff fields

The positional fields are specified in
[GenomeDiff File Format](genomediff-file-format.md#pd-pair-distance-evidence). `PD` shares `DP`'s
two-sided junction shape: `side_1_position` is the last retained base of the left flank (strand
always -1) and `side_2_position` the first retained base of the right flank (strand always +1), so
the two sides bracket exactly the reference bases the event removed, and are adjacent when it removed
none.

Notable `name=value` pairs:

*   **size_shift** *\<int32>* — the estimated shift in mapping distance, in bases. Positive means
    reference sequence is missing from the sample (a deletion of that many bases, which is also the
    separation of the two sides). Negative means sequence was added, and the sides are adjacent: `PD`
    cannot place inserted sequence in the reference, and it cannot tell a pure insertion of *I* bases
    from a replacement of *k* reference bases by *k + I* inserted ones, so it reports the minimal
    reading and lets `size_shift` carry the length.
*   **size_shift_lower**, **size_shift_upper** — profile-likelihood interval for `size_shift`. A call
    whose interval includes zero is rejected with `PAIR_DISTANCE_SIZE`.
*   **score** *\<float>* — the genome-wide E-value score. **This is the test that decides a `PD`
    call.**
*   **seed_z_score** — the rank-sum statistic that opened the candidate region. Sign follows the
    direction of the shift.
*   **position_range** *\<uint32>* — width of the interval of positions the supporting pairs agree
    on. Zero means they pin a single base, which normally happens only when the call snapped onto a
    junction.
*   **shifted_pair_count**, **normal_pair_count**, **ambiguous_pair_count** — of the pairs sampled
    across this point, how many had a distance that decisively supports the call, decisively argues
    against it, and neither. The frequency is computed from the first two.
*   **total_pair_count**, **candidate_covering_count** — the pairs considered, and the peak count
    while the candidate region was open.
*   **distinct_pair_count** — distinct fragment ends among the supporting pairs, so that PCR
    duplicates of one molecule cannot carry a prediction.
*   **repeat_name**, **repeat_name_evidence** — the repeat family the shift's length is consistent
    with, and how that was established. `repeat_name_evidence=junction` means it was read off a
    validated split-read junction whose other side lands in the element; otherwise it was inferred
    from size alone.
*   **repeat_size_candidates**, **repeat_size_difference** — when more than one family fits within
    the estimator's scatter, the honest list rather than a guess.
*   **snapped_to_junction** — present when the coordinates were taken from a validated split-read
    junction inside the interval, and are therefore exact to the base.

## In the HTML report

Accepted items appear on `index.html` under `Unassigned pair distance evidence`; rejected items on
`marginal.html` under `Marginal pair distance evidence`.

`seq id`\
Identifier for the reference sequence. Both sides are always the same sequence.

`position`\
The two bracketing coordinates.

`size`\
`size_shift` — positive for a deletion, negative for an insertion.

`size range`\
The profile-likelihood interval on `size_shift`. An interval spanning zero is a rejection.

`shifted`, `normal`\
Pairs decisively supporting and decisively opposing the call. These are the two terms of `freq`.

`distinct`\
`distinct_pair_count` — how many distinct molecules the support represents.

`freq`, `range`\
Local variant frequency and its 95% confidence interval; the cutoff is applied to the interval. See
[`freq` and `range`](evidence-overview.md#freq-and-range).

`score`\
The genome-wide score. See [`score`](evidence-overview.md#score).

`gene`, `product`\
Annotation at each side.

## Options

`PD` is predicted by default on paired-end data.

`--no-pair-distance-prediction`\
Do not predict pair distance (PD) evidence. `PD` is also skipped automatically on a single-end run or
under `--no-paired-mapping`, and it only considers pairs whose two mates map to the same reference
sequence.

!!! note "`--predict-pair-distance` is deprecated"
    Still accepted so existing command lines keep working, and hidden from `breseq --help`, but it no
    longer switches anything on.

`--pair-distance-seed <int>` (default 3)\
Minimum number of read pairs in the matching distribution tail whose gap covers a position, required
to seed a candidate region.

`--pair-distance-seed-z <float>` (default 0, derived)\
Minimum |z| of the rank-sum statistic required to seed a candidate region. `0` derives it from the
reference length so the genome-wide false-seed rate stays near 1%. Sets the
[candidate seed](#candidate-seed) gate.

`--pair-distance-maximum-span <int>` (default 0, derived)\
Ignore read pairs mapping farther apart than this when seeding. Pairs this long are `DP`'s business,
and an unbounded distance would make the seeding window unbounded. `0` derives it as twice the
paired-mapping distance cutoff.

`--pair-distance-minimum-pairs <int>` (default 3)\
Only accept PD evidence supported by at least this many shifted pairs. Produces
`PAIR_DISTANCE_COUNT`.

`--pair-distance-minimum-distinct <int>` (default 2)\
Only accept PD evidence whose supporting pairs start at at least this many distinct positions.
Produces `PAIR_DISTANCE_DUPLICATES`.

`--pair-distance-minimum-shift <int>` (default 0, derived)\
Only accept PD evidence whose estimated size shift is at least this many bases. `0` derives it from
the width of the paired-mapping distance distribution.

`--pair-distance-tail-quantile <float>` (default 0.05)\
Quantile of the paired-mapping distance distribution defining its tails, for the seed's counting
test.

`--pair-distance-score-cutoff <float>` (default 3)\
Log10 E-value cutoff. `0` = OFF. Sets the [score cutoff](#score-cutoff) gate and produces
`PAIR_DISTANCE_SCORE`.

`--pair-distance-frequency-cutoff <float>` (default: follows `--polymorphism-frequency-cutoff`)\
Only accept PD evidence when the lower 95% confidence bound on its local variant frequency —
shifted pairs divided by shifted plus normal — is at or above this value. `0` = OFF. Produces
`PAIR_DISTANCE_FREQUENCY`.

## Worked examples

Both examples come from `long_ltee_ara_m3_32k_mp2800`, a mate-pair library with a ~2.8 kb insert —
the only test in the tree where the pair machinery sees a kilobase-scale span.

### An insertion, sized but not placed

```text title="tests/long_ltee_ara_m3_32k_mp2800/expected.gd"
PD	1910	.	REL606	546190	-1	REL606	546191	1	ambiguous_pair_count=634
	distinct_pair_count=253	frequency=0.7792	normal_pair_count=217	position_range=21
	repeat_name=IS1	repeat_name_evidence=junction	score=30.8	seed_z_score=-35.15
	shifted_pair_count=766	size_shift=-760	size_shift_lower=-762	size_shift_upper=-760
	gene_name=ybcQ	gene_position=coding (345/384 nt)
```

*(Fields trimmed and wrapped for display.)*

Reading it: `size_shift` is **negative**, so sequence was *added* — the spanning pairs map closer
together than the library says they should. The two sides are adjacent (546190 and 546191), which is
what `PD` always reports for an insertion, because it cannot place the new sequence in reference
coordinates. `size_shift_lower` and `size_shift_upper` are tight and well clear of zero, so the
direction is not in doubt.

`repeat_name_evidence=junction` is the important qualifier: the IS1 identification was **read off a
validated split-read junction** whose other side lands in the element, not inferred from the length
alone. Where `PD` has to infer from size, items instead carry a `repeat_size_candidates` list —
IS150 and IS186 differ by only a hundred bases, which is inside this estimator's scatter, so naming
one of them from size would be a guess.

`position_range` shows the supporting pairs agree on the breakpoint only to within a few tens of
bases. That is why this insertion stays `PD` evidence and is not promoted to a `MOB` mutation, which
would need the base.

### A deletion, corroborated by missing coverage

```text title="tests/long_ltee_ara_m3_32k_mp2800/expected.gd"
PD	1911	.	REL606	547517	-1	REL606	549789	1	ambiguous_pair_count=0
	candidate_covering_count=1096	distinct_pair_count=393	frequency=1.0000
	frequency_lower=0.9973	normal_pair_count=0	position_range=1127	score=87.3
	seed_z_score=56.91	shifted_pair_count=1096	size_shift=2271
	size_shift_lower=2267	size_shift_upper=2276
```

Here `size_shift` is **positive**: the sample is missing reference sequence, and the two sides
bracket exactly the bases removed. `normal_pair_count` is zero — not one pair spanning this point
maps at a normal distance — so `frequency` pins at 1.0 with a tight lower bound, consistent with a
clonal deletion. The `score` is an order of magnitude above the cutoff.

Note `position_range` is over a kilobase here. `PD` is confident *that* sequence is missing and
confident *how much*, while being quite vague about exactly where the boundary falls — a good
illustration of what this statistic does and does not deliver.

!!! tip "PD counts that look too high"
    A run predicting many more `PD` items than its neighbours is not automatically over-predicting.
    Check the `multi-mapped pairs excluded` row of the metrics table: a repeat-rich reference with a long
    insert legitimately produces more real events *and* refuses more ambiguous pairs. On this test,
    most accepted items have a negative `size_shift` matching an annotated IS element's length, and
    the positive-shift ones sit on missing-coverage regions — independent corroboration in both
    directions.

## See also

- [Evidence overview](evidence-overview.md) — how `PD` compares with `JC`, `SC`, `DP` and `MP`
- [DP: Discordant pair evidence](evidence-dp.md) — the per-pair test `PD` complements
- [MC: Missing coverage evidence](evidence-mc.md) — what corroborates a positive `size_shift`
