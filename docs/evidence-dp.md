# DP: Discordant pair evidence

A pair of breakpoints joined by read pairs that are, **individually**, mapped wrongly: in the wrong
orientation, farther apart than the library allows, or with their two mates on two different
reference sequences. It is predicted by default on paired-end data; turn it off with
`--no-discordant-pair-prediction`. It is not currently promoted to any mutation type, so accepted
`DP` items appear as unassigned evidence.

## At a glance

| | |
|---|---|
| Section in `summary.html` | `Paired-End Mapping Distance Information` (the distance fit it depends on) |
| Tables in `summary.html` | `Discordant pair (DP) evidence metrics` + `... gates` |
| Banner in `index.html` | `Unassigned discordant pair evidence` |
| Sort order in the GenomeDiff file | 16 |
| Enabled by | on by default; turn off with `--no-discordant-pair-prediction` |
| Requires | paired reads |
| Promotes to | *(nothing — reported as evidence only)* |
| Rejected items visible in HTML | yes, on `marginal.html` (top 20) |

## The signal

Every read pair carries three facts beyond where its two mates landed: which sequence each landed on,
how far apart they are, and which way each points. A concordant pair agrees with the library — same
sequence, facing inwards, separated by roughly the fragment size. `DP` collects the pairs that do
not, and looks for **breakpoints that many such pairs agree on**.

A pair is discordant when any of three things is true:

- Its mates point the wrong way relative to each other (an inversion signature).
- Its mates map farther apart than the library's distance cutoff (a deletion signature).
- **Its mates map to two different reference sequences** — a plasmid integration, or a translocation
  between two contigs.

That third case is what makes `DP` unique among the evidence types. It is the only one that is
genuinely **two-sided across sequences**: `side_1_seq_id` and `side_2_seq_id` are independent, so a
translocation between two contigs is a single `DP` item rather than two unrelated observations.

!!! note "Cross-sequence pairs have no orientation"
    A pair whose mates land on two different sequences has no meaningful within-sequence orientation
    and no meaningful distance, so it is binned into an orientation slot of its own rather than being
    forced into one of the same-sequence categories. Do not expect a letter-pair orientation on such
    an item.

## What this evidence cannot see

- **Small events, in either direction.** `DP` tests each pair *on its own*, so a pair must clear the
  library's discordance cutoff to count at all — and there is no cutoff on the short side whatsoever.
  An event of a few hundred bases shifts every spanning pair a little, and no individual pair enough.
  That entire band belongs to [`PD`](evidence-pd.md).
- **Base-pair resolution.** A `DP` breakpoint is bracketed by where the discordant pairs' mates fall,
  which is a window the width of the fragment size, not a base.
- **Anything on single-end data.**

Where a `DP` and a [`PD`](evidence-pd.md) item describe the same breakpoint, the `DP` is removed:
`PD` uses the whole pair population where `DP` uses only its tail.

!!! warning "DP is deliberately exempt from the contig-end rule"
    Most evidence types are dropped with `ignore=CONTIG_END` when they sit near the end of a
    reference sequence, because that is where their statistic stops being measurable. `DP` is
    excluded from that rule on purpose. A cross-sequence `DP` sits at two contig ends **by
    construction** — that is what a translocation between two contigs looks like — so applying the
    rule would silently discard every translocation call.

## From reads to evidence

1. **Fit the library.** Measure the paired-mapping distance distribution and the dominant
   orientation, and derive a discordance cutoff from them. All of this is reported in the
   `Paired-End Mapping Distance Information` section of `summary.html`.
2. **Mark discordant pairs.** Each pair is classified against that fit. Cross-sequence pairs are
   marked as such.
3. **Seed candidate regions.** Where at least `--discordant-pair-seed` discordant pairs fall within
   one paired-mapping-distance window, a candidate region is opened.
4. **Pair up the regions.** Two candidate regions sharing at least
   `--discordant-pair-minimum-pairs` read pairs describe a putative junction between them. The
   effective floor here is raised from a background fit, described below.
5. **Count the opposition.** At each side, count the concordant pairs that span the breakpoint. Those
   molecules carry unbroken reference sequence there, so they argue against the call.
6. **Pool siblings.** An insertion creates *two* junctions at one point, and the pairs that would
   have spanned it are divided between them. `--discordant-pair-sibling-window` lets the two be
   judged together, so the weaker one is not rejected for carrying only its share.
7. **Gate.** Apply the frequency and skew tests below.

## Metrics

The `Discordant pair (DP) evidence metrics` table in `summary.html` reports these row for row. Neither
accepts nor rejects anything — they are what the [gates](#gates) below are measured against.

### paired-mapping distance

The library's fitted median, MAD, discordance cutoff and orientation. If the `basis` column notes
that *fragments are shorter than two reads*, most pairs overlap and have no gap for a breakpoint to
fall in, which limits what `DP` can see at all.

### minimum shared read pairs

How many read pairs two candidate regions must share for the junction between them to be examined.
`--discordant-pair-minimum-pairs` sets the floor, but the value in force is often **raised from it**,
and the `basis` column says why: a background of spurious discordant pairs is fitted across the run's
own candidate junctions, and the floor is lifted until the expected number of chance junctions falls
below `--discordant-pair-background-e-value-cutoff`. The fitted null — a negative binomial mean and
size, a Poisson mean, or `not fit` — is reported here.

## Gates

The `Discordant pair (DP) evidence gates` table lists each decision, the rule in force, and the
`reject=` value an item picks up when it fails.

### local frequency

`--discordant-pair-frequency-cutoff`. Rejects as `DISCORDANT_PAIR_FREQUENCY`. The **exact lower
confidence bound** on
`discordant / (discordant + concordant)` pairs spanning the breakpoint, not the point estimate. It
tracks the prediction mode, defaulting to `--polymorphism-frequency-cutoff`. See
[`freq` and `range`](evidence-overview.md#freq-and-range).

### concordant pair skew

`--discordant-pair-skew-cutoff`, default 3.0. Rejects as `CONCORDANT_PAIR_SKEW`. Asks whether the
*shortage* of concordant pairs spanning the breakpoint is itself surprising — if the junction is
real, molecules crossing it should be depleted.

Note the rule reads `> 3.0`, not `≥ 3.0` like every other gate. This is the one test that rejects an
item whose score is **above** its cutoff: a high skew means too many concordant pairs still span the
breakpoint, which argues *against* a junction being there.

This gate is **conditional on having power**. The `basis` column reports the expected number of
concordant pairs spanning a normal position at this coverage, and then either *"so the test can
discriminate"* or *"too few for the test to reject on merit — reported only"*. Below
`--discordant-pair-minimum-crossing` (default 10) the skew is computed and displayed but never
rejects, because a short fragment distribution leaves almost no concordant pair spanning any position
and the test would fire everywhere. In that regime `DP` is accepted or rejected on local frequency
alone.

The null is the empirical crossing distribution where the reference is large enough to support it,
and a negative binomial fit otherwise; the `basis` column says which.

### outcome

Items examined, accepted, rejected by frequency, rejected by skew, and three dispositions unique to
`DP`: `circular-origin artifact`, `dropped with no read pair bridging the placed breakpoints`, and
`folded into another candidate placed at the same breakpoint`.

### A real pair of tables

From `long_ltee_ara_p1_50k_pe101`, a 2x101 paired library of an *E. coli* clone:

**Discordant pair (DP) evidence metrics**

| metric | value | basis |
|---|---|---|
| paired-mapping distance | 232 median, 86 MAD, 956 cutoff, FR | read length 101 |
| minimum shared read pairs | 5 (raised from 2) | spurious-pair background fit to 105 candidate junctions: Poisson mean 4.0000e-01; E-value cutoff 0.050 |

**Discordant pair (DP) evidence gates**

| gate | rule | rejects as | basis |
|---|---|---|---|
| local frequency | ≥ 10.0% | `DISCORDANT_PAIR_FREQUENCY` | exact lower confidence bound on discordant / (discordant + concordant) pairs spanning the breakpoint |
| concordant pair skew | > 3.0 | `CONCORDANT_PAIR_SKEW` | expected concordant pairs spanning a normal position: 23.2 in REL606 ≥ 10.0, so the test can discriminate; null = empirical over 4393233 positions |
| outcome | 98 examined | | 95 accepted, 2 rejected (frequency), 2 rejected (skew), 1 circular-origin artifact |

Two rows repay attention. `minimum shared read pairs` reads **5 (raised from 2)**: the option's
default floor was lifted because a background fitted to this run's own 105 candidate junctions said
that two shared pairs would arise by chance too often. That is the adaptive behaviour described
above, visible in the report — and it belongs in the metrics table rather than the gates table
because it shapes which junctions are examined at all, not which are rejected.

`concordant pair skew` says the expected number of concordant pairs spanning a normal position is
comfortably over the `--discordant-pair-minimum-crossing` floor, **so the test can discriminate** —
this run has the power to use the skew test, and duly rejects two items with it. On a library with a
shorter fragment distribution this basis would instead read *too few for the test to reject on
merit*, and only the frequency gate would be operating.

The single `circular-origin artifact` is the item shown in the first worked example below.

## Why items are rejected

| `reject` value | Shown in the report as | What it means here |
|---|---|---|
| `DISCORDANT_PAIR_FREQUENCY` | Lower confidence bound on discordant pair frequency below cutoff. | Too few of the pairs spanning this breakpoint are discordant, relative to those that cross it normally. On the per-item detail page this generic sentence is replaced with the actual numbers. |
| `CONCORDANT_PAIR_SKEW` | Concordant pair skew score above cutoff. | There are *more* concordant pairs spanning the breakpoint than a real junction there would allow. Only applied when the run has the power for it — see the gate above. |

Items spanning the origin of a circular sequence carry `ignore=CIRCULAR_CHROMOSOME` and are dropped
from the report rather than rejected: position 1 and position *L* are physically adjacent, so pairs
straddling them are not discordant at all. See
[accepted, rejected, ignored](evidence-overview.md#accepted-rejected-ignored).

## GenomeDiff fields

`DP` shares `JC`'s two-sided shape — `side_1_seq_id`, `side_1_position`, `side_1_strand`,
`side_2_seq_id`, `side_2_position`, `side_2_strand` — but with **independent sequence identifiers on
the two sides**. See [GenomeDiff File Format](genomediff-file-format.md#dp-discordant-pair-evidence).

Notable `name=value` pairs:

*   **discordant_count** — read pairs supporting this junction.
*   **distinct_discordant_count** — distinct fragment ends among them, so PCR duplicates of one
    molecule cannot carry a prediction.
*   **candidate_discordant_count** — the count while the candidate region was open, before the
    breakpoints were finally placed.
*   **concordant_count** — pairs spanning the breakpoint normally. The opposition, and the second
    term of the frequency.
*   **expected_concordant_count** — how many concordant pairs *would* be expected to span a normal
    position at this coverage. Compare with `concordant_count`: this is what the skew test tests.
*   **side_1_concordant_count**, **side_2_concordant_count**, **side_1_discordant_count**,
    **side_2_discordant_count** — the same counts measured at each side separately.
*   **side_1_coverage**, **side_2_coverage** — local read-depth coverage at each side. `NA` where the
    side falls in repetitive sequence.
*   **side_1_redundant**, **side_2_redundant** — set when that side maps to more than one place, so
    the coordinate shown is one example.
*   **background_e_value** — the expected number of junctions this well supported arising by chance
    across the run's candidate junctions, from the fitted spurious-pair background.
*   **neg_log10_discordance_p_value** — the concordant pair skew score, compared against
    `--discordant-pair-skew-cutoff`.
*   **pooled_discordant_count** — present when sibling junctions were judged together.
*   **new_junction_coverage** — coverage attributable to the new junction.
*   **frequency**, **frequency_lower**, **frequency_upper** — `discordant / (discordant +
    concordant)` with exact 95% bounds. The cutoff is applied to `frequency_lower`.
*   **ignore** — `CIRCULAR_CHROMOSOME` for an origin-spanning artifact.

## In the HTML report

Accepted items appear on `index.html` under `Unassigned discordant pair evidence`; rejected items on
`marginal.html` under `Marginal discordant pair evidence`, sorted by supporting pair count.

Like `JC`, each `DP` row is **two sub-rows**, one per side of the junction.

`seq id`\
Identifier for the reference sequence at each side. These can differ — that is a translocation.

`position`\
The bracketing coordinate at each side.

`pairs (cov)`\
Concordant pairs spanning the breakpoint, with local coverage.

`disc pairs (cov)`\
Discordant pairs supporting it.

`freq`, `range`\
`discordant / (discordant + concordant)` and its confidence interval; the cutoff is applied to the
interval's lower bound.

`skew`\
`neg_log10_discordance_p_value`. May be reported without being applied — see the
[concordant pair skew](#concordant-pair-skew) gate.

`annotation`, `gene`, `product`\
Annotation at each side.

## Options

`DP` is predicted by default on paired-end data.

`--no-discordant-pair-prediction`\
Do not predict discordant read-pair (DP) evidence. `DP` is also skipped automatically on a single-end
run or under `--no-paired-mapping`, since it needs both mates of a pair.

!!! note "`--predict-discordant-pairs` is deprecated"
    Still accepted so existing command lines keep working, and hidden from `breseq --help`, but it no
    longer switches anything on.

`--discordant-pair-seed <int>` (default 3)\
Minimum discordant read pairs within a paired-mapping-distance window required to seed a candidate
region.

`--discordant-pair-minimum-pairs <int>` (default 2)\
Minimum read pairs shared by two candidate regions for the junction between them to be examined at
all. Sets the floor of the [minimum shared read pairs](#minimum-shared-read-pairs) gate, which the
background fit may raise.

`--discordant-pair-background-e-value-cutoff <float>` (default 0.05)\
Reject evidence whose supporting pair count is expected to arise this many times or more across all
candidate junctions by chance, given the genome-wide background of spurious discordant pairs.
`0` = OFF.

`--discordant-pair-frequency-cutoff <float>` (default: follows `--polymorphism-frequency-cutoff`)\
Only accept when the lower confidence bound on the local variant frequency is above this value.
`0` = OFF. Produces `DISCORDANT_PAIR_FREQUENCY`.

`--discordant-pair-skew-cutoff <float>` (default 3.0)\
Cutoff for the concordant pair skew score. Produces `CONCORDANT_PAIR_SKEW`.

`--discordant-pair-minimum-crossing <float>` (default 10.0)\
Only apply the skew test when at least this many concordant pairs are expected to span a normal
position at the breakpoint's local coverage. Below this the test has no power and does not reject.
`0` = always apply.

`--discordant-pair-sibling-window <int>` (default 20)\
Maximum distance, in bases, between the two breakpoints of one insertion for their `DP` items to be
judged together by the skew test. An insertion creates two junctions at one point and divides the
spanning pairs between them; without this each is tested against the full crossing expectation and
the weaker one is rejected. `0` disables the pooling.

## Worked examples

### An ignored item: the circular chromosome

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
DP	1003	.	REL606	1	1	REL606	4629812	-1	candidate_discordant_count=29
	concordant_count=0.0	discordant_count=29	distinct_discordant_count=29
	expected_concordant_count=23.2	frequency=1.0000	ignore=CIRCULAR_CHROMOSOME
	side_1_gene_name=–/thrL	side_2_gene_name=lasT/–
```

*(Fields trimmed and wrapped for display.)*

This is not a variant, and _breseq_ knows it. The two sides are position 1 and the last base of
REL606 — the two ends of a **circular** chromosome, which are physically adjacent in the cell. Every
fragment that happens to straddle the origin looks, to a linear coordinate system, like a pair mapped
4.6 Mb apart.

The item carries `ignore=CIRCULAR_CHROMOSOME` rather than a `reject=` field, and it is dropped from
the report entirely. That is the right disposition: this is an artifact of how the reference is
*written down*, not a judgement about the data. Note that on the numbers alone it would have been
accepted — the frequency is 1.0 and there is not one concordant pair opposing it.

### A rejected item: too well opposed

```text title="tests/long_ltee_ara_p3_30k_pe150/expected.gd"
DP	2794	.	REL606	2962	1	REL606	7396	1	background_e_value=4.424e-05
	candidate_discordant_count=7	concordant_count=46.0	discordant_count=2
	distinct_discordant_count=2	expected_concordant_count=62.8	frequency=0.0417
	frequency_lower=0.0075	frequency_upper=0.1254	neg_log10_discordance_p_value=7.0
	reject=DISCORDANT_PAIR_FREQUENCY,CONCORDANT_PAIR_SKEW
	side_1_concordant_count=43	side_2_concordant_count=49
```

Two gates fail here, and they are two views of the same problem: **the reference sequence between
these points is demonstrably still intact.**

`concordant_count` counts pairs that span the putative breakpoint normally, and they outnumber the
supporting `discordant_count` by more than twenty to one. So `frequency` is tiny, and
`frequency_lower` falls far below the cutoff — hence `DISCORDANT_PAIR_FREQUENCY`. This is the
frequency test doing exactly what it should: whatever these two pairs represent, it is not present in
any meaningful fraction of the sample.

`CONCORDANT_PAIR_SKEW` says the same thing from the other direction. Compare `concordant_count` with
`expected_concordant_count`: the observed number is close to what an *undisturbed* position at this
coverage would show. A real junction depletes the molecules crossing it; nothing here is depleted.
`neg_log10_discordance_p_value` is well above the 3.0 cutoff.

Note also `distinct_discordant_count` — the two supporting pairs are two distinct molecules, so this
is not a PCR-duplicate artifact. It is simply two pairs that mapped oddly, which in a run with
millions of pairs is unremarkable, and `background_e_value` quantifies exactly that.

This run rejects the large majority of its `DP` items, which is normal for a deeply-sequenced sample:
seeding is deliberately permissive, and the gates are where the work is done.

### When DP correctly predicts nothing

On the `long_ltee_ara_m3_32k_mp2800` mate-pair library, `DP` fires **zero** times, and that is
expected rather than a failure. The distance distribution there is tightly unimodal about a ~2.8 kb
insert with a discordance cutoff well above it, so essentially no pair is an outlier *individually* —
which is precisely the regime `DP` cannot address and [`PD`](evidence-pd.md) exists for. On that same
run `PD` finds a couple of dozen well-corroborated events.

If a change ever made `DP` fire in bulk on such a library, it would be calling library structure, not
variants.

## See also

- [Evidence overview](evidence-overview.md) — how `DP` compares with `JC`, `SC`, `MP` and `PD`
- [PD: Pair distance evidence](evidence-pd.md) — the collective test that covers `DP`'s blind band
- [JC: New junction evidence](evidence-jc.md) — base-pair resolution on the same events
