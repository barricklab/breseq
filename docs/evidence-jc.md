# JC: New junction evidence

Two disjoint pieces of the reference sequence that are **adjacent in the sample**, demonstrated by
reads that align across the join. `JC` is the most demanding evidence type and the most informative:
it is the only one that locates a structural breakpoint to the base *and* identifies both sides. It
is on by default.

## At a glance

| | |
|---|---|
| Section in `summary.html` | `New Junction Evidence` |
| Metrics / gates tables | *(none — see `Junction Candidates Tested`, `Junction Skew Score Calculation` and `Final Junction Predictions`)* |
| Banner in `index.html` | `Unassigned new junction evidence` |
| Sort order in the GenomeDiff file | 12 |
| Enabled by | on by default (disable with `--no-junction-prediction`) |
| Requires | nothing (works on single-end data) |
| Promotes to | `MOB`, `DEL`, `INS`, `AMP`, `CON`, `INT`, `INV` |
| Rejected items visible in HTML | yes, on `marginal.html` (top 10) |

## The signal

If a mobile element inserts itself into a gene, the sample now contains a sequence that reads
"…gene…" up to a point and then "…IS element…" from there on. A read spanning that point cannot align
end-to-end anywhere in the reference — but each of its two halves aligns perfectly, to two places
that are nowhere near each other.

`JC` looks for exactly these **mosaic** alignments, and for junctions supported by enough of them,
distributed well enough, to be believable.

Because both halves align, `JC` recovers everything: the exact base at which the join occurs on each
side, which strand each side continues on, and the identity of both sequences. That is why `JC` is
the only evidence type that can promote to a `MOB` mutation — a mobile element insertion has to be
placed to base-pair resolution and named, and nothing else can do both.

## What this evidence cannot see

- **Sequence absent from the reference.** Both halves must align *somewhere*. An insertion of novel
  sequence has one half with nowhere to go, which is [`MP`](evidence-mp.md)'s signal, or
  [`SC`](evidence-sc.md)'s if the tail is readable but unplaceable.
- **Events whose breakpoints no read happens to span.** A junction needs a read that crosses it with
  enough sequence on both sides to seed two alignments. At low coverage, or with short reads, many
  real junctions are simply never bridged. This is where the pair-based types
  ([`DP`](evidence-dp.md), [`PD`](evidence-pd.md)) earn their keep — they need only that a *fragment*
  spans the point, not a read.
- **Which copy of a repeat was involved.** When a junction side falls in repetitive sequence, the
  coordinate shown is one example among several equally good ones, flagged with `side_n_redundant`.

## From reads to evidence

### Identifying candidate junctions

In a pre-processing step, all read alignments with insertions or deletions of more than 2 bases are
split into their constituent sub-alignments. This tends to be more accurate than looking for such
mutations as [read alignment evidence](evidence-ra.md), because gaps larger than a couple of bases
are hard to align consistently, especially in simple sequence repeats.

Next, for each read with multiple alignments to the reference, all pairs of alignments are tested to
find cases where:

1.  One alignment begins with the first base of the read.
2.  Both alignments together cover a number of bases in the read that is more than 2 bases longer
    than the length covered by any other single alignment.
3.  Both alignments contain at least 5 read bases that do not overlap the other.
4.  One alignment contains at least 10 read bases that do not overlap the other.
5.  There are at most 20 bp unique to the read between matches to the reference.

If a pair of alignments passes, _breseq_ generates the putative sequence of the new junction from the
reference sequence plus any intervening bases unique to the read.

Where the two alignments **overlap**, because some sequence could be assigned to either location,
_breseq_ first trims each alignment to remove portions with mismatched bases or indels, then assigns
as much of the overlap as possible to the side that maps uniquely, or failing that to the side with
the lower reference coordinate. A location is considered non-unique if it has repeat matches. If
`repeat_region` annotation exists in the reference, _breseq_ prefers junctions that exactly overlap
those boundaries.

The candidate junction sequence includes as many flanking reference bases on each end as the longest
read in the data set. (For a 36 bp data set with a five-base overlap, the candidate is 36 × 2 − 5 =
67 bases.)

After processing every read, candidate junctions matching the same reference sequences are combined
and given a **position-hash score**: a count of the number of *different* start position–strand
combinations observed among the supporting reads. This favours junctions supported by reads that are
evenly distributed on each strand and at different offsets from the junction point. Pathological
candidates tend to be supported only by reads that barely overlap the junction and all lie on one
strand.

Candidates are sorted by position-hash score, ties broken by a **minimum-overlap score** that sums,
over all reads, the smaller of each read's two overlaps. Top-scoring candidates are retained until
their cumulative length would exceed 0.1× the reference length or their number would exceed 5000.

### Scoring and accepting junctions

New junctions may also be supported by reads that do not overlap both sides enough to seed alignments
during mapping. To include these, _breseq_ performs a **second alignment step**, mapping all reads to
the candidate junction sequences. For each read it then decides whether its best alignment is to a
junction candidate or to the reference, scoring alignments as matched reference bases minus indel
positions. Alignments covering fewer than 28 bases of the read are discarded.

A position-hash score is recomputed for each candidate from the reads that map best to it. Junctions
are tested in order from most best-alignments to fewest. Reads mapping equally well to the reference
and to one or more junctions are included in these scores.

If a junction is accepted, reads that map equally well elsewhere are assigned to it and removed from
further consideration. Any read still unused after all candidates are tested is assigned to the
reference genome.

For accepted junctions, the ends of the aligning reads are re-added to the alignment database as
split sub-alignments, resolving ambiguous bases so each read base aligns to only one reference base.
These split reads are recognisable in the output by the suffixes **-M1** and **-M2**.

## Gates

`JC` prints no gates table. Its settings appear in three tables in the `New Junction Evidence`
section of `summary.html`.

### Coverage evenness (the "skew" score)

The decisive test. A candidate is accepted if its position-hash score is good enough, measured as
`neg_log10_pos_hash_p_value` — the **skew** score.

In consensus mode the skew is computed against the fitted
[read coverage distribution](methods.md#read-coverage-distribution). _breseq_ tracks what fraction of
position–strand combinations genome-wide have a read starting there, and uses that baseline to
compute the chance of at least one read starting at any given position and strand at a given depth.
The chance of observing the junction's actual position-hash score then follows a binomial with twice
the read length as the number of trials. The negative log10 of that probability is the skew.

A **higher** skew means the observation was less likely — either too few reads, or reads too biased
towards particular start positions. With default settings a junction fails at a skew above 3.0, a
probability below 0.001.

Since version 0.34.0 a saturation correction applies: at very high coverage almost every
position–strand combination has a read starting there, but local shearing and amplification biases
mean some regions never reach that level. Without correction those regions would produce high skew
values and lose real junctions. `--junction-minimum-pr-no-read-start-per-position` (default 0.10)
sets the floor on the fraction of unoccupied start combinations; set it to 0 for pre-0.34.0
behaviour.

In **polymorphism mode** the skew is not useful, because a variant's true coverage is an unknown
fraction of the average. Junctions are assigned a skew of `NT` (Not Tested) and other criteria decide.

### Additional acceptance criteria

Beyond the skew, a junction must:

1.  Be supported by reads mapping to **both strands** of the predicted junction.
2.  Have reads that extend at least **14 bp** into each side of the reference.
3.  Have reads on **each strand** extending at least **9 bp** into each side.
4.  Have reads where the side with the smallest reference overlap extends at least **3 bp** into the
    reference on each side.

### `Junction Candidates Tested`

An `option` / `limit` / `actual` table reporting how much of the available junction evidence was
examined: alignment pairs examined, the coverage-evenness threshold for candidates, how many
candidates were tested, and their total length as a factor of the genome length. Values may read
`NO NO-LIMIT`.

!!! tip "When this table shows the limit was hit"
    If `actual` equals `limit` for *alignment pairs examined*, the run discarded most of its junction
    evidence and real junctions may have been crowded out. This happens on libraries that generate
    enormous numbers of one-off chimeric molecules — mate-pair circularisation junctions are the
    classic case. Each yields a candidate with position-hash 1 that the threshold discards, but they
    flood the budget first. Raising `--junction-alignment-pair-limit` recovers the real junctions.

### `Junction Skew Score Calculation`

Per reference sequence: `pr(no read start)` and the `coverage model` used
(`negative binomial` or `empirical (N positions)`).

### `Final Junction Predictions`

An `option` / `value` table listing the accepted thresholds, including the required coverage evenness
score, the minimum probability assigned that no mapped read will start at a given position and
strand, whether suboptimal matches are allowed, the skew cutoff, and the required overlap into each
uniquely aligned side.

## Why items are rejected

| `reject` value | Shown in the report as | What it means here |
|---|---|---|
| `COVERAGE_EVENNESS_SKEW` | Coverage evenness skew score above cutoff. | The decisive gate. The supporting reads are too few, or too clustered at particular start positions and strands, to be a fair sample of a real junction. Usually the correct verdict on a chimeric molecule or a mapping artifact. |
| `FREQUENCY_CUTOFF` | Frequency below/above cutoff threshold. | The junction is present in too small a fraction of the sample, judged on the confidence bound rather than the point estimate. Also applied during mutation prediction when a `MOB` candidate's junctions disagree about frequency. |
| `BETWEEN_TWO_JUNCTION_ONLY_SEQUENCES` | Between two junction-only reference sequences. | Both sides fall in sequences supplied only as junction targets, so the junction says nothing about the sample's genome. |
| `INDEL_HOMOPOLYMER` | Polymorphic indel expands or contracts a homopolymer stretch. | Applied to `JC` evidence for small indels that could be homopolymer slippage. |
| `POLYMORPHIC_INDEL` | Indel polymorphism suppressed by --polymorphism-no-indels. | Set when `--polymorphism-no-indels` is in force. |

Junctions spanning the origin of a circular sequence carry `ignore=CIRCULAR_CHROMOSOME`, and those
near a contig end `ignore=CONTIG_END`; both are dropped from the report rather than rejected. See
[accepted, rejected, ignored](evidence-overview.md#accepted-rejected-ignored).

## GenomeDiff fields

Positional fields: `side_1_seq_id`, `side_1_position`, `side_1_strand`, `side_2_seq_id`,
`side_2_position`, `side_2_strand`, `overlap` — see
[GenomeDiff File Format](genomediff-file-format.md#jc-new-junction-evidence).

The strands are given as −1 or +1 to indicate how the read leads up to the junction on the first side
and continues after it on the second. The most common type of junction has side 1 strand −1 and side
2 strand +1, which can indicate a deletion.

Notable `name=value` pairs:

*   **pos_hash_score**, **max_pos_hash_score** — the position-hash score and the maximum it could
    have reached. The ratio is the more meaningful quantity.
*   **neg_log10_pos_hash_p_value** — the **skew** score, compared against the cutoff. `NT` in
    polymorphism mode.
*   **min_overlap_score** — the tie-breaking minimum-overlap score.
*   **coverage_plus**, **coverage_minus** — supporting reads by strand. Both must be non-zero.
*   **total_reads**, **total_non_overlap_reads** — reads mapping to the junction, and those excluding
    the overlap region.
*   **max_left**, **max_right**, and their `_plus`/`_minus` variants — the furthest any read extends
    into each side, overall and per strand.
*   **max_min_left**, **max_min_right**, and variants — the acceptance criteria above are expressed in
    these. A zero in `max_min_right` means no read reaches far enough into the right side.
*   **side_1_redundant**, **side_2_redundant** — set when that side maps to more than one place, so
    the coordinate shown is one example. Such sides are highlighted orange in the HTML.
*   **new_junction_frequency**, and its bounds — how much of the sample carries the junction.
*   **alignment_overlap** — bases the two sides have in common.

## In the HTML report

Each `JC` row is **two sub-rows**, one per side. A sub-row highlighted **orange** means that side
maps ambiguously to more than one place, and the coordinate shown is an example.

`* link`\
Links to a results page showing the sequence of the new junction as the reference, and all reads
aligned to the junction.

`? links`\
Links to results pages for each side of the junction, showing the reference sequence at that site and
any reads that aligned better to this original sequence than to the new junction. In some cases (such
as tandem duplications) both the new and old junction can exist in the sample; these pages are how
you check. Reads whose names end in **-M1** or **-M2** mapped better to the new junction.

`seq id`\
Identifiers for the reference sequences involved.

`position`\
Positions of the two sides. Each has an equals sign (=) before or after it representing how the
junction was constructed: the joined pieces approach the given coordinates from the sides with the
equals signs.

`overlap`\
If positive, the number of bp in the junction that could map to either side (generally resolved to
zero by assigning them to one side). If negative, the number of bp unique to the reads crossing the
junction — an insertion relative to the reference.

`reads`\
The total number of reads that map to this junction.

`score`\
The position-hash score in **\<bold angle brackets>** and the minimum-overlap score on the next line.

`skew`\
`neg_log10_pos_hash_p_value`. `NT` in polymorphism mode.

`freq`\
Frequency of the new junction: reads supporting it divided by those supporting it plus those spanning
the original reference sequence at the same breakpoint. `NA` when neither side falls in unique
sequence, so no denominator can be formed.

`range`\
Confidence limits on `freq` — exact (Clopper–Pearson) bounds, taken at the effective depth implied by
how well each read distinguishes the junction from the reference rather than at the raw read count.
This is the interval the frequency cutoffs test, so a junction can be rejected at a frequency above
its cutoff. See [`freq` and `range`](evidence-overview.md#freq-and-range).

`annotation, gene, product`\
Description of the effects of this change on each side. The format is the same as in
[Mutation Display](output.md#mutation-display).

### Junction orientations

<figure>
<img src="../images/jc_side_explanation.png" width="600" />
</figure>

*Figure Credit: Jeff Barrick with additions by Emily Layton*

In the HTML output, equals signs next to the coordinates indicate how the two sides of split reads
supporting a junction are oriented in relation to the reference coordinates joined together in the
sample.

<figure>
<img src="../images/jc_1.png" width="750" />
</figure>

The page reached by the \* link. A partial alignment of reads to the new junction is shown; note the
two joined pieces of reference sequence at the top. This sequence is on the bottom strand of the
reference if start is greater than end. The end of the junction in yellow indicates that end maps to
a repetitive region.

<figure>
<img src="../images/jc_2.png" width="750" />
</figure>

The page reached by one of the ? links. Only a piece of each read maps to this region, ending where
those reads begin matching a disjoint region. The old junction is not supported by any reads and must
no longer exist.

## Options

`--no-junction-prediction`\
Do not predict new sequence junctions.

`--junction-alignment-pair-limit <int>`\
How many alignment pairs to examine when constructing junction candidates. See the tip under
[`Junction Candidates Tested`](#junction-candidates-tested) for when to raise it.

`--junction-minimum-pr-no-read-start-per-position <float>` (default 0.10)\
Floor on the assumed fraction of position–strand combinations with no read starting there, used by
the saturation correction. Set to 0 for pre-0.34.0 behaviour.

Additional junction options appear under `Junction (JC) Evidence Options` in `breseq --help`.

## Worked examples

### An accepted junction

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
JC	182	.	REL606	1	1	REL606	4629812	-1	0	coverage_minus=35	coverage_plus=36
	max_pos_hash_score=200	neg_log10_pos_hash_p_value=0.2	pos_hash_score=53
	side_1_annotate_key=gene	side_1_redundant=0
```

*(Fields trimmed and wrapped for display; a real `.gd` line is one tab-separated row.)*

The two criteria that matter are both visible. `coverage_plus` and `coverage_minus` are large and
nearly equal — the junction is crossed by reads on both strands in similar numbers, which is what
rule 1 of the acceptance criteria demands and what a real join looks like. And
`neg_log10_pos_hash_p_value` is 0.2, far below the 3.0 cutoff: given the depth here, a position-hash
score of 53 out of a possible 200 is entirely unsurprising.

Note the coordinates: position 1 and the last base of REL606. This is the join across the origin of
the circular chromosome — a real junction in the sample, and correctly scored, though not a mutation.

### A rejected junction

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
JC	223	.	REL606	666130	-1	REL606	3697156	-1	0	coverage_minus=8	coverage_plus=8
	max_min_left=21	max_min_left_minus=21	max_min_left_plus=19	max_min_right=0
	max_min_right_minus=0	max_min_right_plus=0	max_pos_hash_score=200
	neg_log10_pos_hash_p_value=3.3	pos_hash_score=13	reject=COVERAGE_EVENNESS_SKEW
	side_1_redundant=1	side_2_redundant=0
```

`neg_log10_pos_hash_p_value` is 3.3, just past the 3.0 cutoff, so the item is rejected for coverage
evenness. Reading the supporting fields shows why that verdict is right rather than marginal.

`pos_hash_score` is 13 out of a possible 200. Sixteen reads support the junction, but they occupy
only thirteen distinct position–strand combinations, so they are piled up rather than spread out —
exactly the pattern the position-hash score exists to detect.

The `max_min_right` fields are the giveaway: **all three are zero**. Not one supporting read reaches
meaningfully into the right side of the junction, on either strand, while `max_min_left` shows they
reach comfortably into the left. A genuine junction is crossed by reads that extend into *both*
sides; these reads all stop at the same place. Combined with `side_1_redundant=1`, marking side 1 as
mapping to several locations, this has the signature of a mapping artifact at a repeat rather than a
real join.

This run rejects only 3 of its 96 junctions. A deeply-sequenced population sample rejects far more —
one of the long tests rejects 145 of 156 — because low-frequency variants and mapping noise both
generate large numbers of weakly-supported candidates.

## See also

- [Evidence overview](evidence-overview.md) — how `JC` compares with `SC`, `DP`, `MP` and `PD`
- [SC: Soft clipping evidence](evidence-sc.md) — the same signal with only one half aligned
- [Mobile element insertions](methods.md#mobile-element-insertions) — how two `JC` items become a `MOB`
- [Read coverage distribution](methods.md#read-coverage-distribution) — the fit the skew score uses
