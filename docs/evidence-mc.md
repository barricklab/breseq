# MC: Missing coverage evidence

A region of the reference where read coverage falls away to nothing, or nearly nothing — the
signature of a deletion. `MC` is on by default and is the principal evidence behind large `DEL`
mutations.

## At a glance

| | |
|---|---|
| Section in `summary.html` | `Coverage Values` (when a user coverage cutoff was set) |
| Metrics / gates tables | *(none — `MC` has no accept/reject step)* |
| Banner in `index.html` | `Unassigned missing coverage evidence` |
| Sort order in the GenomeDiff file | 11 |
| Enabled by | on by default; turn off with `--no-missing-coverage-prediction` |
| Requires | nothing |
| Promotes to | `DEL` |
| Can be rejected | **no** |

## The signal

If a stretch of the reference is absent from the sample, no fragment can originate there, and no read
maps there. `MC` walks the pileup looking for exactly that: regions where **unique** read coverage
drops to zero and stays low.

Unlike the split-read and pair-based types, `MC` never looks at how an individual read aligned. It
sees only depth. That makes it robust — a deletion of any size produces the same unambiguous
signature — and it makes the boundaries approximate, because coverage fades out over a read length
rather than stopping at a base.

Non-zero coverage inside a real deletion is common in practice. A little contaminating DNA from
another sample, or a few high-error reads mapping spuriously, will put a handful of reads in a region
that is genuinely absent. This is why the boundary test uses a coverage *threshold* rather than
requiring literal zero.

## Read coverage distribution

`MC`'s threshold is not a fixed number. It is derived from the coverage distribution fitted to this
run, which is described in full under
[Read coverage distribution](methods.md#read-coverage-distribution) — the same fit that
[`CN`](evidence-cn.md) and the `summary.html` coverage tables use.

In brief: read depth across a genome is overdispersed relative to the Poisson distribution an
idealised random shearing would give, so _breseq_ fits a negative binomial to the depth observed at
unique-only positions, using left-censored data so that genuinely deleted regions do not drag the fit
down. From that fit it calculates the depth below which coverage is surprising.

If the fit fails — very low coverage, a heavily biased library, targeted sequencing — _breseq_ falls
back to a rougher estimate, or may call an entire reference sequence deleted. Two options exist for
this case: `-c`/`--contig-reference` fits all sequences in a file together, for a draft genome whose
contigs are individually too short to fit; `-t`/`--targeted-sequencing` tells _breseq_ to call
mutations regardless of what the coverage distribution looks like, and suppresses deletion calling.

## Seed and extend

Deletion predictions are initiated at every reference position with unique-only coverage of **zero**.
They are extended in each direction, and merged with neighbours, until unique coverage exceeds a
threshold calculated from the fitted distribution for that reference sequence. That cutoff is the
minimum threshold coverage *t* satisfying

$F(t) > 0.05\times\sqrt{L}$,

where *F* is the negative binomial cumulative distribution function with the best-fit mean and size
parameters, and *L* is the reference sequence length. The $\sqrt{L}$ term is a multiple-testing
correction: a longer sequence offers more chances for a low-coverage window to arise by accident, so
it demands a lower threshold.

## Ambiguous boundaries in repeats

A deletion boundary that falls inside a repeat is genuinely ambiguous, and `MC` says so rather than
guessing.

The problem: even when one copy of a repeat is deleted, reads from its surviving copies elsewhere in
the genome still map there, so coverage does not drop. Whether the repeat went with the deletion is
unknowable from depth alone.

_breseq_ handles this in two parts:

- A region of repeat coverage lying **wholly inside** a low-coverage region is assumed deleted along
  with its flanks.
- A region of repeat coverage **overlapping one end** of the prediction makes that end a *range*
  rather than a point. The two limits are the two extreme readings: that the entire contiguous
  repetitive region is missing, and that it is entirely still there. The far limit is found by
  re-running the same extension algorithm on unique coverage *plus normalised repeat coverage*, where
  a read matching *n* places contributes 1/*n* of a read's depth to each.

The `start_range` and `end_range` fields carry the resulting uncertainty; both are zero when the
boundary is unambiguous.

<figure>
<img src="../images/region_coverage_example.png" class="align-center" width="600" height="333" alt="Coverage in a deleted reference region." /><figcaption aria-hidden="true"><strong>Coverage in a deleted reference region.</strong></figcaption>
</figure>

This example shows a region of missing coverage (white background) that extends into a region of
repeat coverage (red line), making the left side end of the missing coverage ambiguous.

## What this evidence cannot see

- **The other side of the event.** `MC` says sequence is absent; it does not say what the two
  surviving flanks are now joined to. That is [`JC`](evidence-jc.md)'s job, and a deletion supported
  by both is far better characterised than one supported by `MC` alone.
- **Base-pair boundaries.** Even outside repeats, the edges are placed where coverage recovers, not
  where the cut was made.
- **Deletions of repeated sequence.** If every copy of a repeat family is present elsewhere, deleting
  one leaves coverage unchanged.
- **Copy number above one.** `MC` only looks downwards. Amplifications are [`CN`](evidence-cn.md).

## Gates

`MC` has none, and prints no gates table. The coverage threshold described above governs where a
region *ends*, not whether it is reported.

The `Coverage Values` table in `summary.html` appears when a user coverage cutoff was set, and lists
per reference sequence the `Calculated average`, `Calculated propagation cutoff`,
`Calculated seed cutoff`, `User defined propagation cutoff` and `User defined seed cutoff`.

## Why items are rejected

They are not. **No `MC` item ever carries a `reject=` field.** `MC` reports the coverage it measured;
the judgement about whether that constitutes a deletion happens later, during
[mutation prediction](methods.md#large-deletions).

## GenomeDiff fields

Positional fields: `seq_id`, `start`, `end`, `start_range`, `end_range` — see
[GenomeDiff File Format](genomediff-file-format.md#mc-missing-coverage-evidence). The region of
missing coverage lies between `[start, start+start_range]` and `[end-end_range, end]`.

Notable `name=value` pairs:

*   **left_outside_cov**, **left_inside_cov** — unique coverage at the last position *outside* the
    region and the first position *inside* it, on the left margin. The contrast between them is the
    evidence: a sharp drop is a clean boundary.
*   **right_inside_cov**, **right_outside_cov** — the same for the right margin.
*   **gene_name**, **gene_product**, **locus_tag** — annotation spanned. For a long region these are
    ranges, and genes only partly covered appear in square brackets.

## In the HTML report

Accepted items appear on `index.html` under `Unassigned missing coverage evidence`, or attached to
the `DEL` mutation they support.

`* links`\
Links to results pages showing the alignment of reads to the left and right margins of the region
with missing coverage.

`÷ link`\
Link to the results page showing a plot of the read coverage in the region of the missing coverage.

`seq id`\
Identifier for the reference sequence where the change is located.

`start, end, size`\
The start and end reference positions and size of the missing coverage. May indicate a range of
positions when one end of the missing coverage is in a repeat region.

`← cov`\
Unique read coverage depth on the left margin. Coverage at the last position outside the region is
shown followed by coverage at the first position inside the region, in brackets.

`→ cov`\
Unique read coverage depth on the right margin. Coverage at the last position inside the region is
shown followed by coverage at the first position outside it.

`gene, description`\
Description of the change's effects for each side. The format of these columns is the same as in
[Mutation Display](output.md#mutation-display).

<figure>
<img src="../images/mc_1.png" width="750" />
</figure>

Read coverage depth around the missing coverage. The white area shows the maximal boundaries of the
predicted range.

The graphed lines are labeled "unique" for reads with only one best match to the reference genome and
"repeat" for multiple equally good matches to repeat sequences (which are down-weighted by how many
matches they have, i.e. a read matching three places contributes 1/3 to the coverage depth at each
matched site). Within each type coverage is graphed separately for reads mapping to the "top" and
"bottom" strands of the reference sequence (i.e., forward and reverse complement matches) to aid in
detecting artifacts, and these sum to the "total" coverage value.

## Options

`MC` is always on and has no `--predict-*` flag. Two options change how the coverage distribution it
depends on is fitted:

`-c <file_path>, --contig-reference <file_path>`\
Fit the coverage distribution of all sequences in this reference file together, as if they were one
chromosome. Use for a draft genome whose contigs are individually too short to fit.

`-t, --targeted-sequencing`\
Call mutations regardless of what the coverage distribution looks like. Deletion mutations are not
called in this case, since low coverage is expected everywhere off-target.

## Worked examples

### A clean deletion boundary

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
MC	166	.	REL606	16975	18001	0	0	gene_name=[nhaA]	gene_product=[nhaA]
	left_inside_cov=2	left_outside_cov=72	locus_tag=[ECB_00018]
	right_inside_cov=1	right_outside_cov=64
```

*(Fields wrapped for display; a real `.gd` line is one tab-separated row.)*

Both `start_range` and `end_range` are **zero**, so neither boundary is ambiguous — no repeat overlaps
either end and the coordinates given are the best available.

The margin coverages are the evidence. On the left, coverage goes from `left_outside_cov=72` to
`left_inside_cov=2` in a single base; on the right it recovers from 1 to 64. That is a coverage cliff,
not a gradual dip. The small non-zero inside values are exactly the residue described above —
a couple of stray or contaminating reads — and are the reason the algorithm extends to a fitted
threshold rather than requiring literal zero.

Compare this region with the [`CN`](evidence-cn.md) item covering the same span, which reports
`copy_number=0` over tiles 17001–18000. The two agree, from independent directions: `CN` from
normalised tile depth, `MC` from the seed-and-extend walk. `MC` gives the sharper boundary.

### An ambiguous boundary in a repeat

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
MC	168	.	REL606	547082	555824	618	0	gene_name=[insB-6]–[ECB_00513]
	left_inside_cov=21	left_outside_cov=22	right_inside_cov=0	right_outside_cov=63
```

`start_range` is 618 while `end_range` is zero, so the two ends are known to very different
precision. The right boundary is sharp — coverage recovers from 0 to 63. The left boundary is not: it
lies somewhere in a 618-base window, and `left_inside_cov=21` against `left_outside_cov=22` shows
why. There is *no coverage cliff on the left at all*.

That is the repeat signature. The gene names confirm it: the region begins in `insB-6`, an IS element
present in many copies across REL606. Reads from the other copies map here regardless of whether this
one was deleted, so depth cannot say where the deletion starts. The 618-base range spans the two
honest extremes — the whole repeat went, or none of it did.

This is the case where `MC` most wants corroboration. A [`JC`](evidence-jc.md) item landing inside
that window would fix the boundary to the base and identify what the surviving flank is joined to;
IS-mediated deletions in this organism usually produce one.

!!! tip "Reading margin coverage"
    The four `*_cov` fields are the quickest sanity check on an `MC` item. A large outside value with
    a near-zero inside value on both margins is a clean deletion. Comparable values across a margin
    mean that boundary is resting on repeat coverage, and the corresponding `_range` field should be
    non-zero to match.

## See also

- [Evidence overview](evidence-overview.md) — how the nine types divide the space
- [Read coverage distribution](methods.md#read-coverage-distribution) — the fit this depends on
- [CN: Copy number evidence](evidence-cn.md) — the tile-based view of the same depth
- [JC: New junction evidence](evidence-jc.md) — what fixes an ambiguous boundary
- [Large deletions](methods.md#large-deletions) — how `MC` becomes a `DEL` mutation
