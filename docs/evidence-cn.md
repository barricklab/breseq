# CN: Copy number evidence

A run of reference tiles whose read depth departs from the genome-wide average by enough to imply a
different number of copies in the sample. `CN` is the only evidence type built purely from **coverage
depth**, with no reference to how any individual read aligned. It is experimental and off by default;
enable it with `--predict-copy-number`.

## At a glance

| | |
|---|---|
| Section in `summary.html` | `Copy Number Variation` |
| Metrics / gates tables | *(none — `CN` has no accept/reject step)* |
| Banner in `index.html` | `Unassigned copy number evidence` |
| Sort order in the GenomeDiff file | 13 |
| Enabled by | `--predict-copy-number` |
| Requires | nothing (works on single-end data) |
| Promotes to | *(nothing — reported as evidence only)* |
| Can be rejected | **no** |

## The signal

If a region of the genome is present in two copies rather than one, twice as many fragments originate
there, and twice as many reads map there. If it is deleted, none do. Copy number is therefore
readable directly off the depth of coverage, without needing a single read to align unusually.

That makes `CN` complementary to every other evidence type in this section. [`JC`](evidence-jc.md),
[`SC`](evidence-sc.md), [`DP`](evidence-dp.md) and [`PD`](evidence-pd.md) all detect the *boundary*
of an event from reads that cross it. `CN` detects the *interior* — it says nothing about where the
boundaries are, but it can see an amplification whose junctions are all inside repeats where no read
maps uniquely.

The reference is divided into fixed-width tiles, the coverage in each is corrected for known biases
and normalised against the genome average, and runs of adjacent tiles sharing a copy number are
merged into a single `CN` item.

## What this evidence cannot see

- **Boundaries.** A `CN` item's `start` and `end` are tile boundaries, not breakpoints. They are
  accurate to the tile size, which is reported on every item.
- **Anything smaller than a tile.** Sub-tile events are invisible.
- **What caused the change.** A doubled region might be a tandem amplification, a duplication
  elsewhere, or a plasmid; coverage alone cannot say. Pair it with `JC` or `DP` to find out.
- **Copy number in a region with no unique coverage.** Repeats accumulate reads from all their
  copies, so depth there does not report the copy number of any one of them.

## Coverage corrections

Raw depth is not proportional to copy number, because two systematic biases act on it first. `CN`
corrects for both, and `summary.html` reports how much each correction achieved.

**GC bias.** Library preparation and amplification favour some base compositions over others, so
depth varies with the GC content of the local sequence independently of copy number.

**Ori-ter replication bias.** In a growing bacterial culture, cells are caught mid-replication, so
sequence near the replication origin is present in more copies *on average across the population*
than sequence near the terminus. The result is a smooth gradient of coverage around the chromosome —
biological, entirely expected, and nothing to do with mutation. Left uncorrected it would make every
region near the origin look amplified.

The `Coverage spread after each correction` table in `summary.html` shows the spread of coverage
uncorrected, after GC correction, and after GC + ori-ter correction. A correction that is working
reduces the spread. The `Ori-ter replication bias fit` table reports the fitted origin and terminus,
the coverage at each, and the peak-to-trough ratio — which is itself a useful measure of how fast the
culture was growing.

Where a fit is not possible the tables say so rather than silently proceeding, with cells reading
`no usable coverage: <reason>`, `no ori-ter bias detected`, or a note beginning `No OTR correction`.

## From reads to evidence

1. **Tile the reference.** Divide each sequence into fixed-width tiles and compute mean coverage in
   each.
2. **Correct for GC bias.** Fit and remove the dependence of coverage on local base composition.
3. **Correct for ori-ter bias.** Fit the replication gradient and remove it.
4. **Normalise.** Divide by the genome-wide average to get `relative_coverage`, where 1.0 is one
   copy.
5. **Assign copy number.** Round the normalised coverage to an integer copy number per tile.
6. **Merge runs.** Combine adjacent tiles sharing a copy number into one item, spanning from the
   first tile's start to the last tile's end.

## Gates

`CN` has none, and prints no gates table. It reports the coverage it measured; there is no hypothesis
being tested, so there is nothing to reject. The relevant diagnostics are the correction-quality
tables in the `Copy Number Variation` section of `summary.html` described above.

## Why items are rejected

They are not. **No `CN` item ever carries a `reject=` field.** If you are looking for why a
copy-number call is absent or implausible, the answer is in the coverage corrections, not in a
rejection reason — check the `Coverage spread after each correction` table first.

## GenomeDiff fields

Positional fields: `seq_id`, `start`, `end`, `copy_number` — see
[GenomeDiff File Format](genomediff-file-format.md#cn-copy-number-evidence).

Notable `name=value` pairs:

*   **relative_coverage** — the corrected, normalised coverage over the item's tiles, where 1.0 means
    one copy. This is the underlying measurement; `copy_number` is it rounded.
*   **tile_size** — the width of the tiles the item was built from, and therefore the resolution of
    its `start` and `end` coordinates.
*   **gene_name**, **gene_product**, **locus_tag** — annotation spanned by the region. For a long
    item these are ranges, and genes only partly covered are shown in square brackets.

## In the HTML report

Accepted items appear on `index.html` under `Unassigned copy number evidence`. `summary.html` carries
the fitting diagnostics and per-sequence copy-number plots.

`seq id`\
Identifier for the reference sequence.

`start`, `end`\
The region's bounds, accurate to `tile size`.

`tile size`\
The tile width used, and hence the coordinate resolution.

`copy number`\
The integer copy number assigned. `0` means the region appears deleted.

`rel cov`\
`relative_coverage` — the corrected coverage relative to the genome average, before rounding. Worth
reading alongside `copy number`: a value of 1.85 rounded to 2 is a different quality of evidence from
one of 1.99.

`gene`, `product`\
Annotation spanned by the region.

## Options

!!! warning "Experimental"
    `CN` prediction is experimental, off by default, and marked HIGHLY EXPERIMENTAL in the command
    line help.

`--predict-copy-number`\
Predict copy number variation evidence.

## Worked examples

### A deleted region

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
CN	278	.	REL606	17001	18000	0	gene_name=[nhaA]	gene_product=[nhaA]
	locus_tag=[ECB_00018]	relative_coverage=0	tile_size=100
```

*(Fields wrapped for display; a real `.gd` line is one tab-separated row.)*

`relative_coverage` is zero — not a single read maps across this kilobase — so the copy number is
zero. The square brackets on `[nhaA]` mean the gene is only partly contained in the region.

Note the coordinates: they are round numbers because they are tile boundaries, and `tile_size` says
they are good to 100 bases. `CN` locates the *interior* of the deletion confidently and its edges
only approximately. To place the breakpoints you want [`MC`](evidence-mc.md), which gives explicit
start and end ranges, or [`JC`](evidence-jc.md) for the exact base.

### A duplicated region

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
CN	282	.	REL606	1969801	1970100	2	gene_name=leuZ–[glyW]
	gene_product=leuZ,cysT,[glyW]	locus_tag=[ECB_t00029]–[ECB_t00031]
	relative_coverage=1.85	tile_size=100
```

Here `relative_coverage` is close to, but not exactly, twice the genome average, and rounds to a copy
number of 2. The gap between 1.85 and 2.0 is ordinary: coverage corrections are imperfect, and the
three-tile region is short enough that sampling noise matters.

The region covers a cluster of tRNA genes. That is worth noticing, because tRNA and rRNA operons are
repetitive, and repeats are exactly where depth-based copy number needs the most care — reads from
several near-identical copies compete for placement, and the coverage attributed to any one of them
is not a clean measurement of its abundance. Corroborate a copy-number call in a repeat region with a
junction before believing it.

!!! tip "Reading an implausible copy-number profile"
    If a report shows broad regions of elevated copy number near one point of the chromosome and
    depressed copy number opposite it, the ori-ter correction has probably failed. Check the
    `Ori-ter replication bias fit` table for a `no ori-ter bias detected` or `No OTR correction`
    note, and the `Coverage spread after each correction` table to see whether the corrections
    actually reduced the spread.

## See also

- [Evidence overview](evidence-overview.md) — how the nine types divide the space
- [MC: Missing coverage evidence](evidence-mc.md) — the boundary-aware view of a deletion
- [JC: New junction evidence](evidence-jc.md) — what identifies the cause of a copy-number change
