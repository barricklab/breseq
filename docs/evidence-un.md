# UN: Unknown base evidence

A contiguous stretch of reference positions where there was insufficient evidence to call any base at
all. `UN` records **absence of information**, not a variant. It is on by default and needs no option.

## At a glance

| | |
|---|---|
| Section in `summary.html` | *(none)* |
| Metrics / gates tables | *(none — `UN` has no accept/reject step)* |
| Banner in `index.html` | *(none — `UN` is not shown as unassigned evidence)* |
| Sort order in the GenomeDiff file | 14 |
| Enabled by | on by default |
| Requires | nothing |
| Promotes to | *(nothing — reported as evidence only)* |
| Can be rejected | **no** |

## The signal

At every reference position _breseq_ asks what base the sample carries. Usually the reads answer
clearly. Sometimes there are no reads at all, or too few, or their base qualities are too poor for
any call to be justified. In that case _breseq_ declines to guess and marks the position **unknown**.

Contiguous runs of unknown positions are merged into a single `UN` item covering `start` to `end`.

## Why this matters

`UN` looks like bookkeeping, and it is the most important bookkeeping in the file.

Consider comparing twenty genomes and finding a mutation in three of them. Did the other seventeen
lack the mutation, or did nobody look? Without `UN` those two situations are indistinguishable —
both appear as an absence in the GenomeDiff file. With `UN`, an absent call inside a `UN` region is
explicitly *unexamined*, while an absent call outside one is a genuine negative.

This is why `UN` items typically outnumber every other kind of evidence in a `.gd` file, often by an
order of magnitude, and why `gdtools` comparison operations pay attention to them. When a
cross-sample comparison shows a mutation present in some samples and absent in others, checking
whether the "absent" samples have `UN` coverage there is the first thing to do.

## What this evidence cannot see

`UN` is a statement about data quality, not about the sample. It does not distinguish a region that
is absent from the sample (which would be [`MC`](evidence-mc.md), and probably a deletion) from a
region that simply sequenced badly. A deleted region will usually produce both an `MC` item and `UN`
items covering the same span; the `MC` is the claim that sequence is missing, the `UN` merely records
that no base could be called there.

## From reads to evidence

1. **Walk the pileup.** At every reference position, evaluate whether the aligned reads support any
   base call.
2. **Mark unknown positions.** Positions where no call is justified — typically from having too few
   reads, or reads whose base qualities do not support a confident call — are marked unknown.
3. **Merge runs.** Contiguous unknown positions are combined into one item spanning `start` to `end`.

`UN` items are also produced by `gdtools MASK`, which converts masked regions into unknown-base
evidence so that masked positions are treated as unexamined rather than as confirmed reference.

## Gates

`UN` has none. It is a record of what could not be determined, so there is no hypothesis to test.

## Why items are rejected

They are not. **No `UN` item ever carries a `reject=` field.**

## GenomeDiff fields

Positional fields: `seq_id`, `start`, `end` — see
[GenomeDiff File Format](genomediff-file-format.md#un-unknown-base-evidence). `UN` carries no
statistics, which is why its lines are the shortest in the file.

## In the HTML report

`UN` has no table of its own in `index.html` — with hundreds or thousands of items per run, listing
them would overwhelm the page. The information is instead visible in two places:

- The **coverage plots** linked from `summary.html`, where unknown regions show as gaps.
- The **GenomeDiff file** itself, `output/output.gd` and `output/annotated.gd`, which is where any
  programmatic comparison should read them from.

## Options

None. `UN` prediction is always on and has no tunable parameters of its own, though anything that
changes what counts as adequate coverage — read filtering, alignment stringency, `--require-match-fraction`
— changes which positions end up unknown.

## Worked example

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
UN	301	.	REL606	15616	15616
UN	302	.	REL606	15618	15619
UN	303	.	REL606	15621	15621
```

Those are whole lines. `UN` items carry no `name=value` pairs at all: the claim is simply that across
this span of REL606, no base could be called.

Single bases and pairs of bases, scattered and interleaved with callable positions, are what
dominates a good run — brief patches of poor quality or low depth. Long `UN` runs are more
interesting:

```text title="tests/long_ltee_ara_p1_50k_pe101/expected.gd"
UN	583	.	REL606	2036717	2053723
```

Seventeen kilobases with nothing callable. A stretch that size is not a data-quality blip; it
coincides with an [`MC`](evidence-mc.md) region and a real deletion. The `UN` item does not make that
claim — it only records that no base could be called — but its length is the hint that sends you to
look.

!!! tip "Comparing samples"
    Before concluding that a mutation found in one sample is absent from another, check whether the
    second sample has a `UN` item covering that position. `gdtools` comparison subcommands use `UN`
    evidence for exactly this purpose, reporting an unexamined position differently from a confirmed
    reference one.

## See also

- [Evidence overview](evidence-overview.md) — how the nine types divide the space
- [MC: Missing coverage evidence](evidence-mc.md) — the claim that sequence is genuinely absent
- [GenomeDiff File Format](genomediff-file-format.md) — the `.gd` line specification
