# Evidence Types

_breseq_ works in two stages. First it gathers **evidence**: local, mechanical observations about how
reads aligned to the reference, each one an assertion that something at a particular place does not
look the way the reference says it should. Then it tries to explain that evidence with **mutational
events** — the `SNP`, `DEL`, `MOB`, `AMP` and friends that appear at the top of `output/index.html`.

This section documents the first stage. There are nine evidence types, and each sees a different
physical signature in the data:

| | Type | Signal | Paired reads only | On by default | Turn off with |
|---|---|---|---|---|---|
| 10 | [RA: Read alignment](evidence-ra.md) | Reads in a pileup disagree with the reference base | no | yes | `--no-read-alignment-prediction` |
| 11 | [MC: Missing coverage](evidence-mc.md) | A stretch of reference has no reads on it | no | yes | `--no-missing-coverage-prediction` |
| 12 | [JC: New junction](evidence-jc.md) | A read splits across two disjoint reference locations | no | yes | `--no-junction-prediction` |
| 13 | [CN: Copy number](evidence-cn.md) | A window's read depth departs from the genome average | no | yes | `--no-copy-number-prediction` |
| 14 | [UN: Unknown base](evidence-un.md) | Too little data to call a base either way | no | yes | *(always on)* |
| 15 | [SC: Soft clipping](evidence-sc.md) | Reads stop aligning part-way and agree on what follows | no | **no** | opt in with `--predict-soft-clipping` |
| 16 | [DP: Discordant pair](evidence-dp.md) | Individual read pairs map in the wrong orientation, too far apart, or on two sequences | yes | yes | `--no-discordant-pair-prediction` |
| 17 | [MP: Missing pair](evidence-mp.md) | Reads pile up facing a point, their mates mapping nowhere | yes | yes | `--no-missing-pair-prediction` |
| 18 | [PD: Pair distance](evidence-pd.md) | Read pairs spanning a point are collectively shifted in mapping distance | yes | yes | `--no-pair-distance-prediction` |

The leading number is the type's **sort order**. It fixes the order evidence appears in a GenomeDiff
file, in `output/index.html`, and in this section's navigation, so all three agree.

**Everything except `SC` is predicted by default.** `DP`, `MP` and `PD` additionally need paired-end
reads, so they are skipped on a single-end run or under `--no-paired-mapping`, and `CN` needs the
separate [CNery](https://github.com/barricklab/CNery) program on your `PATH`.

!!! note "The `--predict-*` flags for CN, DP, MP and PD are deprecated"
    `--predict-copy-number`, `--predict-discordant-pairs`, `--predict-missing-pairs` and
    `--predict-pair-distance` are still accepted so existing command lines keep working, and are
    hidden from `breseq --help`, but they no longer switch anything on — these four types are on
    already. One has a lingering effect: passing `--predict-copy-number` explicitly makes a missing
    CNery program **fatal** rather than a warning.

    `--predict-soft-clipping` is not deprecated. `SC` is still opt-in, and deliberately so: that flag
    also lowers `--require-match-fraction` from 0.9 to 0.5, which changes which alignments are
    accepted genome-wide and so would alter every other evidence type if it were on by default.

## Accepted, rejected, ignored

Every evidence item _breseq_ constructs ends in one of three states. The distinction matters because
two of them are invisible in different ways.

**Accepted.** The item passed every gate. It appears in `output/index.html` — either attached to a
mutation prediction, or in an "unassigned evidence" table if nothing consumed it.

**Rejected.** The item failed at least one gate. It carries a `reject=` field naming every gate it
failed, and it appears on `output/marginal.html` with a row reading *Rejected: &lt;reason&gt;*.
Rejection is a judgement about the data.

!!! note "Only the best rejected items are shown"
    `marginal.html` shows a capped number of rejected items per type — by default the top 20 for
    `RA`, `SC`, `DP`, `MP` and `PD`, and the top 10 for `JC`, each ranked by the sort the table
    header names (frequency, score, or skew). On a run with hundreds of rejected items, the
    `.gd` file is the only complete record. The gates table's `outcome` row always tallies *all* of
    them, so the counts there will exceed what the page displays.

**Ignored.** The item carries an `ignore=` field and is dropped from the report altogether. This is
*not* a judgement about the data — it means the item is an artifact of the reference's shape rather
than of the sample:

| `ignore` value | Meaning |
|---|---|
| `CIRCULAR_CHROMOSOME` | The item spans the origin of a circular sequence, where position 1 and position *L* are physically adjacent. Applies to `JC` and `DP`. |
| `CONTIG_END` | The item sits within a read length or so of the end of a reference sequence, where its statistic stops being measurable. Applies to `RA`, `MC`, `JC`, `SC`, `MP` and `PD`. |
| `masked` | The item falls in a region masked out by `gdtools MASK`. |

`DP` is deliberately exempt from the contig-end rule. A cross-sequence `DP` — a plasmid integration,
or a translocation between two contigs — sits at two contig ends *by construction*, so applying the
rule would silently discard every translocation call.

## "Unassigned" evidence

An accepted item that no mutation prediction consumed is shown on `output/index.html` under its own
banner. The banners, exactly as they appear on the page:

- Unassigned missing coverage evidence
- Unassigned new junction evidence
- Unassigned copy number evidence
- Unassigned soft clipping evidence
- Unassigned discordant pair evidence
- Unassigned missing pair evidence
- Unassigned pair distance evidence

Unassigned is not a failure state. It usually means the evidence is real but _breseq_ could not
work out a single mutational event that explains it — a junction into a repeat family, say, where
the insertion point is clear but the donor copy is not. Unassigned evidence is often the most
interesting material in a report, and is where manual curation starts.

## Which type sees which structural variant

Five types can report structural variation, and they divide the space by what each one physically
requires of a read. This table is the fastest way to work out why an event you know is real was
found by one type and missed by another.

| | [JC](evidence-jc.md) | [SC](evidence-sc.md) | [DP](evidence-dp.md) | [PD](evidence-pd.md) | [MP](evidence-mp.md) |
|---|---|---|---|---|---|
| Needs paired reads | no | no | **yes** | **yes** | **yes** |
| Needs a read split across the breakpoint, both halves mapping | **yes** | no | no | no | no |
| Needs both mates to map | n/a | n/a | **yes** | **yes** | no |
| Needs the new sequence to exist somewhere in the reference | **yes** | **yes** | **yes** | no | no |
| Locates the breakpoint to base-pair resolution | **yes** | **yes** | no | no | no |
| Reports the size of the event | **yes** | no | no | **yes** | no |
| Names the sequence that was inserted | **yes** | partly | **yes** | no | no |
| One-sided (a real event yields two items) | no | **yes** | no | no | **yes** |

Read down the columns rather than across:

- **`JC` is the gold standard and the most demanding.** It needs a single read that spans the
  breakpoint with enough sequence on each side to align both halves. When it fires you get the exact
  base and both sides' identity, which is why it — and only it — promotes to `MOB`.
- **`SC` is `JC` with one half missing.** The read stops aligning part-way; the clipped tail is real
  sequence but is not itself aligned. So you get the exact base but not the other side.
- **`DP` tests each pair on its own.** A pair is discordant when its orientation is wrong, when it
  maps farther apart than the library's cutoff, or when its mates land on two different reference
  sequences. Because it is per-pair, `DP` sees large events and translocations, and it is the only
  type whose two sides carry independent sequence identifiers.
- **`PD` tests the pairs collectively.** It asks whether the pairs whose unsequenced middle spans a
  point are *as a population* shifted longer (sequence missing) or shorter (sequence added). No
  individual pair has to be an outlier, which is exactly the regime `DP` cannot see: events of a few
  hundred bases, in either direction.
- **`MP` is the only type that can see sequence absent from the reference entirely.** A fragment
  crossing into a novel insert puts one mate in reference sequence and the other in the insert, where
  it maps nowhere. `JC` cannot see this (no second half to align), `DP` cannot (needs both mates
  mapped), `SC` cannot say the tail is novel rather than displaced.

### An event of size *S* — who sees it?

Assume a paired library with a few-hundred-base fragment size and all predictors enabled.

- **A single base, or a few.** `RA` alone. None of the structural types apply.
- **Tens of bases.** `RA` for the smallest, `JC` once the event is large enough that reads spanning
  it stop aligning through.
- **A few hundred bases.** `JC` if a read happens to bridge it. Otherwise `PD` — this is the band
  `PD` exists for, and `DP` is blind to it in both directions because no individual pair is unusual.
- **Kilobases.** `JC` at the breakpoint, `DP` from pairs that now straddle it, `PD` from the
  collective shift, and `MC` if the event is a deletion. Convergence of several types on one place is
  the strongest signal a report offers.
- **A novel insertion of any size.** `MP` at each shoulder, plus `SC` if the insert's first bases
  happen to be readable. Nothing else can see it at all.

## Shared vocabulary

Several columns mean the same thing across types. They are defined once here.

### `score`

`SC`, `MP` and `PD` are each decided by a genome-wide score of the same form: **minus the log10 of
the number of items this good expected anywhere in the reference by chance.** So a score of 0 means
one such item is expected per genome, and 3 means one per thousand genomes. Each type's cutoff
defaults to 3.

The important property is that the null is **fitted to the run**, not assumed. _breseq_ measures the
background rate of the relevant artifact across the whole reference — how often mates fail to align,
how often reads are clipped — along with how unevenly that rate is spread, and tests each candidate
against that. A fixed frequency cutoff cannot substitute, because it implicitly assumes the
background is zero, which holds in a simulation and in nothing else.

This is also why the same command line gives wildly different counts on different libraries, and why
each of these types reports its fitted null in a metrics table on `summary.html`. If a run predicts
implausibly many or few, that table is the first thing to read.

`RA`'s `score` is the same *form* — an E-value scaled by the number of reference positions — but the
statistic underneath is a Bayesian model comparison rather than a background rate. See
[RA](evidence-ra.md).

### `freq` and `range`

`freq` is a point estimate of what fraction of the sample carries the variant. `range` is a 95%
confidence interval around it.

**The frequency cutoffs test the interval, not the point estimate.** An item can therefore be
rejected at a frequency that reads as comfortably above its cutoff, because it is the *bound* that
falls below. This is the single most confusing thing in a _breseq_ report, and it is deliberate: a
variant is never rejected merely for having shallow coverage, only for being confidently below the
cutoff.

Which bound is tested depends on what is being asked. In consensus mode the lower bound is tested
("is this confidently the majority allele?"); in polymorphism mode the upper bound is tested for
consensus calls ("can a fixed interpretation still not be ruled out?"). For `RA`, the interval comes
from the fitted allele model, so it widens for low coverage or poor base quality rather than tracking
read count alone. For `JC`, `SC`, `DP`, `MP` and `PD` they are exact (Clopper–Pearson) bounds on the
underlying read or pair counts.

### Distinct-position counts

`SC`, `DP`, `MP` and `PD` all count, alongside their supporting reads, how many **distinct** outer
coordinates or fragment ends those reads occupy (`distinct_read_count`, `distinct_pair_count`,
`distinct_discordant_count`). PCR amplifies one molecule into many identical reads, so a hundred
reads starting at the same coordinate are one observation, not a hundred. Requiring several distinct
positions stops duplicates of a single molecule from carrying a prediction on their own.

The `redundant` flag is related but different: it marks an item whose supporting reads mapped to more
than one place in the reference, so the position shown is one example among several.

## Metrics and gates

Four types — `DP`, `PD`, `SC`, `MP` — print a pair of tables in `output/summary.html`, because their
thresholds are largely derived from the run rather than fixed on the command line, so the options
alone do not explain why a run predicted few or many items.

The two tables separate **what was measured** from **what was decided**:

- *Type* **evidence metrics** — quantities measured or fitted from this run, plus the parameters that
  define them. Nothing here accepts or rejects anything; these are the numbers the decisions are
  taken against.
- *Type* **evidence gates** — the decisions. **Every gate is a rule applied to a metric.**

### The metrics table

`metric`\
What was measured.

`value`\
Its value in this run.

`basis`\
How it was arrived at — what was counted, over what, and which options shaped it. This column is
where the fitted nulls are reported, and it is what makes the table worth reading.

### The gates table

`gate`\
The name of the test, in the order it is applied.

`rule`\
The threshold actually in force. `OFF` means the gate was disabled.

`rejects as`\
**The `reject=` value an item picks up when it fails this gate.** This is the link back from a
*Rejected: …* line in the report to the threshold that produced it, and to the option that sets that
threshold.

`basis`\
Which metric the rule is applied to, and why.

The last row of every gates table is `outcome`, which tallies how many items were examined, accepted,
and rejected by each gate. Those counts reconcile with the evidence tables elsewhere in the report —
and note they count *all* items, including rejected ones the report does not display.

The other five types print no such tables. `RA` and `JC` instead print their settings in the
`Read Alignment Evidence` and `New Junction Evidence` sections of `summary.html`, and `MC`, `CN` and
`UN` have no accept/reject step to report at all.

## All rejection reasons

Every reason any evidence type can carry, with the sentence the report renders for it. Each type's
own page explains what its subset means for that type specifically — the same reason can be
testing quite different things depending on where it appears.

| `reject` value | Types | Shown in the report as |
|---|---|---|
| `SCORE_CUTOFF` | RA, SC | E-value score below prediction cutoff. |
| `FREQUENCY_CUTOFF` | RA, JC | Frequency below/above cutoff threshold. |
| `FREQUENCY_BELOW_CUTOFF` | SC | Lower confidence bound on frequency below cutoff threshold. |
| `FISHER_STRAND` | RA, SC | Biased read strand distribution supporting prediction. |
| `KS_BASE_QUALITY` | RA | Biased base quality scores supporting prediction. |
| `VARIANT_COVERAGE` | RA | Variant not supported by required number of total reads. |
| `TOTAL_COVERAGE` | RA | Genome position does not have required minimum number of aligned reads. |
| `VARIANT_STRAND_COVERAGE` | RA | Variant not supported by required number of reads on each strand. |
| `TOTAL_STRAND_COVERAGE` | RA | Genome position does not have required minimum number of aligned reads on each strand. |
| `INDEL_HOMOPOLYMER` | RA | Polymorphic indel expands or contracts a homopolymer stretch. |
| `SURROUNDING_HOMOPOLYMER` | RA | Polymorphic base substitution creates a homopolymer stretch. |
| `POLYMORPHIC_INDEL` | RA, JC | Indel polymorphism suppressed by --polymorphism-no-indels. |
| `COVERAGE_EVENNESS_SKEW` | JC | Coverage evenness skew score above cutoff. |
| `BETWEEN_TWO_JUNCTION_ONLY_SEQUENCES` | JC | Between two junction-only reference sequences. |
| `CLIPPED_TAIL_CONSENSUS` | SC | Clipped read tails do not agree on a consensus sequence. |
| `LOW_COMPLEXITY_TAIL` | SC | Clipped read tails are a homopolymer or nearly one base; typical of dark-cycle (poly-G) or adapter read-through, not of donor sequence. |
| `NEARBY_BETTER_SOFT_CLIPPING` | SC | Stronger soft-clipping evidence in the same direction within a few bases; probably the same breakpoint. |
| `DISCORDANT_PAIR_FREQUENCY` | DP | Lower confidence bound on discordant pair frequency below cutoff. |
| `CONCORDANT_PAIR_SKEW` | DP | Concordant pair skew score above cutoff. |
| `MISSING_PAIR_SCORE` | MP | Missing pair score below the genome-wide false-positive cutoff. |
| `MISSING_PAIR_COUNT` | MP | Too few reads with unmapped mates support this position. |
| `MISSING_PAIR_DUPLICATES` | MP | Supporting reads with unmapped mates start at too few distinct positions. |
| `MISSING_PAIR_FREQUENCY` | MP | Missing pair local frequency below cutoff. |
| `NEARBY_BETTER_MISSING_PAIR` | MP | A better-scoring missing pair prediction lies within one paired-mapping distance. |
| `PAIR_DISTANCE_SCORE` | PD | Pair distance score below the genome-wide false-positive cutoff. |
| `PAIR_DISTANCE_COUNT` | PD | Too few read pairs with a shifted mapping distance support this position. |
| `PAIR_DISTANCE_DUPLICATES` | PD | Supporting read pairs span too few distinct fragment ends. |
| `PAIR_DISTANCE_FREQUENCY` | PD | Pair distance local frequency below cutoff. |
| `PAIR_DISTANCE_SIZE` | PD | Size of the mapping distance shift is not distinguishable from zero. |
| `PAIR_DISTANCE_INCONSISTENT` | PD | Supporting read pairs cannot all be spanning an event of the estimated size at one position. |

`MC`, `CN` and `UN` are never rejected. They are descriptions of the coverage _breseq_ observed, not
claims that need testing, so there is nothing for a gate to reject.

## Where to go next

- To read a specific evidence table in the HTML report, start from [Output](output.md).
- For the exact field layout of a `.gd` line, see
  [GenomeDiff File Format](genomediff-file-format.md).
- For how evidence becomes a mutation prediction, see
  [Mutation prediction](methods.md#mutation-prediction).
