_breseq_ (pronounced: \\brēz-ˈsēk\\ or *breeze-seq*) is a
computational pipeline specialized for the analysis of short-read
re-sequencing data (e.g. Illumina, 454, IonTorrent, etc.). It uses
reference-based alignment approaches to predict mutations in a sample
relative to an already sequenced genome. _breseq_ is intended for
microbial genomes (\<10 Mb) and re-sequenced samples that are only
slightly diverged from the reference sequence (\<1 mutation per 1000
bp).

_breseq_'s primary advantages over other software programs are that
it can:

1.  Accurately predict new sequence junctions, such as those associated
    with mobile element insertions.
2.  Integrate multiple sources of evidence for genetic changes into
    mutation predictions.
3.  Produce annotated output describing biologically relevant mutational
    events.

_breseq_ version 0.38+ can also be used with long-read sequencing
data (e.g., Nanopore, PacBio). It analyzes long reads by automatically
splitting them into shorter subsequences so that it can use its standard
mapping and analysis methods. By default, pieces cut from the same long read
that are 2000 bases apart are also paired with one another, so that the
discordant pair (DP), missing pair (MP) and pair distance (PD) evidence
normally derived from paired-end libraries can be predicted from long reads
(see `--long-read-pair-distance`; set it to 0 to turn this off). Be aware
that this still discards unique information in long reads that can be used
to identify large-scale structural variation in a genome! We recommend that
you perform _de novo_ assembly and whole-genome comparisons with other
software programs to analyze long-read data for structural variation.

If your Nanopore reads were processed using fast basecalling, we recommend adding
the `--nanopore-fast` option. It rules out false-positive mutation predictions
caused by the high homopolymer indel error rates of this mode: no indel mutations
will be called in homopolymer repeats of 4 or more bases. This option is not
necessary for data processed in high-accuracy basecalling mode.

_breseq_ records the type of reads in each read file: `illumina`,
`nanopore-fast`, `nanopore-hac` (high or super accuracy basecalling),
`pacbio-clr` or `pacbio-hifi`. It detects the type from the read names
when they identify the platform or the basecaller model, and otherwise
guesses it from the read lengths and base qualities. The type is shown in
the read file table of the summary page (marked when it was a guess) and
as the platform (`PL`) of the read group in the output BAM file. To set
the type yourself, prefix the read file with it on the command line:

```
breseq -r reference.gbk nanopore-hac:long_reads.fastq.gz illumina:short_reads_R1.fastq.gz short_reads_R2.fastq.gz
```

The names (and synonyms such as `ont-hac` and `nanopore-sup`) are the ones
[brefito](https://github.com/barricklab/brefito) uses. At present the read
type is only recorded. It does not change how mutations are predicted, and
`--nanopore-fast` is not turned on automatically for `nanopore-fast` reads;
_breseq_ prints a note suggesting it instead.

The platform (Illumina, Nanopore or PacBio) is detected reliably. Whether
Nanopore reads were basecalled in fast or high-accuracy mode is not, unless
the read names include the basecaller model: base qualities differ as much
between libraries as between basecalling modes. _breseq_ still reports its
best guess, but it will not choose any setting from that guess alone.

_breseq_ was initially developed to analyze data from the Lenski
long-term evolution experiment with *E. coli* [(Link to LTEE
Website)](https://the-ltee.org).

However, _breseq_ may be generally useful to researchers who are:

1.  Tracking mutations over time in microbial evolution experiments.
2.  Checking strains for unwanted second-site mutations after genetic
    manipulations.
3.  Identifying mutations that occur during strain improvement or after
    long-term culture of engineered strains.
4.  Discovering what mutations arise in pathogens during infection or
    cause antibiotic resistance.

# Citing _breseq_

Please cite the main _breseq_ publication if you use this software
in your research:

-   Deatherage, D.E., Barrick, J.E. (2014) Identification of mutations
    in laboratory-evolved microbes from next-generation sequencing data
    using *breseq*. *Methods Mol. Biol.* **1151**: 165–188. [Link to
    Full Text](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4239701)

See the [Bibliography](bibliography.md) for a full list of papers that
describe _breseq_. We appreciate you also citing these publications if
they are relevant for your application.
