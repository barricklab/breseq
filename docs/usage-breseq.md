## breseq usage

    breseq -r reference1.gbk [-r reference2.gbk ...] reads1.fastq [reads2.fastq, reads3.fastq ...]

Run the `breseq` mutation prediction pipeline.

Required options:

`-r <file_path>, --reference <file_path>`

Input reference genome sequence files in GenBank, GFF3, or FASTA format.
If there are multiple reference sequences stored in separate files
(e.g., a bacterial genome and a plasmid), this option can be supplied
multiple times.

`reads1.fastq [reads2.fastq, reads3.fastq ...]`

The remaining arguments at the command line are the FASTQ input files of
reads. FASTQ files with base quality scores that are not in [SANGER
format](https://en.wikipedia.org/wiki/FASTQ_format) will be converted.
In addition, reads with \>50% N bases will be removed from the converted
FASTQ file by default. _breseq_ re-calibrates the error rates
for each FASTQ file separately, so data sets that were generated
independently should be stored in different input files.

Commonly used options:

`-h, --help`

Produce help message showing advanced options.

`-n <string>, --name <string>`

Human-readable name of the analysis run for output (DEFAULT=\<none>).

`-j <int>, --num-processors <int>`

Number of processors to use in multithreaded steps (DEFAULT=1).

`-p, --polymorphism-prediction`

Predict polymorphic mutations. Add this option when you are analyzing mixed population (metagenomic) samples.

`--no-copy-number-prediction`, `--no-discordant-pair-prediction`, `--no-missing-pair-prediction`, `--no-pair-distance-prediction`

Copy number (CN), discordant pair (DP), missing pair (MP) and pair distance (PD) evidence are
predicted by default. Each of these flags turns one of them off. The older opt-in flags
(`--predict-copy-number`, `--predict-discordant-pairs`, `--predict-missing-pairs`,
`--predict-pair-distance`) are deprecated: they are still accepted so that existing command lines
keep working, but they no longer do anything.

DP, MP and PD need paired-end reads, so they are skipped for a single-end run or under
`--no-paired-mapping`. CN prediction needs the separate
[CNery](https://github.com/barricklab/CNery) program on your `PATH`; if it is missing, _breseq_
warns and skips CN rather than failing.

Soft clipping (SC) evidence is the exception: it remains opt-in via `--predict-soft-clipping`,
because that option also lowers `--require-match-fraction` from 0.9 to 0.5 and so changes which
read alignments are accepted throughout the analysis.

`--copy-number-resolution <float>`

Spacing of the copy-number grid, in copies (DEFAULT=0.1). **Polymorphism mode only.** In consensus
mode copy number is always called on the integers. Under `-p` it is called on a continuous grid
instead, so a region carried by part of the population is reported at the depth it was measured at —
a `CN` entry may read `copy_number=1.4`. This option sets how finely that grid is spaced; larger
values call fewer, coarser levels. CNery rounds the value to a spacing that divides 1.0 exactly, so
that single copy is always on the grid, and requires it to be above 0.02, the coverage a deleted
region is modeled as retaining.

Note that an `AMP` mutation predicted from a fractional `CN` still carries a whole number of copies
in its `new_copy_number` field — it is a count of copies in a genome, not a measured depth — so a
`CN` reading 2.4 beside an `AMP` reading 3 is expected rather than a discrepancy.

Because a level is only called when it beats the cost of a state change, the grid spacing is also
what keeps ordinary coverage noise from being reported as a copy-number change. It cannot suppress
everything: on a short or heavily fragmented reference, where there is too little data to pin the
single-copy level down, residual mapping and GC bias can still be called as a level slightly off 1.
Treat a `CN` entry near single copy on a small contig with suspicion.

`--no-read-alignment-prediction`, `--no-missing-coverage-prediction`, `--no-junction-prediction`, `--no-homologous-deletion-prediction`

Turn off one kind of evidence prediction. `--no-read-alignment-prediction` also turns off missing
coverage, because both come from the same pileup pass (as do CN, DP, MP and PD, which read what
that pass writes). The older `--skip-RA-MC-prediction`, `--skip-MC-prediction`,
`--skip-JC-prediction` and `--skip-homologous-DEL-prediction` spellings are deprecated: they still
work, but every opt-out is now spelled `--no-X-prediction`.

`--dry-run`

Validate every option, check that the required external programs (`bowtie2`, `gnuplot`,
`samtools`) are installed, and check that every input file exists and every output path can
be written &mdash; then exit without running the pipeline and without creating any files.
Exits with status 0 if everything checks out and non-zero otherwise, so it can gate a real
run: `breseq --dry-run -r reference.gbk reads.fastq && breseq -r reference.gbk reads.fastq`.

Note that the file and folder checks themselves happen on *every* run, not just this one: a
mistyped path fails immediately rather than part way through the analysis. What `--dry-run`
adds is stopping afterwards, and reporting each path it checked.

!!! tip
    For a complete list of options (including many advanced options), please show the full command line help by running `breseq -h` or `breseq --help`.

## Utility subcommands

_breseq_ provides some additional subcommands for further analysis. The subcommands should be used _after_ the main pipeline has been run. It is easiest to run them from within the main output directory of a breseq run, which will include the required `data/reference.fasta` and `data/reference.bam` files, so that you don't have to specify these options on the command line.

### breseq BAM2ALN

Usage:

    breseq BAM2ALN [-b reference.bam -f reference.fasta -o alignment.html -n 200] region1 [region2 region3 ...]

Create an HTML file displaying reads aligned to the specified region or regions.

Commonly used options:

`-h, --help`

Produce help message showing advanced options.

`-b <file_path>, --bam <file_path>`

BAM database file of read alignments (DEFAULT=`data/reference.bam`).

`-f <file_path>, --fasta <file_path>`

FASTA file of reference sequences (DEFAULT=`data/reference.fasta`).

`-o <path>, --output <path>`

Output path. If there is just one region, the name of the output file
(DEFAULT=region1.*). If there are multiple regions, this argument must
be a directory path, and all output files will be output here with names
region1.*, region2.\*, ... (DEFAULT=`.`).

`-r <region> , --region <region>`

Regions to create alignments for. Must be provided as sequence regions
in the format **ACCESSION:START-END**, where **ACCESSION** is a valid
identifier for one of the sequences in the FASTA file, and **START** and
**END** are 1-indexed coordinates of the beginning and end positions.
Any read overlapping these positions will be shown. A separate output
file is created for each region. Regions may be provided at the end of
the command line as unnamed arguments.

`-n <int>, --max-reads <int>`

Maximum number of reads that will be aligned to a region. If there are
more than this many reads, then the reads displayed are randomly chosen
and a warning is added to the output. (DEFAULT=200).

### breseq BAM2COV

Usage:

    breseq BAM2COV [-b reference.bam -f reference.fasta --format PNG -o output.png] region1 [region2 region3 ...]

Create a coverage plot or table for the specified region or regions.

Commonly used options:

`-h, --help`

Produce help message showing advanced options.

`-b <file_path>, --bam <file_path>`

BAM database file of read alignments (DEFAULT=`data/reference.bam`).

`-f <file_path>, --fasta <file_path>`

FASTA file of reference sequences (DEFAULT=`data/reference.fasta`).

`-o <path>, --output <path>`

Output path. If there is just one region, this is the name of the output file and is
used exactly as given (include the desired extension, e.g. `region1.png` or
`region1.tsv`). If there are multiple regions, this argument must be a directory path,
and all output files will be written there with automatically generated names
`region1.<ext>`, `region2.<ext>`, ... where `<ext>` matches the `--format` (DEFAULT=.).

`-r <region>, --region <region>`

Regions to create alignments for. Must be provided as sequence regions
in the format **ACCESSION:START-END**, where **ACCESSION** is a valid
identifier for one of the sequences in the FASTA file, and **START** and
**END** are 1-indexed coordinates of the beginning and end positions.
Any read overlapping these positions will be shown. A separate output
file is created for each region. Regions may be provided at the end of
the command line as unnamed arguments.

`--format <PNG/PDF/SVG/TSV/CSV>`

Format of output. `PNG`, `PDF`, or `SVG` create a coverage **plot**. `TSV`
(tab-separated) or `CSV` (comma-separated) create a coverage **table** instead of a
plot. `TSV` reproduces the output of the deprecated `--table` option. (DEFAULT=PNG).

`-t, --table`

**DEPRECATED** — use `--format TSV` instead. If provided, the output format is set to
`TSV` (a tab-delimited coverage table) and a deprecation warning is printed.

`-1, --total-only`

Only plot/tabulate the total coverage at a position. That is, do not not
output the coverage on each genomic strand.

`--per-read-group`

Repeat every coverage column once per read group (`@RG`) in the BAM file, prefixed
`RG-#_`, where `#` is the read group's index in the BAM header. The normal columns are
still output first and are unchanged, so this only *adds* columns; at each position the
per-read-group values sum to the corresponding total. The same repeat is applied to the
`#`-commented region averages at the end of the table.

breseq writes one read group per read file set, so a pair of paired-end read files shares
a single read group. A BAM file with no read groups produces a single `RG-0` set holding
all of the coverage.

Requires a table output format (`--format TSV` or `--format CSV`); combining it with a
plot format is an error.

`--resolution <int>`

Number of positions to output coverage information for in interval
(0=ALL) (DEFAULT=600).
