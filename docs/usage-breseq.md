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

`--no-junction-prediction`

Do not predict new sequence junctions.

`-p, --polymorphism-prediction`

Predict polymorphic mutations. Add this option when you are analyzing mixed population (metagenomic) samples.

`-x, --nanopore`

Set read alignment and mutation calling options for Nanopore sequencing
data. Important: no indel mutations will be called in homopolymer
repeats of 4 or more bases with this option.

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

Output path. If there is just one region, the name of the output file
(DEFAULT=region1.*). If there are multiple regions, this argument must
be a directory path, and all output files will be output here with names
`region1.*`, `region2.*`, ... (DEFAULT=.).

`-r <region>, --region <region>`

Regions to create alignments for. Must be provided as sequence regions
in the format **ACCESSION:START-END**, where **ACCESSION** is a valid
identifier for one of the sequences in the FASTA file, and **START** and
**END** are 1-indexed coordinates of the beginning and end positions.
Any read overlapping these positions will be shown. A separate output
file is created for each region. Regions may be provided at the end of
the command line as unnamed arguments.

`--format <PNG/PDF>`

Format of output plot: PNG or PDF. (DEFAULT=PNG).

`-t, --table`

Create tab-delimited file of coverage instead of a plot.

`-1, --total-only`

Only plot/tabulate the total coverage at a position. That is, do not not
output the coverage on each genomic strand.

`--resolution <int>`

Number of positions to output coverage information for in interval
(0=ALL) (DEFAULT=600).
