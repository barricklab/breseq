Performs various functions on GenomeDiff format files. Options depend on the
COMMAND supplied. Only a small subset of these commands are described
below. For a full list of _gdtools_ subcommands run it from the
command line with no options.

### gdtools ANNOTATE (or gdtools COMPARE)

Usage:

    gdtools ANNOTATE [-o annotated.html] -r reference.gbk input.1.gd [input.2.gd ... ]

Annotate a file with information about mutations (what genes they
affect, amino acid substitutions, etc.) Default output is to another ,
but an HTML table can be produced with a table of mutations in a single
file or to compare the mutations present in several files. This
subcommand can be called as ANNOTATE or COMPARE. Both have the same
effect.

`-r \<file_path>, --reference=\<file_path>`

Reference sequence files (Genbank, GFF, or FASTA). This option may be
entered multiple times. REQUIRED

`-o \<file_path>, --output=\<file_path>`

File name for the output or HTML. DEFAULT: "annotated.gd" or
"annotated.html".

`-f \<format>,--format \<format>`

Type of output file to generate. See options below (DEFAULT=HTML)

| Format | Description                                                  |
|--------|--------------------------------------------------------------|
| HTML   | Descriptive table viewable in a web browser                  |
| GD     | GenomeDiff with added annotation of mutations                |
| TSV    | Tab-separated values file suitable for input into R or Excel |
| PHYLIP | Alignment file suitable for input into PHYLIP                |
| JSON   | JavaScript object notation file suitable for parsing         |

`input1.gd input2.gd ...  `

Input file(s). This option may be entered multiple times to compare
across files. REQUIRED

!!! warning
    Some advanced attributes for mutations, such as **within** and
    **before**, are ignored when generating compare tables.

### gdtools APPLY

Usage:

    gdtools APPLY [ -o output.gff3 -f GFF3 ] -r reference.gbk input.gd

Apply the mutations described in the input to the reference sequence(s).

A mutation can cross the origin of a circular reference sequence: its position is
near the end of the sequence and its size runs past the last base. The new
sequence is divided at the origin, keeping as many of its bases before the
origin as there were, so a change that does not alter the length (`SUB`, `INV`)
leaves every other coordinate where it was. Features inside an `INV` that
crosses the origin are not carried over to the new sequence (a warning says so).
On a linear sequence, a mutation that runs past the end is cut off there, with a
warning. If the sequence is actually circular, mark it `CIRCULAR` in the
reference file (see [Reference Sequence File Formats](reference-sequence-file-formats.md))
and the warning goes away.

`-r <file_path>, --reference=<file_path>`

Reference sequence files (Genbank, GFF, or FASTA). This option may be
entered multiple times. REQUIRED

`input.gd`

Input file. REQUIRED

`-o <file_path>, --output=<file_path>`

Output file containing the mutated reference genome. DEFAULT: "output.\*"

`-f <output_format>, --format=<output_format>`

Output format. Possible values: `GENBANK`, `FASTA`, or `GFF3`.

`-m <file_path>, --coordinate-map=<file_path>`

Also write a tab-delimited map from the coordinates of the output sequences back to those of the
input reference sequences. Each line is one block of the output that came from one stretch of the
reference: `seq_id`, `applied_start`, `applied_end`, `original_start`, `original_end` and `strand`
(`-` inside an applied inversion, where the original coordinates run backwards). A block of
newly inserted sequence has `.` in the last three columns. Together the blocks cover every base
of the output. This is the map that `breseq --apply-check` uses to report `original_*`
coordinates.

### gdtools CONVERT

Usage:

    gdtools CONVERT -f GVF [ -o output.gvf -a ] -r reference.gbk input.gd

Convert a GenomeDiff file to another format (`GD`, `VCF`, `GVF`, or `JSON`), or a
VCF file to GenomeDiff. `gdtools GD2VCF`, `VCF2GD` and `GD2GVF` are older names
for the same conversions.

`-r <file_path>, --reference=<file_path>`

Reference sequence files (Genbank, GFF, or FASTA). This option may be
entered multiple times. REQUIRED for VCF and GVF output

`-a, --annotate`

Annotate the mutations first. In GVF output this adds `Variant_effect` and the
codon and amino acid attributes to SNPs.

`--gvf-max-sequence-length=<bases>`

In GVF output, a `Reference_seq` or `Variant_seq` longer than this is written as
`~<length>`. Zero means always write the full sequence. DEFAULT: 50

#### How mutations are written as VCF

Each record replaces `REF` with `ALT`, using the same sequences as in the GVF
table below, so it also describes exactly the change that `gdtools APPLY` makes.
Alleles are always written out in full; symbolic alleles (`<DEL>`, `<INS:ME>`)
are not used.

- VCF does not allow an empty allele, so an insertion or a deletion includes the
  base before it (the base after it, at the start of a sequence).
- `INS` entries at one site that differ in `insert_position` are written as one
  record when they are at the same frequency. Separate records with the same
  `POS` would be read as alternatives to one another.
- A mutation that crosses the origin of a circular sequence is written as two
  records, one on each side of it, with the new sequence divided between them
  as `gdtools APPLY` divides it: as many of its bases as there were before the
  origin stay there.
- A mutation `within` another one is omitted with a warning.

#### How mutations are written as GVF

Output follows [Genome Variation Format 1.10](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gvf.md).
Each line says that `Variant_seq` replaces the reference bases from `start` to
`end`, which are `Reference_seq`. It describes exactly the change that
`gdtools APPLY` makes for the same entry.

| GenomeDiff | GVF type | start..end | Reference_seq | Variant_seq |
|---|---|---|---|---|
| SNP | `SNV` | position | the base | new base |
| SUB, same length | `MNV` | replaced bases | those bases | new bases |
| SUB, otherwise | `indel` | replaced bases | those bases | new bases |
| DEL | `deletion` | deleted bases | those bases | `-` |
| INS | `insertion` | position (the site is 3' of it) | `-` | inserted bases |
| AMP, new copy number 2 | `tandem_duplication` | amplified bases | one copy | two copies |
| AMP, otherwise | `copy_number_gain` | amplified bases | one copy | all of the copies |
| MOB | `mobile_element_insertion` | the base before the target site duplication | `-` | new copy of the target site, then the element |
| MOB that deletes target site bases | `indel` | deleted bases | those bases | the element |
| INV | `inversion` | inverted bases | those bases | their reverse complement |
| CON, INT | `substitution` if the same length, otherwise `indel` | replaced bases | those bases | the donor sequence |

- Sequences are always on the reference strand, so the strand column is always
  `+`. The orientation of a mobile element is in `repeat_strand`.
- `-` means no sequence. `~<length>` stands in for a sequence that is longer
  than `--gvf-max-sequence-length`.
- A mutation that is not at 100% frequency lists both alleles, the new one
  first: `Variant_seq=C,A;Variant_freq=0.2500,0.7500;Variant_reads=11:33`.
- The score is the average score of the evidence supporting the mutation. It is
  not Phred scaled.
- `gd_type` and `gd_id` give the type and ID of the GenomeDiff entry. `AMP`
  lines have `copy_number`, `MOB` lines `repeat_name`, `repeat_strand` and
  `duplication_size`, and `CON`/`INT` lines a `Breakpoint_detail` for the donor
  region.
- A feature that crosses the origin of a circular sequence has an `end` past the
  length of the sequence, as in GFF3.
- Several lines at one site (e.g., `INS` entries that differ in
  `insert_position`) are in 5' to 3' order.
- A mutation `within` another one is omitted with a warning. Its position counts
  bases of the new sequence, so it has no reference coordinate.

