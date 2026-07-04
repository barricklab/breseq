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

`-r <file_path>, --reference=<file_path>`

Reference sequence files (Genbank, GFF, or FASTA). This option may be
entered multiple times. REQUIRED

`input.gd`

Input file. REQUIRED

`-o <file_path>, --output=<file_path>`

Output file containing the mutated reference genome. DEFAULT: "output.\*"

`-f <output_format>, --format=<output_format>`

Output format. Possible values: `GENBANK`, `FASTA`, or `GFF3`.

