.. _genomediff-format:

|gdtools| Utility Program
===========================

Performs various functions on |Genome Diff| formatted files. Options depend on the COMMAND supplied. Only a small subset of these commands are described below.
For a full list of |gdtools| subcommands run it from the command line with no options::

  gdtools
  
Command: ANNOTATE/COMPARE
-------------------

Usage::

  gdtools ANNOTATE/COMPARE [-o annotated.html] -r reference.gbk input.1.gd [input.2.gd ... ]

.. program:: gdtools ANNOTATE

Annotate a |Genome Diff| file with information about mutations (what genes they affect, amino acid substitutions, etc.)
Default output is to another |Genome Diff|, but an HTML table can be produced with a table of mutations in a single |Genome Diff| 
file or to compare the mutations present in several |Genome Diff| files.

.. option:: -r <file_path>, --reference=<file_path>

   Reference sequence files (Genbank, GFF, or FASTA). This option may be entered multiple times. REQUIRED

.. option:: -o <file_path>, --output=<file_path>

   File name for the output|Genome Diff| or HTML. DEFAULT: "annotated.gd" or "annotated.html".

.. option:: --html

   Output an HTML table instead of a |Genome Diff| file.

input1.gd input2.gd ...
   Input |Genome Diff| file(s). This option may be entered multiple times to compare across files. REQUIRED

.. WARNING::
   Advanced |Genome Diff| attributes such as **within** and **before** for mutations are ignored when generating compare tables.

Command: APPLY
----------------

Usage::

  gdtools APPLY [ -o output.gff3 -f GFF3 ] -r reference.gbk input.gd

.. program:: gdtools APPLY

Apply the mutations described in the input |Genome Diff| to the reference sequence(s).

.. option:: -r <file_path>, --reference=<file_path>

   Reference sequence files (Genbank, GFF, or FASTA). This option may be entered multiple times. REQUIRED

.. option:: -o <file_path>, --output=<file_path>

   Output file containing the muta. DEFAULT: "output.*"

.. option:: -f <output_format>, --format=<output_format>

   Output format. Possible values are "fasta" or "gff3".

input.gd
   Input |Genome Diff| file. REQUIRED
   
Command: SUBTRACT
-----------------

Usage::

  gdtools SUBTRACT [-o output.gd] input.gd subtract1.gd [subtract2.gd ...]

.. program:: gdtools SUBTRACT

Creates a new |Genome Diff| file of mutations from the input file that are present after removing mutations present in any of the subtracted |Genome Diff| files.

.. option:: -o <file_path>, --output=<file_path> 

   Output |Genome Diff| file. DEFAULT: "output.gd".

input.gd

   Input |Genome Diff| file.

subtract.gd [subtract2.gd ...]

   |Genome Diff| files to subtract from input file.

Command: INTERSECT
------------------

Usage::

  gdtools INTERSECT [-o output.gd] input1.gd input2.gd ...

.. program:: gdtools INTERSECT

Creates a new |Genome Diff| file with mutations that are present in ALL input |Genome Diff| files.

.. option:: -o <file_path>, --output=<file_path> 

   Output |Genome Diff| file. DEFAULT: "output.gd".

input1.gd input2.gd ...

   Input |Genome Diff| files.
