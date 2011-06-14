.. _genomediff-format:

|GD| Utility Program
====================

Performs various functions on genomediff formatted files. Options depend on the COMMAND supplied.

Command: VALIDATE
-----------------

Usage:

   :program:`genomediff` VALIDATE [-o output.gb] input.gd

.. program:: genomdiff VALIDATE

Checks the syntax of a |GD| file and outputs a clean version. Entries are re-numbered consecutively and unused evidence and validation entries are removed form the output |GD|.

.. option:: -o <file_path>, --output=<file_path>

   Output file path. DEFAULT: input.validated.gd

<input.gd>
   Input |GD| file. Only one is allowed. REQUIRED

Command: HTML
---------------

Usage:

   :program:`genomediff` HTML -r reference.gbk [-o output.html] input.gd

.. program:: genomdiff HTML

Create an HTML table of mutations in a |GD| file.

.. option:: -r <file_path>, --reference=<file_path>

   GenBank files for reference sequences. This option may be entered multiple times. REQUIRED

.. option:: -o <file_path>, --output=<file_path>

   Output HTML file containing the mutation table. DEFAULT: "input.html".

<input.gd>
   Input |GD| file. Only one is allowed. REQUIRED


Command: COMPARE
----------------

Usage:

   :program:`genomediff` COMPARE -r reference.gbk [-o output.html] input1.gd [input2.gd ...]

.. program:: genomdiff COMPARE

Create an HTML table comparing mutations from different samples.

.. option:: -r <file_path>, --reference=<file_path>

   GenBank files for reference sequences. This option may be entered multiple times. REQUIRED

.. option:: -o <file_path>, --output=<file_path> 

   Output HTML file containing the comparison table. DEFAULT: "compare.html".

<input1.gd [input2.gd ...]>
   Input |GD| files, one for each sample. REQUIRED
   
Command: SUBTRACT
-----------------

Usage:

   :program:`genomediff` SUBTRACT -o output.gd -1 input1.gd [input2.gd ...]

.. program:: genomdiff SUBTRACT

Create a |GD| containing only mutations that are in |GD| #1 but not in |GD| #2.

.. option:: -o <file_path>, --output=<file_path> 

   Output |GD| file. DEFAULT: "output.gd".

.. option:: -1 <file_path>, -input1=<file_path>

   Input |GD| file #1. REQUIRED

.. option:: -2 <file_path>, -input2=<file_path>

   Input |GD| file #2. REQUIRED
