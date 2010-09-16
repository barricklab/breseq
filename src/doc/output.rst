HTML output
=============

index.html
***************

evidence pages
***************

marginal.html
***************

summary.html
***************

:program:`GenomeDiff` output
=============================

:program:`breseq` outputs its evidence and mutation predictions in a computer-readable :program:`GenomeDiff` text format. 

An example of a :program:`genomediff` file::

   It begins with a header line. Subsequent lines beginning in # are comments or meta data.

format specification
********************

Version line
+++++++++++++++

The first line of the file must define the version::
   
   #=GENOME_DIFF 1.0

Metadata lines
+++++++++++++++

Lines beginning in ``#=<name> <value>`` are interpreted as metadata ``name = value`` pairs. (Thus, the first line is assigning the GENOME_DIFF metadata.) ``<name>`` cannot include whitespace characters. ``<value>`` may include whitespace characters. Lines with the same ``name`` are concatenated. 

Comment lines
++++++++++++++

Subsequent lines beginning in ``#`` are comments. They are skipped by the parser.

Data lines
++++++++++++++++++++++

Data lines describe either a mutation or evidence from an analysis that can potentially support a mutational event. Each line begins with several columns containing information common to all types, then a fixed number, and finally it ends with an arbitrary number of name=value pairs that store optional information.

#. **type** *<string>*

   type of the entry on this line.

#. **id** *<uint32>*

   id of this item. May be set to '+' for manually edited entries.

#. **parent-ids** *<uint32>*
   
   ids of evidence that support this mutation. May be set to '.' or blank.


Valid *mutation* types are: SNP, INS, DEL, SUB, MOB.

Valid *evidence* types are: RA, MC, JC.

Evidence Types
++++++++++++++++++++++

RA
""

RA lines specify Evidence > Read alignment

Line specification:

#. **type** *<string>*

   definition

#. **id** *<uint32>*

   definition 2

#. **parent-id** *<uint32>*

Evidence Types
++++++++++++++++++++++

RA
""" 

*Notes:*

:type: Evidence :: Missing coverage (RA)
:spec: #. **type** (RA) 
       #. **id**
       #. **parent-id**