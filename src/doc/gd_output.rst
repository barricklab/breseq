:program:`GenomeDiff` output
=============================

|breseq| outputs its evidence and mutation predictions in a computer-readable :program:`GenomeDiff` text format. 

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

Lines beginning in ``#=<name> <value>`` are interpreted as metadata ``name = value`` pairs. (Thus, the first line is assigning a metadata item named GENOME_DIFF.) ``<name>`` cannot include whitespace characters. ``<value>`` may include whitespace characters. Lines with the same ``name`` are concatenated with single-spaces added between them. 

Comment lines
++++++++++++++

Subsequent lines beginning with whitespace and ``#`` are comments.

Data lines
++++++++++++++++++++++

Data lines describe either a mutation or evidence from an analysis that can potentially support a mutational event. Each line begins with several columns containing information common to all types, then contains a fixed number of columns depending on the type, and ends with an arbitrary number of name=value pairs that store optional information.

1. **type** *<string>*

   type of the entry on this line.

2. **id** *<uint32>*

   id of this item. May be set to '+' for manually edited entries.

3. **parent-ids** *<uint32>*
   
   ids of evidence that support this mutation. May be set to '.' or left blank.

Valid *mutation* types are: SNP, INS, DEL, SUB, MOB, INV, AMP.

Valid *evidence* types are: RA, MC, JC.

Evidence Types
++++++++++++++++++++++

RA: Read alignment evidence
"""""""""""""""""""""""""""

Line specification:

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **position** *<uint32>*

   position in reference sequence fragment.

6. **insert_position** *<uint32>*

   number of bases inserted after the reference position to get to this base. An value of zero refers to the base. A value of 5 means that this evidence if for the fifth newly inserted column after the reference position.

7. **ref_base** *<char>*

   base in the reference genome.
   
8. **new_base** *<char>*

   new base supported by read alignment evidence.

MC: Missing coverage evidence
"""""""""""""""""""""""""""""

NJ: New junction evidence
"""""""""""""""""""""""""

UN: Unknown base evidence
"""""""""""""""""""""""""



Mutational Event Types
++++++++++++++++++++++

RA
""" 

*Notes:*

:type: Evidence :: Missing coverage (RA)
:spec: #. **type** (RA) 
       #. **id**
       #. **parent-id**