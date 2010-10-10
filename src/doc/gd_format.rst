.. _genomediff-format:

:program:`GenomeDiff` Format
=============================

|breseq| outputs its evidence and mutation predictions in a computer-readable :program:`GenomeDiff` text format. 

An example of a :program:`GenomeDiff` file::

   It begins with a header line. Subsequent lines beginning in # are comments or meta data.

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

Valid *mutation* types are: SNP, SUB, DEL, INS, MOB, AMP, CON, INV.

Valid *evidence* types are: RA, MC, JC, UN.

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

Line specification:

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **start** *<uint32>*

   start position in reference sequence fragment.

6. **end** *<uint32>*

   end position in reference sequence of region.
   
7. **start_range** *<uint32>*

   number of bases to offset *after* the **start position** to define the upper limit of the range where the start of a deletion could be.
   
8. **end_range** *<uint32>*

   number of bases to offset *before* the **end position** to define the lower limit of the range where the start of a deletion could be.
   
Essentially this is evidence of missing coverage between two positions in the ranges [start, start+start_range] [end-end_range, end].


NJ: New junction evidence
"""""""""""""""""""""""""

4. **side_1_seq_id** *<string>*

   id of reference sequence fragment containing side 1 of the junction.

5. **side_1_position** *<uint32>*

   position of side 1 at the junction boundary.
   
6. **side_1_strand** *<1/-1>*

   direction that side 1 continues matching the reference sequence

7. **side_2_seq_id** *<string>*

   id of reference sequence fragment containing side 2 of the junction.
   
8. **side_2_position** *<uint32>*

   position of side 2 at the junction boundary.

9. **side_2_strand** *<1/-1>*

   direction that side 2 continues matching the reference sequence.

9. **overlap** *<uint32>*
   
   Number of bases that the two sides of the new junction have in common.


UN: Unknown base evidence
"""""""""""""""""""""""""

Line specification:

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **start** *<uint32>*

   start position in reference sequence of region.

6. **end** *<uint32>*

   end position in reference sequence of region.

Mutational Event Types
++++++++++++++++++++++

SNP: Base substitution mutation
""""""""""""""""""""""""""""""""

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **position** *<uint32>*

   position in reference sequence fragment.

6. **new_seq** *<char>*

   new base at position

SUB: Multiple base substitution mutation
""""""""""""""""""""""""""""""""""""""""

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **position** *<uint32>*

   position in reference sequence fragment.

6. **size** *<uint32>*

   number of bases *after* the specified reference position to replace with **new_seq**

7. **new_seq** *<string>*

   new base at position


DEL: Deletion mutation
""""""""""""""""""""""

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **position** *<uint32>*

   position in reference sequence fragment.

6. **size** *<uint32>*

   number of bases deleted in reference


INS: Insertion mutation
"""""""""""""""""""""""

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **position** *<uint32>*

   position in reference sequence fragment.

6. **new_seq** *<string>*

   new base inserted *after* the specified rference position

MOB: Mobile element insertion mutation
""""""""""""""""""""""""""""""""""""""

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **position** *<uint32>*

   position in reference sequence fragment.

6. **repeat_name** *<string>*

   name of the mobile element. Should correspond to an annotated **repeat_region** in the reference.

7. **strand** *<1/-1>*

   strand of mobile element insertion.  

8. **duplication_size** *<uint32>*

   number of bases duplicated during insertion, beginning with the specified reference position.
   

AMP: Amplification mutation
"""""""""""""""""""""""""""

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **position** *<uint32>*

   position in reference sequence fragment.

6. **size** *<uint32>*

   number of bases duplicated starting with the specified reference position.

7. **new_copy_number** *<uint32>*

   new number of copies of specified bases. 

CON: Gene conversion mutation
"""""""""""""""""""""""""""""

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **position** *<uint32>*

   position in reference sequence fragment that was the target of gene conversion from another genomic location.

6. **size** *<uint32>*

   number of bases to replace in the reference genome beginning at the specified position.

7. **region** *<sequence:start-end>*

   Region in the reference genome to use as a replacement.

INV: Inversion mutation
"""""""""""""""""""""""

4. **seq_id** *<string>*

   id of reference sequence fragment.

5. **position** *<uint32>*

   position in reference sequence fragment.

6. **size** *<uint32>*

   number of bases in inverted region beginning at the specified reference position.