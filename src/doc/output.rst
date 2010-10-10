Output
======

HTML Archive
************

|breseq| produces the results of an analysis run as a stand-alone HTML archive in the output directory. You can load these files directly in a browser, or copy the directory to a server to allow access via the web.

Important files include:

`output/index.html`
   The main results page. It consists of an upper table showing predicted mutational events and possibly several other tables showing high-quality "orphan" evidence that |breseq| was unable to assigned to mutational events. The format of each row varies depending on the kind of mutational event, as described in :ref:`mutation-display` and :ref:`evidence-display`. 

`output/marginal.html`
   Result page showing evidence for mutations with marginal support. Specifically: (1) RA evidence that supports a mutated base or indel more than the reference sequence, but without sufficient support to pass the cutoff threshold, and (2) JC evidence for a set number of the highest scoring junctions that do not pass tests. The format of these lines is described in :ref:`evidence-display`. 

`output/summary.html`
   Additional information about the read files, reference sequences, analysis settings, and results. Links to plots showing the re-calibrated base error model, coverage distribution for each reference sequence, and coverage across each reference sequence.
   
`output/output.gd`
   Text file of evidence and mutation predictions in computer-readable :ref:`genomediff-format`. This file can be used as input to certain meta-analysis programs to analyze mutations from many samples. 
   
`output/log.txt`
   The original command line used to invoke |breseq|\ .     

.. _mutation-display:   

Mutation Display
++++++++++++++++

Each row displays a predicted mutation in the re-sequenced sample relative to the reference. Examples showing how the format varies depending on the type of mutation are provided in the following sections. 

Column descriptions: 

`evidence`
	Links to the types of evidence that support this particular prediction. See :ref:`evidence-display`.
`seq id`	
	The identifier for the sequence with the mutation.
`position`
	Position in the reference sequence of the mutation. Generally this is where the mutation begins when the mutation affects a range of positions.  
`mutation`
	Description of the mutation. Typically describes how nucleotides are added, substituted, or deleted. May also refer to a mobile element in the genome and how it is inserted at the specified position.
`annotation`
	Description of the mutation. For base substiturions inside genes, indicates the amino acid and codon changes. For other mutations inside genes, gives the coding nucleotides affected. For mutations in intergenic regions gives two relative positions (``+150/-119``) where the number is how many nucleotides from the mutation to the nearest neighboring genes, respectively, and +/- indicate whether the mutation is oriented upstream or downstream of each gene.
`gene`
	Genes affected by the mutation. May be a single gene (``gene``), in an integenic region between genes (``gene1/gene2``), covering a range that completely encompasses several genes (``gene1-geneN``). Brackets around a gene ([gene]) mean that the mutations range ends within that gene.
`description`
	Descriptions of the genes affected by the mutation. Generally these correspond to the genes in the gene column, but if many genes are affected, this field is abbreviated to be a list of all the affected genes with brackets signifying that the mutation begins or ends within the specified gene.
	
All annotation information is taken from the input Genbank files. How informative descriptions of these changes are about how mutations affect genes is dependent on the quality of the reference sequences.

SNP: Single-base substitution
"""""""""""""""""""""""""""""

.. figure:: images/mutation_header.png

.. figure:: images/snp_1.png

SUB: Multiple-base substitution
"""""""""""""""""""""""""""""""

.. figure:: images/mutation_header.png

.. figure:: images/sub_1.png

INS: Sequence insertion 
"""""""""""""""""""""""""""""

.. figure:: images/mutation_header.png

.. figure:: images/ins_1.png

DEL: Sequence deletion 
"""""""""""""""""""""""""""""

For deletion rows, the *position* column gives the first missing reference base and the *mutation* column gives the size of the deletion. Thus, the deleted reference region extends from *position* to *position* + *size* -1.

.. figure:: images/mutation_header.png

.. figure:: images/del_1.png

A single-base deletion at position 139 in an intergenic region at the end of the reference sequence. The deleted nucleotide is located 52 bp downstream of the end of the first gene *nu1* in the genome by. This mutation is supported by :ref:`read-alignment-display` evidence.

.. figure:: images/mutation_header.png

.. figure:: images/del_2.png

A 5,996 bp deletion starting at position 2,338. This deletion begins within the *orf-314* gene and ends past the *ea59* gene. This mutation is supported by :ref:`new-junction-display` and :ref:`missing-coverage-display` evidence.

MOB: Mobile element insertion
"""""""""""""""""""""""""""""

AMP: Sequence amplification
"""""""""""""""""""""""""""""

INV: Chromosomal inversion
"""""""""""""""""""""""""""""

.. _evidence-display:   

Evidence Display
++++++++++++++++

Note that clicking on any evidence link for a mutation prediction will bring up pages with tables showing all items of evidence that |breseq| used to predict the mutational event.

.. _new-junction-display:   

New Junction (JC)
"""""""""""""""""""""""""""""

Each JC row consists of two parts sub-rows, one describing one side of the junction in the reference sequence.

Column descriptions: 

`* link`
    Links to a results page showing the sequence of the new junction as the reference and all reads aligned to the junction.
`? links`
   Links to a results pages for each side of the juncton, that show the original reference sequence at that site and any reads that aligned better to this sequence than to the new junction.  Note that in some cases (such as tandem duplications), it is possible for the new and old junction sequences to still exist in the sample. You can check for this by examining these read alignments. Sequences where the read name has a -M1 or -M2 appended are reads that mapped better to the new junction.
`seq id`	
	The identifiers for the sequences involved in the new junction.
`position`
	Positions in the reference sequence of the two sides of the new junction. Each position has an = before or after it that represents how the junction was constructed. If
`overlap`
    If positive, the number of bp in the junction that could map to either side in the reference sequence. Generally, positive overlap has been resolved to zero by assigning these base pairs to one side of the junction. If negative, the number of bp that are unique to reads mapping across the junction and represent insertions relative to the reference sequence.
`reads`
    The total number of reads that map to this junction.
`score`
    The pos-hash score for the junction in **<bold angle brackets>** and the minimum-overlap score on the next line.
`annotation, gene, product`
	Description of the mutation effects for each side of the junction. The format of these columns is the same as in :ref:`mutation-display`.

Example: 

.. _read-alignment-display:

Read alignment (RA)
"""""""""""""""""""""""""""""

Column descriptions: 

`* link`
    Links to a results page showing the alignment of reads to
`seq id`	
	The identifier for the reference sequence where the evidence is located.    
`position`
   Position in the reference sequence of the single base substitution, insertion, or deletion. Consists of two parts. The first is the reference position, the second is an "insert count" that indicates this is in a column of the alignment that does not exist in the reference sequence (i.e., it is an insertion relative to the reference) and is this many columns past the specified reference position.
`change`	
	The base change.
`freq`	
	Frequency of this base change in the sample. |breseq| currently inly predicts mutations of 0% or 100% frequency.	
`score`
	The log10 ratio of the posterior probability that this position in the sample is the called base to the probability that it is any other base,  minus log10 of the total number of positions in all reference sequences. The higher the score, the more evidence for the mutation.
`cov`
    The number of reads overlapping the mutation. Note that portions of reads that are not aligned (lowercase bases with a white background), ends of reads that have been trimmed because alignments may be ambiguous (lowercase bases with a colored background) and read positions with very low base quality scores that typically indicate sequencing errors (highlighted in yellow) are not counted in this coverage number.
`annotation, gene, product`
	Description of the mutation effects for each side of the junction. The format of these columns is the same as in :ref:`mutation-display`.

Example: 

.. _missing-coverage-display:

Missing coverage (MC)
"""""""""""""""""""""""""""""

Column descriptions: 

`* links`
    Links to results pages showing the alignment of reads to the left and right margins of the region with missing coverage.
`÷ link`	
	Link to the results page showing a plot of the read coverage in the region of the msising coverage.   
`seq id`	
	The identifier for the reference sequence where the evidence is located.   
`start, end, size`
	The reference positions of the missing coverage. May indicate a range of positions when one end of the missing coverage is in a repeat region.
`← cov`
	Unique read coverage depth on the left margin of the region of missing coverage. Coverage at the first position outside the alignment is shown followed by coverage at the first position inside the region of missing coverage in brackets.
`→ cov`
	Unique read coverage depth on the right margin of the region of missing coverage. Coverage at the first position outside the alignment is shown followed by coverage at the first position outside the region.
`gene, description`
	Description of the mutation effects for each side of the junction. The format of these columns is the same as in :ref:`mutation-display`.

Example: 

Processed Data
**************

|breseq| outputs several files that can be used by other tools to further analyze the final *processed* read data.

`data/reference.bam, data/reference.bam.bai`
   The BAM (Binary SAM) formatted database of read alignments to the reference and its index. Along with the *reference.fasta\** files can be used with any :program:`SAMTools` compatible program.
`data/reference.fasta, data/reference.fasta.fai`
   File of all reference sequences and the corresponding index. Along with the *reference.fasta\** files can be used with any :program:`SAMTools` compatible program.
`data/<read_file>.unmatched.fastq`
   These files contain reads from each original file that were not mapped to the reference sequences. This file can be used for de novo assembly to predict if there are novel sequences in the sample.