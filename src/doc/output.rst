Output
======

HTML Archive
************

|breseq| produces the results of an analysis run as a stand-alone HTML archive in the output directory. You can load these files directly in a browser, or copy the directory to a server to allow access via the web.

`output/index.html`
   This is the main results page. It consists of an upper table showing predicted mutational events and possibly several other tables showing high-quality "orphan" evidence that |breseq| was unable to assigned to mutational events. The format of each row varies depending on the kind of mutational event, as described in :ref:`mutation-display` and :ref:`evidence-display`. 

`output/marginal.html`
   This result page shows evidence items that have marginal support. Specifically: (1) RA evidence that supports a mutated base or indel more than the reference sequence, but without sufficient support to pass the cutoff threshold, and (2) JC evidence for a set number of the highest scoring junctions that do not pass tests. The format of these lines is described in :ref:`evidence-display`. 

`output/summary.html`
   This page contains additional information about the read files, reference sequences, analysis settings, and results. It includes plots showing the re-calibrated base error model, coverage distribution for each reference sequence, and  coverage across each reference sequence.

.. _mutation-display:   

Mutation Display
++++++++++++++++

Each row displays a predicted mutation in the re-sequenced sample relative to the reference. Examples showing how the format varies depending on the type of mutation are provided in the following sections. In general the columns contain: 

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

SUB: Multiple-base substitution
"""""""""""""""""""""""""""""""

INS: Sequence insertion 
"""""""""""""""""""""""""""""

DEL: Sequence deletion 
"""""""""""""""""""""""""""""

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

JC: New Junction
"""""""""""""""""""""""""""""

RA: Read alignment
"""""""""""""""""""""""""""""

MC: Missing coverage
"""""""""""""""""""""""""""""

GenomeDiff File
*****************

`output/output.gd`
   |breseq| also outputs its evidence and mutation predictions in a computer-readable file in :ref:`genomediff-format`. These files can be used as input to certain anlysis programs that will compare and tabulate the results of analyzing multiple samples.

Processed Data
**************

|breseq| outputs several files under that can be used by other tools to further analyze the final *processed* read data.

`data/reference.bam, data/reference.bam.bai`
   The BAM (Binary SAM) formatted database of read alignments to the reference and its index. Along with the *reference.fasta\** files can be used with any :program:`SAMTools` compatible program.
`data/reference.fasta, data/reference.fasta.fai`
   File of all reference sequences and the corresponding index. Along with the *reference.fasta\** files can be used with any :program:`SAMTools` compatible program.
`data/<read_file>.unmatched.fastq`
   These files contain reads from each original file that were not mapped to the reference sequences. This file can be used for de novo assembly to predict if there are novel sequences in the sample.