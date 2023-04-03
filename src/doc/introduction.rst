Introduction
==============

|breseq| (pronounced: \\brēz-ˈsēk\\ or *breeze-seq*) is a computational pipeline specialized for the analysis of short-read re-sequencing data (e.g. Illumina, 454, IonTorrent, etc.). It uses reference-based alignment approaches to predict mutations in a sample relative to an already sequenced genome. |breseq| is intended for microbial genomes (<10 Mb) and re-sequenced samples that are only slightly diverged from the reference sequence (<1 mutation per 1000 bp).

|breseq|'s primary advantages over other software programs are that it can:

#. Accurately predict new sequence junctions, such as those associated with mobile element insertions.
#. Integrate multiple sources of evidence for genetic changes into mutation predictions.
#. Produce annotated output describing biologically relevant mutational events.

|breseq| version 0.38+ can also be used with long-read sequencing data (e.g., Nanopore, PacBio). It analyzes long reads by automatically splitting them into shorter subsequences so that it can use its standard mapping and analysis methods. Be aware that this discards unique information in long reads that can be used to identify large-scale structural variation in a genome! We recommend that you perform de novo assembly and whole-genome comparisons with other software programs to analyze long-read data for structural variation.

|breseq| was initially developed to analyze data from the Lenski long-term evolution experiment with *E. coli* `(Link to LTEE Website) <https://the-ltee.org>`_ [Barrick2009a]_ [Barrick2009b]_\ .

However, |breseq| may be generally useful to researchers who are:

#. Tracking mutations over time in microbial evolution experiments.
#. Checking strains for unwanted second-site mutations after genetic manipulations.
#. Identifying mutations that occur during strain improvement or after long-term culture of engineered strains.
#. Discovering what mutations arise in pathogens during infection or cause antibiotic resistance.

Citing |breseq|
---------------

Please cite the main |breseq| publication if you use this software in your research:

* Deatherage, D.E., Barrick, J.E. (2014) Identification of mutations in laboratory-evolved microbes from next-generation sequencing data using *breseq*. *Methods Mol. Biol.* **1151**: 165–188. `Link to Full Text <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4239701>`_

See the :ref:`annotated_bibliography` for a full list of papers that describe |breseq|. We appreciate you also citing these publications if they are relevant for your application.
