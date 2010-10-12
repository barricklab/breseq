Introduction
==============

|breseq| (pronounced: *breeze-seq*) is a computational pipeline for the analysis of short-read re-sequencing data (e.g. 454, Illumina, SOLiD, etc.). It uses reference-based alignment approaches to predict mutations in a sample relative to an already sequenced genome. |breseq| is intended for microbial genomes (<10 Mb) and re-sequenced samples that are only slightly diverged from the reference sequence (<1 mutation per 1000 bp). 

|breseq|'s primary advantages over other existing software programs are that it can:

#. Predict new sequence junctions, such as those associated with mobile element insertions, from single-end read data.
#. Reliably identify short indel mutations by appropriately masking the ends of read alignments.
#. Produce annotated output describing biologically relevant mutational events.

|breseq| has been used to analyze data from the Lenski long-term evolution experiment with *E. coli* [Barrick2009a]_\ [Barrick2009b]_\ .

|breseq| should be generally useful to researchers who are:

#. Following mutations over time in microbial evolution experiments.
#. Checking strains for second-site mutations after genetic manipulations.
#. Identifying mutations that occur during strain improvement or after long-term culture of engineered strains.
#. Discovering what mutations arise in pathogens during infection or cause antibiotic resistance.