Introduction
==============

|breseq| is a computational pipeline for the analysis of short-read re-sequencing data. It maps reads to a reference sequence and predict changes between the sequenced sample and this reference sequence. |breseq| is intended for re-sequencing microbial-sized" genomes (roughly <10 Mb) that are only slightly diverged from the reference sequence (roughly <1 mutation per 1000 bp). It will perform best when the reference sequence is richly annotated. |breseq|'s primary advantages over other existing software programs are that it can predict new new sequence junctions from single-end read data and that it produces annotated output of biologically relevant mutational events.

|breseq| can be used for a variety of applications:

#. Examining mutations over time in microbial evolution experiments.
#. Checking strains for second-site mutations after genetic manipulations.
#. Discovering mutations that occur during strain improvement.
#. Following the molecular evolution of pathogens during infection.


