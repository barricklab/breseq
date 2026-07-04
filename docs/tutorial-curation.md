_breseq_ can predict many mutations correctly, but it doesn't always find all mutations or get all of the predictions entirely correct! Some human intervention is needed to resolve certain complex cases. 

In this tutorial, we'll cover how you can curate the GenomeDiff files output by _breseq_ so that they are suitable for downstream analyses such as measuring the rates and spectra of mutations in evolved genomes and analyzing whether mutations affecting certain genes are more common in certain treatments.

This tutorial also introduces [_brefito_](https://github.com/barricklab/brefito), a set of Snakemake pipelines that can be used to automate these and other analysis steps. Additionally, it covers several utility commands — for analyzing read alignments in more depth, comparing mutations across samples, constructing mutant genomes, and creating phylogenetic trees — that are useful for checking whether your curated mutation files are 100% correct and consistent.

We'll examine resequencing data from clonal isolates that have evolved for 75,000 generations in the [Long-Term Evolution Experiment with _E. coli_ (LTEE)](https://the-ltee.org).

1. [Running _breseq_](tutorial-curation-running-breseq.md)
1. [Editing GenomeDiff files](tutorial-curation-editing-genomediff-files.md)
1. [Validating your predictions](tutorial-curation-validation.md)
1. [Exploring aligned reads](tutorial-curation-exploring-aligned-reads.md)
1. [Automating the curation cycle using _brefito_](tutorial-curation-automating.md)
1. [Common curation cases](tutorial-curation-common-cases.md)
1. [Curating _E. coli_ LTEE genomes](tutorial-curation-ecoli-ltee.md)