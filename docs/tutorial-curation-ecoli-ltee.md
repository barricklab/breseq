[Back to the Main Curation Tutorial Page](tutorial-curation.md)

This tutorial illustrates some advanced curation topics through analyzing
genomes from the [_E. coli_ Long-Term Evolution Experiment (LTEE)](https://the-ltee.org/). 

The long duration of this experiment and the large number of genomes that have been sequenced 
from it lead to some rare (and more challenging) cases occasionally cropping up due to how many
mutations have accumulated. Because of this, we have created some extra utilities to help us curate these genomes.

## Clone and set up the LTEE-Ecoli repository

To get started with the workflow, you need to clone the LTEE-Ecoli GitHub repository.

```bash
git clone https://github.com/barricklab/LTEE-Ecoli.git
``` 

!!! note
    The LTEE-Ecoli repository uses a specific version of _breseq_ (often not the newest one)
    to ensure compatibility between GenomeDiff files and _gdtools_ utility commands.

In this tutorial, we will run the _brefito_ workflow that assists with checking and further curating

## Copy input GenomeDiff files

Create a "curation directory" on a volume where you will have enough free disk space for _breseq_runs.

You MUST name your "curation directory" so that it contains the population designation in the format `Ara+[1-6]` or `Ara-[1-6]` somewhere (for example, `Ara-3-curation` of `LTEE-curation-Ara+4`) or `MAE*`. These are case sensitive. The naming scheme makes it possible for the pipeline to determine the correct ancestor genome to use.

```bash
mkdir Ara-3-curation
cd Ara-3-curation
``` 

Now, create a folder called `gd` within the "curation directory" (`Ara-3-curation` in the example):
```bash
mkdir gd
``` 

Guess what you put here? Correct! The GenomeDiff files directly output by your _breseq_ 
runs of new sample. Let's say those are located under `breseq-references/gd` because you already ran
the _brefito_ `predict-mutations-breseq` workflow in the `Ara-3-curation` directory.
 
Then you could use this command to copy them here:
```bash
cp breseq-references/gd/*.gd gd
```

You should also copy over some already curated GenomdDiff files from the same LTEE population. 
You can find these in the main repository folder `LTEE-clone-curated`. 

If you had that folder side-by-side with `Ara-3-curation` on your disk, you could copy over all of the 
ones from population A-3 this way:
```bash
cp ../LTEE-clone-curated/A-3*.gd gd
```

You will also need a _brefito_ `data.csv` (see the _brefito_ documentation) file in your curation directory. This file tells _brefito_ where to find the read and reference files that it will need for running your analyses. Files like this for each LTEE population will eventually be made available in the LTEE-Ecoli repository. An example of what one should look like is below.

```
sample,type,setting
Ara-2_500gen_763A,illumina-PE,sra://SRR2591050
Ara-2_500gen_763B,illumina-PE,sra://SRR2589073
Ara-2_1000gen_965A,illumina-PE,sra://SRR2584405
Ara-2_1000gen_965B,illumina-PE,sra://SRR2584465
*,reference,https://raw.githubusercontent.com/barricklab/LTEE/7da91974eafac0c5a8f903ae57275795d4395af2/reference/REL606.gbk
```

!!! note
    At first, you might want to only copy over a subset of the LTEE clones that is spread over generations, 
    because the pipeline will take longer the more you use. For final curation you should use them all!

## Initialize the curation directory

Run the LTEE-Ecoli command from your curation directory with the `ltee-ecoli` Conda environment activated.

```bash
brefito curate-LTEE-clones --keep-going
```
!!! note
    Notice that we added the _snakemake_ `--keep-going` flag here. This will cause the pipeline to continue to process
    a many steps as it can for as many genome as possible if one of you GenomeDiff files has a problem due to your edits.

The first time you run the `curate-LTEE-clones` workflow, it will create several directories (`00_header`, `01_curate_add`, `02_curate_remove`) that _are never written over in subsequent invocations of the workflow. You will edit the GenomeDiff files in these directories (and _only_ in these directories)for curating the LTEE mutation lists!

## Edit the header GenomeDiff files

Now edit the newly created `00_header` files for your genomes so they have additional metadata, 
such as `TREATMENT`, `TIME`, `POPULATION`, `CLONE`, and `MUTATOR_STATUS`.

Here's an example you can work from:
```text
#=GENOME_DIFF	1.0
#=TITLE	Ara-5_10000gen_4540A
#=AUTHOR	<yourname>
#=TIME	10000
#=POPULATION	Ara-5
#=TREATMENT	LTEE
#=CLONE	A
#=MUTATOR_STATUS	non-mutator
#=REFSEQ	https://raw.githubusercontent.com/barricklab/LTEE/7da91974eafac0c5a8f903ae57275795d4395af2/reference/REL606.gbk
#=READSEQ	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/001/SRR2589061/SRR2589061_1.fastq.gz
#=READSEQ	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/001/SRR2589061/SRR2589061_2.fastq.gz
```

!!! note
    Those are tabs (not spaces) separating the #=TAG and <value> part of each line.

## Checking your curation of an LTEE population

Now examine some of the other folders created by the _brefito_ `curate-LTEE-clones` workflow.

The ones that are most immediately useful are these:
- `04_final_normalized_gd` - contains the GenomeDiff files that are fully curated. Use these
 with `gdtools APPLY` to test your curation. Once they look complete, use them for your analyses.
- `compare_normalized.html` - the HTML output of `gdtools COMPARE` on the fully curated GenomeDiffs.
- `07_phylogeny/tree.rerooted.tre` - Newick format phylogenetic tree that can be viewed in programs
like [FigTree](https://github.com/rambaut/figtree/releases).

Most of the rest of the files are related to different ways of counting or not counting mutations.

!!! note
    **What's normalization?**
    Some mutations, such as a deletion of an `A` in a run of `AAAAA` could be annotated with 
    multiple positions. `gdtools NORMALIZE` is used to make all of these cases match, which
    is important for deciding whether the same mutation is present in multiple genomes.

!!! note
    **What's masking?**
    Some of the other files and directories mention "masking". What this means is that
    mutations in certain regions of the genome are not counted. Why? If you have a
    mixture of some datasets with longer reads and some with short reads, the short-read
    datasets will miss some mutations in repetitive regions that the others can fine.
    This leads to unequal counting of mutations and disrupts a phylogenetic tree.
    For the LTEE, we use a masking file that asks what regions of the genome
    one would be able to call mutations in with 36-bp reads (the shortest in any dataset).

## Edit the curate add and remove GenomeDiff files

Based on looking at unassigned evidence in your _breseq_ results, you will add mutations 
to the `01_curate_add` GenomeDiff file, as described in the rest of this tutorial. Sometimes you will divide up a mutation predicted by _breseq_ into several mutations that happen one after another such that, in the end, they create the same final change to the genome. Other times _breseq_ might have a false-positive prediction of a mutation due to errors or uncertaintly in the sequencing data. These are the cases when you will need to copy a mutation line from your original GenomeDiff file to the one with the same name in `01_curate_remove`.

The `compare_normalized.html` file is great for finding when later genomes have mutations that
hide earlier mutations or that are described in a different way. If you entered `TIME` and `CLONE`
metadate in your headers, the columns will be sorted by this information.

The directory `07_phylogeny/discrepancies` contains Newick tree files
for all mutations that _don't_ agree with the overall phylogenetic tree. 

Both files are helpful for noticing earlier mutations that were deleted in a later genome 
and need to be marked with `deleted=1`, mutations that were otherwise modified and need to be 
broken down into their constituent parts, and sometimes mutations that were just missed by 
_breseq_ and can be found upon further examination of the
sequencing data (for example with `breseq BAM2COV`, `breseq BAM2ALN`, or IGV).

Inescapably, you'll introduce some typos and other inconsistencies into your `01_curate_add` and `02_curate_remove` GenomeDiff files as you edit them that cause them to fail validation. This will cause steps in the `curate-LTEE-clones` workflow to fail. Look at the `*.log` files referenced in those steps to find the error messages that will tell you which entries have problems.

## Test your curation!

The gold-standard of curation is a "clean apply". What this means is that when we generate the evolved genome from the ancestral genome by applying all of the mutations in the GenomeDiff file and then requery the FASTQ reads against it using _breseq_, there should be _no mutations predicted_!

You can quicky tst your curation using a variation of the _brefito_ predict-mutations-breseq workflow that specifies using the evolved genome files produced as the reference for _breseq_ runs.

```bash
brefito predict-mutations-breseq-mutants --config BRESEQ_OPTIONS="-l 80"
```

Appending `-mutants` to the workflow name makes it use files in the `mutants` directory as the references for the run rather than the file specified in the `data.csv` parameter input file.

This will generate a new directory of _breseq_ output in the `breseq-mutants` directory. Open those file and look for discrepancies that indicate that you have further work to do in curating your GenomeDiff files.