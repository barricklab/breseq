[Back to the Main Curation Tutorial Page](tutorial-curation.md)

In this part of the tutorial, we will automate running the curation cycle on many genomes at once using [_brefito_](https://github.com/barricklab/brefito). 

**What's happening?** 

_brefito_ is a helper script and set of snakemake workflows that make it easier to run these steps on many samples in parallel. Once you set up a control file telling it what to do, it generates and runs all of the commands that you need!

## Installing _brefito_

Follow the [directions for installing _brefito_](https://github.com/barricklab/brefito/wiki/Installation).

Since _brefito_ is actively updated (and we might even fix bugs while doing this tutorial), consider using the **editable install** instructions. This will make it quicker to update your code to the newest version.

## Setting up the data file

To use _brefito_ you will specify a set of samples and the references files and read files that are associated with each sample. This is accomplished by creating a `data.csv` file in the main directory in which you will perform the analysis. For everything in this part of the tutorial we will be describing what you put inside this main directory.

First, read the [_brefito_ overview](https://github.com/barricklab/brefito/wiki/Overview) to learn a bit about the main command and the format of the `data.csv` file. 

A simple template for this file would look like this:
```
sample,type,setting
SAMPLENAME1,illumina-PE,READPATH1
SAMPLENAME2,illumina-PE,READPATH2
...
SAMPLENAMEX,illumina-PE,READPATHX
*,reference,REFERENCEPATH
```
But, we need to change the entries in ALL-CAPS!

The `data.csv` file is in a comma-separated values format. You can edit it in either a text editor or in a spreadsheet program. Be careful to export as a CSV with UTF-8 encoding if you use a spreadsheet program.

For example, here's how we would set up a run of the 75K-generation A and B clones from LTEE population A-1. You can run an analysis on many more than two samples at once using what you will learn here, but at least this example will start to show you how using _brefito_ can simplify your life.

For purposes of this tutorial, we are going to put all of the input **read files** we need in an `input` subdirectory and tell _brefito_ to use these files. If your reads are accessible online (for example, in the SRA or stored on a fileserver that can be accessed publicly or privately through `lftp`) there are other ways to set this up described in the documentation that will save you the step of manually having to copy the files to this location.

Since all of our samples are LTEE genomes and the _E. coli_ REL606 **reference file** is available online, we will give _brefito_ a URL to download that from. (The URL is complex because it is for a specific version from the [LTEE-Ecoli repository](https://github.com/barricklab/LTEE-Ecoli).)

One trick for paired-end Illumina reads is that putting the special wildcard string `{1|2}` inside of the `READPATH` makes it match both files with either a 1 or a 2 there, so you can easily specify both FASTQ files for paired-end reads in one line. It's also possible to put two lines, one with the name of each file, and as long as they have a 1 and corresponding 2 in the filenames and match otherwise, _brefito_ should be able to figure out that they belong together.

So your `data.csv` file might end up looking something like this:
```
sample,type,setting
A-1_75K_clone_A,illumina-PE,input/75K-1A_LTC-0000001_S.R{1|2}.fastq.gz
A-1_75K_clone_B,illumina-PE,input/75K-1B_LTC-0000002_S.R{1|2}.fastq.gz
*,reference,https://raw.githubusercontent.com/barricklab/LTEE/7da91974eafac0c5a8f903ae57275795d4395af2/reference/REL606.gbk
```

To get this to work, you MUST download/copy/move these four files into a new `input` directory that you create within your main directory:
```
75K-1A_LTC-0000001_S.R1.fastq.gz
75K-1A_LTC-0000001_S.R2.fastq.gz
75K-1B_LTC-0000002_S.R1.fastq.gz
75K-1B_LTC-0000002_S.R2.fastq.gz
```

If you had single-end reads, you would specify just `illumina` or `illumina-SE` in place of `illumina-PE`. If you want to combine analyzing FASTQ files from multiple instrument runs (perhaps you didn't get enough coverage the first time), just add extra lines for that sample. You can also specify `nanopore` reads and use these as input. 

!!! tip
    _brefito_ takes care of running the right programs to trim different types of reads and setting _breseq_ options to be appropriate. Usually, its default settings should work, but if you want more control, see the [_brefito_ workflow descriptions](https://github.com/barricklab/brefito/wiki/Workflows) for how you can pass setting to specific programs that are run as part of a pipeline.

You could test whether your `data.csv` is set up correctly by running the _brefito_ command that downloads and copies the input files to its particular locations (Remember to activate your _brefito_ conda environment to run this command!)
```bash
brefito download-data --resources connections=4
```
!!! note
    The `--resources connections=4` part is optional. It tells _brefito_ how many downloads it is allowed to initiate simultaneously. By default, it is 1. Don't set this too high if you are on a shared machine! Something in the range 4-8 can often speed up the downloads by quite a bit, but at some point there are diminishing returns because you only have so much bandwidth in your connection to the internet. If you are following the previous steps, it actually doesn't help to have multiple downloads, because there is just the one reference file that will be downloaded and the read files are already on your machine.

If everything worked, you should now have `reference` and `illumina-reads` subdirectories.

If you get an error, most likely there is something wrong with your paths! Check that the FASTQ files have the right names and are in the `input` directory within the main directory where you are running _brefito_ and have your `data.csv` file.

## Running _breseq_ on all of the samples

Now, running _breseq_ on all samples is as simple as running this command:

```bash
brefito predict-mutations-breseq --config BRESEQ_OPTIONS="-l 80"
```
The part after `BRESEQ_OPTIONS` gets sent to every _breseq_ command, so this limits the nominal coverage read-depth to 80x. (We have a lot of reads for these samples, so it will take a long time to run and a lot of disk space if you don't do this without any real gain in detecting mutations.) If you wanted to also use polymorphism mode, you could add `-p` inside of the double quotes.

!!! caution
    You should not pass the `-j` option for controlling how many threads _breseq_ uses in this way. The number of threads to use is specified by the Snakemake workflow in a different way.

The final output will be in `breseq-references` (_breseq_ run against the files in the `references` directory). 

!!! note
    Notice that a new directory called `illumina-reads-trimmed` is also created because this workflow first trims the reads. If you ever want to _just_ trim the reads, there are different _brefito_ sub-workflows you can use called `trim-illumina-reads` and `trim-nanopore-reads`. You can see all of the workflows by running `brefito --help`.

_brefito_ divides up the _breseq_ output into `data`, `gd`, and `html` subdirectories by sample, so that it is easy to grab just the pieces you need. One thing that is nice about the `gd` output directory is that instead of having to rename the `output.gd` from each _breseq_ run to `SAMPLE.gd` when you are curating them, _brefito_ already does this for you.

If something goes wrong, check the logs in the `log` subdirectory. This is where the output you would normally see going to the terminal is piped (because it would get confusing to have it jumbled together from many _breseq_ runs happening simultaneously).

!!! note
    In the future, we plan to add a say of changing _breseq_ options on a per-sample basis. For now, you should set up different _brefito_ analysis folders for using different options.

## Additional analysis examples

Some common analysis steps can also be automated by _brefito_! 

For example, this 
```bash
brefito compare-mutations-breseq
```
This creates a file `breseq-references/compare_1.html` with the output. If you had sets of samples with different reference sequences, you would have more output files with names ending 2, 3, 4, etc., for each set of samples that shared the same reference sequences.

!!! note
    If you know you want this final output, you can actually run this command FIRST, and _brefito_ will also run all of the required _breseq_ , read trimming, downloading, etc., steps it needs to get to this output!

This command will generate tiled and summary coverage plots for all samples.
```bash
brefito coverage-plots-breseq
```
The output for each sample is under `breseq-references/cov`.

## Curating the GenomeDiff files 

Copy the `breseq-references/gd` directory to a new `genome-diffs`.
```bash
cp -r breseq-references/gd genome-diffs
```

Edit these GenomeDiff files as described in other parts of this tutorial to curate their mutation lists.

When you are ready to check the mutation lists you annotated in your curated GenomeDiff files, it's always good to first quickly check their syntax with this command.

```bash
brefito validate-genome-diff-gdtools
```
The output will be placed in `genome-diffs-validate`. 

You can open each individual text file, or view them all at once by typing this command:
```bash
cat genome-diffs-validate/*
```

Then, you can generate the mutated genomes with this command.

```bash
brefito mutate-genomes-gdtools
```
This will create an output directory called `mutants` with the mutated reference genomes.

Then, re-run _breseq_ against these mutant genomes using this command.

```bash
brefito predict-mutations-breseq-mutants --config BRESEQ_OPTIONS="-l 80"
```

!!! note
    `predict-mutations-breseq` and some other workflows that use reference sequences in the `references` directory by default can be switched over to using other directories of reference files with this `predict-mutations-breseq-*` syntax. Here `*` is `mutants`, but you could create your own directory called whatever and substitute it here.

## Other _brefito_ Functions

### BLAST searches against the reference files

Repeat or near-repeat sequences in your reference genome can make interpreting certain parts of the _breseq_ output difficult (such as **JC** evidence items) when a sequence matches multiple places in the genome.

```
brefito search-blast CACCAAAACGTGCCGAGATGATCCTGTAACCATCATCAGTTGTGAAGTAGTGATTCACGACTTCAAGGCGCTTTTCAAAAGGGTATTTTGGCTTTGACATATTAGGGGCTATTCCATTTCATCGTCCAACAAAATGGGTGCAGTAC A-1_75K_clone_A
```

This creates a long HTML BLAST output file output file ending in `*.html` and a short tabular BLAST output file ending in `*.tsv` inside of the `blast-references` directory. The output of the short tabular file is also output directly so you can look at that immediately.

Instead of providing a DNA sequence at the command line, you can replace that argument with the path to a FASTA file:
```
brefito search-blast input.fasta A-1_75K_clone_A
```

The contents of the file `input.fasta` should be in FASTA format. So, it should look something like this:
```test
>input
CACCAAAACGTGCCGAGATGATCCTGTAACCATCATCAGTTGTGAAGTAGTGATTCACGA
CTTCAAGGCGCTTTTCAAAAGGGTATTTTGGCTTTGACATATTAGGGGCTATTCCATTTC
ATCGTCCAACAAAATGGGTGCAGTAC
```
You can include multiple sequences in this FASTA file.

You can also use the workflow `search-blast-mutants` to BLAST against one of the mutant genomes you generated versus the default of using the references.

!!! warning
    If you leave off the last sample argument (`A-1_75K_clone_A` in the examples above) _brefito_ will perform individual BLAST runs against ALL of your reference sequences. Usually this will be redundant if they all have the same references! If you do want to BLAST against multiple references, but only specific ones, you can continue the command line by adding additional sample names (e.g., `A-1_75K_clone_B`)

**Next:** [Common curation cases](tutorial-curation-common-cases.md)

