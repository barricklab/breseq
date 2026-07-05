In this test drive, we will first download a bacterial genome and FASTQ
files of Illumina reads. Then, we will use _breseq_ to predict
mutations that are present in the resequencing data relative to this
reference genome.

# 1. Download data files

First, create a directory called `test_drive`:

``` bash
$ mkdir test_drive
$ cd test_drive
```

## Reference sequence

_breseq_ prefers the reference sequence in Genbank or GFF3 format.
In this example, the reference sequence is *Escherichia coli* B strain
REL606. The Genbank (Refseq) accession number is: `NC_012967`. You can
search for this sequence at <https://www.ncbi.nlm.nih.gov/> or follow
this [direct link](https://www.ncbi.nlm.nih.gov/nuccore/NC_012967).

Once the sequence is displayed, you will want to select "Show sequence"
from the Display options on the right then click "Update View" and let
the sequence download complete. Finally, use the "Send:" menu to choose
"Complete Record" and Destination: "File" and "Genbank (Full)". It
should start downloading a file called `sequence.gb`. Rename this to
`NC_012967.gbk` after it downloads.

!!! warning
    A common error in using _breseq_ is to download and try to use a
    GenBank file that does not include the DNA sequence of the genome.
    Remember to "Show sequence" from the Display options on the right then
    click "Update View" before downloading to avoid this problem!

If you open the GenBank file that you downloaded and search or scroll way down in a text editor, you should see a section with `ORIGIN` followed by the DNA sequence of the
genome, like this:

``` text
ORIGIN
              1 agcttttcat tctgactgca acgggcaata tgtctctgtg tggattaaaa aaagagtgtc
             61 tgatagcagc ttctgaactg gttacctgcc gtgagtaaat taaaatttta ttgacttagg
            121 tcactaaata ctttaaccaa tataggcata gcgcacagac agataaaaat tacagagtac
            181 acaacatcca tgaaacgcat tagcaccacc attaccacca ccatcaccat taccacaggt
            241 aacggtgcgg gctgacgcgt acaggaaaca cagaaaaaag cccgcacctg acagtgcggg
```

## Read files

We're going to use Illumina genome resequencing data from a strain that
evolved for 20,000 generations in the Lenski long-term evolution experiment
(Barrick et al. 2009). This data is available in the European Nucleotide
Archive (ENA). Go to <https://www.ebi.ac.uk/> and search for the
accession number: `SRR030257`. Then click on the accession number to
open the record and download the two data files using the links in the
'ftp' column.

Move all three of these files into the `test_drive` directory that you
created.

# 2. Run _breseq_

Check to be sure that you have changed into the **test_drive** directory
and that you have all of the input files (and have uncompressed them).

``` bash
$ ls
NC_012967.gbk        SRR030257_1.fastq.gz   SRR030257_2.fastq.gz
```

Now, run breseq:

``` bash
$ breseq -r NC_012967.gbk SRR030257_1.fastq.gz SRR030257_2.fastq.gz
```

The first named argument (-r) is the reference sequence. If you had
multiple reference sequences, you could input multiple ones (e.g., -r
NC_012967.gbk -r plasmid.gbk).

The unnamed arguments at the end of the command line are the read files.
You can input as many as you need and mix FASTQ files from different
sequencing technologies (e.g., Illumina and 454).

For the unnamed read argument(s) but not for the reference sequence 
arguments you can use globs that match multiple files. For example, in
this case, you could use `*.fastq.gz` to match both input FASTQ files.

!!! warning
    Running _breseq_ on a full data set like this is not fast. Depending
    on your computer, this could take several hours. To speed things up, you
    should set the `-j` option to the number of cores on your machine to enable
    multithreaded execution of some steps (e.g., `-j 4` for a quad-core
    machine). If you want to speed this example up, you might also include
    only one of the two input read files on the command line.

# 3. Open _breseq_ output

Open the file **index.html** in the new **output** directory. This
describes the predicted mutations and also evidence for mutations that
_breseq_ could not resolve into mutational events. The tables in
this HTML file are more fully described in the section on
`output-format`.
