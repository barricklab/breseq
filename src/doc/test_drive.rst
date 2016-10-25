.. _test-drive:

Test Drive
==============

In this test drive, we will first download a bacterial genome and FASTQ files of Illumina reads. Then, we will use |breseq| to predict mutations that are present in the re-sequencing data relative to this reference genome.

1. Download data files
---------------------------------

First, create a directory called **test_drive**:

.. code-block:: bash

   $ mkdir test_drive
   $ cd test_drive

Reference sequence
++++++++++++++++++++

|breseq| prefers the reference sequence in Genbank or GFF3 format. In this example, the reference sequence is *Escherichia coli* B strain REL606. The Genbank (Refseq) accession number is: **NC_012967**. You can search for this sequence at http://www.ncbi.nlm.nih.gov/ or follow this `direct link <http://www.ncbi.nlm.nih.gov/nuccore/NC_012967>`_.

Once the sequence is displayed, you will want to select "Show sequence" from the Display options on the right then click "Update View" and let the sequence download complete. Finally, use the "Send:" menu to choose "Complete Record" and Destination: "File" and "Genbank (Full)". It should start downloading a file called "sequence.gb". Rename this to **NC_012967.gbk** after it downloads.

.. warning::
   
   A common error in using |breseq| is to download and try to use a GenBank file that does not include the DNA sequence of the genome. Remember to "Show sequence" from the Display options on the right then click "Update View" before downloading to avoid this problem!

If you open the GenBank file that you downloaded in a text editor, you should see a section with ORIGIN followed by the DNA sequence of the genome, like this:

.. code-block:: text

   ORIGIN
                 1 agcttttcat tctgactgca acgggcaata tgtctctgtg tggattaaaa aaagagtgtc
                61 tgatagcagc ttctgaactg gttacctgcc gtgagtaaat taaaatttta ttgacttagg
               121 tcactaaata ctttaaccaa tataggcata gcgcacagac agataaaaat tacagagtac
               181 acaacatcca tgaaacgcat tagcaccacc attaccacca ccatcaccat taccacaggt
               241 aacggtgcgg gctgacgcgt acaggaaaca cagaaaaaag cccgcacctg acagtgcggg

Read files
++++++++++++++

We're going to use Illumina genome re-sequencing data from a strain that evolved for 20,000 generations in a long-term evolution experiment [Barrick2009a]_. This data is available in the European Nucleotide Archive (ENA). Go to http://www.ebi.ac.uk/ and search for the accession number: **SRR030257**. Then click on the accession number to open the record and download the two data files using the links in the 'ftp' column.

Move all three of these files into the **test_drive** directory that you created.

2. Run |breseq|
-----------------------

Check to be sure that you have changed into the **test_drive** directory and that you have all of the input files (and have uncompressed them).

.. code-block:: bash

   $ ls 
   NC_012967.gbk		SRR030257_1.fastq	SRR030257_2.fastq

Now, run breseq:

.. code-block:: bash

   $ breseq -r NC_012967.gbk SRR030257_1.fastq SRR030257_2.fastq

The first named argument (-r) is the reference sequence. If you had multiple reference sequences, you could input multiple ones (e.g., -r NC_012967.gbk -r plasmid.gbk).

The unnamed arguments at the end of the command line are the read files. You can input as many as you need and mix FASTQ files from different sequencing technologies (e.g., Illumina and 454).

.. warning::
   
   Running |breseq| on a full data set like this is not fast. Depending on your computer, this could take several hours. To speed things up, you can set -j option to the number of cores on your machine to enable multithreaded execution of some steps (e.g., -j 4 for a quad-core machine). If you want to speed this example up, you might also include only one of the two input read files on the command line.


3. Open |breseq| output
----------------------------

Open the file **index.html** in the new **output** directory. This describes the predicted mutations and also evidence for mutations that |breseq| could not resolve into mutational events. The tables in this HTML file are more fully described in the section on :ref:`output-format`.
