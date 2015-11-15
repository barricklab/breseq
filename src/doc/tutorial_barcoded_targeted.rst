.. _tutorial-barcoded-targeted:

Tutorial: Ultra-rare variant detection using consensus reads and targeted sequencing
==========================================================================================

In this exercise, you will analyze targeted sequencing of the initial burst of genetic diversity in a short *E. coli* evolution experiment. This tutorial uses data prepared by a special library preparation technique that adds a "molecular index" to each initial DNA fragment. This enables one to sequence many amplification products from this initial read to achieve lower error rates. In addition pulldowns with biotinylated oligos were used to enrich for only certain genes in the E. coli genome to achieve deeper sequencing of regions that were expected to have beneficial mutations.

.. note:: 
   This tutorial was created for the EMBO Practical Course `Measuring intra-species diversity using high-throughput sequencing <http://events.embo.org/15-htp-sequencing/>`_ held 27–31 July 2015 in Oeras, Portugal.

.. warning::

   If you encounter any |breseq| or |gdtools| errors or crashes in running this tutorial, please `report issues on GitHub <https://github.com/barricklab/breseq/issues>`_.

1. Download data files
---------------------------------

First, create a directory called **sscs_targeted**:

.. code-block:: bash

   $ mkdir tutorial_barcoded_targeted
   $ cd tutorial_barcoded_targeted

Reference sequence
++++++++++++++++++++

The samples sequenced were genomic DNA from populations evolved from a clonal isolate of *Escherichia coli* B strain REL606. We'll contrast some new ways of analyzing this data that require us to use different reference genome setups. 

First, download the entire reference genome.

`Download REL606.gbk via this link <http://barricklab.org/release/breseq_tutorial/REL606.gbk.gz>`_

Then, download these two reference files:

* `on-target.gff3 <http://barricklab.org/release/breseq_tutorial/on-target.gff3.gz>`_
* `off-target.gff3 <http://barricklab.org/release/breseq_tutorial/off-target.gff3.gz>`_

In the first of these (``on-target.gff3``), we've extracted just the 8 target genes with 1400 bp added on each side as three different reference fragments. Each of these reference sequences were created by using ``gdtools APPLY`` to delete the rest of the genome.

In the second of these (``off-target.gff3``), we've masked (via N's) the targeted regions in the whole genome sequence. This reference sequence was created by using ``gdtools MASK``.

This separation of the reference sequence and the 8 targeted regions into different files and (artificial) DNA fragments will makes it easier to calculate certain statistics about enrichment of the target region relative to the rest of the genome and to ignore mutations that occur outside the regions of interest

Read files
++++++++++++++

Paired-end Illumina reads for one population sample taken at generation 139 from a time-course:

* `DED234_GATCAG_L004_R1_001.fastq <http://barricklab.org/release/breseq_tutorial/DED234_GATCAG_L004_R1_001.fastq.gz>`_
* `DED234_GATCAG_L004_R2_001.fastq <http://barricklab.org/release/breseq_tutorial/DED234_GATCAG_L004_R2_001.fastq>`_

2. Generate SSCS Reads
-----------------------

First, we need to pre-process the reads to construct single-strand consensus reads and remove the molecular barcodes. If you have numpy and other Python prerequisites installed, you can do this by downloading this script:

* `Download SSCS_DCS.py <http://barricklab.org/release/breseq_tutorial/SSCS_DCS.py.gz>`_

And then running this command:

.. code-block:: bash

   python SSCS_DCS.py -f1 DED234_GATCAG_L004_R1_001.fastq -f2 DED234_GATCAG_L004_R2_001.fastq -p DED234 -s -d -m 2 --log SSCS_Log

This script will take about 30 minutes to run.

.. warning::

   This script is memory intensive! (16 GB RAM required)

3. Run |breseq| on the consensus reads
-----------------------------------------------------

Use this command to analyze the consensus reads.

.. code-block:: bash

   $ breseq -t -j 8 -o consensus_reads -p --polymorphism-minimum-coverage-each-strand 0 --polymorphism-bias-cutoff 0 --polymorphism-score-cutoff 0 --polymorphism-reject-indel-homopolymer-length 0 --polymorphism-reject-surrounding-homopolymer-length 0 -r on-target.gff3 -s off-target.gff3 DED234_SSCS.fastq

Notice the new **-t** and **-s** options. What are these doing? Look in the help under **Reference File Options**.

Take a look at the |breseq| output. In particular, examine the **summary statistics** to look at the distribution of coverage across the different reference sequence fragments that were used. How effective was the enrichment step at favoring reads from those regions?

4. Run |breseq| on the original reads
-----------------------------------------------------

Your next task is to compare the performance of the consensus reads to the original reads. To make this comparison fair, you need to include the same numbers of original reads and consensus reads. Do this by extracting the same number of lines from the R1 file as exist in the new ``DED234_SSCS.fastq`` file.

.. container:: toggle

   .. container:: header

      **How? I need a hint**

   .. container:: text

      Check out the **wc -l** and **head** unix commands.

You also need to trim the molecular barcodes from these reads (this was done automatically by the ``SSCS_DCS.py`` script for us before).

Then, run a similar |breseq| command to the one above to generate an ``original_reads`` output directory.

Compare the overall |breseq| predictions on each data set (possibly by making a comparison table, as we have in the previous tutorials).

Also, take a look at a part of the |breseq| output that you may not have examined yet. On the **summary statistics** page, click on the link named **errors** in the **Read File Information** table at the top. This graph shows the sequencing error rates in the input reads. Compare the results for the consensus and original reads.
