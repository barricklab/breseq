.. _tutorial-populations:

Tutorial: Population Samples (Polymorphism Mode)
================================================

In this exercise, you will analyze two population (metagenomic) samples using |breseq| to track the frequencies of evolved alleles and changes in genetic diversity in population Ara-3 of the Lenski long-term evolution experiment (LTEE). As discussed in :ref:`tutorial-clones` this population evolved citrate utilization after 31,500 generations.

.. note:: 
   This tutorial was created for the EMBO Practical Course `Measuring intra-species diversity using high-throughput sequencing <http://events.embo.org/15-htp-sequencing/>`_ held 27–31 July 2015 in Oeras, Portugal.

.. warning::

   If you encounter any |breseq| or |gdtools| errors or crashes in running this tutorial, please `report issues on GitHub <https://github.com/barricklab/breseq/issues>`_.

1. Download data files
---------------------------------

First, create a directory called **tutorial_populations**:

.. code-block:: bash

   $ mkdir tutorial_populations
   $ cd tutorial_populations

Reference sequence
++++++++++++++++++++

|breseq| prefers the reference sequence in Genbank or GFF3 format. In this example, the reference sequence is the ancestral *Escherichia coli* B strain REL606. For this tutorial, we are going to use an updated version of the GenBank genome record (accession:**NC_012967**) that contains richer gene annotation information than the version available from NCBI.

* `Download REL606.gbk via this link <http://barricklab.org/release/breseq_tutorial/REL606.gbk.gz>`_

Read files
++++++++++++++

We're going to use Illumina genome re-sequencing data from mixed populations that evolved for up to 40,000 generations in a long-term evolution experiment [Blount2008]_ [Blount2011]_. This data is available in the European Nucleotide Archive (ENA). Go to http://www.ebi.ac.uk/ and search for the accession number: **SRR1721884**. Then click on the accession number to open the record and download the two FASTQ files using the links in the 'ftp' column. 

This particular sample was taken at 20,000 generations from population Ara-3. You'll use it to illustrate running |breseq| in polymorphism mode and the consequences of different filtering options for ruling out false-positives. If you would like to access the entire time-series of population samples from this population check out `SRP051254 <http://www.ebi.ac.uk/ena/data/view/SRP051254>`_.

2. Run |breseq| with default filters
----------------------------------------

Check to be sure that you have changed into the **tutorial_clonal_samples** directory and that you have all of the input files (and have uncompressed them).

.. code-block:: bash

   $ ls 
   NC_012967.gbk   SRR030257_1.fastq   SRR030257_2.fastq

Now, run |breseq| using this command:

.. code-block:: bash

   $ breseq -j 8 -p -o REL8595M-default -r REL606.gbk SRR030257_1.fastq SRR030257_2.fastq

This command is expected to take roughly 30 minutes to an hour to complete.

.. note::

   If you are unable to complete this command, please download the `output for REL8595M-default <http://barricklab.org/release/breseq_tutorial/REL8595M-default.tgz>`_ to continue the tutorial.

Open ``REL8595M-default/output/index.html``. Examine the mutation lines that are highlighted in green, which are predicted to be polymorphic in the mixed population sample (the predicted allele has an estimated maximum likelihood frequency between 0% and 100%). Click on some of the **RA** and **JC** links for these items. How are they different from those that you observed when |breseq| was run in consensus mode on a clone?

3. Run |breseq| with no filters
----------------------------------

Predicting polymorphisms is very prone to false-positives (wrongly predicting genetic variation that is not actually present in a sample). This is, in part, because sequencing data has hotspots and biases that are difficult to adequately capture in a statistical error model (particularly when analyzing only one sample at a time, like in this example). |breseq| has some default options that can be used to filter out variant predictions that look suspicious because of certain biases that they exhibit with respect to how reads align to them versus how they typically align to "normal" positions in the genome. 

In consensus mode, |breseq| results are generally robust to different sequencing coverage depths and types of samples. In polymorphism mode, |breseq| often needs some tuning of parameters and statistical cutoffs depending on characteristics of the input data set in order to not predict either too many (false-positives) or too few (false-negatives) polymorphisms. In addition, it may be necessary to perform more complex analyses of multiple samples or of time courses to gain extra power for discriminating true polymorphisms from errors. Unfortunately, these approaches are outside the scope of a simple |breseq| command!

Bring up the full |breseq| help:

.. code-block:: bash

   $ breseq -h

The relevant options are listed under **Polymorphism (Mixed Population) Options**. Now, we're going to do a |breseq| run in which we disable all of the filters for comparison to the initial run:

.. code-block:: bash

   $ breseq -j 8 -p --polymorphism-reject-indel-homopolymer-length 0 --polymorphism-reject-surrounding-homopolymer-length 0 --polymorphism-bias-cutoff 0 --polymorphism-minimum-coverage-each-strand 0 -o REL8595M-no-filters -r REL606.gbk SRR030257_1.fastq SRR030257_2.fastq

This command is expected to take roughly 30 minutes to an hour to complete.

.. note::

   If you are unable to complete this command, please download the `output for REL8595M-no-filters <http://barricklab.org/release/breseq_tutorial/REL8595M-no-filters.tgz>`_ to continue the tutorial.

4. Compare predictions of mutations
-----------------------------------

Open ``no-filters/output/index.html``. See if you can find examples of mutations that are probably due to different types of sequencing biases by delving into the original |breseq| HTML files that show the read alignments (RA).

You might first want to create a comparison table of the results from the two |breseq| runs.

.. code-block:: bash

   $ cp REL8595M-default/output/output.gd default.gd
   $ cp REL8595M-no-filters/output/output.gd no-filters.gd
   $ gdtools COMPARE -o compare.html -r REL606.gbk default.gd no-filters.gd

Can you find any predictions that look like plausible mutations that were incorrectly rejected by the default filters? 

.. container:: toggle

   .. container:: header

      **Hint**

   .. container:: text

      Look for mutations with intermediate predicted frequencies (closest to 50%).

5. Examine allele frequency time courses
----------------------------------------

Since it would take a long time to create results for all of the mixed population samples, download these |GenomeDiff| output files pre-generated with |breseq| using the default polymorphism filtering options in order to continue the tutorial:

* `Download population_gd <http://barricklab.org/release/breseq_tutorial/population_gd.tgz>`_

If you look at these files, you will also notice that metadata (experiment, population, generation) has been added to these files that enables them to be properly sorted into order.

Make a compare table for all of these files.

.. container:: toggle

   .. container:: header

      **Show me the commands**

   .. code-block:: bash

      $ cd population_gd
      $ gdtools COMPARE -r ../REL606.gbk -o ../time-course.html `ls *.gd`

Open the HTML output file and look at the trajectories of mutations that appear early and later.

Here are a few questions to get you started thinking about the data:

1. What do you notice about the last sample, REL11151 (45,000 generations)? It's located furthest to the right.

.. container:: toggle

   .. container:: header

      **See the answer**

   .. container:: text

      It has many, many, many more predicted mutations than the other samples. This type of result could potentially be due to some pathological characteristic of that particular sequencing dataset. However, in this case it is actually because the population evolved an elevated mutation rate by 36,000 generations which led to an explosion of genetic diversity in the population by 45,000 generations.

2. What other potential problems do you notice with the output?

.. container:: toggle

   .. container:: header

      **For example**

   .. container:: text

      Some mutations may "blink in and out of existence" (be present at one time point and then disappear at later time points only to reappear later). In some cases, this represents actual population dynamics: that lineage may have been close to extinction and then experienced a resurgence as it accumulated additional beneficial mutations. In many cases, however, we know that this is impossible because it does not happen to linked mutations that are present in the same evolved lineage. This type of error is due to improperly rejecting evidence for a polymorphism in one or a few samples.

      One solution to this problem is to adjust the |breseq| options for filtering polymorphism predictions, but this is unlikely to give a clean result for any setting. A second solution is to compile a list of evidence that you force |breseq| to always look at and predict the frequency of using the **--user-evidence-gd** option. If you supply this option, then it will predict and record evidence for that RA or JC item even if it has a frequency of 0% because there are no variants supporting it in a given population sample.
