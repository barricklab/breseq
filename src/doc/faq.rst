Frequently Asked Questions (FAQ)
==================================

.. contents::
   :local:

Does |breseq| use paired-end or mate-paired information?
----------------------------------------------------------
Short answer: No. 

Your reads are mapped in single-end mode even if you are using paired-end or mate-paired data. 
For most microbial genomes, you don't gain much sensitivity (in terms of the number of reference positions at
which there are enough uniquely-mapped reads to call mutations) by doing paired read mapping. Furthermore,
the split-read analysis approach that  |breseq| uses to discover new sequence junctions is more precise 
(finding exact sequence breakpoints) and generally at least as sensitive as predicting structural variation by 
examining read pairs that are mapped with anomalous orientations and insert sizes.

That said, there are definitely cases where this information can be useful, especially if one has data with 
large insert sizes (e.g. an Illumina mate paired library. So, we hope to include some of this functionality 
in the future. It just hasn't been a very high priority.

Why do mutations have a predicted frequency of 100% when I see a mixture in the read pileup?
--------------------------------------------------------------------------------------------
By default, |breseq| is run in **consensus mode** in which it assumes you have a pure clonal
sample of a haploid genome. It therefore uses a statistical model that will only predict 0% 
or 100% for the frequency of each mutation on the main results page. For SNPs and small indels, 
it *does* test a very conservative mixed model that allows intermediate frequencies. If
this model is a better fit to the data, then the mutation will be demoted to the **marginal prediction**
page because |breseq| assumes that it is some kind of artifact in your sample (a sequencing error
hotspot or due to your sample actually being a mixture of two different clones, for example).

When |breseq| is run in **polymorphism mode** (by supplying the **--polymorphism-prediction|-p** option)
it uses a statistical approach that finds the maximum likelihood best frequency over the entire range 0-100%. 
If the prediction of a polymorphism fails some filtering steps (controlled by options) than it is rejected as a
polymorphism and changed to a consensus prediction with 0% or 100% frequency.

Why am I seeing a lot of junction (JC) predictions and no mobile element (MOB) mutations?
---------------------------------------------------------------------------------------------
|breseq| needs to recognize how mobile elements (such as bacterial insertion sequences) are annotated in the
reference genome file. The syntax for doing this is less standardized than it is for annotating genes, but |breseq|
tries to be flexible. For predicting MOB mutations, it will use any annotation item with a major type
of ``repeat_region`` or ``mobile_element``. To determine the name of the mobile element, it will look for a sub-tag in 
this annotation item matching: ``label``, ``note``, ``mobile_element``, ``mobile_element_type``, or ``rpt_family``. For
proper predictions of MOBs by matching together JC predictions, all elements in the same family should have the same name.
  
For example, here are two common ways that IS elements are annotated in GenBank files. Both are recognized by |breseq|.
 
.. code-block:: text

   repeat_region   15387..16731
                   /mobile_element="insertion sequence:IS186A"
   
   mobile_element  321573..322809
                   /mobile_element_type="insertion sequence:IS1236"

If you find another common way that mobile elements are annotated in your input files, feel free to add a feature 
request on GitHub for |breseq| to recognize this format.
