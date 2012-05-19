Usage Summary
==============

:program:`breseq`
------------------

Usage::

  breseq -r reference1.gbk [-r reference2.gbk ...] reads1.fastq [reads2.fastq, reads3.fastq...]

.. program:: breseq

Run the :program:`breseq` pipeline for mutation prediction from genome re-sequencing data.

Required options:

.. option:: -r <file_path>, --reference <file_path> 

   Input reference genome sequence files in GenBank format. If there are multiple reference sequences stored in separate GenBank files (e.g., a bacterial genome and a plasmid), this option can be supplied multiple times.

.. option:: reads1.fastq [reads2.fastq, reads3.fastq...]  

   The remaining arguments at the command line are the FASTQ input files of reads. FASTQ files with base quality scores that are not in `SANGER format <http://en.wikipedia.org/wiki/FASTQ_format>`_ will be converted. In addition, reads with >50% N bases will be removed from the converted FASTQ file by default. |breseq| re-calibrates the error rates for each FASTQ file separately, so data sets that were generated independently should be stored in different input files.

Expert options:

.. option:: --base-quality-cutoff=<int>

   Ignore bases with a quality score lower than this value when calling mutations. This accommodates Illumina formats that use quality scores of 2 to flag bad data. These bases are still used for aligning to the reference genome and are shown highlighted in yellow when drawing alignments, but they do not contribute to read alignment evidence. Default: 3

.. option:: --require-complete-match

   Require the entire read to match for it to be counted as aligned.

.. option:: --required-match-length  

   Require at least this number of nucleotides in the read to match for it to be counted as aligned.
   
.. option:: --predict-polymorphisms

   Identify and predict the frequencies of SNPs and small indels that are polymorphic (appear in only a subpopulation of reads). See :ref:`polymorphism-prediction` for additional options and note that this option is still experimental.

Command: bam2aln
--------------------------

Usage::

   breseq bam2aln [-b input.bam] [-f input.fasta] [-o output/path] region1 [region2 region3 ...]

.. program::`breseq bam2aln`

Creates HTML pileup files displaying reads aligned to each specified region.

Options:

.. option:: -b <file_path>, --bam=<file_path> 

   BAM database file of read alignments. Defaults: reference.bam, data/reference.bam.

.. option:: -f <file_path>, --fasta=<file_path> 

   FASTA file of reference sequences. Defaults: reference.fasta, data/reference.fasta.
   
.. option:: -o <path>, --output=<path> 

   Output path. If there are multiple regions, must be a directory path, and all output files will be output here with names region1.html, region2.html, ... If there is just one region, the output file will be given this name if it is not the name of an already existing directory. Default: current path.
   
.. option:: -n <int>, --max-reads=<int>

   Maximum number of reads that will be aligned to a region. If there are more than this many reads, then the reads displayed are randomly chosen and a warning is added to the output. Default: 1000.

.. option:: region1 [region2, region3, ...]

   Regions to create output for must be provided in the format **FRAGMENT:START-END**, where **FRAGMENT** is a valid identifier for one of the sequences in the FASTA file, and **START** and **END** are 1-indexed coordinates of the beginning and end positions. Any read overlapping these positions will be shown. A separate output file is created for each region.

Command: bam2cov
--------------------------

Usage::

  breseq bam2cov [-b input.bam] [-f input.fasta] [-o output/path] region1 [region2 region3 ...]

.. program::`breseq bam2cov`

Create a coverage plot or table for the specified region or regions.

Options:

.. option:: -b <file_path>, --bam=<file_path> 

   BAM database file of read alignments. Defaults: reference.bam, data/reference.bam.

.. option:: -f <file_path>, --fasta=<file_path> 

   FASTA file of reference sequences. Defaults: reference.fasta, data/reference.fasta.
   
.. option:: -o <path>, --output=<path> 

   Base name of output files. Region specification (seq_id:start-end) appended if there are multiple output files. Default: seq_id:start-end for single regions or the current directory for multiple regions.

.. option:: --plot-format=<plot_format> 

   Format of output plot: PNG or PDF. Default: PNG
   
.. option:: -t, --table

   Create tab delimited file of coverage instead of a plot.

.. option:: -1, --total-only

   Only plot/tabulate total coverage, not per strand coverage.
   
.. option:: --resolution=<int>

  Number of positions to output coverage information for in interval (0=ALL). Default: 600

.. option:: region1 [region2, region3, ...]

   Regions to create output for must be provided in the format **FRAGMENT:START-END**, where **FRAGMENT** is a valid identifier for one of the sequences in the FASTA file, and **START** and **END** are 1-indexed coordinates of the beginning and end positions. A separate output file is created for each region.
