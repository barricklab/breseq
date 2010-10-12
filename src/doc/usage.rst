Usage
==============

:program:`breseq`
------------------

Usage:

   :program:`breseq` -r reference.gbk reads1.fastq [reads2.fastq, reads3.fastq...]

.. program:: breseq

Run the :program:`breseq` pipeline for mutation prediction from genome re-sequencing data.

Required options:

.. option:: -r <file_path>, --reference <file_path> 

   Input reference genome sequence files in GenBank format. If there are multiple reference sequences stored in separate GenBank files (e.g., a bacterial genome and a plasmid), this option can be supplied multiple times.

.. option:: reads1.fastq [reads2.fastq, reads3.fastq...]  

   The remaining arguments at the command line are the FASTQ input files of reads. The FASTQ base quality scores must be in `SANGER format <http://en.wikipedia.org/wiki/FASTQ_format>`_. If you get an error and need to convert your quality scores, see the :ref:`fastq-utils` command. |breseq| re-calibrates the error rates for each FASTQ file separately, so data sets that were generated independently should be stored in different input files.

Expert options:

.. option:: --base-quality-cutoff=<int>

   Ignore bases with a quality score lower than this value when calling mutations. This accommodates Illumina formats that use quality scores of 2 to flag bad data. These bases are still used for aligning to the reference genome and are shown highlighted in yellow when drawing alignments. Default: 3


:program:`bam2aln`
------------------

Usage:

   :program:`bam2aln` [-b input.bam] [-f input.fasta] [-o output/path] region1 [region2 region3 ...]

.. program:: bam2aln

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


:program:`bam2cov`
------------------

Usage:

   :program:`bam2cov` -b input.bam -f input.fasta -o [output/path] region1 [region2, region3, ...]


.. program:: bam2cov

Creates a coverage plot or table for the specified region.
   
Options:

.. option:: -b <file_path>, --fasta <file_path> 

   BAM database file of read alignments. Defaults: reference.bam, data/reference.bam

.. option:: -f <file_path>, --fasta <file_path> 

   FASTA file of reference sequences. Defaults: reference.fasta, data/reference.fasta
   
.. option:: -o <path>, --output <path> 

   Output path. If there are multiple regions, must be a directory path, and all output files will be output here with names region1, region2, ... If there is one region, the output file will be given this name if it is not the name of an already existing directory. Default: current path.

.. option:: region1 [region2, region3, ...]

   Regions to create output for must be provided in the format **FRAGMENT:START-END**, where **FRAGMENT** is a valid identifier for one of the sequences in the FASTA file, and **START** and **END** are 1-indexed coordinates of the beginning and end of the region. A separate output file is created for each region.
   
.. option:: --pdf

   In plot mode, create output plot in PDF format rather than PNG format.

.. option:: -r <int>, --resolution <int>

   In plot mode, maximum mumber of reference positions to plot coverage for within the region. Default: 600.

.. option:: -1, --total_only

   In plot mode, only output the total coverage of unique or repeat read mappings. (Does not break these down into the coverage on each strand of the reference sequence.)

.. option:: -t, --table

   Table mode. Rather than a plot, output a tab-delimited table of the coverage in the specified region to the output file. Also outputs the mean and standard error of the unique coverage within each region to STDOUT.
   

.. _fastq-utils:

:program:`fastq_utils`
-----------------------

Usage:

   :program:`fastq_utils` COMMAND [arguments]

.. program:: fastq_utils

Performs various functions on FASTQ formatted files. Options depend on the COMMAND supplied. There are several different `FASTQ styles <http://en.wikipedia.org/wiki/FASTQ_format>`_ with different base quality score formats.

Command: FORMAT

Usage:

   :program:`fastq_utils` FORMAT [-n 1000|ALL] input.fastq 

Examine reads in a FASTQ file to predict its base quality score format.

.. option:: -n <int>, -n ALL, --num=<int>, --num=ALL

   Number of reads to examine when predicting the format. The keyword 'ALL' means examine every read in input the file.

.. option:: input.fastq

   FASTQ file to examine.

Command: SANGER

Usage:

   :program:`fastq_utils` SANGER -f from_format [-l] input.fastq output.fastq

Convert a FASTQ file to SANGER format.

.. option:: -f <format>, --format=<format>

   Base quality score format of the input FASTQ file. Valid formats are: SANGER, SOLEXA, ILLUMINA_1.3+, ILLUMINA_1.5+. If you are unsure of the format, use the FORMAT command.

.. option:: -l, --list-format

   In the input FASTQ file, quality score lines are white space separated numbers, rather than character strings. 

.. option:: input.fastq

   Input FASTQ file in specified format.

.. option:: output.fastq

   Output FASTQ file in SANGER format.
