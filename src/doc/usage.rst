Usage Summary
==============

:program:`breseq`
------------------

Usage::

  breseq -r reference1.gbk [-r reference2.gbk ...] reads1.fastq [reads2.fastq, reads3.fastq...]

.. program:: breseq

Run the :program:`breseq` mutation prediction pipeline.

Required options:

.. option:: -r <file_path>, --reference <file_path> 

   Input reference genome sequence files in GenBank, GFF3, or FASTA format. If there are multiple reference sequences stored in separate files (e.g., a bacterial genome and a plasmid), this option can be supplied multiple times.

.. option:: reads1.fastq [reads2.fastq, reads3.fastq...]  

   The remaining arguments at the command line are the FASTQ input files of reads. FASTQ files with base quality scores that are not in `SANGER format <http://en.wikipedia.org/wiki/FASTQ_format>`_ will be converted. In addition, reads with >50% N bases will be removed from the converted FASTQ file by default. |breseq| re-calibrates the error rates for each FASTQ file separately, so data sets that were generated independently should be stored in different input files.

Commonly used options:

.. option:: -h, --help

   Produce help message showing advanced options.

.. option:: -n <string>, --name <string>

   Human-readable name of the analysis run for output (DEFAULT=<none>).

.. option:: -j <int>, --num-processors <int>

   Number of processors to use in multithreaded steps (DEFAULT=1).

.. option:: --no-junction-prediction

   Do not predict new sequence junctions.
   
.. option:: -p, --predict-polymorphisms

   Predict polymorphic (mixed) mutations.

For a complete list of options, please see the command line help (by using the -h option).

Command: bam2aln
--------------------------

Usage::

  breseq BAM2ALN [-b reference.bam -f reference.fasta -o alignment.html -n 200] region1 [region2 region3 ...]

.. program:: breseq_bam2aln

Display reads aligned to the specified region or regions.

Commonly used options:

.. option:: -b <file_path>, --bam <file_path> 

   BAM database file of read alignments (DEFAULT=data/reference.bam).

.. option:: -f <file_path>, --fasta <file_path> 

   FASTA file of reference sequences (DEFAULT=data/reference.fasta).

.. option:: -o <path>, --output <path> 

   Output path. If there is just one region, the name of the output file (DEFAULT=region1.*). If there are multiple regions, this argument must be a directory path, and all output files will be output here with names region1.*, region2.*, ... (DEFAULT=.).

.. option:: -r <region> , --region <region>

   Regions to create alignments for. Must be provided as sequence regions in the format **ACCESSION:START-END**, where **ACCESSION** is a valid identifier for one of the sequences in the FASTA file, and **START** and **END** are 1-indexed coordinates of the beginning and end positions. Any read overlapping these positions will be shown. A separate output file is created for each region. Regions may be provided at the end of the command line as unnamed arguments.

.. option:: -n <int>, --max-reads <int>

   Maximum number of reads that will be aligned to a region. If there are more than this many reads, then the reads displayed are randomly chosen and a warning is added to the output. (DEFAULT=200).


Command: bam2cov
--------------------------

Usage::

  breseq BAM2COV [-b reference.bam -f reference.fasta --format PNG -o output.png] region1 [region2 region3 ...]

.. program:: breseq_bam2cov

Create a coverage plot or table for the specified region or regions.

Commonly used options:

.. option:: -b <file_path>, --bam <file_path> 

   BAM database file of read alignments (DEFAULT=data/reference.bam).

.. option:: -f <file_path>, --fasta <file_path> 

   FASTA file of reference sequences (DEFAULT=data/reference.fasta).
   
.. option:: -o <path>, --output <path> 

   Output path. If there is just one region, the name of the output file (DEFAULT=region1.*). If there are multiple regions, this argument must be a directory path, and all output files will be output here with names region1.*, region2.*, ... (DEFAULT=.).

.. option:: -r <region>, --region <region>

   Regions to create alignments for. Must be provided as sequence regions in the format **ACCESSION:START-END**, where **ACCESSION** is a valid identifier for one of the sequences in the FASTA file, and **START** and **END** are 1-indexed coordinates of the beginning and end positions. Any read overlapping these positions will be shown. A separate output file is created for each region. Regions may be provided at the end of the command line as unnamed arguments.

.. option:: --format <PNG/PDF> 

   Format of output plot: PNG or PDF. (DEFAULT=PNG).
   
.. option:: -t, --table

   Create tab-delimited file of coverage instead of a plot.

.. option:: -1, --total-only

   Only plot/tabulate the total coverage at a position. That is, do not not output the coverage on each genomic strand.
   
.. option:: --resolution <int>

  Number of positions to output coverage information for in interval (0=ALL) (DEFAULT=600).
