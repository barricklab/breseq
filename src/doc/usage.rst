Usage
==============

:program:`breseq`
------------------

Usage:

   :program:`breseq` -r reference.gbk reads1.fastq [reads2.fastq, reads3.fastq...]

.. program:: breseq

Runs the :program:`breseq` pipeline for mutation prediction from short-read genome re-sequencing data.

Options:

.. option:: -r <file_path>, --reference <file_path> 

   Input reference genome sequence files in GenBank format. If there are multiple reference sequences stored in separate GenBank files (i.e., a bacterial genome and a plasmid),  this option can be supplied multiple times.

.. option:: reads1.fastq [reads2.fastq, reads3.fastq...]  

   The remaining arguments at the command line are the FASTQ input files of reads. Any FASTQ quality score style (e.g., Sanger, Solexa, Illumina 1.5+) is accepted, :program:`breseq` internally re-calibrates the error rates. It does this for each FASTQ file separately, so data sets that were generated independently should be stored in separate input files.

Other tools
------------------

:program:`bam2aln`
*********************

Usage:

   :program:`bam2aln` [-b input.bam] [-f input.fasta] [-o output/path] region1 [region2 region3 ...]

.. program:: bam2aln

Creates HTML pileup files displaying reads aligned to each specified region.

Options:

.. option:: -b <file_path>, --bam=<file_path> 

   BAM database file of read alignments. Defaults: ``reference.bam``, ``data/reference.bam``.

.. option:: -f <file_path>, --fasta=<file_path> 

   FASTA file of reference sequences. Defaults: ``reference.fasta``, ``data/reference.fasta``.
   
.. option:: -o <path>, --output=<path> 

   Output path. If there are multiple regions, must be a directory path, and all output files will be output with names region1.html, region2.html, ... If there is one region, the output file will be given this name if it is not the name of an already existing directory. Default: current path.
   
.. option:: -n <int>, --max-reads=<int>

   Maximum number of reads that will be aligned to a region. If there are more than this many reads then the reads displayed are randomly chosen displayed and a warning is added to the output. Default: 1000.

.. option:: region1 [region2, region3, ...]

   Regions to create output for must be provided in the format ``FRAGMENT:START-END``, where ``FRAGMENT`` is a valid identifier for one of the sequences in the FASTA file, and ``START`` and ``END`` are 1-indexed coordinates of the beginning and end of the alignment. A separate output file is created for each region.


:program:`bam2cov`
******************

Usage:

   :program:`bam2cov` -b input.bam -f input.fasta -o [output/path] region1 [region2, region3, ...]


.. program:: bam2cov

Creates a coverage table or image for the specified region.
   
Options:

.. option:: -b <file_path>, --fasta <file_path> 

   BAM database file of read alignments. Defaults: ``reference.bam``, ``data/reference.bam``

.. option:: -f <file_path>, --fasta <file_path> 

   FASTA file of reference sequences. Defaults: ``reference.fasta``, ``data/reference.fasta``
   
.. option:: -o <path>, --output <path> 

   Output path. If provided, all output files are created in this location. Default: current path.

.. option:: region1 [region2, region3, ...]

   Regions to create output for must be provided in the format ``FRAGMENT:START-END``, where ``FRAGMENT`` is a valid identifier for one of the sequences in the FASTA file, and ``START`` and ``END`` are 1-indexed coordinates of the beginning and end of the alignment. A separate output file is created for each region.