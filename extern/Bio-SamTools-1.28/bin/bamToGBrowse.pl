#!/usr/bin/perl

use strict;
use Bio::DB::Sam;
use Bio::DB::Sam::SamToGBrowse;

my $has_bigwig = eval {require Bio::DB::BigFile;1};

@ARGV >= 1 or die <<USAGE;
Usage: $0 <directory_path> [<ref.fa>]

Given the path to a directory and the fasta file for the reference
sequence, do the following:

 1. Find all SAM files in the indicated directory and convert them
    into BAM. These must end in one of the extensions ".sam" or ".sam.gz".
    A series of <base>.bam files will be created. This step will be
    skipped if SAM files are absent or the corresponding BAM files
    are already present.

 2. Sort the newly created BAM files, creating a series of files named
    <base>_sorted.bam.

 3. Index BAM files that need indexing. This step will look for
      files named <base>_sorted.bam

 4. Create a set of BigWig files representing the coverage graph. These
      will be named <base>.bw.

 5. Create a skeletal GBrowse config file named "gbrowse.conf" that
    serves as a starting point for viewing these files. Previous versions
    of this file will be appended to.

You can prepopulate the directory with sorted and indexed BAM files,
in which case the script will not attempt to resort or reindex them.
Unless the Fasta file is explicitly provided, this script will look in
the designated directory for ONE .fa file to use.

Note that you will need temporary space in the directory equivalent to
the size of the largest SAM file being processed.

This script takes a long time to run and uses significant amounts of
RAM when generating the coverage graphs. To improve both performance
and memory consumption, you can install the following C binaries:

  bedGraphToBigWig  -- download from http://hgdownload.cse.ucsc.edu/admin/exe
                       or build from source downloadable from
                       http;//hgdownload.cse.ucsc.edu/admin/jksrc.zip
  genomeCoverageBed -- download from http://code.google.com/p/bedtools/

Place these executables into your path, e.g. /usr/local/bin or ~/bin.
USAGE
    ;

my($dir,$fasta) = @ARGV;
my $converter = Bio::DB::Sam::SamToGBrowse->new($dir,$fasta,'verbose');
$converter->run();

exit 0;

