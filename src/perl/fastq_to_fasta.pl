#!/usr/bin/env perl -w

###
# Pod Documentation
###

=head1 NAME

fastq_to_fasta.pl

=head1 SYNOPSIS

Usage: fastq_to_fasta.pl input.fastq output.fasta

Convert a FASTQ formatted file to FASTA. Important: This script only 
works with FASTQ files where the sequence and qualities occupy only a 
single line each.

=head1 DESCRIPTION

=over

=back

=head1 AUTHOR

Jeffrey Barrick
<jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2008.  All rights reserved.

=cut

###
# End Pod Documentation
###


use strict;
use Data::Dumper;

use FindBin;
use lib $FindBin::Bin;

use Bio::Location::Simple;
use Bio::DB::Flat;
use Bio::SeqIO;

use FastqLite;

#Get options
use Getopt::Long;
use Pod::Usage;
my ($help, $man);
my $options;
my ($strip_names, $paired, $number_per_file);
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'strip|s' => \$options->{strip_names},
	'paired|p' => \$options->{paired},
	'number-per-file|n=s' => \$options->{number_per_file},
	
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my ( $input, $output ) = @ARGV;

my $reads = FastqLite::split_fastq_to_fasta($input, $output, $options);
print "$reads reads converted.\n"
