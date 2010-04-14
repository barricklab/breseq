#!/usr/bin/env perl -w

###
# Pod Documentation
###

=head1 NAME

fasta_quals_to_fastq.pl

=head1 SYNOPSIS

Usage: fasta_list_to_fastq.pl input.list_fastq output.fastq

Convert a FASTQ file where the quality line is a list of numbers separated by spaces
into a fastq file where the qualities are given as a character string.
and a quality score file into a FASTQ file. 

=head1 DESCRIPTION

=over

=back

=head1 AUTHOR

Jeffrey Barrick
<barrick@msu.edu>

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


#print "$input_1 $input_2 $output\n";

#Get options
use Getopt::Long;
use Pod::Usage;
use FastqLite;
my ($help, $man);
my $quality_score_format = 'phred';
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'quality-score-format|q=s' => \$quality_score_format,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(-exitstatus => 0, -verbose => 2) if (scalar @ARGV != 2);
my ( $input_file, $output_file) = @ARGV;


print "Using quality score format: $quality_score_format\n";
my $in = FastqLite->new( -file => "<$input_file", -list_format => 1, -quality_format => $quality_score_format);
my $output_fh;
open $output_fh, ">$output_file" or die "Could not open $output_file for output\n";
while (my $read_seq = $in->next_seq) 
{
	FastqLite::write_seq($output_fh, $read_seq);
}
