#!/usr/bin/env perl -w

###
# Pod Documentation
###

=head1 NAME

sam2html.pl

=head1 SYNOPSIS

Usage: sam2html.pl [-o output] database.bam reference.fasta REL606:1-100

Create an HTML alignment of a region. Region must not be too big!

=head1 DESCRIPTION

=over

=back

=head1 AUTHOR

Jeffrey Barrick
<barrick@msu.edu>

=head1 COPYRIGHT

Copyright 2009.  All rights reserved.

=cut

###
# End Pod Documentation
###


use strict;
use Data::Dumper;

use FindBin;
use lib $FindBin::Bin;

use CGI qw/:standard/;

use Breseq::AlignmentOutput;

#Get options
use Getopt::Long;
use Pod::Usage;
my ($help, $man);
my ($input_file, $start, $end);
my $output_file = 'output.html';

GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'output|o=s' => \$output_file,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my ($bam_path, $fasta_path, @regions) = @ARGV;

die "No regions defined.\n" if (scalar @regions == 0);
die "No BAM file defined.\n" if (!$bam_path);
die "No fasta file defined.\n" if (!$fasta_path);

open OUT, ">$output_file";
print OUT start_html(-title => "bam2html");

foreach my $region (@regions)
{
	my $ao =Breseq::AlignmentOutput->new;
	print OUT $ao->html_alignment($bam_path, $fasta_path, $region);
}
print OUT end_html;
close OUT;
