#!/usr/bin/perl -w

###
# Pod Documentation
###

=head1 NAME

fasta_quals_to_fastq.pl

=head1 SYNOPSIS

Usage: fasta_quals_to_fastq.pl input.fna input.qual > new.fastq

Convert a FASTA file and a quality score file into a FASTQ file. 
Output is to STDOUT. Assumes the quality score file has the same sequence
order and names as the FASTA file, and following each entry a list of integers
that give PHRED style quality scores for each base.

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
my ($help, $man);
my ($input_file, $start, $end);
my $strip_names = 0;
GetOptions(
	'help|?' => \$help, 'man' => \$man,
	'strip-names|s' => \$strip_names,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(-exitstatus => 0, -verbose => 2) if (scalar @ARGV != 2);
my ( $input_1, $input_2 ) = @ARGV;

# Solexa quality scores are -10 log10(P/(1-P))
# $Q = ord($q) - 64;

# PHRED (454) quality scores are -10 log10(P) 
# $Q = ord($q) - 33;

#Open Input/output

open INPUT_1, "<$input_1";
open INPUT_2, "<$input_2";

#my $in_1 = Bio::SeqIO->new( -file => "<$input_1", -format => "fastq");
#my $in_2 = Bio::SeqIO->new( -file => "<$input_2", -format => "fastq");
#my $out  = Bio::SeqIO->new( -file => ">$output", -format => "fastq");

#my $seq_1 = $in_1->next_seq;
#my $seq_2 = $in_2->next_seq;

my $line_1 = <INPUT_1>;
my $line_2 = <INPUT_2>;

my $num_seq = 0;
while ($line_1)
{
	#
	$num_seq++;
	chomp $line_1;
	my $read_name = $line_1;
	$read_name =~ s/^>/@/;
	if ($strip_names)
	{
		$read_name = "\@$num_seq\n";
	}

	##collect sequence
	my $sequence;
	
	while ($line_1 = <INPUT_1>)
	{
		last if ($line_1 =~ m/^>/);
		chomp $line_1;
		$sequence .= $line_1;
	}
	
	##collect qualities
	my @quality_list;
	
	while ($line_2 = <INPUT_2>)
	{
		last if ($line_2 =~ m/^>/);
		chomp $line_2;
		push @quality_list, split /\s+/, $line_2;
	}
	
	#output...
	print "$read_name\n$sequence\n+\n";
#	print +(phred_quality_to_char(0));
	foreach my $quality (@quality_list)
	{
		print +(phred_quality_to_char($quality));
	}
	print "\n";
	
	
	##print in fastq format
	
	#$seq_1 = $in_1->next_seq;
	#$seq_2 = $in_2->next_seq;
	
	#$out->write_seq($seq_1);
	#$out->write_seq($seq_2);
}

#my $index_m1 = $index-1;
#system "cat $input_1\_$index_m1 $input_1\_$index > temp; mv temp > $input_1\_$index_m1";
#system "cat $input_2\_$index_m1 $input_2\_$index > temp; mv temp > $input_2\_$index_m1";



sub phred_to_solexa_quality
{
	my ($pq) = @_;
	my $e = exp(-$pq * log(10) / 10);
		
	my $sq = -10 * log ($e / (1-$e)) / log (10);
	return int floor $sq;
}

sub phred_quality_to_char
{
	my ($q) = @_;
	return chr($q+33);
}