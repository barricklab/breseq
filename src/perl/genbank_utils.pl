#!/usr/bin/env perl

###
# Pod Documentation
###

=head1 NAME

genbank_utils.pl

=head1 SYNOPSIS

Usage: genbank_utils.pl SPLIT -i input.gbk -o output.gbk

Perform various operations on GenBank files that preserve 

=head1 DESCRIPTION

=over

=head1 COMMAND: SPLIT

=over

=item B<-i> <file path> 

Input GanBank file. Output files will be named in the form <input>_#.

=item B<-n> <integer>

Number of equal files to split GenBank into.

=back

=back

=head1 AUTHOR

Jeffrey Barrick
<jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2010.  All rights reserved.

=cut

##
## Note that really we should use the average unique coverage as a poisson mean
## and fit this form of generalized linear model to the observation in each file
## simultaneously to get a per-position offset to the counts. Then we can calculate
## a true expected 95% confidence interval. 
##

###
# End Pod Documentation
###
use strict;

use File::Path;
use FindBin;
use lib $FindBin::Bin;
use Data::Dumper;
use POSIX;

#Get options
use Getopt::Long;
use Pod::Usage;

my $command = shift @ARGV;

if ($command eq 'SPLIT')
{
	do_split();
}
elsif ("\U$command" eq 'MUTATE')
{
	do_mutate();
}
elsif ("\U$command" eq 'FASTA')
{
	do_convert_to_fasta();
}

sub do_split
{
	use Bio::SeqIO;
	my ($help, $man);
	my $verbose;
	my $input;
	my $n = 2;

	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'number|n=s' => \$n,
	);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;

	my $input = shift @ARGV;
	
	my $ref_i = Bio::SeqIO->new( -file => "<$input");
	my $ref_seq = $ref_i->next_seq;
	my @features = $ref_seq->get_SeqFeatures();

	for (my $i=0; $i<$n; $i++) 
	{
		my $start = floor($i* $ref_seq->length / $n) + 1;
		my $end = floor(($i+1)* $ref_seq->length / $n);

		print "Fragment $i: start: $start end: $end\n";
	
		#this strips all features...
		my $fragment = $ref_seq->trunc($start, $end);
		$fragment->display_id($fragment->display_id . "-$i");
		$fragment->accession($fragment->accession . "-$i");
	
		my $fragment_loc = Bio::Location::Simple->new(-start => $start, -end => $end, -strand => 0 );
		foreach my $f (@features)
		{
		
			if ($fragment_loc->overlaps($f->location))
			{
				#print $f->start . " " . $f->end . "\n";
			
				if ($f->location->start < $fragment_loc->start)
				{
					$f->location->start_pos_type('BEFORE');
					$f->location->start($fragment_loc->start);
				}
			
				if ($f->location->end > $fragment_loc->end)
				{
					$f->location->end_pos_type('AFTER');
					$f->location->end($fragment_loc->end);
				}			
			
				$f->location->start( $f->location->start + 1 - $fragment_loc->start );
				$f->location->end( $f->location->end + 1 - $fragment_loc->start );
						
				$fragment->add_SeqFeature($f);
			}
		}
	
	
		my $output = $input;
		$output =~ s/\.(.+?)$/-$i.$1/;
	
		my $ref_o = Bio::SeqIO->new( -file => ">$output", -format => "GenBank");
		$ref_o->write_seq($fragment);
	}
}

## Changes fasta file
sub do_mutate
{
	use Bio::SeqIO;
	
	my ($help, $man);
	my $verbose;
	my $output = "output.fna";

	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'output|o=s' => \$output,
	);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;

	my $input = shift @ARGV;
	my $change = shift @ARGV;

	my $ref_i = Bio::SeqIO->new( -file => "<$input");
	my $ref_seq = $ref_i->next_seq;
	my $seq = $ref_seq->seq;
	
	my $ref_o = Bio::SeqIO->new( -file => ">$output");
	
	my ($pos, $len, $new) = split "-", $change;
	
	print "Position = $pos\n";
	print "Length   = $len\n";
	print "New      = $new\n";

	$new =~ s/\.//g;
	substr($seq, $pos-1, $len) = $new;
	
	$ref_seq->seq($seq);
	$ref_o->write_seq($ref_seq);	
}

## Changes fasta file
sub do_convert_to_fasta
{
	use Bio::SeqIO;
	
	my ($help, $man);
	my $verbose;
	my $output = "output.fna";

	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'output|o=s' => \$output,
	);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;

	my $input = shift @ARGV;

	my $ref_i = Bio::SeqIO->new( -file => "<$input");
	my $ref_seq = $ref_i->next_seq;
	
	my $ref_o = Bio::SeqIO->new( -file => ">$output");
	$ref_o->write_seq($ref_seq);	
}