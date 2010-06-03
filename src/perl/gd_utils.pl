#!/usr/bin/env perl

###
# Pod Documentation
###

=head1 NAME

gd_utils.pl

=head1 SYNOPSIS

Usage: gd_utils.pl COMMAND [options]

=head1 DESCRIPTION

Perform various operations on "genome diff" files.

=over

=head1 COMMAND: FILTER

=over

=item B<-r|--reference> file path

Genbank file or files that will be used for annotating mutations.

=back

=item B<-i|--input-path> directory path

File containing columns that can be translated into the required information about SNP mutations. May
provide multiple times to create multiple output files.

=item B<-g|--genome-info> file path

Names or correspondance table. Connects filenames to the generation, line, and clone designation of a genome. 

=item B<--output-file> genbank sequence id

Tab-delimited output file. DEFAULT = 'output.tab'.

=item B<--ignore-mutation> seq_id:position:ref_base>mut_base

Ignore this SNP mutation when counting mutations. Any of '/', ':', '_', and '>' can be used to separate the items.

=item B<--seq-id> genbank sequence id

Force all mutations to be on this sequence id, regardless of other information. OPTIONAL.

=back

=head1 AUTHOR

Jeffrey Barrick <jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2009.  All rights reserved.

=cut

###
# End Pod Documentation
###

our $VERSION = '0';

#### Standard Perl Modules ####
use strict;
use Data::Dumper;
use File::Path;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil floor);

#### Bio Perl ####
use Bio::SeqIO;

#allows us to load modules and call other scripts in same directory
use FindBin;
use lib $FindBin::Bin;
$ENV{PATH} = "$ENV{PATH}:" . $FindBin::Bin;

use Breseq::GenomeDiff;

##First argument is always the command
my $command = shift @ARGV;

#Could auto-complete abbreviated commands here...

if ("\U$command" eq 'FILTER')
{
	do_filter();
}
elsif ("\U$command" eq 'ANNOTATE')
{
	do_annotate();
}
elsif ("\U$command" eq 'UNION')
{
	do_union();
}
elsif ("\U$command" eq 'INTERSECTION')
{
	do_intersecton();
}
elsif ("\U$command" eq 'SUBTRACT')
{
	do_subtract();
}
else
{
	die "Unknown command";
}


sub do_filter
{	
	my ($help, $man, $verbose);
	my $output = 'output.gd'; 
	my $removed = 'removed.gd';
	my $input = ();
	
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'verbose|v' => \$verbose,
		'output|o=s' => \$output,
		'input|i=s' => \$input,
		'removed|r=s' => \$removed,
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;

	my @input_conditions = @ARGV;

	our @conditions = ();
	foreach my $ic (@input_conditions)
	{		
		my $c;
	 	if (!($ic =~ m/^(.+?)\s*?(==|<=|>=|<|>|!=|\s+eq\s+|\s+ne\s+|\s+lt\s+|\s+gt\s+|\s+le\s+|\s+ge\s+)\s*?(.+)$/))
		{
			print STDERR "Could not interpret filter: $ic\n";
			next;
		}
		
		$c->{key} = $1;
		$c->{value} = $3;
		$c->{comparison} = $2;
		push @conditions, $c;
		
		print STDERR "Condition: $c->{key} | $c->{comparison} | $c->{value}\n";
		
	}
	
	my $gd = Breseq::GenomeDiff->new(-FILE_NAME => $input);
		
	### screen out polymorphism predictions at this step
	sub mutation_filter
	{
		my ($mut) = @_;
		
		
		my $accept = 1;
		
		foreach my $c (@conditions)
		{						
			next if (!defined $mut->{$c->{key}});
			
			#print STDERR "return (\"$mut->{$c->{key}}\" $c->{comparison} \"$c->{value}\")\n";
			
			my $res = 1;
			if ($c->{comparison} =~ m/(>)/)
			{
				$res = eval "return ($mut->{$c->{key}} $c->{comparison} $c->{value})";
			}
			else
			{
				$res = eval "return (\"$mut->{$c->{key}}\" $c->{comparison} \"$c->{value}\")";
			}
		 	$accept = ($accept && $res);
		}
		
		return $accept;
	}	
		
	my $removed_gd = $gd->filter(\&mutation_filter);
	$removed_gd->write($removed) if ($removed);	
	$gd->write($output);
}

sub do_annotate
{	
	use Breseq::Settings;
	use Breseq::ReferenceSequence;
	use Breseq::Output;
	
	my $settings = Breseq::Settings->new_annotate;
	
	## load information about reference sequences from GenBank files
	my $ref_seq_info = Breseq::ReferenceSequence::load_ref_seq_info($settings);
	
	my $gd = Breseq::GenomeDiff->new();
	foreach my $gd_file (@{$settings->{input_genome_diffs}})
	{
		my $mutation_info = $gd->read($gd_file);		
	}	
		
	##
	# Annotate mutations
	##
	print STDERR "Annotating mutations...\n";
	Breseq::ReferenceSequence::annotate_mutations($ref_seq_info, $gd);
	#print Dumper(\@mutations); ##DEBUG	

	##
	# Plot coverage of genome and large deletions
	##
	print STDERR "Drawing coverage graphs...\n";
	Breseq::Output::draw_coverage($settings, $ref_seq_info, $gd);
	
	##
	# Creare evidence files containing alignments and coverage plots
	##
	if (!$settings->{no_alignment_generation})
	{
		Breseq::Output::create_evidence_files($settings, $gd);
	}

		###
		# Sort the junctions by unique coordinates or by their scores
		###

=comment	
		if ($settings->{sort_junctions_by_score})
		{
			#don't show junctions supported by only a few reads
			@hybrids = grep {$_->{total_non_overlap_reads} > 2} @hybrids;

			#sort and truncate list
			@hybrids = sort { -($a->{score} <=> $b->{score}) || ($a->{total_reads} <=> $a->{total_reads}) } @hybrids;
			my $last = $settings->{max_junctions_to_print};
			$last = scalar @hybrids if (scalar @hybrids < $last);
			$#hybrids = $last-1;	
		}
=cut

	###
	## HTML output
	###	

	print STDERR "Creating index HTML table...\n";	

	my $summary = {};
	my $sequence_conversion_summary_file_name = $settings->file_name('sequence_conversion_summary_file_name');	
	$summary->{sequence_conversion} = Storable::retrieve($sequence_conversion_summary_file_name);

	my $index_html_file_name = $settings->file_name('index_html_file_name');	
	Breseq::Output::html_index($index_html_file_name, $settings, $summary, $ref_seq_info, $gd);	
}

sub do_union
{
	my ($help, $man, $verbose);
	my $output = 'output.gd'; 
	my @input_gd_files = ();
	
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'verbose|v' => \$verbose,
		'output|o=s' => \$output,
		'input|i=s' => \@input_gd_files,
			
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;

	my @gds;
	foreach my $gd_file (@input_gd_files)
	{
		my $gd = GenomeDiff->new(-FILE_NAME => $gd_file);
		push @gds, $gd;
	}

	my $output_gd = GenomeDiff::union(\@gds);
	$output_gd->write($output);
}

sub do_intersection
{	
	my ($help, $man, $verbose);
	my $output = 'output.gd'; 
	my @input_gd_files = ();
	
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'verbose|v' => \$verbose,
		'output|o=s' => \$output,
		'input|i=s' => \@input_gd_files,
			
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;

	my @gds;
	foreach my $gd_file (@input_gd_files)
	{
		my $gd = GenomeDiff->new(-FILE_NAME => $gd_file);
		push @gds, $gd;
	}

	my $output_gd = GenomeDiff::intersection(\@gds);
	$output_gd->write($output);
}

sub do_subtract
{
	my ($help, $man, $verbose);
	my @input1_gd_files = ();
	my @input2_gd_files = ();
	my $output = 'output.gd'; 
	GetOptions(
		'help|?' => \$help, 'man' => \$man,
		'verbose|v' => \$verbose,
		'output|o=s' => \$output,
		'input1|1=s' => \@input1_gd_files,
		'input2|2=s' => \@input2_gd_files,
	) or pod2usage(2);
	pod2usage(1) if $help;
	pod2usage(-exitstatus => 0, -verbose => 2) if $man;
	
	my @gds1;
	foreach my $gd1_file (@input1_gd_files)
	{
		my $gd1 = Breseq::GenomeDiff->new(-FILE_NAME => $gd1_file);
		push @gds1, $gd1;
	}
	
	my @gds2;
	foreach my $gd2_file (@input2_gd_files)
	{
		my $gd2 = Breseq::GenomeDiff->new(-FILE_NAME => $gd2_file);
		push @gds2, $gd2;
	}	
	
	die if ((!@gds1) || (!@gds2));
	
	my $output_gd = Breseq::GenomeDiff::subtract(\@gds1, \@gds2);
	$output_gd->write($output);
}

