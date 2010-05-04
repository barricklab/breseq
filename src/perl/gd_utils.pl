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
		
	my $removed_gd = $gd->filter_mutations(\&mutation_filter);
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

	#load all genome diff files
	my (@mutations, @deletions, @unknowns, @hybrids);
	
	foreach my $gd_file (@{$settings->{input_genome_diffs}})
	{
		my $mutation_info = Breseq::Output::read_genome_diff($gd_file);

		push @mutations, @{$mutation_info->{mutations}};
		push @deletions, @{$mutation_info->{deletions}};
		push @unknowns, @{$mutation_info->{unknowns}};
		push @hybrids, @{$mutation_info->{hybrids}};
	}	
	
	###
	## Gather together insertions and deletions that occur next to each other
	## ...unless we are predicting polymorphisms
	###

	if (!$settings->{polymorphism_prediction})
	{
		my $lc;
		for (my $i=0; $i<scalar @mutations; $i++)
		{
			my $c = $mutations[$i];	 

			#if the same position and both are tandem insertion
			if ($lc && ($lc->{start} == $c->{start}) && ($lc->{insert_end}+1 == $c->{insert_start}))
			{
				$lc->{insert_end} = $c->{insert_start};
				$lc->{new_seq} .= $c->{new_seq};
				$lc->{quality} .= ",$c->{quality}";
				$lc->{total_coverage_string} .= ",$c->{total_coverage_string}";
				$lc->{best_coverage_string} .= ",$c->{best_coverage_string}";
				splice @mutations, $i, 1;
				$i--;
				next;
			}

			#if the positions are next to each other and both are deletions...
			if ($lc && ($lc->{end}+1 == $c->{start}) && ($lc->{new_seq} eq '.') && ($c->{new_seq} eq '.') )
			{
				$lc->{end} = $c->{start};
				$lc->{ref_seq} .= $c->{ref_seq};
				$lc->{quality} .= ",$c->{quality}";
				$lc->{total_coverage_string} .= ",$c->{total_coverage_string}";
				$lc->{best_coverage_string} .= ",$c->{best_coverage_string}";
				splice @mutations, $i, 1;
				$i--;
				next;
			}

			$lc = $c;
		}
	}
	
	##
	# Annotate mutations and deletions
	sub mutation_annotation {}	
	##
	print STDERR "Annotating within-read mutations...\n";
	Breseq::ReferenceSequence::annotate_mutations(undef, undef, $ref_seq_info, \@mutations);
	#print Dumper(\@mutations); ##DEBUG

	print STDERR "Annotating deletions...\n";
	Breseq::ReferenceSequence::annotate_deletions(undef, undef, $ref_seq_info, \@deletions);
	#print Dumper(\@deletions); ##DEUG		

	##
	# Write text output files
	sub text_output {}
	##
	print STDERR "Creating text output files...\n";
#	Breseq::Output::save_text_mutation_file($settings->file_name('mutations_text_file_name') , \@mutations);
#	Breseq::Output::save_text_deletion_file($settings->file_name('deletions_text_file_name') , \@deletions);
#	Breseq::Output::save_text_unknown_file($settings->file_name('unknowns_text_file_name') , \@unknowns);
	
	##
	# Plot coverage of genome and large deletions
	sub plot_coverage {}
	##
	print STDERR "Drawing coverage graphs...\n";
	$settings->create_path('coverage_graph_path');
	my $coverage_graph_path = $settings->file_name('coverage_graph_path');	
	my $deletions_text_file_name = $settings->file_name('deletions_text_file_name');
	Breseq::Output::save_text_deletion_file($deletions_text_file_name, \@deletions);
	
	foreach my $seq_id (@{$ref_seq_info->{seq_ids}})
	{
		my $this_plot_coverage_done_file_name = $settings->file_name('plot_coverage_done_file_name', {'@'=>$seq_id});

		if (!-e $this_plot_coverage_done_file_name)
		{
			my $this_complete_coverage_text_file_name = $settings->file_name('complete_coverage_text_file_name', {'@'=>$seq_id});			
			my $res = Breseq::Shared::system("graph_coverage.pl -t $coverage_graph_path -p $settings->{coverage_graph_path} -i $deletions_text_file_name -c $this_complete_coverage_text_file_name --seq_id=$seq_id");				
			die if ($res);
			open DONE, ">$this_plot_coverage_done_file_name";
			close DONE;
		}
		else
		{
			print STDERR "Drawing coverage graphs already complete.\n";
		}

		#need to assign link names even if coverage was already drawn
		my $i=1;
		my @this_deletions = grep {$_->{seq_id} eq $seq_id} @deletions if ($seq_id);
		foreach my $del (@this_deletions)
		{
			$del->{coverage_graph_link} = "$settings->{local_coverage_graph_path}/$seq_id\.$i\.pdf";
			$i++;
		}
	}
#	$settings->remove_path('deletions_text_file_name');

	sub html_output {}

	##
	# Output SNPs
	##

	### make alignments first because we fill in a link field used by the summary file
	my @composite_list; ##

	### screen out polymorphism predictions at this step
	if ($settings->{polymorphism_prediction})
	{
		foreach my $c (@mutations)
		{
			my $accept = 1;

			if ($c->{polymorphism})
			{
				my $polymorphism = $c->{polymorphism};
				#print Dumper($polymorphism);

				$accept = 0 if ($polymorphism->{log10_e_value} < $settings->{polymorphism_log10_e_value_cutoff});
				$accept = 0 if ($polymorphism->{fisher_strand_p_value} < $settings->{polymorphism_fisher_strand_p_value_cutoff});
				$accept = 0 if ($polymorphism->{fraction} < $settings->{polymorphism_fraction_cutoff});
				$accept = 0 if ($polymorphism->{fraction} > 1-$settings->{polymorphism_fraction_cutoff});	
			}	

			print "ACCEPT: $accept\n";

			push @composite_list, $c if ($accept);
		}
		@mutations = @composite_list;
	}
	else
	{
		push @composite_list, @mutations;
	}

	# 
	### look, we have to invent intervals for the deletions so each has an upstream and downstream.
	foreach my $deletion (@deletions)
	{
	 	#ok, this has a circular reference, which may possibly be a very bad thing.
	 	$deletion->{upstream_interval} = { start => $deletion->{start}-1, end => $deletion->{start}-1, 
	 		deletion => $deletion, seq_id => $deletion->{seq_id} };

	 	$deletion->{downstream_interval} = { start => $deletion->{end}+1, end => $deletion->{end}+1, 
 			deletion => $deletion, seq_id => $deletion->{seq_id} };
		push @composite_list, $deletion->{upstream_interval};
		push @composite_list, $deletion->{downstream_interval};
	}	

	### look, we have to invent intervals for the non-rearranged versions of the hybrid reads as well.
	### proper alignments can be made as part of the composite list

	###
	# Sort the junctions by unique coordinates or by their scores
	###
	
	if ($settings->{sort_junctions_by_score})
	{
		#sort and truncate list
		@hybrids = sort { -($a->{score} <=> $b->{score}) || ($a->{total_reads} <=> $a->{total_reads}) } @hybrids;
		my $last = $settings->{max_junctions_to_print};
		$last = scalar @hybrids if (scalar @hybrids < $last);
		$#hybrids = $last-1;	
	}

	print STDERR "Annotating rearrangements...\n";
	Breseq::ReferenceSequence::annotate_rearrangements($settings, undef, $ref_seq_info, \@hybrids);

	foreach my $hybrid (@hybrids)
	{
		#print STDERR Dumper($hybrid);	
		push @composite_list, $hybrid;	
	 	push @composite_list, $hybrid->{interval_1};
	 	push @composite_list, $hybrid->{interval_2};
	}
	
	

	# ### first name all the files/links, so backlinks work

	my $reference_bam_file_name = $settings->file_name('reference_bam_file_name');
	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');

	foreach my $c (@composite_list)
	{
	 	my $html_alignment_file_name = "$c->{seq_id}_$c->{start}_$c->{end}_alignment.html";
	 	$c->{link} = "$settings->{local_alignment_path}/$html_alignment_file_name";
	 	$c->{file_name} = "$settings->{alignment_path}/$html_alignment_file_name";
		$c->{bam_path} = $reference_bam_file_name;
		$c->{fasta_path} = $reference_fasta_file_name;
#		print Dumper($c);
	}

	## hybrids use different BAM files for making the alignments!!!
	my $junction_bam_file_name = $settings->file_name('junction_bam_file_name');
	my $junction_fasta_file_name = $settings->file_name('candidate_junction_fasta_file_name');

	foreach my $c (@hybrids)
	{	
		$c->{bam_path} = $junction_bam_file_name;
		$c->{fasta_path} = $junction_fasta_file_name;
		
		## rename junctions so they don't clobber normal alignments
		my $html_alignment_file_name = "JCT_$c->{seq_id}_$c->{start}_$c->{end}_alignment.html";
	 	$c->{link} = "$settings->{local_alignment_path}/$html_alignment_file_name";
	 	$c->{file_name} = "$settings->{alignment_path}/$html_alignment_file_name";
	
		foreach my $int ('interval_1', 'interval_2')
		{
			$html_alignment_file_name = "JCT_$c->{$int}->{seq_id}_$c->{$int}->{start}_$c->{$int}->{end}_alignment.html";
			$c->{$int}->{link} = "$settings->{local_alignment_path}/$html_alignment_file_name";
		 	$c->{$int}->{file_name} = "$settings->{alignment_path}/$html_alignment_file_name";
		}
	}

	### now create alignment files
	if (!$settings->{no_alignment_generation})
	{
		$settings->create_path('alignment_path');

		print STDERR "Creating alignment HTML files...\n";
		foreach my $c (@composite_list) # , @hybrids)
		{			
			print STDERR "Creating alignment file: $c->{link}\n";
			Breseq::Output::html_alignment_file($settings, $c);		
		}
	}


	###
	## HTML output
	###	

	print STDERR "Creating full HTML table...\n";	
	my $mutation_file_name = $settings->file_name('mutations_html_file_name');
	$mutation_file_name = $settings->{html_mutation_file} if ($settings->{html_mutation_file});
	Breseq::Output::html_full_table($mutation_file_name, $settings, $ref_seq_info, \@mutations, \@deletions, \@hybrids);
	
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
		my $gd1 = GenomeDiff->new(-FILE_NAME => $gd1_file);
		push @gds1, $gd1;
	}

	my @gds2;
	foreach my $gd2_file (@input2_gd_files)
	{
		my $gd2 = GenomeDiff->new(-FILE_NAME => $gd2_file);
		push @gds2, $gd2;
	}
	
	die if ((!@gds1) || (!@gds2));
	
	my $output_gd = GenomeDiff::subtract(\@gds1, \@gds2);
	$output_gd->write($output);
}

