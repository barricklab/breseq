###
# Pod Documentation
###

=head1 NAME

Breseq::CoverageDistribution

=head1 SYNOPSIS

Perl modules used internally by breseq.

=head1 AUTHOR

Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010 Michigan State University

breseq is free software; you can redistribute it and/or modify it under the terms the 
GNU General Public License as published by the Free Software Foundation; either 
version 1, or (at your option) any later version.

=cut

###
# End Pod Documentation
###

use strict;

package Breseq::CoverageDistribution;
use vars qw(@ISA);

use Breseq::Shared;
use Data::Dumper;
use FindBin;

sub new
{
	my($caller, $args) = @_;
	
	my $class = ref($caller) || $caller;
	my $self = $args;
	bless ($self, $class);

	$self->{path} = '.' if (!defined $self->{path});
	$self->{r_script} = $FindBin::Bin . "/../lib/perl5/Breseq/coverage_distribution.r";
	
	return $self; 
}


sub fit
{
	my (	$self,
			$distribution_file, 
			$plot_file, 
			$deletion_propagation_pr_cutoff, 
			$junction_coverage_pr_cutoff, 
			$junction_accept_pr_cutoff, 
			$junction_keep_pr_cutoff, 
			$junction_max_score
	) = @_;
		
	my $log_file_name = "$self->{path}/$$.r.log";
	my $command = "R --vanilla < $self->{r_script} > $log_file_name";
	$command .= " distribution_file=$distribution_file plot_file=$plot_file";
	$command .= " deletion_propagation_pr_cutoff=$deletion_propagation_pr_cutoff";
	$command .= " junction_coverage_pr_cutoff=$junction_coverage_pr_cutoff";
	$command .= " junction_accept_pr_cutoff=$junction_accept_pr_cutoff";
	$command .= " junction_keep_pr_cutoff=$junction_keep_pr_cutoff";
	$command .= " junction_max_score=$junction_max_score";

	Breseq::Shared::system($command);

	open ROUT, "<$log_file_name" or die;
	my @lines = <ROUT>;
	close ROUT;
	unlink $log_file_name;
	
	chomp @lines;
	@lines = grep s/^\[1\]\s+//, @lines;
#	print STDERR Dumper(@lines);
	
	return(@lines);
}

## helper functions

sub analyze_unique_coverage_distributions
{
	my ($settings, $summary, $ref_seq_info, $plot_key, $distribution_key) = @_;

	foreach my $seq_id (@{$ref_seq_info->{seq_ids}})
	{               
		analyze_unique_coverage_distribution($settings, $seq_id, $summary, $plot_key, $distribution_key);
	}
}
sub analyze_unique_coverage_distribution
{
	my ($settings, $seq_id, $summary, $plot_key, $distribution_key) = @_;

	##initialize summary information
	$summary->{unique_coverage}->{$seq_id}->{nbinom_size_parameter} = 'ND';
	$summary->{unique_coverage}->{$seq_id}->{nbinom_mean_parameter} = 'ND';
	$summary->{unique_coverage}->{$seq_id}->{nbinom_prob_parameter} = 'ND'; 
	$summary->{unique_coverage}->{$seq_id}->{average} = 1;
	$summary->{unique_coverage}->{$seq_id}->{variance} = 'ND';
	$summary->{unique_coverage}->{$seq_id}->{dispersion} = 'ND';
	$summary->{unique_coverage}->{$seq_id}->{deletion_coverage_propagation_cutoff} = 5;

	my $unique_only_coverage_plot_file_name = $settings->file_name($plot_key, {'@'=>$seq_id});
	my $unique_only_coverage_distribution_file_name = $settings->file_name($distribution_key, {'@'=>$seq_id});

	### Define various coverage thresholds...	
	my $sequence_length = $summary->{sequence_conversion}->{reference_sequences}->{$seq_id}->{length};

	### DELETION PROPAGATION CUTOFF
	## One-tailed test p=0.05, Bonferroni correction
	# my $del_propagation_pr_cutoff = 0.05 / $sequence_length;

	## One-tailed test p=0.01, no Bonferroni correction
	#my $del_propagation_pr_cutoff = 0.01;

	## We really want somewhere between these two, try this...
	my $deletion_propagation_pr_cutoff = 0.05 / sqrt($sequence_length);

	### NEW JUNCTION COVERAGE CUTOFFS
	## Arbitrary value that seems to work....
	my $junction_coverage_pr_cutoff = 1/$sequence_length; # *0.05

	## We really want somewhere between these two, try this...
	my $junction_accept_pr_cutoff = 0.01;
	my $junction_keep_pr_cutoff = 0.01 / sqrt($sequence_length);
	my $junction_max_score = int(2 * $summary->{sequence_conversion}->{avg_read_length});

	my $dist = Breseq::CoverageDistribution->new({});
	my @lines = $dist->fit($unique_only_coverage_distribution_file_name, $unique_only_coverage_plot_file_name,
	        $deletion_propagation_pr_cutoff, $junction_coverage_pr_cutoff, $junction_accept_pr_cutoff, $junction_keep_pr_cutoff, $junction_max_score);

	#First two lines are negative binomial parameters.
	#Next three lines are average, standard deviation, and index of overdispersion

	#Put these into summary
	$summary->{unique_coverage}->{$seq_id}->{nbinom_size_parameter} = $lines[0];
	$summary->{unique_coverage}->{$seq_id}->{nbinom_mean_parameter} = $lines[1]; 
	#Calculated by formula, prob = size/(size + mu)
	$summary->{unique_coverage}->{$seq_id}->{nbinom_prob_parameter} = $lines[0] / ($lines[0] + $lines[1]); 
	$summary->{unique_coverage}->{$seq_id}->{average} = $lines[2];
	$summary->{unique_coverage}->{$seq_id}->{variance} = $lines[3];
	$summary->{unique_coverage}->{$seq_id}->{dispersion} = $lines[4];

	$summary->{unique_coverage}->{$seq_id}->{deletion_coverage_propagation_cutoff} = $lines[5];
	$summary->{unique_coverage}->{$seq_id}->{junction_coverage_cutoff} = $lines[6];
	$summary->{unique_coverage}->{$seq_id}->{junction_accept_score_cutoff} = $lines[7];
	$summary->{unique_coverage}->{$seq_id}->{junction_keep_score_cutoff} = $lines[8];
}


return 1;