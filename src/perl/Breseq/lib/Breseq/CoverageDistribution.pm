###
# Pod Documentation
###

=head1 NAME
Breseq::CoverageDistribution.pm

=head1 SYNOPSIS

Module for fitting and plotting the distribution of read-depth coverage.

=head1 DESCRIPTION

=head1 AUTHOR

Jeffrey Barrick
<jeffrey@barricklab.org>

=head1 COPYRIGHT

Copyright 2009.  All rights reserved.

=cut

###
# End Pod Documentation
###

use strict;

package Breseq::CoverageDistribution;

use vars qw(@ISA);
use Bio::Root::Root;
@ISA = qw( Bio::Root::RootI );

use Bio::DB::Sam;
use Breseq::Fastq;
use Breseq::Shared;
use Data::Dumper;
use FindBin;

sub new
{
	my($caller,@args) = @_;
	my $class = ref($caller) || $caller;
	my $self = new Bio::Root::Root($caller, @args);

	bless ($self, $class);
	($self->{path}) = $self->Bio::Root::RootI::_rearrange([qw(PATH)], @args);
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

return 1;