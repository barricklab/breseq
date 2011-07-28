###
# Pod Documentation
###

=head1 NAME

Breseq::CoverageOutput

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

package Breseq::CoverageOutput;
use vars qw(@ISA);

use FindBin;

use Breseq::Shared;
use Data::Dumper;


sub new
{
	my($caller,$args) = @_;
	my $class = ref($caller) || $caller;

	my $self = $args;
	bless ($self, $class);
	
	# These should be defined in $args
	#($self->{fasta}, $self->{bam}, $self->{path})	
	
	$self->{path} = '.' if (!defined $self->{path});
	$self->{r_script} = $FindBin::Bin . "/../lib/perl5/Breseq/plot_coverage.r";
	
	$self->{bam_path} = $self->{bam};
	$self->{fasta_path} = $self->{fasta};
			
	return $self; 
}


sub plot_coverage
{
	my ($self, $region, $output, $options) = @_;

	## Set up defaults
	$options = {} if (!defined $options);
	$options->{pdf} = 0 if (!defined $options->{pdf});
	$options->{resolution} = 600 if (!defined $options->{resolution});
	$options->{total_only} = 0 if (!defined $options->{total_only});
	$options->{shaded_flanking} = 0 if (!defined $options->{shaded_flanking}); #how much gray on each end
	die if (!$options->{reference_length});
	
	#print STDERR Dumper($options);
	
	my ($seq_id, $start, $end, $insert_start, $insert_end);
	($seq_id, $start, $end, $insert_start, $insert_end, $region)  = Breseq::Shared::region_to_coords($region, $options->{reference_length});
	$self->throw("Invalid region $region") if (!$seq_id || !$start || !$end);
		
	## extend the region and re-check
	my $extended_start = $start - $options->{shaded_flanking};
	$extended_start = 0 if ($extended_start < 0);
	my $extended_end = $end + $options->{shaded_flanking};
	
	#call to region_to_coords will fix the end if it has been extended too far
	
	
	my $extended_region = "$seq_id:$extended_start-$extended_end";
	my ($extended_seq_id, $extended_insert_start, $extended_insert_end);
	($extended_seq_id, $extended_start, $extended_end, $extended_insert_start, $extended_insert_end, $extended_region)  = Breseq::Shared::region_to_coords($extended_region, $options->{reference_length});
			
	my $size = $extended_end - $extended_start + 1;
	
	$output = $region if (!defined $output);
	$output =~ s/:/_/g;
	$output .= ($options->{pdf}) ? ".pdf" : ".png";
	
	## if no guidance provided about how much to downsample, aim for 1E2-1E3 total points.
	my $downsample = $options->{downsample};
	my $resolution = $options->{resolution};
	if (!defined $downsample)
	{
		$downsample = int($size/$resolution);
		$downsample = 1 if ($downsample < 1);
	}
	
#	print STDERR "Downsample: $downsample\n";
	
	my $tmp_coverage = "$self->{path}/$$.coverage.tab";

	my $bin_path = $FindBin::Bin;
	my $cmdline = "$bin_path/cbreseq TABULATE_COVERAGE --bam $self->{bam_path} --fasta $self->{fasta_path} --region $extended_region --output $tmp_coverage --downsample $downsample";
	Breseq::Shared::system($cmdline);

	my $log_file_name = "$self->{path}/$$.r.log";
	Breseq::Shared::system("R --vanilla in_file=$tmp_coverage out_file=$output pdf_output=$options->{pdf} total_only=$options->{total_only} window_start=$start window_end=$end < $self->{r_script} > $log_file_name");
	
	#clean up
	unlink $tmp_coverage;
	unlink $log_file_name;
}

return 1;

