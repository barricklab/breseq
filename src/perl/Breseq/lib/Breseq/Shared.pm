###
# Pod Documentation
###

=head1 NAME

Breseq::Shared

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

package Breseq::Shared;
use strict;

require Exporter;
our @ISA = qw( Exporter );
our @EXPORT = qw( );

use Data::Dumper;


sub system
{
	my ($command, $silent, $continue) = @_;	
	print STDERR "[system] $command\n" if (!$silent);
	my $res = CORE::system $command;
	print STDERR "Error: $!\nResult code: $res\n" if ($res);
	die if ($res && !$continue);
	return $res;
}

sub capture_system
{
	my ($command, $silent, $continue) = @_;	
	print STDERR "[system] $command\n" if (!$silent);
	$command .= " > $$.script_output";
	my $res = CORE::system $command;
	print STDERR "Error: $!\nResult code: $res\n" if ($res);
	die if ($res && !$continue);
	
	open IN, "$$.script_output";
	my @lines = <IN>;
	close IN;
	unlink("$$.script_output");
	return join("", @lines);
}


sub region_to_coords
{
	my ($region, $reference_length) = @_;
	die if (!defined $region);
	
	my ($seq_id, $start, $end);
	my ($insert_start, $insert_end) = (0, 0);
	#syntax that includes insert counts
	# e.g. NC_001416:4566.1-4566.1
		
	#strip commas!
	$region =~ s/,//g;
	if ($region =~ m/(.+)\:(\d+)\.(\d+)-(\d+)\.(\d+)/)
	{	
		($seq_id, $start, $insert_start, $end, $insert_end) = ($1, $2, $3, $4, $5);		
	}
	elsif ($region =~ m/(.+)\:(\d+)(\.\.|\-)(\d+)/)
	{
		($seq_id, $start, $end) = ($1, $2, $4);
	}
	else
	{
		($seq_id, $start, $end) = split /:|\.\.|\-/, $region;
	}
		
	die "Deprecated usage!" if (ref(\$reference_length) ne "SCALAR");
	
	($start, $end) = (1, $reference_length) if (!defined $start && !defined $end);
	$end = $start if (!defined $end);

	##check the start and end for sanity....	
	$start = 1 if ($start < 1); 
	$end = 1 if ($end < 1); 
	$end = $reference_length if ($end > $reference_length); 
	$start = $reference_length if ($start > $reference_length); 

#	die "Problem parsing region: \'$region\'\n" if ($start > $end);

	## return cleaned up region
	$region = "$seq_id:$start-$end";
	return ($seq_id, $start, $end, $insert_start, $insert_end, $region);
}

sub revcom
{
	my ($seq) = @_;
	$seq =~ tr/ATCG/TAGC/;
	return reverse $seq;
}

sub polymorphism_statistics
{
	our ($settings, $summary, $ref_seq_info) = @_;

	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
	my @seq_ids = @{$ref_seq_info->{seq_ids}};
       
	## some local variable lookups for convenience
	my $total_ref_length = 0;
	foreach my $seq_id (@seq_ids)
	{
		$total_ref_length += length $ref_seq_info->{ref_strings}->{$seq_id};
	}
	my $log10_ref_length = log($total_ref_length) / log(10);

	##
	## Replacement for below
	##
	## ToDo: This should really make a different column for each input read set.
	##
	my $coverage_fn = $settings->file_name('unique_only_coverage_distribution_file_name', {'@'=>""});
	my $outputdir = `dirname $coverage_fn`;
	chomp $outputdir; $outputdir .= '/';
	my $count_file_name = $outputdir .= "error_counts.tab";


	open COUNT, "<$count_file_name";
	my $count_header_line = <COUNT>;
	$count_header_line = <COUNT>;
	chomp $count_header_line;
	my @count_header_list = split /\t/, $count_header_line; 
	my $quality_column;
	my $count_column;
	for (my $i=0; $i<scalar @count_header_list; $i++)
	{
	          $quality_column = $i if ($count_header_list[$i] eq 'quality');
	          $count_column = $i if ($count_header_list[$i] eq 'count');
	}

	#print "$count_column $quality_column\n";

	my @quality_count_list;
	while (my $line = <COUNT>)
	{
	        chomp $line;
	        #print "$_\n";
	        my @line_list = split /\t/, $line; 
	        $quality_count_list[$line_list[$quality_column]] += $line_list[$count_column];
	}
	close COUNT;

	my $genome_error_counts_file_name = $settings->file_name('genome_error_counts_file_name');

	open GEC, ">$genome_error_counts_file_name";
	for (my $i=1; $i<scalar @quality_count_list; $i++)
	{
	        my $val = 0;
	        $val = $quality_count_list[$i] if (defined $quality_count_list[$i]);
	        print GEC "$val\n";
	}
	close GEC;

	my $polymorphism_statistics_input_file_name = $settings->file_name('polymorphism_statistics_input_file_name');
	my $polymorphism_statistics_output_file_name = $settings->file_name('polymorphism_statistics_output_file_name');

	### Load the older GenomeDiff and add new fields
	my $ra_mc_genome_diff_file_name = $settings->file_name('ra_mc_genome_diff_file_name');
	my $gd = GenomeDiff->new( {in => $ra_mc_genome_diff_file_name} );

	my $polymorphism_statistics_r_script_file_name = $settings->file_name('polymorphism_statistics_r_script_file_name');
	my $polymorphism_statistics_r_script_log_file_name = $settings->file_name('polymorphism_statistics_r_script_log_file_name');
	my $total_reference_length = $summary->{sequence_conversion}->{total_reference_sequence_length};

	Breseq::Shared::system("R --vanilla total_length=$total_reference_length in_file=$polymorphism_statistics_input_file_name out_file=$polymorphism_statistics_output_file_name qual_file=$genome_error_counts_file_name < $polymorphism_statistics_r_script_file_name > $polymorphism_statistics_r_script_log_file_name");

	## Read R file and add new results corresponding to all columns
	open ROUT, "<$polymorphism_statistics_output_file_name" or die "Could not find file: $polymorphism_statistics_output_file_name";
	my $header = <ROUT>;
	chomp $header;
	my @header_list = split /\t/, $header;

	my $new_gd = GenomeDiff->new();
	foreach my $mut ($gd->list)
	{
		## lines only exist for RA evidence
		if ($mut->{type} ne 'RA')
		{
			$new_gd->add($mut);
			next;
		}

		## lines only exist for polymorphisms
		if (($mut->{frequency} == 1) || ($mut->{frequency} == 0))
		{
			$new_gd->add($mut);
			next;
		}

		my $line = <ROUT>;
		chomp $line;
		my @line_list = split /\t/, $line;

		for (my $i=0; $i< scalar @header_list; $i++)
		{
			$mut->{$header_list[$i]} = $line_list[$i];
			die "Incorrect number of items on line:\n$line" if (!defined $line_list[$i]);
		}

		## Evalue cutoff again (in case we are only running this part)
		GenomeDiff::add_reject_reason($mut, "EVALUE") if ($mut->{polymorphism_quality} < $settings->{polymorphism_log10_e_value_cutoff});

		## Frequency cutoff
		if ( ($mut->{frequency} < $settings->{polymorphism_frequency_cutoff}   )
		  || ($mut->{frequency} > 1-$settings->{polymorphism_frequency_cutoff} ) )
		{
			GenomeDiff::add_reject_reason($mut, "POLYMORPHISM_FREQUENCY_CUTOFF");                                   
		}

		## Minimum coverage on both strands
		my $polymorphism_coverage_limit_both_bases = $settings->{polymorphism_coverage_both_strands};
		my $passed = 1;
		my ($top,$bot) = split /\//, $mut->{ref_cov};
		$passed &&= $top >= $polymorphism_coverage_limit_both_bases;
		$passed &&= $bot >= $polymorphism_coverage_limit_both_bases;            
		($top,$bot) = split /\//, $mut->{new_cov};
		$passed &&= $top >= $polymorphism_coverage_limit_both_bases;
		$passed &&= $bot >= $polymorphism_coverage_limit_both_bases;

		GenomeDiff::add_reject_reason($mut, "POLYMORPHISM_STRAND") if (!$passed);
		GenomeDiff::add_reject_reason($mut, "KS_QUALITY_P_VALUE") if ($mut->{ks_quality_p_value} < $settings->{polymorphism_bias_p_value_cutoff});
		GenomeDiff::add_reject_reason($mut, "FISHER_STRAND_P_VALUE") if ($mut->{fisher_strand_p_value} < $settings->{polymorphism_bias_p_value_cutoff});

		###### Optionally, ignore if in a homopolymer stretch
		if (defined $settings->{polymorphism_reject_homopolymer_length})
		{
			my $test_length = 20;
			my $seq_id = $mut->{seq_id};
			my $end_pos = $mut->{position};
			my $start_pos = $end_pos - $test_length + 1;
			$start_pos = 1 if ($start_pos < 1);
			my $length = length $ref_seq_info->{ref_strings}->{$seq_id};     
			my $bases = substr($ref_seq_info->{ref_strings}->{$seq_id}, $start_pos-1, ($end_pos - $start_pos + 1) );

			#print Dumper($mut);
			#print "$bases\n";

			my $same_base_length = 0;
			my $first_base = substr($bases, $end_pos-$start_pos, 1);
			for (my $i=$end_pos; $i>=$start_pos; $i--)
			{
				my $this_base = substr($bases, $i-$start_pos, 1);
				last if ($first_base ne $this_base);
				$same_base_length++;
			}

			#print "$same_base_length\n";
			if ($same_base_length >= $settings->{polymorphism_reject_homopolymer_length})
			{
				GenomeDiff::add_reject_reason($mut, "HOMOPOLYMER_STRETCH");
			}
		}

		if ($mut->{reject} && ($mut->{polymorphism_quality} > $settings->{mutation_log10_e_value_cutoff}) && ($mut->{frequency} > 0.5) )
		{
			#print Dumper($mut);
			$mut->{frequency} = 1;
			delete $mut->{reject};

			## FIX -- need to re-evaluate whether it would have been accepted as a normal mutation 
			## This is NOT the right quality being used here. Need a separate quality for consensus call and polymorphism call!
			GenomeDiff::add_reject_reason($mut, "EVALUE") if ($mut->{polymorphism_quality} < $settings->{mutation_log10_e_value_cutoff});
		}

		$new_gd->add($mut);

		## END EXPERIMENTAL
	}

	### Write out the file which now has much more data
	my $polymorphism_statistics_ra_mc_genome_diff_file_name = $settings->file_name('polymorphism_statistics_ra_mc_genome_diff_file_name');
	$new_gd->write($polymorphism_statistics_ra_mc_genome_diff_file_name);
        
}

return 1;

