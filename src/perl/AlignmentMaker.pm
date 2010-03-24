###
# Pod Documentation
###

=head1 NAME

AlignmentMaker.pm

=head1 SYNOPSIS

Module for making predictions from hybrid reads.

=head1 DESCRIPTION

=head1 AUTHOR

Jeffrey Barrick
<jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2008.  All rights reserved.

=cut

###
# End Pod Documentation
###

package AlignmentMaker;
use strict;
use CGI qw/:standard *table *Tr *code start_b end_b/;

require Exporter;
our @ISA = qw( Exporter );
our @EXPORT = qw();

use Bio::Root::Root;

use MummerDelta;
use FastqLite;
use BreseqUtility;

use Data::Dumper;

=head2 add_aligned_reads_to_intervals

 Title   : add_aligned_reads_to_intervals
 Usage   : $diff = add_aligned_reads_to_intervals( -file_name => 'sequence.HCDiff.txt' );
 Function: This is a data preparation step, aligned reads are added as {aligned_reads} to each 
           interval hash item which must have {start} and {end} tags.
 Returns : 

=cut

our $maximum_to_align = 1000;

sub add_aligned_reads_to_intervals
{
	my ($settings, $summary, $interval_list_ref) = @_;
	my $max = $settings->{alignment_read_limit};

	## create empty arrays and check to be sure all required fields are there
	my $num_intervals = scalar @$interval_list_ref;
	print STDERR "  Loading reads aligned to $num_intervals intervals...\n";
	INTERVAL: foreach my $interval (@$interval_list_ref)
	{
		$interval->{aligned_reads} = [] if (!defined $interval->{aligned_reads});
	
		$interval->{original_start} = $interval->{start} if (!defined $interval->{original_start});
		$interval->{original_end} = $interval->{end} if (!defined $interval->{original_end});

		if ( (!defined $interval->{start}) || (!defined $interval->{end}) || (!defined $interval->{seq_id}) )
		{
			print STDERR "Required information missing for interval.\n";
			print Dumper($interval);
			die;
		}
	}
	
	## go through each original fastq file
	## using the same encoding of read names as we did
	## originally, so that we can match up and restore the original names
	READ_FILE: foreach my $read_file ($settings->read_files)
	{		
		my $fastq_file_name = $settings->read_file_to_fastq_file_name($read_file);
		$fastq_file_name = $settings->file_name('trimmed_fastq_file_name', {'#'=>$read_file}) if ($settings->{trim_reads});
		my $fastq_file_index = $settings->read_file_to_fastq_file_index($read_file);

		#open the fastq file 
		my $in = FastqLite->new( -file => "<$fastq_file_name");
		my $sequence_index = 1;
		
		my $seq = $in->next_seq;
		my $encoded_read_name = "$fastq_file_index:$sequence_index";
		
		my $i = -1;
		INDEX: while (1) 
		{		
			$i++;
			my $this_mummer_delta_file_name = $settings->file_name('corrected_mummer_delta_file_name', {'#' => $read_file, '=' => $i});	
			last INDEX if (!-e $this_mummer_delta_file_name);	
			my $md = MummerDelta->new(-file_name => "$this_mummer_delta_file_name");
			
			my @match_list = $md->get_next_list;
			READ: while (scalar @match_list > 0)
			{			
				next READ if (scalar @match_list > 1); ##don't add redundants		
				while ($seq && ($encoded_read_name ne $match_list[0]->{query}))
				{
					last READ if (defined $max && ($sequence_index >= $max));
					$seq = $in->next_seq;
					$sequence_index++;
					$encoded_read_name = "$fastq_file_index:$sequence_index";
					
					if ($sequence_index % 10000 == 0)
					{
						print STDERR "    $sequence_index reads processed.\n";
					}
				}
						
				MATCH: foreach my $m (@match_list)
				{			
					$m->{read_fastq} = $seq;
					#record all of the intervals that it overlaps;
					my @overlap_list;
					
					## replace the encoded read name with the actual name!!
					$m->{query} = $seq->{id}; 
					
					INTERVAL: foreach my $interval (@$interval_list_ref)
					{
						### there is a limit to how many reads we will align for practical purposes
						next INTERVAL if (scalar @{$interval->{aligned_reads}} >= $maximum_to_align);	
						#print STDERR "$interval->{reference} ne $m->{reference}\n";
							
						## nice try, but the match is to the wrong reference sequence	
						next if ($interval->{seq_id} ne $m->{reference});
				
						if ( ($interval->{original_start} >= $m->{reference_start}) && ($interval->{original_start} <= $m->{reference_end})
						  || ($interval->{original_end} >= $m->{reference_start}) && ($interval->{original_end} <= $m->{reference_end}) )
						{
							push @{$interval->{aligned_reads}}, $m;
							#print STDERR Dumper($m);
						}
					}
				} #end MATCH
			} continue {
				@match_list = $md->get_next_list;
			}
		} #end INDEX
	} #end READ_FILE
}

## already have the matches
sub add_read_fastq_to_intervals
{
	my ($settings, $summary, $interval_list_ref) = @_;
	print STDERR "  Loading reads aligned to hybrid intervals...\n";

	## Part 1: Go through all of the matches and mark which need 'read_fastq' assigned

	my %need_sequences;
	
	foreach my $interval (@$interval_list_ref)
	{
		foreach my $a (@{$interval->{aligned_reads}})
		{
			push @{$need_sequences{$a->{query}}}, $a;
		}
	}
	
	## Part 2: Go through fastq file
	
	## go through each original fastq file
	## using the same encoding of read names as we did
	## originally, so that we can match up and restore the original names
	READ_FILE: foreach my $read_file ($settings->read_files)
	{
		my $fastq_file_name = $settings->read_file_to_fastq_file_name($read_file);
		$fastq_file_name = $settings->file_name('trimmed_fastq_file_name', {'#'=>$read_file}) if ($settings->{trim_reads});
		my $fastq_file_index = $settings->read_file_to_fastq_file_index($read_file);

		#open the fastq file 
		my $in = FastqLite->new( -file => "<$fastq_file_name");
		my $sequence_index = 1;
		
		while (my $seq = $in->next_seq)
		{
			if ($sequence_index % 10000 == 0)
			{
				print STDERR "    $sequence_index reads processed.\n";
			}
			
			my $encoded_read_name = "$fastq_file_index:$sequence_index";
			$sequence_index++;
			
			foreach my $m (@{$need_sequences{$encoded_read_name}})
			{
				$m->{query} = $seq->{id};
				$m->{read_fastq} = $seq;
			}
		}
	}	
}

=head2 text_alignment

 Title   : text_alignment
 Usage   : text_alignment( $interval );
 Function: 
 Returns : text version of alignment
 
=cut

our $end_gap_char = '-';
our $internal_gap_char = '.';

sub text_alignment
{
	my $debug = 0;
	my ($interval) = @_; # ref_seq_ref field is to a string, not a Bio::Seq
	
	my $alignment;
	@{$alignment->{indels}} = ();
	
	if ( (!defined $interval->{aligned_reads}) || (scalar @{$interval->{aligned_reads}} == 0) )
	{
		print STDERR "WARNING!\nNo reads aligned to interval: start = $interval->{start}, end = $interval->{end}\n";
		return $alignment;
	}
	
	foreach my $m (@{$interval->{aligned_reads}})
	{
		#print Dumper($m);
	
		#note, deep copy necessary because match pointer may be used by other intervals
		my $new_line;
		
		$new_line->{strand} = $m->{strand};
		$new_line->{reference_end} = $m->{reference_end};
		$new_line->{reference_start} = $m->{reference_start};
		$new_line->{query_start} = $m->{query_start};
		$new_line->{query_end} = $m->{query_end};
		$new_line->{query_length} = $m->{query_length};
		$new_line->{mismatches} = $m->{mismatches};
		$new_line->{query} = $m->{query};
		@{$new_line->{original_indels}} = @{$m->{indels}};
		@{$new_line->{indels}} = ();
		
		$new_line->{aligned_seq} = ( ($m->{strand} == -1) ? FastqLite::revcom($m->{read_fastq}->{seq}) : $m->{read_fastq}->{seq});
		
		@{$new_line->{quals}} = FastqLite::quals($m->{read_fastq});
		@{$new_line->{quals}} = reverse @{$new_line->{quals}} if ($m->{strand} == -1);
		
		#negative indels indicate an insertion in this sequence relative to the reference and other sequence
		#  --> need to insert gaps in OTHER sequences if they don't already have one BEFORE this genomic location
		#positive indels indicate a deletion in this sequence relative to the reference and other sequences
		#  --> only need to remember these for inserting gaps in THIS sequence
		
		#correct indels to refer to genome positions 
		my $on_position = $new_line->{reference_start};
		foreach my $indel (@{$m->{indels}})
		{
			my $sign = ($indel > 0) ? +1 : -1;
			$on_position += abs ($indel-1);	
			
			## deletions in reference sequence are BEFORE previous position
			$on_position-- if ($indel < 0);
			
			#deletions in reference
			$on_position -= 1 if ($sign == -1);
					
			push @{$new_line->{indels}}, $on_position * $sign;
		}
 
		#all indel lists are sorted by reference genome position
		my @new_indels = (); #in new read, not in alignment
		my @missing_indels = (); #already in alignment
		my @existing_indels = (); #in read and already in alignment
		
		#update
		my @read_indels;
		push @read_indels, @{$new_line->{indels}};
		my @alignment_indels;
		push @alignment_indels, @{$alignment->{indels}};
		
		my $next_read_indel = shift @read_indels;
		my $next_alignment_indel = shift @alignment_indels;
		
		#deletion in reference (if negative) or deletion in query (if positive)
		#go through in order and classify as new in read, new in alignment, or existing
		while ((defined $next_read_indel) || (defined $next_alignment_indel))
		{
			if ($debug)
			{	
				print STDERR "Testing...";
				print STDERR (defined $next_read_indel ? " $next_read_indel": "none" );
				print STDERR (defined $next_alignment_indel ? " $next_alignment_indel": "none" );
				print STDERR "\n";
			}
		
			#read list ran out, remaining are in alignment, but not read
			if (!defined $next_read_indel)
			{
				push @missing_indels,  $next_alignment_indel;
				$next_alignment_indel = shift @alignment_indels;
			}
			#alignment list ran out, remaining are new in read
			elsif (!defined $next_alignment_indel)
			{
				push @new_indels,  $next_read_indel;
				$next_read_indel = shift @read_indels;
			}
			#exists in both read and alignment
			elsif ( $next_read_indel == $next_alignment_indel)
			{
				push @existing_indels, $next_read_indel;
				$next_read_indel = shift @read_indels;
				$next_alignment_indel = shift @alignment_indels;
			}
			#exists only in read
			elsif ( (abs($next_read_indel) < abs($next_alignment_indel)) 
				|| ( (abs($next_read_indel) == abs($next_alignment_indel)) && ($next_read_indel < $next_alignment_indel) ) )
			{
				push @new_indels,  $next_read_indel;
				$next_read_indel = shift @read_indels;
			}
			#exists only in alignment
			else
			{
				push @missing_indels,  $next_alignment_indel;
				$next_alignment_indel = shift @alignment_indels;
			}
		}	
		
		print "both read and alignment @existing_indels\n" if ($debug);
		print "read, not alignment @new_indels\n" if ($debug);
		print "alignment, not read @missing_indels\n" if ($debug);
			

		#remember, we must align the whole read (and display which part is aligned)		
		$new_line->{extra_start} = ($m->{strand} == +1) ? $m->{query_start} - 1 : $m->{query_length} - $m->{query_end};
		$new_line->{read_start} = $new_line->{reference_start} - $new_line->{extra_start};		
		$alignment->{start} = $new_line->{read_start} if (!defined $alignment->{start});
		
		#alignment start changed, add gaps to all sequences
		if ($new_line->{read_start} < $alignment->{start})
		{
			foreach my $line (@{$alignment->{lines}})
			{
				$line->{aligned_seq} = ($end_gap_char x ($alignment->{start} - $new_line->{read_start})) . $line->{aligned_seq}
			}
			$alignment->{start} = $new_line->{read_start};
		}

		#remember, we must align the whole read (and display which part is aligned)
		$new_line->{extra_end} = ($m->{strand} == +1) ? $m->{query_length} - $m->{query_end} : $m->{query_start} - 1;
		$new_line->{read_end} = $new_line->{reference_end} + $new_line->{extra_end};		
		$alignment->{end} = $new_line->{read_end} if (!defined $alignment->{end});
		
		#alignment end changed, add gaps to all sequences
		if ($alignment->{end} < $new_line->{read_end})
		{
			foreach my $line (@{$alignment->{lines}})
			{
				$line->{aligned_seq} .= ($end_gap_char x ($new_line->{read_end} - $alignment->{end}));
			}
			$alignment->{end} = $new_line->{read_end};
		}
		
		## add gaps to this line
		#keep positive in this sequence and negative in other sequences!
		my @insert_read_indels = grep {$_>0} @existing_indels; #alignment and read
		push @insert_read_indels, grep {$_>0} @new_indels;     #read, not alignment
		push @insert_read_indels, grep {$_<0} @missing_indels; #alignment, not read
		@insert_read_indels = sort { (abs($a) <=> abs($b)) || ($a <=> $b) } @insert_read_indels;
		
		#if there is a negative in read, then we need to downshift where inserts from others occur 
		my @offset_read_indels = grep {$_<0} @new_indels;
		push @offset_read_indels, grep {$_<0} @existing_indels;
		my $offset_read_indel = shift @offset_read_indels;
		
		my $gaps_added = 0;
		print "Insert read indels @insert_read_indels\n" if ($debug);
		foreach my $ref_pos (@insert_read_indels)
		{
			while (defined $offset_read_indel && (abs($offset_read_indel) < abs($ref_pos)))
			{
				$gaps_added++;
				$offset_read_indel = shift @offset_read_indels;
			}
		
			my $insert_position = $gaps_added + abs($ref_pos) - $new_line->{read_start} - 1;
			if ($insert_position < 1)
			{
				substr $new_line->{aligned_seq}, 0, 0, $end_gap_char;
			}
			elsif ($insert_position >= length $new_line->{aligned_seq})
			{
				substr $new_line->{aligned_seq}, length($new_line->{aligned_seq}), 0, $end_gap_char;
			}
			else
			{
				#CHANGED, removed -1 to position
				substr $new_line->{aligned_seq}, $gaps_added + abs($ref_pos) - $new_line->{read_start}, 0, $internal_gap_char;
			}
			$gaps_added++;
		}
		
		## add gaps to other sequences (negatives in new_indels)
		my @insert_alignment_indels = grep {$_<0} @new_indels;
		@insert_alignment_indels = sort { (abs($a) <=> abs($b)) || ($a <=> $b) } @insert_alignment_indels;
		
		my @existing_alignment_indels;
		push @existing_alignment_indels, grep {$_<0} @missing_indels;
		push @existing_alignment_indels, grep {$_<0} @existing_indels;
		@existing_alignment_indels = sort { (abs($a) <=> abs($b)) || ($a <=> $b) } @existing_alignment_indels;

		print "Insert alignment indels @insert_alignment_indels\n" if ($debug);

		my $existing_alignment_indel = shift @existing_alignment_indels;
		$gaps_added = 0;		
		foreach my $ref_pos (@insert_alignment_indels)
		{
			#these will be in the SAME column for all other entries
			#which we can figure out from existing gaps and alignment start, end

			#count how many of existing indels are less than this one and add them
			
			while (defined $existing_alignment_indel && (abs($existing_alignment_indel) < abs($ref_pos)))
			{
				$gaps_added++;
				$existing_alignment_indel = shift @existing_alignment_indels;
			}

			foreach my $line (@{$alignment->{lines}})
			{			
				my $insert_char = ( ($line->{read_start} <= abs($ref_pos)) && ($line->{read_end} >= abs($ref_pos)) )
					? $internal_gap_char : $end_gap_char;
								#CHANGED, removed -1 to position
	
				my $insert_pos = $gaps_added + abs($ref_pos) - $alignment->{start};	
				$insert_pos = 0 if ($insert_pos < 0);
				$insert_pos = length $line->{aligned_seq} if ($insert_pos > length $line->{aligned_seq});
				substr $line->{aligned_seq}, $insert_pos, 0, $insert_char;
			}
		}		
		
		#add end gaps to this line		
 		$new_line->{aligned_seq} = ($end_gap_char x ($new_line->{read_start} - $alignment->{start}) )
			.  $new_line->{aligned_seq}
			. ($end_gap_char x ($alignment->{end} - $new_line->{read_end}) );
			
		
		push @{$alignment->{lines}}, $new_line;
		
		# keep all indels and re-sort
		push @{$alignment->{indels}}, @new_indels; 
		@{$alignment->{indels}} = sort { (abs($a) <=> abs($b)) || ($a <=> $b) } @{$alignment->{indels}};


		print Dumper($alignment) if ($debug);
		foreach my $line (@{$alignment->{lines}})
		{
			print "$line->{aligned_seq}\n" if ($debug);
		}
	}
	
	## REFERENCE SEQUENCE(s) and CHANGE LINE
	##
	## Get reference sequence spanning aligned region
	
	add_change_sequence($alignment, $interval);
	
	## take care of bounds of alignment passing outside of the reference sequence	
	if (defined $interval->{ref_seq_ref})
	{
		## all aligned reads should have the same reference sequence name...
		## NORMAL MODE
		if (!$interval->{alignment_reference_info_list})
		{
			add_reference_sequence($alignment, $interval);
		}
		## MODE ALLOWING SPLIT REFERENCE SEQUENCES
		else 
		{
			foreach my $alignment_reference_info (@{$interval->{alignment_reference_info_list}})
			{
				add_reference_sequence($alignment, $interval, $alignment_reference_info);
			}
		}
	}
		
	## not strictly necessary, stitch together lines
	$alignment->{text} = '';
	foreach my $line (@{$alignment->{lines}})
	{
		$alignment->{text} .= $line->{aligned_seq} . "\n";
	}
	
	return $alignment;
}

sub add_change_sequence
{
	my ($alignment, $interval) = @_;
	
	##mark the original start if we have shifter to later in a gene on the opposite strand
	my $start = (defined $interval->{original_start}) ? $interval->{original_start} : $interval->{start};
	my $end = (defined $interval->{original_end}) ? $interval->{original_end} : $interval->{end};

	$alignment->{aligned_change_seq} = '';
	return if ($interval->{alignment_empty_change_line});
	
	$alignment->{aligned_change_seq} = ' ' x ($alignment->{end} - $alignment->{start} + 1);
	#could use different symbols for whether start < end or start > end
	foreach my $star ($start..$end)
	{
		substr $alignment->{aligned_change_seq}, $star - $alignment->{start}, 1, '|';
	}
	foreach my $star ($end..$start)
	{
		substr $alignment->{aligned_change_seq}, $star - $alignment->{start}, 1, '|';
	}
	
	
	$alignment->{aligned_change_seq} = add_gaps_to_sequence($alignment, $alignment->{aligned_change_seq}, ' ');
}

sub add_reference_sequence
{
	my ($alignment, $interval, $alignment_reference_info) = @_;

	my $ref_seq_ref = $interval->{ref_seq_ref};
	my $ref_seq_length = length $$ref_seq_ref;
		
	my $ref_seq_start = $alignment->{start};
	my $ref_seq_end = $alignment->{end};

	my $add_start = ($ref_seq_start - 1);
	$ref_seq_start = 1 if ($add_start < 0);

	my $add_end = ($ref_seq_end - $ref_seq_length);
	$ref_seq_end = $ref_seq_length if ($add_end > 0);

#print Dumper($interval->{aligned_reads});
#print Dumper($interval->{seq_id});

	my $item;
	$item->{seq_id} = $interval->{aligned_reads}->[0]->{reference};
	$item->{start} = $ref_seq_start;
	$item->{end} = $ref_seq_end;
	$item->{label} = "$item->{seq_id}/$item->{start}-$item->{end}";
	$item->{aligned_seq} = substr $$ref_seq_ref, $ref_seq_start-1, $ref_seq_end - $ref_seq_start + 1;
	
	if ($alignment_reference_info->{truncate_start})
	{				
		my $delete_size = $alignment_reference_info->{truncate_start} - $item->{start};
		if ($delete_size > 0)
		{
			substr $item->{aligned_seq}, 0, $delete_size, ($end_gap_char x $delete_size);
			$item->{start} = $alignment_reference_info->{truncate_start};
		}
		my $new_size = $item->{end} - $item->{start} + 1;
		my $renumber_start = $alignment_reference_info->{ghost_start};
		my $renumber_end = $renumber_start + $alignment_reference_info->{ghost_strand} * ($new_size - 1);
		my $actual_ref = $item->{seq_id};
		$actual_ref = $alignment_reference_info->{ghost_seq_id} if (defined $alignment_reference_info->{ghost_seq_id});
		
		$item->{label} = "$actual_ref/$renumber_start-$renumber_end";
	}
	
	if ($alignment_reference_info->{truncate_end})
	{
		my $delete_size = $item->{end} - $alignment_reference_info->{truncate_end};
		if ($delete_size > 0)
		{
			substr $item->{aligned_seq}, -$delete_size, $delete_size, ($end_gap_char x $delete_size);
			$item->{end} = $alignment_reference_info->{truncate_end};
		}
		
		my $new_size = $item->{end} - $item->{start} + 1;
		my $renumber_end = $alignment_reference_info->{ghost_end};
		my $renumber_start = $renumber_end + $alignment_reference_info->{ghost_strand} * ($new_size - 1);
		my $actual_ref = $item->{seq_id};
		$actual_ref = $alignment_reference_info->{ghost_seq_id} if (defined $alignment_reference_info->{ghost_seq_id});
		
		$item->{label} = "$actual_ref/$renumber_start-$renumber_end";
	}

	$item->{aligned_seq} = ($end_gap_char x -$add_start) . $item->{aligned_seq} if ($add_start < 0);
	$item->{aligned_seq} = $item->{aligned_seq} . ($end_gap_char x $add_end)  if ($add_end > 0);
	
	$item->{aligned_seq} = add_gaps_to_sequence($alignment, $item->{aligned_seq});
		
	push @{$alignment->{reference_list}}, $item;
}

sub add_gaps_to_sequence
{
	my ($alignment, $sequence, $gap_char) = @_;
	$gap_char = $internal_gap_char if (!defined $gap_char);
	
	my $new_sequence = $sequence;
	my $gaps_added = 0;
	
	#print STDERR Dumper($alignment, $sequence, $gap_char);
	
	
	foreach my $ref_pos (@{$alignment->{indels}})
	{
		next if ($ref_pos > 0);
		##negative insertions (deletions in reference occur BEFORE the reference position)
		substr $new_sequence, $gaps_added + abs($ref_pos) - $alignment->{start}, 0, $gap_char;
		$gaps_added++;
	}	
	return $new_sequence;
}

=head2 set_quality_range

 Title   : set_quality_range
 Usage   : set_quality_range( min, max );
 Function: 
 Returns : set minimum and maximum qualities
 
=cut

##so these parameters can be set from outside
our $max_quality;
our $min_quality;
our $quality_type;
our @qual_cutoffs;

sub set_quality_range
{
	my ($_quality_type, $_min_quality, $_max_quality, $_quality_cdf_ref, $verbose) = @_;

	$quality_type = $_quality_type;
	$min_quality = $_min_quality;
	$max_quality = $_max_quality;
	
	#exclude max and min scores
	my @cutoff_levels = (0.05, 0.5, 0.95);
	my $not_ends_pr = 1-($_quality_cdf_ref->[-1] - $_quality_cdf_ref->[-2] + $_quality_cdf_ref->[0]);
	
	#print STDERR "$not_ends_pr\n";
	if ($not_ends_pr > 0 )
	{
		#find 5%, 50%, 5% thresholds from cdf
		@qual_cutoffs = ();
		for (my $i=0; $i<scalar @$_quality_cdf_ref-1; $i++)
		{
			my $test_pr = ($_quality_cdf_ref->[$i]-$_quality_cdf_ref->[0])/$not_ends_pr;
			
			#print STDERR "$i $test_pr\n";
			for (my $j=0; $j<scalar @cutoff_levels; $j++)
			{
				if ( (!defined $qual_cutoffs[$j]) && ($test_pr >$cutoff_levels[$j]) )
				{
					$qual_cutoffs[$j] = $min_quality+$i;
				}
			}
		}
	}
		
	## assign any that weren't
	for (my $j=0; $j<scalar @cutoff_levels; $j++)
	{
		$qual_cutoffs[$j] = $min_quality + (scalar @$_quality_cdf_ref) if (!defined $qual_cutoffs[$j]) ;
	}
	
	#@qual_cutoffs = ( $min_quality, $min_quality + int(0.5*($max_quality - $min_quality)), $max_quality );	
	print STDERR "quality range cutoffs @qual_cutoffs\n" if ($verbose);
}


=head2 html_alignment

 Title   : html_alignment
 Usage   : html_alignment( $interval );
 Function: 
 Returns : html version of alignment
 
=cut

our $base_colors_hash = {
	'G' => [ "rgb(210,210,210)", "rgb(140,140,140)", "rgb(70,70,70)", "rgb(0,0,0)" ],
	'C' => [ "rgb(120,120,255)", "rgb(60,60,255)", "rgb(0,0,255)", "rgb(0,0,150)" ],
	'A' => [ "rgb(255,180,180)", "rgb(255,100,100)", "rgb(255,20,20)", "rgb(200,0,0)" ],
	'T' => [ "rgb(180,255,180)", "rgb(100,255,100)", "rgb(20,255,20)", "rgb(0,200,0)" ],
};


sub html_alignment
{
	my ($interval) = @_;
	
	## everything is built from the text alignment
	my $alignment = text_alignment($interval);
	
	## htmlize the change sequence	
	if (defined $alignment->{aligned_change_seq})
	{
		$alignment->{aligned_html_change_seq} = $alignment->{aligned_change_seq};
		$alignment->{aligned_html_change_seq} =~ s/\|/&darr;/g;
		$alignment->{aligned_html_change_seq} =~ s/ /&nbsp;/g;
		$alignment->{aligned_html_change_seq} = code($alignment->{aligned_html_change_seq});
	}
	
	## htmlize the reference sequence
	foreach my $r (@{$alignment->{reference_list}})
	{		
		$r->{aligned_html_seq} = htmlize_aligned_sequence($r->{aligned_seq});
	}
	
	## go through each sequence and HTMLize it.
	for my $line (@{$alignment->{lines}})
	{	
		$line->{aligned_html_seq} = htmlize_aligned_sequence($line->{aligned_seq}, $line->{strand}, $line->{query_start}, $line->{query_end}, $line->{query_length}, $line->{quals}, \@qual_cutoffs);
	}	
	
	## create legend information
	$alignment->{legend} = htmlize_aligned_sequence('ATCG', +1, 1, 4, 4, [0,0,0,0], \@qual_cutoffs);
	
	foreach (my $i=0; $i<scalar @qual_cutoffs; $i++)
	{
		my $cutoff = $qual_cutoffs[$i];
		$alignment->{legend} .= " &lt; $cutoff &le ";
		$alignment->{legend} .= htmlize_aligned_sequence('ATCG', +1, 1, 4, 4, [$qual_cutoffs[$i], $qual_cutoffs[$i], $qual_cutoffs[$i], $qual_cutoffs[$i]], \@qual_cutoffs);
	}
	
	return $alignment;
}

sub htmlize_aligned_sequence
{
	my ($text_seq, $strand, $start_match, $end_match, $query_length, $quals_ref, $qual_cutoffs_ref ) = @_;
	my $html_seq;
	$html_seq .= start_code;
	my $ungapped_char_count = 0;
	
	#query start and $query end need to be adjusted if on opposite strand
	if (defined ($strand) && ($strand == -1))
	{
		($start_match, $end_match) = ($query_length - $end_match+1, $query_length - $start_match+1);
	}
		
	foreach my $c (split //, $text_seq)
	{
		## make dashed ungapped
		if ($c eq $end_gap_char)
		{
			$html_seq .= "&#8209;";
		}
		elsif ($c eq $internal_gap_char)
		{
			$html_seq .= $c;
			##nothing
		}
		## add color to characters to signify quality
		else
		{
			##color by base and quality if provided
			my $qual;
			my $color;
			if (defined $quals_ref && ($c =~ m/[ATCG]/))
			{
				$qual = $quals_ref->[$ungapped_char_count];

				my $color_num = 0;
				
				while ((defined $qual_cutoffs_ref->[$color_num]) && ($qual >= $qual_cutoffs_ref->[$color_num]))
				{
					$color_num++;
				}	
				
				$color = $base_colors_hash->{$c}->[$color_num];
			}
			## no base quality provided -- assume BEST
			else
			{
				$color = $base_colors_hash->{$c}->[-1];
			}
		
			## outside of matched portion of read
			if ( (!defined $color)
			  || ((defined $start_match) && ($ungapped_char_count + 1 < $start_match))
			  || ((defined $end_match)   && ($ungapped_char_count + 1 > $end_match)) )
			{
				$html_seq .= font({-style => "color: rgb(0,0,0); background-color: rgb(255,255,255);"}, $c);
			}
			else
			{
				#print STDERR "$c $qual $color\n";
				$html_seq .= font({-style => "color: rgb(255,255,255); background-color:$color;"}, $c);
			}
			
			$ungapped_char_count++;

		}
	}
	$html_seq .= end_code;

	return $html_seq;
}

return 1;

