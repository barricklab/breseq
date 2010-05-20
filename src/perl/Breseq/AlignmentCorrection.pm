###
# Pod Documentation
###

=head1 NAME

MummerDelta.pm

=head1 SYNOPSIS

Module for reading and writing MUMmer delta files.

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

package Breseq::AlignmentCorrection;
require Exporter;

use strict;

our @ISA = qw( Exporter );
our @EXPORT = qw();
our $VERSION = '0.01';

use Data::Dumper;


#compare matches to candidate junctions with matches to original genome
sub correct_alignments
{
	my $verbose = 0;
	my ($settings, $summary, $ref_seq_info) = @_;
	my $gene_list_hash_ref = $ref_seq_info->{gene_lists};
	my $is_list_hash_ref = $ref_seq_info->{is_lists};
	my $flanking_length = $settings->{max_read_length};
				
	## for now we just use mapping qualities from ssaha2, but could load ref sequences this way
	my $reference_faidx_file_name = $settings->file_name('reference_fasta_file_name');
	my $reference_fai = Bio::DB::Sam::Fai->load($reference_faidx_file_name);
		
	my $candidate_junction_file_name = $settings->file_name('candidate_junction_fasta_file_name');
	## if there were no candidate junctions (file is empty) then we seg fault if we try to use samtools on it...
	$settings->{no_junction_prediction} = 1 if ( (!-e $candidate_junction_file_name) || (-s $candidate_junction_file_name == 0) );
	my	$candidate_junction_fai = Bio::DB::Sam::Fai->load($candidate_junction_file_name) if (!$settings->{no_junction_prediction});		
	
	my $minimum_best_score = 0;
	my $minimum_best_score_difference = 0;
	
	my %matched_junction;
	my %degenerate_matches;
	my $i = 0;
	
	my $resolved_reference_sam_file_name = $settings->file_name('resolved_reference_sam_file_name');
	my $RREF;
	open $RREF, ">$resolved_reference_sam_file_name" or die;
	
	my $resolved_junction_sam_file_name = $settings->file_name('resolved_junction_sam_file_name');
	my $RCJ;
	open $RCJ, ">$resolved_junction_sam_file_name" or die;
	
	my $reference_header;
	my $candidate_junction_header;
	foreach my $read_struct ($settings->read_structures)
	{	
		my $read_file = $read_struct->{base_name};	
		print STDERR "  READ FILE::$read_file\n";

		## setup the summary
		my $s;
		$s->{unmatched_reads} = 0;

		## Traverse the original fastq files to keep track of order
		## b/c some matches may exist in only one or the other file
		
		my @in_fastq;
		my @fastq_file_name;
		my @fastq_file_index;
		my @out_unmatched_fastq;
		my $i=0;
		for (my $i=0; $i < scalar @{$read_struct->{base_names}}; $i++)
		{
			my $this_read_file = $read_struct->{base_names}->[$i];			
			$fastq_file_name[$i] = $settings->read_file_to_fastq_file_name($this_read_file);	
			$fastq_file_index[$i] = $settings->read_file_to_fastq_file_index($this_read_file);		
				
			$in_fastq[$i] = Breseq::Fastq->new(-file => $fastq_file_name[$i]);
			
			if ($settings->{unmatched_reads})
			{				
				my $unmatched_file_name = $settings->file_name('unmatched_read_file_name', {'#'=>$this_read_file});
				$out_unmatched_fastq[$i] = Breseq::Fastq->new(-file => ">$unmatched_file_name");
			}
		}		

		my $reference_sam_file_name = $settings->file_name('reference_sam_file_name', {'#'=>$read_file});
		my $reference_tam = Bio::DB::Tam->open($reference_sam_file_name) or die "Could not open $reference_sam_file_name";
		my $reference_faidx_file_name = $settings->file_name('reference_faidx_file_name');
		$reference_header = $reference_tam->header_read2($reference_faidx_file_name) or throw("Error reading reference fasta index file: $reference_faidx_file_name");		
		my $reference_al;
		my $last_reference_alignment;
		
		
		my $candidate_junction_tam;
		if (!$settings->{no_junction_prediction})
		{
			my $candidate_junction_sam_file_name = $settings->file_name('candidate_junction_sam_file_name', {'#'=>$read_file});
			$candidate_junction_tam = Bio::DB::Tam->open($candidate_junction_sam_file_name) or die " Could not open candidate junction SAM file\n";
			my $candidate_junction_faidx_file_name = $settings->file_name('candidate_junction_faidx_file_name');
			$candidate_junction_header = $candidate_junction_tam->header_read2($candidate_junction_faidx_file_name) or die("Error reading reference fasta index file: $candidate_junction_faidx_file_name");			
		}
		
		my $candidate_junction_al;
		my $last_candidate_junction_alignment;		
				
		#proceed through all of the alignments
		if (!$settings->{no_junction_prediction})
		{
			($candidate_junction_al, $last_candidate_junction_alignment) 
				= Breseq::Shared::tam_next_read_alignments($candidate_junction_tam, $candidate_junction_header, $last_candidate_junction_alignment);		
		}
		
		($reference_al, $last_reference_alignment) 
			= Breseq::Shared::tam_next_read_alignments($reference_tam, $reference_header, $last_reference_alignment);		
		
		my $f = 0;
		READ: while (my $seq = $in_fastq[$f]->next_seq)
		{			
			$i++;
			last if ($settings->{alignment_read_limit} && ($i > $settings->{alignment_read_limit}));
			print STDERR "    READS:$i\n" if ($i % 10000 == 0);
			
			#Does this read have candidate junction matches?
			my $this_candidate_junction_al = [];
			if (($candidate_junction_al) && ($candidate_junction_al->[0]->qname =~ m/^$seq->{id}/))
			{
				$this_candidate_junction_al = $candidate_junction_al;
				($candidate_junction_al, $last_candidate_junction_alignment) 
					= Breseq::Shared::tam_next_read_alignments($candidate_junction_tam, $candidate_junction_header, $last_candidate_junction_alignment);
			}
			
			#Does this read have reference sequence matches?
			my $this_reference_al = [];
			if (($reference_al) && ($reference_al->[0]->qname =~ m/^$seq->{id}/))
			{
				$this_reference_al = $reference_al;
				($reference_al, $last_reference_alignment) 
					= Breseq::Shared::tam_next_read_alignments($reference_tam, $reference_header, $last_reference_alignment);
			}		
			
			## Nothing to be done if there were no matches to either
			if ((@$this_candidate_junction_al == 0) && (@$this_reference_al == 0))
			{
				$s->{unmatched_reads}++;
				## if not, then write to unmatched read file
				if ($settings->{unmatched_reads}) 
				{
					$out_unmatched_fastq[$f]->write_seq($seq);
				}
				next READ;
			}
						
			###			
			## Matches to candidate junctions may not overlap the junction.
			##
			## Reduce this list to those that overlap ANY PART of the junction.
			## Keep the ones that only match through the overlap, but do not propagate
			## to the opposite side, separate because they are only additional
			## evidence for predicted junctions and NOT support for the new junction
			## on their own. (They will also match the original reference genome equally well).
			###

			my $this_overlap_only_al = []; #reads that just map to the overlap region of a junction
			
			for (my $i=0; $i<scalar @$this_candidate_junction_al; $i++)
			{
				my $a = $this_candidate_junction_al->[$i];
				next if ($a->unmapped);
				
				my $junction_id = $candidate_junction_header->target_name()->[$a->tid];
				my $scj = Breseq::Shared::junction_name_split($junction_id);
				my $overlap = $scj->{overlap};

				## find the start and end coordinates of the overlap
				my ($junction_start, $junction_end);
				
				$junction_start = $flanking_length + 1;
				$junction_end = $flanking_length + abs($overlap);

				## If it didn't overlap the junction at all, remove it
				if ( ($a->start > $junction_end) || ($a->end < $junction_start) )
				{
					splice  @$this_candidate_junction_al, $i, 1;
					$i--;
				}
				else
				{	
					## If it overlaps ONLY the overlap part of the junction
					## then keep it separate in a special list...
					## Note that this only applies if overlap is negative
					## (if positive then the overlap part is unique sequence)
					if ( ($overlap > 0) && ($a->end <= $junction_end) || ($a->start >= $junction_start) )
					{						
						my $overlap_only_match = splice  @$this_candidate_junction_al, $i, 1;
						push @$this_overlap_only_al, $overlap_only_match;
						$i--;
					}
				}		
			}

			###			
			## Determine if the read has a better match to a candidate junction
			## or to the reference sequence.
			###

			### There are three possible kinds of reads at this point
			##
			## 1: Read has a best match to the reference genome
			## 2: Read has a best match (or multiple best matches) to junctions
			## 3: Read has an equivalent match to the reference genome
			##      and goes into the overlap part of a junction condidate
			###
			
			my $best_candidate_junction_score = -9999;
			my $best_reference_score = -9999;

			if (@$this_candidate_junction_al)
			{
				my $ca = $this_candidate_junction_al->[0];
				$best_candidate_junction_score = $ca->aux_get("AS");
			}
			elsif (@$this_overlap_only_al)
			{
				my $ca = $this_overlap_only_al->[0];
				$best_candidate_junction_score = $ca->aux_get("AS");
			}
			
			if (@$this_reference_al)
			{
				my $ra = $this_reference_al->[0];
				$best_reference_score = $ra->aux_get("AS");
			}


			### The best match we found to the reference was no better than the best to the
			### candidate junction. This read potentially supports the candidate junction.
			###
			### ONLY allow EQUAL matches through if they match the overlap only, otherwise
			### you can get predictions of new junctions with all reads supporting them
			### actually mapping perfectly to the reference.
			my $best_match_to_reference = ($best_candidate_junction_score < $best_reference_score);			
			if  (!$best_match_to_reference)
			{	
				my $equal_match_to_reference = ($best_candidate_junction_score == $best_reference_score);				
				
				## We flag whether the best hit was overlap only.
				my $best_candidate_junction_is_overlap_only = 0;
				my @this_dominant_candidate_junction_al;
				
				## only accept equal matches if they are to the overlap only 
				## -- preventing them from adding to the score in cases where a normal read match
				## to certain parts of the reference also creates junctions 
				if (@$this_candidate_junction_al && !$equal_match_to_reference)
				{
					@this_dominant_candidate_junction_al = _alignment_list_to_dominant_best($minimum_best_score, undef, @$this_candidate_junction_al);
				}
				elsif (@$this_overlap_only_al)
				{
					$best_candidate_junction_is_overlap_only = 1;
					@this_dominant_candidate_junction_al = _alignment_list_to_dominant_best($minimum_best_score, undef, @$this_overlap_only_al);
				}
				
				$best_match_to_reference = 1 if (scalar @this_dominant_candidate_junction_al == 0);
				
				if (!$best_match_to_reference)
				{
				
					my $item = {
						reference_alignments => $this_reference_al, 					# reference sequence alignments
						dominant_alignments => \@this_dominant_candidate_junction_al,   #the BEST candidate junction alignments
						dominant_alignment_is_overlap_only => $best_candidate_junction_is_overlap_only,
						fastq_file_index => $fastq_file_index[$f],						#index of the fastq file this read came from
					};
					#print Dumper($item);
	
					## Just one best hit, we put this in the hash that is used to predict junctions first
					if (scalar @this_dominant_candidate_junction_al == 1)
					{
						my $a = $this_dominant_candidate_junction_al[0];
						my $junction_id = $candidate_junction_header->target_name()->[$a->tid];
						#print "$junction_id\n";
						push @{$matched_junction{$junction_id}}, $item;
					}					
					## Multiple equivalent matches to junctions, ones with most hits later will win these matches
					## these matches are added BEFORE scoring a junction prediction
					else 
					{					
						foreach my $a (@this_dominant_candidate_junction_al)
						{	
							my $junction_id = $candidate_junction_header->target_name()->[$a->tid];
							my $read_name = $a->qname;
							$degenerate_matches{$junction_id}->{$read_name} = $item;
						}
					}
				}
			}
			
			### best match is to the reference, record in that SAM file.
			if ($best_match_to_reference && (scalar @$this_reference_al > 0))
			{
				_write_reference_matches($minimum_best_score, $minimum_best_score_difference, $reference_fai, $ref_seq_info, $RREF, $reference_header, $fastq_file_index[$f], @$this_reference_al);
			}

		} continue {
			$f++;
			$f %= scalar @in_fastq;
		}
	
		## save statistics
		$summary->{alignment_correction}->{$read_file} = $s;
	}	
	
	###			
	## Determine which junctions are real, prefer ones with most matches
	###
	
	my @rejected_keys = ();
	my %junction_test_info;
	
	## This keeps track of what overlap only reads we have already written as reference reads so
	## they don't appear multiple times. Better record keeping might obviate the need for this.
	my %written_overlap_only_reads;
	
	## first deal with ones with unique matches
	my @sorted_keys = sort {-(scalar @{$matched_junction{$a}} <=> scalar @{$matched_junction{$b}})} keys %matched_junction;
			
	my @new_keys;
	foreach my $key (@sorted_keys)
	{
		my ($failed, $has_non_overlap_only) = _test_junction($key, \%matched_junction, \%degenerate_matches, \%junction_test_info, $minimum_best_score, $minimum_best_score_difference, $reference_fai, $ref_seq_info, $RREF, $reference_header, $RCJ, $candidate_junction_header, \%written_overlap_only_reads);

		if (!$failed && !$has_non_overlap_only)
		{
			push @new_keys, $key; 
		}
#		else
		elsif ($has_non_overlap_only)
		{
			push @rejected_keys, $key;
		}
	}
	
	
	## next deal with ones with degenerate matches only	
	@sorted_keys = sort {-(scalar keys %{$degenerate_matches{$a}} <=> scalar keys %{$degenerate_matches{$b}})} keys %degenerate_matches;
	while (@sorted_keys)
	{
		my $key = shift @sorted_keys;
		
		print "Trying degenerate $key...\n" if ($verbose);
		
		next if (!defined $degenerate_matches{$key}); #they can be removed 
		
		my ($failed, $has_non_overlap_only) = _test_junction($key, \%matched_junction, \%degenerate_matches, \%junction_test_info, $minimum_best_score, $minimum_best_score_difference, $reference_fai, $ref_seq_info, $RREF, $reference_header, $RCJ, $candidate_junction_header);
		if (!$failed && !$has_non_overlap_only)
		{
			@sorted_keys = sort {-(scalar keys %{$degenerate_matches{$a}} <=> scalar keys %{$degenerate_matches{$b}})} keys %degenerate_matches;
			push @new_keys, $key;
		}
#		else
		elsif ($has_non_overlap_only)
		{
			push @rejected_keys, $key;
		}
	}
	@sorted_keys = @new_keys;
	
	#print successful ones out
	print "Successful hybrids\n" if ($verbose);

	#re-sort since some have gained reads from degenerate matches
	@sorted_keys = sort {-(scalar @{$matched_junction{$a}} <=> scalar @{$matched_junction{$b}})} @sorted_keys;
	@rejected_keys = sort {-(scalar @{$matched_junction{$a}} <=> scalar @{$matched_junction{$b}})} @rejected_keys;

	my @hybrid_predictions;
 	foreach my $key (@sorted_keys)
 	{	  
		print "$key\n" if ($verbose);
		my $item = _junction_to_hybrid_list_item($key, $ref_seq_info, scalar @{$matched_junction{$key}}, $junction_test_info{$key});	
		push @hybrid_predictions, $item;
		
		## create matches from UNIQUE sides of each match to reference genome
		## this fixes, for example appearing to not have any coverage at the origin of a circular DNA fragment
		### currently we do not add coverage to the IS element (which we would have to know all copies of to do...)
	
		if ( $settings->{add_split_junction_sides} )
		{
			foreach my $match (@{$matched_junction{$key}})
			{
				my $a = $match->{dominant_alignments}->[0];
				my $fastq_file_index = $match->{fastq_file_index};
				foreach my $side (1,2)
				{
					my $side_key = 'interval_' . $side;
	
                    ## Do not count for coverage if it is redundant!!
					if (!$item->{$side_key}->{redundant})
					{	
						##write out match corresponding to this part to SAM file
						
						##NOTE: this will not trim reads that extend to very near the end of the CJ sequence correctly!						
						my $trim = _trim_ambiguous_ends($a, $candidate_junction_header, $candidate_junction_fai);						
						Breseq::Shared::tam_write_moved_alignment(
							$RREF, 
							$a, 
							$fastq_file_index, 
							$item->{$side_key}->{seq_id}, 
							$item->{$side_key}->{start}, 
							$item->{$side_key}->{strand},
							$item->{$side_key}->{overlap}, 
							$side, 
							$item->{flanking_left}, 
							$item->{alignment_overlap}, 
							$trim
						);		
					}
				}
			}
		}
	}
	
	my @rejected_hybrid_predictions = ();
		
	foreach my $key (@rejected_keys)
	{
		print "$key\n" if ($verbose);
		my $item = _junction_to_hybrid_list_item($key, $ref_seq_info, scalar @{$matched_junction{$key}}, $junction_test_info{$key});
		$item->{marginal}=1;
		push @rejected_hybrid_predictions, $item;
	}
	push @hybrid_predictions, @rejected_hybrid_predictions;

 	return @hybrid_predictions;	
}

sub _write_reference_matches
{
	my ($minimum_best_score, $minimum_best_score_difference, $reference_fai, $ref_seq_info, $RREF, $reference_header, $fastq_file_index, @reference_al) = @_;
	
	@reference_al = _alignment_list_to_dominant_best($minimum_best_score, $minimum_best_score_difference, @reference_al);
	
	my @trims;
	foreach my $a (@reference_al)
	{
	#	push @trims, _trim_ambiguous_ends($a, $reference_header, $reference_fai);
		push @trims, _trim_ambiguous_ends($a, $reference_header, $reference_fai, $ref_seq_info); #slightly faster than using fai	
	}
	
	Breseq::Shared::tam_write_read_alignments($RREF, $reference_header, $fastq_file_index, \@reference_al, \@trims);
}

sub _alignment_begins_with_match_in_read
{
	my ($a) = @_;
	return 0 if ($a->unmapped);
	
	my $ca = $a->cigar_array;
	return ($ca->[-1]->[0] eq 'M') if ($a->reversed);
	return ($ca->[0]->[0] eq 'M');
}

sub _alignment_length_on_query
{
	my ($a) = @_;
	
	return 0 if ($a->unmapped);
	
	my $ca = $a->cigar_array;	
	my $start = 1;
	$start += $ca->[0]->[1] if ($ca->[0]->[0] eq 'S');
	my $end = $a->query->length;
	$end -= $ca->[-1]->[1] if ($ca->[-1]->[0] eq 'S');
	
	($start, $end) = ($a->query->length - $start + 1, $a->query->length - $end + 1) if ($a->reversed);
	($start, $end) = ($end, $start) if ($a->reversed);
	
	return ($end-$start+1);
}

sub _alignment_list_to_dominant_best
{
	###MOVE TO SETTINGS
	my $minimum_match_length = 28; 
	
	my ($minimum_best_score, $minimum_best_score_difference, @al) = @_;
	return () if (scalar @al <= 0);
	
	## require a minimum length of the read to be mapped
 	@al = grep {_alignment_length_on_query($_) >= $minimum_match_length} @al;
	return () if (scalar @al == 0);
	
	my $best_score = $al[0]->aux_get("AS");
	
	## no scores meet minimum
	return () if (defined $minimum_best_score && ($best_score < $minimum_best_score));
	
	## how many reads share the best score?
	my $last_best = 0;
	while (($last_best+1 < scalar @al) && ($al[$last_best+1]->aux_get("AS") == $best_score))
	{
		$last_best++;
	}
	
	## no scores meet minimum difference between best and next best
	if (defined $minimum_best_score_difference && (scalar @al > $last_best+1))
	{
		my $second_best_score = $al[$last_best+1]->aux_get("AS");
		return () if ($second_best_score + $minimum_best_score_difference >= $best_score)
	}
	
	return splice @al, 0, $last_best+1;
}


sub _trim_ambiguous_ends
{
	my $verbose = 0;
	my ($a, $header, $fai, $ref_seq_info) = @_;
	
#	print $a->qname . "\n";	
#	$verbose = ($a->qname eq '20AOWAAXX-Lenski:6:6:817:62');
	
	# has two keys: 'left' and 'right' which are how far to inset in REFERENCE
	my $trims;
	
	#which reference sequence?
	my $seq_id = $header->target_name->[$a->tid];
	my ($ref_strings, $ref_seq_length);
	
	## using $fai is more compatible, and must currently be used for junctions
	## using $ref_seq_info is slightly quicker, and currently used for the reference sequence
	if (defined $ref_seq_info)
	{		
		$ref_strings = $ref_seq_info->{ref_strings};	
		$ref_seq_length = length($ref_strings->{$seq_id});
	}
	## >>> transition to not using ref_seq_info
	else
	{
		$ref_seq_length  = $header->target_len->[$a->tid];
	}
	
	#create sequence snippets that we need to pay attention to ends of sequence
	my $expand_by = 36;
	my $expand_left = ($a->start-1 < $expand_by) ? $a->start-1 : $expand_by;
	my $expand_right = ($ref_seq_length - $a->end < $expand_by) ? $ref_seq_length-$a->end : $expand_by;
	
	
	my $expanded_ref_string = '';	
	if (defined $ref_strings)
	{
		$expanded_ref_string = substr $ref_strings->{$seq_id}, $a->start-$expand_left-1, ($a->end+$expand_right) - ($a->start-$expand_left) + 1;
	}
	## >>> transition to not using ref_seq_info
	else
	{
		my $expanded_ref_range = $seq_id . ':' . ($a->start-$expand_left) . '-' . ($a->end+$expand_right);
		$expanded_ref_string = $fai->fetch($expanded_ref_range);		
	}

	my $ref_string;
	if (defined $ref_strings)
	{
		$ref_string = substr $ref_strings->{$seq_id}, $a->start-1, $a->end - $a->start + 1;
	}	
	## >>> transition to not using ref_seq_info	
	else
	{
		my $ref_range = $seq_id . ':' . $a->start . '..' . $a->end;
		$ref_string = $fai->fetch( $seq_id . ':' . $a->start . '-' . $a->end );
	}
		
	my ($q_start, $q_end) = Breseq::Shared::alignment_query_start_end($a, { 'no_reverse' => 1} );
	my $q_length = $a->l_qseq;
	my $qry_string = substr $a->qseq, $q_start-1, $q_end - $q_start + 1;
	my $full_qry_string = $a->qseq;
	
	#take maximum of inset for query and reference
	my ($left_ref_inset, $right_ref_inset) = _ambiguous_end_offsets_from_sequence($ref_string);
	
	#add UNALIGNED bases at te end of reads
	$left_ref_inset += $q_start - 1;
	$right_ref_inset += $q_length - $q_end;
	
	
	## save a little time if qry and ref sequences are identical.
	my ($left_qry_inset, $right_qry_inset) = (0,0);
	if ($ref_string ne $qry_string)
	{
		($left_qry_inset, $right_qry_inset) = _ambiguous_end_offsets_from_sequence($qry_string);
		#add UNALIGNED bases at te end of reads
		$left_qry_inset += $q_start - 1;
		$right_qry_inset += $q_length - $q_end;
		
	}
	
	my ($left_full_qry_inset, $right_full_qry_inset) = _ambiguous_end_offsets_from_sequence($full_qry_string);
	my ($left_ref_expanded_inset, $right_ref_expanded_inset) 
		= _ambiguous_end_offsets_from_expanded_sequence($expand_left, $expand_right, $expanded_ref_string);
	
	if ($verbose)
	{
		print "Whole Read: $ref_string\n";
		print "Qry Start, End: $q_start, $q_end\n";	
		print "Ref: $ref_string\n";
		print "Ref insets: $left_ref_inset, $right_ref_inset\n";
		print "Qry: $qry_string\n";
		print "Qry insets: $left_qry_inset, $right_qry_inset\n";
		print "Full Qry: $full_qry_string\n";
		print "Full Qry insets: $left_full_qry_inset, $right_full_qry_inset\n";
		print "Expanded: $expanded_ref_string\n";
		print "Expand: $expand_left, $expand_right\n";
		print "Expanded Ref insets: $left_ref_expanded_inset, $right_ref_expanded_inset\n";		
	}


	$left_qry_inset = ($left_qry_inset > $left_full_qry_inset) ? $left_qry_inset : $left_full_qry_inset;
	$right_qry_inset = ($right_qry_inset > $right_full_qry_inset) ? $right_qry_inset : $right_full_qry_inset;

	$left_ref_inset = ($left_ref_inset > $left_ref_expanded_inset) ? $left_ref_inset : $left_ref_expanded_inset;
	$right_ref_inset = ($right_ref_inset > $right_ref_expanded_inset) ? $right_ref_inset : $right_ref_expanded_inset;

	if ($verbose)
	{
		print "Ref insets: $left_ref_inset, $right_ref_inset\n";
		print "Full Qry insets: $left_full_qry_inset, $right_full_qry_inset\n";
		print "Expanded Ref insets: $left_ref_expanded_inset, $right_ref_expanded_inset\n";		
	}

	##
	# Correct insets in the ref sequence to read coordinates, which must pay attention to gaps
	##
	
	#if no gaps, then just need to take the greater one
	my $left_ref_count = ($left_ref_inset > $left_qry_inset) ? $left_ref_inset : $left_qry_inset;
	my $left_qry_count = $left_ref_count;
	my $right_ref_count = ($right_ref_inset > $right_qry_inset) ? $right_ref_inset : $right_qry_inset;
	my $right_qry_count = $right_ref_count;

	if ($verbose)
	{
		print "Qry Count: $left_qry_count, $right_qry_count\n";
	}
	
		
	# my $cigar = $a->cigar_array;
	# ##skip soft padding
	# shift @$cigar if ($cigar->[0]->[0] eq 'S');
	# 
	# my $ungapped_left_qry_inset = $left_ref_inset;
	# foreach (my $i=0; $i<$left_ref_count; $i++)
	# {
	# 	$i;
	# }

	return ( {'L'=>$left_qry_count, 'R'=>$right_qry_count} );
}

sub _ambiguous_end_offsets_from_expanded_sequence
{
	my ($expand_left, $expand_right, $ref_string) = @_;
	
	my $left_inset = 0;
	my $right_inset = 0;
	
	{ #left side

		#maximum size to check is $expand_by
		my $test_left_inset = 0;
		while ($test_left_inset < $expand_left)
		{
			my $found_left_inset = $test_left_inset;
			my $test_length = $expand_left-$test_left_inset;

			my $match_found = 0;
			my $test_end_string;
			while ( ($test_length > 0) && !$match_found)
			{	
				$test_end_string = substr $ref_string, $found_left_inset, $test_length;
				#force removal of end nucleotides if no higher order repeats found
				$found_left_inset-- if ($test_length == 1);			
				my $test_interior_string = substr $ref_string, $found_left_inset+$test_length, $test_length;
				
				while ($test_end_string eq $test_interior_string)
				{
					$match_found = 1;
					$found_left_inset += $test_length;
					$test_interior_string = substr $ref_string, $found_left_inset+$test_length, $test_length;
				}
				$found_left_inset += $test_length if ($match_found);
				$test_length--;
			} 
				
			#test partial matches (that continue part of repeat further)
			 #note: already starts at one less than actual repeat size
			my $test_partial_size = $test_length; 
			while ($test_partial_size > 0)
			{
				my $test_partial_end_string = substr $test_end_string, 0, $test_partial_size;
				my $test_interior_string = substr $ref_string, $found_left_inset, $test_partial_size;
				if ($test_partial_end_string eq $test_interior_string)
				{
					$found_left_inset += $test_partial_size;
					last;
				}
				$test_partial_size--;
			}
		
			$left_inset = $found_left_inset-$expand_left if ($found_left_inset-$expand_left > $left_inset);
			$test_left_inset++;
		}
	}
	
	{ #right side

		#maximum size to check is $expand_by
		my $test_right_inset = 0;
		while ($test_right_inset < $expand_right)
		{
			my $found_right_inset = $test_right_inset;
			my $test_length = $expand_right-$test_right_inset;

			my $match_found = 0;
			my $test_end_string;
			while ( ($test_length > 0) && !$match_found)
			{	
				$test_end_string = substr $ref_string, $found_right_inset, $test_length;
				#force removal of end nucleotides if no higher order repeats found
				$found_right_inset-- if ($test_length == 1);			
				my $test_interior_string = substr $ref_string, $found_right_inset+$test_length, $test_length;
				
				while ($test_end_string eq $test_interior_string)
				{
					$match_found = 1;
					$found_right_inset += $test_length;
					$test_interior_string = substr $ref_string, $found_right_inset+$test_length, $test_length;
				}
				$found_right_inset += $test_length if ($match_found);
				$test_length--;
			} 
				
			#test partial matches (that continue part of repeat further)
			 #note: already starts at one less than actual repeat size
			my $test_partial_size = $test_length; 
			while ($test_partial_size > 0)
			{
				my $test_partial_end_string = substr $test_end_string, 0, $test_partial_size;
				my $test_interior_string = substr $ref_string, $found_right_inset, $test_partial_size;
				if ($test_partial_end_string eq $test_interior_string)
				{
					$found_right_inset += $test_partial_size;
					last;
				}
				$test_partial_size--;
			}
		
			$right_inset = $found_right_inset-$expand_right if ($found_right_inset-$expand_right > $right_inset);
			$test_right_inset++;
		}
	}

	
	
	return ($left_inset, $right_inset);
}


sub _ambiguous_end_offsets_from_sequence
{
	my ($ref_string) = @_;
	
	#test longest substrings
	my $left_inset = 0;
	my $right_inset = 0;
	
	{ #left side
		my $test_length = int((length $ref_string) / 2);
		my $match_found = 0;
		my $test_end_string;
		while ( ($test_length > 0) && !$match_found)
		{	
			$test_end_string = substr $ref_string, $left_inset, $test_length;
			#force removal of end nucleotides if no higher order repeats found
			$left_inset-- if ($test_length == 1);			
			my $test_interior_string = substr $ref_string, $left_inset+$test_length, $test_length;
			
			while ($test_end_string eq $test_interior_string)
			{
				$match_found = 1;
				$left_inset += $test_length;
				$test_interior_string = substr $ref_string, $left_inset+$test_length, $test_length;
			}
			$left_inset += $test_length if ($match_found);
			$test_length--;
		} 
		
		#test partial matches (that continue part of repeat further)
		 #note: already starts at one less than actual repeat size
		my $test_partial_size = $test_length; 
		while ($test_partial_size > 0)
		{
			my $test_partial_end_string = substr $test_end_string, 0, $test_partial_size;
			my $test_interior_string = substr $ref_string, $left_inset, $test_partial_size;
			if ($test_partial_end_string eq $test_interior_string)
			{
				$left_inset += $test_partial_size;
				last;
			}
			$test_partial_size--;
		}
		
	}
	
	{ #right side
		my $test_length = int((length $ref_string) / 2);
		my $match_found = 0;
		my $test_end_string;
		while ( ($test_length > 0) && !$match_found)
		{	
			$test_end_string = substr $ref_string, -$right_inset-$test_length, $test_length;
			#force removal of end nucleotides if no higher order repeats found
			$right_inset-- if ($test_length == 1);			
			my $test_interior_string = substr $ref_string, -$right_inset-2 * $test_length, $test_length;
						
			while ($test_end_string eq $test_interior_string)
			{
				$match_found = 1;
				$right_inset += $test_length;
				$test_interior_string = substr $ref_string, -$right_inset-2 * $test_length, $test_length;
			}
			$right_inset += $test_length if ($match_found);
			
			$test_length--;
		} 
		
		#test partial matches (that continue part of repeat further)
		 #note: already starts at one less than actual repeat size
		my $test_partial_size = $test_length; 
		while ($test_partial_size > 0)
		{
			my $test_partial_end_string = substr $test_end_string, -$test_partial_size, $test_partial_size;
			my $test_interior_string = substr $ref_string, -$right_inset-$test_partial_size, $test_partial_size;
			if ($test_partial_end_string eq $test_interior_string)
			{
				$right_inset += $test_partial_size;
				last;
			}
			$test_partial_size--;
		}
	}
	
	return ($left_inset, $right_inset);
}


sub _test_junction
{
	my ($key, $matched_junction_ref, $degenerate_matches_ref, $junction_test_info_ref, $minimum_best_score, $minimum_best_score_difference, $reference_fai, $ref_seq_info, $RREF, $reference_header, $RCJ, $candidate_junction_header, $written_overlap_only_reads) = @_;
	
	my $test_info;
	my @unique_matches = ();
	@unique_matches = @{$matched_junction_ref->{$key}} if (defined $matched_junction_ref->{$key});
	my @degenerate_matches = ();
	@degenerate_matches = map { $degenerate_matches_ref->{$key}->{$_} } sort keys %{$degenerate_matches_ref->{$key}} if (defined $degenerate_matches_ref->{$key});		
		
	my $failed = 0;
#	print "Junction Candidate: $key Unique Matches: " . (scalar @unique_matches) . " Degenerate Matches: " . (scalar @degenerate_matches) . "\n";

	#### TEST 1: There should be a minimal number of reads supporting the junction
#	my $minimum_number_of_reads_for_junction = 3;	
#	$failed = 1 if (scalar @unique_matches < $minimum_number_of_reads_for_junction);

	#### TEST 2: There should be at least one read that goes a certain number of bp into the nonoverlap sequence on each side of the junction

	my $max_left_per_strand = { '0'=> 0, '1'=>0 };
	my $max_right_per_strand = { '0'=> 0, '1'=>0 };
	my $max_min_left_per_strand = { '0'=> 0, '1'=>0 };
	my $max_min_right_per_strand = { '0'=> 0, '1'=>0 };	
	my $count_per_strand = { '0'=> 0, '1'=>0 };
	my $total_non_overlap_reads = 0;
	my $score = 0;
	
	## is there at least one read that isn't overlap only?
	## displaying ones where it doesn't as marginals looks really confusing
	my $has_non_overlap_only = 1; 
	
	
	### we also need to count degenerate matches b/c sometimes ambiguity unfairly penalizes real reads...
	READ: foreach my $item (@unique_matches, @degenerate_matches)
	{
		## we don't want to count matches that do not extend through the overlap region
		next READ if ($item->{dominant_alignment_is_overlap_only});
		$total_non_overlap_reads++;
		$has_non_overlap_only = 0;
		
		#If there were no degenerate matches, then we could just take the
		#one and only match in the 'dominant_alignments' array
		#my $a = $item->{dominant_alignments}->[0]; 
		
		## as it is, we must be sure we are looking at the one that matches
		my $a;
		DOMINANT_ALIGNMENT: foreach my $candidate_a (@{$item->{dominant_alignments}})
		{
			my $junction_id = $candidate_junction_header->target_name()->[$candidate_a->tid];
			if ($key eq $junction_id)
			{
				$a = $candidate_a;
				last DOMINANT_ALIGNMENT;
			}
		}
		die if (!defined $a);
		
		my $junction_id = $candidate_junction_header->target_name()->[$a->tid];
		my $scj = Breseq::Shared::junction_name_split($junction_id);
		my $overlap = $scj->{overlap};
		my $flanking_left = $scj->{flanking_left};

		my $rev_key = ($a->reversed ? '1' : '0');

		$count_per_strand->{$rev_key}++;

		##The left side goes exactly up to the flanking length
		my $this_left = $flanking_left;
		$this_left = $this_left - $a->start+1;

		#The right side starts after moving past any overlap (negative or positive)
		my $this_right = $flanking_left+1;
		$this_right += abs($overlap);
		$this_right = $a->end - $this_right+1;

#		print "  " . $a->start . "-" . $a->end . " " . $overlap . " " . $rev_key . "\n";
#		print "  " . $item->{alignments}->[0]->qname . " LEFT: $this_left RIGHT: $this_right\n";

		## Update:
		### Score = the minimum unique match length on a side
		### Max_Min = the maximum of the minimum length match sides
		### Max = the maximum match on a side
		### Note that the max and min filtering is really a kind of poor man's KS test
		###   if we implemented that with a certain coverage cutoff it would be a
		###   more principled way of doing things... 
		if ($this_left < $this_right)
		{
			$score += $this_left;
			$max_min_left_per_strand->{$rev_key} = $this_left if ($max_min_left_per_strand->{$rev_key} < $this_left);
			
		}
		else
		{
			$score += $this_right;
			$max_min_right_per_strand->{$rev_key} = $this_right if ($max_min_right_per_strand->{$rev_key} < $this_right);
		}
	
		$max_left_per_strand->{$rev_key} = $this_left if ($max_left_per_strand->{$rev_key} < $this_left);
		$max_right_per_strand->{$rev_key} = $this_right if ($max_right_per_strand->{$rev_key} < $this_right);
		
	}				
	
	my $max_left = ($max_left_per_strand->{'0'} > $max_left_per_strand->{'1'}) ? $max_left_per_strand->{'0'} : $max_left_per_strand->{'1'};
	my $max_right = ($max_right_per_strand->{'0'} > $max_right_per_strand->{'1'}) ? $max_right_per_strand->{'0'} : $max_right_per_strand->{'1'};

	my $max_min_left = ($max_min_left_per_strand->{'0'} > $max_min_left_per_strand->{'1'}) ? $max_min_left_per_strand->{'0'} : $max_min_left_per_strand->{'1'};
	my $max_min_right = ($max_min_right_per_strand->{'0'} > $max_min_right_per_strand->{'1'}) ? $max_min_right_per_strand->{'0'} : $max_min_right_per_strand->{'1'};
	
	use Math::CDF;	
	## it is possible for the counts on each strand to be zero because all matches are overlap only
	my $strand_p_value = Math::CDF::pbinom($count_per_strand->{0}, $count_per_strand->{0}+$count_per_strand->{1}, 0.5);
	if (defined $strand_p_value)
	{
		$strand_p_value = 1-$strand_p_value if ($strand_p_value > 0.5);
		$strand_p_value = sprintf "%.1e", $strand_p_value; #round immediately
	}
	else
	{
		$strand_p_value = 'NA';
	}
	
	$test_info = {
		max_left => $max_left,
		max_left_minus => $max_left_per_strand->{0},
		max_left_plus => $max_left_per_strand->{1},
		max_right => $max_right,
		max_right_minus => $max_right_per_strand->{0},
		max_right_plus =>$max_right_per_strand->{1},
		max_min_right => $max_min_right,
		max_min_right_minus => $max_min_right_per_strand->{0},
		max_min_right_plus =>$max_min_right_per_strand->{1},		
		max_min_left => $max_min_left,
		max_min_left_minus => $max_min_left_per_strand->{0},
		max_min_left_plus =>$max_min_left_per_strand->{1},		
		coverage_minus => $count_per_strand->{0},
		coverage_plus => $count_per_strand->{1},
		strand_p_value => $strand_p_value,
		total_non_overlap_reads => $total_non_overlap_reads,
		score => $score,
	};
		
	## These parameters still need additional testing
	my $alignment_on_each_side_cutoff = 16; #14
	my $alignment_on_each_side_cutoff_per_strand = 13; #9
	my $alignment_on_each_side_min_cutoff = 5;
		
	$failed = 	   ($max_left < $alignment_on_each_side_cutoff) 
				|| ($max_right < $alignment_on_each_side_cutoff)
	       		|| ($max_left_per_strand->{'0'} < $alignment_on_each_side_cutoff_per_strand) 
				|| ($max_left_per_strand->{'1'} < $alignment_on_each_side_cutoff_per_strand)
	       		|| ($max_right_per_strand->{'0'} < $alignment_on_each_side_cutoff_per_strand) 
				|| ($max_right_per_strand->{'1'} < $alignment_on_each_side_cutoff_per_strand)
	       		|| ($max_right_per_strand->{'0'} < $alignment_on_each_side_cutoff_per_strand) 
				|| ($max_right_per_strand->{'1'} < $alignment_on_each_side_cutoff_per_strand)
				|| ($max_min_left < $alignment_on_each_side_min_cutoff)
				|| ($max_min_right < $alignment_on_each_side_min_cutoff)
	;

	#### TEST X: Overlap should not be biased such that one side of the junction often has more of the
	####         read overlapping it than the other. Use a sign test.
	####
	####   >>>>  But it's really not this simple since reads might be biased by the sequences they match???	
	
	## to implement
	
	
	### If we passed all the tests, or we were only testing degenerate junctions
	### add degenerate matches and make them unavailable for other junctions	
	### degenerate matches is a hash of junction_ids of read_names
	if (!$failed) # || !scalar (@unique_matches == 0))
	{
		if (defined $degenerate_matches_ref->{$key})
		{
			foreach my $read_name (keys %{$degenerate_matches_ref->{$key}})
			{
				my $degenerate_match = $degenerate_matches_ref->{$key}->{$read_name};
				my $matched_alignment;
			
				#purge all references to this from the degenerate match hash
				#so that they will not be counted for other junctions
				foreach my $a (@{$degenerate_match->{dominant_alignments}})
				{
					my $junction_id = $candidate_junction_header->target_name()->[$a->tid];
					$matched_alignment = $a if ($key eq $junction_id);
					delete $degenerate_matches_ref->{$junction_id}->{$read_name};
					
					if (scalar keys %{$degenerate_matches_ref->{$junction_id}} == 0)
					{
						delete $degenerate_matches_ref->{$junction_id};
					}
				}
			
				## add to the matched_junction
				@{$degenerate_match->{dominant_alignments}} = ($matched_alignment);
				push @{$matched_junction_ref->{$key}}, $degenerate_match;
			}
		}
	}	
	
	## Write out the matches to the proper SAM file(s) depending on whether the junction succeeded or failed
	foreach my $junction_read (@{$matched_junction_ref->{$key}})
	{
		my $fastq_file_index = $junction_read->{fastq_file_index};
		
		## Write matches to reference sequences if we failed
		if ($failed)
		{		
			my $this_reference_al = $junction_read->{reference_alignments};
			my $read_name;
			$read_name = $this_reference_al->[0]->qname if ((defined $this_reference_al) && (scalar @$this_reference_al > 0));
			
			## hashing by read name here makes sure that we only write each read once in the reference file...
			if (!$junction_read->{dominant_alignment_is_overlap_only} || ($read_name && !$written_overlap_only_reads->{$read_name}))
			{
				_write_reference_matches($minimum_best_score, $minimum_best_score_difference, $reference_fai, $ref_seq_info, $RREF, $reference_header, $fastq_file_index, @$this_reference_al);
				$written_overlap_only_reads->{$read_name} = 1;
			}
		}
		
		## REGARDLESS: write matches to the candidate junction SAM file 
		if (!$has_non_overlap_only) 
		{
			my @this_dominant_candidate_junction_al = @{$junction_read->{dominant_alignments}}; 
			Breseq::Shared::tam_write_read_alignments($RCJ, $candidate_junction_header, $fastq_file_index, \@this_dominant_candidate_junction_al);
			$RCJ->flush();
		}
	}
		
	$junction_test_info_ref->{$key} = $test_info;
	return ($failed, $has_non_overlap_only);
}

sub _junction_to_hybrid_list_item
{
	my ($key, $ref_seq_info, $total_reads, $test_info) = @_;
	
	## split the key to an item with information about the junction
	my $item = Breseq::Shared::junction_name_split($key);
	$item->{key} = $key;
	$item->{test_info} = $test_info;
	
	## overlap may be adjusted below... this messes up making the alignment
	## 'alignment_overlap' is the original one that applies to the candidate junction BAM file
	## 'overlap' is a version where overlap has been resolved if possible for adding sides of the
	##    alignment
	$item->{alignment_overlap} = $item->{overlap};			
				
	##
	## Create three intervals for making alignments and html output
	##    it would be nice to flag whether we think the original junction is still there
	##
 	
	## (1) The first alignment has information relative to the junction candidates
	if ($item->{alignment_overlap} == 0)
	{
		$item->{start} = $item->{flanking_left};
		$item->{end} = $item->{flanking_left}+1;			
		$item->{mark} = '/';
	}
	elsif ($item->{alignment_overlap} > 0)
	{
		$item->{start} = $item->{flanking_left}+1;
		$item->{end} = $item->{flanking_left}+$item->{alignment_overlap};
		$item->{mark} = '|';
	}
	else ## ($item->{overlap} < 0)
	{
		$item->{start} = $item->{flanking_left}+1;
		$item->{end} = $item->{flanking_left}-$item->{alignment_overlap};
		$item->{mark} = '*';
	}
	$item->{seq_id} = $key;
	
	## These are now automatically created from splitting the junction key
	## (2) Alignment for the reference sequence that disagrees with the junction #1.
	## (3) Alignment for the reference sequence that disagrees with the junction #2.

	$item->{total_reads} = $total_reads;
	
	## Correct for overlapping IS elements
	
	###
	# IS insertion overlap correction
	#
	# For these the coordinates may have been offset incorrectly initially (because both sides of the junction may look unique)
	# The goal is to offset through positive overlap to get as close as possible to the ends of the IS
	###

	foreach my $key ('interval_1', 'interval_2')
	{
		## Determine IS elements
		## Is it within an IS or near the boundary of an IS in the direction leading up to the junction?			
		if (my $is = Breseq::ReferenceSequence::find_closest_is_element($item->{$key}, $ref_seq_info->{is_lists}->{$item->{$key}->{seq_id}}, 200, $item->{$key}->{strand}))
		{
			$item->{$key}->{is}->{gene} = $is->{gene};
			$item->{$key}->{is}->{interval} = ($is->{strand} == +1) ? "$is->{start}-$is->{end}" : "$is->{end}-$is->{start}"; 
			$item->{$key}->{is}->{product} = $is->{product};
		}
	}
	
	## use of $j is historical due to a move of this part of the function from annotate_rearrangements
	my $j = $item;
	sub add_is_coords_from_interval
	{
		my ($c) = @_;
		return if (!defined $c->{is}); 
		
		my ($is_start, $is_end) = split /-/, $c->{is}->{interval};
		$c->{is}->{strand} = ($is_start < $is_end) ? +1 : -1; 
		$c->{is}->{start} = ($is_start < $is_end) ? $is_start : $is_end; 
		$c->{is}->{end} = ($is_start < $is_end) ? $is_end : $is_start;
	}
	
	add_is_coords_from_interval($j->{interval_1});
	add_is_coords_from_interval($j->{interval_2});
	
	$j->{interval_1}->{read_side} = -1;
	$j->{interval_2}->{read_side} = +1;
		
	## Determine which side of the junction is the IS and which is unique
	## these point to the correct initial interval...
	if (defined $j->{interval_1}->{is} && !defined $j->{interval_2}->{is})
	{
		if (abs($j->{interval_1}->{is}->{start} - $j->{interval_1}->{start}) <= 20)
		{
			$j->{is_interval} = $j->{interval_1};
			$j->{is_interval}->{is}->{side_key} = 'start';
		}
		elsif (abs($j->{interval_1}->{is}->{end} - $j->{interval_1}->{start}) <= 20 )
		{
			$j->{is_interval} = $j->{interval_1};
			$j->{is_interval}->{is}->{side_key} = 'end';
		}
		$j->{unique_interval} = $j->{interval_2};
	}
	
	if (!defined $j->{interval_1}->{is} && defined $j->{interval_2}->{is})
	{
		if (abs($j->{interval_2}->{is}->{start} - $j->{interval_2}->{start}) <= 20)
		{
			$j->{is_interval} = $j->{interval_2};
			$j->{is_interval}->{is}->{side_key} = 'start';
		}
		elsif (abs($j->{interval_2}->{is}->{end} - $j->{interval_2}->{start}) <= 20 )
		{
			$j->{is_interval} = $j->{interval_2};
			$j->{is_interval}->{is}->{side_key} = 'end';
		}
		$j->{unique_interval} = $j->{interval_1};
	}
	
	## both were IS! -- define as redundant here
	if (defined $j->{interval_1}->{is} && defined $j->{interval_2}->{is})
	{
		$j->{interval_1}->{redundant} = 1;
		$j->{interval_2}->{redundant} = 1;
	}
	
	#by default, overlap is included on both sides of the junction (possibly changed below)
	$item->{interval_1}->{overlap} = 0;
	$item->{interval_2}->{overlap} = 0;
		
	## Resolve redundant overlap
	if ($item->{overlap} > 0)
	{
		$item->{interval_1}->{overlap} = $item->{overlap};
		$item->{interval_2}->{overlap} = $item->{overlap};
		
		## If there was in IS, resolve overlap so it goes to the edge of the IS element
		if (defined $j->{is_interval})
		{			
			### first, adjust the repetitive sequence boundary to get as close to the IS as possible
			my $move_dist = $j->{is_interval}->{strand} * ($j->{is_interval}->{is}->{$j->{is_interval}->{is}->{side_key}} - $j->{is_interval}->{start});
			$move_dist = 0 if ($move_dist < 0);
			$move_dist = $j->{overlap} if ($move_dist > $j->{overlap});
			$j->{is_interval}->{start} += $j->{is_interval}->{strand} * $move_dist;
			$j->{overlap} -= $move_dist;
			$j->{is_interval}->{overlap} -= $move_dist;
			$j->{is_interval}->{end} = $j->{is_interval}->{start};
	
			### second, adjust the unique sequence side with any remaining overlap
			$j->{unique_interval}->{start} += $j->{unique_interval}->{strand} * $j->{overlap};	
			$j->{unique_interval}->{end} = $j->{unique_interval}->{start};
			$j->{unique_interval}->{overlap} -= $j->{overlap};
			
			$j->{overlap} = 0;
			$j->{is_interval}->{redundant} = 1;
		}
	
		### Should this be "elsif" or "if"?
		### If there is no IS element and both sides are unique,
		### then give overlap to first side. This gives proper support for junctions.
		### and ensures we don't count this coverage twice.
		elsif (!$item->{interval_1}->{redundant} && !$item->{interval_2}->{redundant})
		{
			my $strand_direction = ($item->{interval_2}->{strand} > 0) ? +1 : -1;
	
			$item->{interval_2}->{start} += $item->{overlap} * $strand_direction;
			$item->{interval_2}->{end} += $item->{overlap} * $strand_direction;
			$item->{interval_2}->{overlap} = 0;
			$item->{overlap} = 0;			
		}
		
		#print STDERR Dumper($item);
	}
	
	### Note: Other adjustments to overlap can happen at the later annotation stage
	### and they will not affect coverage for calling deletions or mutations
	### because they will be in redundantly matched sides of junctions
	### However, this can cause confusing alignments
	return $item;
}


return 1;
