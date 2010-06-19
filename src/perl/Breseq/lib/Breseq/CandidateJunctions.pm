###
# Pod Documentation
###

=head1 NAME

CandidateJunction.pm

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

package Breseq::CandidateJunction;
use strict;

require Exporter;
our @ISA = qw( Exporter );
our @EXPORT = qw( $join_string );

use Bio::DB::Sam;

use Breseq::Shared;
use Breseq::Fastq;

use Data::Dumper;

### constants

=head2 identify_candidate_junctions

 Title   : identify_candidate_junctions
 Usage   : identify_candidate_junctions( );
 Function:
 Returns : 
 
=cut

sub identify_candidate_junctions
{
	my $verbose = 0;
	our ($settings, $summary, $ref_seq_info) = @_;

	#set up files and local variables from settings
	my $ref_strings = $ref_seq_info->{ref_strings};
	
	my $reference_faidx_file_name = $settings->file_name('reference_faidx_file_name');
	my $fai = Bio::DB::Sam::Fai->load($settings->file_name('reference_fasta_file_name'));
	
	my $candidate_junction_fasta_file_name = $settings->file_name('candidate_junction_fasta_file_name');
	my $out = Bio::SeqIO->new(-file => ">$candidate_junction_fasta_file_name", -format => 'fasta');
		
	#hash by junction sequence concatenated with name of counts
	my $candidate_junctions = {};
	
	### summary data for this step
	my $hcs;
	
	my @read_files = $settings->read_files;	
	my $i = 0;

	BASE_FILE: foreach my $read_struct ($settings->read_structures)
	{	
		my $read_file = $read_struct->{base_name};	
		print STDERR "  READ FILE::$read_file\n";
		
		## Zero out summary information
		my $s;
		$s->{reads_with_unique_matches} = 0;
		$s->{reads_with_non_unique_matches} = 0;
		$s->{reads_with_redundant_matches} = 0;
		$s->{reads_with_possible_hybrid_matches} = 0;
		$s->{reads_with_more_multiple_matches} = 0;
		$s->{reads_with_only_unwanted_matches} = 0;		
		$s->{reads_with_best_unwanted_matches} = 0;		
#		$s->{reads_with_no_matches} = 0; #calculate as total minus number found
		$s->{reads_with_ambiguous_hybrids} = 0;
		
		my $reference_sam_file_name = $settings->file_name('reference_sam_file_name', {'#'=>$read_file});
		my $tam = Bio::DB::Tam->open($reference_sam_file_name) or die("Could not open reference same file: $reference_sam_file_name");		
		my $header = $tam->header_read2($reference_faidx_file_name) or die("Error reading reference fasta index file: $reference_faidx_file_name");		
		my $last_alignment;
		my $al;
		
		ALIGNMENT_LIST: while (1)
		{		
			($al, $last_alignment) = Breseq::Shared::tam_next_read_alignments($tam, $header, $last_alignment);
			last ALIGNMENT_LIST if (!$al);
					
			$i++;
			print STDERR "    ALIGNED READ:$i\n" if ($i % 10000 == 0);

			last ALIGNMENT_LIST if (defined $settings->{candidate_junction_read_limit} && ($i > $settings->{candidate_junction_read_limit}));
			_alignments_to_candidate_junctions($settings, $summary, $ref_seq_info, $candidate_junctions, $fai, $header, @$al);
		}	
		
		$hcs->{read_file}->{$read_file} = $s;
	}		
		
		
	###	
	## Combine hash into a list, one item for each unique sequence (also joining reverse complements)
	###		
	
	my @combined_candidate_junctions;
	my $handled_seq;
	my %ids_to_print;
	
	## We should really sort the candidate junctions on something (seq_id then position?) so that we get 
	## consistent results!!
		
	JUNCTION_SEQ: foreach my $junction_seq (keys %$candidate_junctions)
	{	
		next if ($handled_seq->{$junction_seq});
		
		my @combined_candidate_junction_list = ();
		my $rc_junction_seq = Breseq::Fastq::revcom($junction_seq);
				
		JUNCTION_ID: foreach my $junction_id (keys %{$candidate_junctions->{$junction_seq}} )
		{
			## add redundancy to the $junction_id
			my $cj = $candidate_junctions->{$junction_seq}->{$junction_id};
			$junction_id .= $Breseq::Shared::junction_name_separator . $cj->{r1} . $Breseq::Shared::junction_name_separator . $cj->{r2};
			push @combined_candidate_junction_list, { id=>$junction_id, score=>$cj->{score}, seq=>$junction_seq, rc_seq=>$rc_junction_seq };
		}
		$handled_seq->{$junction_seq}++;
		
		## add the reverse complement
		if (defined $candidate_junctions->{$rc_junction_seq})
		{
			JUNCTION_ID: foreach my $junction_id (keys %{$candidate_junctions->{$rc_junction_seq}} )
			{
				## add redundancy to the $junction_id
				my $cj = $candidate_junctions->{$junction_seq}->{$junction_id};
				
				### IS IT CORRECT TO REVERSE THE REDUNDANCY HERE???
				$junction_id .= $Breseq::Shared::junction_name_separator . $cj->{r2} . $Breseq::Shared::junction_name_separator . $cj->{r1};
				push @combined_candidate_junction_list, { id=>$junction_id, score=>$cj->{score}, seq=>$rc_junction_seq, rc_seq=>$junction_seq };
			}
			$handled_seq->{$rc_junction_seq}++;
		}
		@combined_candidate_junction_list = sort { -($a->{score} <=> $b->{score}) } @combined_candidate_junction_list;
		
		##save just the best one
		push @combined_candidate_junctions, $combined_candidate_junction_list[0];
		
		##make sure it isn't a duplicate junction id -- how is this even possible?
		die "Duplicate junction id is being printed: $combined_candidate_junction_list[0]->{id}" if ($ids_to_print{$combined_candidate_junction_list[0]->{id}});
		$ids_to_print{$combined_candidate_junction_list[0]->{id}}++;
	}
	
	@combined_candidate_junctions = sort {-($a->{score} <=> $b->{score}) || (length($a->{seq}) <=> length($b->{seq}))} @combined_candidate_junctions;
	print Dumper(@combined_candidate_junctions) if ($verbose);
	
	###
	## Limit the number of candidate junctions that we print by:
	##   (1) A maximum number of candidate junctions
	##   (2) A maximum length of the sequences in candidate junctions
	###
	
	print STDERR "  Taking top candidate junctions..." . "\n";
	
	## adding up the lengths might be too time-consuming to be worth it...
	my $total_cumulative_cj_length = 0;
	my $total_candidate_junction_number = scalar @combined_candidate_junctions;
	foreach my $c (@combined_candidate_junctions)
	{
		$total_cumulative_cj_length += length $c->{seq};
	}
		
	my @duplicate_sequences;
	my $cumulative_cj_length = 0;
	my $lowest_accepted_score = 'NA';
	
	## Right now we limit the candidate junctions to have a length no longer than the reference sequence.
	my $cj_length_limit = int($summary->{sequence_conversion}->{total_reference_sequence_length} * $settings->{maximum_candidate_junction_length_factor});
	my $maximum_candidate_junctions = $settings->{maximum_candidate_junctions};
	my $minimum_candidate_junctions = $settings->{minimum_candidate_junctions};

	print STDERR sprintf ("  Minimum number to keep: %7d \n", $minimum_candidate_junctions);
	print STDERR sprintf ("  Maximum number to keep: %7d \n", $maximum_candidate_junctions);
	print STDERR sprintf ("  Maximum length to keep: %7d bases\n", $cj_length_limit);
	
	print STDERR "    Initial: Number = " . $total_candidate_junction_number . ", Cumulative Length = " . $total_cumulative_cj_length . " bases\n";


	if ((defined $settings->{maximum_candidate_junctions}) && (@combined_candidate_junctions > 0))
	{	
		my @remaining_ids = ();
		my @list_in_waiting = ();
		my $add_cj_length = 0;
		my $num_duplicates = 0;
	
		my $i = 0;
		my $current_score = $combined_candidate_junctions[$i]->{score};
	
		## Check to make sure that adding everything from the last iteration doesn't put us over any limits...
		my $new_number = scalar(@remaining_ids) + scalar(@list_in_waiting);
		my $new_length = $cumulative_cj_length + $add_cj_length;
		while (	( $new_number <= $minimum_candidate_junctions ) || (($new_length <= $cj_length_limit) && ($new_number <= $maximum_candidate_junctions)) )
		{			
			## OK, add everything from the last iteration
			$cumulative_cj_length += $add_cj_length;
			push @remaining_ids, @list_in_waiting;
			$lowest_accepted_score = $current_score;
			
			## Zero out what we will add
			$add_cj_length = 0;
			@list_in_waiting = ();
			$num_duplicates = 0;
			
			## Check to make sure we haven't exhausted the list
			last if ($i >= scalar @combined_candidate_junctions);

			$current_score = $combined_candidate_junctions[$i]->{score};
			CAND: while (($i < scalar @combined_candidate_junctions) && ($combined_candidate_junctions[$i]->{score} == $current_score))
			{
				my $c = $combined_candidate_junctions[$i];
				foreach my $dup (@duplicate_sequences)
				{
					#is this a subsequence of one already in the list?
					if ( ($dup =~ m/$c->{seq}/) || ($dup =~ m/$c->{rc_seq}/) || ($c->{seq} =~ m/$dup/) || ($c->{rc_seq} =~ m/$dup/) ) 
					{
						$num_duplicates++;
						print "Dup $dup\n$c->{seq}\n" if ($verbose);
						next CAND;
					}
				}
				push @duplicate_sequences, $c->{seq};
				push @list_in_waiting, $c;
				$add_cj_length += length $c->{seq};
				
			} continue {
				$i++;
			}
			
			$new_number = scalar(@remaining_ids) + scalar(@list_in_waiting);
			$new_length = $cumulative_cj_length + $add_cj_length;
			
			print STDERR sprintf("      Testing Score %5d: Number = %4d, Length = %6d, Duplicate = %4d\n", $current_score, (scalar @list_in_waiting), $add_cj_length, $num_duplicates);
			
		}
		@combined_candidate_junctions = @remaining_ids;
	}
	
	my $accepted_candidate_junction_number = scalar @combined_candidate_junctions;
	print STDERR "    Accepted: Number = $accepted_candidate_junction_number, Score >= $lowest_accepted_score, Cumulative Length = $cumulative_cj_length bases\n";
	
	## Save summary statistics
	$hcs->{total}->{number} = $total_candidate_junction_number;	
	$hcs->{total}->{length} = $total_cumulative_cj_length;
	
	$hcs->{accepted}->{number} = $accepted_candidate_junction_number;	
	$hcs->{accepted}->{length} = $cumulative_cj_length;
	$hcs->{accepted}->{score_cutoff} = $lowest_accepted_score;
	
	###
	## Print out the candidate junctions, sorted by the lower coordinate, higher coord, then number
	###
	sub by_ref_seq_coord
	{
		my $acj = Breseq::Shared::junction_name_split($a->{id});
		my $bcj = Breseq::Shared::junction_name_split($b->{id});
		return (($ref_seq_info->{seq_order}->{$acj->{side_1}->{seq_id}} <=> $ref_seq_info->{seq_order}->{$bcj->{side_1}->{seq_id}}) ||  ($acj->{side_1}->{position} <=> $bcj->{side_1}->{position})); 
	}
	
	@combined_candidate_junctions = sort by_ref_seq_coord @combined_candidate_junctions;
	
	foreach my $junction (@combined_candidate_junctions)
	{
		#print Dumper($ids_to_print);
		
		my $seq = Bio::Seq->new(
			-display_id => $junction->{id}, -seq => $junction->{seq});
		$out->write_seq($seq);
	}
	
	## create SAM faidx
	my $samtools = $settings->ctool('samtools');
	Breseq::Shared::system("$samtools faidx $candidate_junction_fasta_file_name") if (scalar @combined_candidate_junctions > 0);
	
	$summary->{candidate_junction} = $hcs;
}

sub _alignments_to_candidate_junctions
{
	my ($settings, $summary, $ref_seq_info, $candidate_junctions, $fai, $header, @al) = @_;

	my $verbose = 0;

	### TO DO:
	### We don't want to predict a new junction candidates from hits that are completely contained
	### in other perfect alignments. This can lead to problems later for repetitive regions.
	### Furthermore, it's just a waste of time... 
		
	### We are only concerned with combining intervals that are different in the reference they match
	#print "Before removing same-reference alignments: " . (scalar @al) . "\n" if ($verbose);
	my ($al_ref, $r_ref) = _unique_alignments_with_redundancy($fai, $header, @al);
	#print "After removing same-reference alignments: " . (scalar @al) . "\n" if ($verbose);

	return if (scalar @$al_ref == 1);
				
	### Now is our chance to decide which groups of matches are compatible,
	### to change their boundaries and to come up with a final list.		
	
	### For keeping track of how many times unique reference sequences (ignoring overlap regions)
	### were used to construct a junction. We must mark redundant sides AFTER correcting for overlap.
	my %redundant_junction_sides;
	my @junctions;	
		
	### Try adding together each pair of matches to make a longer match
	A1: for (my $i=0; $i<scalar @$al_ref; $i++)
	{		
		my $a1 = $al_ref->[$i];
		my ($a1_start, $a1_end) = Breseq::Shared::alignment_query_start_end($a1);	
						
		next A1 if ($a1_start == 0); #this is only true if read has no matches
		my $a1_length = $a1_end - $a1_start + 1;

		A2: for (my $j=$i+1; $j<scalar @$al_ref; $j++)
		{
			my $a2 = $al_ref->[$j];					
			my ($a2_start, $a2_end) = Breseq::Shared::alignment_query_start_end($a2);
			next A2 if ($a2_start == 0); #this is only true if read has no matches
			my $a2_length = $a2_end - $a2_start + 1;
			
			#print join(" ", $a1_start, $a1_end, $a2_start, $a2_end) . "\n" if ($verbose);
						
			my $union_start = ($a1_start < $a2_start) ? $a1_start : $a2_start;
			my $union_end = ($a1_end > $a2_end) ? $a1_end : $a2_end;			
			my $union_length = $union_end - $union_start + 1;

			my $intersection_start = ($a1_start > $a2_start) ? $a1_start : $a2_start;
			my $intersection_end = ($a1_end < $a2_end) ? $a1_end : $a2_end;  
			my $intersection_length = $intersection_end - $intersection_start + 1;
		
			#intersection can be a negative number
			$union_length += $intersection_length if ($intersection_length<0);
			my $intersection_length_nonzero = ($intersection_length > 0) ? $intersection_length : 0;			
						
			#print "$i $j $union_start-$union_end $union_length $intersection_start-$intersection_end $intersection_length\n" if ($verbose);		

			## If the intersection is just one of the hits, then there is no point...
			## At least one side should have this much unique sequence.
			my $extra_length = 10;
			
			if ( (($a1_length >= $intersection_length_nonzero+$extra_length) && ($a2_length > $intersection_length_nonzero))
			  || (($a1_length > $intersection_length_nonzero) && ($a2_length >= $intersection_length_nonzero+$extra_length)) )
			# We used to make sure BOTH sides had at least this much extra length			
			#if (($a1_length > $intersection_length_nonzero + $extra_length) && ($a2_length > $intersection_length_nonzero + $extra_length))
			{
				my $r1 = $r_ref->[$i];
				my $r2 = $r_ref->[$j]; 
				
				## we pass back and forth the redundancies in case they switch sides
				my ($junction_seq_string, $side_1_ref_seq, $side_2_ref_seq, $junction_coord_1, $junction_coord_2, @junction_id_list);
				
				($junction_seq_string, $side_1_ref_seq, $side_2_ref_seq, $junction_coord_1, $junction_coord_2, $r1, $r2, @junction_id_list)
					= _alignments_to_candidate_junction($settings, $summary, $ref_seq_info, $fai, $header, $a1, $a2, $r1, $r2);
							
				# a value of zero gets added if they were unique, >0 if redundant b/c they matched same reference sequence
				$redundant_junction_sides{$side_1_ref_seq}->{$junction_coord_1} += $r1-1;
				$redundant_junction_sides{$side_2_ref_seq}->{$junction_coord_2} += $r2-1;
				
				# also add to the reverse complement, because we can't be sure of the strandedness
				# (alternately we could reverse complement if the first base was an A or C, for example
				## it seems like there could possibly be some cross-talk between sides of junctions here, that
				## could snarl things up, but I'm not sure?
				## TO DO: I'm too tired of this section to do it now, but the correct sequence strand could be 
				## decided (keeping track of all the reversals) in _alignments_to_candidate_junction
				$redundant_junction_sides{Breseq::Fastq::revcom($side_1_ref_seq)}->{$junction_coord_1} += $r1-1;
				$redundant_junction_sides{Breseq::Fastq::revcom($side_2_ref_seq)}->{$junction_coord_2} += $r2-1;
				
				## Add a score based on the current read. We want to favor junctions with reads that overlap each side quite a bit
				## so we add the minimum that the read extends into each side of the candidate junction (not counting the overlap).
				my $a1_length_diff = ($a1_length - $intersection_length_nonzero);
				my $a2_length_diff = ($a2_length - $intersection_length_nonzero);
				
				my $score = ($a1_length_diff < $a2_length_diff) ? $a1_length_diff : $a2_length_diff;
				#print "score: $score\n";
				
				### Old score squared both sides (giving advantage to very uneven overlap on each side)
				#$candidate_junctions->{$junction_seq_string}->{$junction_id} += ($a1_length - $union_length) ** 2 + ($a2_length - $union_length) ** 2);

				push @junctions, { 'list' => \@junction_id_list, 'string' => $junction_seq_string, 'score' => $score, 'side_1_ref_seq' => $side_1_ref_seq, 'side_2_ref_seq' => $side_2_ref_seq};
			}
		}
	}	
	
	return if (scalar @junctions == 0);
	print $al_ref->[0]->qname . "\n"  if ($verbose);
	
	
	#only now that we've looked through everything can we determine whether the reference sequence matched 
	#on a side was unique, after correcting for overlap
	print Dumper(\%redundant_junction_sides) if ($verbose);
	
	foreach my $jct (@junctions)
	{	
		print Dumper($jct) if ($verbose);
		
		my @junction_id_list = @{$jct->{list}};
		my $junction_seq_string = $jct->{string};
		my $score = $jct->{score};
		
		my $side_1_ref_seq = $jct->{side_1_ref_seq};
		my $side_2_ref_seq = $jct->{side_2_ref_seq};

		## these are the total number of redundant matches to that side of the junction
		# the only way to be unique is to have at most one coordinate corresponding to that sequence
		# and not have that reference sequence redundantly matched (first combination of alignment coords)
		my $total_r1 = 0;
		foreach my $key (keys %{$redundant_junction_sides{$side_1_ref_seq}})
		{
			$total_r1 += $redundant_junction_sides{$side_1_ref_seq}->{$key} + 1;
		}
		
		my $total_r2 = 0;
		foreach my $key (keys %{$redundant_junction_sides{$side_2_ref_seq}})
		{
			$total_r2 += $redundant_junction_sides{$side_2_ref_seq}->{$key} + 1;
		}
		
		my $junction_id = Breseq::Shared::junction_name_join(@junction_id_list);
		print "$junction_id\n" if ($verbose);
		
		## initialize candidate junction if it didn't exist
		## they are redundant by default, until proven otherwise
		my $cj = $candidate_junctions->{$junction_seq_string}->{$junction_id};
		
		if (!defined $cj)
		{
			##redundancies of each side
			$cj->{r1} = 1;
			$cj->{r2} = 1;
			
			##maximum nonoverlapping match size on each side
			$cj->{L1} = 0;
			$cj->{L2} = 0;
			
			$candidate_junctions->{$junction_seq_string}->{$junction_id} = $cj;
		}
		
		## Update score of junction and the redundancy of each side
		$cj->{score} += $score;
		
		print "Totals: $total_r1, $total_r2\n" if ($verbose);
		print "Redundancy (before): $cj->{r1} ($cj->{L1}) $cj->{r2} ($cj->{L2})\n" if ($verbose);		
		my $side_1_ref_match_length = length $side_1_ref_seq;
		my $side_2_ref_match_length = length $side_2_ref_seq;
		
		
		## Longer reads into a side trump the redundancy of shorter reads for two reasons
		##   (1) obviously, the longer the read, the better the chance it is unique
		##   (2) subtly, if the read barely has enough to match both sides of the junction,
		##       there are situations where you can predict uniqueness because a short match
		##       maps one place only with the overlap included, and the non-overlapping part is
		##       unique, but once you have a longer match, you see that the non-overlapping
		##       part was really redundant.
		if ($side_1_ref_match_length > $cj->{L1})
		{
			$cj->{L1} = $side_1_ref_match_length;
			$cj->{r1} = ($total_r1 > 1) ? 1 : 0;
		}
		
		if ($side_2_ref_match_length > $cj->{L2})
		{
			$cj->{L2} = $side_1_ref_match_length;
			$cj->{r2} = ($total_r2 > 1) ? 1 : 0;;
		}		
		print "Redundancy (after): $cj->{r1} ($cj->{L1}) $cj->{r2} ($cj->{L2})\n" if ($verbose);		


	}	
	
}


### May dramatically help 454 data analysis, longer reads will lead to 
### more spurious candidate junctions.
### experimental: this should pay attention to how many mismatches there are!@
sub _uncontained_alignments
{
	my (@al) = @_;
	
	## assumes longer alignments are given first
	A1: for (my $i=0; $i<scalar @al; $i++)
	{		
		my $a1 = $al[$i];
		my ($a1_start, $a1_end) = Breseq::Shared::alignment_query_start_end($a1);			
		next A1 if ($a1_start == 0); #this is only true if read has no matches
		my $a1_length = $a1_end - $a1_start + 1;

		A2: for (my $j=$i+1; $j<scalar @al; $j++)
		{
			my $a2 = $al[$i];
			my ($a2_start, $a2_end) = Breseq::Shared::alignment_query_start_end($a2);			
			my $a2_length = $a2_end - $a2_start + 1;
			
			if ( ($a2_start >= $a1_start) && ($a2_start <= $a1_end ) && ($a1_length > $a2_length) )
			{
				splice @al, $j, 1;
				$j--;
			}
		}
	}

	return @al;
}


## Combines alignments where reads match to exactly
## the same reference sequence, so that there are
## fewer pairs to consider when making candidate junction
## sequences.
sub _unique_alignments_with_redundancy
{
	my ($fai, $header, @al) = @_;
	my @new_al;
	my $matched_reference_sequence_hash;
	my @redundancy;
	foreach my $a (@al)
	{
		my $strand = 1 - 2 * $a->reversed;
		my $interval = $header->target_name()->[$a->tid] . ":" . $a->start . "-" . $a->end;
		my $dna = $fai->fetch($interval);
		$dna = Breseq::Fastq::revcom($dna) if ($strand == -1);
				
		my $found = $matched_reference_sequence_hash->{$dna};
		if ($found)
		{
			$redundancy[$found]++;
			next;
		}
		
		push @new_al, $a;
		$matched_reference_sequence_hash->{$dna} = scalar(@redundancy);
		push @redundancy, 1;
	}
	
	return (\@new_al, \@redundancy);
}

=head2 alignments_to_candidate_junction

 Title   : alignments_to_candidate_junction
 Usage   : alignments_to_candidate_junction( );
 Function: Creates sequences that reproduce the junctions predicted in the list of hybrid reads.
 Returns : 
 
=cut

sub _alignments_to_candidate_junction
{
	my ($settings, $summary, $ref_seq_info, $fai, $header, $a1, $a2, $redundancy_1, $redundancy_2) = @_;
		
	my $verbose = 0;
	#$verbose = $a1->qname eq "30K88AAXX_LenskiSet2:8:3:1374:1537";
	#$verbose = ($a1->start == 1) || ($a2->start == 1);	
		
	## set up local settings
	my $flanking_length = $settings->{max_read_length};
	my $reference_sequence_string_hash_ref = $ref_seq_info->{ref_strings};

	### Method
	###
	### Hash junctions by a key showing the inner coordinate of the read.
	### and the direction propagating from that position in the reference sequence.
	### Prefer the unique or lower coordinate side of the read for main hash.
	###
	### REL606__1__1__REL606__4629812__0__0__
	### means the junction sequence is 36-1 + 4629812-4629777 from the reference sequence
	###
	### On the LEFT side:  0 means this is highest coord of alignment, junction seq begins at lower coord
	###                    1 means this is lowest coord of alignment, junction seq begins at higher coord
	### On the RIGHT side: 0 means this is highest coord of alignment, junction seq continues to lower coord
	###                    1 means this is lowest coord of alignment, junction seq continues to higher coord
	###
	### Note that there are more fields now than shown in this example...
	###
	### Need the junction key to include the offset to get to the junction within the read for cases where
	### the junction is near the end of the sequence... test case exists in JEB574.
				
	my %printed_keys;
	my $i = 0;

	my $read_id = $a1->query->name;
	
	#First, sort matches by their order in the query
	my ($q1, $q2) = ($a1, $a2);
	my ($q1_start, $q1_end) = Breseq::Shared::alignment_query_start_end($q1);
	my ($q2_start, $q2_end) = Breseq::Shared::alignment_query_start_end($q2);
		
	print "$q1_start, $q1_end, $q2_start, $q2_end\n" if ($verbose);
	# Reverse the coordinates to be consistently such that 1 refers to lowest...
	if ($q2_start < $q1_start)
	{
		($q1, $q2, $q1_start, $q2_start, $q1_end, $q2_end) = ($q2, $q1, $q2_start, $q1_start, $q2_end, $q1_end);
		($redundancy_1, $redundancy_2) = ($redundancy_2, $redundancy_1);	
	}
	#create hash key and store information about the location of this hit
	my $hash_strand_1 = ($q1->reversed) ? +1 : 0;
	my $hash_seq_id_1 = $header->target_name()->[$q1->tid];

	my $hash_strand_2 = ($q2->reversed) ? 0 : +1;
	my $hash_seq_id_2 = $header->target_name()->[$q2->tid];

	#how much overlap is there between the two matches?
	#positive if the middle sequence can match either side of the read
	#negative if there is sequence in the read NOT matched on either side 
	my $overlap = -($q2_start - $q1_end - 1);

	###
	## OVERLAP MISMATCH CORRECTION
	###
	# If there are mismatches in one or the other read in the overlap region
	# then we need to adjust the coordinates. Why? All sequences that we
	# retrieve are from the reference sequence and there are two choices
	# for where to extract this non-necessarily identical sequence!
	
	## save these as variables, because we may have to adjust them
	my ($r1_start, $r1_end) = ($q1->start, $q1->end);
	my ($r2_start, $r2_end) = ($q2->start, $q2->end);

	if ($verbose)
	{
		my $ref_seq_matched_1 = substr $reference_sequence_string_hash_ref->{$hash_seq_id_1},  $r1_start-1, $r1_end - $r1_start + 1;
		my $ref_seq_matched_2 = substr $reference_sequence_string_hash_ref->{$hash_seq_id_2},  $r2_start-1, $r2_end - $r2_start + 1;

		print "Alignment #1\n";
		print "qpos: $q1_start-$q1_end rpos: $r1_start-$r1_end reversed: " . $q1->reversed . "\n"; 
		print $q1->qseq . "\n" . $ref_seq_matched_1 . "\n";
		print Dumper($q1->cigar_array);

		print "Alignment #2\n";
		print "qpos: $q2_start-$q2_end rpos: $r2_start-$r2_end reversed: " . $q2->reversed . "\n"; 
		print $q2->qseq . "\n" . $ref_seq_matched_2 . "\n";
		print Dumper($q2->cigar_array);
	}

	##debug print information
	print "=== overlap: " . $overlap . "\n" if ($verbose);

	## Adjust the overlap in cases where there is a mismatch within the overlap region
	if ($overlap > 0)
	{
		my ($q1_move, $r1_move) = _num_matches_from_end($q1, $reference_sequence_string_hash_ref->{$hash_seq_id_1}, -1, $overlap);
		my ($q2_move, $r2_move) = _num_matches_from_end($q2, $reference_sequence_string_hash_ref->{$hash_seq_id_2}, +1, $overlap);

		if ($q1_move)
		{
			print "ALIGNMENT #1 OVERLAP MISMATCH: $q1_move, $r1_move\n" if ($verbose);
			#change where it ENDS
			$q1_end -= $q1_move;
			if ($q1->reversed)
			{
				$r1_start += $r1_move;
			}
			else
			{
				$r1_end -= $r1_move;
			}
		}
		
		if ($q2_move)
		{
			print "ALIGNMENT #2 OVERLAP MISMATCH: $q2_move, $r2_move\n" if ($verbose);
			#change where it STARTS
			$q2_start += $q2_move;
			if (!$q2->reversed)
			{
				$r2_start += $r2_move;
			}
			else
			{
				$r2_end -= $r2_move;
			}	
		}
		
		#re-calculate the overlap
		$overlap = -($q2_start - $q1_end - 1);
		print "=== overlap corrected for mismatches $overlap\n" if ($verbose);
	}
	
	## create hash coords AFTER overlap adjustment
	my $hash_coord_1 =  ($hash_strand_1 == +1) ? $r1_start : $r1_end;
	my $hash_coord_2 = ($hash_strand_2 == +1) ? $r2_start : $r2_end;
	
	if ($verbose)
	{
		my $ref_seq_m_1 = substr $reference_sequence_string_hash_ref->{$hash_seq_id_1},  $r1_start-1, $r1_end - $r1_start + 1;
		my $ref_seq_m_2 = substr $reference_sequence_string_hash_ref->{$hash_seq_id_2},  $r2_start-1, $r2_end - $r2_start + 1;
		
		print "Alignment #1\n";
		print "qpos: $q1_start-$q1_end rpos: $r1_start-$r1_end reversed: " . $q1->reversed . "\n"; 
		print $q1->qseq . "\n" . $ref_seq_m_1 . "\n";
		print Dumper($q1->cigar_array);

		print "Alignment #2\n";
		print "qpos: $q2_start-$q2_end rpos: $r2_start-$r2_end reversed: " . $q2->reversed . "\n"; 
		print $q2->qseq . "\n" . $ref_seq_m_2 . "\n";
		print Dumper($q2->cigar_array);
	}
	
	
	## Calculate an offset that only applies if the overlap is positive (sequence is shared between the two ends)
	my $overlap_offset = 0;
	$overlap_offset = $overlap if ($overlap > 0);
	print "Overlap offset: $overlap_offset\n" if ($verbose);
	
	## record what parts of the reference sequence were actually matched on each side
	## this is to determine whether that side was redundant or unique in the reference

	my ($ref_seq_matched_1, $ref_seq_matched_2);
	my $ref_seq_matched_length_1 = $r1_end - $r1_start + 1;
	my $ref_seq_matched_length_2 = $r2_end - $r2_start + 1;

	### create the sequence of the candidate junction 	
	my $junction_seq_string = '';

	##first end
	my $flanking_left = $flanking_length;
	if ($hash_strand_1 == 0) #alignment is not reversed
	{
		## $start_pos is in 1-based coordinates
		my $start_pos = $hash_coord_1 - ($flanking_left - 1) - $overlap_offset;
		if ($start_pos < 1)
		{
			print "START POS 1: $start_pos < 0\n" if ($verbose);
			$flanking_left += ($start_pos-1);
			$start_pos = 1;
		}
		my $add_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_1},  $start_pos-1, $flanking_left+$overlap_offset;
		print "1F: $add_seq\n" if ($verbose);
		$junction_seq_string .= $add_seq;
		
		$ref_seq_matched_1 = substr $add_seq, -$overlap_offset-$ref_seq_matched_length_1, $ref_seq_matched_length_1;
	}
	else 
	{	
		## $end_pos is in 1-based coordinates
		my $end_pos = $hash_coord_1 + ($flanking_left - 1) + $overlap_offset;
		if ($end_pos > length $reference_sequence_string_hash_ref->{$hash_seq_id_1})
		{
			print "END POS 1: ($end_pos < length\n" if ($verbose);
			$flanking_left -= ($end_pos - length $reference_sequence_string_hash_ref->{$hash_seq_id_1});
			$end_pos = length $reference_sequence_string_hash_ref->{$hash_seq_id_1};
		}	
		
		my $add_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_1},  $end_pos - ($flanking_left+$overlap_offset), $flanking_left+$overlap_offset;
		$add_seq = Breseq::Fastq::revcom($add_seq);
		print "1R: $add_seq\n" if ($verbose);
		$junction_seq_string .= $add_seq;
		
		$ref_seq_matched_1 = substr $add_seq, -$overlap_offset-$ref_seq_matched_length_1, $ref_seq_matched_length_1;
	}
	
	## Add any unique junction sequence that was only in the read
	## and NOT present in the reference genome		
	my $unique_read_seq_string = '';
	$unique_read_seq_string = substr $q1->qseq, $q1_end, -$overlap if ($overlap < 0);
	$junction_seq_string .= $unique_read_seq_string;
	
	print "+U: $unique_read_seq_string\n" if ($verbose);
		
	##second end
	my $flanking_right = $flanking_length;
	if ($hash_strand_2 == +1) #alignment is not reversed 
	{
		## $end_pos is in 1-based coordinates
		my $end_pos = $hash_coord_2 + ($flanking_right - 1) + $overlap_offset;
		if ($end_pos > length $reference_sequence_string_hash_ref->{$hash_seq_id_2})
		{
			print "END POS 2: ($end_pos < length\n" if ($verbose);
			$flanking_right -= ($end_pos - length $reference_sequence_string_hash_ref->{$hash_seq_id_2});
			$end_pos = length $reference_sequence_string_hash_ref->{$hash_seq_id_2};
		}	
		my $add_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_2},  $end_pos - $flanking_right, $flanking_right;
		print "2F: $add_seq\n" if ($verbose);
		$junction_seq_string .= $add_seq;
		
		$ref_seq_matched_2 = substr $add_seq, -$ref_seq_matched_length_2, $ref_seq_matched_length_2;
	}
	else # ($m_2->{hash_strand} * $m_2->{read_side} == -1)
	{
		## $start_pos is in 1-based coordinates
		my $start_pos = $hash_coord_2 - ($flanking_right - 1) - $overlap_offset;
		if ($start_pos < 1)
		{
			print "START POS 2: $start_pos < 0\n" if ($verbose);
			$flanking_right += ($start_pos-1);
			$start_pos = 1;
		}
		my $add_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_2},  $start_pos-1, $flanking_right;
		$add_seq = Breseq::Fastq::revcom($add_seq);
		print "2R: $add_seq\n" if ($verbose);
		$junction_seq_string .= $add_seq;
		
		$ref_seq_matched_2 = substr $add_seq, -$ref_seq_matched_length_2, $ref_seq_matched_length_2;
	}	
	print "3: $junction_seq_string\n" if ($verbose);
	
	
	#print "$read_id\n";
	#print "$junction_id\n";
	#print "$junction_seq_string\n";
	
#	my $hash_redundancy_1 = ($r1 > 1) ? 1 : 0;
#	my $hash_redundancy_2 = ($r2 > 1) ? 1 : 0;

	## create hash coords after this adjustment
	if ($hash_strand_1 == 0)
	{
		$r1_end -= $overlap_offset;
	}
	else #reversed
	{
		$r1_start += $overlap_offset;
	}

	if ($hash_strand_2 == 0)
	{
		$r2_end -= $overlap_offset;
	}
	else #reversed
	{
		$r2_start += $overlap_offset;
	}

	##matched reference sequence
	$ref_seq_matched_1 = substr $reference_sequence_string_hash_ref->{$hash_seq_id_1},  $r1_start-1, $r1_end - $r1_start + 1;
	$ref_seq_matched_2 = substr $reference_sequence_string_hash_ref->{$hash_seq_id_2},  $r2_start-1, $r2_end - $r2_start + 1;

	
	#want to be sure that lowest ref coord is always first for consistency 
	if ( ($hash_seq_id_2 lt $hash_seq_id_1) || (($hash_seq_id_1 eq $hash_seq_id_2) && ($hash_coord_2 < $hash_coord_1)) )
	{
		($hash_coord_1, $hash_coord_2) = ($hash_coord_2, $hash_coord_1);
		($hash_strand_1, $hash_strand_2) = ($hash_strand_2, $hash_strand_1);
		($hash_seq_id_1, $hash_seq_id_2) = ($hash_seq_id_2, $hash_seq_id_1);
		($redundancy_1, $redundancy_2) = ($redundancy_2, $redundancy_1);
		($flanking_left, $flanking_right) = ($flanking_right, $flanking_left);
		($ref_seq_matched_1, $ref_seq_matched_2) = ($ref_seq_matched_2, $ref_seq_matched_1);
			
		$junction_seq_string = Breseq::Fastq::revcom($junction_seq_string);
		$unique_read_seq_string = Breseq::Fastq::revcom($unique_read_seq_string);
	}
			
	my @junction_id_list = (
		$hash_seq_id_1, 
		$hash_coord_1, 
		$hash_strand_1, 
		$hash_seq_id_2, 
		$hash_coord_2, 
		$hash_strand_2, 
		$overlap, 
		$unique_read_seq_string, 
		$flanking_left, 
		$flanking_right
	);	
	
	my $junction_coord_1 = join "::", $hash_seq_id_1, $hash_coord_1, $hash_strand_1;
	my $junction_coord_2 = join "::", $hash_seq_id_2, $hash_coord_2, $hash_strand_2;
		
	my $junction_id = Breseq::Shared::junction_name_join(@junction_id_list);
	print "JUNCTION ID: " . $junction_id . "\n" if ($verbose);
	
	die "Junction sequence not found: $junction_id " . $q1->qname . " " . $a2->qname  if (!$junction_seq_string);
	die "Incorrect length for $junction_seq_string: $junction_id " . $q1->qname . " " . $a2->qname if (length $junction_seq_string != $flanking_left + $flanking_right + abs($overlap));
		
	return ($junction_seq_string, $ref_seq_matched_1, $ref_seq_matched_2, $junction_coord_1, $junction_coord_2, $redundancy_1, $redundancy_2, @junction_id_list);
}


## find the maximum number of positions one can go from the end without finding a mismatch
## returns undef if this is a perfect match. (it doesn't matter whether returned numbers 
## are in query or in reference, since they will be the same.) If there are mismatches,
## we want the one that is closest to the desired overlap without going over...
##
## $dir can be:
##    -1 , which implies start at the highest read coord and go to lowest
##    +1 , which implies start at the lowest read coord and go to highest

sub _num_matches_from_end
{
	my $verbose = 0;
	my ($a, $refseq_str, $dir, $overlap) = @_;
	
	my $reversed = $a->reversed;
	my ($q_start, $q_end) = Breseq::Shared::alignment_query_start_end($a);
	my $q_length = $q_end - $q_start + 1;
	my ($q_seq_start, $q_seq_end) = ($q_start, $q_end);
	($q_seq_start, $q_seq_end) = ($a->l_qseq - $q_seq_end + 1, $a->l_qseq - $q_seq_start + 1) if ($reversed);
	my $q_str = substr $a->qseq, $q_seq_start-1, $q_seq_end-$q_seq_start+1;
	my @q_array = split //, $q_str;
	
	my ($r_start, $r_end, $r_length) = ($a->start, $a->end, $a->length);
	my $r_str = substr $refseq_str, $r_start-1, $r_length;
	my @r_array = split //, $r_str;
	
	if ($verbose)
	{
		print "==== " . $a->qname . "\n";
		print "direction: $dir\n";
		print $a->qseq . "\n";
		print "$q_start-$q_end $q_seq_start-$q_seq_end $reversed\n";
		print "$q_str\n";
		print "$r_start-$r_end\n";
		print "$r_str\n";
	}
	
	$dir = -1 * $dir if ($reversed);
	@q_array = reverse @q_array if ($dir == -1);
	@r_array = reverse @r_array if ($dir == -1);
	
	my $cigar_array_ref = $a->cigar_array;
	##remove soft padding
	shift @$cigar_array_ref if ($cigar_array_ref->[0]->[0] eq 'S');
	pop @$cigar_array_ref if ($cigar_array_ref->[-1]->[0] eq 'S');
	@$cigar_array_ref = reverse @$cigar_array_ref if ($dir == -1);

	my $qry_mismatch_pos = undef;
	my $ref_mismatch_pos = undef;
	
	my $r_pos = 0;
	my $q_pos = 0;
	POS: while (($q_pos < $overlap) && ($r_pos < $r_length) && ($q_pos < $q_length))
	{
		#get rid of previous items...
		shift @$cigar_array_ref if ($cigar_array_ref->[0]->[1]==0);
		$cigar_array_ref->[0]->[1]--;

		#handle indels
		if ($cigar_array_ref->[0]->[0] eq 'I')
		{
			$r_pos++;
			$ref_mismatch_pos = $r_pos;
			$qry_mismatch_pos = $q_pos;
		}
		elsif ($cigar_array_ref->[0]->[0] eq 'D')
		{
			$q_pos++;
			$ref_mismatch_pos = $r_pos;
			$qry_mismatch_pos = $q_pos;
		}
		else
		{
			if ($q_array[$q_pos] ne $r_array[$r_pos])
			{
				$ref_mismatch_pos = $r_pos;
				$qry_mismatch_pos = $q_pos;
			}
			$r_pos++;
			$q_pos++;
		}
		
		print "$r_pos $q_pos\n" if ($verbose);		
	}
	
	print "INTERESTING!!!\n" if ($verbose && $qry_mismatch_pos);
	
	## make 1 indexed...
	$qry_mismatch_pos++ if (defined $qry_mismatch_pos);
	$ref_mismatch_pos++ if (defined $ref_mismatch_pos);

	return ($qry_mismatch_pos, $ref_mismatch_pos);
}



return 1;

