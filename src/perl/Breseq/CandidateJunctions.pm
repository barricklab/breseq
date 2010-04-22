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
		
		$hcs->{$read_file} = $s;
	}		
		
		
	###	
	## Combine hash into a list, one item for each unique sequence (also joining reverse complements)
	###		
	
	my @combined_candidate_junctions;
	my $handled_seq;
	my %ids_to_print;
		
	JUNCTION_SEQ: foreach my $junction_seq (keys %$candidate_junctions)
	{	
		next if ($handled_seq->{$junction_seq});
		
		my @combined_candidate_junction_list = ();
		my $rc_junction_seq = Breseq::Fastq::revcom($junction_seq);
		
		JUNCTION_ID: foreach my $junction_id (keys %{$candidate_junctions->{$junction_seq}} )
		{
			push @combined_candidate_junction_list, { id=>$junction_id, score=>$candidate_junctions->{$junction_seq}->{$junction_id}, seq=>$junction_seq, rc_seq=>$rc_junction_seq };
		}
		$handled_seq->{$junction_seq}++;
		
		## add the reverse complement
		if (defined $candidate_junctions->{$rc_junction_seq})
		{
			JUNCTION_ID: foreach my $junction_id (keys %{$candidate_junctions->{$rc_junction_seq}} )
			{
				push @combined_candidate_junction_list, { id=>$junction_id, score=>$candidate_junctions->{$rc_junction_seq}->{$junction_id}, seq=>$rc_junction_seq, rc_seq=>$junction_seq };
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
	## Limit the number of candidate junctions that we print (and record the cutoff)
	##   eliminate subsequences as we go.
	###
	
#	print STDERR "  Before taking top candidates:" . (scalar @combined_candidate_junctions) . "\n";
	
	my @duplicate_sequences;
	if ((defined $settings->{maximum_candidate_junctions}) && (@combined_candidate_junctions > 0))
	{	
		my $cutoff_reads_for_candidate_junction;
		my @remaining_ids = ();
		my @list_in_waiting = ();
	
		my $i = 0;
		my $current_reads_for_candidate_junction = $combined_candidate_junctions[$i]->{score};
	
		while (scalar(@remaining_ids) + scalar(@list_in_waiting) <= $settings->{maximum_candidate_junctions})
		{			
			push @remaining_ids, @list_in_waiting;
			@list_in_waiting = ();
			$cutoff_reads_for_candidate_junction = $current_reads_for_candidate_junction;
			last if ($i >= scalar @combined_candidate_junctions);

			$current_reads_for_candidate_junction = $combined_candidate_junctions[$i]->{score};
			CAND: while (($i < scalar @combined_candidate_junctions) && ($current_reads_for_candidate_junction == $combined_candidate_junctions[$i]->{score}))
			{
				##Check here 
				my $c = $combined_candidate_junctions[$i];
				foreach my $dup (@duplicate_sequences)
				{
					#is this a subsequence of one already in the list?
					if ( ($dup =~ m/$c->{seq}/) || ($dup =~ m/$c->{rc_seq}/) || ($c->{seq} =~ m/$dup/) || ($c->{rc_seq} =~ m/$dup/) ) 
					{
						print "Dup $dup\n$c->{seq}\n" if ($verbose);
						next CAND;
					}
				}
				
				push @duplicate_sequences, $c->{seq};
				push @list_in_waiting, $c;
			} continue {
				$i++;
			}
			
#			print "    Number with junction score $current_reads_for_candidate_junction = " . (scalar @list_in_waiting) . "\n";
		}
		@combined_candidate_junctions = @remaining_ids;
		$summary->{cutoff_reads_for_candidate_junction} = $cutoff_reads_for_candidate_junction;
		
#		print STDERR "  Junction score cutoff for printing: $cutoff_reads_for_candidate_junction\n";
	}
	
#	print STDERR "  After taking top candidates:" . (scalar @combined_candidate_junctions) . "\n";
	
	
	###
	## Print out the candidate junctions, sorted by the lower coordinate, higher coord, then number
	###
	sub by_ref_seq_coord
	{
		my @al = Breseq::Shared::junction_name_split($a->{id});
		my @bl = Breseq::Shared::junction_name_split($b->{id});
		return (($ref_seq_info->{seq_order}->{$al[0]} <=> $ref_seq_info->{seq_order}->{$bl[0]}) ||  ($al[1] <=> $bl[1])); 
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
	Breseq::Shared::system("samtools faidx $candidate_junction_fasta_file_name") if (scalar @combined_candidate_junctions > 0);
	
	$summary->{candidate_junction} = $hcs;
}

sub _alignments_to_candidate_junctions
{
	my ($settings, $summary, $ref_seq_info, $candidate_junctions, $fai, $header, @al) = @_;

	my $verbose = 0;
		
	### We are only concerned with different matches.
	#print "Before removing same-reference alignments: " . (scalar @al) . "\n" if ($verbose);
	
	my ($al_ref, $r_ref) = _unique_alignments_with_redundancy($fai, $header, @al);
	#print "After removing same-reference alignments: " . (scalar @al) . "\n" if ($verbose);

	return if (scalar @$al_ref == 1);
				
	### Now is our chance to decide which groups of matches are compatible,
	### to change their boundaries and to come up with a final list.		
	
	### for keeping track of redundancy
	my %redundant_junction_sides;	
	
	### Try adding together each pair of matches to make a longer match
	for (my $i=0; $i<scalar @$al_ref; $i++)
	{		
		for (my $j=$i+1; $j<scalar @$al_ref; $j++)
		{
			my $a1 = $al_ref->[$i];
			my $a2 = $al_ref->[$j];					
					
			my ($a1_start, $a1_end) = Breseq::Shared::alignment_query_start_end($a1);			
			my ($a2_start, $a2_end) = Breseq::Shared::alignment_query_start_end($a2);
		
			#the read may be there but have no match!
			next if (($a1_start == 0) || ($a2_start == 0));
			
			my $a1_length = $a1_end - $a1_start + 1;
			my $a2_length = $a2_end - $a2_start + 1;
			
			print join(" ", $a1_start, $a1_end, $a2_start, $a2_end) . "\n" if ($verbose);
						
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
		
			#limit on how different they should be
			my $extra_length = 10;
			my $extra_length_2 = 15;
			
			## if the intersection is one of the hits, then well... there is no point...
			if ( (($a1_length >= $intersection_length_nonzero+$extra_length) && ($a2_length > $intersection_length_nonzero))
			  || (($a1_length > $intersection_length_nonzero) && ($a2_length >= $intersection_length_nonzero+$extra_length)) )
#			if (($a1_length > $intersection_length_nonzero + $extra_length) && ($a2_length > $intersection_length_nonzero + $extra_length))
			{
				my $r1 = $r_ref->[$i];
				my $r2 = $r_ref->[$j]; 
				my ($junction_id, $junction_seq_string, $q1, $q2) = _alignments_to_candidate_junction($settings, $summary, $ref_seq_info, $fai, $header, $a1, $a2, $r1, $r2);

				#add to number of observations			
				$candidate_junctions->{$junction_seq_string}->{$junction_id} += ($a1_length - $union_length) ** 2 + ($a2_length - $union_length) ** 2;
			}
		}
	}		
}

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
	my ($settings, $summary, $ref_seq_info, $fai, $header, $a1, $a2, $r1, $r2) = @_;
		
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
	($q1, $q2, $q1_start, $q2_start, $q1_end, $q2_end) = ($q2, $q1, $q2_start, $q1_start, $q2_end, $q1_end) if ($q2_start < $q1_start);		
	
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
	#If there are mismatches in one or the other read in the overlap region
	#then we need to adjust the coordinates. Why? All sequences that we
	#retrieve are from the reference sequence and there are two choices
	#for where to extract this non-necessarily identical sequence
	
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
		print "=== new overlap $overlap\n" if ($verbose);
	
	}
	
	## create hash coords after this adjustment
	my $hash_coord_1 =  ($hash_strand_1 == +1) ? $r1_start : $r1_end;
	my $hash_coord_2 = ($hash_strand_2 == +1) ? $r2_start : $r2_end;

	##matched reference sequence
	my $ref_seq_matched_1 = substr $reference_sequence_string_hash_ref->{$hash_seq_id_1},  $r1_start-1, $r1_end - $r1_start + 1;
	my $ref_seq_matched_2 = substr $reference_sequence_string_hash_ref->{$hash_seq_id_2},  $r2_start-1, $r2_end - $r2_start + 1;
	
	if ($verbose)
	{
		print "Alignment #1\n";
		print "qpos: $q1_start-$q1_end rpos: $r1_start-$r1_end reversed: " . $q1->reversed . "\n"; 
		print $q1->qseq . "\n" . $ref_seq_matched_1 . "\n";
		print Dumper($q1->cigar_array);

		print "Alignment #2\n";
		print "qpos: $q2_start-$q2_end rpos: $r2_start-$r2_end reversed: " . $q2->reversed . "\n"; 
		print $q2->qseq . "\n" . $ref_seq_matched_2 . "\n";
		print Dumper($q2->cigar_array);
	}
	
	
	## this is the sequence NOT present in the reference genome		
	my $unique_read_seq_string = '';
	$unique_read_seq_string = substr $q1->qseq, $q1_end, -$overlap if ($overlap < 0);
	
	### create the sequence of the candidate junction 	
	my $junction_seq_string = '';

	my $first_query_match_seq = substr $q1->qseq, $q1_start-1, $q1_end - $q1_start + 1;
	if ($q1->reversed)
	{
		$first_query_match_seq = Breseq::Fastq::revcom($first_query_match_seq);
	}
	
	## this only applies if the overlap is positive (sequence is shared between the two ends)
	my $overlap_offset = 0;
	$overlap_offset = $overlap if ($overlap > 0);

	print "Overlap offset: $overlap_offset\n" if ($verbose);

	##first end
	my $first_overlap_seq = '';
	if ($hash_strand_1 == 0) #alignment is not reversed
	{
		my $add_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_1},  $hash_coord_1-($flanking_length-1)-1-$overlap_offset, $flanking_length+$overlap_offset;
		print "1F: $add_seq\n" if ($verbose);
		$junction_seq_string .= $add_seq;
		$first_overlap_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_1},  $hash_coord_1-$overlap, $overlap if ($overlap > 0);		
	}
	else 
	{
		my $add_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_1},  $hash_coord_1-1, $flanking_length+$overlap_offset;
		$add_seq = Breseq::Fastq::revcom($add_seq);
		print "1R: $add_seq\n" if ($verbose);
		$junction_seq_string .= $add_seq;

		if ($verbose && ($overlap > 0))
		{
			$first_overlap_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_1},  $hash_coord_1-1, $overlap;
			$first_overlap_seq = Breseq::Fastq::revcom($first_overlap_seq);
		}		
	}
	
	print "read: $first_query_match_seq\n" if ($verbose);
	print "$hash_strand_1 $overlap $first_overlap_seq\n" if ($verbose);
	$junction_seq_string .= $unique_read_seq_string;
	print "+: $unique_read_seq_string\n" if ($verbose);
		
	##second end
	my $second_overlap_seq = '';
	if ($hash_strand_2 == +1) #alignment is not reversed 
	{
		print "Size: " . ($flanking_length-$overlap_offset) . "\n" if ($verbose);
		
		my $add_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_2},  $hash_coord_2+$overlap_offset-1, $flanking_length;
		print "2F: $add_seq\n" if ($verbose);

		$second_overlap_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_2},  $hash_coord_2-1, $overlap if ($overlap > 0);

		$junction_seq_string .= $add_seq;
	}
	else # ($m_2->{hash_strand} * $m_2->{read_side} == -1)
	{
		my $add_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_2},  $hash_coord_2-($flanking_length-1)-$overlap_offset-1, $flanking_length;
		$add_seq = Breseq::Fastq::revcom($add_seq);
		
		print "2R: $add_seq\n" if ($verbose);

		$junction_seq_string .= $add_seq;
		
		if ($verbose && ($overlap > 0))
		{
			$second_overlap_seq = substr $reference_sequence_string_hash_ref->{$hash_seq_id_2},  $hash_coord_2-$overlap, $overlap;
			$second_overlap_seq = Breseq::Fastq::revcom($second_overlap_seq);
		}
				
	}
	print "$hash_strand_2 $overlap $second_overlap_seq\n" if ($verbose);
	
	print "3: $junction_seq_string\n" if ($verbose);
	
	
	#print "$read_id\n";
	#print "$junction_id\n";
	#print "$junction_seq_string\n";
	
	my $hash_redundancy_1 = ($r1 > 1) ? 1 : 0;
	my $hash_redundancy_2 = ($r2 > 1) ? 1 : 0;
	
	#want to be sure that lowest ref coord is always first for consistency 
	if ( ($hash_seq_id_2 lt $hash_seq_id_1) || (($hash_seq_id_1 eq $hash_seq_id_2) && ($hash_coord_2 < $hash_coord_1)) )
	{
		($hash_coord_1, $hash_coord_2) = ($hash_coord_2, $hash_coord_1);
		($hash_strand_1, $hash_strand_2) = ($hash_strand_2, $hash_strand_1);
		($hash_seq_id_1, $hash_seq_id_2) = ($hash_seq_id_2, $hash_seq_id_1);
		($hash_redundancy_1, $hash_redundancy_2) = ($hash_redundancy_2, $hash_redundancy_1);

		$junction_seq_string = Breseq::Fastq::revcom($junction_seq_string);
		$unique_read_seq_string = Breseq::Fastq::revcom($unique_read_seq_string);
	}
		
	my $junction_id = Breseq::Shared::junction_name_join($hash_seq_id_1, $hash_coord_1, $hash_strand_1, $hash_redundancy_1, $hash_seq_id_2, $hash_coord_2, $hash_strand_2, $hash_redundancy_2, $overlap, $unique_read_seq_string, $flanking_length);	
	print "JUNCTION ID: $junction_id\n" if ($verbose);
	
	die "Junction sequence not found: $junction_id " . $q1->qname . " " . $a2->qname  if (!$junction_seq_string);
## this guard accidentally catches junctions near the ends of sequences...
#	die "Incorrect length for $junction_seq_string: $junction_id " . $q1->qname . " " . $a2->qname if (length $junction_seq_string != 2*$flanking_length + abs($overlap));
		
	return ($junction_id, $junction_seq_string, $q1, $q2);
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

