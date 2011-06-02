/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#include "breseq/common.h"

#include "breseq/candidate_junctions.h"

using namespace std;

namespace breseq {

	// Private

	CandidateJunction::CandidateJunction() {}

	void CandidateJunction::_alignments_to_candidate_junction(map_t settings, map_t summary, map_t ref_seq_info, map_t fai, map_t header, map_t a1, map_t a2, map_t redundancy_1, map_t redundancy_2) {}

	void CandidateJunction::_alignments_to_candidate_junctions(Settings settings, Summary summary,  const cReferenceSequences& ref_seq_info, map_t candidate_junctions, faidx_t* fai, bam_header_t* header, vector<bam1_t*> al_ref)
	{
		bool verbose = false;

		if (verbose)
		{
			cout << endl << "###########################" << endl;
			cout << bam1_qname(al_ref[0]);
			cout << endl << "###########################" << endl;
		}

		// Must still have multiple matches to support a new junction.
		if (al_ref.size() <= 1)
			return;

		// Now is our chance to decide which groups of matches are compatible,
		// to change their boundaries and to come up with a final list.

		// For keeping track of how many times unique reference sequences (ignoring overlap regions)
		// were used to construct a junction. We must mark redundant sides AFTER correcting for overlap.
//		my %redundant_junction_sides;
//		my @junctions;

		if (verbose)
		{
			cout << bam1_qname(al_ref[0]) << endl;
			cout << "Total matches: " << al_ref.size() << endl;
		}

		vector<bam1_t*> list1, list2;

	  	// Try only pairs where one match starts at the beginning of the read >>>
		// This saves a number of comparisons and gets rid of a lot of bad matches.
		int32_t max_union_length = 0;

		for (int32_t i = 0; i < al_ref.size(); i++)
		{
			bam1_t* a = al_ref[i];

			int32_t a_start, a_end;
			alignment_query_start_end(a, a_start, a_end);

			if (verbose) cout << "(" << a_start << ", " << a_end << ")" << endl;

			int32_t length = a_end - a_start + 1;
			if (max_union_length < length)
				max_union_length = length;

			if (a_start == 1)
				list1.push_back(a);
			else
				list2.push_back(a);
		}

		// Alternately try all pairs
		//	list1 = al_ref;
		//	list2 = al_ref;

		// Pairs must CLEAR the maximum length of any one alignment by a certain amount
		max_union_length += settings.required_extra_pair_total_length;

		// The first match in this category is the longest
		if (verbose)
		{
			cout << "  List1: " << list1.size() << endl;
			cout << "  List2: " << list2.size() << endl;
		}

		vector<PassedPair> passed_pair_list;

		// Try adding together each pair of matches to make a junction, by looking at read coordinates
		for (int32_t i = 0; i < list1.size(); i++)
		{
			bam1_t* a1 = list1[i];

			int32_t a1_start, a1_end;
			alignment_query_start_end(a1, a1_start, a1_end);

			for (int32_t j = 0; j < list2.size(); j++)
			{
				bam1_t* a2 = list1[i];

				int32_t a2_start, a2_end;
				alignment_query_start_end(a2, a2_start, a2_end);

				// if either end is the same, it is going to fail.
				// Note: we already checked for the start, when we split into the two lists!

				// don't allow it to succeed no matter what if both of these reads are encompassed by a longer read from list one
				//#next A2 if ($a2_end <= $encompass_end);

				// Check a slew of guards to prevent predicting too many junctions to handle.
				int32_t a1_unique_length, a2_unique_length, union_length;
				bool passed = _check_read_pair_requirements(settings, a1_start, a1_end, a2_start, a2_end, a1_unique_length, a2_unique_length, union_length);

				if (passed && (union_length >= max_union_length))
				{
					//destroy any contained matches -- possibly leaving some leeway
					if (union_length > max_union_length)
					{
						max_union_length = union_length;
						for (vector<PassedPair>::iterator it = passed_pair_list.end() - 1; it >= passed_pair_list.begin(); it--)
							if ((*it).union_length < max_union_length)
								passed_pair_list.erase(it);
					}

					PassedPair passed_pair = { a1 = a1, a2 = a2, union_length = union_length, a1_unique_length = a1_unique_length, a2_unique_length = a2_unique_length };
					passed_pair_list.push_back(passed_pair);
				}

			}
		}
	}

	bool CandidateJunction::_check_read_pair_requirements(Settings settings, int32_t a1_start, int32_t a1_end, int32_t a2_start, int32_t a2_end, int32_t& a1_unique_length, int32_t& a2_unique_length, int32_t& union_length)
	{
		bool verbose = false;

		if (verbose)
			cout << "=== Match1: " << a1_start << "-" << a1_end << "   Match2: " << a2_start << "-" << a2_end << endl;

		// 0. Require one match to start at the beginning of the read
		if (a1_start != 1 && a2_start != 1)
			return false;
		// Already checked when two lists were constructed: TEST AND REMOVE

		int32_t a1_length = a1_end - a1_start + 1;
		int32_t a2_length = a2_end - a2_start + 1;

		int32_t union_start = (a1_start < a2_start) ? a1_start : a2_start;
		int32_t union_end = (a1_end > a2_end) ? a1_end : a2_end;
		union_length = union_end - union_start + 1;

		int32_t intersection_start = (a1_start > a2_start) ? a1_start : a2_start;
		int32_t intersection_end = (a1_end < a2_end) ? a1_end : a2_end;
		int32_t intersection_length = intersection_end - intersection_start + 1;

		if (intersection_length < 0)
			union_length += intersection_length;
		int32_t intersection_length_positive = (intersection_length > 0) ? intersection_length : 0;
		int32_t intersection_length_negative = (intersection_length < 0) ? -1 * intersection_length : 0;

		//
		// CHECKS
		//

		if (verbose)
			cout << "    Union: " << union_length << "   Intersection: " << intersection_length << endl;

		//// 1. Require maximum negative overlap (inserted unique sequence length) to be less than some value
		if (intersection_length_negative > settings.maximum_inserted_junction_sequence_length)
			return false;

		//// 2. Require both ends to extend a certain minimum length outside of the overlap
		if (a1_length < intersection_length_positive + settings.required_both_unique_length_per_side)
			return false;

		if (a2_length < intersection_length_positive + settings.required_both_unique_length_per_side);
			return false;

		//// 3. Require one end to extend a higher minimum length outside of the overlap
		if ((a1_length < intersection_length_positive + settings.required_one_unique_length_per_side)
				  && (a2_length < intersection_length_positive + settings.required_one_unique_length_per_side))
			return false;

		//// 4. Require both matches together to cover a minimum part of the read.
		if (union_length < settings.required_match_length)
			return false;

		//// Add a score based on the current read. We want to favor junctions with reads that overlap each side quite a bit
		//// so we add the minimum that the read extends into each side of the candidate junction (not counting the overlap).
		a1_unique_length = (a1_length - intersection_length_positive);
		a2_unique_length = (a2_length - intersection_length_positive);

		if (verbose)
			cout << "    Unique1: " << a1_unique_length << "   Unique2: " << a2_unique_length << endl;

		return true;
	}

	void CandidateJunction::_entire_read_matches(map_t a) {}
	void CandidateJunction::_num_matches_from_end(map_t a, map_t refseq_str, map_t dir, map_t overlap) {}
	void CandidateJunction::_split_indel_alignments(Settings settings, Summary summary, bam_header_t* header, ofstream& PSAM, int32_t min_indel_split_len, vector<bam1_t*> al_ref) {}
	void CandidateJunction::_by_ref_seq_coord(map_t a, map_t b, map_t ref_seq_info) {}
	void CandidateJunction::_by_score_unique_coord(map_t a, map_t b) {}
	void CandidateJunction::_tam_write_split_alignment(map_t fh, map_t header, map_t min_indel_split_len, map_t a) {}

	// Public

	/*! Preprocesses alignments
	 */
	void CandidateJunction::preprocess_alignments(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info)
	{
		cout << "Preprocessing alignments." << endl;

		// get the cutoff for splitting alignments with indels
		boost::optional<int32_t> min_indel_split_len = settings.preprocess_junction_min_indel_split_length;

		// includes best matches as they are
		string preprocess_junction_best_sam_file_name = settings.preprocess_junction_best_sam_file_name;
		ofstream BSAM;
		BSAM.open(preprocess_junction_best_sam_file_name.c_str());
		assert(BSAM.is_open());

		string reference_faidx_file_name = settings.reference_fasta_file_name;
		faidx_t* reference_fai = fai_load(reference_faidx_file_name.c_str());
		for (int32_t index = 0; index < settings.read_structures.size(); index++)
		{
			Settings::ReadStructure read_struct = settings.read_structures[index];
			cerr << "  READ FILE::" << read_struct.base_name << endl;

			string reference_sam_file_name = settings.reference_sam_file_name;
			string reference_faidx_file_name = settings.reference_faidx_file_name;

			tamFile tam = sam_open(reference_sam_file_name.c_str()); // or die("Could not open reference same file: $reference_sam_file_name");
			bam_header_t* header = sam_header_read2(reference_faidx_file_name.c_str()); // or die("Error reading reference fasta index file: $reference_faidx_file_name");

			// includes all matches, and splits long indels
			string preprocess_junction_split_sam_file_name = settings.preprocess_junction_split_sam_file_name;
			ofstream PSAM;
			PSAM.open(preprocess_junction_split_sam_file_name.c_str());
			assert(PSAM.is_open());

			bam1_t* last_alignment;
			vector<bam1_t*> al_ref;
			int32_t i = 0;
			while (true)
			{
				// resolve_alignments.c
				al_ref = tam_next_read_alignments(tam, header, last_alignment, false);

				if (al_ref.size() == 0)
					break;

				if (++i % 10000 == 0)
					cerr << "    ALIGNED READ:" << i << endl;

				// for testing...
				if (settings.candidate_junction_read_limit != 0 && i > settings.candidate_junction_read_limit) break;

				// write split alignments
				if (min_indel_split_len)
					_split_indel_alignments(settings, summary, header, PSAM, min_indel_split_len.get(), al_ref);

				// write best alignments
				if (settings.candidate_junction_score_method.compare("POS_HASH") == 0)
				{
					int32_t best_score = _eligible_read_alignments(settings, header, reference_fai, ref_seq_info, al_ref);
					tam_write_read_alignments(BSAM, header, 0, al_ref, boost::optional<vector<Trim> >());
				}
			}
			sam_close(tam);
			PSAM.close();
		}
		BSAM.close();
	}

	/*! Predicts candidate junctions
	 */
	void CandidateJunction::identify_candidate_junctions(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info)
	{
		int32_t verbose = 0;

		// set up some options that are global to this module
		/*
		$required_both_unique_length_per_side = $settings->{required_both_unique_length_per_side};
		$required_one_unique_length_per_side = $settings->{required_one_unique_length_per_side};
		$maximum_inserted_junction_sequence_length = $settings->{maximum_inserted_junction_sequence_length};
	 	$required_match_length = $settings->{required_match_length};
		$required_extra_pair_total_length = $settings->{required_extra_pair_total_length};
		*/

		// set up files and local variables from settings
		//my $ref_strings = ref_seq_info[*].m_fasta_sequence.m_sequence;

		string reference_faidx_file_name = settings.reference_faidx_file_name;
		faidx_t* fai = fai_load(settings.reference_fasta_file_name.c_str());

		string candidate_junction_fasta_file_name = settings.candidate_junction_fasta_file_name;

		cFastaFile out(candidate_junction_fasta_file_name, ios_base::out);

		// hash by junction sequence concatenated with name of counts
		map_t candidate_junctions;

		// summary data for this step
//		my $hcs;
//
//		my @read_files = $settings->read_files;
		int32_t i = 0;

		for (int32_t j = 0; j < settings.read_structures.size(); j++)
		{
			Settings::ReadStructure read_struct = settings.read_structures[j];

			string read_file = read_struct.base_name;
			cerr << "  READ FILE::" << read_file << endl;

			// Zero out summary information
			map<string, int32_t> s;
			s["reads_with_unique_matches"] = 0;
			s["reads_with_non_unique_matches"] = 0;
			s["reads_with_redundant_matches"] = 0;
			s["reads_with_possible_hybrid_matches"] = 0;
			s["reads_with_more_multiple_matches"] = 0;
			s["reads_with_only_unwanted_matches"] = 0;
			s["reads_with_best_unwanted_matches"] = 0;
			//s["reads_with_no_matches"] = 0; #calculate as total minus number found
			s["reads_with_ambiguous_hybrids"] = 0;

			// Decide which input SAM file we are using...
			string reference_sam_file_name = settings.preprocess_junction_split_sam_file_name;
			if (settings.candidate_junction_score_method.compare("POS_HASH") != 0)
				reference_sam_file_name = settings.reference_sam_file_name;

			tamFile tam = sam_open(reference_sam_file_name.c_str()); // or die("Could not open reference same file: $reference_sam_file_name");
			bam_header_t* header = sam_header_read2(reference_faidx_file_name.c_str()); // or die("Error reading reference fasta index file: $reference_faidx_file_name");
			bam1_t* last_alignment;
			vector<bam1_t*> al_ref;

			while (true)
			{
				// resolve_alignments.c
				al_ref = tam_next_read_alignments(tam, header, last_alignment, false);

				if (al_ref.size() == 0)
					break;

				if (++i % 10000 == 0)
					cerr << "    ALIGNED READ:" << i << " CANDIDATE JUNCTIONS:" << candidate_junctions.size() << endl;

				// for testing...
				if (settings.candidate_junction_read_limit != 0 && i > settings.candidate_junction_read_limit)
					break;

				_alignments_to_candidate_junctions(settings, summary, ref_seq_info, candidate_junctions, fai, header, al_ref);
			}
//
//			$hcs->{read_file}->{$read_file} = $s;
		}

		//
		// Calculate pos_hash score for each candidate junction now that we have the complete list of read match support
		//

//		JUNCTION_SEQ: foreach my $junction_seq (sort keys %$candidate_junctions)
//		{
//			JUNCTION_ID: foreach my $junction_id (keys %{$candidate_junctions->{$junction_seq}} )
//			{
//				my $cj = $candidate_junctions->{$junction_seq}->{$junction_id};
//				$cj->{pos_hash_score} = scalar keys %{$cj->{read_begin_hash}};
//
//				print ">>>$junction_id\n" if ($verbose);
//				print Dumper($cj->{read_begin_hash}) if ($verbose);
//			}
//		}
//
		//
		// Combine hash into a list, one item for each unique sequence (also joining reverse complements)
		//
//
//		my %observed_pos_hash_score_distribution;
//		my %observed_min_overlap_score_distribution;
//
//		my @combined_candidate_junctions;
//		my $handled_seq;
//		my %ids_to_print;
//
//		## Sort here is to get reproducible ordering.
//		JUNCTION_SEQ: foreach my $junction_seq (sort keys %$candidate_junctions)
//		{
//			print "Handling $junction_seq\n" if ($verbose);
//
//			## We may have already done the reverse complement
//			next if ($handled_seq->{$junction_seq});
//
//			my @combined_candidate_junction_list = ();  ## holds all junctions with the same seq
//			my $rc_junction_seq = Breseq::Fastq::revcom($junction_seq);
//
//			JUNCTION_ID: foreach my $junction_id (keys %{$candidate_junctions->{$junction_seq}} )
//			{
//				## add redundancy to the $junction_id
//				my $cj = $candidate_junctions->{$junction_seq}->{$junction_id};
//				$junction_id .= $Breseq::Shared::junction_name_separator . $cj->{r1} . $Breseq::Shared::junction_name_separator . $cj->{r2};
//				push @combined_candidate_junction_list, { id=>$junction_id, pos_hash_score=>$cj->{pos_hash_score}, min_overlap_score=>$cj->{min_overlap_score}, seq=>$junction_seq, rc_seq=>$rc_junction_seq };
//			}
//			$handled_seq->{$junction_seq}++;
//
//			## add the reverse complement
//			if (defined $candidate_junctions->{$rc_junction_seq})
//			{
//				JUNCTION_ID: foreach my $junction_id (keys %{$candidate_junctions->{$rc_junction_seq}} )
//				{
//					## add redundancy to the $junction_id (reversed)
//					my $cj = $candidate_junctions->{$rc_junction_seq}->{$junction_id};
//
//					$junction_id .= $Breseq::Shared::junction_name_separator . $cj->{r2} . $Breseq::Shared::junction_name_separator . $cj->{r1};
//					push @combined_candidate_junction_list, { id=>$junction_id, pos_hash_score=>$cj->{pos_hash_score}, min_overlap_score=>$cj->{min_overlap_score}, seq=>$rc_junction_seq, rc_seq=>$junction_seq };
//				}
//				$handled_seq->{$rc_junction_seq}++;
//			}
//
//			## Sort by unique coordinate, then redundant (or second unique) coordinate to get reliable ordering for output
//			sub by_score_unique_coord
//			{
//				my $a_item = Breseq::Shared::junction_name_split($a->{id});
//				my $a_uc = ($a_item->{side_1}->{redundant} != 0) ? $a_item->{side_2}->{position} : $a_item->{side_1}->{position};
//				my $a_rc = ($a_item->{side_1}->{redundant} != 0) ? $a_item->{side_1}->{position} : $a_item->{side_2}->{position};
//
//				my $b_item = Breseq::Shared::junction_name_split($b->{id});
//				my $b_uc = ($b_item->{side_1}->{redundant} != 0) ? $b_item->{side_2}->{position} : $b_item->{side_1}->{position};
//				my $b_rc = ($b_item->{side_1}->{redundant} != 0) ? $b_item->{side_1}->{position} : $b_item->{side_2}->{position};
//
//				return ( -($a->{pos_hash_score} <=> $b->{pos_hash_score}) || -($a->{min_overlap_score} <=> $b->{min_overlap_score}) || ($a_uc <=> $b_uc) || ($a_rc <=> $b_rc));
//			}
//			@combined_candidate_junction_list = sort by_score_unique_coord @combined_candidate_junction_list;
//			my $best_candidate_junction = $combined_candidate_junction_list[0];
//
//			## Save the score in the distribution
//			Breseq::Shared::add_score_to_distribution(\%observed_pos_hash_score_distribution, $best_candidate_junction->{pos_hash_score});
//			Breseq::Shared::add_score_to_distribution(\%observed_min_overlap_score_distribution, $best_candidate_junction->{min_overlap_score});
//
//			## Check minimum requirements
//			next JUNCTION_SEQ if ($best_candidate_junction->{pos_hash_score} < $settings->{minimum_candidate_junction_pos_hash_score});
//			next JUNCTION_SEQ if ($best_candidate_junction->{min_overlap_score} < $settings->{minimum_candidate_junction_min_overlap_score});
//
//			## Make sure it isn't a duplicate junction id -- this should NEVER happen and causes downstream problem.
//			## <--- Begin sanity check
//			if ($ids_to_print{$best_candidate_junction->{id}})
//			{
//				print STDERR "Attempt to create junction candidate with duplicate id: $combined_candidate_junction_list[0]->{id}\n";
//				print STDERR "==Existing junction==\n";
//				print STDERR Dumper($ids_to_print{$best_candidate_junction->{id}});
//				print STDERR "==New junction==\n";
//				print STDERR Dumper($best_candidate_junction);
//
//				die if ($best_candidate_junction->{seq} ne $ids_to_print{$best_candidate_junction->{id}}->{seq});
//				next JUNCTION_SEQ;
//			}
//			$ids_to_print{$best_candidate_junction->{id}} = $best_candidate_junction;
//			## <--- End sanity check
//
//			push @combined_candidate_junctions, $best_candidate_junction;
//		}
//
//		@combined_candidate_junctions = sort {-($a->{pos_hash_score} <=> $b->{pos_hash_score}) || -($a->{min_overlap_score} <=> $b->{min_overlap_score}) || (length($a->{seq}) <=> length($b->{seq}))} @combined_candidate_junctions;
//		print Dumper(@combined_candidate_junctions) if ($verbose);
//
//		###
//		## Limit the number of candidate junctions that we print by:
//		##   (1) A maximum number of candidate junctions
//		##   (2) A maximum length of the sequences in candidate junctions
//		###
//
		cerr << "  Taking top candidate junctions..." << endl;
//
//		## adding up the lengths might be too time-consuming to be worth it...
//		my $total_cumulative_cj_length = 0;
//		my $total_candidate_junction_number = scalar @combined_candidate_junctions;
//		foreach my $c (@combined_candidate_junctions)
//		{
//			$total_cumulative_cj_length += length $c->{seq};
//		}
//
//		my @duplicate_sequences;
//		my $cumulative_cj_length = 0;
//		my $lowest_accepted_pos_hash_score = 'NA';
//		my $lowest_accepted_min_overlap_score = 'NA';
//
//		## Right now we limit the candidate junctions to have a length no longer than the reference sequence.
//		my $cj_length_limit = int($summary->{sequence_conversion}->{total_reference_sequence_length} * $settings->{maximum_candidate_junction_length_factor});
//		my $maximum_candidate_junctions = $settings->{maximum_candidate_junctions};
//		my $minimum_candidate_junctions = $settings->{minimum_candidate_junctions};
//
//		print STDERR sprintf ("  Minimum number to keep: %7d \n", $minimum_candidate_junctions);
//		print STDERR sprintf ("  Maximum number to keep: %7d \n", $maximum_candidate_junctions);
//		print STDERR sprintf ("  Maximum length to keep: %7d bases\n", $cj_length_limit);
//
//		print STDERR "    Initial: Number = " . $total_candidate_junction_number . ", Cumulative Length = " . $total_cumulative_cj_length . " bases\n";
//
//		if ((defined $settings->{maximum_candidate_junctions}) && (@combined_candidate_junctions > 0))
//		{
//			my @remaining_ids = ();
//			my @list_in_waiting = ();
//			my $add_cj_length = 0;
//			my $num_duplicates = 0;
//
//			my $i = 0;
//			my $current_pos_hash_score = $combined_candidate_junctions[$i]->{pos_hash_score};
//			my $current_min_overlap_score = $combined_candidate_junctions[$i]->{min_overlap_score};
//
//			## Check to make sure that adding everything from the last iteration doesn't put us over any limits...
//			my $new_number = scalar(@remaining_ids) + scalar(@list_in_waiting);
//			my $new_length = $cumulative_cj_length + $add_cj_length;
//			while (	( $new_number <= $minimum_candidate_junctions ) || (($new_length <= $cj_length_limit) && ($new_number <= $maximum_candidate_junctions)) )
//			{
//				## OK, add everything from the last iteration
//				$cumulative_cj_length += $add_cj_length;
//				push @remaining_ids, @list_in_waiting;
//				$lowest_accepted_pos_hash_score = $current_pos_hash_score;
//				$lowest_accepted_min_overlap_score = $current_min_overlap_score;
//
//				## Zero out what we will add
//				$add_cj_length = 0;
//				@list_in_waiting = ();
//				$num_duplicates = 0;
//
//				## Check to make sure we haven't exhausted the list
//				last if ($i >= scalar @combined_candidate_junctions);
//
//				$current_pos_hash_score = $combined_candidate_junctions[$i]->{pos_hash_score};
//				$current_min_overlap_score = $combined_candidate_junctions[$i]->{min_overlap_score};
//				CANDIDATE: while (
//					    ($i < scalar @combined_candidate_junctions)
//				     && ($combined_candidate_junctions[$i]->{pos_hash_score} == $current_pos_hash_score)
//				     && ($combined_candidate_junctions[$i]->{min_overlap_score} == $current_min_overlap_score)
//					)
//				{
//					my $c = $combined_candidate_junctions[$i];
//					push @list_in_waiting, $c;
//					$add_cj_length += length $c->{seq};
//
//				} continue {
//					$i++;
//				}
//
//				$new_number = scalar(@remaining_ids) + scalar(@list_in_waiting);
//				$new_length = $cumulative_cj_length + $add_cj_length;
//
//				print STDERR sprintf("      Testing Pos Hash Score = %4d, Min Overlap Score = %4d, Num = %6d, Length = %6d\n", $current_pos_hash_score, $current_min_overlap_score, (scalar @list_in_waiting), $add_cj_length);
//
//			}
//			@combined_candidate_junctions = @remaining_ids;
//		}
//
//		my $accepted_candidate_junction_number = scalar @combined_candidate_junctions;
//		print STDERR "    Accepted: Number = $accepted_candidate_junction_number, Pos Hash Score >= $lowest_accepted_pos_hash_score, Min Overlap Score >= $lowest_accepted_min_overlap_score, Cumulative Length = $cumulative_cj_length bases\n";
//
//		## Save summary statistics
//		$hcs->{total}->{number} = $total_candidate_junction_number;
//		$hcs->{total}->{length} = $total_cumulative_cj_length;
//
//		$hcs->{accepted}->{number} = $accepted_candidate_junction_number;
//		$hcs->{accepted}->{length} = $cumulative_cj_length;
//		$hcs->{accepted}->{pos_hash_score_cutoff} = $lowest_accepted_pos_hash_score;
//		$hcs->{accepted}->{min_overlap_score_cutoff} = $lowest_accepted_min_overlap_score;
//
//		$hcs->{pos_hash_score_distribution} = \%observed_pos_hash_score_distribution;
//		$hcs->{min_overlap_score_distribution} = \%observed_min_overlap_score_distribution;
//
//		###
//		## Print out the candidate junctions, sorted by the lower coordinate, higher coord, then number
//		###
//		sub by_ref_seq_coord
//		{
//			my $acj = Breseq::Shared::junction_name_split($a->{id});
//			my $bcj = Breseq::Shared::junction_name_split($b->{id});
//			return (($ref_seq_info->{seq_order}->{$acj->{side_1}->{seq_id}} <=> $ref_seq_info->{seq_order}->{$bcj->{side_1}->{seq_id}}) ||  ($acj->{side_1}->{position} <=> $bcj->{side_1}->{position}));
//		}
//
//		@combined_candidate_junctions = sort by_ref_seq_coord @combined_candidate_junctions;
//
//		foreach my $junction (@combined_candidate_junctions)
//		{
//			#print Dumper($ids_to_print);
//
//			my $seq = Bio::Seq->new(
//				-display_id => $junction->{id}, -seq => $junction->{seq});
//			$out->write_seq($seq);
//		}
//
//		## create SAM faidx
//		my $samtools = $settings->ctool('samtools');
//		Breseq::Shared::system("$samtools faidx $candidate_junction_fasta_file_name") if (scalar @combined_candidate_junctions > 0);
//
//		$summary->{candidate_junction} = $hcs;
	}
  
} // namespace breseq
