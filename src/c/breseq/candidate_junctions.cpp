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

	CandidateJunctions::CandidateJunctions() {}

	bool CandidateJunctions::_alignments_to_candidate_junction(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info, faidx_t* fai, bam_header_t* header, bam1_t* a1, bam1_t* a2,
															  int32_t& redundancy_1, int32_t& redundancy_2, string& junction_seq_string, string& ref_seq_matched_1, string& ref_seq_matched_2, string& junction_coord_1, string& junction_coord_2, int32_t& read_begin_coord, JunctionList& junction_id_list)
	{
		bool verbose = false;

		// set up local settings
		int32_t flanking_length = settings.max_read_length;
//		my $reference_sequence_string_hash_ref = $ref_seq_info->{ref_strings};

		// Method
		//
		// Hash junctions by a key showing the inner coordinate of the read.
		// and the direction propagating from that position in the reference sequence.
		// Prefer the unique or lower coordinate side of the read for main hash.
		//
		// REL606__1__1__REL606__4629812__0__0__
		// means the junction sequence is 36-1 + 4629812-4629777 from the reference sequence
		//
		// On the LEFT side:  0 means this is highest coord of alignment, junction seq begins at lower coord
		//                    1 means this is lowest coord of alignment, junction seq begins at higher coord
		// On the RIGHT side: 0 means this is highest coord of alignment, junction seq continues to lower coord
		//                    1 means this is lowest coord of alignment, junction seq continues to higher coord
		//
		// Note that there are more fields now than shown in this example...
		//
		// Need the junction key to include the offset to get to the junction within the read for cases where
		// the junction is near the end of the sequence... test case exists in JEB574.

//		my %printed_keys;
		int32_t i = 0;

		string read_id = bam1_qname(a1);

		// First, sort matches by their order in the query
		bam1_t* q1 = a1;
		bam1_t* q2 = a2;
		int32_t q1_start, q1_end;
		int32_t q2_start, q2_end;
		alignment_query_start_end(q1, q1_start, q1_end);
		alignment_query_start_end(q2, q2_start, q2_end);

		if (verbose)
			cout << q1_start << ", " << q1_end << ", " << q2_start << ", " << q2_end << endl;

		// Reverse the coordinates to be consistently such that 1 refers to lowest...
		if (q2_start < q1_start)
		{
			swap(q1_start, q2_start);
			swap(q1_end, q2_end);
			swap(redundancy_1, redundancy_2);
			swap(q1, q2);
		}

		// create hash key and store information about the location of this hit

		cReferenceSequences req_seq_info_copy = ref_seq_info;
		uint32_t seq_id;

		bool hash_strand_1 = is_reversed(q1);
		string hash_seq_id_1 = header->target_name[q1->core.tid];
		seq_id = req_seq_info_copy.seq_id_to_index(hash_seq_id_1);
		string ref_seq_1 = ref_seq_info[seq_id].m_fasta_sequence.m_sequence;

		bool hash_strand_2 = !is_reversed(q2);
		string hash_seq_id_2 = header->target_name[q2->core.tid];
		seq_id = req_seq_info_copy.seq_id_to_index(hash_seq_id_2);
		string ref_seq_2 = ref_seq_info[seq_id].m_fasta_sequence.m_sequence;

		// how much overlap is there between the two matches?
		// positive if the middle sequence can match either side of the read
		// negative if there is sequence in the read NOT matched on either side
		int32_t overlap = -1 * (q2_start - q1_end - 1);

		//
		// OVERLAP MISMATCH CORRECTION
		//
		// If there are mismatches in one or the other read in the overlap region
		// then we need to adjust the coordinates. Why? All sequences that we
		// retrieve are from the reference sequence and there are two choices
		// for where to extract this non-necessarily identical sequence!

		// save these as variables, because we may have to adjust them
		int32_t r1_start = q1_start;
		int32_t r1_end = q1_end;
		int32_t r2_start = q2_start;
		int32_t r2_end = q2_end;

		if (verbose)
		{
			ref_seq_matched_1 = ref_seq_1.substr(r1_start - 1, r1_end - r1_start + 1);
			ref_seq_matched_2 = ref_seq_2.substr(r2_start - 1, r2_end - r2_start + 1);

			cout << "==============> Initial Matches" << endl;
			cout << "Alignment #1" << endl;
			cout << "qpos: " << q1_start << "-" << q1_end << " rpos: " << r1_start << "-" << r1_end << " reversed: " << is_reversed(q1) << endl;
			cout << bama_qseq(q1) << endl << ref_seq_matched_1 << endl;
//			print Dumper($q1->cigar_array);

			cout << "Alignment #2" << endl;
			cout << "qpos: " << q2_start << "-" << q2_end << " rpos: " << r2_start << "-" << r2_end << " reversed: " << is_reversed(q2) << endl;
			cout << bama_qseq(q2) << endl << ref_seq_matched_2 << endl;
//			print Dumper($q2->cigar_array);
			cout << "<==============" << endl;

			// debug print information
			cout << "=== overlap: " << overlap << endl;
		}

		// Adjust the overlap in cases where there is a mismatch within the overlap region
		if (overlap > 0)
		{
			int32_t q1_move, q2_move, r1_move, r2_move;
			_num_matches_from_end(q1, ref_seq_1, false, overlap, q1_move, r1_move);
			_num_matches_from_end(q2, ref_seq_2, true, overlap, q2_move, r2_move);

			if (q1_move >= 0)
			{
				if (verbose)
					cout << "ALIGNMENT #1 OVERLAP MISMATCH: " << q1_move << ", " << r1_move << endl;
				// change where it ENDS
				q1_end -= q1_move;
				if (is_reversed(q1))
					r1_start += r1_move;
				else
					r1_end -= r1_move;
			}

			if (q2_move >= 0)
			{
				if (verbose)
					cout << "ALIGNMENT #2 OVERLAP MISMATCH: " << q2_move << ", " << r2_move << endl;
				// change where it STARTS
				q2_start += q2_move;
				if (!is_reversed(q2))
					r2_start += r2_move;
				else
					r2_end -= r2_move;
			}

			// We need to re-check that they have the required amount of unique length if they moved
			if (q1_move >= 0 || q2_move >= 0)
			{
				if (verbose)
					cout << "Rejecting alignment because there is not not enough overlap" << endl;
				int32_t dont_care;
				bool passed = _check_read_pair_requirements(settings, q1_start, q1_end, q2_start, q2_end, dont_care, dont_care, dont_care);
				if (!passed) return false;
			}

			//re-calculate the overlap
			overlap = -1 * (q2_start - q1_end - 1);
			if (verbose)
				cout << "=== overlap corrected for mismatches " << overlap << endl;
		}

		// create hash coords AFTER overlap adjustment
		int32_t hash_coord_1 = (hash_strand_1) ? r1_start : r1_end;
		int32_t hash_coord_2 = (hash_strand_2) ? r2_start : r2_end;

		// these are the positions of the beginning and end of the read, across the junction
		// query 1 is the start of the read, which is why we hash by this coordinate
		// (it is less likely to be shifted by a nucleotide or two by base errors)
		read_begin_coord = (hash_strand_1) ? r1_end : r1_start;

		if (verbose)
		{
			string ref_seq_matched_1 = ref_seq_1.substr(r1_start - 1, r1_end - r1_start + 1);
			string ref_seq_matched_2 = ref_seq_2.substr(r2_start - 1, r2_end - r2_start + 1);

			cout << "==============> Initial Matches" << endl;
			cout << "Alignment #1" << endl;
			cout << "qpos: " << q1_start << "-" << q1_end << " rpos: " << r1_start << "-" << r1_end << " reversed: " << is_reversed(q1) << endl;
			cout << bama_qseq(q1) << endl << ref_seq_matched_1 << endl;
//			print Dumper($q1->cigar_array);

			cout << "Alignment #2" << endl;
			cout << "qpos: " << q2_start << "-" << q2_end << " rpos: " << r2_start << "-" << r2_end << " reversed: " << is_reversed(q2) << endl;
			cout << bama_qseq(q2) << endl << ref_seq_matched_2 << endl;
//			print Dumper($q2->cigar_array);
			cout << "<==============" << endl;
		}

		// Calculate an offset that only applies if the overlap is positive (sequence is shared between the two ends)
		int32_t overlap_offset = (overlap > 0) ? overlap : 0;
		if (verbose)
			cout << "Overlap offset: " << overlap_offset << endl;

		// record what parts of the reference sequence were actually matched on each side
		// this is to determine whether that side was redundant or unique in the reference

		int32_t ref_seq_matched_length_1 = r1_end - r1_start + 1;
		int32_t ref_seq_matched_length_2 = r2_end - r2_start + 1;

		// create the sequence of the candidate junction
		junction_seq_string = "";

		// first end
		int32_t flanking_left = flanking_length;
		if (!hash_strand_1) // alignment is not reversed
		{
			// start_pos is in 1-based coordinates
			int32_t start_pos = hash_coord_1 - (flanking_left - 1) - overlap_offset;
			if (start_pos < 1)
			{
				if (verbose)
					cout << "START POS 1: " << start_pos << " < 0" << endl;
				flanking_left += start_pos - 1;
				start_pos = 1;
			}
			string add_seq = ref_seq_1.substr(start_pos - 1, flanking_left + overlap_offset);
			if (verbose) cout << "1F: " << add_seq << endl;
			junction_seq_string += add_seq;

			ref_seq_matched_1 = add_seq.substr(add_seq.size() - overlap_offset - ref_seq_matched_length_1, ref_seq_matched_length_1);
		}
		else
		{
			// end_pos is in 1-based coordinates
			int32_t end_pos = hash_coord_1 + (flanking_left - 1) + overlap_offset;
			if (end_pos > ref_seq_1.size())
			{
				if (verbose) cout << "END POS 1: (" << end_pos << " < length" << endl;
				flanking_left -= end_pos - ref_seq_1.size();
				end_pos = ref_seq_1.size();
			}

			string add_seq = ref_seq_1.substr(end_pos - (flanking_left + overlap_offset), flanking_left + overlap_offset);
			add_seq = reverse_complement(add_seq);
			if (verbose) cout << "1R: " << add_seq << endl;
			junction_seq_string += add_seq;

			ref_seq_matched_1 = add_seq.substr(add_seq.size() - overlap_offset - ref_seq_matched_length_1, ref_seq_matched_length_1);
		}

		// Add any unique junction sequence that was only in the read
		// and NOT present in the reference genome
		string unique_read_seq_string = "";
		if (overlap < 0)
			unique_read_seq_string = bama_qseq(q1).substr(q1_end, -1 * overlap);
		junction_seq_string += unique_read_seq_string;

		if (verbose) cout << "+U: " << unique_read_seq_string << endl;

		// second end
		int32_t flanking_right = flanking_length;
		if (hash_strand_2) //alignment is not reversed
		{
			// end_pos is in 1-based coordinates
			int32_t end_pos = hash_coord_2 + (flanking_right - 1) + overlap_offset;
			if (end_pos > ref_seq_2.size())
			{
				if (verbose)
					cout << "END POS 2: (" << end_pos << " < length" << endl;
				flanking_right -= (end_pos - ref_seq_2.size());
				end_pos = ref_seq_2.size();
			}
			string add_seq = ref_seq_2.substr(end_pos - flanking_right, flanking_right);
			if (verbose) cout << "2F: " << add_seq << endl;
			junction_seq_string += add_seq;

			ref_seq_matched_2 = add_seq.substr(add_seq.size() - ref_seq_matched_length_2, ref_seq_matched_length_2);
		}
		else // ($m_2->{hash_strand} * $m_2->{read_side} == -1)
		{
			// start_pos is in 1-based coordinates
			int32_t start_pos = hash_coord_2 - (flanking_right - 1) - overlap_offset;
			if (start_pos < 1)
			{
				if (verbose) cout << "START POS 2: " << start_pos << " < 0" << endl;
				flanking_right += start_pos - 1;
				start_pos = 1;
			}
			string add_seq = ref_seq_2.substr(start_pos - 1, flanking_right);
			add_seq = reverse_complement(add_seq);
			if (verbose) cout << "2R: " << add_seq << endl;
			junction_seq_string += add_seq;

			ref_seq_matched_2 = add_seq.substr(add_seq.size() - ref_seq_matched_length_2, ref_seq_matched_length_2);
		}

		if (verbose) cout << "3: " << junction_seq_string << endl;

		// create hash coords after this adjustment
		if (!hash_strand_1)
			r1_end -= overlap_offset;
		else //reversed
			r1_start += overlap_offset;

		if (!hash_strand_2)
			r2_end -= overlap_offset;
		else //reversed
			r2_start += overlap_offset;

		// matched reference sequence
		ref_seq_matched_1 = ref_seq_1.substr(r1_start - 1, r1_end - r1_start + 1);
		ref_seq_matched_2 = ref_seq_2.substr(r2_start - 1, r2_end - r2_start + 1);

		// want to be sure that lowest ref coord is always first for consistency
		if ( hash_seq_id_1.compare(hash_seq_id_2) > 0 || (hash_seq_id_1.compare(hash_seq_id_2) == 0 && hash_coord_2 < hash_coord_1) )
		{
			swap(hash_coord_1, hash_coord_2);
			swap(hash_strand_1, hash_strand_2);
			swap(hash_seq_id_1, hash_seq_id_2);
			swap(redundancy_1, redundancy_2);
			swap(flanking_left, flanking_right);
			swap(ref_seq_matched_1, ref_seq_matched_2);

			junction_seq_string = reverse_complement(junction_seq_string);
			unique_read_seq_string = reverse_complement(unique_read_seq_string);
		}

		JunctionList::Side side_1 = {
			hash_seq_id_1,	// seq_id
			hash_coord_1,	// position
			hash_strand_1,	// strand
			-1				// redundant: set to magic (uninitialized) value
		}, side_2 = {
			hash_seq_id_2,	// seq_id
			hash_coord_2,	// position
			hash_strand_2,	// strand
			-1				// redundant: set to magic (uninitialized) value
		};
		JunctionList new_junction =
		{
			side_1,
			side_2,
			overlap, 				// alignment_overlap
			unique_read_seq_string,	// unique_read_sequence
			flanking_left,
			flanking_right
		}; junction_id_list = new_junction;

		string join1[] = { hash_seq_id_1, boost::lexical_cast<string>(hash_coord_1), boost::lexical_cast<string>(hash_strand_1) };
		junction_coord_1 = join(join1, "::");
		string join2[] = { hash_seq_id_2, boost::lexical_cast<string>(hash_coord_2), boost::lexical_cast<string>(hash_strand_2) };
		junction_coord_2 = join(join2, "::");

		string junction_id = junction_name_join(junction_id_list);
		if (verbose)
		{
			cout << "READ ID: " << bam1_qname(a1) << endl;
			cout << "JUNCTION ID: " << junction_id << endl;
		}

		assert(junction_seq_string.size() > 0); // die "Junction sequence not found: $junction_id " . $q1->qname . " " . $a2->qname  if (!$junction_seq_string);
		assert(junction_seq_string.size() != flanking_left + flanking_right + abs(overlap)); // die "Incorrect length for $junction_seq_string: $junction_id " . $q1->qname . " " . $a2->qname if (length $junction_seq_string != $flanking_left + $flanking_right + abs($overlap));

		return true;
	}

	void CandidateJunctions::_alignments_to_candidate_junctions(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info, map<string, map<string, CandidateJunction> >& candidate_junctions, faidx_t* fai, bam_header_t* header, vector<bam1_t*> al_ref)
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
		map<string, map<string, int32_t> > redundant_junction_sides;
		vector<JunctionListContainer> junctions;

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

		for (int32_t i = 0; i < passed_pair_list.size(); i++)
		{
			PassedPair pp = passed_pair_list[i];

			// localize variables
			bam1_t* a1 = pp.a1;
			bam1_t* a2 = pp.a2;
			int32_t a1_unique_length = pp.a1_unique_length;
			int32_t a2_unique_length = pp.a2_unique_length;

			// start redundancy at 1
			int32_t r1 = 1;
			int32_t r2 = 1;

			// we pass back and forth the redundancies in case they switch sides
			string junction_seq_string;
			string side_1_ref_seq, side_2_ref_seq;
			string junction_coord_1, junction_coord_2;
			int32_t read_begin_coord;
			JunctionList junction_id_list;

			bool passed = _alignments_to_candidate_junction(settings, summary, ref_seq_info, fai, header, a1, a2,
															r1, r2, junction_seq_string, side_1_ref_seq, side_2_ref_seq, junction_coord_1, junction_coord_2, read_begin_coord, junction_id_list);
			if (!passed) continue;

			// a value of zero gets added if they were unique, >0 if redundant b/c they matched same reference sequence
			redundant_junction_sides[side_1_ref_seq][junction_coord_1] += r1 - 1;
			redundant_junction_sides[side_2_ref_seq][junction_coord_2] += r2 - 1;

			// also add to the reverse complement, because we can't be sure of the strandedness
			// (alternately we could reverse complement if the first base was an A or C, for example
			// it seems like there could possibly be some cross-talk between sides of junctions here, that
			/// could snarl things up, but I'm not sure?
			/// TO DO: I'm too tired of this section to do it now, but the correct sequence strand could be
			/// decided (keeping track of all the reversals) in _alignments_to_candidate_junction
			redundant_junction_sides[reverse_complement(side_1_ref_seq)][junction_coord_1] += r1 - 1;
			redundant_junction_sides[reverse_complement(side_2_ref_seq)][junction_coord_2] += r2 - 1;

			int32_t min_overlap_score = min(a1_unique_length, a2_unique_length);

			JunctionListContainer junction =
			{
				junction_id_list,		// list
				junction_seq_string,	// str
				min_overlap_score,
				read_begin_coord,
				side_1_ref_seq,
				side_2_ref_seq
			};
			junctions.push_back(junction);
		}

		if (verbose)
			cout << "  Junctions: " << junctions.size() << endl;
		//#	print "JACKPOT!!!\n" if (scalar @junctions > 5);

		// Done if everything already ruled out...
		if (junctions.size() == 0) return;

		if (verbose)
			cout << bam1_qname(al_ref[0]) << endl;

		// only now that we've looked through everything can we determine whether the reference sequence matched
		// on a side was unique, after correcting for overlap
		for (int32_t i = 0; i < junctions.size(); i++)
		{
			JunctionListContainer jct = junctions[i];
			JunctionList junction_id_list = jct.list;
			string junction_seq_string = jct.str;
			int32_t min_overlap_score = jct.min_overlap_score;
			int32_t read_begin_coord = jct.read_begin_coord;

			string side_1_ref_seq = jct.side_1_ref_seq;
			string side_2_ref_seq = jct.side_2_ref_seq;

			// these are the total number of redundant matches to that side of the junction
			// the only way to be unique is to have at most one coordinate corresponding to that sequence
			// and not have that reference sequence redundantly matched (first combination of alignment coords)

			uint32_t total_r1 = 0, total_r2 = 0;
			map<string, int32_t>::iterator it;
			map<string, int32_t> side = redundant_junction_sides[side_1_ref_seq];

			for (it = side.begin(); it != side.end(); it++)
				total_r1 += (*it).second + 1;

			side = redundant_junction_sides[side_2_ref_seq];

			for (it = side.begin(); it != side.end(); it++)
				total_r2 += (*it).second + 1;

			string junction_id = junction_name_join(junction_id_list);
			if (verbose) cout << junction_id << endl;

			// initialize candidate junction if it didn't exist
			// they are redundant by default, until proven otherwise
			CandidateJunction cj;

			if (candidate_junctions.count(junction_seq_string) > 0 && candidate_junctions[junction_seq_string].count(junction_id) > 0)
			{
				cj = candidate_junctions[junction_seq_string][junction_id];
			}
			else
			{
				//redundancies of each side
				cj.r1 = cj.r2 = 1;

				//maximum nonoverlapping match size on each side
				cj.L1 = cj.L2 = 0;
			}

			// Update score of junction and the redundancy of each side
			cj.min_overlap_score += min_overlap_score;
			cj.read_begin_hash[read_begin_coord]++;

			if (verbose)
			{
				cout << "Totals: " << total_r1 << ", " << total_r2 << endl;
				cout << "Redundancy (before): " << cj.r1 << " (" << cj.L1 << ") " << cj.r2 << " (" << cj.L2 << ")" << endl;
			}
			int32_t side_1_ref_match_length = side_1_ref_seq.size();
			int32_t side_2_ref_match_length = side_2_ref_seq.size();

			// Longer reads into a side trump the redundancy of shorter reads for two reasons
			//   (1) obviously, the longer the read, the better the chance it is unique
			//   (2) subtly, if the read barely has enough to match both sides of the junction,
			//       there are situations where you can predict uniqueness because a short match
			//       maps one place only with the overlap included, and the non-overlapping part is
			//       unique, but once you have a longer match, you see that the non-overlapping
			//       part was really redundant.
			if (side_1_ref_match_length > cj.L1)
			{
				cj.L1 = side_1_ref_match_length;
				cj.r1 = (total_r1 > 1) ? 1 : 0;
			}

			if (side_2_ref_match_length > cj.L2)
			{
				cj.L2 = side_1_ref_match_length; //TODO: Confirm this shouldn't be side_2
				cj.r2 = (total_r1 > 1) ? 1 : 0;
			}

			if (verbose)
				cout << "Redundancy (after): " << cj.r1 << " (" << cj.L1 << ") " << cj.r2 << " (" << cj.L2 << ")" << endl;

			candidate_junctions[junction_seq_string][junction_id] = cj;
		}
	}

	bool CandidateJunctions::_check_read_pair_requirements(Settings settings, int32_t a1_start, int32_t a1_end, int32_t a2_start, int32_t a2_end, int32_t& a1_unique_length, int32_t& a2_unique_length, int32_t& union_length)
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

	void CandidateJunctions::_entire_read_matches(map_t a) {}

	void CandidateJunctions::_num_matches_from_end(bam1_t* a, string refseq_str, bool dir, int32_t overlap, int32_t& qry_mismatch_pos, int32_t& ref_mismatch_pos)
	{
		bool verbose = false;

		// Read
		// start, end, length  of the matching sequence (in top strand coordinates)
		int32_t q2_start, q2_end;

		int32_t q_seq_start, q_seq_end;
		alignment_query_start_end(a, q_seq_start, q_seq_end);
		int32_t q_length = alignment_query_length(a);

		bool reversed = is_reversed(a);
		// sequence of the matching part of the query (top genomic strand)
		string a_qseq = bama_qseq(a);
		string q_str = a_qseq.substr(q_seq_start - 1, q_seq_end - q_seq_start + 1);
		vector<string> q_array;
		boost::split(q_array, q_str, boost::is_any_of("/"));

		// Reference
		// start, end, length of match in reference
		int32_t r_start = a->core.pos + 1;
		int32_t r_end = bam_calend(&a->core, bam1_cigar(a));
		int32_t r_length = r_end - r_start + 1;

		// sequence of match in reference (top genomic strand)
		string r_str = refseq_str.substr(r_start - 1, r_length);
		vector<string> r_array;
		boost::split(r_array, r_str, boost::is_any_of("/"));

		if (verbose)
		{
			cout << "====> Num Matches from End" << endl;
			cout << bam1_qname(a) << endl;
			cout << "direction: " << dir << endl;
			cout << "Read sequence: " << a_qseq << endl;
			cout << "Read Match coords: " << q_seq_start << "-" << q_seq_end << " " << reversed << endl;
			cout << "Read Match sequence: " << q_str << endl;
			cout << "Reference Match coords: " << r_start << "-" << r_end << endl;
			cout << "Reference Match sequence: " << r_str << endl;
		}

		if (reversed) dir = !dir;
		if (!dir)
		{
			reverse(q_array.begin(), q_array.end());
			reverse(r_array.begin(), r_array.end());
		}

		uint32_t* cigar_array = bam1_cigar(a);
		vector<uint32_t> cigar_array_ref;
		for (int32_t i = 0; i <= a->core.n_cigar; i++)
			cigar_array_ref.push_back(cigar_array[i]);

		char op_0 = cigar_array_ref[0] & BAM_CIGAR_MASK;
		char op_last = cigar_array_ref[cigar_array_ref.size()-1] & BAM_CIGAR_MASK;

		//remove soft padding
		if (op_0 == 'S') cigar_array_ref.erase(cigar_array_ref.begin());
		if (op_last == 'S') cigar_array_ref.erase(cigar_array_ref.end());
		if (!dir) reverse(cigar_array_ref.begin(), cigar_array_ref.end());

		qry_mismatch_pos = -1;
		ref_mismatch_pos = -1;

		int32_t r_pos = 0;
		int32_t q_pos = 0;
		while ((q_pos < overlap) && (r_pos < r_length) && (q_pos < q_length))
		{
			// get rid of previous items...
			uint32_t len_0 = cigar_array_ref[0] >> BAM_CIGAR_SHIFT;
			if (len_0 == 0) cigar_array_ref.erase(cigar_array_ref.begin());

			len_0 = cigar_array_ref[0] >> BAM_CIGAR_SHIFT;
			len_0--;
			len_0 <<= BAM_CIGAR_SHIFT;
			cigar_array_ref[0] &= (len_0 | BAM_CIGAR_MASK);

			// handle indels
			op_0 = cigar_array_ref[0] & BAM_CIGAR_MASK;
			if (op_0 == 'I')
			{
				q_pos++;
				ref_mismatch_pos = r_pos;
				qry_mismatch_pos = q_pos;
			}
			else if (op_0 == 'D')
			{
				r_pos++;
				ref_mismatch_pos = r_pos;
				qry_mismatch_pos = q_pos;
			}
			else
			{
				if (q_array[q_pos] != r_array[r_pos])
				{
					ref_mismatch_pos = r_pos;
					qry_mismatch_pos = q_pos;
				}
				r_pos++;
				q_pos++;
			}

			if (verbose)
				cout << r_pos << " " << q_pos << endl;
		}

		// make 1 indexed...
		if (qry_mismatch_pos >= 0) qry_mismatch_pos++;
		if (ref_mismatch_pos >= 0) ref_mismatch_pos++;

		if (verbose)
		{
			if (qry_mismatch_pos >= 0)
				cout << "Query Mismatch At = " << qry_mismatch_pos << " | Overlap = " << overlap << endl;
			cout << "<====" << endl;
		}
	}

	void CandidateJunctions::_split_indel_alignments(Settings settings, Summary summary, bam_header_t* header, ofstream& PSAM, int32_t min_indel_split_len, vector<bam1_t*> al_ref) {}
	void CandidateJunctions::_by_ref_seq_coord(map_t a, map_t b, map_t ref_seq_info) {}
	void CandidateJunctions::_by_score_unique_coord(map_t a, map_t b) {}
	void CandidateJunctions::_tam_write_split_alignment(map_t fh, map_t header, map_t min_indel_split_len, map_t a) {}

	// Public

	/*! Preprocesses alignments
	 */
	void CandidateJunctions::preprocess_alignments(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info)
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
	void CandidateJunctions::identify_candidate_junctions(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info)
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
		map<string, map<string, CandidateJunction> > candidate_junctions;

		// summary data for this step
		struct SummaryData
		{
			struct Total
			{
				int32_t number;
				int32_t length;
			} total;

			struct Accepted
			{
				int32_t number;
				int32_t length;
				int32_t pos_hash_score_cutoff;
				int32_t min_overlap_score_cutoff;
			} accepted;

			int32_t pos_hash_score_distribution;
			int32_t min_overlap_score_distribution;

			map<string, map<string, int32_t> > read_file;
		} hcs;

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

			hcs.read_file[read_file] = s;
		}

		//
		// Calculate pos_hash score for each candidate junction now that we have the complete list of read match support
		//
		for (map<string, map<string, CandidateJunction> >::iterator outer_it = candidate_junctions.begin(); outer_it != candidate_junctions.end(); outer_it++)
		{
			for (map<string, CandidateJunction>::iterator it = (*outer_it).second.begin(); it != (*outer_it).second.end(); it++)
			{
				string junction_id = (*it).first;
				CandidateJunction cj = (*it).second;
				cj.pos_hash_score = cj.read_begin_hash.size();

				if (verbose)
				{
					cout << ">>>" << junction_id << endl; // prints junction id
					//print Dumper($cj->{read_begin_hash});
				}

			}
		}

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
