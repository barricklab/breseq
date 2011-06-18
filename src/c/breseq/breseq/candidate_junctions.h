/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_CANDIDATE_JUNCTIONS_H_
#define _BRESEQ_CANDIDATE_JUNCTIONS_H_

#include "resolve_alignments.h"
#include "alignment.h"

#include "common.h"

using namespace std;

namespace breseq {

	typedef std::map<string, string> map_t;

	class CandidateJunctions
	{
		public:

		struct CandidateJunction
		{
			int32_t r1;
			int32_t r2;
			int32_t L1;
			int32_t L2;
			int32_t min_overlap_score;
			int32_t pos_hash_score;
			map<int32_t, int32_t> read_begin_hash;

			struct Sorter {
				bool operator() (const string& a, const string& b) const { return (a < b); }
			};
		};

		struct Junction
		{
			string hash_seq_id_1;
			int32_t hash_coord_1;
			bool hash_strand_1;
			string hash_seq_id_2;
			int32_t hash_coord_2;
			bool hash_strand_2;
			int32_t overlap;
			string unique_read_seq_string;
			uint32_t flanking_left;
			uint32_t flanking_right;
		};

		struct PassedPair
		{
			alignment	a1;
			alignment a2;
			int32_t union_length;
			int32_t a1_unique_length;
			int32_t a2_unique_length;
		};

		/*! Preprocesses alignments
		 */
		static void preprocess_alignments(const Settings& settings, Summary& summary, const cReferenceSequences& ref_seq_info);

		/*! Predicts candidate junctions
		 */
		static void identify_candidate_junctions(const Settings& settings, Summary& summary, const cReferenceSequences& ref_seq_info);

		private:

		CandidateJunctions();

		struct JunctionListContainer
		{
			JunctionList list;
			string str;
			int32_t min_overlap_score;
			int32_t read_begin_coord;
			string side_1_ref_seq;
			string side_2_ref_seq;
		};
		struct CombinedCandidateJunction
		{
			string id;
			int32_t pos_hash_score;
			int32_t min_overlap_score;
			string seq;
			string rc_seq;

			// Sort by unique coordinate, then redundant (or second unique) coordinate to get reliable ordering for output
			static bool sort_by_score_unique_coord(const CombinedCandidateJunction& a, const CombinedCandidateJunction &b)
			{
				JunctionList a_item = junction_name_split(a.id);
				int32_t a_uc = a_item.side_1.position;
				int32_t a_rc = a_item.side_2.position;
				if (a_item.side_1.redundant != 0) swap(a_uc, a_rc);

				JunctionList b_item = junction_name_split(b.id);
				int32_t b_uc = b_item.side_1.position;
				int32_t b_rc = b_item.side_2.position;
				if (b_item.side_1.redundant != 0) swap(b_uc, b_rc);

				if (b.pos_hash_score != a.pos_hash_score)
					return (b.pos_hash_score < a.pos_hash_score);
				else if (b.min_overlap_score != a.min_overlap_score)
					return (b.min_overlap_score < a.min_overlap_score);
				else if (a_uc != b_uc)
					return (a_uc < b_uc);
				else
					return (a_rc < b_rc);
			}

			static bool sort_by_scores_and_seq_length(const CombinedCandidateJunction& a, const CombinedCandidateJunction &b)
			{
				if (b.pos_hash_score != a.pos_hash_score)
					return (b.pos_hash_score < a.pos_hash_score);
				else if (b.min_overlap_score != a.min_overlap_score)
					return (b.min_overlap_score < a.min_overlap_score);
				else
					return (a.seq.size() < b.seq.size());
			}

			static bool sort_by_ref_seq_coord(const CombinedCandidateJunction& a, const CombinedCandidateJunction &b)
			{
				JunctionList acj = junction_name_split(a.id);
				JunctionList bcj = junction_name_split(b.id);
				//TODO: Uncomment this code after supplying it with a ref_seq_info with a seq_order field
				/*if (ref_seq_info.seq_order[acj.side_1.seq_id] != ref_seq_info.seq_order[bcj.side_1.seq_id])
					return (ref_seq_info.seq_order[acj.side_1.seq_id] < ref_seq_info.seq_order[bcj.side_1.seq_id]);
				else*/
					return (acj.side_1.position < bcj.side_1.position);
			}
		};

		static bool _alignments_to_candidate_junction(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info, faidx_t* fai, bam_header_t* header, alignment a1, alignment a2,
														int32_t& redundancy_1, int32_t& redundancy_2, string& junction_seq_string, string& ref_seq_matched_1, string& ref_seq_matched_2, string& junction_coord_1, string& junction_coord_2, int32_t& read_begin_coord, JunctionList& junction_id_list);
		static void _alignments_to_candidate_junctions(Settings settings, Summary summary,  const cReferenceSequences& ref_seq_info, map<string, map<string, CandidateJunction>, CandidateJunction::Sorter>& candidate_junctions, faidx_t* fai, bam_header_t* header, vector<alignment> al_ref);
		static bool _check_read_pair_requirements(Settings settings, int32_t a1_start, int32_t a1_end, int32_t a2_start, int32_t a2_end, int32_t& a1_unique_length, int32_t& a2_unique_length, int32_t& union_length);
        
		static void _num_matches_from_end(alignment a, string refseq_str, bool dir, int32_t overlap, int32_t& qry_mismatch_pos, int32_t& ref_mismatch_pos);
		static void _split_indel_alignments(const Settings& settings, Summary& summary, tam_file& PSAM, int32_t min_indel_split_len, alignment_list& alignments);
		static void _by_ref_seq_coord(map_t a, map_t b, map_t ref_seq_info);
		static void _by_score_unique_coord(map_t a, map_t b);

	}; // class CandidateJunction

} // namespace breseq

#endif
