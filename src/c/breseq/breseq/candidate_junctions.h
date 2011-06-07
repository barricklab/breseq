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

#include "breseq/resolve_alignments.h"

#include "breseq/common.h"

using namespace std;

namespace breseq {

	typedef std::map<string, boost::any> map_t;

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
			int32_t flanking_left;
			int32_t flanking_right;
		};

		struct PassedPair
		{
			bam1_t*	a1;
			bam1_t* a2;
			int32_t union_length;
			int32_t a1_unique_length;
			int32_t a2_unique_length;
		};

		/*! Preprocesses alignments
		 */
		static void preprocess_alignments(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info);

		/*! Predicts candidate junctions
		 */
		static void identify_candidate_junctions(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info);

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

		static bool _alignments_to_candidate_junction(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info, faidx_t* fai, bam_header_t* header, bam1_t* a1, bam1_t* a2,
														int32_t& redundancy_1, int32_t& redundancy_2, string& junction_seq_string, string& ref_seq_matched_1, string& ref_seq_matched_2, string& junction_coord_1, string& junction_coord_2, int32_t& read_begin_coord, JunctionList& junction_id_list);
		static void _alignments_to_candidate_junctions(Settings settings, Summary summary,  const cReferenceSequences& ref_seq_info, map<string, map<string, CandidateJunction> >& candidate_junctions, faidx_t* fai, bam_header_t* header, vector<bam1_t*> al_ref);
		static bool _check_read_pair_requirements(Settings settings, int32_t a1_start, int32_t a1_end, int32_t a2_start, int32_t a2_end, int32_t& a1_unique_length, int32_t& a2_unique_length, int32_t& union_length);
		static void _entire_read_matches(map_t a);
		static void _num_matches_from_end(bam1_t* a, string refseq_str, bool dir, int32_t overlap, int32_t& qry_mismatch_pos, int32_t& ref_mismatch_pos);
		static void _split_indel_alignments(Settings settings, Summary summary, bam_header_t* header, ofstream& PSAM, int32_t min_indel_split_len, vector<bam1_t*> al_ref);
		static void _by_ref_seq_coord(map_t a, map_t b, map_t ref_seq_info);
		static void _by_score_unique_coord(map_t a, map_t b);
		static void _tam_write_split_alignment(map_t fh, map_t header, map_t min_indel_split_len, map_t a);

	}; // class CandidateJunction

} // namespace breseq

#endif
