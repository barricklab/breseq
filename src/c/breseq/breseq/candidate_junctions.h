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

#include <list>
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/variant.hpp>

#include "breseq/annotated_sequence.h"
#include "breseq/common.h"

using namespace std;

namespace breseq {

	typedef std::string key_t; //!< Diff entry keys.
	typedef boost::variant<char,uint8_t,uint32_t,int,double,std::string,std::pair<int,int> > value_t; //!< Diff entry values.
	typedef std::map<key_t, value_t> map_t; //!< Diff entry key-value map.

	class CandidateJunction
	{
		public:

		struct Settings
		{
			// Fields

			string candidate_junction_fasta_file_name;
			string candidate_junction_faidx_file_name;
			string candidate_junction_sam_file_name;
			string candidate_junction_score_method;
			string jc_genome_diff_file_name;
			string preprocess_junction_best_sam_file_name;
			string preprocess_junction_split_sam_file_name;
			string reference_fasta_file_name;
			string reference_faidx_file_name;
			string reference_sam_file_name;
			string resolved_reference_sam_file_name;
			string resolved_junction_sam_file_name;
			string unmatched_read_file_name;

			bool no_junction_prediction;
			bool unmatched_reads;
			bool add_split_junction_sides;
			bool require_complete_match;

			int32_t alignment_read_limit;
			int32_t candidate_junction_read_limit;
			int32_t max_read_length;
			int32_t maximum_read_mismatches;
			int32_t required_match_length;

			boost::optional<int32_t> preprocess_junction_min_indel_split_length;

			struct ReadStructure
			{
				string base_name;
			};
			vector<ReadStructure> read_structures;

			// Constructors

			Settings();
		};

		struct Summary
		{
			// Fields

			struct AlignmentCorrection
			{
				string read_file;
				struct NewJunction
				{
					int32_t observed_min_overlap_score_distribution;
					int32_t accepted_min_overlap_score_distribution;
					int32_t observed_pos_hash_score_distribution;
					int32_t accepted_pos_hash_score_distribution;
				};
				list<NewJunction> new_junctions;

			} alignment_correction;

			struct PreprocessCoverage
			{
				int32_t junction_accept_score_cutoff_1;
				int32_t junction_accept_score_cutoff_2;
			} preprocess_coverage;

			// Constructors

			Summary();
		};

		/*! Preprocesses alignments
		 */
		static void preprocess_alignments(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info);

		/*! Predicts candidate junctions
		 */
		static void identify_candidate_junctions(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info);


		private:

		CandidateJunction();

		static void _alignments_to_candidate_junction(map_t settings, map_t summary, map_t ref_seq_info, map_t fai, map_t header, map_t a1, map_t a2, map_t redundancy_1, map_t redundancy_2);
		static void _alignments_to_candidate_junctions(map_t settings, map_t summary, map_t ref_seq_info, map_t candidate_junctions, map_t fai, map_t header, map_t al_ref);
		static void _check_read_pair_requirements(map_t a1_start, map_t a1_end, map_t a2_start, map_t a2_end);
		static void _entire_read_matches(map_t a);
		static void _num_matches_from_end(map_t a, map_t refseq_str, map_t dir, map_t overlap);
		static void _split_indel_alignments(Settings settings, Summary summary, bam_header_t* header, ofstream& PSAM, int32_t min_indel_split_len, vector<bam1_t*> al_ref);
		static void _by_ref_seq_coord(map_t a, map_t b, map_t ref_seq_info);
		static void _by_score_unique_coord(map_t a, map_t b);
		static void _tam_write_split_alignment(map_t fh, map_t header, map_t min_indel_split_len, map_t a);

	}; // class CandidateJunction

} // namespace breseq

#endif
