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
	void CandidateJunction::_alignments_to_candidate_junctions(map_t settings, map_t summary, map_t ref_seq_info, map_t candidate_junctions, map_t fai, map_t header, map_t al_ref) {}
	void CandidateJunction::_check_read_pair_requirements(map_t a1_start, map_t a1_end, map_t a2_start, map_t a2_end) {}
	void CandidateJunction::_entire_read_matches(map_t a) {}
	void CandidateJunction::_num_matches_from_end(map_t a, map_t refseq_str, map_t dir, map_t overlap) {}
	void CandidateJunction::_split_indel_alignments(Settings settings, Summary summary, bam_header_t* header, ofstream& PSAM, int32_t min_indel_split_len, vector<bam1_t*> al_ref) {}
	void CandidateJunction::_by_ref_seq_coord(map_t a, map_t b, map_t ref_seq_info) {}
	void CandidateJunction::_by_score_unique_coord(map_t a, map_t b) {}
	void CandidateJunction::_tam_write_split_alignment(map_t fh, map_t header, map_t min_indel_split_len, map_t a) {}

	// Public

	CandidateJunction::Settings::Settings(){}

	CandidateJunction::Summary::Summary(){}

	/*! Preprocesses alignments
	 */
	void CandidateJunction::preprocess_alignments(Settings settings, Summary summary, const cReferenceSequences& ref_seq_info)
	{
		string read_file = "???";
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
			cerr << "  READ FILE::" << read_file << endl;

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
					int32_t best_score;
//					($best_score, @$al_ref) = Breseq::AlignmentCorrection::_eligible_read_alignments($settings, $header, $reference_fai, $ref_seq_info, @$al_ref);
					tam_write_read_alignments(BSAM, header, 0, al_ref);
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
		cout << "Identifying candidate junctions." << endl;
	}
  
} // namespace breseq
