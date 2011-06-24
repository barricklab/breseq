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

#ifndef _BRESEQ_RESOLVE_ALIGNMENTS_H_
#define _BRESEQ_RESOLVE_ALIGNMENTS_H_

#include "common.h"

#include "annotated_sequence.h"

using namespace std;

namespace breseq {

	struct MatchedJunction
	{
		vector<alignment> reference_alignments;
		vector<alignment> junction_alignments;
		uint32_t fastq_file_index;
		int32_t mapping_quality_difference;
		uint32_t degenerate_count;
	};

	struct VectorSize {
		string junction_id; uint32_t size;
		static bool sort_by_size(const VectorSize& lhs, const VectorSize& rhs) {
			return lhs.size < rhs.size;
		}
	};

	uint32_t _eligible_read_alignments(const Settings& settings, const cReferenceSequences& ref_seq_info, vector<alignment>& alignments);
	bool _test_read_alignment_requirements(const Settings& settings, const cReferenceSequences& ref_seq_info, const alignment& a);
	bool _alignment_overlaps_junction(const vector<JunctionInfo>& junction_info_list, alignment in_a);
	bool _test_junction(const Settings& settings, /*const map<string, uint32_t>& summary_info,*/ const string& junction_seq_id, map<string, vector<MatchedJunction> >& matched_junction_ref, map<string, map<string, MatchedJunction> >& degenerate_matches_ref, map<string, CandidateJunction>& junction_test_info_ref, cReferenceSequences& ref_seq_info, ifstream& RREF, ifstream& RCJ, tam_file& reference_tam, tam_file& junction_tam, bool& has_non_overlap_alignment);
	Trim _trim_ambiguous_ends(const alignment& a, const tam_file& tam, cReferenceSequences& ref_seq_info);
	void _write_reference_matches(const Settings& settings, cReferenceSequences& ref_seq_info, alignment_list& reference_alignments, tam_file& reference_tam, uint32_t fastq_file_index);
	vector<string> get_sorted_junction_ids(map<string, vector<MatchedJunction> >& map, const vector<string>& keys);
	vector<string> get_sorted_junction_ids(map<string, map<string, MatchedJunction> >& map, const vector<string>& keys);


	void resolve_alignments(
						  const bool junction_prediction,
						  const string &reference_fasta,
						  const string &junction_fasta,
						  const string &reference_sam_path,
						  const string &junction_sam_path,
						  const string &resolved_path,
						  const string &data_path,
						  const string &features_file,
						  const cReadFiles &read_files,
						  const uint32_t max_read_length,
						  const uint32_t alignment_read_limit
					  );
	
}

#endif
