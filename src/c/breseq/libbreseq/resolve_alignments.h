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
#include "genome_diff.h"
#include "calculate_trims.h"
#include "candidate_junctions.h"


using namespace std;

namespace breseq {

  class CandidateJunction
	{
  public:
		int32_t pos_hash_score;
		map<uint32_t, uint32_t> read_begin_hash;
    
		struct TestInfo {
			int32_t max_left;
			int32_t max_left_minus;
			int32_t max_left_plus;
			int32_t max_right;
			int32_t max_right_minus;
			int32_t max_right_plus;
			int32_t max_min_right;
			int32_t max_min_right_minus;
			int32_t max_min_right_plus;
			int32_t max_min_left;
			int32_t max_min_left_minus;
			int32_t max_min_left_plus;
			uint32_t coverage_minus;
			uint32_t coverage_plus;
			uint32_t total_non_overlap_reads;
			uint32_t pos_hash_score;
      bool redundant_1;
      bool redundant_2;
		} test_info;
    
		struct Sorter {      
      bool operator() (const string& lhs, const string& rhs) const
      { 
        // if (lhs.size() != rhs.size()) 
        //  return lhs.size()<rhs.size();
        return  lhs < rhs;
        
        //return lhs.compare(rhs);
      }
		};
	};

  
	class MatchedJunction
	{
  public:
    MatchedJunction() {}
    
    MatchedJunction(
                    const alignment_list& _reference_alignments,
                    const alignment_list&  _junction_alignments,
                    uint32_t _fastq_file_index,
                    int32_t _mapping_quality_difference,
                    uint32_t _degenerate_count
                    )
          :reference_alignments(_reference_alignments)
          ,junction_alignments(_junction_alignments)
          ,fastq_file_index(_fastq_file_index)
          ,mapping_quality_difference(_mapping_quality_difference)
          ,degenerate_count(_degenerate_count)
    { }
    
		alignment_list reference_alignments;
		alignment_list junction_alignments;
		uint32_t fastq_file_index;
		int32_t mapping_quality_difference;
		uint32_t degenerate_count;
	};

	struct VectorSize {
		string junction_id; uint32_t size; uint32_t size2;
    VectorSize(string _junction_id, uint32_t _size, uint32_t _size2) 
      : junction_id(_junction_id), size(_size), size2(_size2) {}
		static bool sort_reverse_by_size(const VectorSize& lhs, const VectorSize& rhs) {
			if (lhs.size != rhs.size) return lhs.size > rhs.size;
      if (lhs.size2 != rhs.size2) return lhs.size2 > rhs.size2;
      return lhs.junction_id < rhs.junction_id;
		}
	};

	bool _alignment_overlaps_junction(const vector<JunctionInfo>& junction_info_list, const alignment_wrapper& in_a);
	diff_entry _junction_to_hybrid_list_item(const string& key, cReferenceSequences& ref_seq_info, CandidateJunction& test_info);
	bool _test_junction(const Settings& settings, Summary& summary, const string& junction_seq_id, map<string, vector<MatchedJunction> >& matched_junction_ref, map<string, map<string, MatchedJunction> >& degenerate_matches_ref, map<string, CandidateJunction>& junction_test_info_ref, cReferenceSequences& ref_seq_info, const SequenceTrimsList& trims_list, tam_file& reference_tam, tam_file& junction_tam, bool& has_non_overlap_alignment);
	Trims _trim_ambiguous_ends(const alignment_wrapper& a, const SequenceTrimsList& trim_list);
  void read_trims(SequenceTrimsList& trims, const cReferenceSequences& ref_seqs, const string &in_trims_file_name ); 

	void _write_reference_matches(const Settings& settings, cReferenceSequences& ref_seq_info, const SequenceTrimsList& trim_list, alignment_list& reference_alignments, tam_file& reference_tam, uint32_t fastq_file_index);

	vector<string> get_sorted_junction_ids(map<string, vector<MatchedJunction> >& unique_map, map<string, map<string, MatchedJunction> >& degenerate_map, const vector<string>& keys);


	void resolve_alignments(
                          Settings& settings,
                          Summary& summary,
                          cReferenceSequences& ref_seq_info,
                          const bool junction_prediction,
                          cReadFiles &read_files,
                          const uint32_t max_read_length,
                          const uint32_t alignment_read_limit
                          );

}

#endif
