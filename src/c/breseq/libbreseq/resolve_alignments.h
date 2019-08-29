/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2017 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_RESOLVE_ALIGNMENTS_H_
#define _BRESEQ_RESOLVE_ALIGNMENTS_H_

#include "common.h"

#include "reference_sequence.h"
#include "genome_diff.h"
#include "calculate_trims.h"
#include "candidate_junctions.h"
#include "pileup_base.h"


using namespace std;

namespace breseq {

  class ResolveJunctionInfo : public JunctionInfo
  {
  public:
    // Extended properties for resolve_alignments.cpp
		string key;
		int32_t overlap;
		uint32_t unique_side;
		uint32_t is_side;
    
    ResolveJunctionInfo() : JunctionInfo() {};
    
    ResolveJunctionInfo(const string& junction_name) 
    : JunctionInfo(junction_name)
    , key(junction_name)
    , overlap(0)
    , is_side(0)
    {
    }
  };
  
  class JunctionTestInfo {
  public:
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
    uint32_t max_pos_hash_score;
    uint32_t unique_matches_size;
    uint32_t repeat_matches_size;
    bool has_reads_with_both_different_start_and_end;
    bool redundant_1;
    bool redundant_2;
    string junction_id;
    double neg_log10_pos_hash_p_value;
    uint32_t side_1_continuation;
    uint32_t side_2_continuation;
    vector<string> reject_reasons;
    
    bool operator <(const JunctionTestInfo& _in) const
    {
      // sort by pos_hash_score, unique_matches_size, repeat_matches_size
      
      if (this->pos_hash_score != _in.pos_hash_score) {
        return (this->pos_hash_score < _in.pos_hash_score);
      }
      
      if (this->unique_matches_size != _in.unique_matches_size) {
        return (this->unique_matches_size < _in.unique_matches_size);
      }
      
      return (this->repeat_matches_size < _in.repeat_matches_size);
    }
  };
  
    
  typedef map<string, JunctionTestInfo> JunctionTestInfoMap;
  
	class JunctionMatch
	{
  public:
    alignment_list reference_alignments;
		alignment_list junction_alignments;
		uint32_t fastq_file_index;
		int32_t mapping_quality_difference;
		uint32_t degenerate_count;
    
    JunctionMatch() {}
    
    JunctionMatch(
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
  
  typedef counted_ptr<JunctionMatch> JunctionMatchPtr;
  typedef map<string, vector<JunctionMatchPtr> > UniqueJunctionMatchMap;      // Map of junction_id to list of MatchedJunctions
  typedef map<string, map<string, JunctionMatchPtr> > RepeatJunctionMatchMap; // Map of junction_id to read_name to MatchedJunction
  
  
  class PosHashProbabilityTable {
  public:
        
    struct Parameters {
      double negative_binomial_size;
      double negative_binomial_prob;
      double chance_per_pos_strand_no_read_start;
      double average_coverage;
    };
    
    map<string, Parameters> param;
    uint32_t average_read_length;
    map<string, map<uint32_t, map<uint32_t, double> > > probability_table;
    
    PosHashProbabilityTable(Summary& summary);
    
    double probability(string& seq_id, uint32_t pos_hash_score, uint32_t overlap);
  };

  
  class pos_hash_p_value_table : public vector<vector<double> > {
  public:
    pos_hash_p_value_table() {};
    
    pos_hash_p_value_table(const string& in_file_name)
    {      
      // file may not have been created if fitting failed
      if (!file_exists(in_file_name.c_str())) return;
      
      ifstream in_file(in_file_name.c_str());
      ASSERT(!in_file.fail(), "Could not open file: " + in_file_name);
      
      string line; 
      while (!in_file.eof()) {
        getline(in_file, line);
        if (line == "") continue;
        vector<string> line_list = split(line, "\t");
        
        vector<double> converted_line;
        for (vector<string>::iterator it=line_list.begin(); it!=line_list.end(); it++) {
          converted_line.push_back(from_string<double>(*it));
        }
        
        (*this).push_back(converted_line);
      }
    }
  };
  
  void calculate_continuation(
                              ResolveJunctionInfo& rji, 
                              cReferenceSequences& ref_seq_info, 
                              cReferenceSequences& junction_ref_seq_info, 
                              uint32_t& side_1_continuation,
                              uint32_t& side_2_continuation
                              );
  
  void resolve_alignments(
                          Settings& settings,
                          Summary& summary,
                          cReferenceSequences& ref_seq_info,
                          const bool junction_prediction,
                          cReadFiles &read_files
                          );
  
  void load_junction_alignments(
                                const Settings& settings, 
                                Summary& summary, 
                                cReadFiles& read_files, 
                                cReferenceSequences& ref_seq_info,
                                cReferenceSequences& junction_ref_seq_info,
                                SequenceTrimsList& trims_list,
                                map<string,uint32_t>& all_junction_ids,
                                bool junction_prediction,
                                const vector<ResolveJunctionInfo>& junction_info_list,
                                UniqueJunctionMatchMap& unique_junction_match_map,
                                RepeatJunctionMatchMap& repeat_junction_match_map,
                                tam_file& resolved_reference_tam
                                );
  
  void load_sam_only_alignments(
                           const Settings& settings, 
                           Summary& summary, 
                           cReadFiles& read_files, 
                           cReferenceSequences& ref_seq_info,
                           SequenceTrimsList& trims_list,
                           tam_file& resolved_reference_tam
                                );
  
  bool alignment_overlaps_junction(const vector<ResolveJunctionInfo>& junction_info_list, const alignment_wrapper& in_a);

  
  void score_junction(
                      const Settings& settings, 
                      Summary& summary, 
                      const string& junction_id, 
                      UniqueJunctionMatchMap& matched_junction_ref, 
                      RepeatJunctionMatchMap& degenerate_matches_ref, 
                      tam_file& junction_tam,
                      JunctionTestInfo& junction_test_info, 
                      vector<ResolveJunctionInfo>& junction_info_list,
                      cReferenceSequences& ref_seq_info, 
                      cReferenceSequences& junction_ref_seq_info
                      );
  
  void resolve_junction(
                        const Settings& settings,
                        Summary& summary,
                        cReferenceSequences& ref_seq_info,
                        cReferenceSequences& junction_ref_seq_info,
                        const SequenceTrimsList& trims_list,
                        const string& junction_id,
                        UniqueJunctionMatchMap& unique_junction_match_map,
                        RepeatJunctionMatchMap& repeat_junction_match_map,
                        tam_file& reference_tam,
                        tam_file& junction_tam,
                        bool failed,
                        bool has_non_overlap_alignment
                        );

	void _write_reference_matches(
                                const Settings& settings, 
                                Summary& summary,
                                cReferenceSequences& ref_seq_info, 
                                const SequenceTrimsList& trim_list, 
                                alignment_list& reference_alignments, 
                                tam_file& reference_tam, 
                                uint32_t fastq_file_index
                                );
                                

	vector<string> get_sorted_junction_ids(
                                         UniqueJunctionMatchMap& unique_map, 
                                         RepeatJunctionMatchMap& degenerate_map, 
                                         const vector<string>& keys
                                         );

  
  cDiffEntry junction_to_diff_entry(
                                    const string& key, 
                                    cReferenceSequences& ref_seq_info, 
                                    JunctionTestInfo& test_info
                                    );

  void  assign_one_junction_read_counts(
                                        const Settings& settings,
                                        Summary& summary,
                                        cDiffEntry& j,
                                        int32_t require_overlap = 0
                                        );
  
  void  assign_junction_read_counts(
                                    const Settings& settings,
                                    Summary& summary,
                                    cGenomeDiff& gd
                                    );
  
  void  assign_junction_read_coverage(
                                      const Settings& settings,
                                      Summary& summary,
                                      cGenomeDiff& gd
                                      );

  
  // Pileup class for fetching reads that align across from start to end
  class junction_read_counter : pileup_base {
  public:
    junction_read_counter(const string& bam, const string& fasta, bool verbose)
      : pileup_base(bam, fasta), _verbose(verbose) {};
    
    uint32_t count(
                   const string& seq_id, 
                   const int32_t start, 
                   const int32_t end, 
                   const map<string,bool> ignore_read_names, 
                   map<string,bool>& counted_read_names
                   );
    
    virtual void fetch_callback ( const alignment_wrapper& a );
    
  protected:
    uint32_t _count;
    int32_t _start;
    int32_t _end;
    map<string,bool> _ignore_read_names, _counted_read_names;
    bool _verbose;
  };

}

#endif
