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


#ifndef _BRESEQ_SETTINGS_H_
#define _BRESEQ_SETTINGS_H_

#include <map>
#include <string>
#include <vector>
#include <list>
#include <stdint.h>

using namespace std;

namespace breseq {

  class cReferenceSequences;

  // We need to be able to group read files for two reasons
  // 1) They may be paired-end, so we want to map them together
  // 2) They may have the same error rates, so we want to treat them together for error analysis

  struct cReadFile {
  public:
    string m_fastq_file_name;
    string m_base_name;
    uint32_t m_paired_end_group;    // indicated what file contains paired reads
    uint32_t m_error_group;         // indicates what other read files have the same error rates
    uint32_t m_id;                  // index used to refer to this fastq file in BAM
  };
  
  typedef vector<vector<cReadFile> > cReadFileGroup;
  
  
  class cReadFiles : public vector<cReadFile> {
    
  protected:
    
  public:
    cReadFiles() {};
    cReadFiles(const vector<string>& read_file_names);
    ~cReadFiles() {};
    
    void Init(const vector<string>& read_file_names);
    
  };
  
  struct Settings
	{
    
    // Set up defaults here
    Settings(const string& run_path="");
    
		// Fields
		map<string, bool> installed;
		string bin_path;
    
    // Path to files....
    string sequence_conversion_path;
    string reference_trim_file_name;
    
    string reference_alignment_path;
		string reference_sam_file_name;
    
    string candidate_junction_path;
		string candidate_junction_fasta_file_name;
		string candidate_junction_faidx_file_name;
		string jc_genome_diff_file_name;
		string preprocess_junction_best_sam_file_name;
		string preprocess_junction_split_sam_file_name;

    string candidate_junction_alignment_path;
    string candidate_junction_sam_file_name; 
    
    string alignment_correction_path;
		string resolved_reference_sam_file_name;
		string resolved_junction_sam_file_name;
    
    string data_path;
		string reference_fasta_file_name;
		string reference_faidx_file_name;
    string reference_features_file_name;
		string unmatched_read_file_name;
    
    string output_path;
    string local_evidence_path;
    string evidence_path;
    string evidence_genome_diff_file_name;
    string final_genome_diff_file_name;
    
    // Options...
    
		bool no_junction_prediction;
    string candidate_junction_score_method;
		bool unmatched_reads;
		bool add_split_junction_sides;
		bool require_complete_match;
    
		uint32_t alignment_read_limit;
		uint32_t candidate_junction_read_limit;
		uint32_t minimum_candidate_junction_pos_hash_score;
		uint32_t minimum_candidate_junction_min_overlap_score;
		int32_t minimum_candidate_junctions;
		int32_t maximum_candidate_junctions;
		double maximum_candidate_junction_length_factor;
		int32_t maximum_read_length;
		int32_t maximum_inserted_junction_sequence_length;
		int32_t max_read_mismatches;
		int32_t required_both_unique_length_per_side;
		int32_t required_one_unique_length_per_side;
		int32_t required_extra_pair_total_length;
		int32_t required_match_length;
    
		int32_t preprocess_junction_min_indel_split_length;
    
		cReadFiles read_structures;

		struct Coverage {
			int32_t junction_coverage_cutoff;
		};
		map<string,Coverage> unique_coverage;
    
    //!@GRC Setting options needed for HTML outputs
    string print_run_name; //!< need to set default to "unnamed"
    bool hide_circular_genome_junctions;
  	bool polymorphism_prediction;
		bool lenski_format;
		bool no_evidence;
		bool shade_frequencies;
		bool no_header;
    
		// Utility function to substitute specific details into a generic file name
		static string file_name(const string& file_name_key, const string& substitute = "", const string& with = "")
		{
      string s(file_name_key);
        
      if (substitute.size() > 0)
      {
        size_t pos = s.find(substitute);
        if (pos != string::npos)
        {
          s.replace(pos, 1, with);
        }
      }
      
			return s;
		}
    
    // assumes things are in our path for now
    string ctool(string tool_name) const
		{
      //			my ($self, $tool_name, $allow_fail) = @_;
      //
      //			if (!$self->{installed}->{$tool_name})
      //			{
      //				if ($allow_fail)
      //				{
      //					$self->warn("Executable \"$tool_name\" not found in breseq bin path\"$self->{bin_path}\".");
      //					return undef; # couldn't find it, but it's not an error.
      //				}
      //				else
      //				{
      //					$self->throw("Executable \"$tool_name\" not found in breseq bin path\"$self->{bin_path}\".");
      //				}
      //			}
      
			return  tool_name;
		}
	};

  struct Summary
	{
		// Fields
    
		struct AlignmentCorrection
		{
			map<string, map<string, int32_t> > read_file;
			struct NewJunction
			{
				map<int32_t,int32_t> observed_min_overlap_score_distribution;
				map<int32_t,int32_t> accepted_min_overlap_score_distribution;
				map<int32_t,int32_t> observed_pos_hash_score_distribution;
				map<int32_t,int32_t> accepted_pos_hash_score_distribution;
			} new_junctions;
      
		} alignment_correction;
    
		struct Coverage {
			int32_t junction_accept_score_cutoff;
		};
		map<string,Coverage> preprocess_coverage;
		map<string,Coverage> unique_coverage;

		struct CandidateJunctionSummaryData
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
      
			map<int32_t, int32_t> pos_hash_score_distribution;
			map<int32_t, int32_t> min_overlap_score_distribution;
      
			map<string, map<string, int32_t> > read_file;
		} candidate_junction;
    
		struct SequenceConversion
		{
			uint32_t total_reference_sequence_length;
			uint32_t max_read_length;
			cReferenceSequences* reference_sequences;
		} sequence_conversion;
	};

  

} // breseq namespace

#endif
