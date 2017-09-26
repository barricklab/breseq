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


#ifndef _BRESEQ_SUMMARY_H_
#define _BRESEQ_SUMMARY_H_

#include "common.h"
#include "storable.h"

namespace breseq{
  
  class ReadFileSummary : public JSONStorable<ReadFileSummary> {
  public:
    uint64_t num_unmatched_reads;
    uint64_t num_unmatched_bases;
    uint64_t num_total_reads;
    uint64_t num_total_bases;
    
    
    ReadFileSummary()
    : num_unmatched_reads(0)
    , num_unmatched_bases(0)
    , num_total_reads(0)
    , num_total_bases(0)
    {}
  };

  class PosHashScoreDistribution : public map<int32_t, int32_t> {
  public:
    PosHashScoreDistribution() {}
    
    inline void add_score(int32_t score)
    {
      if (this->count(score) == 0)
        (*this)[score] = 1; // Initialize value
      else
        (*this)[score]++;
    }
  };
  
  class AlignmentResolutionSummary : public JSONStorable<AlignmentResolutionSummary>
  {
  public:
    map<string,ReadFileSummary> read_file;
    uint64_t total_unmatched_reads;
    uint64_t total_unmatched_bases;
    uint64_t total_reads;
    uint64_t total_bases;
    int32_t max_sam_base_quality_score; // @JEB only filled in for aligned_sam_mode
    
    PosHashScoreDistribution observed_pos_hash_score_distribution;
    PosHashScoreDistribution accepted_pos_hash_score_distribution;
    
    double total_reads_mapped_to_references;
    uint64_t total_bases_mapped_to_references;
    vector<double> reads_mapped_to_references;
    
    AlignmentResolutionSummary()
    : total_unmatched_reads(0)
    , total_unmatched_bases(0)
    , total_reads(0)
    , total_bases(0)
    , total_reads_mapped_to_references(0.0)
    , total_bases_mapped_to_references(0)
    {}
    
  };
  
  
  class CoverageSummary : public JSONStorable<CoverageSummary>
  {
  public:
    double deletion_coverage_propagation_cutoff;
    double deletion_coverage_seed_cutoff;
    double nbinom_size_parameter;
    double nbinom_mean_parameter;
    double nbinom_prob_parameter;
    double nbinom_variance;
    double nbinom_dispersion;
    double average;
    double variance;
    double dispersion;
    
    CoverageSummary()
    : deletion_coverage_propagation_cutoff(0.0)
    , deletion_coverage_seed_cutoff(0.0)
    , nbinom_size_parameter(0.0)
    , nbinom_mean_parameter(0.0)
    , nbinom_prob_parameter(0.0)
    , nbinom_variance(0.0)
    , nbinom_dispersion(0.0)
    , average(0.0)
    , variance(0.0)
    , dispersion(0.0)
    {}
  };
  
  class AnalyzeFastqSummary : public JSONStorable<AnalyzeFastqSummary>
  {
  public:
    uint32_t min_read_length;
    uint32_t max_read_length;
    double avg_read_length;
    uint64_t original_reads;
    uint64_t too_short_filtered_reads;
    uint64_t same_base_filtered_reads;
    uint64_t N_filtered_reads;
    uint64_t num_reads;         // original_reads = homopolymer_filtered_reads + N_filtered_reads + num_reads
    uint32_t min_quality_score;
    uint32_t max_quality_score;
    uint64_t original_num_bases;
    uint64_t num_bases;
    string original_qual_format;
    string quality_format;
    string converted_fastq_name;
    
    AnalyzeFastqSummary()
    : min_read_length(0)
    , max_read_length(0)
    , avg_read_length(0.0)
    , original_reads(0)
    , too_short_filtered_reads(0)
    , same_base_filtered_reads(0)
    , N_filtered_reads(0)
    , num_reads(0)
    , min_quality_score(0)
    , max_quality_score(0)
    , original_num_bases(0)
    , num_bases(0)
    { }
    
    AnalyzeFastqSummary(
                 uint32_t _min_read_length,
                 uint32_t _max_read_length,
                 double _avg_read_length,
                 uint64_t _original_reads,
                 uint64_t _too_short_filtered_reads,
                 uint64_t _same_base_filtered_reads,
                 uint64_t _N_filtered_reads,
                 uint64_t _num_reads,
                 uint32_t _min_quality_score,
                 uint32_t _max_quality_score,
                 uint64_t _original_num_bases,
                 uint64_t _num_bases,
                 const string& _original_qual_format,
                 const string& _quality_format,
                 const string& _converted_fastq_name
                 )
    : min_read_length(_min_read_length)
    , max_read_length(_max_read_length)
    , avg_read_length(_avg_read_length)
    , original_reads(_original_reads)
    , too_short_filtered_reads(_too_short_filtered_reads)
    , same_base_filtered_reads(_same_base_filtered_reads)
    , N_filtered_reads(_N_filtered_reads)
    , num_reads(_num_reads)
    , min_quality_score(_min_quality_score)
    , max_quality_score(_max_quality_score)
    , original_num_bases(_original_num_bases)
    , num_bases(_num_bases)
    , original_qual_format(_original_qual_format)
    , quality_format(_quality_format)
    , converted_fastq_name(_converted_fastq_name)
    { }
    
  };
  
  class PreprocessAlignmentsSummary : public JSONStorable<PreprocessAlignmentsSummary>
  {
  public:
    int64_t aligned_reads;
    int64_t alignments;
    int64_t alignments_split_on_indels;
    int64_t reads_with_alignments_split_on_indels;
    int64_t split_alignments;
    int64_t reads_with_split_alignments;
    
    PreprocessAlignmentsSummary()
    : aligned_reads(0)
    , alignments(0)
    , alignments_split_on_indels(0)
    , reads_with_alignments_split_on_indels(0)
    , split_alignments(0)
    , reads_with_split_alignments(0)
    {}
  };
  
  class CandidateJunctionSummary : public JSONStorable<CandidateJunctionSummary>
  {
  public:

    int64_t total_number;
    int64_t total_length;

    int64_t accepted_number;
    int64_t accepted_length;
    int32_t accepted_pos_hash_score_cutoff;
    
    uint64_t passed_alignment_pairs_considered;
    PosHashScoreDistribution pos_hash_score_distribution;
    
    CandidateJunctionSummary() : passed_alignment_pairs_considered(0) {}
  };
  
  class SequenceConversionSummary : public JSONStorable<SequenceConversionSummary>
  {
  public:
    float avg_read_length;
    uint32_t max_qual;
    uint64_t num_reads;
    uint64_t original_num_reads;
    uint64_t num_bases;
    uint64_t original_num_bases;
    map<string, AnalyzeFastqSummary> reads;
    uint64_t total_reference_sequence_length;
    uint32_t max_read_length;
    uint32_t min_read_length;
    
    SequenceConversionSummary()
    : avg_read_length(0.0)
    , max_qual(255)
    , num_reads(0)
    , original_num_reads(0)
    , num_bases(0)
    , original_num_bases(0)
    , total_reference_sequence_length(0)
    , max_read_length(0)
    , min_read_length(0)
    { }
  };
  
  class ErrorCountSummary
  {
  public:
    double no_pos_hash_per_position_pr;
    
    ErrorCountSummary() {}
  };
  
  class CoverageSummaries : public map<string, CoverageSummary>, public JSONStorable<CoverageSummaries>
  {
  public:
    CoverageSummaries() {}
  };
  
  class ErrorCountSummaries : public map<string, ErrorCountSummary>, public JSONStorable<ErrorCountSummaries>
  {
  public:
    ErrorCountSummaries() {}
  };
  
	class Summary : public JSONStorable<Summary>
	{
	public:
    AlignmentResolutionSummary alignment_resolution;
		CoverageSummaries preprocess_coverage;
		CoverageSummaries unique_coverage;
    PreprocessAlignmentsSummary preprocess_alignments;
    CandidateJunctionSummary candidate_junction;
    SequenceConversionSummary sequence_conversion;
    ErrorCountSummaries preprocess_error_count;
    
    Summary() {}
	};
  
  // ReadFileSummary
  void to_json(json& j, const ReadFileSummary& s);
  void from_json(const json& j, ReadFileSummary& s);
  
  // PosHashScoreDistribution
  void to_json(json& j, const PosHashScoreDistribution& s);
  void from_json(const json& j, PosHashScoreDistribution& s);
  
  // AlignmentResolutionSummary
  void to_json(json& j, const AlignmentResolutionSummary& s);
  void from_json(const json& j, AlignmentResolutionSummary& s);
  
  // CoverageSummary
  void to_json(json& j, const CoverageSummary& s);
  void from_json(const json& j, CoverageSummary& s);
  
  // AnalyzeFastqSummary
  void to_json(json& j, const AnalyzeFastqSummary& s);
  void from_json(const json& j, AnalyzeFastqSummary& s);
  
  // PreprocessAlignmentsSummary
  void to_json(json& j, const PreprocessAlignmentsSummary& s);
  void from_json(const json& j, PreprocessAlignmentsSummary& s);
  
  // CandidateJunctionSummary
  void to_json(json& j, const CandidateJunctionSummary& s);
  void from_json(const json& j, CandidateJunctionSummary& s);

  // SequenceConversionSummmary
  void to_json(json& j, const SequenceConversionSummary& s);
  void from_json(const json& j, SequenceConversionSummary& s);
  
  // ErrorCountSummary
  void to_json(json& j, const ErrorCountSummary& s);
  void from_json(const json& j, ErrorCountSummary& s);
  
  // CoverageSummaries
  void to_json(json& j, const CoverageSummaries& s);
  void from_json(const json& j, CoverageSummaries& s);
  
  // ErrorCountSummaries
  void to_json(json& j, const ErrorCountSummaries& s);
  void from_json(const json& j, ErrorCountSummaries& s);
  
  // Summary
  void to_json(json& j, const Summary& s);
  void from_json(const json& j, Summary& s);
  
} // breseq namespace

#endif
