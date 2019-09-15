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
  
  class cReferenceSequences;
  class cAnnotatedSequence;
  class Settings;
  class cReferenceSequenceSettings;
  
  // PRIVATE Summaries
  //
  // Files for internal use only. Not documented.
  
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
  
  
  class AlignmentResolutionReferenceSummary : public JSONStorable<AlignmentResolutionReferenceSummary>
  {
  public:
    // These amounts are divided by the number of different equal mappings
    double reads_mapped_to_reference;
    double bases_mapped_to_reference;
    
    AlignmentResolutionReferenceSummary()
    : reads_mapped_to_reference(0)
    , bases_mapped_to_reference(0)
    {}
    
    virtual ~AlignmentResolutionReferenceSummary() {}

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
    
    uint64_t total_reads_mapped_to_references;
    uint64_t total_bases_mapped_to_references;
    
    // This is stored as a vector rather than as a map for speed in adding them up
    vector<AlignmentResolutionReferenceSummary> reference;
    
    AlignmentResolutionSummary()
    : total_unmatched_reads(0)
    , total_unmatched_bases(0)
    , total_reads(0)
    , total_bases(0)
    , max_sam_base_quality_score(0)
    , total_reads_mapped_to_references(0)
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
    uint32_t read_length_min;
    uint32_t read_length_max;
    double read_length_avg;
    // original_reads = homopolymer_filtered_reads + num_filtered_too_many_N_reads + num_filtered_coverage_limit_reads + num_reads
    uint64_t num_original_reads;
    uint64_t num_filtered_too_short_reads;
    uint64_t num_filtered_same_base_reads;
    uint64_t num_filtered_too_many_N_reads;
    uint64_t num_filtered_coverage_limit_reads;
    uint64_t num_reads;
    uint32_t min_quality_score;
    uint32_t max_quality_score;
    uint64_t num_original_bases;
    uint64_t num_bases;
    string quality_format_original;
    string quality_format;
    string converted_fastq_name;
    
    AnalyzeFastqSummary()
    : read_length_min(0)
    , read_length_max(0)
    , read_length_avg(0.0)
    , num_original_reads(0)
    , num_filtered_too_short_reads(0)
    , num_filtered_same_base_reads(0)
    , num_filtered_too_many_N_reads(0)
    , num_filtered_coverage_limit_reads(0)
    , num_reads(0)
    , min_quality_score(0)
    , max_quality_score(0)
    , num_original_bases(0)
    , num_bases(0)
    { }
    
    AnalyzeFastqSummary(
                 uint32_t _read_length_min,
                 uint32_t _read_length_max,
                 double _read_length_avg,
                 uint64_t _num_original_reads,
                 uint64_t _num_filtered_too_short_reads,
                 uint64_t _num_filtered_same_base_reads,
                 uint64_t _num_filtered_too_many_N_reads,
                 uint64_t _num_filtered_coverage_limit_reads,
                 uint64_t _num_reads,
                 uint32_t _min_quality_score,
                 uint32_t _max_quality_score,
                 uint64_t _num_original_bases,
                 uint64_t _num_bases,
                 const string& _quality_format_original,
                 const string& _quality_format,
                 const string& _converted_fastq_name
                 )
    : read_length_min(_read_length_min)
    , read_length_max(_read_length_max)
    , read_length_avg(_read_length_avg)
    , num_original_reads(_num_original_reads)
    , num_filtered_too_short_reads(_num_filtered_too_short_reads)
    , num_filtered_same_base_reads(_num_filtered_same_base_reads)
    , num_filtered_too_many_N_reads(_num_filtered_too_many_N_reads)
    , num_filtered_coverage_limit_reads(_num_filtered_coverage_limit_reads)
    , num_reads(_num_reads)
    , min_quality_score(_min_quality_score)
    , max_quality_score(_max_quality_score)
    , num_original_bases(_num_original_bases)
    , num_bases(_num_bases)
    , quality_format_original(_quality_format_original)
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
    float read_length_avg;
    uint32_t max_qual;
    uint64_t num_reads;
    uint64_t num_original_reads;
    uint64_t num_bases;
    uint64_t num_original_bases;
    map<string, AnalyzeFastqSummary> reads;
    uint64_t total_reference_sequence_length;
    uint32_t read_length_max;
    uint32_t read_length_min;
    
    SequenceConversionSummary()
    : read_length_avg(0.0)
    , max_qual(255)
    , num_reads(0)
    , num_original_reads(0)
    , num_bases(0)
    , num_original_bases(0)
    , total_reference_sequence_length(0)
    , read_length_max(0)
    , read_length_min(0)
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

  // AlignmentResolutionReferenceSummary
  void to_json(json& j, const AlignmentResolutionReferenceSummary& s);
  void from_json(const json& j, AlignmentResolutionReferenceSummary& s);
  
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
  
  
  // PUBLIC Summaries
  //
  // Special summary objects refactored to collate all information about read files / reference files
  // together and be output as a file for users. Should have constructors that take other summary objects
  // or information about the run to let them transfer over the fields.
  
  class PublicReadFileSummary : public JSONStorable<PublicReadFileSummary> {
  public:
    uint32_t read_length_min;
    uint32_t read_length_max;
    double read_length_avg;
    // num_original_reads = homopolymer_filtered_reads + num_filtered_too_many_N_reads + num_reads
    uint64_t num_original_reads;
    uint64_t num_filtered_too_short_reads;
    uint64_t num_filtered_same_base_reads;
    uint64_t num_filtered_too_many_N_reads;
    uint64_t num_reads;
    uint32_t min_quality_score;
    uint32_t max_quality_score;
    uint64_t num_original_bases;
    uint64_t num_bases;
    string quality_format_original;
    string quality_format;
    
    uint64_t num_aligned_reads;
    uint64_t num_aligned_bases;
    double fraction_aligned_reads;
    double fraction_aligned_bases;
    
    PublicReadFileSummary() {}
    PublicReadFileSummary(const ReadFileSummary &rfs, const AnalyzeFastqSummary &afs);
  };
  
  class PublicReadFileSummaries : public map<string, PublicReadFileSummary>, public JSONStorable<PublicReadFileSummaries> {
  public:
    PublicReadFileSummaries() {}
  };
  
  class PublicReadSummary : public JSONStorable<PublicReadSummary> {
  public:
    
    uint64_t total_reads;
    uint64_t total_bases;
    uint64_t total_aligned_reads;
    uint64_t total_aligned_bases;
    double total_fraction_aligned_reads;
    double total_fraction_aligned_bases;
    
    PublicReadFileSummaries read_file;
    
    PublicReadSummary() {}
    PublicReadSummary(const Summary &s);
  };
  
  class PublicReferenceSummary : public JSONStorable<PublicReferenceSummary> {
  public:
    
    uint64_t length;
    uint64_t num_features;
    uint64_t num_genes;
    uint64_t num_repeats;
    
    double num_reads_mapped_to_reference;
    double num_bases_mapped_to_reference;
    
    double coverage_deletion_coverage_propagation_cutoff;
    double coverage_deletion_coverage_seed_cutoff;
    double coverage_nbinom_size_parameter;
    double coverage_nbinom_mean_parameter;
    double coverage_nbinom_prob_parameter;
    double coverage_nbinom_variance;
    double coverage_nbinom_dispersion;
    double coverage_average;
    double coverage_variance;
    double coverage_dispersion;
    
    int32_t coverage_group;
    bool junction_only;
    
    PublicReferenceSummary() {}
    PublicReferenceSummary(
                           const AlignmentResolutionReferenceSummary &arrs,
                           const CoverageSummary &cs,
                           const cAnnotatedSequence &r,
                           const cReferenceSequenceSettings &rss
                           );
  };

  class PublicReferenceSummaries : public map<string, PublicReferenceSummary>, public JSONStorable<PublicReferenceSummaries> {
  public:
    PublicReferenceSummaries() {}
  };
  
  class PublicReferencesSummary : public JSONStorable<PublicReferencesSummary> {
  public:

    uint64_t total_length;
    uint64_t total_features;
    uint64_t total_genes;
    uint64_t total_repeats;
    
    PublicReferenceSummaries reference;
    
    PublicReferencesSummary() {}
    PublicReferencesSummary(const Summary &s, const cReferenceSequences& r, const cReferenceSequenceSettings &rss);
  };
  
  class PublicOptionsSummary : public JSONStorable<PublicOptionsSummary> {
  public:
    
    //! Settings: Workflow
    string custom_run_name;
    int32_t num_processors;
    bool skip_read_filtering;
    bool skip_junction_prediction;
    bool skip_mutation_prediction;
    bool skip_deletion_prediction;
    bool skip_alignment_or_plot_generation;

    //! Settings: Read File
    bool aligned_sam_mode;
    double  read_file_coverage_fold_limit;
    uint32_t read_file_read_length_min;
    double read_file_max_same_base_fraction;
    double read_file_max_N_fraction;
    
    //! Settings: Read Alignment
    string bowtie2_scoring;
    string bowtie2_stage1;
    string bowtie2_stage2;
    string bowtie2_junction;
    uint64_t bowtie2_junction_maximum_alignments_to_consider_per_read;
    uint64_t bowtie2_genome_maximum_alignments_to_consider_per_read;
    uint32_t minimum_mapping_quality;
    uint32_t require_match_length;
    double   require_match_fraction;
    int32_t  maximum_read_mismatches;
    
    //! Settings: Candidate Junction
    int32_t  preprocess_junction_min_indel_split_length;
    int32_t required_both_unique_length_per_side;
    double   required_both_unique_length_per_side_fraction;
    int32_t required_one_unique_length_per_side;
    uint32_t unmatched_end_minimum_read_length;
    double   unmatched_end_length_factor;
    uint32_t maximum_junction_sequence_insertion_length;
    uint32_t maximum_junction_sequence_overlap_length;
    double maximum_junction_sequence_negative_overlap_length_fraction;
    uint32_t maximum_junction_sequence_negative_overlap_length_minimum;
    double maximum_junction_sequence_positive_overlap_length_fraction;
    uint32_t maximum_junction_sequence_positive_overlap_length_minimum;
    uint32_t highly_redundant_junction_ignore_passed_pair_limit;
    uint64_t maximum_junction_sequence_passed_alignment_pairs_to_consider;
    uint32_t minimum_candidate_junction_pos_hash_score;
    uint32_t minimum_candidate_junctions;
    uint32_t maximum_candidate_junctions;
    double maximum_candidate_junction_length_factor;

    //! Settings: Alignment Resolution
    bool add_split_junction_sides;
    uint32_t minimum_alignment_resolution_pos_hash_score;
    double minimum_pr_no_read_start_per_position;
    int32_t junction_minimum_side_match;
    double junction_pos_hash_neg_log10_p_value_cutoff;
    
    //! Settings: Mutation Identification
    string user_evidence_genome_diff_file_name;
    uint32_t base_quality_cutoff;
    uint32_t quality_score_trim;
    double deletion_coverage_propagation_cutoff;
    double deletion_coverage_seed_cutoff;
    bool polymorphism_prediction;
    bool mixed_base_prediction;
    bool targeted_sequencing;
    bool print_mutation_identification_per_position_file;
    double mutation_log10_e_value_cutoff;
    uint32_t consensus_minimum_variant_coverage_each_strand;
    uint32_t consensus_minimum_total_coverage_each_strand;
    uint32_t consensus_minimum_variant_coverage;
    uint32_t consensus_minimum_total_coverage;
    uint32_t consensus_reject_indel_homopolymer_length;
    uint32_t consensus_reject_surrounding_homopolymer_length;
    
    double consensus_frequency_cutoff;
    double polymorphism_log10_e_value_cutoff;
    double polymorphism_bias_p_value_cutoff;
    double polymorphism_frequency_cutoff;
    uint32_t polymorphism_minimum_variant_coverage_each_strand;
    uint32_t polymorphism_minimum_total_coverage_each_strand;
    uint32_t polymorphism_minimum_variant_coverage;
    uint32_t polymorphism_minimum_total_coverage;
    uint32_t polymorphism_reject_indel_homopolymer_length;
    uint32_t polymorphism_reject_surrounding_homopolymer_length;
    bool no_indel_polymorphisms;

    //! Settings: Mutation Prediction
    int32_t size_cutoff_AMP_becomes_INS_DEL_mutation;
    
    //! Settings: Output
    uint32_t max_displayed_reads;
    bool no_javascript;
    string header_genome_diff_file_name;
    uint32_t max_nucleotides_to_show_in_tables;
    uint32_t max_rejected_read_alignment_evidence_to_show;
    uint32_t max_rejected_junction_evidence_to_show;
    bool hide_circular_genome_junctions;

    PublicOptionsSummary() {}
    PublicOptionsSummary(const Settings &t);
  };
  
  class PublicSummary : public JSONStorable<PublicSummary>
  {
  public:
    PublicReadSummary reads;
    PublicReferencesSummary references;
    PublicOptionsSummary options;
    
    PublicSummary() {}
    PublicSummary(const Summary &s, const Settings &t, const cReferenceSequences &r);
  };
  
  // PublicReadFileSummary
  void to_json(json& j, const PublicReadFileSummary& s);
  void from_json(const json& j, PublicReadFileSummary& s);
  
  // PublicReadFileSummaries
  void to_json(json& j, const PublicReadFileSummaries& s);
  void from_json(const json& j, PublicReadFileSummaries& s);

  // PublicReadSummary
  void to_json(json& j, const PublicReadSummary& s);
  void from_json(const json& j, PublicReadSummary& s);

  // PublicReferenceSummary
  void to_json(json& j, const PublicReferenceSummary& s);
  void from_json(const json& j, PublicReferenceSummary& s);

  // PublicReferenceSummaries
  void to_json(json& j, const PublicReferenceSummaries& s);
  void from_json(const json& j, PublicReferenceSummaries& s);
  
  // PublicReferencesSummary
  void to_json(json& j, const PublicReferencesSummary& s);
  void from_json(const json& j, PublicReferencesSummary& s);

  // PublicOptionsSummary
  void to_json(json& j, const PublicOptionsSummary& s);
  void from_json(const json& j, PublicOptionsSummary& s);

  // PublicSummary
  void to_json(json& j, const PublicSummary& s);
  void from_json(const json& j, PublicSummary& s);

} // breseq namespace

#endif
