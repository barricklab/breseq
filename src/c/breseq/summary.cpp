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

#include "libbreseq/summary.h"
#include "libbreseq/settings.h"
#include "libbreseq/reference_sequence.h"

using namespace std;
using nlohmann::json;

namespace breseq
{
 
// ReadFileSummary
void to_json(json& j, const ReadFileSummary& s)
{
  j = json{
    {"num_unmatched_reads", s.num_unmatched_reads},
    {"num_unmatched_bases", s.num_unmatched_bases},
    {"num_total_reads", s.num_total_reads},
    {"num_total_bases", s.num_total_bases},
  };
}

void from_json(const json& j, ReadFileSummary& s)
{
  s.num_unmatched_reads = j.at("num_unmatched_reads").get<uint64_t>();
  s.num_unmatched_bases = j.at("num_unmatched_bases").get<uint64_t>();
  s.num_total_reads = j.at("num_total_reads").get<uint64_t>();
  s.num_total_bases = j.at("num_total_bases").get<uint64_t>();
}
  
// PosHashScoreDistribution
void to_json(json& j, const PosHashScoreDistribution& s)
{
  for (PosHashScoreDistribution::const_iterator it = s.begin(); it != s.end(); it++) {
    j[to_string<int32_t>(it->first)] = it->second;
  }
}

void from_json(const json& j, PosHashScoreDistribution& s)
{
  for (json::const_iterator it = j.begin(); it != j.end(); ++it) {
    s[from_string<int32_t>(it.key())] = it.value();
  }
}
 
// AlignmentResolutionReferenceSummary
void to_json(json& j, const AlignmentResolutionReferenceSummary& s)
{
  j = json{
    {"reads_mapped_to_reference", s.reads_mapped_to_reference},
    {"bases_mapped_to_reference", s.bases_mapped_to_reference},
  };
}

void from_json(const json& j, AlignmentResolutionReferenceSummary& s)
{
  s.reads_mapped_to_reference = j.at("reads_mapped_to_reference").get<double>();
  s.bases_mapped_to_reference = j.at("bases_mapped_to_reference").get<double>();
}
  
// AlignmentResolutionSummary
void to_json(json& j, const AlignmentResolutionSummary& s)
{
  j = json{
    {"read_file", s.read_file},
    {"total_unmatched_reads", s.total_unmatched_reads},
    {"total_unmatched_bases", s.total_unmatched_bases},
    {"total_reads", s.total_reads},
    {"total_bases", s.total_bases},
    {"max_sam_base_quality_score", s.max_sam_base_quality_score},
    {"observed_pos_hash_score_distribution", s.observed_pos_hash_score_distribution},
    {"accepted_pos_hash_score_distribution", s.accepted_pos_hash_score_distribution},
    {"total_reads_mapped_to_references", s.total_reads_mapped_to_references},
    {"total_bases_mapped_to_references", s.total_bases_mapped_to_references},
    {"reference", s.reference},
  };
}

void from_json(const json& j, AlignmentResolutionSummary& s)
{
  s.read_file = j.at("read_file").get<map<string,ReadFileSummary> >();
  s.total_unmatched_reads = j.at("total_unmatched_reads").get<uint64_t>();
  s.total_unmatched_bases = j.at("total_unmatched_bases").get<uint64_t>();
  s.total_reads = j.at("total_reads").get<uint64_t>();
  s.total_bases = j.at("total_bases").get<uint64_t>();
  s.max_sam_base_quality_score = j.at("max_sam_base_quality_score").get<int32_t>();
  s.observed_pos_hash_score_distribution = j.at("observed_pos_hash_score_distribution").get<PosHashScoreDistribution>();
  s.accepted_pos_hash_score_distribution = j.at("accepted_pos_hash_score_distribution").get<PosHashScoreDistribution>();
  s.total_reads_mapped_to_references = j.at("total_reads_mapped_to_references").get<uint64_t>();
  s.total_bases_mapped_to_references = j.at("total_bases_mapped_to_references").get<uint64_t>();
  s.reference = j.at("reference").get<vector<AlignmentResolutionReferenceSummary> >();
}
  
// CoverageSummary
void to_json(json& j, const CoverageSummary& s)
{
  j = json{
    {"deletion_coverage_propagation_cutoff", s.deletion_coverage_propagation_cutoff},
    {"deletion_coverage_seed_cutoff", s.deletion_coverage_seed_cutoff},
    {"nbinom_size_parameter", s.nbinom_size_parameter},
    {"nbinom_mean_parameter", s.nbinom_mean_parameter},
    {"nbinom_prob_parameter", s.nbinom_prob_parameter},
    {"nbinom_variance", s.nbinom_variance},
    {"nbinom_dispersion", s.nbinom_dispersion},
    {"average", s.average},
    {"variance", s.variance},
    {"dispersion", s.dispersion},
  };
}

void from_json(const json& j, CoverageSummary& s)
{
  s.deletion_coverage_propagation_cutoff = j.at("deletion_coverage_propagation_cutoff").get<double>();
  s.deletion_coverage_seed_cutoff = j.at("deletion_coverage_seed_cutoff").get<double>();
  s.nbinom_size_parameter = j.at("nbinom_size_parameter").get<double>();
  s.nbinom_mean_parameter = j.at("nbinom_mean_parameter").get<double>();
  s.nbinom_prob_parameter = j.at("nbinom_prob_parameter").get<double>();
  s.nbinom_variance = j.at("nbinom_variance").get<double>();
  s.nbinom_dispersion = j.at("nbinom_dispersion").get<double>();
  s.average = j.at("average").get<double>();
  s.variance = j.at("variance").get<double>();
  s.dispersion = j.at("dispersion").get<double>();
}
  
// AnalyzeFastqSummary
void to_json(json& j, const AnalyzeFastqSummary& s)
{
  j = json{
    {"read_length_min", s.read_length_min},
    {"read_length_max", s.read_length_max},
    {"read_length_avg", s.read_length_avg},
    {"num_original_reads", s.num_original_reads},
    {"num_filtered_too_short_reads", s.num_filtered_too_short_reads},
    {"num_filtered_same_base_reads", s.num_filtered_same_base_reads},
    {"num_filtered_too_many_N_reads", s.num_filtered_too_many_N_reads},
    {"num_filtered_coverage_limit_reads", s.num_filtered_coverage_limit_reads},
    {"num_reads", s.num_reads},
    {"min_quality_score", s.min_quality_score},
    {"max_quality_score", s.max_quality_score},
    {"num_original_bases", s.num_original_bases},
    {"num_bases", s.num_bases},
    {"quality_format_original", s.quality_format_original},
    {"quality_format", s.quality_format},
    {"converted_fastq_name", s.converted_fastq_name},
  };
}

void from_json(const json& j, AnalyzeFastqSummary& s)
{
  s.read_length_min = j.at("read_length_min").get<uint32_t>();
  s.read_length_max = j.at("read_length_max").get<uint32_t>();
  s.read_length_avg = j.at("read_length_avg").get<double>();
  s.num_original_reads = j.at("num_original_reads").get<uint64_t>();
  s.num_filtered_too_short_reads = j.at("num_filtered_too_short_reads").get<uint64_t>();
  s.num_filtered_same_base_reads = j.at("num_filtered_same_base_reads").get<uint64_t>();
  s.num_filtered_too_many_N_reads = j.at("num_filtered_too_many_N_reads").get<uint64_t>();
  s.num_filtered_coverage_limit_reads = j.at("num_filtered_coverage_limit_reads").get<uint64_t>();
  s.num_reads = j.at("num_reads").get<uint64_t>();
  s.min_quality_score = j.at("min_quality_score").get<uint32_t>();
  s.max_quality_score = j.at("max_quality_score").get<uint32_t>();
  s.num_original_bases = j.at("num_original_bases").get<uint64_t>();
  s.num_bases = j.at("num_bases").get<uint64_t>();
  s.quality_format_original = j.at("quality_format_original").get<string>();
  s.quality_format = j.at("quality_format").get<string>();
  s.converted_fastq_name = j.at("converted_fastq_name").get<string>();
}
  
// PreprocessAlignmentsSummary
void to_json(json& j, const PreprocessAlignmentsSummary& s)
{
  j = json{
    {"aligned_reads", s.aligned_reads},
    {"alignments", s.alignments},
    {"alignments_split_on_indels", s.alignments_split_on_indels},
    {"reads_with_alignments_split_on_indels", s.reads_with_alignments_split_on_indels},
    {"split_alignments", s.split_alignments},
    {"reads_with_split_alignments", s.reads_with_split_alignments},
  };
}

void from_json(const json& j, PreprocessAlignmentsSummary& s)
{
  s.aligned_reads = j.at("aligned_reads").get<uint64_t>();
  s.alignments = j.at("alignments").get<uint64_t>();
  s.alignments_split_on_indels = j.at("alignments_split_on_indels").get<uint64_t>();
  s.reads_with_alignments_split_on_indels = j.at("reads_with_alignments_split_on_indels").get<uint64_t>();
  s.split_alignments = j.at("split_alignments").get<uint64_t>();
  s.reads_with_split_alignments = j.at("reads_with_split_alignments").get<uint64_t>();
}

// CandidateJunctionSummary
void to_json(json& j, const CandidateJunctionSummary& s)
{
  j = json{
    {"total_number", s.total_number},
    {"total_length", s.total_length},
    {"accepted_number", s.accepted_number},
    {"accepted_length", s.accepted_length},
    {"accepted_pos_hash_score_cutoff", s.accepted_pos_hash_score_cutoff},
    {"passed_alignment_pairs_considered", s.passed_alignment_pairs_considered},
    {"pos_hash_score_distribution", s.pos_hash_score_distribution},
  };
}

void from_json(const json& j, CandidateJunctionSummary& s)
{
  s.total_number = j.at("total_number").get<uint64_t>();
  s.total_length = j.at("total_length").get<uint64_t>();
  s.accepted_number = j.at("accepted_number").get<uint64_t>();
  s.accepted_length = j.at("accepted_length").get<uint64_t>();
  s.accepted_pos_hash_score_cutoff = j.at("accepted_pos_hash_score_cutoff").get<uint64_t>();
  s.passed_alignment_pairs_considered = j.at("passed_alignment_pairs_considered").get<uint64_t>();
  s.pos_hash_score_distribution = j.at("pos_hash_score_distribution").get<PosHashScoreDistribution>();
}

// SequenceConversionSummmary
void to_json(json& j, const SequenceConversionSummary& s)
{
  j = json {
    {"read_length_avg", s.read_length_avg},
    {"max_qual", s.max_qual},
    {"num_reads", s.num_reads},
    {"num_original_reads", s.num_original_reads},
    {"num_bases", s.num_bases},
    {"num_original_bases", s.num_original_bases},
    {"reads", s.reads},
    {"total_reference_sequence_length", s.total_reference_sequence_length},
    {"read_length_max", s.read_length_max},
    {"read_length_min", s.read_length_min},
  };
}

void from_json(const json& j, SequenceConversionSummary& s)
{
  s.read_length_avg = j.at("read_length_avg").get<float>();
  s.max_qual = j.at("max_qual").get<uint32_t>();
  s.num_reads = j.at("num_reads").get<uint64_t>();
  s.num_original_reads = j.at("num_original_reads").get<uint64_t>();
  s.num_bases = j.at("num_bases").get<uint64_t>();
  s.num_original_bases = j.at("num_original_bases").get<uint64_t>();
  s.reads = j.at("reads").get<map<string, AnalyzeFastqSummary> >();
  s.total_reference_sequence_length = j.at("total_reference_sequence_length").get<uint64_t>();
  s.read_length_max = j.at("read_length_max").get<uint32_t>();
  s.read_length_min = j.at("read_length_min").get<uint32_t>();
}

// cErrorCountSummary
void to_json(json& j, const ErrorCountSummary& s)
{
  j = json{
    {"no_pos_hash_per_position_pr", s.no_pos_hash_per_position_pr}
  };
}

void from_json(const json& j, ErrorCountSummary& s)
{
  s.no_pos_hash_per_position_pr = j.at("no_pos_hash_per_position_pr").get<double>();
}
  
void to_json(json& j, const CoverageSummaries& s)
{
  for (CoverageSummaries::const_iterator it = s.begin(); it != s.end(); it++) {
    j[it->first] = it->second;
  }
}

void from_json(const json& j, CoverageSummaries& s)
{
  for (json::const_iterator it = j.begin(); it != j.end(); ++it) {
    s[it.key()] = it.value();
  }
}
  
void to_json(json& j, const ErrorCountSummaries& s)
{
  for (ErrorCountSummaries::const_iterator it = s.begin(); it != s.end(); it++) {
    j[it->first] = it->second;
  }
}

void from_json(const json& j, ErrorCountSummaries& s)
{
  for (json::const_iterator it = j.begin(); it != j.end(); ++it) {
    s[it.key()] = it.value();
  }
}

// Summary
void to_json(json& j, const Summary& s)
{
  j = json{
    {"sequence_conversion", s.sequence_conversion},
    {"candidate_junction", s.candidate_junction},
    {"alignment_resolution", s.alignment_resolution},
    {"preprocess_coverage", s.preprocess_coverage},
    {"unique_coverage", s.unique_coverage},
    {"preprocess_error_count", s.preprocess_error_count},
    {"preprocess_alignments", s.preprocess_alignments},
  };
}

void from_json(const json& j, Summary& s)
{
  s.alignment_resolution = j.at("alignment_resolution").get<AlignmentResolutionSummary>();
  s.preprocess_coverage = j.at("preprocess_coverage").get<CoverageSummaries>();
  s.unique_coverage = j.at("unique_coverage").get<CoverageSummaries>();
  s.preprocess_alignments = j.at("preprocess_alignments").get<PreprocessAlignmentsSummary>();
  s.candidate_junction = j.at("candidate_junction").get<CandidateJunctionSummary>();
  s.sequence_conversion = j.at("sequence_conversion").get<SequenceConversionSummary>();
  s.preprocess_error_count = j.at("preprocess_error_count").get<ErrorCountSummaries>();
}
  
PublicReadFileSummary::PublicReadFileSummary(const ReadFileSummary &rfs, const AnalyzeFastqSummary &afs)
{
  read_length_min = afs.read_length_min;
  read_length_max = afs.read_length_max;
  read_length_avg = afs.read_length_avg;
  num_original_reads = afs.num_original_reads;
  num_filtered_too_short_reads = afs.num_filtered_too_short_reads;
  num_filtered_same_base_reads = afs.num_filtered_same_base_reads;
  num_filtered_too_many_N_reads = afs.num_filtered_too_many_N_reads;
  num_reads = afs.num_reads;
  min_quality_score = afs.min_quality_score;
  max_quality_score = afs.max_quality_score;
  num_original_bases = afs.num_original_bases;
  num_bases = afs.num_bases;
  quality_format_original = afs.quality_format_original;
  quality_format = afs.quality_format;
  
  num_aligned_reads = rfs.num_total_reads - rfs.num_unmatched_reads;
  num_aligned_bases = rfs.num_total_bases - rfs.num_unmatched_bases;
  fraction_aligned_reads = static_cast<double>(num_aligned_reads) / static_cast<double>(rfs.num_total_reads);
  fraction_aligned_bases = static_cast<double>(num_aligned_bases) / static_cast<double>(rfs.num_total_bases);
}

  
PublicReadSummary::PublicReadSummary(const Summary &s)
{
  
  for(map<string,AnalyzeFastqSummary>::const_iterator it = s.sequence_conversion.reads.begin();
      it != s.sequence_conversion.reads.end(); it++) {
    
    string seq_key = it->first;
    this->read_file.emplace(seq_key, PublicReadFileSummary(s.alignment_resolution.read_file.at(seq_key), it->second));
  }
  
  total_reads = s.alignment_resolution.total_reads;
  total_bases = s.alignment_resolution.total_bases;
  total_aligned_reads = s.alignment_resolution.total_reads - s.alignment_resolution.total_unmatched_reads;
  total_aligned_bases = s.alignment_resolution.total_bases - s.alignment_resolution.total_unmatched_bases;
  total_fraction_aligned_reads = static_cast<double>(total_aligned_reads) / static_cast<double>(s.alignment_resolution.total_reads);
  total_fraction_aligned_bases = static_cast<double>(total_aligned_bases) / static_cast<double>(s.alignment_resolution.total_bases);
}

  
PublicReferenceSummary::PublicReferenceSummary(
                                               const AlignmentResolutionReferenceSummary &arrs,
                                               const CoverageSummary &cs,
                                               const cAnnotatedSequence &r,
                                               const cReferenceSequenceSettings &rss
                                               )
{
  length = r.get_sequence_length();
  num_features = r.m_features.size();
  num_genes = r.m_genes.size();
  num_repeats = r.m_repeats.size();
  
  num_reads_mapped_to_reference = arrs.reads_mapped_to_reference;
  num_bases_mapped_to_reference = arrs.bases_mapped_to_reference;
  
  coverage_deletion_coverage_propagation_cutoff = cs.deletion_coverage_propagation_cutoff;
  coverage_deletion_coverage_seed_cutoff = cs.deletion_coverage_seed_cutoff;
  coverage_nbinom_size_parameter = cs.nbinom_size_parameter;
  coverage_nbinom_mean_parameter = cs.nbinom_mean_parameter;
  coverage_nbinom_prob_parameter = cs.nbinom_prob_parameter;
  coverage_nbinom_variance = cs.nbinom_variance;
  coverage_nbinom_dispersion = cs.nbinom_dispersion;
  coverage_average = cs.average;
  coverage_variance = cs.variance;
  coverage_dispersion = cs.dispersion;
  
  coverage_group = rss.m_seq_id_to_coverage_group_map.at(r.m_seq_id);
  junction_only = rss.m_junction_only_seq_id_set.count(r.m_seq_id);
}
  
PublicReferencesSummary::PublicReferencesSummary(
                                                 const Summary &s,
                                                 const cReferenceSequences& r,
                                                 const cReferenceSequenceSettings &rss
                                                 )
{
  total_length = r.get_total_length();
  total_features = 0;
  total_genes = 0;
  total_repeats = 0;
  
  for(size_t i = 0; i != s.alignment_resolution.reference.size(); i++) {
    const cAnnotatedSequence& seq = r.at(i);
    const string& seq_id = seq.m_seq_id;
    
    total_features += seq.m_features.size();
    total_genes += seq.m_genes.size();
    total_repeats += seq.m_repeats.size();
    
    this->reference.emplace(
                            seq_id,
                            PublicReferenceSummary(
                                                   s.alignment_resolution.reference.at(i),
                                                   s.unique_coverage.at(seq_id),
                                                   seq,
                                                   rss
                                                   )
                            );
  }
}
  
PublicOptionsSummary::PublicOptionsSummary(const Settings &t)
{
  //! Settings: Workflow
  custom_run_name = t.custom_run_name;
  num_processors = t.num_processors;
  skip_read_filtering = t.skip_read_filtering;
  skip_new_junction_prediction = t.skip_new_junction_prediction;
  skip_read_alignment_and_missing_coverage_prediction = t.skip_read_alignment_and_missing_coverage_prediction;
  skip_missing_coverage_prediction = t.skip_missing_coverage_prediction;
  skip_alignment_or_plot_generation = t.skip_alignment_or_plot_generation;
  
  //! Settings: Read File
  aligned_sam_mode = t.aligned_sam_mode;
  read_file_coverage_fold_limit = t.read_file_coverage_fold_limit;
  read_file_read_length_min = t.read_file_read_length_min;
  read_file_max_same_base_fraction = t.read_file_max_same_base_fraction;
  read_file_max_N_fraction = t.read_file_max_N_fraction;
  
  //! Settings: Read Alignment
  bowtie2_scoring = t.bowtie2_scoring;
  bowtie2_stage1 = t.bowtie2_stage1;
  bowtie2_stage2 = t.bowtie2_stage2;
  bowtie2_junction = t.bowtie2_junction;
  bowtie2_junction_maximum_alignments_to_consider_per_read = t.bowtie2_junction_maximum_alignments_to_consider_per_read;
  bowtie2_genome_maximum_alignments_to_consider_per_read = t.bowtie2_genome_maximum_alignments_to_consider_per_read;
  minimum_mapping_quality = t.minimum_mapping_quality;
  require_match_length = t.require_match_length;
  require_match_fraction = t.require_match_fraction;
  maximum_read_mismatches = t.maximum_read_mismatches;
  
  //! Settings: Candidate Junction
  preprocess_junction_min_indel_split_length = t.preprocess_junction_min_indel_split_length;
  required_both_unique_length_per_side = t.required_both_unique_length_per_side;
  required_both_unique_length_per_side_fraction = t.required_both_unique_length_per_side_fraction;
  required_one_unique_length_per_side = t.required_one_unique_length_per_side;
  unmatched_end_minimum_read_length = t.unmatched_end_minimum_read_length;
  unmatched_end_length_factor = t.unmatched_end_length_factor;
  maximum_junction_sequence_insertion_length = t.maximum_junction_sequence_insertion_length;
  maximum_junction_sequence_overlap_length = t.maximum_junction_sequence_overlap_length;
  maximum_junction_sequence_negative_overlap_length_fraction = t.maximum_junction_sequence_negative_overlap_length_fraction;
  maximum_junction_sequence_negative_overlap_length_minimum = t.maximum_junction_sequence_negative_overlap_length_minimum;
  maximum_junction_sequence_positive_overlap_length_fraction = t.maximum_junction_sequence_positive_overlap_length_fraction;
  maximum_junction_sequence_positive_overlap_length_minimum = t.maximum_junction_sequence_positive_overlap_length_minimum;
  highly_redundant_junction_ignore_passed_pair_limit = t.highly_redundant_junction_ignore_passed_pair_limit;
  maximum_junction_sequence_passed_alignment_pairs_to_consider = t.maximum_junction_sequence_passed_alignment_pairs_to_consider;
  minimum_candidate_junction_pos_hash_score = t.minimum_candidate_junction_pos_hash_score;
  minimum_candidate_junctions = t.minimum_candidate_junctions;
  maximum_candidate_junctions = t.maximum_candidate_junctions;
  maximum_candidate_junction_length_factor = t.maximum_candidate_junction_length_factor;
  
  //! Settings: Alignment Resolution
  add_split_junction_sides = t.add_split_junction_sides;
  minimum_alignment_resolution_pos_hash_score = t.minimum_alignment_resolution_pos_hash_score;
  minimum_pr_no_read_start_per_position = t.minimum_pr_no_read_start_per_position;
  junction_minimum_side_match = t.junction_minimum_side_match;
  junction_pos_hash_neg_log10_p_value_cutoff = t.junction_pos_hash_neg_log10_p_value_cutoff;
  
  //! Settings: Mutation Identification
  user_evidence_genome_diff_file_name = t.user_evidence_genome_diff_file_name;
  base_quality_cutoff = t.base_quality_cutoff;
  quality_score_trim = t.quality_score_trim;
  deletion_coverage_propagation_cutoff = t.deletion_coverage_propagation_cutoff;
  deletion_coverage_seed_cutoff = t.deletion_coverage_seed_cutoff;
  polymorphism_prediction = t.polymorphism_prediction;
  mixed_base_prediction = t.mixed_base_prediction;
  targeted_sequencing = t.targeted_sequencing;
  print_mutation_identification_per_position_file = t.print_mutation_identification_per_position_file;
  mutation_log10_e_value_cutoff = t.mutation_log10_e_value_cutoff;
  consensus_minimum_variant_coverage_each_strand = t.consensus_minimum_variant_coverage_each_strand;
  consensus_minimum_total_coverage_each_strand = t.consensus_minimum_total_coverage_each_strand;
  consensus_minimum_variant_coverage = t.consensus_minimum_variant_coverage;
  consensus_minimum_total_coverage = t.consensus_minimum_total_coverage;
  consensus_frequency_cutoff = t.consensus_frequency_cutoff;
  polymorphism_log10_e_value_cutoff = t.polymorphism_log10_e_value_cutoff;
  polymorphism_bias_p_value_cutoff = t.polymorphism_bias_p_value_cutoff;
  polymorphism_frequency_cutoff = t.polymorphism_frequency_cutoff;
  polymorphism_minimum_variant_coverage_each_strand = t.polymorphism_minimum_variant_coverage_each_strand;
  polymorphism_minimum_total_coverage_each_strand = t.polymorphism_minimum_total_coverage_each_strand;
  polymorphism_minimum_variant_coverage = t.polymorphism_minimum_variant_coverage;
  polymorphism_minimum_total_coverage = t.polymorphism_minimum_total_coverage;
  polymorphism_reject_indel_homopolymer_length = t.polymorphism_reject_indel_homopolymer_length;
  polymorphism_reject_surrounding_homopolymer_length = t.polymorphism_reject_surrounding_homopolymer_length;
  no_indel_polymorphisms = t.no_indel_polymorphisms;
  
  //! Settings: Mutation Prediction
  size_cutoff_AMP_becomes_INS_DEL_mutation = t.size_cutoff_AMP_becomes_INS_DEL_mutation;
  
  //! Settings: Output
  max_displayed_reads = t.max_displayed_reads;
  no_javascript = t.no_javascript;
  header_genome_diff_file_name = t.header_genome_diff_file_name;
  max_nucleotides_to_show_in_tables = t.max_nucleotides_to_show_in_tables;
  max_rejected_read_alignment_evidence_to_show = t.max_rejected_read_alignment_evidence_to_show;
  max_rejected_junction_evidence_to_show = t.max_rejected_junction_evidence_to_show;
  hide_circular_genome_junctions = t.hide_circular_genome_junctions;
}
  
  PublicSummary::PublicSummary(const Summary &s, const Settings &t, const cReferenceSequences &r)
  : reads(s)
  , references(s, r, t.refseq_settings)
  , options(t)
{ }
  
// PublicReadFileSummary
void to_json(json& j, const PublicReadFileSummary& s)
{
  j = json{
    {"read_length_min", s.read_length_min},
    {"read_length_max", s.read_length_max},
    {"read_length_avg", s.read_length_avg},
    {"num_original_reads", s.num_original_reads},
    {"num_filtered_too_short_reads", s.num_filtered_too_short_reads},
    {"num_filtered_same_base_reads", s.num_filtered_same_base_reads},
    {"num_filtered_too_many_N_reads", s.num_filtered_too_many_N_reads},
    
    {"num_reads", s.num_reads},
    {"min_quality_score", s.min_quality_score},
    {"max_quality_score", s.max_quality_score},
    {"num_original_bases", s.num_original_bases},
    {"num_bases", s.num_bases},
    {"quality_format_original", s.quality_format_original},
    {"quality_format", s.quality_format},
    
    {"num_aligned_reads", s.num_aligned_reads},
    {"num_aligned_bases", s.num_aligned_bases},
    {"fraction_aligned_reads", s.fraction_aligned_reads},
    {"fraction_aligned_bases", s.fraction_aligned_bases},
  };
}
  
void from_json(const json& j, PublicReadFileSummary& s)
{
  s.read_length_min = j.at("read_length_min").get<uint32_t>();
  s.read_length_max = j.at("read_length_max").get<uint32_t>();
  s.read_length_avg = j.at("read_length_avg").get<double>();
  s.num_original_reads = j.at("num_original_reads").get<uint64_t>();
  s.num_filtered_too_short_reads = j.at("num_filtered_too_short_reads").get<uint64_t>();
  s.num_filtered_same_base_reads = j.at("num_filtered_same_base_reads").get<uint64_t>();
  s.num_filtered_too_many_N_reads = j.at("num_filtered_too_many_N_reads").get<uint64_t>();
  
  s.num_reads = j.at("num_reads").get<uint64_t>();
  s.min_quality_score = j.at("min_quality_score").get<uint32_t>();
  s.max_quality_score = j.at("max_quality_score").get<uint32_t>();
  s.num_original_bases = j.at("num_original_bases").get<uint64_t>();
  s.num_bases = j.at("num_bases").get<uint64_t>();
  s.quality_format_original = j.at("quality_format_original").get<string>();
  s.quality_format = j.at("quality_format").get<string>();

  s.num_aligned_reads = j.at("num_aligned_reads").get<uint64_t>();
  s.num_aligned_bases = j.at("num_aligned_bases").get<uint64_t>();
  s.fraction_aligned_reads = j.at("fraction_aligned_reads").get<double>();
  s.fraction_aligned_bases = j.at("fraction_aligned_bases").get<double>();
}
  
void to_json(json& j, const PublicReadFileSummaries& s)
{
  for (PublicReadFileSummaries::const_iterator it = s.begin(); it != s.end(); it++) {
    j[it->first] = it->second;
  }
}

void from_json(const json& j, PublicReadFileSummaries& s)
{
  for (json::const_iterator it = j.begin(); it != j.end(); ++it) {
    s[it.key()] = it.value();
  }
}
  
// PublicReadSummary
void to_json(json& j, const PublicReadSummary& s)
{
  
  j = json{
    {"total_reads", s.total_reads},
    {"total_bases", s.total_bases},
    {"total_aligned_reads", s.total_aligned_reads},
    {"total_aligned_bases", s.total_aligned_bases},
    {"total_fraction_aligned_reads", s.total_fraction_aligned_reads},
    {"total_fraction_aligned_bases", s.total_fraction_aligned_bases},
    {"read_file", s.read_file},
  };
  
}

void from_json(const json& j, PublicReadSummary& s)
{
  s.total_reads = j.at("total_reads").get<uint64_t>();
  s.total_bases = j.at("total_bases").get<uint64_t>();
  s.total_aligned_reads = j.at("total_aligned_reads").get<uint64_t>();
  s.total_aligned_bases = j.at("total_aligned_bases").get<uint64_t>();
  s.total_fraction_aligned_reads = j.at("total_fraction_aligned_reads").get<double>();
  s.total_fraction_aligned_bases = j.at("total_fraction_aligned_bases").get<double>();
  s.read_file = j.at("read_file").get<PublicReadFileSummaries>();
}
  
// PublicReferenceFileSummary
void to_json(json& j, const PublicReferenceSummary& s)
{
  j = json{
    {"length", s.length},
    {"num_features", s.num_features},
    {"num_genes", s.num_genes},
    {"num_repeats", s.num_repeats},
    
    {"num_reads_mapped_to_reference", s.num_reads_mapped_to_reference},
    {"num_bases_mapped_to_reference", s.num_bases_mapped_to_reference},
    
    {"coverage_deletion_coverage_propagation_cutoff", s.coverage_deletion_coverage_propagation_cutoff},
    {"coverage_deletion_coverage_seed_cutoff", s.coverage_deletion_coverage_seed_cutoff},
    {"coverage_nbinom_size_parameter", s.coverage_nbinom_size_parameter},
    {"coverage_nbinom_mean_parameter", s.coverage_nbinom_mean_parameter},
    {"coverage_nbinom_prob_parameter", s.coverage_nbinom_prob_parameter},
    {"coverage_nbinom_variance", s.coverage_nbinom_variance},
    {"coverage_nbinom_dispersion", s.coverage_nbinom_dispersion},
    {"coverage_average", s.coverage_average},
    {"coverage_variance", s.coverage_variance},
    {"coverage_dispersion", s.coverage_dispersion},
    
    {"coverage_group", s.coverage_group},
    {"junction_only", s.junction_only},
  };
}
  
void from_json(const json& j, PublicReferenceSummary& s)
{
  s.length = j.at("length").get<uint64_t>();
  s.num_features = j.at("num_features").get<uint64_t>();
  s.num_genes = j.at("num_genes").get<uint64_t>();
  s.num_repeats = j.at("num_repeats").get<uint64_t>();
  
  s.num_reads_mapped_to_reference = j.at("num_reads_mapped_to_reference").get<double>();
  s.num_bases_mapped_to_reference = j.at("num_bases_mapped_to_reference").get<double>();
  
  s.coverage_deletion_coverage_propagation_cutoff = j.at("coverage_deletion_coverage_propagation_cutoff").get<double>();
  s.coverage_deletion_coverage_seed_cutoff = j.at("coverage_deletion_coverage_seed_cutoff").get<double>();
  s.coverage_nbinom_size_parameter = j.at("coverage_nbinom_size_parameter").get<double>();
  s.coverage_nbinom_mean_parameter = j.at("coverage_nbinom_mean_parameter").get<double>();
  s.coverage_nbinom_prob_parameter = j.at("coverage_nbinom_prob_parameter").get<double>();
  s.coverage_nbinom_variance = j.at("coverage_nbinom_variance").get<double>();
  s.coverage_nbinom_dispersion = j.at("coverage_nbinom_dispersion").get<double>();
  s.coverage_average = j.at("coverage_average").get<double>();
  s.coverage_variance = j.at("coverage_variance").get<double>();
  s.coverage_dispersion = j.at("coverage_dispersion").get<double>();
  
  
  s.coverage_group = j.at("coverage_group").get<int32_t>();
  s.junction_only = j.at("junction_only").get<bool>();
}
  
// PublicReferenceSummaries
void to_json(json& j, const PublicReferenceSummaries& s)
{
  for (PublicReferenceSummaries::const_iterator it = s.begin(); it != s.end(); it++) {
    j[it->first] = it->second;
  }
}
  
void from_json(const json& j, PublicReferenceSummaries& s)
{
  for (json::const_iterator it = j.begin(); it != j.end(); ++it) {
    s[it.key()] = it.value();
  }
}

// PublicReferencesSummary
void to_json(json& j, const PublicReferencesSummary& s)
{
  j = json{
    {"total_length", s.total_length},
    {"total_features", s.total_features},
    {"total_genes", s.total_genes},
    {"total_repeats", s.total_repeats},
    
    {"reference", s.reference},
  };
}
  
void from_json(const json& j, PublicReferencesSummary& s)
{
  s.total_length = j.at("total_length").get<uint64_t>();
  s.total_features = j.at("total_features").get<uint64_t>();
  s.total_genes = j.at("total_genes").get<uint64_t>();
  s.total_repeats = j.at("total_repeats").get<uint64_t>();
  
  s.reference = j.at("reference").get<PublicReferenceSummaries>();
}


// PublicOptionsSummary
void to_json(json& j, const PublicOptionsSummary& s)
{
  
  j = json{
    //! Settings: Workflow
    {"workflow", json{
      {"custom_run_name", s.custom_run_name},
      {"num_processors", s.num_processors},
      {"skip_read_filtering", s.skip_read_filtering},
      {"skip_new_junction_prediction", s.skip_new_junction_prediction},
      {"skip_read_alignment_and_missing_coverage_prediction", s.skip_read_alignment_and_missing_coverage_prediction},
      {"skip_missing_coverage_prediction", s.skip_missing_coverage_prediction},
      {"skip_alignment_or_plot_generation", s.skip_alignment_or_plot_generation},
      }
    },
    
    //! Settings: Read File
    {"read_file", json{
      {"aligned_sam_mode", s.aligned_sam_mode},
      {"read_file_coverage_fold_limit", s.read_file_coverage_fold_limit},
      {"read_file_read_length_min", s.read_file_read_length_min},
      {"read_file_max_same_base_fraction", s.read_file_max_same_base_fraction},
      {"read_file_max_N_fraction", s.read_file_max_N_fraction},
      }
    },
    
    //! Settings: Read Alignment
    {"read_alignment", json{
      {"bowtie2_scoring", s.bowtie2_scoring},
      {"bowtie2_stage1", s.bowtie2_stage1},
      {"bowtie2_stage2", s.bowtie2_stage2},
      {"bowtie2_junction", s.bowtie2_junction},
      {"bowtie2_junction_maximum_alignments_to_consider_per_read", s.bowtie2_junction_maximum_alignments_to_consider_per_read},
      {"bowtie2_genome_maximum_alignments_to_consider_per_read", s.bowtie2_genome_maximum_alignments_to_consider_per_read},
      {"minimum_mapping_quality", s.minimum_mapping_quality},
      {"require_match_length", s.require_match_length},
      {"require_match_fraction", s.require_match_fraction},
      {"maximum_read_mismatches", s.maximum_read_mismatches},
      }
    },
    
    //! Settings: Candidate Junction
    {"candidate_junction", json{
      {"preprocess_junction_min_indel_split_length", s.preprocess_junction_min_indel_split_length},
      {"required_both_unique_length_per_side", s.required_both_unique_length_per_side},
      {"required_both_unique_length_per_side_fraction", s.required_both_unique_length_per_side_fraction},
      {"required_one_unique_length_per_side", s.required_one_unique_length_per_side},
      {"unmatched_end_minimum_read_length", s.unmatched_end_minimum_read_length},
      {"unmatched_end_length_factor", s.unmatched_end_length_factor},
      {"maximum_junction_sequence_insertion_length", s.maximum_junction_sequence_insertion_length},
      {"maximum_junction_sequence_overlap_length", s.maximum_junction_sequence_overlap_length},
      {"maximum_junction_sequence_negative_overlap_length_fraction", s.maximum_junction_sequence_negative_overlap_length_fraction},
      {"maximum_junction_sequence_negative_overlap_length_minimum", s.maximum_junction_sequence_negative_overlap_length_minimum},
      {"maximum_junction_sequence_positive_overlap_length_fraction", s.maximum_junction_sequence_positive_overlap_length_fraction},
      {"maximum_junction_sequence_positive_overlap_length_minimum", s.maximum_junction_sequence_positive_overlap_length_minimum},
      {"highly_redundant_junction_ignore_passed_pair_limit", s.highly_redundant_junction_ignore_passed_pair_limit},
      {"maximum_junction_sequence_passed_alignment_pairs_to_consider", s.maximum_junction_sequence_passed_alignment_pairs_to_consider},
      {"minimum_candidate_junction_pos_hash_score", s.minimum_candidate_junction_pos_hash_score},
      {"minimum_candidate_junctions", s.minimum_candidate_junctions},
      {"maximum_candidate_junctions", s.maximum_candidate_junctions},
      {"maximum_candidate_junction_length_factor", s.maximum_candidate_junction_length_factor},
      }
    },
    
    //! Settings: Alignment Resolution
    {"alignment_resolution", json{
      {"add_split_junction_sides", s.add_split_junction_sides},
      {"minimum_alignment_resolution_pos_hash_score", s.minimum_alignment_resolution_pos_hash_score},
      {"minimum_pr_no_read_start_per_position", s.minimum_pr_no_read_start_per_position},
      {"junction_minimum_side_match", s.junction_minimum_side_match},
      {"junction_pos_hash_neg_log10_p_value_cutoff", s.junction_pos_hash_neg_log10_p_value_cutoff},
      }
    },
    
    //! Settings: Mutation Identification
    {"mutation_identification", json{
      {"user_evidence_genome_diff_file_name", s.user_evidence_genome_diff_file_name},
      {"base_quality_cutoff", s.base_quality_cutoff},
      {"quality_score_trim", s.quality_score_trim},
      {"deletion_coverage_propagation_cutoff", s.deletion_coverage_propagation_cutoff},
      {"deletion_coverage_seed_cutoff", s.deletion_coverage_seed_cutoff},
      {"polymorphism_prediction", s.polymorphism_prediction},
      {"mixed_base_prediction", s.mixed_base_prediction},
      {"targeted_sequencing", s.targeted_sequencing},
      {"print_mutation_identification_per_position_file", s.print_mutation_identification_per_position_file},
      {"mutation_log10_e_value_cutoff", s.mutation_log10_e_value_cutoff},
      {"consensus_minimum_variant_coverage_each_strand", s.consensus_minimum_variant_coverage_each_strand},
      {"consensus_minimum_total_coverage_each_strand", s.consensus_minimum_total_coverage_each_strand},
      {"consensus_minimum_variant_coverage", s.consensus_minimum_variant_coverage},
      {"consensus_minimum_total_coverage", s.consensus_minimum_total_coverage},
      {"consensus_frequency_cutoff", s.consensus_frequency_cutoff},
      {"consensus_reject_indel_homopolymer_length", s.consensus_reject_indel_homopolymer_length},
      {"consensus_reject_surrounding_homopolymer_length", s.consensus_reject_surrounding_homopolymer_length},
      {"polymorphism_log10_e_value_cutoff", s.polymorphism_log10_e_value_cutoff},
      {"polymorphism_bias_p_value_cutoff", s.polymorphism_bias_p_value_cutoff},
      {"polymorphism_frequency_cutoff", s.polymorphism_frequency_cutoff},
      {"polymorphism_minimum_variant_coverage_each_strand", s.polymorphism_minimum_variant_coverage_each_strand},
      {"polymorphism_minimum_total_coverage_each_strand", s.polymorphism_minimum_total_coverage_each_strand},
      {"polymorphism_minimum_variant_coverage", s.polymorphism_minimum_variant_coverage},
      {"polymorphism_minimum_total_coverage", s.polymorphism_minimum_total_coverage},      
      {"polymorphism_reject_indel_homopolymer_length", s.polymorphism_reject_indel_homopolymer_length},
      {"polymorphism_reject_surrounding_homopolymer_length", s.polymorphism_reject_surrounding_homopolymer_length},
      {"no_indel_polymorphisms", s.no_indel_polymorphisms},
      }
    },
    
    //! Settings: Mutation Prediction
    {"mutation_prediction", json{
      {"size_cutoff_AMP_becomes_INS_DEL_mutation", s.size_cutoff_AMP_becomes_INS_DEL_mutation},
      }
    },
    
    //! Settings: Output
    {"output", json{
      {"max_displayed_reads", s.max_displayed_reads},
      {"no_javascript", s.no_javascript},
      {"header_genome_diff_file_name", s.header_genome_diff_file_name},
      {"max_nucleotides_to_show_in_tables", s.max_nucleotides_to_show_in_tables},
      {"max_rejected_read_alignment_evidence_to_show", s.max_rejected_read_alignment_evidence_to_show},
      {"max_rejected_junction_evidence_to_show", s.max_rejected_junction_evidence_to_show},
      {"hide_circular_genome_junctions", s.hide_circular_genome_junctions},
      }
    },
  };
}
  
void from_json(const json& j, PublicOptionsSummary& s)
{
  (void)j;
  (void)s;
  ERROR("Not implemented");
}


// PublicSummary
void to_json(json& j, const PublicSummary& s)
{
  j = json{
    {"reads", s.reads},
    {"references", s.references},
    {"options", s.options},
  };
}

  
void from_json(const json& j, PublicSummary& s)
{
  s.reads = j.at("reads").get<PublicReadSummary>();
  s.references = j.at("references").get<PublicReferencesSummary>();
  s.options = j.at("options").get<PublicOptionsSummary>();
}

  
}
