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
    {"reads_mapped_to_references", s.reads_mapped_to_references},
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
  s.reads_mapped_to_references = j.at("reads_mapped_to_references").get<vector<double> >();
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
    {"min_read_length", s.min_read_length},
    {"max_read_length", s.max_read_length},
    {"avg_read_length", s.avg_read_length},
    {"original_reads", s.original_reads},
    {"too_short_filtered_reads", s.too_short_filtered_reads},
    {"same_base_filtered_reads", s.same_base_filtered_reads},
    {"N_filtered_reads", s.N_filtered_reads},
    {"num_reads", s.num_reads},
    {"min_quality_score", s.min_quality_score},
    {"max_quality_score", s.max_quality_score},
    {"original_num_bases", s.original_num_bases},
    {"num_bases", s.num_bases},
    {"original_qual_format", s.original_qual_format},
    {"quality_format", s.quality_format},
    {"converted_fastq_name", s.converted_fastq_name},
  };
}

void from_json(const json& j, AnalyzeFastqSummary& s)
{
  s.min_read_length = j.at("min_read_length").get<uint32_t>();
  s.max_read_length = j.at("max_read_length").get<uint32_t>();
  s.avg_read_length = j.at("avg_read_length").get<double>();
  s.original_reads = j.at("original_reads").get<uint64_t>();
  s.too_short_filtered_reads = j.at("too_short_filtered_reads").get<uint64_t>();
  s.same_base_filtered_reads = j.at("same_base_filtered_reads").get<uint64_t>();
  s.N_filtered_reads = j.at("N_filtered_reads").get<uint64_t>();
  s.num_reads = j.at("num_reads").get<uint64_t>();
  s.min_quality_score = j.at("min_quality_score").get<uint32_t>();
  s.max_quality_score = j.at("max_quality_score").get<uint32_t>();
  s.original_num_bases = j.at("original_num_bases").get<uint64_t>();
  s.num_bases = j.at("num_bases").get<uint64_t>();
  s.original_qual_format = j.at("original_qual_format").get<string>();
  s.quality_format = j.at("quality_format").get<string>();
  s.converted_fastq_name = j.at("converted_fastq_name").get<string>();
}
  
// PreprocessAlignmentsSummary
void to_json(json& j, const PreprocessAlignmentsSummary& s)
{
  j = json{
    {"aligned_reads", s.aligned_reads},
    {"aligned_reads", s.alignments},
    {"aligned_reads", s.alignments_split_on_indels},
    {"aligned_reads", s.reads_with_alignments_split_on_indels},
    {"aligned_reads", s.split_alignments},
    {"aligned_reads", s.reads_with_split_alignments},
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
    {"avg_read_length", s.avg_read_length},
    {"max_qual", s.max_qual},
    {"num_reads", s.num_reads},
    {"original_num_reads", s.original_num_reads},
    {"num_bases", s.num_bases},
    {"original_num_bases", s.original_num_bases},
    {"reads", s.reads},
    {"total_reference_sequence_length", s.total_reference_sequence_length},
    {"max_read_length", s.max_read_length},
    {"min_read_length", s.min_read_length},
  };
}

void from_json(const json& j, SequenceConversionSummary& s)
{
  s.avg_read_length = j.at("avg_read_length").get<float>();
  s.max_qual = j.at("max_qual").get<uint32_t>();
  s.num_reads = j.at("num_reads").get<uint64_t>();
  s.original_num_reads = j.at("original_num_reads").get<uint64_t>();
  s.num_bases = j.at("num_bases").get<uint64_t>();
  s.original_num_bases = j.at("original_num_bases").get<uint64_t>();
  s.reads = j.at("reads").get<map<string, AnalyzeFastqSummary> >();
  s.total_reference_sequence_length = j.at("total_reference_sequence_length").get<uint64_t>();
  s.max_read_length = j.at("max_read_length").get<uint32_t>();
  s.min_read_length = j.at("min_read_length").get<uint32_t>();
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

  
}
