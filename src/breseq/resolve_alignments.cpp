/*****************************************************************************

 AUTHORS

   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors

 LICENSE AND COPYRIGHT

   Copyright (c) 2008-2010 Michigan State University
   Copyright (c) 2011-2025 The University of Texas at Austin
   Copyright (c) 2025-     Michigan State University

   breseq is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2, or (at your option) any later version.

   SPDX-License-Identifier: GPL-2.0-or-later

*****************************************************************************/

#include "resolve_alignments.h"

#include "genome_diff.h"
#include "fastq.h"
#include "fasta.h"
#include "alignment.h"
#include "identify_mutations.h"
#include "reference_sequence.h"
#include "stats.h"
#include "output.h"

using namespace std;

namespace breseq {

uint32_t qmissing (double tail_value, double pr_missing)
{
  int32_t missing = 0;
  double test_pr = 1;
  double pr_no_cov = pow(pr_missing, 2);
  while (test_pr > tail_value) {
    missing++;
    test_pr *= pr_no_cov;
  }
  return missing;
}
  
// Continuation is how many bases match the exact same past where the junction is in the reference
// We need to count this for cases of short duplications and deletions in tandem repeats to not
// penalize them when we count the "evenness" score and relative coverage.

// Graphical Explanation of the Types of Continuation:
//
// Junction:    ++++++++++ &&&&&&&&&&  side_1_junction_seq        side_2_junction_seq
// Reference 1: ++++++++++ ??????????  side_1_reference_left_seq  side_1_reference_right_seq
// Reference 2: ^^^^^^^^^^ &&&&&&&&&&  side_2_reference_left_seq  side_2_reference_right_seq
//
// ++++++++++   left_junction_seq
// &&&&&&&&&&   right_junction_seq
// ^^^^^^^^^^   left_reference_sequence
// ??????????   right_reference_sequence

// side_1_continuation  +++ = ^^^ (left side)
// side_2_continuation  &&& = ??? (right side)

void calculate_continuation(
                            ResolveJunctionInfo& rji, 
                            cReferenceSequences& ref_seq_info, 
                            cReferenceSequences& junction_ref_seq_info, 
                            uint32_t& side_1_continuation,
                            uint32_t& side_2_continuation
                            )
{
  bool verbose = false;
  
  // At this point we are an object that has been initialized with junction information
  ASSERT(rji.sides[0].seq_id != "", "Uninitialized variable");
  
  //the first part of the junction
  
  int32_t left_compare_length = rji.flanking_left + max(rji.overlap, 0);
  string left_junction_seq = junction_ref_seq_info.get_sequence_1(rji.key, 1, left_compare_length) + rji.unique_read_sequence;
  
  int32_t right_compare_length = rji.flanking_right + max(rji.overlap, 0);
  int32_t junction_sequence_length = junction_ref_seq_info[rji.key].get_sequence_length();
  string right_junction_seq = rji.unique_read_sequence + junction_ref_seq_info.get_sequence_1(rji.key, junction_sequence_length - right_compare_length + 1, junction_sequence_length);

  if (verbose) {
    cout << ">>>" << rji.key << endl;
    cout << "LEFT | OVERLAP | RIGHT = " << left_compare_length << " " << rji.overlap << " " << right_compare_length << endl;
    cout << "JUNCTION LEFT   :: " << left_junction_seq << endl;
    cout << "JUNCTION RIGHT  :: " << right_junction_seq << endl;
  }

  // the same length of reference sequence extending past the matched part of the junction
  string right_reference_seq;
  int32_t seq_length_0 = ref_seq_info[rji.sides[0].seq_id].get_sequence_length();

  if (rji.sides[0].strand == -1) {
  
    int32_t start_pos = min(seq_length_0, rji.sides[0].position + 1);
    int32_t end_pos = min(seq_length_0, rji.sides[0].position + right_compare_length);
    
    if (rji.sides[0].position + 1 <= seq_length_0) {
      right_reference_seq = ref_seq_info.get_sequence_1(rji.sides[0].seq_id, start_pos, end_pos);
    }
  }  else if (rji.sides[0].strand == +1) {
  
    int32_t start_pos = max(1, rji.sides[0].position - left_compare_length);
    int32_t end_pos = max(1, rji.sides[0].position - 1);
    
    if (rji.sides[0].position - 1 >= 1) {
      right_reference_seq = ref_seq_info.get_sequence_1(rji.sides[0].seq_id, start_pos, end_pos);
      right_reference_seq = reverse_complement(right_reference_seq);
    }
    
  }
  
  string left_reference_seq;
  int32_t seq_length_1 = ref_seq_info[rji.sides[1].seq_id].get_sequence_length();

  if (rji.sides[1].strand == +1) {
    
    int32_t start_pos = max(1, rji.sides[1].position - right_compare_length);
    int32_t end_pos = max(1, rji.sides[1].position - 1);
    
    if (rji.sides[1].position - 1 >= 1) {
      left_reference_seq = ref_seq_info.get_sequence_1(rji.sides[1].seq_id, start_pos, end_pos);
    }
    
  }  else if (rji.sides[1].strand == -1) {
    
    int32_t start_pos = min(seq_length_1, rji.sides[1].position + 1);
    int32_t end_pos = min(seq_length_1, rji.sides[1].position + right_compare_length);
    
    if (rji.sides[1].position + 1 <= seq_length_1) {
      left_reference_seq = ref_seq_info.get_sequence_1(rji.sides[1].seq_id, start_pos, end_pos);
      left_reference_seq = reverse_complement(left_reference_seq);
    }
  }
  
  if (verbose) {
    cout << "REFERENCE LEFT  :: " << left_reference_seq << endl;
    cout << "REFERENCE RIGHT :: " << right_reference_seq << endl;
  }
  
  // Compare left side going backwards
  side_2_continuation = 0;
  for (uint32_t i = 0; i < left_reference_seq.size() ; i++) {
    if (left_reference_seq[left_reference_seq.size() - i - 1] != left_junction_seq[left_junction_seq.size() - i - 1] )
        break;
    side_2_continuation++;
  }
  
  // Compare right side going forwards
  side_1_continuation = 0;
  for (uint32_t i = 0; i < right_reference_seq.size() ; i++) {
    if (right_reference_seq[i] != right_junction_seq[i] )
        break;
    side_1_continuation++;
  }
  
  if (verbose) {
    cout << "SIDE 1 CONTINUATION (LEFT-TO-RIGHT) :: " << side_1_continuation << endl;
    cout << "SIDE 2 CONTINUATION (RIGHT-TO-LEFT) :: " << side_2_continuation << endl;
  }
  
}
  
PosHashProbabilityTable::PosHashProbabilityTable(Summary& summary, const Settings& settings)
{
  average_read_length = static_cast<int32_t>(round(summary.sequence_conversion.read_length_avg));
  for (map<string,CoverageSummary>::iterator it=summary.preprocess_coverage.begin();
       it != summary.preprocess_coverage.end(); it++) {
    
    string seq_id = it->first;
    CoverageSummary& cov = it->second;
    
    double no_pos_hash_per_position_pr = summary.preprocess_error_count[seq_id].no_pos_hash_per_position_pr;
    if (no_pos_hash_per_position_pr < settings.minimum_pr_no_read_start_per_position) {
      no_pos_hash_per_position_pr = settings.minimum_pr_no_read_start_per_position;
    }
    
    Parameters p = {
          cov.nbinom_size_parameter,
          cov.nbinom_prob_parameter,
          no_pos_hash_per_position_pr,
          cov.nbinom_mean_parameter
    };
    
    param[it->first] = p;
  }
}

double PosHashProbabilityTable::probability(string& seq_id, uint32_t pos_hash_score, uint32_t max_pos_hash_score)
{
  bool verbose = false;
  
  // Use hashed value  if it exists.
  if (
      probability_table.count(seq_id) 
      && probability_table[seq_id].count(pos_hash_score) 
      && probability_table[seq_id][pos_hash_score].count(max_pos_hash_score)
      )
    return probability_table[seq_id][pos_hash_score][max_pos_hash_score];
  
  Parameters& p = param[seq_id];

  // no coverage was fit for this fragment, allow it to pass for this fragment no matter what by passing high value
  // @JEB we may want to disallow these junctions instead, by passing a negative value here?
  if (p.negative_binomial_size == 0) {
    if (verbose) cout << "No coverage fit for fragment: " << seq_id << endl;
    return 999999;
  }
  // Calculate this entry in the table -- 
  uint32_t max_coverage = static_cast<uint32_t>(round(10*p.average_coverage));
        
  double pr = 0;
  
  if (verbose) {
    cout << "Calculating: seq_id " << seq_id << " pos_hash_score " << pos_hash_score << " max_pos_hash_score " << max_pos_hash_score << endl;
    cout << "Coverage from 1 to " << max_coverage << endl;
    cout << "Negative Binomial Fit: Size = " << p.negative_binomial_size << " Prob = " << p.negative_binomial_prob << endl;
  }
  
  // This calculation can be incredibly slow when there is very high coverage 100,000s
  // estimate using this many equal sized bins.
  const   int32_t target_num_bins_for_estimate = 10000;
  int32_t bin_size = trunc(static_cast<double>(max_coverage)/target_num_bins_for_estimate)+1;
  
  for (uint32_t this_coverage=bin_size; this_coverage<= max_coverage; this_coverage+=bin_size) {
    
    // This calculation takes care of the bin size, we are getting the
    // full probabiliy all the way in the range from this_coverage to this_coverage + bin_size - 1
    double this_cov_pr =  nbdtr(this_coverage, p.negative_binomial_size, p.negative_binomial_prob)
                        - nbdtr(this_coverage-bin_size, p.negative_binomial_size, p.negative_binomial_prob);

    // This calculation uses the middle coverage value in the bin as an estimate for the
    // probability across the entire bin
    double this_coverage_middle = this_coverage + (bin_size-1) / 2;
    double this_ratio_of_coverage_to_average = this_coverage_middle / static_cast<double>(p.average_coverage);
    double this_chance_per_pos_strand_read_start = 1 - pow(p.chance_per_pos_strand_no_read_start, this_ratio_of_coverage_to_average);

    double this_pos_hash_pr = 0;
    for (uint32_t i=0; i <= pos_hash_score; i++) {

      //chance of getting pos_hash_score or lower
      this_pos_hash_pr += binomial(this_chance_per_pos_strand_read_start, max_pos_hash_score, i);
    }
        
    double this_pr = this_cov_pr*this_pos_hash_pr;
    pr += this_pr;
    
    if (verbose) {
      cout << "  Cov: " << this_coverage << " Cov Pr: " << this_cov_pr << " Pos Hash Pr: " << this_pos_hash_pr << " Read Start Pr: " << this_chance_per_pos_strand_read_start << endl;
      cout << "  This Pr: " << this_pr << " Cumulative Pr: " << pr << endl;
    }    
  }
  double log_pr = -log(pr)/log(10);
  probability_table[seq_id][pos_hash_score][max_pos_hash_score] = log_pr;
  
  if (verbose) {
    cout << "  -Log10 Cumulative Pr: " << log_pr << endl;
  }
  
  return log_pr;
}

  
// Compares matches to candidate junctions with matches to original genome
void resolve_alignments(
                        Settings& settings,
                        Summary& summary,
                        cReferenceSequences& ref_seq_info,
                        bool junction_prediction,
                        cReadFileSets& read_files
                        ) 
{    
	bool verbose = false;
  
  // local variables for convenience
  int32_t read_length_avg = static_cast<int32_t>(round(summary.sequence_conversion.read_length_avg));
  
  // Initialize out summary so we can add to it
  summary.alignment_resolution.reference.resize(ref_seq_info.size(), AlignmentResolutionReferenceSummary());
  
  // Load the reference sequence trims, for writing resolved alignments
  SequenceTrimsList trims_list;
  read_trims(trims_list, ref_seq_info, settings.reference_trim_file_name);
  
  // Junction sequence trims are loaded below, once junction_ref_seq_info is available
  SequenceTrimsList junction_trims_list;

	// ####
	// ##	Junction sequences
	// ####

	//## if there were no candidate junctions (file is empty) then we seg fault if we try to use samtools on it...
  cReferenceSequences junction_ref_seq_info(false);
	if (junction_prediction
			&& !file_exists(settings.candidate_junction_fasta_file_name.c_str())
			&& !file_empty(settings.candidate_junction_fasta_file_name.c_str())
		)
  {
		junction_prediction = 0;
  }
    
	vector<ResolveJunctionInfo> junction_info_list;
    

  //## Preload all of the information about junctions
  //## so that we only have to split the names once
  
  if (junction_prediction) {
    junction_ref_seq_info.LoadFiles(make_vector<string>(settings.candidate_junction_fasta_file_name));

    // Load the candidate junction sequence trims, for writing resolved alignments
    read_trims(junction_trims_list, junction_ref_seq_info, settings.candidate_junction_trim_file_name);

		string junction_sam_file_name = settings.file_name(settings.candidate_junction_sam_file_name, "#", read_files[0].m_files[0].m_base_name);
		bam_file junction_tam(junction_sam_file_name, settings.candidate_junction_fasta_file_name, ios::in);

    junction_info_list.resize(junction_tam.bam_header->n_targets);
    
		for (int i = 0; i < junction_tam.bam_header->n_targets; i++) {
			junction_info_list[i] = ResolveJunctionInfo(junction_tam.bam_header->target_name[i]);
		}
    if (verbose) cout << "Number of candidate junctions: " << junction_info_list.size() << endl;
  }

	// ####
	// ## Output files
	// ####

	cGenomeDiff gd;
    
  bam_file resolved_reference_tam(settings.resolved_reference_sam_file_name, settings.reference_fasta_file_name, ios::out);
  bam_file resolved_junction_tam(settings.resolved_junction_sam_file_name, settings.candidate_junction_fasta_file_name, ios::out);
  
  settings.track_intermediate_file(settings.bam_done_file_name, settings.resolved_reference_sam_file_name);
  settings.track_intermediate_file(settings.bam_done_file_name, settings.resolved_junction_sam_file_name);
  
  UniqueJunctionMatchMap unique_junction_match_map;    // Map of junction_id to MatchedJunction
	RepeatJunctionMatchMap repeat_junction_match_map;  // Map of junction_id to read_name to MatchedJunction

  // stores all junction ids that we have encountered
  map<string,uint32_t> all_junction_ids;
  
  if (!settings.aligned_sam_mode) {
  
    load_junction_alignments(
                            settings, 
                            summary, 
                            read_files, 
                            ref_seq_info,
                            junction_ref_seq_info,
                            trims_list,
                            all_junction_ids,
                            junction_prediction,
                            junction_info_list,
                            unique_junction_match_map,
                            repeat_junction_match_map,
                            resolved_reference_tam
                            );
  } else {
    
    load_sam_only_alignments(
                             settings, 
                             summary, 
                             read_files, 
                             ref_seq_info,
                             trims_list,
                             resolved_reference_tam
                             );
  }
  
  // Be sure to add user defined junctions
  if (junction_prediction) {
    if (settings.user_evidence_genome_diff_file_name != "") {
      cFastaFile ff(settings.candidate_junction_fasta_file_name, ios::in);
      cFastaSequence sequence;
      while( ff.read_sequence(sequence) ) {
        all_junction_ids[sequence.get_name()]++;
      }
    }
  }
  
  if (verbose)
  {
    cout << "Total junction ids: " << all_junction_ids.size() << endl;
  }
  
	////
	// Determine which junctions are real, prefer ones with most matches
	////
  
	list<JunctionTestInfo> junction_test_info_list; // scoring information about junctions
	list<JunctionTestInfo> passed_junction_test_info_list;
	list<JunctionTestInfo> rejected_junction_test_info_list;
  
  ////
  // Score all of the matches.
  ////
  
  vector<string> junction_ids = get_keys(all_junction_ids);
  for(vector<string>::iterator it=junction_ids.begin(); it != junction_ids.end(); it++) {
    const string& junction_id = *it;
    JunctionTestInfo junction_test_info;
    score_junction(
                   settings, 
                   summary, 
                   junction_id, 
                   unique_junction_match_map, 
                   repeat_junction_match_map,
                   resolved_junction_tam,
                   junction_test_info,
                   junction_info_list,
                   ref_seq_info,
                   junction_ref_seq_info
                   );
    junction_test_info_list.push_back(junction_test_info);

    
    if (verbose) 
    {
      cout << "Scoring Junction: " << junction_id << endl;
      cout << "  Pos hash score: " << junction_test_info.pos_hash_score << endl;
      cout << "  Number of unique matches: " << unique_junction_match_map[junction_id].size() << endl;
      size_t num_degenerate_matches = repeat_junction_match_map.count(junction_id) ? repeat_junction_match_map[junction_id].size() : 0;
      cout << "  Number of degenerate matches: " << num_degenerate_matches << endl;
    }
/*    
 // Immediately reject and resolve junctions that have no overlap matches??
 // @JEB needs testing...

    if (test_info.total_non_overlap_reads > 0) {
      junction_test_info_list.push_back(test_info);
    }
    else {
      resolve_junction(
                       settings,
                       summary,
                       ref_seq_info,
                       junction_ref_seq_info,
                       trims_list,
                       junction_trims_list,
                       test_info.junction_id,
                       unique_junction_match_map,
                       repeat_junction_match_map,
                       *reference_tam,
                       *junction_tam,
                       true, // it failed
                       test_info.total_non_overlap_reads > 0 // no non_overlap_alignments
                       ); 
      all_junction_ids.erase(test_info.junction_id);
    }
*/
  }
  
	///
	// Candidate junctions with unique matches
	///

  if (verbose) {
    cout << ">>>>> BEFORE SORT >>>>>" << endl;
    for(list<JunctionTestInfo>::iterator it = junction_test_info_list.begin(); it != junction_test_info_list.end(); it++) {
      JunctionTestInfo& junction_test_info = *it;
      string key = junction_test_info.junction_id;
      cout << key << endl;
      cout << "  Pos hash score: " << junction_test_info.pos_hash_score << endl;
      cout << "  Unique matches: " << junction_test_info.unique_matches_size << endl;
      cout << "  Repeat matches: " << junction_test_info.repeat_matches_size << endl;
    }
  }
  
	//sort junction ids from lowest to highest pos_hash score
  junction_test_info_list.sort();
  //junction_test_info_list.reverse();

  
  if (verbose) {
    cout << ">>>>> AFTER SORT >>>>>" << endl;
    for(list<JunctionTestInfo>::iterator it = junction_test_info_list.begin(); it != junction_test_info_list.end(); it++) {
      JunctionTestInfo& junction_test_info = *it;
      string key = junction_test_info.junction_id;
      cout << endl << key << endl;
      cout << "  Pos hash score: " << junction_test_info.pos_hash_score << endl;
      cout << "  Unique matches: " << junction_test_info.unique_matches_size << endl;
      cout << "  Repeat matches: " << junction_test_info.repeat_matches_size << endl;
    }
  }
  
  PosHashProbabilityTable pos_hash_p_value_calculator(summary, settings);
  
  while(!junction_test_info_list.empty() ) {
    
    JunctionTestInfo& junction_test_info = junction_test_info_list.back();
    const string& junction_id = junction_test_info.junction_id;
    
    // We need to re-score because repeat matches may have been taken by better junctions
    score_junction(
                   settings, 
                   summary, 
                   junction_id, 
                   unique_junction_match_map, 
                   repeat_junction_match_map,
                   resolved_junction_tam,
                   junction_test_info,
                   junction_info_list,
                   ref_seq_info,
                   junction_ref_seq_info
                   );
    
    
    ResolveJunctionInfo junction_info(junction_id);
    
    // Test the best-scoring junction.
    bool failed = false;  
    bool record = true;
    
    ///////////////////////////////////////////////////////////////
    // Tests related to the pos hash score (coverage evenness)
    
    // One can actually have a passing E-value score with a pos_hash score of zero under some circumstances...
    if (junction_test_info.pos_hash_score == 0) {
      junction_test_info.reject_reasons.push_back("COVERAGE_EVENNESS_SKEW");
    }
    
    // Check that it met the minimum pos hash score criterion
    if (junction_test_info.pos_hash_score < settings.minimum_alignment_resolution_pos_hash_score) {
      junction_test_info.reject_reasons.push_back("COVERAGE_EVENNESS_SKEW");
    }
    
    // We don't even record junctions in the genome diff if they don't meet the minimum_alignment_resolution_pos_hash_score
    // criterion, or if they just have no matches, unless they are user-defined junction
    record = (junction_test_info.reject_reasons.size() == 0) || junction_info.user_defined;
    
    
    // If both are on a junction-only sequence then don't count it -- EVEN IF USER DEFINED
    // This is purposeful after re-accepting user_defined junctions
    if ( settings.junction_only_seq_id_set().count(junction_info.sides[0].seq_id) && settings.junction_only_seq_id_set().count(junction_info.sides[1].seq_id) ) {
      junction_test_info.reject_reasons.push_back("BETWEEN_TWO_JUNCTION_ONLY_SEQUENCES");
    }
    
    // If no two reads have different start and end values from each other, regardless of strand, then fail (Chuck it!)
    if (!junction_test_info.has_reads_with_both_different_start_and_end) {
      junction_test_info.reject_reasons.push_back("COVERAGE_EVENNESS_SKEW");
    }
    
      junction_test_info.neg_log10_pos_hash_p_value = -1;
    
    // In targeted mode it can be rejected for the reasons above, but we don't have a coverage distribution so we can't reject it on that basis
    // Zero junction_pos_hash_neg_log10_p_value_cutoff means calculating this p-value is off
    if (!settings.targeted_sequencing && settings.junction_pos_hash_neg_log10_p_value_cutoff) {
      
      double neg_log10_p_value_1 = pos_hash_p_value_calculator.probability(junction_info.sides[0].seq_id, junction_test_info.pos_hash_score, junction_test_info.max_pos_hash_score);
      double neg_log10_p_value_2 = pos_hash_p_value_calculator.probability(junction_info.sides[1].seq_id, junction_test_info.pos_hash_score, junction_test_info.max_pos_hash_score);
      
      // Take the *least* significantly below pos_hash cutoff
      // Why this way? Consider an element moving from a plasmid to a genome
      // we don't want to penalize it with requiring coverage typical of the plasmid.
      junction_test_info.neg_log10_pos_hash_p_value = min(neg_log10_p_value_1, neg_log10_p_value_2);
      
      if (junction_test_info.neg_log10_pos_hash_p_value > settings.junction_pos_hash_neg_log10_p_value_cutoff) {
        junction_test_info.reject_reasons.push_back("COVERAGE_EVENNESS_SKEW");
      }
    }
        
    if (verbose) 
    {
      cout << "Testing Junction: " << junction_id << endl;
      cout << "  " << (junction_test_info.reject_reasons.size() ? "FAILED" : "SUCCESS") << endl;
      cout << "  Pos hash score: " << junction_test_info.pos_hash_score << endl;
      cout << "  Neg log10 pos hash p-value: " << junction_test_info.neg_log10_pos_hash_p_value
           << " [ " << settings.junction_pos_hash_neg_log10_p_value_cutoff << " ]" << endl;
      cout << "  Number of unique matches: " << unique_junction_match_map[junction_id].size() << endl;
      size_t num_degenerate_matches = repeat_junction_match_map.count(junction_id) ? repeat_junction_match_map[junction_id].size() : 0;
      cout << "  Number of degenerate matches: " << num_degenerate_matches << endl;
      cout << "  Number of total_non_overlap reads: " << junction_test_info.total_non_overlap_reads  << endl;
    }
    
    // Resolve junction no matter what. DO NOT break out of this loop before getting here.
    resolve_junction(
                     settings,
                     summary,
                     ref_seq_info,
                     junction_ref_seq_info,
                     trims_list,
                     junction_trims_list,
                     junction_test_info.junction_id,
                     unique_junction_match_map,
                     repeat_junction_match_map,
                     resolved_reference_tam,
                     resolved_junction_tam,
                     junction_test_info.reject_reasons.size(),
                     junction_test_info.total_non_overlap_reads > 0
                     ); 
    
    // However, we record only some junctions...
    if (record) {
      if (!junction_test_info.reject_reasons.size())
        passed_junction_test_info_list.push_back(junction_test_info);
      else
        rejected_junction_test_info_list.push_back(junction_test_info);
    }
    
    junction_test_info_list.pop_back();
    
    // @JEB 2019-06-12 List is now sorted after each one accepted to ensure best results
    junction_test_info_list.sort();
  }
    
  PosHashScoreDistribution accepted_pos_hash_score_distribution;
  for(list<JunctionTestInfo>::iterator it = passed_junction_test_info_list.begin(); it != passed_junction_test_info_list.end(); it++)
  {
    JunctionTestInfo& junction_test_info = *it;
		string key = junction_test_info.junction_id;
		if (verbose) cout << key << endl;
		cDiffEntry item = junction_to_diff_entry(key, ref_seq_info, junction_test_info);
		gd.add(item);

		// save the score in the distribution
		accepted_pos_hash_score_distribution.add_score(junction_test_info.pos_hash_score);

		// Create matches from UNIQUE sides of each match to reference genome
		// this fixes, for example appearing to not have any coverage at the origin of a circular DNA fragment
		// We do NOT add coverage to REDUNDANT sides because we don't know which copy.
		if (!settings.add_split_junction_sides) continue;

		for (uint32_t j = 0; j < unique_junction_match_map[key].size(); j++)
		{
			JunctionMatch& match = *unique_junction_match_map[key][j];
			bam_alignment& a = *(match.junction_alignments.front().get());
			uint32_t fastq_file_index = match.fastq_file_index;
      
      // at this point, all degeneracy should have been removed!
      assert(match.junction_alignments.size() == 1);

			for (uint32_t side = 1; side <= 2; side++)
			{
				string side_key = "side_" + to_string(side);

				// Do not count for coverage or write it out if it is redundant!!
				if (from_string<int32_t>(item[side_key + "_redundant"])) continue;
        
				// Write out match corresponding to this part to SAM file. write_moved_alignment()
				// computes the split read's XL/XR trims from the real-genome trims_list at its final
				// (post soft-clip) coordinates -- leaving the junction/middle side untrimmed.
				resolved_reference_tam.write_moved_alignment(
					a,
          resolved_junction_tam.target_name(a),
					fastq_file_index,
					item[side_key + "_seq_id"],
					from_string<int32_t>(item[side_key + "_position"]),
					from_string<int32_t>(item[side_key + "_strand"]),
					from_string<int32_t>(item[side_key + "_overlap"]),
					side,
					from_string<int32_t>(item["flanking_left"]),
					from_string<int32_t>(item["alignment_overlap"]),
          match.junction_alignments,
					&trims_list,
          &ref_seq_info,
          true
				);
			}
		}
	}

  for(list<JunctionTestInfo>::iterator it = rejected_junction_test_info_list.begin(); it != rejected_junction_test_info_list.end(); it++)
  {
    JunctionTestInfo& junction_test_info = *it;
		string key = junction_test_info.junction_id;
		cDiffEntry item = junction_to_diff_entry(key, ref_seq_info, junction_test_info);
    
    // Copy over the reject reasons
    for (vector<string>::iterator itr = junction_test_info.reject_reasons.begin(); itr != junction_test_info.reject_reasons.end(); itr++) {
      item.add_reject_reason(*itr);
    }
		gd.add(item);
	}

  // Save summary statistics
	summary.alignment_resolution.accepted_pos_hash_score_distribution = accepted_pos_hash_score_distribution;
  
  // Write the genome diff file
	gd.write(settings.jc_genome_diff_file_name);
  settings.track_intermediate_file(settings.output_done_file_name, settings.jc_genome_diff_file_name);

}
    
// Per-mate result of resolving one read's alignments against the reference and (optionally)
// candidate junctions -- lifted out of the old single-file loop body so it can be called once
// per mate when processing a paired read file set in lockstep.
struct MateResolution
{
  bool mapped_anywhere = false;
  alignment_list this_reference_alignments;
  alignment_list this_junction_alignments;
  uint32_t best_reference_score = 0;
  uint32_t best_junction_score = 0;
  int32_t mapping_quality_difference = 0; // best_junction_score - best_reference_score
};

static MateResolution resolve_one_mate(
                                       const Settings& settings,
                                       cReferenceSequences& ref_seq_info,
                                       cReferenceSequences& junction_ref_seq_info,
                                       const vector<ResolveJunctionInfo>& junction_info_list,
                                       bool junction_prediction,
                                       cFastqSequence& seq,
                                       alignment_list& reference_alignments, // peek buffer, mutated
                                       alignment_list& junction_alignments,  // peek buffer, mutated
                                       bam_file& reference_tam,
                                       bam_file* junction_tam,
                                       ReadFileSummary& read_file_summary_info,
                                       cFastqFile* unmapped_fastq
                                       )
{
  MateResolution m;

  read_file_summary_info.num_total_reads++;
  read_file_summary_info.num_total_bases += seq.length();

  // Does this read have eligible reference sequence matches?
  if ((reference_alignments.size() > 0) && (seq.m_name == reference_alignments.front()->read_name()))
  {
    m.this_reference_alignments = reference_alignments;
    reference_tam.read_alignments(reference_alignments, false);
    m.best_reference_score = eligible_read_alignments(settings, ref_seq_info, m.this_reference_alignments);
  }

  // Does this read have eligible candidate junction matches?
  if (junction_prediction && (junction_alignments.size() > 0) && (seq.m_name == junction_alignments.front()->read_name()))
  {
    m.this_junction_alignments = junction_alignments;
    junction_tam->read_alignments(junction_alignments, false);

    ///
    // Matches to candidate junctions MUST overlap the junction.
    //
    // Reduce this list to those that overlap ANY PART of the junction.
    // Alignments that extend only into the overlap region, are only additional
    //  evidence for predicted junctions and NOT support for a new junction on
    // their own. (They will also match the original reference genome equally well).
    // ... but this last point only if overlap >=0 for the junction
    ///

    for (alignment_list::iterator it = m.this_junction_alignments.begin(); it != m.this_junction_alignments.end(); )
    {
      if (!alignment_overlaps_junction(junction_info_list, it->get()))
        it = m.this_junction_alignments.erase(it);
      else
        it++;
    }

    m.best_junction_score = eligible_read_alignments(settings, junction_ref_seq_info, m.this_junction_alignments, settings.junction_allow_suboptimal_matches, m.best_reference_score);
  }

  // Nothing to be done if there were no eligible matches to either
  // Record in the unmatched FASTQ data file
  if ((m.this_junction_alignments.size() == 0) && (m.this_reference_alignments.size() == 0))
  {
    read_file_summary_info.num_unmapped_reads++;
    read_file_summary_info.num_unmapped_read_bases += seq.length();

    if (unmapped_fastq) {
      unmapped_fastq->write_sequence(seq);
    }
    m.mapped_anywhere = false;
  }
  else
  {
    m.mapped_anywhere = true;
  }

  // if < 0, then the best match is to the reference
  m.mapping_quality_difference = static_cast<int32_t>(m.best_junction_score) - static_cast<int32_t>(m.best_reference_score);

  return m;
}

// Dispatches one mate's already-resolved alignments exactly as the original single-file loop
// did: best match to reference is written immediately; best match to a candidate junction is
// deferred into the unique/repeat junction match maps for later adjudication. No-op if the mate
// had no alignments anywhere.
static void dispatch_mate_result(
                                 const Settings& settings,
                                 Summary& summary,
                                 cReferenceSequences& ref_seq_info,
                                 const SequenceTrimsList& trims_list,
                                 bam_file& resolved_reference_tam,
                                 bam_file* junction_tam,
                                 uint32_t fastq_file_index,
                                 const string& read_name,
                                 MateResolution& m,
                                 map<string,uint32_t>& all_junction_ids,
                                 UniqueJunctionMatchMap& unique_junction_match_map,
                                 RepeatJunctionMatchMap& repeat_junction_match_map
                                 )
{
  if (!m.mapped_anywhere) return;

  // best match is to the reference, record in that SAM file.
  if (m.mapping_quality_difference <= 0)
  {
    _write_reference_matches(settings, summary, ref_seq_info, trims_list, m.this_reference_alignments, resolved_reference_tam, fastq_file_index);
  }
  else
  {
    JunctionMatchPtr junction_match_ptr(
                                        new JunctionMatch(
                                                          m.this_reference_alignments,   // reference sequence alignments
                                                          m.this_junction_alignments,    // the BEST candidate junction alignments
                                                          fastq_file_index,              // index of the fastq file this read came from
                                                          m.mapping_quality_difference,  // difference between reference junction alignments (in # mismatches)
                                                          0                              //
                                                          )
                                        );

    ////
    // Just one best hit to candidate junctions, that is better than every match to the reference
    ////
    if ((m.this_junction_alignments.size() == 1) && (m.mapping_quality_difference > 0))
    {
      bam_alignment& a = *(m.this_junction_alignments.front().get());
      string junction_id = junction_tam->bam_header->target_name[a.reference_target_id()];
      unique_junction_match_map[junction_id].push_back( junction_match_ptr );
      all_junction_ids[junction_id]++;
    }
    ////
    // Multiple equivalent matches to junctions and reference, ones with most hits later will win these repeat matches
    // If mapping_quality_difference > 0, then they will count for scoring
    ////
    else
    {
      junction_match_ptr->degenerate_count = m.this_junction_alignments.size(); // mark as degenerate
      for(alignment_list::iterator it=m.this_junction_alignments.begin(); it!=m.this_junction_alignments.end(); it++)
      {
        bam_alignment& a = *(it->get());
        string junction_id = junction_tam->bam_header->target_name[a.reference_target_id()];
        repeat_junction_match_map[junction_id][read_name] = junction_match_ptr;
        all_junction_ids[junction_id]++;
      }
    }
  }
}

// Result of classifying a read pair's reference alignments against the expected majority
// orientation/distance cutoff -- generalizes best_pair_orientation_and_distance
// (candidate_junctions.cpp) to track ALL same-tid combinations (not just the smallest-distance
// one), so multi-mapped mates can be downselected to just the alignments that participate in at
// least one concordant combination with the other mate.
struct ConcordancePairing
{
  bool any_same_tid_combo_exists = false;
  bool any_concordant_combo_exists = false;
  bam_alignment* best_a = NULL;  // smallest-distance same-tid combo's mate-1 alignment
  bam_alignment* best_b = NULL;  // ...and mate-2 alignment
  string best_orientation;
  int64_t best_distance = -1;
  set<bam_alignment*> keep_mate1; // alignments participating in >=1 concordant combo
  set<bam_alignment*> keep_mate2;
};

// Coordinate-order orientation (F/R per mate, RR folded to FF) + outer-span distance for a
// specific pair of alignments -- shared by classify_pair's combination search and by the final
// BAM pair-info marking step (which must recompute fresh from the actual surviving alignments,
// not reuse classify_pair's "best" combo -- see mark_pair_info call site for why).
struct OrientationDistance { string orientation; int64_t distance; };

static OrientationDistance compute_orientation_and_distance(bam_alignment* a, bam_alignment* b)
{
  uint32_t a_start = a->reference_start_1();
  uint32_t a_end = a->reference_end_1();
  uint32_t b_start = b->reference_start_1();
  uint32_t b_end = b->reference_end_1();

  bool a_is_lower = (a_start <= b_start);
  bam_alignment* lower_alignment = a_is_lower ? a : b;
  bam_alignment* higher_alignment = a_is_lower ? b : a;

  string orientation;
  orientation += lower_alignment->reversed() ? 'R' : 'F';
  orientation += higher_alignment->reversed() ? 'R' : 'F';
  if (orientation == "RR") orientation = "FF";

  int64_t distance = static_cast<int64_t>(max(a_end, b_end)) - static_cast<int64_t>(min(a_start, b_start));

  return OrientationDistance{orientation, distance};
}

// Sets standard BAM pairing fields (BAM_FPAIRED, RNEXT/PNEXT/TLEN via mtid/mpos/isize,
// BAM_FMREVERSE) on both mates of an unambiguously-resolved pair, plus BAM_FPROPER_PAIR when
// breseq called the pair concordant (readable back later via that flag -- no separate custom
// concordance tag needed) and an XO:Z: orientation tag (no standard-field equivalent).
static void mark_pair_info(bam_alignment* a1, bam_alignment* a2, bool same_tid,
                            const string& orientation, int64_t distance, bool is_concordant)
{
  bam_alignment* mates[2] = {a1, a2};
  for (int i = 0; i < 2; i++)
  {
    bam_alignment* self = mates[i];
    bam_alignment* mate = mates[1 - i];
    bool self_is_leftmost = self->reference_start_1() <= mate->reference_start_1();

    self->core.flag |= BAM_FPAIRED;
    if (is_concordant) self->core.flag |= BAM_FPROPER_PAIR;
    if (mate->reversed()) self->core.flag |= BAM_FMREVERSE;
    self->core.mtid = mate->reference_target_id();
    self->core.mpos = mate->core.pos; // 0-based, matches core.mpos convention
    self->core.isize = same_tid ? (self_is_leftmost ? distance : -distance) : 0; // 0: undefined across references

    self->aux_set("XP", 'Z', orientation.size() + 1, (void*)orientation.c_str()); // 'Z' length includes null terminator; XO is already used by bowtie2 (gap opens)
  }
}

static ConcordancePairing classify_pair(
                                        alignment_list& mate1_alignments,
                                        alignment_list& mate2_alignments,
                                        const string& majority_orientation,
                                        double distance_cutoff
                                        )
{
  ConcordancePairing result;
  int64_t best_dist_seen = -1;

  for (alignment_list::iterator ita = mate1_alignments.begin(); ita != mate1_alignments.end(); ita++)
  {
    bam_alignment* pa = ita->get();
    for (alignment_list::iterator itb = mate2_alignments.begin(); itb != mate2_alignments.end(); itb++)
    {
      bam_alignment* pb = itb->get();
      if (pa->reference_target_id() != pb->reference_target_id()) continue;

      result.any_same_tid_combo_exists = true;

      OrientationDistance od = compute_orientation_and_distance(pa, pb);
      const string& orientation = od.orientation;
      int64_t distance = od.distance;

      if ((best_dist_seen < 0) || (distance < best_dist_seen))
      {
        best_dist_seen = distance;
        result.best_a = pa;
        result.best_b = pb;
        result.best_orientation = orientation;
        result.best_distance = distance;
      }

      if ((orientation == majority_orientation) && (distance <= distance_cutoff))
      {
        result.any_concordant_combo_exists = true;
        result.keep_mate1.insert(pa);
        result.keep_mate2.insert(pb);
      }
    }
  }

  return result;
}

// Erases every alignment from the list whose pointer isn't in keep -- same erase-in-list-loop
// idiom as eligible_read_alignments (candidate_junctions.cpp).
static void downselect_to_kept(alignment_list& alignments, const set<bam_alignment*>& keep)
{
  for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); )
  {
    if (keep.count(it->get()) == 0)
      it = alignments.erase(it);
    else
      it++;
  }
}

// Strand-aware 5' anchor position of an alignment: reference_start_1() if forward,
// reference_end_1() if reversed (reference_stranded_bounds_1 already computes exactly this).
static uint32_t stranded_anchor_position(bam_alignment* a)
{
  uint32_t start, end;
  a->reference_stranded_bounds_1(start, end);
  return start;
}

static void write_discordant_pair_row(
                                      ofstream& out,
                                      const string& read_number,
                                      const string& orientation,
                                      const string& seq_id,
                                      const string& start_1,
                                      const string& start_2,
                                      int64_t distance
                                      )
{
  out << read_number << "," << orientation << "," << seq_id << "," << start_1 << "," << start_2 << "," << distance << endl;
}

void load_junction_alignments(
                              const Settings& settings,
                              Summary& summary,
                              cReadFileSets& read_files,
                              cReferenceSequences& ref_seq_info,
                              cReferenceSequences& junction_ref_seq_info,
                              SequenceTrimsList& trims_list,
                              map<string,uint32_t>& all_junction_ids,
                              bool junction_prediction,
                              const vector<ResolveJunctionInfo>& junction_info_list,
                              UniqueJunctionMatchMap& unique_junction_match_map,
                              RepeatJunctionMatchMap& repeat_junction_match_map,
                              bam_file& resolved_reference_tam
                              )
{
  bool verbose = false;
  uint32_t reads_processed = 0;

  cFastqFile * unmapped_fastq = NULL;
  if (settings.output_unmapped_reads) {
    string unmapped_read_file_name = settings.unmapped_reads_fastq_file_name;
    unmapped_fastq = new cFastqFile(unmapped_read_file_name, ios::out);
  }

  cFastqQualityConverter fqc("SANGER", "SANGER");
  uint32_t fastq_file_index = 0;

  for (cReadFileSets::iterator rfs_it = read_files.begin(); rfs_it != read_files.end(); rfs_it++)
  {
    const cReadFileSet& rfs = *rfs_it;

    if (!rfs.is_paired())
    {
      ///
      //  UNPAIRED: the exact original single-file logic, now scoped to one read file set.
      ///
      bam_file* reference_tam = NULL;
      bam_file* junction_tam = NULL;

      const cReadFile& rf = rfs.m_files[0];
      string fastq_file_name = read_files.base_name_to_read_file_name(rf.m_base_name);

      end_progress_line();
      cerr << "  READ FILE:" << rf.m_base_name << endl;

      ReadFileSummary read_file_summary_info;

      // Traverse the original fastq files to keep track of order
      // b/c some matches may exist in only one or the other file

      cFastqFile in_fastq(fastq_file_name, ios::in);

      string reference_sam_file_name = settings.file_name(settings.reference_sam_file_name, "#", rf.m_base_name);
      string reference_fasta = settings.reference_fasta_file_name;

      reference_tam = new bam_file(reference_sam_file_name, settings.reference_fasta_file_name, ios::in);

      if (junction_prediction)
      {
        string junction_sam_file_name = settings.file_name(settings.candidate_junction_sam_file_name, "#", rf.m_base_name);
        junction_tam = new bam_file(junction_sam_file_name, settings.candidate_junction_fasta_file_name, ios::in);
      }

      alignment_list junction_alignments;

      //proceed through all of the alignments
      if (junction_prediction)
        junction_tam->read_alignments(junction_alignments, false);

      alignment_list reference_alignments;
      reference_tam->read_alignments(reference_alignments, false);

      ///
      //  Test each read for its matches to the reference and candidate junctions
      ///

      cFastqSequence seq;
      while (in_fastq.read_sequence(seq, fqc)) // READ
      {
        if ((settings.resolve_alignment_read_limit) && (reads_processed >= settings.resolve_alignment_read_limit))
          break; // to next file

        reads_processed++;
        read_file_summary_info.num_total_reads++;
        read_file_summary_info.num_total_bases+=seq.length();

        if (reads_processed % 10000 == 0) {
          ostringstream progress_message;
          progress_message << "    READS:" << setw(12) << right << reads_processed;
          print_progress_line(progress_message.str());
        }

        if (verbose)
          cerr << "===> Read: " << seq.m_name << endl;

        uint32_t best_junction_score = 0;
        uint32_t best_reference_score = 0;

        // Does this read have eligible reference sequence matches?
        alignment_list this_reference_alignments;
        if ((reference_alignments.size() > 0) && (seq.m_name == reference_alignments.front()->read_name()))
        {
          this_reference_alignments = reference_alignments;
          reference_tam->read_alignments(reference_alignments, false);

          if (verbose) {
            cerr << " Before Overlap Reference alignments = " << this_reference_alignments.size() << endl;
          }
          best_reference_score = eligible_read_alignments(settings, ref_seq_info, this_reference_alignments);
        }

        // Does this read have eligible candidate junction matches?
        alignment_list this_junction_alignments;

        if ((junction_alignments.size() > 0) && (seq.m_name == junction_alignments.front()->read_name()))
        {

          this_junction_alignments = junction_alignments;
          junction_tam->read_alignments(junction_alignments, false);

          if (verbose) {
            cerr << " Before Overlap Junction alignments = " << this_junction_alignments.size() << endl;
          }

          ///
          // Matches to candidate junctions MUST overlap the junction.
          //
          // Reduce this list to those that overlap ANY PART of the junction.
          // Alignments that extend only into the overlap region, are only additional
          //  evidence for predicted junctions and NOT support for a new junction on
          // their own. (They will also match the original reference genome equally well).
          // ... but this last point only if overlap >=0 for the junction
          ///

          for (alignment_list::iterator it = this_junction_alignments.begin(); it != this_junction_alignments.end(); )
          {
            if (!alignment_overlaps_junction(junction_info_list, it->get()))
              it = this_junction_alignments.erase(it);
            else
              it++;
          }

          best_junction_score = eligible_read_alignments(settings, junction_ref_seq_info, this_junction_alignments, settings.junction_allow_suboptimal_matches, best_reference_score);
        }

        // Nothing to be done if there were no eligible matches to either
        // Record in the unmatched FASTQ data file
        if ((this_junction_alignments.size() == 0) && (this_reference_alignments.size() == 0))
        {
          read_file_summary_info.num_unmapped_reads++;
          read_file_summary_info.num_unmapped_read_bases+=seq.length();

          if (unmapped_fastq) {
            unmapped_fastq->write_sequence(seq);
          }
        }

        ///
        // Determine if the read has a better match to a candidate junction
        // or to the reference sequence.
        ///

        /// There are three possible kinds of reads at this point
        //
        // 1: Read has a best match to the reference genome
        // --> Write this match and we are done
        // 2: Read has a best match (or multiple best matches) to junctions
        // --> Keep an item that describes these matches
        // 3: Read has an equivalent match to the reference genome
        //      and goes into the overlap part of a junction condidate
        // --> Keep an item that is not used during scoring
        ///

        // if < 0, then the best match is to the reference
        int32_t mapping_quality_difference = best_junction_score - best_reference_score;

        if (verbose)
        {
          cerr << " Best junction score: " << best_junction_score << endl;
          cerr << " Best reference score: " << best_reference_score << endl;
          cerr << " Mapping quality difference: " << mapping_quality_difference << endl;
          cerr << " Final Reference alignments = " << this_reference_alignments.size() << endl;
          cerr << " Final Candidate junction alignments = " << this_junction_alignments.size() << endl;
        }

        if ((this_junction_alignments.size() == 0) && (this_reference_alignments.size() == 0))
          continue;

        ///
        // The best match we found to the reference was no better than the best to the
        // candidate junction. This read potentially supports the candidate junction.
        //
        // ONLY allow EQUAL matches through if they match the overlap only, otherwise
        // you can get predictions of new junctions with all reads supporting them
        // actually mapping perfectly to the reference.
        ///

        // best match is to the reference, record in that SAM file.
        if (mapping_quality_difference <= 0)
        {
          if (verbose)
            cout << "Best alignment to reference. MQD: " << mapping_quality_difference << endl;

          _write_reference_matches(settings, summary, ref_seq_info, trims_list, this_reference_alignments, resolved_reference_tam, fastq_file_index);
        }
        else
        {
          if (verbose)
            cout << "Best alignment is to candidate junction. MQD: " << mapping_quality_difference << endl;

          JunctionMatchPtr junction_match_ptr(
                                              new JunctionMatch(
                                                                this_reference_alignments,    // reference sequence alignments
                                                                this_junction_alignments,     // the BEST candidate junction alignments
                                                                fastq_file_index,             // index of the fastq file this read came from
                                                                mapping_quality_difference,   // difference between reference junction alignments (in # mismatches)
                                                                0                             //
                                                                )
                                              );

          ////
          // Just one best hit to candidate junctions, that is better than every match to the reference
          ////
          if ((this_junction_alignments.size() == 1) && (mapping_quality_difference > 0))
          {
            bam_alignment& a = *(this_junction_alignments.front().get());
            string junction_id = junction_tam->bam_header->target_name[a.reference_target_id()];
            unique_junction_match_map[junction_id].push_back( junction_match_ptr );
            all_junction_ids[junction_id]++;
          }
          ////
          // Multiple equivalent matches to junctions and reference, ones with most hits later will win these repeat matches
          // If mapping_quality_difference > 0, then they will count for scoring
          ////
          else
          {
            if (verbose)
              cout << "this_junction_alignments: " << this_junction_alignments.size() << endl;

            junction_match_ptr->degenerate_count = this_junction_alignments.size(); // mark as degenerate
            for(alignment_list::iterator it=this_junction_alignments.begin(); it!=this_junction_alignments.end(); it++)
            {
              bam_alignment& a = *(it->get());
              string junction_id = junction_tam->bam_header->target_name[a.reference_target_id()];
              repeat_junction_match_map[junction_id][seq.m_name] = junction_match_ptr;
              all_junction_ids[junction_id]++;
            }
          }
        } // READ
      } // End loop through every $read_struct

      {
        ostringstream progress_message;
        progress_message << "    READS:" << setw(12) << right << reads_processed;
        print_progress_line(progress_message.str());
      }
      end_progress_line();

      // save statistics
      summary.alignment_resolution.read_file[rf.m_base_name] = read_file_summary_info;
      summary.alignment_resolution.total_unmapped_reads += read_file_summary_info.num_unmapped_reads;
      summary.alignment_resolution.total_unmapped_read_bases += read_file_summary_info.num_unmapped_read_bases;

      summary.alignment_resolution.total_reads += read_file_summary_info.num_total_reads;
      summary.alignment_resolution.total_bases += read_file_summary_info.num_total_bases;

      if (junction_tam != NULL) delete junction_tam;
      if (reference_tam != NULL) delete reference_tam;

      fastq_file_index++;
      continue;
    }

    ///
    //  PAIRED: stream both mates in lockstep so we can classify each read pair's orientation
    //  and distance against this read file set's majority orientation/distance cutoff.
    ///

    const cReadFile& rf1 = rfs.m_files[0];
    const cReadFile& rf2 = rfs.m_files[1];
    uint32_t fastq_file_index_1 = fastq_file_index;
    uint32_t fastq_file_index_2 = fastq_file_index + 1;
    fastq_file_index += 2;

    end_progress_line();
    cerr << "  READ FILE SET:" << rfs.m_base_name << endl;

    ReadFileSummary read_file_summary_info_1;
    ReadFileSummary read_file_summary_info_2;

    string fastq_file_name_1 = read_files.base_name_to_read_file_name(rf1.m_base_name);
    string fastq_file_name_2 = read_files.base_name_to_read_file_name(rf2.m_base_name);

    cFastqFile in_fastq_1(fastq_file_name_1, ios::in);
    cFastqFile in_fastq_2(fastq_file_name_2, ios::in);

    string reference_sam_file_name_1 = settings.file_name(settings.reference_sam_file_name, "#", rf1.m_base_name);
    string reference_sam_file_name_2 = settings.file_name(settings.reference_sam_file_name, "#", rf2.m_base_name);

    bam_file* reference_tam_1 = new bam_file(reference_sam_file_name_1, settings.reference_fasta_file_name, ios::in);
    bam_file* reference_tam_2 = new bam_file(reference_sam_file_name_2, settings.reference_fasta_file_name, ios::in);

    bam_file* junction_tam_1 = NULL;
    bam_file* junction_tam_2 = NULL;
    if (junction_prediction)
    {
      string junction_sam_file_name_1 = settings.file_name(settings.candidate_junction_sam_file_name, "#", rf1.m_base_name);
      string junction_sam_file_name_2 = settings.file_name(settings.candidate_junction_sam_file_name, "#", rf2.m_base_name);
      junction_tam_1 = new bam_file(junction_sam_file_name_1, settings.candidate_junction_fasta_file_name, ios::in);
      junction_tam_2 = new bam_file(junction_sam_file_name_2, settings.candidate_junction_fasta_file_name, ios::in);
    }

    alignment_list junction_alignments_1, junction_alignments_2;
    if (junction_prediction)
    {
      junction_tam_1->read_alignments(junction_alignments_1, false);
      junction_tam_2->read_alignments(junction_alignments_2, false);
    }

    alignment_list reference_alignments_1, reference_alignments_2;
    reference_tam_1->read_alignments(reference_alignments_1, false);
    reference_tam_2->read_alignments(reference_alignments_2, false);

    // Majority orientation / distance cutoff computed earlier in the pipeline (candidate-junction
    // preprocessing) for this read file set.
    const string majority_orientation = summary.paired_mapping_distance_distribution[rfs.m_base_name].majority_orientation;
    const double distance_cutoff = summary.paired_mapping_distance_distribution[rfs.m_base_name].distance_cutoff;

    string discordant_pairs_file_name = Settings::file_name(settings.discordant_pairs_file_name, "#", rfs.m_base_name);
    ofstream discordant_csv_out(discordant_pairs_file_name.c_str());
    ASSERT(discordant_csv_out.is_open(), "Could not write to file: " + discordant_pairs_file_name);
    discordant_csv_out << "read_number,orientation,seq_id,start_1,start_2,distance" << endl;

    cFastqSequence seq1, seq2;
    while (in_fastq_1.read_sequence(seq1, fqc) && in_fastq_2.read_sequence(seq2, fqc)) // READ PAIR
    {
      // Checked once per PAIR (not per mate) so a limit hit can't truncate one mate but not
      // the other.
      if ((settings.resolve_alignment_read_limit) && (reads_processed >= settings.resolve_alignment_read_limit))
        break; // to next read file set

      reads_processed += 2;

      if (reads_processed % 10000 == 0) {
        ostringstream progress_message;
        progress_message << "    READS:" << setw(12) << right << reads_processed;
        print_progress_line(progress_message.str());
      }

      MateResolution m1 = resolve_one_mate(settings, ref_seq_info, junction_ref_seq_info, junction_info_list, junction_prediction,
                                            seq1, reference_alignments_1, junction_alignments_1,
                                            *reference_tam_1, junction_tam_1, read_file_summary_info_1, unmapped_fastq);
      MateResolution m2 = resolve_one_mate(settings, ref_seq_info, junction_ref_seq_info, junction_info_list, junction_prediction,
                                            seq2, reference_alignments_2, junction_alignments_2,
                                            *reference_tam_2, junction_tam_2, read_file_summary_info_2, unmapped_fastq);

      // A mate counts as "a reference match" only if it actually has reference alignments to
      // use -- mapping_quality_difference <= 0 alone isn't sufficient, since both scores can be
      // tied at 0 (e.g. no reference alignments at all, but some junction alignments survive
      // eligible_read_alignments with their score clamped to 0), which would otherwise leave
      // this_reference_alignments empty despite passing the "reference is best" check.
      bool m1_is_reference_match = m1.mapped_anywhere && (m1.mapping_quality_difference <= 0) && (m1.this_reference_alignments.size() > 0);
      bool m2_is_reference_match = m2.mapped_anywhere && (m2.mapping_quality_difference <= 0) && (m2.this_reference_alignments.size() > 0);

      bool both_reference_best = m1_is_reference_match && m2_is_reference_match;

      if (both_reference_best)
      {
        ConcordancePairing pairing = classify_pair(m1.this_reference_alignments, m2.this_reference_alignments,
                                                    majority_orientation, distance_cutoff);

        if (pairing.any_concordant_combo_exists)
        {
          // Concordant: downselect each mate's alignments to just the ones participating in a
          // concordant combination (a no-op if the mate was already unique/already concordant).
          downselect_to_kept(m1.this_reference_alignments, pairing.keep_mate1);
          downselect_to_kept(m2.this_reference_alignments, pairing.keep_mate2);
        }
        else if (pairing.any_same_tid_combo_exists)
        {
          // Discordant: closest same-tid combination still fails the orientation/cutoff test.
          // Both mates are still written to resolved_reference_sam_file_name below (with their
          // original, non-downselected alignment lists) -- only the CSV logging differs here.
          string seq_id = reference_tam_1->bam_header->target_name[pairing.best_a->reference_target_id()];
          write_discordant_pair_row(discordant_csv_out, seq1.m_name, pairing.best_orientation, seq_id,
                                     to_string(stranded_anchor_position(pairing.best_a)),
                                     to_string(stranded_anchor_position(pairing.best_b)),
                                     pairing.best_distance);
        }
        else
        {
          // Discordant: mates map to different reference sequences entirely -- no orientation
          // or distance is meaningful, use each mate's own primary alignment for diagnostics.
          // Both mates are still written to resolved_reference_sam_file_name below.
          bam_alignment* a1 = m1.this_reference_alignments.front().get();
          bam_alignment* a2 = m2.this_reference_alignments.front().get();
          string seq_id = reference_tam_1->bam_header->target_name[a1->reference_target_id()];
          write_discordant_pair_row(discordant_csv_out, seq1.m_name, "NA", seq_id,
                                     to_string(stranded_anchor_position(a1)),
                                     to_string(stranded_anchor_position(a2)),
                                     -1);
        }

        // Whenever both mates end up with exactly one alignment (whether concordant or
        // discordant), mark them as paired in the BAM using standard fields. Recompute
        // orientation/distance fresh from these specific surviving alignments rather than
        // reusing pairing.best_a/best_b, which can point at a different alignment than the one
        // that survived downselection.
        if ((m1.this_reference_alignments.size() == 1) && (m2.this_reference_alignments.size() == 1))
        {
          bam_alignment* a1 = m1.this_reference_alignments.front().get();
          bam_alignment* a2 = m2.this_reference_alignments.front().get();
          bool same_tid = (a1->reference_target_id() == a2->reference_target_id());
          OrientationDistance od = same_tid ? compute_orientation_and_distance(a1, a2) : OrientationDistance{"NA", 0};
          mark_pair_info(a1, a2, same_tid, od.orientation, od.distance, pairing.any_concordant_combo_exists);
        }

        _write_reference_matches(settings, summary, ref_seq_info, trims_list, m1.this_reference_alignments, resolved_reference_tam, fastq_file_index_1);
        _write_reference_matches(settings, summary, ref_seq_info, trims_list, m2.this_reference_alignments, resolved_reference_tam, fastq_file_index_2);
      }
      else
      {
        // Either mate is junction-best, or one/both mates are fully unmapped -- these fall
        // through to the existing, unmodified per-mate handling. A singleton mapping (one mate
        // reference-best, the other fully unmapped) is still logged as discordant for
        // visibility, but the mapped mate's alignment is still written normally.
        bool m1_singleton_reference = m1_is_reference_match && !m2.mapped_anywhere;
        bool m2_singleton_reference = m2_is_reference_match && !m1.mapped_anywhere;

        if (m1_singleton_reference || m2_singleton_reference)
        {
          bam_alignment* mapped_alignment = m1_singleton_reference
            ? m1.this_reference_alignments.front().get()
            : m2.this_reference_alignments.front().get();
          string seq_id = reference_tam_1->bam_header->target_name[mapped_alignment->reference_target_id()];
          string start_1 = m1_singleton_reference ? to_string(stranded_anchor_position(mapped_alignment)) : "";
          string start_2 = m2_singleton_reference ? to_string(stranded_anchor_position(mapped_alignment)) : "";
          write_discordant_pair_row(discordant_csv_out, seq1.m_name, "NA", seq_id, start_1, start_2, -1);
        }

        dispatch_mate_result(settings, summary, ref_seq_info, trims_list, resolved_reference_tam, junction_tam_1, fastq_file_index_1, seq1.m_name, m1, all_junction_ids, unique_junction_match_map, repeat_junction_match_map);
        dispatch_mate_result(settings, summary, ref_seq_info, trims_list, resolved_reference_tam, junction_tam_2, fastq_file_index_2, seq2.m_name, m2, all_junction_ids, unique_junction_match_map, repeat_junction_match_map);
      }
    } // End loop through every read pair

    {
      ostringstream progress_message;
      progress_message << "    READS:" << setw(12) << right << reads_processed;
      print_progress_line(progress_message.str());
    }
    end_progress_line();

    // save statistics
    summary.alignment_resolution.read_file[rf1.m_base_name] = read_file_summary_info_1;
    summary.alignment_resolution.read_file[rf2.m_base_name] = read_file_summary_info_2;

    summary.alignment_resolution.total_unmapped_reads += read_file_summary_info_1.num_unmapped_reads + read_file_summary_info_2.num_unmapped_reads;
    summary.alignment_resolution.total_unmapped_read_bases += read_file_summary_info_1.num_unmapped_read_bases + read_file_summary_info_2.num_unmapped_read_bases;

    summary.alignment_resolution.total_reads += read_file_summary_info_1.num_total_reads + read_file_summary_info_2.num_total_reads;
    summary.alignment_resolution.total_bases += read_file_summary_info_1.num_total_bases + read_file_summary_info_2.num_total_bases;

    discordant_csv_out.close();

    if (junction_tam_1 != NULL) delete junction_tam_1;
    if (junction_tam_2 != NULL) delete junction_tam_2;
    delete reference_tam_1;
    delete reference_tam_2;

  } // End of Read File Set loop

  if (unmapped_fastq != NULL) delete unmapped_fastq;
}
  
  
void load_sam_only_alignments(
                              const Settings& settings, 
                              Summary& summary, 
                              cReadFileSets& read_files, 
                              cReferenceSequences& ref_seq_info,
                              SequenceTrimsList& trims_list,
                              bam_file& resolved_reference_tam
                              )
{
  
  uint32_t reads_processed = 0;
  summary.alignment_resolution.max_sam_base_quality_score = 0;

  bam_file* reference_tam = NULL;

  // One gzipped unmatched read file produced
  cFastqFile * out_unmapped_fastq = NULL;
  if (settings.output_unmapped_reads) {
    string unmapped_read_file_name = settings.unmapped_reads_fastq_file_name;
    out_unmapped_fastq = new cFastqFile(unmapped_read_file_name, ios::out);
  }

  for (cReadFileSets::iterator rfs_it = read_files.begin(); rfs_it != read_files.end(); rfs_it++) {
    if (rfs_it->is_paired()) {
      cerr << "  NOTE: Pairing-aware alignment resolution (discordant-pair detection and downselection)" << endl;
      cerr << "        is not applied in --aligned-sam mode. Read file set: " << rfs_it->m_base_name << endl;
    }
  }

  vector<cReadFile> flat_sam_read_files = read_files.flat_files();
  for (uint32_t sam_file_index = 0; sam_file_index < flat_sam_read_files.size(); sam_file_index++)
  {

    ReadFileSummary read_file_summary_info;

    const cReadFile& rf = flat_sam_read_files[sam_file_index];
    string reference_sam_file_name = read_files.base_name_to_read_file_name(rf.m_base_name);

    end_progress_line();
    cerr << "  READ FILE:" << rf.m_base_name << endl;
        
    reference_tam = new bam_file(reference_sam_file_name, settings.reference_fasta_file_name, ios::in); 
    
    
    ///
    //  Test each read for its matches to the reference and candidate junctions
    ///
        
    alignment_list this_reference_alignments;
    while (reference_tam->read_alignments(this_reference_alignments, false)) {
            
      // Grab this value before eligible alignments may remove all alignments
      uint32_t this_read_length = this_reference_alignments.front().get()->read_length();
      
      read_file_summary_info.num_total_reads++;
      read_file_summary_info.num_total_bases+=this_read_length;
      
      if ((settings.resolve_alignment_read_limit) && (reads_processed >= settings.resolve_alignment_read_limit))
        break; // to next file
      
      reads_processed++;
      if (reads_processed % 10000 == 0) {
        ostringstream progress_message;
        progress_message << "    READS:" << setw(12) << right << reads_processed;
        print_progress_line(progress_message.str());
      }
      
      // Does this read have eligible reference sequence matches? (junction matches are not possible)
      uint32_t best_reference_score = eligible_read_alignments(settings, ref_seq_info, this_reference_alignments);
      
      // if < 0, then the best match is to the reference
      // if == 0, then it is unmapped
      if ( (this_reference_alignments.size() == 0) || (best_reference_score == 0) ) {
 
          read_file_summary_info.num_unmapped_reads++;
          read_file_summary_info.num_unmapped_read_bases+=this_read_length;
          
        if (out_unmapped_fastq) {
          cFastqSequence seq;
          seq.m_sequence = this_reference_alignments.front().get()->read_char_sequence();
          seq.m_qualities = this_reference_alignments.front().get()->read_base_quality_char_string();
          out_unmapped_fastq->write_sequence(seq);
        }
        
        continue;
      }
      
      uint8_t* base_qualities = this_reference_alignments.front()->read_base_quality_bam_sequence();
      for(uint32_t base_index=0; base_index<this_reference_alignments.front()->read_length(); base_index++) {
        summary.alignment_resolution.max_sam_base_quality_score = max(summary.alignment_resolution.max_sam_base_quality_score, static_cast<int32_t>(base_qualities[base_index]));
      }
      
      // best match is to the reference, record in that SAM file.
      _write_reference_matches(settings, summary, ref_seq_info, trims_list, this_reference_alignments, resolved_reference_tam, sam_file_index);
      
    } // End loop through every read in file

    {
      ostringstream progress_message;
      progress_message << "    READS:" << setw(12) << right << reads_processed;
      print_progress_line(progress_message.str());
    }
    end_progress_line();

    summary.alignment_resolution.read_file[flat_sam_read_files[sam_file_index].m_base_name] = read_file_summary_info;
    
    summary.alignment_resolution.total_unmapped_reads += read_file_summary_info.num_unmapped_reads;
    summary.alignment_resolution.total_unmapped_read_bases += read_file_summary_info.num_unmapped_read_bases;
    
    summary.alignment_resolution.total_reads += read_file_summary_info.num_total_reads;
    summary.alignment_resolution.total_bases += read_file_summary_info.num_total_bases;
    
    // safe only because we know they are always or never used
    if (reference_tam != NULL) delete reference_tam;
    
  } // End of Read File loop
}

/*! Tests whether a read alignment to a candidate junction extends across
 *  the entire junction (past overlapping or unique sequence) and thus
 *  can be used as real evidence for the junction.
 */
bool alignment_overlaps_junction(const vector<ResolveJunctionInfo>& junction_info_list, const alignment_wrapper& a)
{
  // unmapped reads don't overlap the junction
  if (a.unmapped()) return false;
  
  int32_t tid = a.reference_target_id();
  const JunctionInfo& this_junction_info = junction_info_list[a.reference_target_id()];
  
  // DESCRIPTION OF LOGIC
  //
  // for overlap == 0
  //   junction is coords [flanking_left + 1, flanking_left - 0]
  // for overlap < 0 (unique sequence in read)
  //   junction is coords [flanking_left + 1, flanking_left - overlap]
  // for overlap > 0 (overlapping sequence matching to each side)
  //   junction is coords [flanking_left + 1, flanking_left + overlap]
  
  uint32_t junction_start, junction_end;
  
  junction_start = this_junction_info.flanking_left + 1;
  junction_end = this_junction_info.flanking_left + abs(this_junction_info.alignment_overlap);

  //## If it didn't overlap the junction at all
  //## Check coordinates in the "reference" (the JUNCTION sequence)
  if (a.reference_start_1() > junction_end) return false;
  if (a.reference_end_1() < junction_start) return false;
  
  return true;
}


void _write_reference_matches(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, const SequenceTrimsList& trims_list, alignment_list& reference_alignments, bam_file& reference_tam, uint32_t fastq_file_index)
{
  (void)settings;
	// Nice try, no alignments
	if (reference_alignments.size() == 0) return;

  double redundancy_corrected_count = 1.0 / static_cast<double>(reference_alignments.size());
  for(alignment_list::iterator it=reference_alignments.begin(); it!=reference_alignments.end(); it++)
  {
    summary.alignment_resolution.reference[(*it)->reference_target_id()].reads_mapped_to_reference  +=redundancy_corrected_count;
    summary.alignment_resolution.reference[(*it)->reference_target_id()].bases_mapped_to_reference  +=redundancy_corrected_count * (*it)->query_match_length();
  }
  summary.alignment_resolution.total_reads_mapped_to_references+=1;
  summary.alignment_resolution.total_bases_mapped_to_references+=reference_alignments.front()->query_match_length();

  // write_alignments() computes the XL/XR trims itself from trims_list (recomputing them for
  // any read whose ends it soft-clips), so we no longer precompute a Trims vector here.
	reference_tam.write_alignments((int32_t)fastq_file_index, reference_alignments, &trims_list, &ref_seq_info, true, true);
}
  
/*! Calculates various statistics about reads overlapping a junction
 */
void score_junction(
                    const Settings& settings, 
                    Summary& summary, 
                    const string& junction_id, 
                    UniqueJunctionMatchMap& unique_junction_match_map, 
                    RepeatJunctionMatchMap& repeat_junction_match_map, 
                    bam_file& resolved_junction_tam, 
                    JunctionTestInfo& junction_test_info, 
                    vector<ResolveJunctionInfo>& junction_info_list,
                    cReferenceSequences& ref_seq_info,
                    cReferenceSequences& junction_ref_seq_info
                  )
{
  bool verbose = false;
  
	if (verbose) cout << "Testing " << junction_id << endl;
  
	// There are two kinds of matches to a candidate junction:
  
	// (1) Reads that uniquely map to one candidate junction (but any number of times to reference)
	vector<JunctionMatchPtr>* unique_matches = NULL;
	if (unique_junction_match_map.count(junction_id)) {
		unique_matches = &(unique_junction_match_map[junction_id]);
  }
  
	// (2) Reads that uniquely map equally well to more than one candidate junction (and any number of times to reference)
	map<string, JunctionMatchPtr>* repeat_matches = NULL;
	if (repeat_junction_match_map.count(junction_id))
		repeat_matches = &(repeat_junction_match_map[junction_id]);
  
	// FAI target id -- there is no easy way to get this short of loading the entire array and going through them...
	// Debatable about whether we save more string comparisons by doing this here or each time
  
	// @JEB: optimization: it would be more efficient to hash junctions by target id rather than alignment junction names!!
	uint32_t junction_tid;
	for (junction_tid = 0; junction_tid < static_cast<uint32_t>(resolved_junction_tam.bam_header->n_targets); junction_tid++)
		if (resolved_junction_tam.bam_header->target_name[junction_tid] == junction_id) break;
	ASSERT(junction_tid < static_cast<uint32_t>(resolved_junction_tam.bam_header->n_targets), "Junction target id out of range.");
  
  uint32_t unique_matches_size = (unique_matches) ? unique_matches->size() : 0;
  uint32_t repeat_matches_size = (repeat_matches) ? repeat_matches->size() : 0;
  
	if (verbose) {
		cout << "Testing Junction Candidate: " << junction_id << endl;
		cout << "Unique Matches: " << unique_matches_size << " Degenerate Matches: " << repeat_matches_size << endl;
	}
  
	//// TEST 1: Reads that go a certain number of bp into the nonoverlap sequence on each side of the junction on each strand
	map<bool,int32_t> max_left_per_strand = make_map<bool,int32_t>(true,0)(false,0);
	map<bool,int32_t> max_right_per_strand = make_map<bool,int32_t>(true,0)(false,0);
	map<bool,int32_t> max_min_left_per_strand = make_map<bool,int32_t>(true,0)(false,0);
	map<bool,int32_t> max_min_right_per_strand = make_map<bool,int32_t>(true,0)(false,0);
	map<bool,uint32_t> count_per_strand = make_map<bool,uint32_t>(true,0)(false,0);
	uint32_t total_non_overlap_reads = 0;
	map<int32_t,bool> pos_hash;
  vector<pair<uint32_t,uint32_t> > start_end_check_list;
  bool has_reads_with_both_different_start_and_end(false);
  uint32_t pos_hash_count(0);
  
	// basic information about the junction
	ResolveJunctionInfo scj(junction_id);
	int32_t alignment_overlap = scj.alignment_overlap;
	int32_t flanking_left = scj.flanking_left;
    
	// We also need to count degenerate matches b/c sometimes ambiguity unfairly penalizes real junctions...
	vector<JunctionMatchPtr> items;
  if (unique_matches)
    for (vector<JunctionMatchPtr>::iterator it = unique_matches->begin(); it != unique_matches->end(); it++)
      items.push_back(*it);
  if (repeat_matches)
    for (map<string, JunctionMatchPtr>::iterator it = repeat_matches->begin(); it != repeat_matches->end(); it++)
      items.push_back(it->second);
  
	for (uint32_t i = 0; i < items.size(); i++) // READ (loops over unique_matches, degenerate_matches)
	{
		JunctionMatchPtr& item = items[i];
    
    if (verbose) cout << "  " << item->junction_alignments.front()->read_name() << endl;

    //! Do not count reads that map the reference equally well toward the score.
		if (item->mapping_quality_difference == 0) {
      if (verbose) cout << "    X Degenerate" << endl;
      continue; 
    }
    
    // Determine which alignment we are working with.
    
		// If there were no degenerate matches, then we could just take the
		// one and only match in the 'junction_alignments' array.
		// As it is, we must be sure we are looking at the one that matches
		alignment_wrapper* a = NULL;
    for(alignment_list::iterator it=item->junction_alignments.begin(); it!=item->junction_alignments.end(); it++)
		{
			alignment_wrapper* candidate_a = &(**it);
			if (candidate_a->reference_target_id() == junction_tid) {
				a = candidate_a;
				break;
			}
		}
		assert(a != NULL);

    // Only count alignments tied for best
    uint32_t is_best;
    bool tag_found = a->aux_get_i(kBreseqBestAlignmentScoreBAMTag, is_best);
    ASSERT(tag_found, "Did not find tag " + string(kBreseqBestAlignmentScoreBAMTag) + " for alignment");
    if (!is_best) 
      continue;
    
    ///
    // CHECK to be sure that this read overlaps the junction.
    // this_left and this_right are how far it extends into each side (past any overlap)
    ///
    
    // The left side goes exactly up to the flanking length
    uint32_t reference_start_1 = a->reference_start_1(); //settings.base_quality_cutoff);
    if (reference_start_1 == UNDEFINED_UINT32) continue;
		int32_t this_left = flanking_left + 1 - reference_start_1;
    
		// The right side starts after moving past any overlap (negative or positive)
    uint32_t reference_end_1 = a->reference_end_1(); //settings.base_quality_cutoff);
    if (reference_end_1 == UNDEFINED_UINT32) continue;
		int32_t this_right = reference_end_1 - flanking_left - abs(alignment_overlap);
    
    // doesn't span. The default for required_both_unique_length_per_side is 1, implying any overlap here is fine.
    if ((this_right < settings.required_both_unique_length_per_side) || (this_left < settings.required_both_unique_length_per_side)) {
      if (verbose) cout << "    X Does not span junction" << endl;
     continue; 
    }
    
    ////
    // CHECK that alignment starts at the first base of the query
    // and covers a certain amount of the read before counting toward the pos_hash score
    ////
    
    if (a->query_stranded_start_1() != 1) {
      if (verbose) cout << "    X First read base does not match" << endl;
      continue; 
    }
    
    if (a->query_stranded_end_1() < settings.required_junction_read_end_min_coordinate(a->read_length())) {
      if (verbose) cout << "    X End of read does not match as far as required" << endl;
      continue;
    }
    
    ////
    // COUNT reads that overlap both sides toward statistics other than pos_hash
    ////
    
    total_non_overlap_reads++;
    
		bool rev_key = a->reversed();
		count_per_strand[rev_key]++;
    
    ////
    // COUNT reads that overlap both sides toward the pos_hash_score
    ////
    
    // Note that reference here is the junction's sequence, not the reference genome sequence!
    uint32_t stranded_reference_start, stranded_reference_end;
    a->reference_stranded_bounds_1(stranded_reference_start, stranded_reference_end);
    
    if (verbose)
			cout << "  " << item->junction_alignments.front()->read_name() << ' ' << static_cast<int32_t>(rev_key) << ' ' << stranded_reference_start << endl;
    
    
    if (!pos_hash.count(stranded_reference_start))
    {
      pos_hash[stranded_reference_start] = true;
      pos_hash_count++;
    }
    if (verbose) cout << "    Y pos_hash: " << stranded_reference_start << endl;
    
    if (!has_reads_with_both_different_start_and_end) {
      for(vector<pair<uint32_t,uint32_t> >::iterator it=start_end_check_list.begin(); it != start_end_check_list.end();it++) {
        if ( (it->first != reference_start_1) && (it->second != reference_end_1) ) {
          has_reads_with_both_different_start_and_end = true;
          break;
        }
      }
      
      if (!has_reads_with_both_different_start_and_end) {
        start_end_check_list.push_back(std::pair<uint32_t,uint32_t>(reference_start_1,reference_end_1));
      }
    }
      
    // Update:
    // Max_Min = the maximum of the minimum length match sides
    // Max = the maximum match on a side
    // Note that the max and min filtering is really a kind of poor man's KS test
    //   if we implemented that with a certain coverage cutoff it would be a
    //   more principled way of doing things...
    if (this_left < this_right) {
      if (max_min_left_per_strand[rev_key] < this_left)
        max_min_left_per_strand[rev_key] = this_left;
    }
    else
    {
      if (max_min_right_per_strand[rev_key] < this_right)
        max_min_right_per_strand[rev_key] = this_right;
    }
    
    max_left_per_strand[rev_key] = max(this_left, max_left_per_strand[rev_key]);
    max_right_per_strand[rev_key] = max(this_right, max_right_per_strand[rev_key]);
	}
  
	int32_t max_left = max(max_left_per_strand[false], max_left_per_strand[true]);
	int32_t max_right = max(max_right_per_strand[false], max_right_per_strand[true]);
  
	int32_t max_min_left = max(max_min_left_per_strand[false], max_min_left_per_strand[true]);
	int32_t max_min_right = max(max_min_right_per_strand[false], max_min_right_per_strand[true]);
  
  // UPDATE REDUNDANCY
  // =================
  // If all the matches were to repeats, then at least one of the sides of this junction
  // needs to be marked as newly degenerate. Even through the whole junction sequences were
  // unique, no reads extended far enough to disambiguate between them
  
  // @JEB: TODO broken? needs to construct a unique list of all junction_ids supported by all degenerate matches
  
  bool redundant[2] = {false, false};

  map<uint32_t,bool> repeat_junction_tid_map;

  if (!unique_matches && repeat_matches) {

    // READ (loops over repeats matches)
    for (uint32_t i = 0; i < items.size(); i++) 
    {
      JunctionMatchPtr& item = items[i];
      
      for(alignment_list::iterator it=item->junction_alignments.begin(); it!=item->junction_alignments.end(); it++)
      {
        repeat_junction_tid_map[(*it)->reference_target_id()] = true;
      }
    }
    
    vector<uint32_t> repeat_junction_tid_list = get_keys(repeat_junction_tid_map);
    
    vector<uint32_t>::iterator it = repeat_junction_tid_list.begin();
    ResolveJunctionInfo main = junction_info_list[*it];
    for (it++; it != repeat_junction_tid_list.end(); it++)
    {
      ResolveJunctionInfo test = junction_info_list[*it]; // the junction key
      
      for (uint32_t best_side_index=0; best_side_index<=1; best_side_index++)
      {
        bool found = false;
        for (uint32_t test_side_index=0; test_side_index<=1; test_side_index++)
        {
          if (main.sides[best_side_index] == test.sides[test_side_index])
            found=true;
        }
        
        // didn't find this side -- mark as redundant
        if (!found)
        {
          redundant[best_side_index] = true;
          //if (verbose) cout << "Marking side " << best_side_index << " as redundant." << endl;
        }
      }
    }
  }
  
  ////
  // Calculate the maximum possible pos_hash_score
  ////
  
  uint32_t side_1_continuation, side_2_continuation;
  calculate_continuation(scj, ref_seq_info, junction_ref_seq_info, side_1_continuation, side_2_continuation);
  
  uint32_t read_length_avg = static_cast<uint32_t>(round(summary.sequence_conversion.read_length_avg));
  
  // For  junctions the number of start sites where reads crossing the
  // new junction sequence could occur is reduced by overlap (positive or negative):
  int32_t overlap_positions_max_pos_hash_score_reduction = abs(scj.alignment_overlap);
  int32_t max_pos_hash_score = 2 * (read_length_avg - 1 - 2 * (settings.required_both_unique_length_per_side - 1) - overlap_positions_max_pos_hash_score_reduction - side_1_continuation - side_2_continuation);
  
  //@JEB Actually, this can happen if the read lengths vary, so we better not rule it out as an error!
  /*
  ASSERT(pos_hash_count <= max_pos_hash_score, 
         "Pos hash score (" + to_string(pos_hash_count) + ") is greater than calculated maximum (" + to_string(max_pos_hash_score) 
         + ").\nFor junction: " + junction_id );
  */
  
	// Save the test info about this junction.
	JunctionTestInfo this_junction_test_info = {
		max_left,                           //max_left
		max_left_per_strand[false],         //max_left_minus
		max_left_per_strand[true],          //max_left_plus
		max_right,                          //max_right
		max_right_per_strand[false],        //max_right_minus
		max_right_per_strand[true],         //max_right_plus
		max_min_right,                      //max_min_right
		max_min_right_per_strand[false],    //max_min_right_minus
		max_min_right_per_strand[true],     //max_min_right_plus
		max_min_left,                       //max_min_left
		max_min_left_per_strand[false],     //max_min_left_minus
		max_min_left_per_strand[true],      //max_min_left_plus
    count_per_strand[false],            //coverage_minus
    count_per_strand[true],             //coverage_plus
		total_non_overlap_reads,            //total_non_overlap_reads
    pos_hash_count,                     //pos_hash_score
    static_cast<uint32_t>(max(0,max_pos_hash_score)),          //max_pos_hash_score
    unique_matches_size,                //unique_matches_size
    repeat_matches_size,                //repeat_matches_size
    has_reads_with_both_different_start_and_end, //2 diff reads have both different start and different end
    redundant[0],
    redundant[1],
    junction_id,
    999999999.99,
    side_1_continuation,
    side_2_continuation,
    vector<string>()
	};
  
  junction_test_info = this_junction_test_info;
}
  
/*! deals with the reads corresponding to a successful or failed junction
 */
void resolve_junction(
                      const Settings& settings,
                      Summary& summary,
                      cReferenceSequences& ref_seq_info,
                      cReferenceSequences& junction_ref_seq_info,
                      const SequenceTrimsList& trims_list,
                      const SequenceTrimsList& junction_trims_list,
                      const string& junction_id,
                      UniqueJunctionMatchMap& unique_junction_match_map,
                      RepeatJunctionMatchMap& repeat_junction_match_map,
                      bam_file& resolved_reference_tam,
                      bam_file& resolved_junction_tam,
                      bool failed,
                      bool has_non_overlap_alignment
                      )
{
  (void) summary;
  bool verbose = false;
  
  // There are two kinds of matches to a candidate junction:
  
	// (1) Reads that uniquely map to one candidate junction (but any number of times to reference)
	vector<JunctionMatchPtr>* unique_matches = NULL;
	if (unique_junction_match_map.count(junction_id)) {
		unique_matches = &(unique_junction_match_map[junction_id]);
  }
  
	// (2) Reads that uniquely map equally well to more than one candidate junction (and any number of times to reference)
	map<string, JunctionMatchPtr>* repeat_matches = NULL;
	if (repeat_junction_match_map.count(junction_id))
		repeat_matches = &(repeat_junction_match_map[junction_id]);
  
  // @JEB: optimization: it would be more efficient to hash junctions by target id rather than alignment junction names!!
	uint32_t junction_tid;
	for (junction_tid = 0; junction_tid < static_cast<uint32_t>(resolved_junction_tam.bam_header->n_targets); junction_tid++)
		if (resolved_junction_tam.bam_header->target_name[junction_tid] == junction_id) break;
	ASSERT(junction_tid < static_cast<uint32_t>(resolved_junction_tam.bam_header->n_targets), "Junction target id out of range.");

  
  // DEGENERATE JUNCTION MATCHES
	// ===========================
	// Determine the fate of degenerate reads that map to this junction
  
	if (repeat_matches)
	{
    
    // Figure out whether each side is redundantly matched here...
    
		for (map<string, JunctionMatchPtr>::iterator it = repeat_matches->begin(); it != repeat_matches->end(); it++)
		{      
			JunctionMatchPtr& repeat_match_ptr = it->second;
      JunctionMatch& repeat_match = *repeat_match_ptr;
      
			uint32_t fastq_file_index = repeat_match.fastq_file_index;
      
			// Success for this candidate junction...
			// purge all references to this from the degenerate match hash
			// so that they will not be counted for other junctions
			if (!failed)
			{
        if (verbose) cout << "Unique matches before size: " << (unique_matches ? unique_matches->size() : 0) << endl;
        
				// Purge all references to this read from the degenerate match hash
        // so that it cannot be counted for any other junction
        
        counted_ptr<bam_alignment> matched_alignment(NULL);
        for (alignment_list::iterator it2=repeat_match.junction_alignments.begin(); it2 !=repeat_match.junction_alignments.end(); )
				{          
          // we make a copy and then increment, in case the current iterator value will be erased
					counted_ptr<bam_alignment> a = *it2; it2++;
          string test_junction_seq_id = resolved_junction_tam.target_name(*a);
          
          //this is the one for the current candidate junction
          if (a->reference_target_id() == junction_tid)
          {
            matched_alignment = a;
          }
          else
          {
            size_t deleted = repeat_junction_match_map[test_junction_seq_id].erase(a->read_name());
          }
          
          if (repeat_junction_match_map[test_junction_seq_id].size() == 0)
          {
            repeat_junction_match_map.erase(test_junction_seq_id);
          }
        }
        
				assert(matched_alignment.get() != NULL);
				repeat_match.junction_alignments.clear();
        repeat_match.junction_alignments.push_back(matched_alignment);
        
        // We need to add this degenerately matched read to the other ones supporting this junction
        // Create empty list if necessary...
        if (unique_junction_match_map.count(junction_id) == 0) {
          unique_junction_match_map.insert( pair<string, vector<JunctionMatchPtr> >(junction_id, vector<JunctionMatchPtr>()) );
          unique_matches = &(unique_junction_match_map[junction_id]);
        }
				unique_matches->push_back(repeat_match_ptr);
        
        if (verbose) cout << "Unique matches after size: " << (unique_matches ? unique_matches->size() : 0) << endl;
        
			}
      
			// Failure for this candidate junction...
			// Remove just the degenerate hits to this candidate junction
			// Once all have failed, then we need to add the reference alignments (if any)!
			else
			{
				repeat_match.degenerate_count--;
        
				if (verbose) cout << "New Degenerate match count: " << repeat_match.degenerate_count << endl;
        
				// This degenerate match missed on all opportunities,
				// we should add it to the reference sequence
				if (repeat_match.degenerate_count == 0)
				{
					alignment_list& this_reference_al = repeat_match.reference_alignments;
					_write_reference_matches(settings, summary, ref_seq_info, trims_list, this_reference_al, resolved_reference_tam, fastq_file_index);
				}
        
        counted_ptr<bam_alignment> matched_alignment(NULL);
        for (alignment_list::iterator it2=repeat_match.junction_alignments.begin(); it2 !=repeat_match.junction_alignments.end(); it2++)
				{
					counted_ptr<bam_alignment>& candidate_a = *it2; //this is the one for the current candidate junction
					if (candidate_a->reference_target_id() == junction_tid)
          {
						matched_alignment = candidate_a;
            break;
          }
        }
        
        // Write alignment to SAM file for candidate junctions regardless of success...
        // Note that successful ones get written below, because they were pushed to the other list
        assert(matched_alignment.get() != NULL);
        if (has_non_overlap_alignment) {
          alignment_list alignments;
          alignments.read_base_quality_char_string = repeat_match.junction_alignments.read_base_quality_char_string;
          alignments.read_base_quality_char_string_reversed = repeat_match.junction_alignments.read_base_quality_char_string_reversed;

          alignments.push_back(matched_alignment);
          resolved_junction_tam.write_alignments(fastq_file_index, alignments, &junction_trims_list, &junction_ref_seq_info, true);
        }
			}
		}
    
		// We are completely done with degenerate matches to this junction id.
		// Deleting them here means that we will never go through this loop with them again
		// and is necessary for not doubly writing them.
		repeat_junction_match_map.erase(junction_id);
	}
  
	// UNIQUE JUNCTION MATCHES
	// =======================
  // If there were no unique matches to begin with, we may have created this entry by processing degenerate junctions...
  
  
  if (unique_matches)
  {
    if (verbose) cout << "Printing size:" << (unique_matches ? unique_matches->size() : 0) << endl;
    
    for (uint32_t i = 0; i< unique_matches->size(); i++)
    {
      JunctionMatch& item = *((*unique_matches)[i]);
      // Write out the matches to the proper SAM file(s) depending on whether the junction succeeded or failed
      uint32_t fastq_file_index = item.fastq_file_index;
      
      // ONLY if we failed: write matches to reference sequences
      if (failed)
      {
        alignment_list this_reference_al = item.reference_alignments;
        _write_reference_matches(settings, summary, ref_seq_info, trims_list, this_reference_al, resolved_reference_tam, fastq_file_index);
      }
      
      // REGARDLESS of success: write matches to the candidate junction SAM file
      resolved_junction_tam.write_alignments(fastq_file_index, item.junction_alignments, &junction_trims_list, &junction_ref_seq_info, true);
    }
  }

}

cDiffEntry junction_to_diff_entry(
                                         const string& key, 
                                         cReferenceSequences& ref_seq_info,
                                         JunctionTestInfo& test_info
                                         )
{
	// split the key to an item with information about the junction
	ResolveJunctionInfo jc(key);

	jc.key = key;

	// overlap may be adjusted below... this messes up making the alignment
	// 'alignment_overlap' is the original one that applies to the candidate junction BAM file
	// 'overlap' is a version where overlap has been resolved if possible for adding sides of the
	//    alignment
	jc.overlap = jc.alignment_overlap;

	// Redundancy is loaded from the key, but we doubly enforce it when IS elements are involved.

	// Correct for overlapping IS elements

	///
	// IS insertion overlap correction
	//
	// For these the coordinates may have been offset incorrectly initially (because both sides of the junction may look unique)
	// The goal is to offset through positive overlap to get as close as possible to the ends of the IS
	///

	cFeatureLocation* repeat_ptr(NULL);
	for (int32_t i = 0; i <= 1; i++)
	{
		// Determine IS elements
		// Is it within an IS or near the boundary of an IS in the direction leading up to the junction?
    int32_t max_distance_to_repeat = 20;
		repeat_ptr = cReferenceSequences::find_closest_repeat_region_boundary(jc.sides[i].position, ref_seq_info[jc.sides[i].seq_id].m_repeats, max_distance_to_repeat, jc.sides[i].strand);
		if (repeat_ptr != NULL) {
			jc.sides[i].is = repeat_ptr;
			jc.sides[i].is_interval = (repeat_ptr->get_strand() == 1)
        ? to_string(repeat_ptr->get_start_1()) + "-" + to_string(repeat_ptr->get_end_1())
        : to_string(repeat_ptr->get_end_1()) + "-" + to_string(repeat_ptr->get_start_1());
		}
	}

	// Determine which side of the junction is the IS and which is unique
	// these point to the correct initial interval...
	jc.is_side = UNDEFINED_UINT32;
	if (jc.sides[0].is && !jc.sides[1].is)
	{
		if (abs(static_cast<int32_t>(jc.sides[0].is->get_start_1()) - static_cast<int32_t>(jc.sides[0].position)) <= 20)
		{
			jc.is_side = 0;
			jc.sides[jc.is_side].is_side_key = "start";
		}
		else if (abs(static_cast<int32_t>(jc.sides[0].is->get_end_1()) - static_cast<int32_t>(jc.sides[0].position)) <= 20 )
		{
			jc.is_side = 0;
			jc.sides[jc.is_side].is_side_key = "end";
		}
		jc.unique_side = 1;
	}

	else if (!jc.sides[0].is && jc.sides[1].is)
	{
		if (abs(static_cast<int32_t>(jc.sides[1].is->get_start_1()) - static_cast<int32_t>(jc.sides[1].position)) <= 20)
		{
			jc.is_side = 1;
			jc.sides[jc.is_side].is_side_key = "start";
		}
		else if (abs(static_cast<int32_t>(jc.sides[1].is->get_end_1()) - static_cast<int32_t>(jc.sides[1].position)) <= 20 )
		{
			jc.is_side = 1;
			jc.sides[jc.is_side].is_side_key = "end";
		}
		jc.unique_side = 0;
	}
  
	// By default, overlap is included on both sides of the junction (possibly changed below)
	jc.sides[0].overlap = 0;
	jc.sides[1].overlap = 0;

	// Resolve redundant overlap
	if (jc.overlap > 0)
	{
		jc.sides[0].overlap = jc.overlap;
		jc.sides[1].overlap = jc.overlap;

		// If there was in IS, resolve overlap so it goes to the edge of the IS element
		if (jc.is_side != UNDEFINED_UINT32)
		{
			// first, adjust the repetitive sequence boundary to get as close to the IS as possible
      assert(jc.sides[jc.is_side].is_side_key.size() > 0);
			int32_t move_dist = jc.sides[jc.is_side].strand * (static_cast<int32_t>((jc.sides[jc.is_side].is_side_key == "start" 
          ? jc.sides[jc.is_side].is->get_start_1() : jc.sides[jc.is_side].is->get_end_1())) - jc.sides[jc.is_side].position);

			if (move_dist < 0) move_dist = 0;
			if (move_dist > jc.overlap) move_dist = jc.overlap ;

			jc.sides[jc.is_side].position += jc.sides[jc.is_side].strand * move_dist;
			jc.overlap -= move_dist;
			jc.sides[jc.is_side].overlap -= move_dist;

			// second, adjust the unique sequence side with any remaining overlap
			jc.sides[jc.unique_side].position += jc.sides[jc.unique_side].strand * jc.overlap;
			jc.sides[jc.unique_side].overlap -= jc.overlap;

			jc.overlap = 0;
		}
		/// If there is no IS element and
		//    (1) both sides are unique
		// OR (2) only the second side is redundant,
		// OR (3) both sides are redundant
		/// then give overlap to first side.
		/// This gives proper support for junctions.
		/// and ensures we don't count this coverage twice.
		else if ((!jc.sides[0].redundant) || (jc.sides[0].redundant && jc.sides[1].redundant) )
		{
			int32_t strand_direction = (jc.sides[1].strand > 0 ? 1 : -1);
			jc.sides[1].position += jc.overlap * strand_direction;
			jc.sides[1].overlap = 0;
			jc.overlap = 0;
		}
		else  // side_1 was redundant, give overlap to side_2
		{
			int32_t strand_direction = (jc.sides[0].strand > 0 ? 1 : -1);
			jc.sides[0].position += jc.overlap * strand_direction;
			jc.sides[0].overlap = 0;
			jc.overlap = 0;
		}

	}
  
	// Flatten things to only what we want to keep
  cDiffEntry item(JC);
	item
		("side_1_seq_id", jc.sides[0].seq_id)
		("side_1_position", to_string(jc.sides[0].position))
		("side_1_redundant", to_string(jc.sides[0].redundant))
		("side_1_strand", to_string(jc.sides[0].strand))
		("side_1_overlap", to_string(jc.sides[0].overlap))

		("side_2_seq_id", jc.sides[1].seq_id)
		("side_2_position", to_string(jc.sides[1].position))
		("side_2_redundant", to_string(jc.sides[1].redundant))
		("side_2_strand", to_string(jc.sides[1].strand))
		("side_2_overlap", to_string(jc.sides[1].overlap))

		("key", jc.key)
		("alignment_overlap", to_string(jc.alignment_overlap))
		("overlap", to_string(jc.overlap))
		("flanking_left", to_string(jc.flanking_left))
		("flanking_right", to_string(jc.flanking_right))

    ("unique_read_sequence", to_string(jc.unique_read_sequence))
	;
  
  if (jc.user_defined) item[USER_DEFINED] = "1";
  
  //	## may want to take only selected of these fields in the future.
  
  item
  ("max_left", to_string(test_info.max_left))
  ("max_left_minus", to_string(test_info.max_left_minus))
  ("max_left_plus", to_string(test_info.max_left_plus))
  ("max_right", to_string(test_info.max_right))
  ("max_right_minus", to_string(test_info.max_right_minus))
  ("max_right_plus", to_string(test_info.max_right_plus))
  ("max_min_right", to_string(test_info.max_min_right))
  ("max_min_right_minus", to_string(test_info.max_min_right_minus))
  ("max_min_right_plus", to_string(test_info.max_min_right_plus))
  ("max_min_left", to_string(test_info.max_min_left))
  ("max_min_left_minus", to_string(test_info.max_min_left_minus))
  ("max_min_left_plus", to_string(test_info.max_min_left_plus))
  ("coverage_minus", to_string(test_info.coverage_minus))
  ("coverage_plus", to_string(test_info.coverage_plus))
  ("total_non_overlap_reads", to_string(test_info.total_non_overlap_reads))
  ("pos_hash_score", to_string(test_info.pos_hash_score))
  ("max_pos_hash_score", to_string(test_info.max_pos_hash_score))
  ("neg_log10_pos_hash_p_value", test_info.neg_log10_pos_hash_p_value == -1 ? "NT" : to_string(test_info.neg_log10_pos_hash_p_value, 1, false))
  ("side_1_continuation", to_string<uint32_t>(test_info.side_1_continuation))
  ("side_2_continuation", to_string<uint32_t>(test_info.side_2_continuation))
  ;

	/// Note: Other adjustments to overlap can happen at the later annotation stage
	/// and they will not affect coverage for calling deletions or mutations
	/// because they will be in REDUNDANTLY matched sides of junctions
	return item;
}
  

//sort junction ids based on size of vector contained in map
vector<string> get_sorted_junction_ids(
                                       UniqueJunctionMatchMap& unique_map, 
                                       RepeatJunctionMatchMap& degenerate_map, 
                                       const vector<string>& keys
                                       )
{
  bool verbose = false;
  
  vector<VectorSize> vector_sizes;
  for (uint32_t i = 0; i < keys.size(); i++)
  {
    const string& junction_id = keys[i];
    
    // may or may not exist
    uint32_t num_degenerate_matches = (degenerate_map.count(junction_id)) ? degenerate_map[junction_id].size() : 0;
    uint32_t num_unique_matches = unique_map[junction_id].size();
    
    if (verbose) {
      cout << "Pre-sort Junction: " << junction_id << endl;
      cout << "  Number of unique matches: " << num_unique_matches << endl;   
      
      for(vector<JunctionMatchPtr>::iterator it=unique_map[junction_id].begin(); it!= unique_map[junction_id].end(); it++) {
        JunctionMatch& m = **it;
        cout << "    " << m.junction_alignments.begin()->get()->read_name() << endl;
      }
      cout << "  Number of degenerate matches: " << num_degenerate_matches << endl;
      
      for(map<string, JunctionMatchPtr>::iterator it=degenerate_map[junction_id].begin(); it!= degenerate_map[junction_id].end(); it++) {
        JunctionMatch& m = *it->second;
        cout << "    " << m.junction_alignments.begin()->get()->read_name() << endl;
      }

    }
    
    VectorSize info(keys[i], num_unique_matches+num_degenerate_matches, num_degenerate_matches);
    vector_sizes.push_back(info);
  }
  sort(vector_sizes.begin(), vector_sizes.end(), VectorSize::sort_reverse_by_size);
  
  if (verbose) {
    cout << "SORTED--->" << endl;
  }
  
  vector<string> sorted_junction_ids;
  for (uint32_t i = 0; i < vector_sizes.size(); i++) {
    sorted_junction_ids.push_back(vector_sizes[i].junction_id);
    
    cout << "Post-sort Junction: " << vector_sizes[i].junction_id << endl;
    cout << "  Number of unique matches: " << vector_sizes[i].size << endl;   
    cout << "  Number of degenerate matches: " << vector_sizes[i].size2 << endl;
  }
  return sorted_junction_ids;
}

void  assign_one_junction_read_counts(
                                  const Settings& settings,
                                  Summary& summary,
                                  cDiffEntry& j,
                                  const counted_ptr<junction_read_counter>& reference_jrc,
                                  const counted_ptr<junction_read_counter>& junction_jrc
                                  )
{
  
  bool verbose = false;
  bool debug_output = settings.junction_debug;
  
  uint32_t read_length_avg = summary.sequence_conversion.read_length_avg;
 
  fstream ofile;
  if (settings.junction_debug) {
    ofile.open(settings.junction_debug_file_name.c_str(), ios_base::out | ios_base::app);
    ASSERT(ofile.good(), "Could not open file " + settings.junction_debug_file_name);
    ofile << j << endl;
  }
 
  map<string,bool> empty_read_names;
  map<string,bool> junction_read_names;
  
  int32_t start, end;
  
  uint32_t side_1_continuation = from_string<uint32_t>(j["side_1_continuation"]);
  uint32_t side_2_continuation = from_string<uint32_t>(j["side_2_continuation"]);
  
  // Must be positive
  int32_t extra_stranded_require_overlap = 0;
  if (j.entry_exists("read_count_offset")) {
    extra_stranded_require_overlap = max(from_string<int32_t>(j["read_count_offset"]), 0);
  }
  
  // This is for requiring a certain number of bases (at least one) to match past the normal point
  // where a read could be uniquely assigned to the junction (or a side)
  int32_t minimum_side_match_correction = settings.junction_minimum_side_match - 1;

  
  // Print out the junction we are processing
  if (verbose) cerr << endl << "ASSIGNING READ COUNTS TO JUNCTION" << endl << j << endl << endl; 
  
  
  if (verbose) {
    cerr << "==SIDE 1==" << endl;
    cerr << "  position:" << j[SIDE_1_POSITION] << endl;
    cerr << "  strand:" << j[SIDE_1_STRAND] << endl;
    cerr << "  flanking:" << j["flanking_left"] << endl;
    cerr << "  continuation:" << j["side_1_continuation"] << endl;
    cerr << "  overlap:" << j[SIDE_1_OVERLAP] << endl;

    cerr << "==OVERLAP== " << j[ALIGNMENT_OVERLAP] << endl;
    
    cerr << "==SIDE 2==" << endl;
    cerr << "  position:" << j[SIDE_2_POSITION] << endl;
    cerr << "  strand:" << j[SIDE_2_STRAND] << endl;
    cerr << "  flanking:" << j["flanking_right"] << endl;
    cerr << "  continuation:" << j["side_2_continuation"] << endl;
    cerr << "  overlap:" << j[SIDE_2_OVERLAP] << endl;
  }
  
  if (settings.junction_debug) {
    ofile << "==SIDE 1==" << endl;
    ofile << "  position:" << j[SIDE_1_POSITION] << endl;
    ofile << "  strand:" << j[SIDE_1_STRAND] << endl;
    ofile << "  flanking:" << j["flanking_left"] << endl;
    ofile << "  continuation:" << j["side_1_continuation"] << endl;
    ofile << "  overlap:" << j[SIDE_1_OVERLAP] << endl;

    ofile << "==OVERLAP== " << j[ALIGNMENT_OVERLAP] << endl;
    
    ofile << "==SIDE 2==" << endl;
    ofile << "  position:" << j[SIDE_2_POSITION] << endl;
    ofile << "  strand:" << j[SIDE_2_STRAND] << endl;
    ofile << "  flanking:" << j["flanking_right"] << endl;
    ofile << "  continuation:" << j["side_2_continuation"] << endl;
    ofile << "  overlap:" << j[SIDE_2_OVERLAP] << endl;
  }
  
  // New Junction
  start = from_string<uint32_t>(j["flanking_left"]) - side_1_continuation;
  end = start + abs(from_string<int32_t>(j[ALIGNMENT_OVERLAP])) + 1 + side_2_continuation;
  start -= minimum_side_match_correction;
  end += minimum_side_match_correction;
  
  int32_t alignment_overlap = from_string<int32_t>(j[ALIGNMENT_OVERLAP]);
  int32_t non_negative_alignment_overlap = alignment_overlap;
  non_negative_alignment_overlap = max(0, non_negative_alignment_overlap);
  
  if (debug_output) {
    j["junction_start_pos_for_counting"] = to_string(start);
    j["junction_end_pos_for_counting"] = to_string(end);
  }
  j["junction_possible_overlap_registers_before_trimming"] = to_string(read_length_avg - abs(end - start));
  j["junction_possible_overlap_registers"] = (junction_jrc.get() != NULL) ? to_string(junction_jrc->count_confident_overlap_registers(j["key"], start, end, read_length_avg)) : "0";

  if (settings.junction_debug) ofile << "JUNCTION: start " << start << " end " << end << endl;
  if (verbose) cerr << "JUNCTION: start " << start << " end " << end << endl;

  j[NEW_JUNCTION_READ_COUNT] = (junction_jrc.get() != NULL) ? to_string(junction_jrc->count(j["key"], start, end, empty_read_names, junction_read_names)) : "0";
  
  if (settings.junction_debug) {
    ofile << "JUNCTION" << endl;
    for (map<string,bool>::iterator it = junction_read_names.begin(); it != junction_read_names.end(); it++) {
      ofile << it->first << endl;
    }
  }
  
  // @JEB 2022-10-25 +GGG JC mutation not assigned correct frequency
  // Is this a junction that adds bases between two adjacent nucleotides?
  // in this case we need to add the continuation on each side to the opposite side
  // ---------|GGGGGG---------   +GGG at |
  // ----------GGGGGG            side_1_continuation = 6   DON'T COUNT
  //           GGGGGG---------   side_2_continuation = 0   DON'T COUNT
  
  bool adjacent_junction = false;
  if (j[SIDE_1_SEQ_ID] == j[SIDE_2_SEQ_ID]) {
    if (from_string<int32_t>(j[SIDE_1_POSITION]) + 1 == from_string<int32_t>(j[SIDE_2_POSITION])) {
      adjacent_junction = true;
    }
  }
    
  // New side 1
  if ( (j[SIDE_1_REDUNDANT] != "1") && (j["side_1_annotate_key"] != "repeat") ) {
    int32_t side_1_strand = from_string<int32_t>(j[SIDE_1_STRAND]);
    start = from_string<uint32_t>(j[SIDE_1_POSITION]);
    end = start;
    int32_t overlap_correction = non_negative_alignment_overlap - from_string<int32_t>(j[SIDE_1_OVERLAP]);
      

    if (side_1_strand == +1) {
      start = start - 1;
      start -= overlap_correction;
      end += extra_stranded_require_overlap;
      start -= minimum_side_match_correction;
      end += minimum_side_match_correction;
      start -= side_1_continuation;
      if (adjacent_junction) {
        end += side_2_continuation;
      }
    } else {
      end = end + 1;
      start -= extra_stranded_require_overlap;
      end += overlap_correction;
      start -= minimum_side_match_correction;
      end += minimum_side_match_correction;
      end += side_1_continuation;
      if (adjacent_junction) {
        start -= side_2_continuation;
      }
    }
    
    if (settings.junction_debug) ofile << "SIDE 1: start " << start << " end " << end << endl;       
    if (verbose) cerr << "SIDE 1: start " << start << " end " << end << endl;
    if (debug_output) {
      j["side_1_start_pos_for_counting"] = to_string(start);
      j["side_1_end_pos_for_counting"] = to_string(end);
    }
    
    j["side_1_possible_overlap_registers_before_trimming"] = to_string(read_length_avg - abs(end - start));
    j["side_1_possible_overlap_registers"] = (reference_jrc.get() != NULL) ? to_string(reference_jrc->count_confident_overlap_registers(j[SIDE_1_SEQ_ID], start, end, read_length_avg)) : "0";

    j[SIDE_1_READ_COUNT] = (reference_jrc.get() != NULL) ? to_string(reference_jrc->count(j[SIDE_1_SEQ_ID], start, end, junction_read_names, empty_read_names)) : "0";
    
    if (settings.junction_debug) {
      ofile << "SIDE_1" << endl;
      for (map<string,bool>::iterator it = empty_read_names.begin(); it != empty_read_names.end(); it++) {
        ofile << it->first << endl;
      }
    }
  
  } else {
    j[SIDE_1_READ_COUNT] = "NA";
    j["side_1_possible_overlap_registers"] = "NA";
    j["side_1_possible_overlap_registers_before_trimming"] = "NA";
  }

  // New side 2
  int32_t side_2_possibilities = 0;

  if ( (j[SIDE_2_REDUNDANT] != "1") && (j["side_2_annotate_key"] != "repeat") ) {
    
    int32_t side_2_strand = from_string<int32_t>(j[SIDE_2_STRAND]);
    start = from_string<uint32_t>(j[SIDE_2_POSITION]);
    end = start;
    int32_t overlap_correction = non_negative_alignment_overlap - from_string<int32_t>(j[SIDE_2_OVERLAP]);
    
    if (side_2_strand == +1) {
      start = start - 1;
      start -= overlap_correction;
      end += extra_stranded_require_overlap;
      start -= minimum_side_match_correction;
      end += minimum_side_match_correction;
      start -= side_2_continuation;
      if (adjacent_junction) {
        end += side_1_continuation;
      }
    } else {
      end = end + 1;
      start -= extra_stranded_require_overlap;
      end += overlap_correction;
      start -= minimum_side_match_correction;
      end += minimum_side_match_correction;
      end += side_2_continuation;
      if (adjacent_junction) {
        start -= side_2_continuation;
      }
    }
    
    if (settings.junction_debug) ofile << "SIDE 2: start " << start << " end " << end << endl;
    if (verbose) cerr << "SIDE 2: start " << start << " end " << end << endl;
    if (debug_output) {
      j["side_2_start_pos_for_counting"] = to_string(start);
      j["side_2_end_pos_for_counting"] = to_string(end);
    }
    
    j["side_2_possible_overlap_registers_before_trimming"] = to_string(read_length_avg - abs(end - start));
    j["side_2_possible_overlap_registers"] = (reference_jrc.get() != NULL) ? to_string(reference_jrc->count_confident_overlap_registers(j[SIDE_2_SEQ_ID], start, end, read_length_avg)) : "0";

    j[SIDE_2_READ_COUNT] = (reference_jrc.get() != NULL) ? to_string(reference_jrc->count(j[SIDE_2_SEQ_ID], start, end, junction_read_names, empty_read_names)) : "0";

    if (settings.junction_debug) {
      ofile << "SIDE_2" << endl;
      for (map<string,bool>::iterator it = empty_read_names.begin(); it != empty_read_names.end(); it++) {
        ofile << it->first << endl;
      }
    }

  } else {
    j[SIDE_2_READ_COUNT] = "NA";
    j["side_2_possible_overlap_registers"] = "NA";
    j["side_2_possible_overlap_registers_before_trimming"] = "NA";
  }
  
  if (verbose) {
    cerr << "==Possibilities for overlapped reads==" << endl;
    cerr << "Side 1  : " << j["side_1_possible_overlap_registers"] << endl;
    cerr << "Junction: " << j["junction_possible_overlap_registers"] << endl;
    cerr << "Side 2  : " << j["side_2_possible_overlap_registers"] << endl;
  }
  if (settings.junction_debug) {
    ofile << "==Possibilities for overlapped reads==" << endl;
    ofile << "Side 1  : " << j["side_1_possible_overlap_registers"] << endl;
    ofile << "Junction: " << j["junction_possible_overlap_registers"] << endl;
    ofile << "Side 2  : " << j["side_2_possible_overlap_registers"] << endl;
  }
  
  //@ded calculate frequency of each junction. //
  double a, b;

  // A possible_overlap_registers of zero (the candidate-junction pre-filter in
  // candidate_junctions.cpp should normally have already rejected these, but this is a
  // defensive backstop) means there is no valid denominator to normalize by -- treat it the
  // same as a read count of "NA" (no information), rather than dividing by zero.
  uint32_t junction_possible_overlap_registers = from_string<uint32_t>(j["junction_possible_overlap_registers"]);
  uint32_t side_1_possible_overlap_registers = (j[SIDE_1_READ_COUNT] != "NA") ? from_string<uint32_t>(j["side_1_possible_overlap_registers"]) : 0;
  uint32_t side_2_possible_overlap_registers = (j[SIDE_2_READ_COUNT] != "NA") ? from_string<uint32_t>(j["side_2_possible_overlap_registers"]) : 0;

  double c = 0;
  bool have_c = (junction_possible_overlap_registers != 0);
  if (have_c) {
    c = from_string<uint32_t>(j[NEW_JUNCTION_READ_COUNT]);
    c /= static_cast<double>(junction_possible_overlap_registers);
  }

  // @JEB: divide the side X counts by 2, if both were counted
  // or by 1 if that side of the alignment was ambiguous (for edges of IS-elements)
  double d = 2;
  if ((j[SIDE_1_READ_COUNT] == "NA") || (side_1_possible_overlap_registers == 0)) {
    a = 0; //"NA" in read count sets value to 1 not 0
    d--;
  } else {
    a = from_string<uint32_t>(j[SIDE_1_READ_COUNT]);
    a /= static_cast<double>(side_1_possible_overlap_registers);
  }

  if ((j[SIDE_2_READ_COUNT] == "NA") || (side_2_possible_overlap_registers == 0)) {
    b = 0;
    d--;
  } else {
    b = from_string<uint32_t>(j[SIDE_2_READ_COUNT]);
    b /= static_cast<double>(side_2_possible_overlap_registers);
  }

  // We cannot assign a frequency if the denominator is zero
  if ((d == 0) || !have_c) {
    j[POLYMORPHISM_FREQUENCY] = "NA";
    j[FREQUENCY] = j[POLYMORPHISM_FREQUENCY];

    j[PREDICTION] = "unknown";
  } else {
    double new_junction_frequency_value = c /(c + ((a+b)/d) );
    j[POLYMORPHISM_FREQUENCY] = to_string(new_junction_frequency_value, settings.polymorphism_precision_places, true);
    j[FREQUENCY] = j[POLYMORPHISM_FREQUENCY];

    // Determine what kind of prediction we are
    
    // We may have added FREQUENCY_CUTOFF previously, so clear it.
    j.remove_reject_reason("FREQUENCY_CUTOFF");
    
    // Order of precedence depends on mode
    if (settings.polymorphism_prediction) {
      
      // Above 1-cutoff, we reject unless we change to consensus
      if (new_junction_frequency_value > 1.0 - settings.polymorphism_frequency_cutoff - settings.polymorphism_precision_decimal) {
        
        if (new_junction_frequency_value >= settings.consensus_frequency_cutoff - settings.polymorphism_precision_decimal) {
          j[PREDICTION] = "consensus";
          j[FREQUENCY] = "1";
        } else {
          j.add_reject_reason("FREQUENCY_CUTOFF");
          j[PREDICTION] = "polymorphism";
        }
      }
      // Below the cutoff, just reject
      else if (new_junction_frequency_value < settings.polymorphism_frequency_cutoff - settings.polymorphism_precision_decimal) {
        j.add_reject_reason("FREQUENCY_CUTOFF");
        j[PREDICTION] = "polymorphism";
      } else {
        j[PREDICTION] = "polymorphism";
      }
      
    } else {
      if (new_junction_frequency_value >=  settings.consensus_frequency_cutoff - settings.polymorphism_precision_decimal) {
        j[PREDICTION] = "consensus";
        j[FREQUENCY] = "1";
      } else if ((new_junction_frequency_value >= settings.polymorphism_frequency_cutoff - settings.polymorphism_precision_decimal) && (new_junction_frequency_value <= 1 - settings.polymorphism_frequency_cutoff - settings.polymorphism_precision_decimal)) {
        j[PREDICTION] = "polymorphism";
      } else {
        j[PREDICTION] = "polymorphism";
        j.add_reject_reason("FREQUENCY_CUTOFF");
      }
    }
  }
  
  // Finally assign average coverages based on fragments

  if ((j[SIDE_1_READ_COUNT] == "NA") || (side_1_possible_overlap_registers == 0)) {
      j[SIDE_1_COVERAGE] = "NA";
  }
  else {
    double side_1_correction = static_cast<double>(side_1_possible_overlap_registers) / read_length_avg;
    j[SIDE_1_COVERAGE] = to_string(from_string<double>(j[SIDE_1_READ_COUNT]) / summary.unique_coverage[j[SIDE_1_SEQ_ID]].average / side_1_correction, 2);
  }

  if ((j[SIDE_2_READ_COUNT] == "NA") || (side_2_possible_overlap_registers == 0)) {
    j[SIDE_2_COVERAGE] = "NA";
  }
  else {
    double side_2_correction = static_cast<double>(side_2_possible_overlap_registers) / read_length_avg;
    j[SIDE_2_COVERAGE] = to_string(from_string<double>(j[SIDE_2_READ_COUNT]) / summary.unique_coverage[j[SIDE_2_SEQ_ID]].average / side_2_correction, 2);
  }

  //corrects for overlap making it less likely for a read to span
  if (!have_c) {
    j[NEW_JUNCTION_COVERAGE] = "NA";
  } else {
    double overlap_correction = static_cast<double>(junction_possible_overlap_registers) / read_length_avg;
    double new_junction_average_read_count = (summary.unique_coverage[j[SIDE_1_SEQ_ID]].average + summary.unique_coverage[j[SIDE_2_SEQ_ID]].average) / 2;
    j[NEW_JUNCTION_COVERAGE] = to_string(from_string<double>(j[NEW_JUNCTION_READ_COUNT]) / new_junction_average_read_count / overlap_correction, 2);
  }

}
  
  
  
void  assign_junction_read_counts(
                                  const Settings& settings,
                                  Summary& summary,
                                  cGenomeDiff& gd
                                  )
{
  // Could be added as a parameter to reduce problems due to one-base mismatches
  static int32_t require_overlap = 0;
  
  diff_entry_list_t jc = gd.get_list(make_vector<gd_entry_type>(JC));

  if (jc.size() == 0) return;
  // Next calls can fail if there are no junctions (and therefore no FASTA file of junctions).
  
  // Create the file (to overwrite previous versions since we use append later)
  if (settings.junction_debug) {
    fstream ofile(settings.junction_debug_file_name.c_str(), ios_base::out);
    ASSERT(ofile.good(), "Could not open file " + settings.junction_debug_file_name);
  }
  
  
  // Fetch all the junction reads supporting
  // Keep track of how well they match the reference versus the putative new junctions.
  // right now this is in terms of mismatches (adding unmatched read bases as mismatches)
  
  // Create read counters in a way that they will automatically be cleaned up
  counted_ptr<junction_read_counter> reference_jrc(NULL);
  if (file_exists(settings.reference_bam_file_name.c_str())) {
    reference_jrc = counted_ptr<junction_read_counter>(new junction_read_counter(settings.reference_bam_file_name, settings.reference_fasta_file_name, settings.verbose, settings.reference_trim_file_name));
  }

  counted_ptr<junction_read_counter> junction_jrc(NULL);
  if (file_exists(settings.junction_bam_file_name.c_str()) && file_exists(settings.candidate_junction_fasta_file_name.c_str())) {
    junction_jrc = counted_ptr<junction_read_counter>(new junction_read_counter(settings.junction_bam_file_name, settings.candidate_junction_fasta_file_name, settings.verbose, settings.candidate_junction_trim_file_name));
  }
  
  for (diff_entry_list_t::iterator it = jc.begin(); it != jc.end(); it++) {
    cDiffEntry& j = **it;
    assign_one_junction_read_counts(settings, summary, j, reference_jrc, junction_jrc);
  }
  
}

  
void  assign_junction_read_coverage(
                                  const Settings& settings,
                                  Summary& summary,
                                  cGenomeDiff& gd
                                  )
{
  (void) settings;
  diff_entry_list_t jc = gd.get_list(make_vector<gd_entry_type>(JC));
  
  if (jc.size() == 0) return;
  // Next calls can fail if there are no junctions (and therefore no FASTA file of junctions).
  
    
  for (diff_entry_list_t::iterator it = jc.begin(); it != jc.end(); it++)
  {
    cDiffEntry& j = **it;
      
    // This sections just normalizes read counts to the average coverage of the correct sequence fragment
    
    double side_1_correction = (summary.sequence_conversion.read_length_avg - 1 - abs(from_string<double>(j["alignment_overlap"])) - from_string<double>(j["continuation_left"])) / (summary.sequence_conversion.read_length_avg - 1);
    
    if (j[SIDE_1_READ_COUNT] == "NA")
      j[SIDE_1_COVERAGE] = "NA";
    else
      j[SIDE_1_COVERAGE] = to_string(from_string<double>(j[SIDE_1_READ_COUNT]) / summary.unique_coverage[j[SIDE_1_SEQ_ID]].average / side_1_correction, 2);
    
    double side_2_correction = (summary.sequence_conversion.read_length_avg - 1 - abs(from_string<double>(j["alignment_overlap"])) - from_string<double>(j["continuation_right"])) / (summary.sequence_conversion.read_length_avg - 1);
    
    if (j[SIDE_2_READ_COUNT] == "NA")
      j[SIDE_2_COVERAGE] = "NA";
    else
      j[SIDE_2_COVERAGE] = to_string(from_string<double>(j[SIDE_2_READ_COUNT]) / summary.unique_coverage[j[SIDE_2_SEQ_ID]].average, 2);
    
    //corrects for overlap making it less likely for a read to span
    double overlap_correction = (summary.sequence_conversion.read_length_avg - 1 - abs(from_string<double>(j["alignment_overlap"])) - from_string<double>(j["continuation_left"]) - from_string<double>(j["continuation_right"])) / (summary.sequence_conversion.read_length_avg - 1);
    double new_junction_average_read_count = (summary.unique_coverage[j[SIDE_1_SEQ_ID]].average + summary.unique_coverage[j[SIDE_2_SEQ_ID]].average) / 2;
    
    j[NEW_JUNCTION_COVERAGE] = to_string(from_string<double>(j[NEW_JUNCTION_READ_COUNT]) / new_junction_average_read_count / overlap_correction, 2);
    
  }
}
  
  
junction_read_counter::junction_read_counter(const string& bam, const string& fasta, bool verbose, const string& trim_file_name)
  : pileup_base(bam, fasta), _verbose(verbose)
{
  _trims_list.resize(num_targets());
  for (uint32_t tid = 0; tid < num_targets(); tid++) {
    string this_file_name = Settings::file_name(trim_file_name, "@", target_name(tid));
    _trims_list[tid].ReadFile(this_file_name, target_length(tid));
  }
}

uint32_t junction_read_counter::count_confident_overlap_registers(
                                      const string& seq_id,
                                      const int32_t window_start,
                                      const int32_t window_end,
                                      const uint32_t read_length_avg
                                      ) const
{
  int32_t tid = seq_id_to_target_id(seq_id);
  if (tid < 0) return 0;

  return breseq::count_confident_overlap_registers(_trims_list[tid], target_length(static_cast<uint32_t>(tid)), window_start, window_end, read_length_avg);
}

uint32_t junction_read_counter::count(
                                      const string& seq_id,
                                      const int32_t start,
                                      const int32_t end,
                                      const map<string,bool> ignore_read_names,
                                      map<string,bool>& counted_read_names
                                      )
{
  _verbose = false; //for checking
  _ignore_read_names = ignore_read_names;
  _counted_read_names.clear();
  
  
  _count = 0;
  _start = start;
  _end = end;
  
  // it's possible that we will be sent negative values (by design)
  if (_start < 1) {
    return _count;
  }
  
    
  // The reads we overlap need to overlap both the start and the end to count.
  // We can retrieve them all by extracting things that overlap the start 
  // (and then checking if they also overlap the end).
  string region = seq_id + ":" + to_string(start) + "-" + to_string(start);
  
  if (_verbose) cerr << "junction_read_counter::count " << seq_id << ":" << start << "-" << end << endl;
  
  do_fetch(region);
  
  if (_verbose) cerr << "COUNT: " << _count << endl;

  
  counted_read_names = _counted_read_names;
  return _count;
}
  
void junction_read_counter::fetch_callback ( const alignment_wrapper& a )
{
  // The target_id will always be right.
  // Just check to be sure the start and end of the alignment go across the desired start and end.
  
  if (_verbose) cout << "  " << a.read_name();
  
  // Store the scores in a hash that can be resolved to see whether the read would have gone to the junction
  // or the reference . We can count. 
  // For certain kinds of junctions, we need to know how far they are identical on the end to properly normalize the others for that overlap.

  
  // read is to be ignored
  if (_ignore_read_names.count(a.read_name())
      ||a.read_name().find("-M1") != string::npos
      ||a.read_name().find("-M2") != string::npos
      )   
  {
    if (_verbose) cout << "  IGNORED" << endl;
    
    //ERROR("Read being counted twice for junction.");
    return;
  }
  
  // Don't count redundant
  if (a.redundancy() > 1) {
    if (_verbose) cout << "  REDUNDANT" << endl;
    return;
  }
  
  uint32_t q_start, q_end;
  a.reference_bounds_1(q_start, q_end);

  // Shrink the reference span inward by any trimmed (untrustworthy, e.g. repeat-ambiguous)
  // bases at the ends of the read, so a read that only *appears* to reach a coordinate
  // because of a repeat isn't counted as confidently spanning it. trim_left()/trim_right()
  // read the XL/XR tags (see get_alignment_trims() in calculate_trims.cpp), which are counted
  // from the absolute ends of the read and already include any existing soft-clip offset.
  uint32_t left_extra = a.trim_left();
  left_extra = (left_extra > a.query_start_0()) ? (left_extra - a.query_start_0()) : 0;
  uint32_t right_extra = a.trim_right();
  right_extra = (right_extra > (a.read_length() - a.query_end_1())) ? (right_extra - (a.read_length() - a.query_end_1())) : 0;

  q_start += left_extra;
  q_end = (right_extra < q_end) ? q_end - right_extra : 0;

  if (q_start > q_end) {
    if (_verbose) cout << "  TRIMMED AWAY" << endl;
    return;
  }

  if (_verbose) cout << "  " << q_start << "-" << q_end << "  " << _start << "-" << _end;

  if ((static_cast<int32_t>(q_start) <= _start) && (static_cast<int32_t>(q_end) >= _end)) {
    if (_verbose) cout << "  COUNTED" << endl;
    _count++;
  } else {
    if (_verbose) cout << "  NO OVERLAP" << endl;
    return;
  }
    
  // record that we counted this read
  _counted_read_names[a.read_name()] = true;
}
  
  
} // namespace breseq

