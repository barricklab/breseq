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
#include "coverage_distribution.h"

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

  // Cache each coverage group's empirical unique-coverage histogram (read once per group) plus the
  // empirical-vs-nbinom decision. The histogram comes from the preprocess coverage tab
  // (candidate_junction_path/@.unique_only_coverage_distribution.tab), which is kept alive to this
  // stage (see breseq_cmdline.cpp -- re-keyed to alignment_correction_done_file_name).
  map<uint32_t, coverage_histogram_t> group_hist;

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
          cov.nbinom_mean_parameter,
          false, // use_empirical (set below)
          0,     // deletion_floor
          0.0,   // coverage_hist_total
          vector<double>()
    };

    // Load this seq_id's coverage-group histogram and decide empirical vs nbinom. Only meaningful when
    // there is a fitted coverage distribution (average > 0); otherwise probability() short-circuits.
    if (cov.nbinom_mean_parameter > 0) {
      uint32_t group = settings.seq_id_to_coverage_group(seq_id);
      if (!group_hist.count(group)) {
        string hist_file = settings.file_name(settings.coverage_junction_distribution_file_name, "@", to_string<uint32_t>(group));
        if (file_exists(hist_file.c_str())) {
          group_hist[group] = read_coverage_histogram(hist_file);
        } else {
          group_hist[group] = coverage_histogram_t(); // empty -> stays on nbinom
        }
      }
      const coverage_histogram_t& hist = group_hist[group];
      if (!hist.n.empty()) {
        double N = 0.0; int32_t deletion_floor = 0;
        bool use_emp = use_empirical_pos_hash_coverage(hist, cov.nbinom_mean_parameter,
                                                       settings.junction_pos_hash_neg_log10_p_value_cutoff,
                                                       N, deletion_floor);
        p.use_empirical = use_emp;
        p.deletion_floor = deletion_floor;
        p.coverage_hist_total = N;
        if (use_emp) p.coverage_hist = hist.n;
      }
    }

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
    cout << "Coverage model: " << (p.use_empirical ? "EMPIRICAL" : "negative binomial") << endl;
    cout << "Negative Binomial Fit: Size = " << p.negative_binomial_size << " Prob = " << p.negative_binomial_prob << endl;
  }

  // Marginalize over local coverage. The coverage WEIGHT is either the empirical unique-coverage
  // histogram (faithful low-coverage tail; used for large references) or the fitted negative binomial
  // (parametric tail with no 1/N floor; used for small references). Everything else -- the per-position
  // read-start chance and the binomial pos_hash CDF -- is identical for both.
  if (p.use_empirical) {

    // Bin the histogram like the nbinom path so very high coverage stays bounded.
    int32_t hist_max = static_cast<int32_t>(p.coverage_hist.size()) - 1;
    const int32_t target_num_bins_for_estimate = 10000;
    int32_t bin_size = trunc(static_cast<double>(hist_max)/target_num_bins_for_estimate)+1;

    for (int32_t bin_start = p.deletion_floor + 1; bin_start <= hist_max; bin_start += bin_size) {
      int32_t bin_end = min(bin_start + bin_size - 1, hist_max);
      double bin_count = 0;
      for (int32_t c = bin_start; c <= bin_end; c++) bin_count += p.coverage_hist[c];
      if (bin_count <= 0) continue;

      double this_cov_pr = bin_count / p.coverage_hist_total;
      double this_coverage_middle = (bin_start + bin_end) / 2.0;
      double this_ratio_of_coverage_to_average = this_coverage_middle / static_cast<double>(p.average_coverage);
      double this_chance_per_pos_strand_read_start = 1 - pow(p.chance_per_pos_strand_no_read_start, this_ratio_of_coverage_to_average);

      double this_pos_hash_pr = 0;
      for (uint32_t i=0; i <= pos_hash_score; i++) {
        //chance of getting pos_hash_score or lower
        this_pos_hash_pr += binomial(this_chance_per_pos_strand_read_start, max_pos_hash_score, i);
      }

      pr += this_cov_pr * this_pos_hash_pr;
    }

  } else {

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
  }
  double log_pr = -log(pr)/log(10);
  probability_table[seq_id][pos_hash_score][max_pos_hash_score] = log_pr;
  
  if (verbose) {
    cout << "  -Log10 Cumulative Pr: " << log_pr << endl;
  }
  
  return log_pr;
}


// Defined later (after the static resolve helpers it uses): votes the specific IS copy per unique
// locus and writes the held pairs into the still-open resolved reference SAM.
void write_held_discordant_pairs(Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info,
                                 const SequenceTrimsList& trims_list, bam_file& resolved_reference_tam,
                                 vector<HeldDiscordantPair>& held_discordant_pairs);

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
    
  // These two become data/reference.bam and 06_bam/junction.bam (samtools sort preserves @RG), so
  // this is where read groups reach the output a user -- or CNery, via bam2cov -- actually reads.
  const cReadGroupList read_groups = settings.read_groups();
  bam_file resolved_reference_tam(settings.resolved_reference_sam_file_name, settings.reference_fasta_file_name, ios::out, read_groups);
  bam_file resolved_junction_tam(settings.resolved_junction_sam_file_name, settings.candidate_junction_fasta_file_name, ios::out, read_groups);
  
  settings.track_intermediate_file(settings.bam_done_file_name, settings.resolved_reference_sam_file_name);
  settings.track_intermediate_file(settings.bam_done_file_name, settings.resolved_junction_sam_file_name);
  
  UniqueJunctionMatchMap unique_junction_match_map;    // Map of junction_id to MatchedJunction
	RepeatJunctionMatchMap repeat_junction_match_map;  // Map of junction_id to read_name to MatchedJunction

  // Ambiguous-IS-side discordant pairs held aside during the streaming pass; their specific copy is
  // chosen by a per-locus vote and written back in the post-streaming merge-back (see below).
  vector<HeldDiscordantPair> held_discordant_pairs;

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
                            resolved_reference_tam,
                            held_discordant_pairs
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
	// Failing (marginal) junctions are resolved in a second pass, after every passing
	// junction has claimed its reads, so that failing-junction read assignment cannot
	// affect the passing junctions.
	list<JunctionTestInfo> deferred_failed_junction_test_info_list;
  
  ////
  // Score all of the matches.
  ////
  
  // MEASUREMENT ONLY (--junction-debug): create/truncate the per-read scoring table once here,
  // since score_junction() appends to it per candidate.
  if (settings.junction_debug) {
    ofstream mapq_debug_header(settings.junction_mapq_debug_file_name.c_str(), ios_base::out);
    ASSERT(mapq_debug_header.good(), "Could not open file " + settings.junction_mapq_debug_file_name);
    mapq_debug_header << "junction_id\tread_name\tdelta_this_candidate\tdelta_best_candidate"
                         "\tbest_reference_score\tas_this_candidate\tis_best_tag\tdegenerate_count"
                         "\tthis_left\tthis_right\tread_length\tn_ref_placements"
                         "\tlog10_odds_for_junction\tweight" << endl;
  }

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

  // Runs the pass/fail guard battery on a (freshly scored) junction: appends any
  // reject reasons to junction_test_info.reject_reasons, sets its neg_log10_pos_hash_p_value,
  // and returns whether the junction should be recorded in the genome diff at all.
  // Used for the passing junctions in the first pass and again for the failing
  // (marginal) junctions in the second pass, so both are judged by identical criteria --
  // there, on the basis of the reads each junction keeps after deduplication.
  auto test_junction_guards = [&](JunctionTestInfo& junction_test_info, ResolveJunctionInfo& junction_info) -> bool {

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

    return record;
  };

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

    // Test the best-scoring junction. (record is recomputed for deferred failing
    // junctions in the second pass, so it is not stored here.)
    test_junction_guards(junction_test_info, junction_info);

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
    
    if (!junction_test_info.reject_reasons.size()) {
      // PASSED: resolve now. This claims the junction's reads and removes them from all
      // inferior candidates (including failing ones), exactly as before.
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
                       false, // passed
                       junction_test_info.total_non_overlap_reads > 0
                       );
      // record is always true here (no reject reasons)
      passed_junction_test_info_list.push_back(junction_test_info);
    }
    else {
      // FAILED: defer resolution to the second pass. Deferring leaves this junction's
      // reads in the maps (just as the old failing branch did, which never purged), so
      // every passing junction still sees the identical map state -- passing results are
      // unchanged.
      deferred_failed_junction_test_info_list.push_back(junction_test_info);
    }

    junction_test_info_list.pop_back();

    // @JEB 2019-06-12 List is now sorted after each one accepted to ensure best results
    junction_test_info_list.sort();
  }

  ///
  // Resolve failing (marginal) junctions, best first.
  //
  // Now that every passing junction has claimed its reads, resolve the failing junctions
  // the same winner-take-all way: the best remaining marginal junction claims each shared
  // read and it is removed from inferior marginal candidates. Because this happens only
  // after the passing pass, it cannot change any passing-junction assignment.
  ///
  deferred_failed_junction_test_info_list.sort();

  while (!deferred_failed_junction_test_info_list.empty()) {

    JunctionTestInfo& junction_test_info = deferred_failed_junction_test_info_list.back();
    const string& junction_id = junction_test_info.junction_id;

    // Re-score against the current maps so pos_hash_score, read counts, and coverage
    // reflect only the reads this marginal junction actually keeps (after passing junctions
    // and better marginals took theirs). score_junction() overwrites the whole struct,
    // including reject_reasons.
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

    // Re-run the pass/fail guards on the deduplicated scores. A marginal junction's support
    // can only shrink here (reads only leave), so it can never turn into a passing junction;
    // but one whose deduplicated pos_hash_score now falls below the recording threshold is
    // dropped entirely (record == false) instead of being shown as a marginal supported by
    // reads that really belong to a better junction.
    ResolveJunctionInfo junction_info(junction_id);
    bool record = test_junction_guards(junction_test_info, junction_info);

    // Resolve with the winner-take-all failing behavior: claim reads for this junction and
    // remove them from inferior marginal candidates. Done even when the junction is not
    // recorded, so its reads are still deduplicated away from the other marginals.
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
                     junction_test_info.reject_reasons.size(), // failed
                     junction_test_info.total_non_overlap_reads > 0
                     );

    if (record) {
      rejected_junction_test_info_list.push_back(junction_test_info);
    }

    deferred_failed_junction_test_info_list.pop_back();

    // Re-sort so the next best-remaining marginal (whose score may have dropped after this
    // one claimed shared reads) is processed next.
    deferred_failed_junction_test_info_list.sort();
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
			const read_group_ref& rg = match.rg;
      
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
					rg,
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

  // Post-streaming DP merge-back: vote the specific IS copy per unique locus and write the held
  // ambiguous discordant pairs into the still-open resolved reference SAM (see the helper below).
  write_held_discordant_pairs(settings, summary, ref_seq_info, trims_list, resolved_reference_tam, held_discordant_pairs);

  // Save summary statistics
	summary.alignment_resolution.accepted_pos_hash_score_distribution = accepted_pos_hash_score_distribution;
  
  // Write the genome diff file
	gd.write(settings.jc_genome_diff_file_name);
  settings.track_intermediate_file(settings.output_done_file_name, settings.jc_genome_diff_file_name);

}
    
// Re-weight one window's reads at a given junction frequency and return the two evidence sums.
// The per-read prior odds AGAINST the junction are (1-f)/f -- the frequency itself -- which is
// the only self-consistent choice: the question "did THIS read come from the junction" has prior
// exactly the fraction of molecules that carry it.
//
// A fixed prior cannot serve here. The intuition that most candidate junctions are false is a
// statement about a JUNCTION, not about a read, so applying it per read charges it once for every
// read of support and deflates a junction's frequency in proportion to its own depth. Measured on
// 35 bp polymorphism data, a fixed prior of 512:1 lost 22 of 118 accepted junctions and no value
// of the constant avoided it -- the error is in where the prior is applied, not in its size.
// Suppressing spurious junctions is left to the junction-level tests that already exist
// (pos_hash_score and the confidence bounds below).
static void reweight_window(const vector<junction_read_counter::read_evidence_t>& evidence,
                            double f, double& sum_own, double& sum_other, double& sum_own_sq,
                            bool own_is_junction)
{
  sum_own = 0.0; sum_other = 0.0; sum_own_sq = 0.0;
  // Guard the odds at the extremes so a frequency pinned at 0 or 1 cannot produce inf/NaN; at
  // these magnitudes the weights are saturated anyway.
  double prior_odds_against_junction = (f <= 1e-12) ? 1e12 : ((f >= 1.0 - 1e-12) ? 1e-12 : (1.0 - f) / f);

  for (size_t i = 0; i < evidence.size(); i++) {
    double w_own = 1.0;
    if (evidence[i].have_odds) {
      // The ratio is ALWAYS stored in favour of the junction, on both windows. fetch_callback folds
      // the reference window's sign flip in when it records the value, so flipping again here
      // would turn every decisively-reference read into a decisively-junction one -- which drives
      // every junction to frequency 1. Only the ANSWER is flipped, to name this window's side.
      double w_junction = junction_read_weight(evidence[i].log10_odds_for_junction, evidence[i].n_ref_placements,
                                               prior_odds_against_junction);
      w_own = own_is_junction ? w_junction : (1.0 - w_junction);
    }
    sum_own += w_own;
    sum_own_sq += w_own * w_own;
    sum_other += 1.0 - w_own;
  }
}

// Best (least negative) quality-aware log-likelihood over a read's alignments to one reference.
// kNoHypothesisLogLikelihood means the read had no alignment there at all, which is different from
// having a bad one and must not be conflated with a numeric likelihood -- 0 would read as a
// PERFECT match, the opposite of what is meant.
static const int32_t kNoHypothesisLogLikelihood = INT32_MIN;

// The likelihoods below are computed and stamped UNCONDITIONALLY, including under
// --junction-no-read-weighting, which does not read them.
//
// Skipping the work when the weighting is off looks free -- it is the only cost this feature adds to
// a stage that walks every alignment of every read, and it saves ~10% of that stage. It is not free.
// The tags are written here, during alignment resolution, but consumed much later in the Output
// stage, and re-running Output alone is the normal way to compare weighting configurations cheaply.
// Make the tags conditional and an Output-only re-run WITH weighting, on top of a resolution stage
// produced WITHOUT it, finds no tags and silently emits unweighted frequencies -- the failure this
// model has already hit three times, and one that no test can catch because an unweighted read
// conserves perfectly.
//
// The saving also applies only to a diagnostic flag: with the weighting on, which is the default,
// the work happens either way. So the trade was a speedup on a path nobody takes for a silent wrong
// answer on one they do. Keeping this unconditional makes "the tags are present" an invariant of any
// resolution stage, which is why no runtime check for their absence is needed.

// Fills per_alignment in list order. Kept as its own pass so the CIGAR walk runs once per
// alignment rather than once to find the best and again to stamp each one.
static void _hypothesis_log_likelihoods(alignment_list& alignments, const cReferenceSequences& ref_seq_info,
                                        vector<int32_t>& per_alignment)
{
  per_alignment.clear();
  per_alignment.reserve(alignments.size());
  for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++) {
    per_alignment.push_back(alignment_log10_likelihood(*(it->get()), ref_seq_info));
  }
}

// Best log-likelihood the read achieves under one hypothesis, over the alignments as they stood
// BEFORE eligible_read_alignments() pruned them -- which is the same list best_reference_score and
// best_junction_score are computed from, and it matters.
//
// eligible_read_alignments() ends with test_read_alignment_requirements(), a minimum-mapped-length
// rule that can empty the list while still returning a nonzero best score. The reads it empties are
// exactly the ones this model exists to weigh: a read spanning a breakpoint matches the reference
// over only part of its length, so its reference alignment is short and gets dropped. Reading the
// competing likelihood off the pruned list therefore found nothing for those reads and fell back to
// full weight -- leaving the weighting inert on its entire target population while every test
// stayed green, because an unweighted read conserves perfectly. Measured on Ara-2_60K before the
// fix: 8,831 of 13,029 junction-supporting reads had no competing likelihood on record.
//
// The length requirement is a routing and QC rule, not a claim about likelihood. A discarded
// alignment is still the best explanation the reference genome can offer for the read, and that is
// what the ratio needs.
static int32_t _best_hypothesis_log_likelihood(const alignment_list& alignments, const cReferenceSequences& ref_seq_info)
{
  int32_t best = kNoHypothesisLogLikelihood;
  for (alignment_list::const_iterator it = alignments.begin(); it != alignments.end(); it++) {
    if (it->get()->unmapped()) continue;   // no coordinates to score against
    int32_t ll = alignment_log10_likelihood(*(it->get()), ref_seq_info);
    if ((best == kNoHypothesisLogLikelihood) || (ll > best)) best = ll;
  }
  return best;
}

// Stamp each alignment with both terms of its junction-vs-reference likelihood ratio -- its own
// log-likelihood under the reference it is aligned to, and the best the same read achieved under
// the COMPETING hypothesis -- so the comparison survives into the BAMs and can be recovered during
// read counting in the Output stage (see kBreseqOwnHypothesisLogLikelihoodBAMTag, settings.cpp).
// This is the only point in the pipeline where both hypotheses are in scope for the same read:
// reads that go on to support a passing junction are written only to junction.bam, so their
// reference alignments -- and with them any way to recompute the reference term -- are dropped.
//
// When the read had no alignment at all under the competing hypothesis, neither tag is written and
// the read falls back to full weight toward the BAM it sits in. Stamping a placeholder instead
// would be read downstream as a real likelihood.
static void _stamp_hypothesis_log_likelihoods(alignment_list& alignments,
                                              const vector<int32_t>& own_log_likelihoods,
                                              int32_t other_log_likelihood)
{
  if (other_log_likelihood == kNoHypothesisLogLikelihood) return;

  size_t i = 0;
  for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++, i++) {
    int32_t own_log_likelihood = own_log_likelihoods[i];
    it->get()->aux_set(kBreseqOwnHypothesisLogLikelihoodBAMTag, 'i', sizeof(int32_t), (uint8_t*)&own_log_likelihood);
    it->get()->aux_set(kBreseqOtherHypothesisLogLikelihoodBAMTag, 'i', sizeof(int32_t), (uint8_t*)&other_log_likelihood);
  }
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

  // Alignments as they stood before eligible_read_alignments() pruned them. The likelihood ratio is
  // read off these, for the reason spelled out on _best_hypothesis_log_likelihood(); copying is
  // cheap because the list holds counted pointers.
  alignment_list all_reference_alignments, all_junction_alignments;

  // Does this read have eligible reference sequence matches?
  if ((reference_alignments.size() > 0) && (seq.m_name == reference_alignments.front()->read_name()))
  {
    m.this_reference_alignments = reference_alignments;
    reference_tam.read_alignments(reference_alignments, false);
    all_reference_alignments = m.this_reference_alignments;
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

    // min_match_score of 0, NOT best_reference_score. Passing the reference score here erases
    // every junction alignment that fails to beat the reference, which destroys the score margin
    // for exactly the reads that end up counted on a junction's reference SIDES: they are left
    // with either an exact tie (difference 0) or an empty list, and an empty list makes
    // eligible_read_alignments return 0 so the difference becomes a meaningless
    // -best_reference_score. Keeping the true best junction score makes the difference a real
    // signed number on both sides of the comparison.
    //
    // This cannot change which reads route to junctions. Routing tests difference > 0, i.e. best
    // junction score > best reference score; any alignment meeting that already cleared the old
    // floor, so its score is unchanged. Only reads that still route to the reference see a
    // different (now meaningful) value, and their junction alignments are never written.
    //
    // Snapshot AFTER the overlap filter above and before the pruning, so the likelihood ratio's
    // junction term is taken over the same set best_junction_score is: an alignment that misses the
    // junction entirely is not evidence for it.
    all_junction_alignments = m.this_junction_alignments;
    m.best_junction_score = eligible_read_alignments(settings, junction_ref_seq_info, m.this_junction_alignments, settings.junction_allow_suboptimal_matches, 0);
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

  // Both likelihoods must be computed before either is stamped: each list needs the OTHER list's
  // best value, and stamping does not change it, but reading it after a partial stamp would invite
  // exactly the kind of ordering bug that silently produces a self-referential ratio.
  int32_t best_reference_log_likelihood = _best_hypothesis_log_likelihood(all_reference_alignments, ref_seq_info);
  int32_t best_junction_log_likelihood  = _best_hypothesis_log_likelihood(all_junction_alignments, junction_ref_seq_info);

  vector<int32_t> reference_log_likelihoods, junction_log_likelihoods;
  _hypothesis_log_likelihoods(m.this_reference_alignments, ref_seq_info, reference_log_likelihoods);
  _hypothesis_log_likelihoods(m.this_junction_alignments, junction_ref_seq_info, junction_log_likelihoods);

  _stamp_hypothesis_log_likelihoods(m.this_reference_alignments, reference_log_likelihoods, best_junction_log_likelihood);
  _stamp_hypothesis_log_likelihoods(m.this_junction_alignments, junction_log_likelihoods, best_reference_log_likelihood);

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
                                 const read_group_ref& rg,
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
    _write_reference_matches(settings, summary, ref_seq_info, trims_list, m.this_reference_alignments, resolved_reference_tam, rg);
  }
  else
  {
    JunctionMatchPtr junction_match_ptr(
                                        new JunctionMatch(
                                                          m.this_reference_alignments,   // reference sequence alignments
                                                          m.this_junction_alignments,    // the BEST candidate junction alignments
                                                          rg,              // index of the fastq file this read came from
                                                          m.mapping_quality_difference,  // difference between reference junction alignments (in # mismatches)
                                                          0,                             //
                                                          static_cast<int32_t>(m.best_reference_score)
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

  // Order by each mate's 5' end (forward -> start, reverse -> end), NOT by leftmost mapped
  // coordinate. Leftmost coordinate ties for short/overlapping fragments and would mislabel
  // normal FR fragments as RF; the 5' end gives the correct read geometry regardless of overlap.
  uint32_t a_5p = a->reversed() ? a_end : a_start;
  uint32_t b_5p = b->reversed() ? b_end : b_start;
  bool a_is_lower = (a_5p <= b_5p);
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

// Marks the MAPPED mate of a SINGLETON pair -- one mate reference-best, the other with no alignment
// anywhere (neither to the reference nor to any candidate junction). Without this a singleton is
// indistinguishable in the BAM from several unrelated cases, because mark_pair_info above is the only
// thing that ever sets BAM_FPAIRED and it runs only when BOTH mates resolved to the reference.
//
// Standard BAM semantics: BAM_FPAIRED | BAM_FMUNMAP, with RNEXT/PNEXT set to the read's own position
// (the samtools convention for a mate-unmapped record) and TLEN 0. Setting mtid is NOT optional:
// bam_file::write_alignments branches on (mtid == tid) and an untouched mtid of -1 reads back as
// 4294967295, indexing target_name[] out of bounds. BAM_FMREVERSE is deliberately left unset (it is
// undefined when 0x8 is set), and there is no XP tag -- with no second alignment there is no pair
// orientation to record.
//
// Every alignment in the list is marked, not just the first: write_alignments emits all of them, and
// mate-unmapped is a property of the READ, not of which copy it was placed at.
//
// This is what Missing Pair (MP) evidence keys on: a mate that maps only across a candidate junction
// still has mapped_anywhere == true, so a marked singleton means the mate's sequence is absent from
// the reference AND from every candidate junction -- i.e. genuinely novel.
static void mark_mate_unmapped(alignment_list& alignments)
{
  for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++)
  {
    bam_alignment* a = it->get();
    a->core.flag |= BAM_FPAIRED | BAM_FMUNMAP;
    a->core.mtid  = a->core.tid;
    a->core.mpos  = a->core.pos;
    a->core.isize = 0;
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

// Deterministic representative of a (possibly multiply-mapped) mate: the alignment with the lowest
// (reference target id, reference start). Used to break the tie for a redundant discordant pair so
// every read from the same multicopy element canonicalizes to ONE copy (see the tie-break block in
// resolve_reference_pair). A size-1 list returns its sole alignment.
static bam_alignment* pick_lowest_alignment(alignment_list& alignments)
{
  bam_alignment* best = NULL;
  for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++)
  {
    bam_alignment* a = it->get();
    if (best == NULL
        || a->reference_target_id() < best->reference_target_id()
        || (a->reference_target_id() == best->reference_target_id()
            && a->reference_start_1() < best->reference_start_1()))
      best = a;
  }
  return best;
}

// Repeat-copy identifier "seq_id:start-end" for the copy an alignment sits in (interior or within a
// short distance of a boundary), or "" if the alignment is not in an annotated repeat.
static string dp_copy_id_for_alignment(cReferenceSequences& ref_seq_info, bam_alignment* a)
{
  string seq_id = ref_seq_info[a->reference_target_id()].m_seq_id;
  int32_t md = 50;
  cFeatureLocation* f = cReferenceSequences::find_closest_repeat_region_boundary(
      a->reference_start_1(), ref_seq_info[seq_id].m_repeats, md, a->reversed() ? -1 : +1, true);
  if (f == NULL) return "";
  return seq_id + ":" + to_string(f->get_start_1()) + "-" + to_string(f->get_end_1());
}

#ifdef BRESEQ_WRITE_DISCORDANT_PAIRS_CSV
// Strand-aware 5' anchor position of an alignment: reference_start_1() if forward,
// reference_end_1() if reversed (reference_stranded_bounds_1 already computes exactly this).
static uint32_t stranded_anchor_position(bam_alignment* a)
{
  uint32_t start, end;
  a->reference_stranded_bounds_1(start, end);
  return start;
}
#endif

#ifdef BRESEQ_WRITE_DISCORDANT_PAIRS_CSV
// Single-character strand code used by the legacy discordant-pairs CSV: forward -> 'F', reverse -> 'R'.
static inline char strand_char(const bam_alignment* a) { return a->reversed() ? 'R' : 'F'; }

static void write_discordant_pair_row(
                                      ofstream& out,
                                      const string& read_number,
                                      const string& orientation,
                                      const string& seq_id_1,
                                      const string& start_1,
                                      const string& seq_id_2,
                                      const string& start_2,
                                      int64_t distance
                                      )
{
  out << read_number << "," << orientation << "," << seq_id_1 << "," << start_1 << "," << seq_id_2 << "," << start_2 << "," << distance << endl;
}
#endif

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
                              bam_file& resolved_reference_tam,
                              vector<HeldDiscordantPair>& held_discordant_pairs
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
  // One read group per read file SET, numbered over the sets that survived conversion and
  // --limit-fold-coverage filtering -- the same list Settings::read_groups() emits @RG lines for,
  // so this index and the header agree. Both mates of a paired set share it.
  uint32_t set_index = 0;

  // Per-reference-position difference array counting concordant read pairs whose inner gap spans each
  // position, pooled across all paired read-file sets. Indexed by BAM target id (consistent across sets
  // -- same reference). Lazily allocated on the first paired set; fit to a negative binomial after the
  // loop as the null distribution for scoring Discordant Pair evidence.
  vector<vector<int32_t> > crossing_diff;
  vector<string> crossing_seq_id;  // target id -> seq_id (captured so it survives header deletion)

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

      // A one-file set is a read group with no mates, so the writers leave its mate bits alone.
      const read_group_ref rg(set_index, false, false);

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

        // Alignments as they stood before eligible_read_alignments() pruned them; the likelihood
        // ratio is read off these. See _best_hypothesis_log_likelihood().
        alignment_list all_reference_alignments, all_junction_alignments;

        // Does this read have eligible reference sequence matches?
        alignment_list this_reference_alignments;
        if ((reference_alignments.size() > 0) && (seq.m_name == reference_alignments.front()->read_name()))
        {
          this_reference_alignments = reference_alignments;
          reference_tam->read_alignments(reference_alignments, false);

          if (verbose) {
            cerr << " Before Overlap Reference alignments = " << this_reference_alignments.size() << endl;
          }
          all_reference_alignments = this_reference_alignments;
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

          // min_match_score of 0, not best_reference_score -- see the matching comment in
          // resolve_one_mate() for why the reference score must not be used as a floor here.
          // Snapshot after the overlap filter, before the pruning, as resolve_one_mate() does.
          all_junction_alignments = this_junction_alignments;
          best_junction_score = eligible_read_alignments(settings, junction_ref_seq_info, this_junction_alignments, settings.junction_allow_suboptimal_matches, 0);
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

        // See the matching block in resolve_one_mate(): compute both terms of the likelihood ratio
        // before stamping either.
        int32_t best_reference_log_likelihood = _best_hypothesis_log_likelihood(all_reference_alignments, ref_seq_info);
        int32_t best_junction_log_likelihood  = _best_hypothesis_log_likelihood(all_junction_alignments, junction_ref_seq_info);

        vector<int32_t> reference_log_likelihoods, junction_log_likelihoods;
        _hypothesis_log_likelihoods(this_reference_alignments, ref_seq_info, reference_log_likelihoods);
        _hypothesis_log_likelihoods(this_junction_alignments, junction_ref_seq_info, junction_log_likelihoods);

        _stamp_hypothesis_log_likelihoods(this_reference_alignments, reference_log_likelihoods, best_junction_log_likelihood);
        _stamp_hypothesis_log_likelihoods(this_junction_alignments, junction_log_likelihoods, best_reference_log_likelihood);

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

          _write_reference_matches(settings, summary, ref_seq_info, trims_list, this_reference_alignments, resolved_reference_tam, rg);
        }
        else
        {
          if (verbose)
            cout << "Best alignment is to candidate junction. MQD: " << mapping_quality_difference << endl;

          JunctionMatchPtr junction_match_ptr(
                                              new JunctionMatch(
                                                                this_reference_alignments,    // reference sequence alignments
                                                                this_junction_alignments,     // the BEST candidate junction alignments
                                                                rg,             // index of the fastq file this read came from
                                                                mapping_quality_difference,   // difference between reference junction alignments (in # mismatches)
                                                                0,                            //
                                                                static_cast<int32_t>(best_reference_score)
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

      set_index++;
      continue;
    }

    ///
    //  PAIRED: stream both mates in lockstep so we can classify each read pair's orientation
    //  and distance against this read file set's majority orientation/distance cutoff.
    ///

    const cReadFile& rf1 = rfs.m_files[0];
    const cReadFile& rf2 = rfs.m_files[1];
    // Both mates land in the SAME read group and are told apart by BAM_FREAD1/BAM_FREAD2.
    const read_group_ref rg_1(set_index, false, true);
    const read_group_ref rg_2(set_index, true,  true);
    set_index++;

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

    // Lazily allocate the per-position concordant-pair-crossing accumulator on the first paired set
    // (the reference header is identical across sets, so target ids stay consistent).
    if (crossing_diff.empty()) {
      int32_t nt = reference_tam_1->bam_header->n_targets;
      crossing_diff.resize(nt);
      crossing_seq_id.resize(nt);
      for (int32_t t = 0; t < nt; t++) {
        string sid = reference_tam_1->bam_header->target_name[t];
        crossing_seq_id[t] = sid;
        crossing_diff[t].assign(ref_seq_info[sid].get_sequence_length() + 2, 0);
      }
    }

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
    // preprocessing = the PRELIMINARY fit) for this read file set. This is what we use to assign the
    // concordant/discordant BAM flags below.
    const string majority_orientation = summary.preliminary_paired_mapping_distance_distribution[rfs.m_base_name].majority_orientation;
    const double distance_cutoff = summary.preliminary_paired_mapping_distance_distribution[rfs.m_base_name].distance_cutoff;

    // Count read pairs that receive BAM pair flags here (the FINAL pass), split concordant/discordant.
    uint64_t mapped_pairs = 0;
    uint64_t concordant_pairs = 0;

#ifdef BRESEQ_WRITE_DISCORDANT_PAIRS_CSV
    string discordant_pairs_file_name = Settings::file_name(settings.discordant_pairs_file_name, "#", rfs.m_base_name);
    ofstream discordant_csv_out(discordant_pairs_file_name.c_str());
    ASSERT(discordant_csv_out.is_open(), "Could not write to file: " + discordant_pairs_file_name);
    discordant_csv_out << "read_number,orientation,seq_id_1,start_1,seq_id_2,start_2,distance" << endl;
#endif

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
#ifdef BRESEQ_WRITE_DISCORDANT_PAIRS_CSV
        else if (pairing.any_same_tid_combo_exists)
        {
          // Discordant: closest same-tid combination still fails the orientation/cutoff test.
          // Both mates are still written to resolved_reference_sam_file_name below (with their
          // original, non-downselected alignment lists) -- only the CSV logging differs here.
          // Orientation letters are written per-mate in start_1/start_2 order (not folded) so
          // the discordant-pairs plot can color each point by its own mate's true strand.
          string seq_id = reference_tam_1->bam_header->target_name[pairing.best_a->reference_target_id()];
          string orientation = string() + strand_char(pairing.best_a) + strand_char(pairing.best_b);
          write_discordant_pair_row(discordant_csv_out, seq1.m_name, orientation, seq_id,
                                     to_string(stranded_anchor_position(pairing.best_a)),
                                     seq_id,
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
          string seq_id_1 = reference_tam_1->bam_header->target_name[a1->reference_target_id()];
          string seq_id_2 = reference_tam_1->bam_header->target_name[a2->reference_target_id()];
          string orientation = string() + strand_char(a1) + strand_char(a2);
          write_discordant_pair_row(discordant_csv_out, seq1.m_name, orientation, seq_id_1,
                                     to_string(stranded_anchor_position(a1)),
                                     seq_id_2,
                                     to_string(stranded_anchor_position(a2)),
                                     -1);
        }
#endif

        // Whenever both mates end up with exactly one alignment (whether concordant or
        // discordant), mark them as paired in the BAM using standard fields. Recompute
        // orientation/distance fresh from these specific surviving alignments rather than
        // reusing pairing.best_a/best_b, which can point at a different alignment than the one
        // that survived downselection.
        bool held_aside = false;   // set when an ambiguous discordant pair is deferred (see below)
        if ((m1.this_reference_alignments.size() == 1) && (m2.this_reference_alignments.size() == 1))
        {
          bam_alignment* a1 = m1.this_reference_alignments.front().get();
          bam_alignment* a2 = m2.this_reference_alignments.front().get();
          bool same_tid = (a1->reference_target_id() == a2->reference_target_id());
          OrientationDistance od = same_tid ? compute_orientation_and_distance(a1, a2) : OrientationDistance{"NA", 0};
          mark_pair_info(a1, a2, same_tid, od.orientation, od.distance, pairing.any_concordant_combo_exists);

          // Count this flag-assigned pair for the final summary (concordant vs discordant).
          ++mapped_pairs;
          if (pairing.any_concordant_combo_exists) ++concordant_pairs;

          // Tally the concordant pair's inner gap into the per-position crossing accumulator (the null
          // distribution for DP scoring): +1 over [left_read.end+1, right_read.start-1].
          if (pairing.any_concordant_combo_exists && same_tid && !crossing_diff.empty()) {
            int32_t tid = a1->reference_target_id();
            int32_t s1 = a1->reference_start_1(), e1 = a1->reference_end_1();
            int32_t s2 = a2->reference_start_1(), e2 = a2->reference_end_1();
            int32_t gap_lo = ((s1 <= s2) ? e1 : e2) + 1;   // end of the earlier-starting read + 1
            int32_t gap_hi = ((s1 <= s2) ? s2 : s1) - 1;   // start of the later-starting read - 1
            vector<int32_t>& d = crossing_diff[tid];
            if (gap_hi >= gap_lo && gap_lo >= 1 && static_cast<size_t>(gap_hi) + 1 < d.size()) {
              d[gap_lo]++;
              d[gap_hi + 1]--;
            }
          }
        }
        else if (!pairing.any_concordant_combo_exists &&
                 ((m1.this_reference_alignments.size() > 1) || (m2.this_reference_alignments.size() > 1)))
        {
          // Redundant discordant pair: no concordant combination exists and at least one mate maps
          // to multiple copies (e.g. a multicopy IS element). Which copy the discordant-pair flags
          // land on determines which IS copy DP evidence reports, so we do NOT decide it per-pair
          // here (the lowest copy is often a sequence variant). When exactly one mate is unique, HOLD
          // THE PAIR ASIDE and let a global per-locus vote in the post-streaming merge-back pick the
          // read-supported copy. Both mates' full alignment lists are kept (X1>1 coverage preserved);
          // they are written there, not inline. Not counted in mapped_pairs/concordant_pairs.
          bool m1_unique = (m1.this_reference_alignments.size() == 1);
          bool m2_unique = (m2.this_reference_alignments.size() == 1);
          if (m1_unique != m2_unique) {
            HeldDiscordantPair h;
            MateResolution& um = m1_unique ? m1 : m2;
            MateResolution& rm = m1_unique ? m2 : m1;
            h.unique_alignments = um.this_reference_alignments;
            h.redundant_alignments = rm.this_reference_alignments;
            h.unique_rg = m1_unique ? rg_1 : rg_2;
            h.redundant_rg = m1_unique ? rg_2 : rg_1;
            bam_alignment* ua = um.this_reference_alignments.front().get();
            h.unique_seq_id = ref_seq_info[ua->reference_target_id()].m_seq_id;
            h.unique_position = ua->reference_start_1();
            h.window = distance_cutoff;
            held_discordant_pairs.push_back(h);
            held_aside = true;
          } else {
            // Both mates redundant (IS-to-IS): no unique anchor to vote by; fall back to marking each
            // mate's lowest copy inline (written below with the rest).
            bam_alignment* a1 = pick_lowest_alignment(m1.this_reference_alignments);
            bam_alignment* a2 = pick_lowest_alignment(m2.this_reference_alignments);
            bool same_tid = (a1->reference_target_id() == a2->reference_target_id());
            OrientationDistance od = same_tid ? compute_orientation_and_distance(a1, a2) : OrientationDistance{"NA", 0};
            mark_pair_info(a1, a2, same_tid, od.orientation, od.distance, /*is_concordant=*/false);
          }
        }

        // Held pairs are written in the post-streaming merge-back (after the copy vote); all others
        // are written inline here exactly as before.
        if (!held_aside) {
          _write_reference_matches(settings, summary, ref_seq_info, trims_list, m1.this_reference_alignments, resolved_reference_tam, rg_1);
          _write_reference_matches(settings, summary, ref_seq_info, trims_list, m2.this_reference_alignments, resolved_reference_tam, rg_2);
        }
      }
      else
      {
        // Either mate is junction-best, or one/both mates are fully unmapped -- these fall
        // through to the existing, unmodified per-mate handling. A singleton mapping (one mate
        // reference-best, the other fully unmapped) has its mapped mate's alignments flagged
        // BAM_FPAIRED|BAM_FMUNMAP below (see mark_mate_unmapped), and is still logged as
        // discordant for visibility; either way the alignment itself is written normally.
        bool m1_singleton_reference = m1_is_reference_match && !m2.mapped_anywhere;
        bool m2_singleton_reference = m2_is_reference_match && !m1.mapped_anywhere;

#ifdef BRESEQ_WRITE_DISCORDANT_PAIRS_CSV
        if (m1_singleton_reference || m2_singleton_reference)
        {
          bam_alignment* mapped_alignment = m1_singleton_reference
            ? m1.this_reference_alignments.front().get()
            : m2.this_reference_alignments.front().get();
          string seq_id = reference_tam_1->bam_header->target_name[mapped_alignment->reference_target_id()];
          // Singleton: only one mate maps. Write its single strand letter and put its position in
          // the first (seq_id_1/start_1) slot; the second slot is left empty.
          string orientation = string(1, strand_char(mapped_alignment));
          write_discordant_pair_row(discordant_csv_out, seq1.m_name, orientation, seq_id,
                                     to_string(stranded_anchor_position(mapped_alignment)),
                                     "", "", -1);
        }
#endif

        // Mark before dispatching: dispatch_mate_result routes a reference-best mate to
        // _write_reference_matches -> bam_file::write_alignments, which preserves these flags
        // (fix_flags strips only 0x100; 0x80 is then set from the read group). A junction-best
        // mate is never marked here, so the
        // flags can never reach write_moved_alignment's RNEXT/PNEXT handling.
        if (m1_singleton_reference) mark_mate_unmapped(m1.this_reference_alignments);
        if (m2_singleton_reference) mark_mate_unmapped(m2.this_reference_alignments);

        dispatch_mate_result(settings, summary, ref_seq_info, trims_list, resolved_reference_tam, junction_tam_1, rg_1, seq1.m_name, m1, all_junction_ids, unique_junction_match_map, repeat_junction_match_map);
        dispatch_mate_result(settings, summary, ref_seq_info, trims_list, resolved_reference_tam, junction_tam_2, rg_2, seq2.m_name, m2, all_junction_ids, unique_junction_match_map, repeat_junction_match_map);
      }
    } // End loop through every read pair

    {
      ostringstream progress_message;
      progress_message << "    READS:" << setw(12) << right << reads_processed;
      print_progress_line(progress_message.str());
    }
    end_progress_line();

    // Final paired-mapping-distance record: carry the fit fields from the preliminary (stage-03)
    // record and attach the actual concordant/mapped pair counts from this flag-assignment pass.
    {
      PairedMappingDistanceDistributionSummary final_pmdd = summary.preliminary_paired_mapping_distance_distribution[rfs.m_base_name];
      final_pmdd.mapped_pairs = static_cast<double>(mapped_pairs);
      final_pmdd.concordant_pairs = static_cast<double>(concordant_pairs);
      summary.paired_mapping_distance_distribution[rfs.m_base_name] = final_pmdd;
    }

    // save statistics
    summary.alignment_resolution.read_file[rf1.m_base_name] = read_file_summary_info_1;
    summary.alignment_resolution.read_file[rf2.m_base_name] = read_file_summary_info_2;

    summary.alignment_resolution.total_unmapped_reads += read_file_summary_info_1.num_unmapped_reads + read_file_summary_info_2.num_unmapped_reads;
    summary.alignment_resolution.total_unmapped_read_bases += read_file_summary_info_1.num_unmapped_read_bases + read_file_summary_info_2.num_unmapped_read_bases;

    summary.alignment_resolution.total_reads += read_file_summary_info_1.num_total_reads + read_file_summary_info_2.num_total_reads;
    summary.alignment_resolution.total_bases += read_file_summary_info_1.num_total_bases + read_file_summary_info_2.num_total_bases;

#ifdef BRESEQ_WRITE_DISCORDANT_PAIRS_CSV
    discordant_csv_out.close();
#endif

    if (junction_tam_1 != NULL) delete junction_tam_1;
    if (junction_tam_2 != NULL) delete junction_tam_2;
    delete reference_tam_1;
    delete reference_tam_2;

  } // End of Read File Set loop

  // Write each reference sequence's INTERIOR concordant-pair crossing histogram to a tab file (the null
  // distribution used to score DP evidence). Interior = >= distance_cutoff (max insert) from each end --
  // the only positions a pair can fully span (near-end positions are crossing-depleted). Prefix-sum the
  // difference array into per-position counts. Distributions are CSV intermediates, not summary JSON.
  double dist_cutoff = 0.0;
  for (PairedMappingDistanceDistributionSummaries::const_iterator it = summary.preliminary_paired_mapping_distance_distribution.begin();
       it != summary.preliminary_paired_mapping_distance_distribution.end(); it++)
    if (it->second.distance_cutoff > dist_cutoff) dist_cutoff = it->second.distance_cutoff;
  int32_t margin = static_cast<int32_t>(dist_cutoff);

  for (size_t t = 0; t < crossing_diff.size(); t++) {
    vector<int32_t>& d = crossing_diff[t];
    for (size_t pos = 1; pos < d.size(); pos++) d[pos] += d[pos - 1];
    int32_t seqlen = static_cast<int32_t>(d.size()) - 2;
    int32_t lo = margin + 1, hi = seqlen - margin;
    map<int32_t, uint64_t> hist;
    for (int32_t pos = lo; pos <= hi; pos++) {
      int32_t c = d[pos]; if (c < 0) c = 0;
      hist[c]++;
    }
    string fn = Settings::file_name(settings.concordant_pair_crossing_distribution_file_name, "#", crossing_seq_id[t]);
    ofstream out(fn.c_str());
    out << "crossing\tcount" << endl;
    int32_t maxc = hist.empty() ? 0 : hist.rbegin()->first;
    for (int32_t c = 0; c <= maxc; c++) out << c << "\t" << (hist.count(c) ? hist[c] : 0) << endl;
    out.close();
  }

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

  // Resolution is flattened per file here, but read groups are still per SET, as everywhere else --
  // so a flat file index has to be mapped back to its set plus which mate it is. Using the flat
  // index directly as a group index would name a group that does not exist whenever a set is paired.
  vector<read_group_ref> flat_read_groups;
  {
    uint32_t on_set_index = 0;
    for (cReadFileSets::iterator rfs_it = read_files.begin(); rfs_it != read_files.end(); rfs_it++, on_set_index++)
      for (size_t m = 0; m < rfs_it->m_files.size(); m++)
        flat_read_groups.push_back(read_group_ref(on_set_index, (m == 1), rfs_it->is_paired()));
  }

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
      _write_reference_matches(settings, summary, ref_seq_info, trims_list, this_reference_alignments, resolved_reference_tam, flat_read_groups[sam_file_index]);
      
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


void _write_reference_matches(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, const SequenceTrimsList& trims_list, alignment_list& reference_alignments, bam_file& reference_tam, const read_group_ref& rg)
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
	reference_tam.write_alignments(rg, reference_alignments, &trims_list, &ref_seq_info, true, true);
}

// Post-streaming DP merge-back. Held ambiguous discordant pairs (one unique mate + one IS mate that
// maps to several near-identical copies) are clustered by unique locus; per cluster the IS copy
// compatible with the MOST pairs wins the vote (ties -> lowest coordinate). Each pair is then written
// to the resolved reference SAM: ALL best-score copies of the redundant mate (redundant coverage
// preserved, X1>1) plus the unique mate, with the discordant-pair flags placed on the redundant
// alignment sitting in the chosen copy (or the pair's own lowest copy when it has none there).
void write_held_discordant_pairs(Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info,
                                 const SequenceTrimsList& trims_list, bam_file& resolved_reference_tam,
                                 vector<HeldDiscordantPair>& held_discordant_pairs)
{
  if (held_discordant_pairs.empty()) return;

  sort(held_discordant_pairs.begin(), held_discordant_pairs.end(),
       [](const HeldDiscordantPair& a, const HeldDiscordantPair& b) {
         if (a.unique_seq_id != b.unique_seq_id) return a.unique_seq_id < b.unique_seq_id;
         return a.unique_position < b.unique_position;
       });

  size_t i = 0;
  while (i < held_discordant_pairs.size())
  {
    // Cluster held pairs by unique locus (new cluster when the gap exceeds the running window).
    size_t j = i + 1;
    int32_t cluster_end = held_discordant_pairs[i].unique_position;
    double w = max(1.0, held_discordant_pairs[i].window);
    while (j < held_discordant_pairs.size()
           && held_discordant_pairs[j].unique_seq_id == held_discordant_pairs[i].unique_seq_id
           && (held_discordant_pairs[j].unique_position - cluster_end) <= w) {
      cluster_end = held_discordant_pairs[j].unique_position;
      w = max(w, held_discordant_pairs[j].window);
      j++;
    }

    // Vote: each pair contributes its DISTINCT candidate IS copies; the copy compatible with the most
    // pairs wins (ties -> lowest coordinate). Identical copies share a sequence, so the choice among
    // them does not affect the downstream MOB mob_region comparison.
    map<string,int> copy_votes;
    for (size_t k = i; k < j; k++) {
      set<string> pair_copies;
      for (alignment_list::iterator it = held_discordant_pairs[k].redundant_alignments.begin();
           it != held_discordant_pairs[k].redundant_alignments.end(); it++) {
        string c = dp_copy_id_for_alignment(ref_seq_info, it->get());
        if (!c.empty()) pair_copies.insert(c);
      }
      for (set<string>::iterator c = pair_copies.begin(); c != pair_copies.end(); c++)
        copy_votes[*c]++;
    }
    string chosen_copy; int best_votes = -1; int32_t best_start = 0;
    for (map<string,int>::iterator v = copy_votes.begin(); v != copy_votes.end(); v++) {
      size_t colon = v->first.rfind(':');
      size_t dash = (colon == string::npos) ? string::npos : v->first.find('-', colon);
      int32_t start = (colon != string::npos && dash != string::npos)
                      ? from_string<int32_t>(v->first.substr(colon + 1, dash - colon - 1)) : 0;
      if (v->second > best_votes || (v->second == best_votes && start < best_start)) {
        best_votes = v->second; best_start = start; chosen_copy = v->first;
      }
    }

    for (size_t k = i; k < j; k++) {
      HeldDiscordantPair& h = held_discordant_pairs[k];
      bam_alignment* ua = h.unique_alignments.front().get();
      bam_alignment* chosen = NULL;
      if (!chosen_copy.empty()) {
        for (alignment_list::iterator it = h.redundant_alignments.begin(); it != h.redundant_alignments.end(); it++) {
          if (dp_copy_id_for_alignment(ref_seq_info, it->get()) == chosen_copy) { chosen = it->get(); break; }
        }
      }
      if (chosen == NULL) chosen = pick_lowest_alignment(h.redundant_alignments);

      bool same_tid = (ua->reference_target_id() == chosen->reference_target_id());
      OrientationDistance od = same_tid ? compute_orientation_and_distance(ua, chosen) : OrientationDistance{"NA", 0};
      mark_pair_info(ua, chosen, same_tid, od.orientation, od.distance, /*is_concordant=*/false);

      _write_reference_matches(settings, summary, ref_seq_info, trims_list, h.unique_alignments, resolved_reference_tam, h.unique_rg);
      _write_reference_matches(settings, summary, ref_seq_info, trims_list, h.redundant_alignments, resolved_reference_tam, h.redundant_rg);
    }

    i = j;
  }
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
  
  // MEASUREMENT ONLY (--junction-debug). One row per (this candidate, read), recording the
  // read's score against THIS candidate next to the best-over-all-candidates difference the
  // pipeline currently uses. The two differ whenever a read aligns better to some OTHER candidate,
  // which is the case the current argmax mishandles: near-identical candidates are alternative
  // descriptions of one event, not competing origins. Writes nothing when the flag is off.
  counted_ptr<ofstream> mapq_debug(NULL);
  if (settings.junction_debug) {
    mapq_debug = counted_ptr<ofstream>(new ofstream(settings.junction_mapq_debug_file_name.c_str(), ios_base::app));
  }

	for (uint32_t i = 0; i < items.size(); i++) // READ (loops over unique_matches, degenerate_matches)
	{
		JunctionMatchPtr& item = items[i];

    if (verbose) cout << "  " << item->junction_alignments.front()->read_name() << endl;

    if (mapq_debug.get() != NULL) {
      // Deliberately a separate, non-asserting lookup of the alignment to this candidate: the
      // real one below runs after guards that would skip some of the rows we want to see.
      const alignment_wrapper* dbg_a = NULL;
      for (alignment_list::iterator it = item->junction_alignments.begin(); it != item->junction_alignments.end(); it++) {
        if ((*it)->reference_target_id() == junction_tid) { dbg_a = &(**it); break; }
      }
      if (dbg_a != NULL) {
        uint32_t dbg_as = 0, dbg_is_best = 0;
        dbg_a->aux_get_i(kBreseqAlignmentScoreBAMTag, dbg_as);
        dbg_a->aux_get_i(kBreseqBestAlignmentScoreBAMTag, dbg_is_best);
        int32_t delta_this = static_cast<int32_t>(dbg_as) - item->best_reference_score;
        // The weight column reports what the model actually uses: the quality-aware log-likelihood
        // ratio stamped on this alignment, not the score delta in the column beside it. They can
        // disagree, and that disagreement is the point of the table.
        //
        // "NA" when the read carries no competing likelihood, which is NOT the same as a ratio of
        // zero. Printing the missing case as 0 is how a bug that left the weighting inert on 8,831
        // of 13,029 reads read as "these reads are simply undecidable" for a full measurement cycle.
        uint32_t dbg_own_ll_raw = 0, dbg_other_ll_raw = 0;
        double dbg_log10_odds = 0.0;
        bool dbg_have_odds = dbg_a->aux_get_i(kBreseqOwnHypothesisLogLikelihoodBAMTag, dbg_own_ll_raw)
                          && dbg_a->aux_get_i(kBreseqOtherHypothesisLogLikelihoodBAMTag, dbg_other_ll_raw);
        if (dbg_have_odds) {
          dbg_log10_odds = static_cast<double>(static_cast<int32_t>(dbg_own_ll_raw) - static_cast<int32_t>(dbg_other_ll_raw))
                         / kAlignmentLogLikelihoodScale;
        }
        string dbg_odds_field = dbg_have_odds ? to_string(dbg_log10_odds) : string("NA");
        string dbg_weight_field = dbg_have_odds
                                ? to_string(junction_read_weight(dbg_log10_odds, dbg_a->redundancy(), 1.0))
                                : string("NA");
        *mapq_debug << junction_id
                    << "\t" << dbg_a->read_name()
                    << "\t" << delta_this                            // Delta against THIS candidate
                    << "\t" << item->mapping_quality_difference      // Delta against the BEST candidate
                    << "\t" << item->best_reference_score
                    << "\t" << dbg_as
                    << "\t" << dbg_is_best
                    << "\t" << item->degenerate_count
                    << "\t" << (flanking_left + 1 - static_cast<int32_t>(dbg_a->reference_start_1()))   // this_left
                    << "\t" << (static_cast<int32_t>(dbg_a->reference_end_1()) - flanking_left - abs(alignment_overlap)) // this_right
                    << "\t" << dbg_a->read_length()
                    << "\t" << dbg_a->redundancy()
                    << "\t" << dbg_odds_field
                    << "\t" << dbg_weight_field
                    << endl;
      }
    }

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
    
    // Measured over CONFIDENT bases only. this_left/this_right are how far the read reaches past
    // the breakpoint, and that reach IS the evidence that it spans the junction -- so a tail of
    // miscalled bases must not be allowed to supply it. The cutoff makes these accessors walk in
    // past any base at or below it, the same definition the RA caller and the alignment viewer use.
    uint32_t reference_start_1 = a->reference_start_1(settings.base_quality_cutoff);
    if (reference_start_1 == UNDEFINED_UINT32) continue;
		int32_t this_left = flanking_left + 1 - reference_start_1;

		// The right side starts after moving past any overlap (negative or positive)
    uint32_t reference_end_1 = a->reference_end_1(settings.base_quality_cutoff);
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
      
			const read_group_ref& rg = repeat_match.rg;
      
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
			// This is the best remaining (failing/marginal) junction for this read, so it
			// claims the read exclusively -- exactly mirroring the passing branch above.
			// The read is purged from every inferior candidate so that different marginal
			// junctions are not supported by the same reads. Because failing junctions are
			// resolved only after all passing junctions have already claimed their reads,
			// this never steals a read from a passing junction.
			else
			{
        // Purge all references to this read from the degenerate match hash
        // so that it cannot be counted or written for any inferior junction,
        // and collapse to the single alignment for THIS junction.
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
            repeat_junction_match_map[test_junction_seq_id].erase(a->read_name());
          }

          if (repeat_junction_match_map[test_junction_seq_id].size() == 0)
          {
            repeat_junction_match_map.erase(test_junction_seq_id);
          }
        }

        assert(matched_alignment.get() != NULL);
        repeat_match.junction_alignments.clear();
        repeat_match.junction_alignments.push_back(matched_alignment);

        // This read's best candidate junction is a rejected one, so its true home is the
        // reference genome. (Any read that also matched a passing junction would already
        // have been claimed by it in the earlier pass and removed from here, so a read
        // reaching this point failed on every real junction.) Write it back once.
        _write_reference_matches(settings, summary, ref_seq_info, trims_list, repeat_match.reference_alignments, resolved_reference_tam, rg);

        // Write the read's alignment to the (marginal) junction SAM file once.
        if (has_non_overlap_alignment) {
          alignment_list alignments;
          alignments.read_base_quality_char_string = repeat_match.junction_alignments.read_base_quality_char_string;
          alignments.read_base_quality_char_string_reversed = repeat_match.junction_alignments.read_base_quality_char_string_reversed;

          alignments.push_back(matched_alignment);
          resolved_junction_tam.write_alignments(rg, alignments, &junction_trims_list, &junction_ref_seq_info, true);
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
      const read_group_ref& rg = item.rg;
      
      // ONLY if we failed: write matches to reference sequences
      if (failed)
      {
        alignment_list this_reference_al = item.reference_alignments;
        _write_reference_matches(settings, summary, ref_seq_info, trims_list, this_reference_al, resolved_reference_tam, rg);
      }
      
      // REGARDLESS of success: write matches to the candidate junction SAM file
      resolved_junction_tam.write_alignments(rg, item.junction_alignments, &junction_trims_list, &junction_ref_seq_info, true);
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

	// Normalize each redundant side onto the lowest-coordinate CONSENSUS copy of its repeat family, so
	// a redundant JC side always reports the consensus IS. DP evidence carries the specific (possibly
	// variant) copy and can override the MOB's IS via mob_region during prediction. MOB prediction is
	// otherwise unaffected (it keys on the repeat family name and the unique side, not the copy).
	for (uint32_t i = 0; i < 2; i++) {
		if (jc.sides[i].redundant) {
			ref_seq_info.canonicalize_redundant_side_to_consensus_copy(jc.sides[i].seq_id, jc.sides[i].position, jc.sides[i].strand);
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

  // Weighted counterparts of the raw counts above. The raw fields stay integers -- they are
  // user-facing and cGenomeDiff::read_counts_for_entry() parses them with from_string<uint32_t> --
  // so the weighted values live alongside rather than replacing them.
  // Four sums, each staying in the window that produced it: the junction window's support for the
  // junction and for the reference, and each side window's support for the reference and for the
  // junction. Every window's two sums add to its own raw count.
  double junction_weighted = 0.0, junction_reference_weighted = 0.0, junction_weight_sq = 0.0;
  double side_weighted[2] = {0.0, 0.0};
  double side_weight_sq[2] = {0.0, 0.0};
  double side_junction_weighted[2] = {0.0, 0.0};
  // Per-read evidence, captured per window so the mixture below can re-weight in memory. The
  // reference counter is reused for both sides, so each side's list must be copied out before the
  // next count() overwrites it.
  vector<junction_read_counter::read_evidence_t> junction_evidence, side_evidence[2];
  if (junction_jrc.get() != NULL) {
    junction_weighted = junction_jrc->sum_weight();
    junction_reference_weighted = junction_jrc->sum_complement();
    junction_weight_sq = junction_jrc->sum_weight_sq();
    junction_evidence = junction_jrc->counted_read_evidence();
  }

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

    // This side's own weighted evidence. The junction reads' reference-hypothesis residual is
    // added later, once we know which sides are actually usable -- see the consolidation below.
    if (reference_jrc.get() != NULL) {
      side_weighted[0] = reference_jrc->sum_weight();
      side_weight_sq[0] = reference_jrc->sum_weight_sq();
      side_junction_weighted[0] = reference_jrc->sum_complement();
      side_evidence[0] = reference_jrc->counted_read_evidence();
    }
    
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

    if (reference_jrc.get() != NULL) {
      side_weighted[1] = reference_jrc->sum_weight();
      side_weight_sq[1] = reference_jrc->sum_weight_sq();
      side_junction_weighted[1] = reference_jrc->sum_complement();
      side_evidence[1] = reference_jrc->counted_read_evidence();
    }

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
    c = junction_weighted;
    c /= static_cast<double>(junction_possible_overlap_registers);
  }

  // @JEB: divide the side X counts by 2, if both were counted
  // or by 1 if that side of the alignment was ambiguous (for edges of IS-elements)
  double d = 2;
  if ((j[SIDE_1_READ_COUNT] == "NA") || (side_1_possible_overlap_registers == 0)) {
    a = 0; //"NA" in read count sets value to 1 not 0
    d--;
  } else {
    a = side_weighted[0];
    a /= static_cast<double>(side_1_possible_overlap_registers);
  }

  if ((j[SIDE_2_READ_COUNT] == "NA") || (side_2_possible_overlap_registers == 0)) {
    b = 0;
    d--;
  } else {
    b = side_weighted[1];
    b /= static_cast<double>(side_2_possible_overlap_registers);
  }

  // Fit the frequency and the per-read assignments together (EM on a two-component mixture).
  //
  // The E-step asks, for each read, "did this come from the junction?" using the current frequency
  // as the prior; the M-step recomputes the frequency from those soft assignments. Both steps run
  // in the rate space defined below, so the register normalization is inside the loop rather than
  // applied to a converged count. Reads are re-weighted from the stored evidence, not re-fetched.
  //
  // This is what removes the arbitrary constant: there is no prior to tune, because the only
  // self-consistent per-read prior is the frequency being estimated. Convergence is fast (the
  // mixture is one-dimensional and the likelihood is well behaved); the iteration cap is a
  // backstop, not a tuning knob.
  if (settings.junction_weight_reads && have_c && (d > 0)) {
    const uint32_t k_max_iterations = 50;
    const double k_tolerance = 1e-6;
    double f_fit = (a + b) / d;
    f_fit = (c + f_fit > 0.0) ? c / (c + f_fit) : 0.0;   // start from the unweighted estimate

    for (uint32_t iter = 0; iter < k_max_iterations; iter++) {
      double j_own, j_other, j_own_sq;
      reweight_window(junction_evidence, f_fit, j_own, j_other, j_own_sq, true);

      double s_own[2] = {0.0, 0.0}, s_other[2] = {0.0, 0.0}, s_own_sq[2] = {0.0, 0.0};
      for (uint32_t k = 0; k < 2; k++) {
        reweight_window(side_evidence[k], f_fit, s_own[k], s_other[k], s_own_sq[k], false);
      }

      // Rates, each in its own window's register basis (see the note below on why nothing moves
      // between windows).
      double junction_rate = j_own / static_cast<double>(junction_possible_overlap_registers);
      double reference_rate = j_other / static_cast<double>(junction_possible_overlap_registers);
      double side_ref = 0.0, side_junc = 0.0;
      if ((j[SIDE_1_READ_COUNT] != "NA") && (side_1_possible_overlap_registers != 0)) {
        side_ref  += s_own[0]   / static_cast<double>(side_1_possible_overlap_registers);
        side_junc += s_other[0] / static_cast<double>(side_1_possible_overlap_registers);
      }
      if ((j[SIDE_2_READ_COUNT] != "NA") && (side_2_possible_overlap_registers != 0)) {
        side_ref  += s_own[1]   / static_cast<double>(side_2_possible_overlap_registers);
        side_junc += s_other[1] / static_cast<double>(side_2_possible_overlap_registers);
      }
      junction_rate += side_junc / d;
      reference_rate += side_ref / d;

      double f_new = (junction_rate + reference_rate > 0.0)
                     ? junction_rate / (junction_rate + reference_rate) : 0.0;
      bool converged = (fabs(f_new - f_fit) < k_tolerance);
      f_fit = f_new;

      if (converged || (iter + 1 == k_max_iterations)) {
        junction_weighted = j_own;
        junction_reference_weighted = j_other;
        junction_weight_sq = j_own_sq;
        for (uint32_t k = 0; k < 2; k++) {
          side_weighted[k] = s_own[k];
          side_junction_weighted[k] = s_other[k];
          side_weight_sq[k] = s_own_sq[k];
        }
        // Recompute the per-window rates from the converged assignments.
        c = junction_weighted / static_cast<double>(junction_possible_overlap_registers);
        if ((j[SIDE_1_READ_COUNT] != "NA") && (side_1_possible_overlap_registers != 0)) {
          a = side_weighted[0] / static_cast<double>(side_1_possible_overlap_registers);
        }
        if ((j[SIDE_2_READ_COUNT] != "NA") && (side_2_possible_overlap_registers != 0)) {
          b = side_weighted[1] / static_cast<double>(side_2_possible_overlap_registers);
        }
        j["junction_mixture_iterations"] = to_string(iter + 1);
        break;
      }
    }
  }

  // Record the split, and check that each window conserved its own reads.
  //
  // Conservation is PER WINDOW, not across windows. An earlier version moved a junction read's
  // reference-hypothesis share over to a side window, which is wrong twice over: each window is
  // normalized by its own register count, so a count moved across changes basis; and the sides
  // enter the frequency as a mean, so weight w removed from the junction reappears as only w/2.
  // Keeping every read's whole 1.0 inside the window that observed it makes the invariant exact
  // and lets the frequency combine the pieces in rate space, where they are comparable.
  bool side_usable[2] = { (j[SIDE_1_READ_COUNT] != "NA"), (j[SIDE_2_READ_COUNT] != "NA") };
  {
    j["new_junction_weighted_read_count"] = (junction_jrc.get() != NULL) ? to_string(junction_weighted, 2) : "0";
    // The junction window's own reference evidence: reads that span the breakpoint but do not
    // discriminate. It belongs in the frequency's denominator at the JUNCTION register basis,
    // which is why it is reported separately rather than folded into a side.
    j["new_junction_reference_weighted_read_count"] = (junction_jrc.get() != NULL) ? to_string(junction_reference_weighted, 2) : "0";
    j["side_1_weighted_read_count"] = side_usable[0] ? to_string(side_weighted[0], 2) : "NA";
    j["side_2_weighted_read_count"] = side_usable[1] ? to_string(side_weighted[1], 2) : "NA";

    double junction_raw = from_string<double>(j[NEW_JUNCTION_READ_COUNT]);
    ASSERT(fabs((junction_weighted + junction_reference_weighted) - junction_raw) < 1e-6 * (1.0 + junction_raw),
           "Junction window weighting did not conserve reads for " + j["key"] + ": raw "
           + to_string(junction_raw, 6) + " vs split "
           + to_string(junction_weighted + junction_reference_weighted, 6));

    for (uint32_t k = 0; k < 2; k++) {
      if (!side_usable[k]) continue;
      double side_raw = from_string<double>(j[(k == 0) ? SIDE_1_READ_COUNT : SIDE_2_READ_COUNT]);
      ASSERT(fabs((side_weighted[k] + side_junction_weighted[k]) - side_raw) < 1e-6 * (1.0 + side_raw),
             "Side " + to_string(k + 1) + " window weighting did not conserve reads for " + j["key"]
             + ": raw " + to_string(side_raw, 6) + " vs split "
             + to_string(side_weighted[k] + side_junction_weighted[k], 6));
    }
  }

  // Kish effective sample size of the weighted evidence, (sum w)^2 / (sum w^2), pooled over the
  // junction and the usable sides with the same side-averaging the frequency uses. Fractional
  // weights are not binomial trials, so this is the honest "n" for a confidence bound built on
  // them: it equals the raw count when every weight is 1 and falls as the weights spread out.
  // Recorded now; the bounds start using it in a later change.
  {
    double n_sides_for_eff = 0.0, side_w = 0.0, side_w_sq = 0.0;
    if (j[SIDE_1_READ_COUNT] != "NA") { side_w += side_weighted[0]; side_w_sq += side_weight_sq[0]; n_sides_for_eff++; }
    if (j[SIDE_2_READ_COUNT] != "NA") { side_w += side_weighted[1]; side_w_sq += side_weight_sq[1]; n_sides_for_eff++; }
    double total_w = junction_weighted + ((n_sides_for_eff > 0.0) ? side_w / n_sides_for_eff : 0.0);
    double total_w_sq = junction_weight_sq + ((n_sides_for_eff > 0.0) ? side_w_sq / n_sides_for_eff : 0.0);
    j["junction_effective_depth"] = (total_w_sq > 0.0) ? to_string(total_w * total_w / total_w_sq, 2) : "0";
  }

  // The cross terms: each window sees evidence for BOTH hypotheses, so the junction window
  // contributes to the reference rate and the side windows contribute to the junction rate. Every
  // term is divided by its own window's register count, because that is the number of chances that
  // window had to observe the read -- this is why the split is not done by moving counts around.
  //
  // The read sets are disjoint (a read counted at the junction is excluded from the sides by
  // name), so these are complementary observation channels and their rates ADD. The two side
  // windows are instead two measurements of the same local reference coverage, so they are
  // averaged -- which is exactly what the existing (a+b)/d already did.
  //
  // With weighting off every cross term is zero and this collapses to the original
  // c / (c + (a+b)/d) identically; that exact reduction is what makes the change safe to land
  // before the prior has been calibrated.
  double c_from_sides = 0.0, r_from_sides = 0.0;
  if (d > 0) {
    if ((j[SIDE_1_READ_COUNT] != "NA") && (side_1_possible_overlap_registers != 0)) {
      c_from_sides += side_junction_weighted[0] / static_cast<double>(side_1_possible_overlap_registers);
    }
    if ((j[SIDE_2_READ_COUNT] != "NA") && (side_2_possible_overlap_registers != 0)) {
      c_from_sides += side_junction_weighted[1] / static_cast<double>(side_2_possible_overlap_registers);
    }
    c_from_sides /= d;
  }
  if (have_c) {
    r_from_sides = junction_reference_weighted / static_cast<double>(junction_possible_overlap_registers);
  }

  // We cannot assign a frequency if the denominator is zero
  if ((d == 0) || !have_c) {
    // "NA" rather than omission, matching the sentinel this entry already uses for an undefined
    // side count. A frequency that could not be assigned has no interval either.
    j[FREQUENCY] = "NA";
    j[FREQUENCY_LOWER] = "NA";
    j[FREQUENCY_UPPER] = "NA";

    j[PREDICTION] = "unknown";
  } else {
    double junction_rate = c + c_from_sides;
    double reference_rate = ((a + b) / d) + r_from_sides;
    double new_junction_frequency_value = (junction_rate + reference_rate > 0.0)
                                          ? junction_rate / (junction_rate + reference_rate) : 0.0;
    j[FREQUENCY] = to_string(new_junction_frequency_value, settings.polymorphism_precision_places, true);

    // Determine what kind of prediction we are
    
    // We may have added FREQUENCY_CUTOFF previously, so clear it.
    j.remove_reject_reason("FREQUENCY_CUTOFF");
    
    // Exact (Clopper-Pearson) confidence bounds on the junction frequency, so the cutoffs below are
    // confidence statements rather than point-estimate comparisons -- a junction is never rejected
    // merely for having few reads, only for being confidently on the wrong side of a cutoff.
    //
    // The interval is CENTRED on the frequency this entry actually reports, and given the width
    // implied by an effective depth. This mirrors RA_frequency_bounds() (identify_mutations.cpp),
    // and for the same reason: a bound has to describe the number shown, or "rejected at frequency
    // 0.43" carrying an interval computed from 0.21 is unreadable.
    //
    // It used to use the RAW integer counts instead, on the grounds that a binomial needs genuine
    // trials, with the observation that raw and normalized frequencies then agreed to a median of
    // 0.0000 across the goldens. That is true for long reads at high depth and false where it
    // matters: on 35 bp polymorphism data the two differ by a median of 0.015 and up to 0.16, and
    // once reads are weighted by how well they discriminate the gap grows to 0.83. The "genuine
    // trials" property is restored instead by the effective depth below, which equals the read
    // count when every read is decisive and shrinks as the weights spread out.
    double freq_effective_depth = from_string<double>(j["junction_effective_depth"]);
    double freq_successes = new_junction_frequency_value * freq_effective_depth;
    double freq_lower = binomial_frequency_lower_bound(freq_successes, freq_effective_depth);
    double freq_upper = binomial_frequency_upper_bound(freq_successes, freq_effective_depth);

    // Recorded, not just used and discarded: the interval that decides the call below is then the
    // one the .gd and the report show. Same formatting as the frequency it brackets.
    j[FREQUENCY_LOWER] = to_string(freq_lower, settings.polymorphism_precision_places, true);
    j[FREQUENCY_UPPER] = to_string(freq_upper, settings.polymorphism_precision_places, true);

    // Two bounds decide all three outcomes. Note the asymmetry is deliberate: the LOWER bound gates
    // "is there a junction here at all", while the UPPER bound gates "is it fixed rather than
    // polymorphic". Testing the high side with the lower bound instead would require
    // ln(alpha)/ln(cutoff) reads even at 100% frequency -- 14 at a cutoff of 0.8 -- and would drop
    // real junctions for want of depth.
    bool below_floor = (settings.polymorphism_frequency_cutoff > 0.0)
                       && (freq_lower < settings.polymorphism_frequency_cutoff);
    bool is_consensus = (settings.consensus_frequency_cutoff > 0.0)
                        && (freq_upper >= settings.consensus_frequency_cutoff);

    if (below_floor) {
      // Not confidently a junction at all. Rejected in both modes.
      j[PREDICTION] = "polymorphism";
      j.add_reject_reason("FREQUENCY_CUTOFF");
    } else if (is_consensus) {
      // [frequency] keeps the fitted estimate; [prediction] carries the snap to 1, which
      // mutation_predictor re-applies via cDiffEntry::mutation_frequency() when it builds the
      // mutation. Overwriting the frequency here would leave the bounds bracketing a number the
      // entry no longer reports.
      j[PREDICTION] = "consensus";
    } else {
      j[PREDICTION] = "polymorphism";
    }

    // DIAGNOSTIC, not used for the call: what the old raw-count interval would have decided.
    // Kept so a call that moves can be attributed to the re-centring rather than to the weighting.
    if (settings.junction_debug) {
      double raw_new = from_string<double>(j[NEW_JUNCTION_READ_COUNT]);
      double raw_sides = 0.0; double raw_side_n = 0.0;
      if (j[SIDE_1_READ_COUNT] != "NA") { raw_sides += from_string<double>(j[SIDE_1_READ_COUNT]); raw_side_n++; }
      if (j[SIDE_2_READ_COUNT] != "NA") { raw_sides += from_string<double>(j[SIDE_2_READ_COUNT]); raw_side_n++; }
      double raw_total = raw_new + ((raw_side_n > 0.0) ? (raw_sides / raw_side_n) : 0.0);
      double raw_lower = binomial_frequency_lower_bound(raw_new, raw_total);
      double raw_upper = binomial_frequency_upper_bound(raw_new, raw_total);
      string old_call = "polymorphism";
      if ((settings.polymorphism_frequency_cutoff > 0.0) && (raw_lower < settings.polymorphism_frequency_cutoff)) {
        old_call = "rejected";
      } else if ((settings.consensus_frequency_cutoff > 0.0) && (raw_upper >= settings.consensus_frequency_cutoff)) {
        old_call = "consensus";
      }
      j["junction_rawbounds_prediction"] = old_call;
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
  
  
  
// Build the pair of read counters used to tally junction evidence, configured identically for
// every caller. Either may come back NULL when its BAM is absent, which callers already handle.
//
// This exists because there is more than one place that counts junction reads, and they must agree
// on how a read is split between the junction and the reference and on which bases are trusted.
// When they disagreed, the later caller silently overwrote the earlier one's results -- and no
// assertion fired, because an unweighted count is internally consistent.
void make_junction_read_counters(const Settings& settings,
                                 counted_ptr<junction_read_counter>& reference_jrc,
                                 counted_ptr<junction_read_counter>& junction_jrc)
{
  reference_jrc = counted_ptr<junction_read_counter>(NULL);
  junction_jrc = counted_ptr<junction_read_counter>(NULL);

  if (file_exists(settings.reference_bam_file_name.c_str())) {
    reference_jrc = counted_ptr<junction_read_counter>(new junction_read_counter(settings.reference_bam_file_name, settings.reference_fasta_file_name, settings.verbose, settings.reference_trim_file_name));
  }
  if (file_exists(settings.junction_bam_file_name.c_str()) && file_exists(settings.candidate_junction_fasta_file_name.c_str())) {
    junction_jrc = counted_ptr<junction_read_counter>(new junction_read_counter(settings.junction_bam_file_name, settings.candidate_junction_fasta_file_name, settings.verbose, settings.candidate_junction_trim_file_name));
  }

  // Opposite delta signs: in the junction BAM a record's own score belongs to the junction
  // hypothesis, in the reference BAM it belongs to the reference. See set_weighting().
  if (reference_jrc.get() != NULL) {
    reference_jrc->set_weighting(settings.junction_weight_reads, -1);
    reference_jrc->set_base_quality_cutoff(settings.base_quality_cutoff);
  }
  if (junction_jrc.get() != NULL) {
    junction_jrc->set_weighting(settings.junction_weight_reads, +1);
    junction_jrc->set_base_quality_cutoff(settings.base_quality_cutoff);
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
  
  counted_ptr<junction_read_counter> reference_jrc, junction_jrc;
  make_junction_read_counters(settings, reference_jrc, junction_jrc);
  
  for (diff_entry_list_t::iterator it = jc.begin(); it != jc.end(); it++) {
    cDiffEntry& j = **it;
    assign_one_junction_read_counts(settings, summary, j, reference_jrc, junction_jrc);
  }
  
}

  
junction_read_counter::junction_read_counter(const string& bam, const string& fasta, bool verbose, const string& trim_file_name)
  : pileup_base(bam, fasta), _verbose(verbose)
  , _weighting_enabled(false)  // off until set_weighting() -- every read counts 1
  , _base_quality_cutoff(0)   // raw bounds until set_base_quality_cutoff()
  , _delta_sign(1)
  , _sum_weight(0.0), _sum_weight_sq(0.0), _sum_complement(0.0)
{

  _trims_list.resize(num_targets());
  for (uint32_t tid = 0; tid < num_targets(); tid++) {
    string this_file_name = Settings::file_name(trim_file_name, "@", target_name(tid));
    _trims_list[tid].ReadFile(this_file_name, target_length(tid));
  }
}

void junction_read_counter::set_weighting(bool enabled, int32_t delta_sign)
{
  _weighting_enabled = enabled;
  _delta_sign = delta_sign;
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
  _sum_weight = 0.0;
  _sum_weight_sq = 0.0;
  _sum_complement = 0.0;
  _read_evidence.clear();
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
  // Confident bounds: a read must reach the window with real bases, not a miscalled tail.
  q_start = a.reference_start_1(_base_quality_cutoff);
  q_end   = a.reference_end_1(_base_quality_cutoff);
  if ((q_start == UNDEFINED_UINT32) || (q_end == UNDEFINED_UINT32)) {
    if (_verbose) cout << "  NO CONFIDENT BASES" << endl;
    return;
  }

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

  // Split this read between the two hypotheses. w is the posterior that it came from the junction;
  // 1-w is the matching evidence for the reference. Every counted read contributes exactly 1.0 in
  // total, so the frequency's numerator and denominator stay on the same footing -- weighting one
  // side only would shrink the denominator and inflate every frequency.
  //
  // A read missing either log-likelihood tag has no competing hypothesis on record (it aligned to
  // only one of the two references), so it falls back to w = 1 toward the BAM it sits in.
  //
  // X8/X7, not AS: AS is breseq's rescored ALIGNMENT SCORE, which is quality-blind -- every bowtie2
  // pass runs --ignore-quals and alignment_score() reads only the CIGAR and the reference, so a
  // mismatch at Q40 and one at Q2 move it identically. The X8/X7 pair carries the quality-aware
  // log-likelihood of each hypothesis instead (alignment_log10_likelihood(), stamped during
  // resolution), which is what makes the ratio below a real likelihood ratio. Routing still uses
  // AS; only this weighting reads the new tags.
  //
  // Watch the tag names: an earlier version read X5 here, which is set in memory during resolution
  // but never serialized, so it silently found nothing and fell back to weight 1 -- and an
  // unweighted read conserves perfectly, so no assertion could catch it.
  //
  // w is the weight toward THIS window's own hypothesis, so the fallback of 1.0 means "count this
  // read exactly as before" whichever BAM we are reading. Deriving it by flipping a junction-side
  // weight would break that: with weighting disabled the junction weight is 1.0, and flipping it
  // would zero every reference side and force the frequency to 1.
  double w = 1.0;
  double read_log10_odds = 0.0;
  bool have_odds = false;
  uint32_t own_ll_raw, other_ll_raw;
  if (a.aux_get_i(kBreseqOwnHypothesisLogLikelihoodBAMTag, own_ll_raw)
      && a.aux_get_i(kBreseqOtherHypothesisLogLikelihoodBAMTag, other_ll_raw)) {
    // aux_get_i hands back a uint32_t; these tags are signed (a log-likelihood is <= 0), so the
    // round trip through int32_t is what recovers the value.
    int32_t own_ll = static_cast<int32_t>(own_ll_raw);
    int32_t other_ll = static_cast<int32_t>(other_ll_raw);
    read_log10_odds = _delta_sign * static_cast<double>(own_ll - other_ll) / kAlignmentLogLikelihoodScale;
    have_odds = true;
    if (_weighting_enabled) {
      double w_junction = junction_read_weight(read_log10_odds, a.redundancy(), 1.0);
      w = (_delta_sign > 0) ? w_junction : (1.0 - w_junction);
    }
  }

  // Both halves stay in THIS window. A read's evidence is not moved to another window, because
  // each window is normalized by its own number of confident overlap registers and a count moved
  // across would silently change register basis. The caller combines the four sums in rate space.
  _sum_weight += w;
  _sum_weight_sq += w * w;
  _sum_complement += 1.0 - w;

  // Keep the raw evidence so the caller can re-weight this read at a different prior without
  // fetching it again.
  read_evidence_t ev;
  ev.log10_odds_for_junction = read_log10_odds;
  ev.have_odds = have_odds;                         // false = no competing alignment on record
  ev.n_ref_placements = a.redundancy();
  _read_evidence.push_back(ev);

  // record that we counted this read
  _counted_read_names[a.read_name()] = true;
}
  
  
} // namespace breseq

