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

#include "libbreseq/resolve_alignments.h"

#include "libbreseq/genome_diff.h"
#include "libbreseq/fastq.h"
#include "libbreseq/fasta.h"
#include "libbreseq/alignment.h"
#include "libbreseq/annotated_sequence.h"
#include "libbreseq/chisquare.h"

using namespace std;

namespace breseq {
    
  
double combination (int32_t num, int32_t choose)
{
  double log_result = 0.0;
  for (int32_t i=num; i > choose; i--)
  {
    log_result += log(i);
  }
  for (int32_t i=2; i <= num - choose; i++)
  {
    log_result -= log(i);
  }
 
  return exp(log_result);
}
  
double binomial (double pr_success, int32_t num_trials, int32_t num_successes)
{
  double ret_val = combination(num_trials, num_successes) * pow(pr_success, num_successes) * pow(1-pr_success, num_trials-num_successes);
  return ret_val;
}
  
uint32_t qbinomial (double tail_value, int32_t num_trials, double pr_success)
{
  ASSERT((tail_value >= 0) && (tail_value < 1), "probability out of range");
  double cumulative_pr = 0.0;

  int32_t num_successes;
  for (num_successes=0; num_successes < num_trials; num_successes++) {
    cumulative_pr += binomial(pr_success, num_trials, num_successes);
    if (cumulative_pr > tail_value) break;
  }
  
  return num_successes;
}
  
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
  
PosHashProbabilityTable::PosHashProbabilityTable(Summary& summary)
{
  average_read_length = summary.sequence_conversion.avg_read_length;
  for (map<string,Coverage>::iterator it=summary.preprocess_coverage.begin();
       it != summary.preprocess_coverage.end(); it++) {
    
    string seq_id = it->first;
    Coverage& cov = it->second;
    
    Parameters p = {
          cov.nbinom_size_parameter,
          cov.nbinom_prob_parameter,
          summary.error_count[seq_id].no_pos_hash_per_position_pr,
          cov.nbinom_mean_parameter
    };
    
    param[it->first] = p;
  }
}

double PosHashProbabilityTable::probability(string& seq_id, uint32_t pos_hash_score, uint32_t overlap)
{
  bool verbose = false;
  if (
      probability_table.count(seq_id) 
      && probability_table[seq_id].count(pos_hash_score) 
      && probability_table[seq_id][pos_hash_score].count(overlap)
      )
    return probability_table[seq_id][pos_hash_score][overlap];
  
  Parameters& p = param[seq_id];

  // no coverage was fit for this fragment.
  if (p.negative_binomial_size == 0) return 999999;
  
  uint32_t junction_max_score = 2 * average_read_length;
  uint32_t max_coverage = 10*p.average_coverage;
        
  double pr = 0;
  
  if (verbose)
    cout << "Calculating: seq_id " << seq_id << " pos_hash_score " << pos_hash_score << " overlap " << overlap << endl;
      
  for (uint32_t this_coverage=1; this_coverage<= max_coverage; this_coverage++) {
    
    double this_cov_pr =  nbdtr(this_coverage, p.negative_binomial_size, p.negative_binomial_prob)
                        - nbdtr(this_coverage-1, p.negative_binomial_size, p.negative_binomial_prob); 
    
    if (verbose)
      cout << "NB: " << this_coverage << " " << p.negative_binomial_size << " " << p.negative_binomial_prob << endl;

    double this_chance_per_pos_strand_read_start = 1 - pow(p.chance_per_pos_strand_no_read_start, (this_coverage / static_cast<double>(p.average_coverage)));

    double this_pos_hash_pr = 0;
    for (uint32_t i=0; i <= pos_hash_score; i++) {

      //chance of getting pos_hash_score or lower
      this_pos_hash_pr += binomial(this_chance_per_pos_strand_read_start, junction_max_score - 2 * overlap, i);
    }
    
    if (verbose)
      cout << this_coverage << " " << this_cov_pr << " " << this_pos_hash_pr << " " << this_chance_per_pos_strand_read_start << endl;
    
    double this_pr = this_cov_pr*this_pos_hash_pr;
    pr += this_pr;
    
    if (verbose)
      cout << "  " << pr << endl;      
  }
  double log_pr = -log(pr)/log(10);
  probability_table[seq_id][pos_hash_score][overlap] = log_pr;
  return log_pr;
}

  
// Compares matches to candidate junctions with matches to original genome
void resolve_alignments(
                        const Settings& settings,
                        Summary& summary,
                        cReferenceSequences& ref_seq_info,
                        bool junction_prediction,
                        cReadFiles& read_files
                        ) 
{    
	bool verbose = false;

  // @JEB eventually deprecate this
  map<string,pos_hash_p_value_table> p_value_table_map;
  if (settings.use_r_junction_p_value_table_check) {      

    for(vector<cAnnotatedSequence>::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); it++) {
      const string& seq_id = it->m_seq_id;
      string coverage_junction_pos_hash_p_value_file_name = Settings::file_name(settings.coverage_junction_pos_hash_p_value_file_name, 
                                                                                "@", seq_id);
      pos_hash_p_value_table p_value_table(coverage_junction_pos_hash_p_value_file_name);
      p_value_table_map[seq_id] = p_value_table;
    }
  }
  
  //@JEB removed
  //calculate_cutoffs(settings, summary, ref_seq_info);
  
  // local variables for convenience
  map<string,int32_t>& distance_cutoffs = summary.alignment_resolution.distance_cutoffs;
  storable_map<string, storable_vector<int32_t> >& pos_hash_cutoffs = summary.alignment_resolution.pos_hash_cutoffs;
  int32_t avg_read_length = round(summary.sequence_conversion.avg_read_length);

  
  // Load the reference sequence trims, for writing resolved alignments
  SequenceTrimsList trims_list;
  read_trims(trims_list, ref_seq_info, settings.reference_trim_file_name);
  
  // Create junction trims directly from FASTA
  SequenceTrimsList junction_trims_list;
  
  if (junction_prediction)
  {  
    cFastaFile ff(settings.candidate_junction_fasta_file_name, ios::in);
    cFastaSequence fs;
    while (ff.read_sequence(fs))
    {
      SequenceTrims st(fs.m_sequence);
      junction_trims_list.push_back(st);
    }
  }
    
	// ####
	// ##	Junction sequences
	// ####

	//## if there were no candidate junctions (file is empty) then we seg fault if we try to use samtools on it...
  cReferenceSequences junction_ref_seq_info;
	if (junction_prediction
			&& !file_exists(settings.candidate_junction_fasta_file_name.c_str())
			&& !file_empty(settings.candidate_junction_fasta_file_name.c_str())
		)
  {
		junction_prediction = 0;
  }
    
	vector<JunctionInfo> junction_info_list;
    

  //## Preload all of the information about junctions
  //## so that we only have to split the names once
  
  if (junction_prediction) {
    junction_ref_seq_info.ReadFASTA(settings.candidate_junction_fasta_file_name);
		string junction_sam_file_name = settings.file_name(settings.candidate_junction_sam_file_name, "#", read_files[0].m_base_name);
		tam_file junction_tam(junction_sam_file_name, settings.candidate_junction_fasta_file_name, ios::in);

    junction_info_list.resize(junction_tam.bam_header->n_targets);
    
		for (int i = 0; i < junction_tam.bam_header->n_targets; i++) {
			junction_info_list[i] = JunctionInfo(junction_tam.bam_header->target_name[i]);
		}
    if (verbose) cout << "Number of candidate junctions: " << junction_info_list.size() << endl;
	}

	//####
	//##	Output files
	//####

	genome_diff gd;
    
  tam_file resolved_reference_tam(settings.resolved_reference_sam_file_name, settings.reference_fasta_file_name, ios::out);
  tam_file resolved_junction_tam(settings.resolved_junction_sam_file_name, settings.candidate_junction_fasta_file_name, ios::out);
    
  UniqueJunctionMatchMap unique_junction_match_map;    // Map of junction_id to MatchedJunction
	RepeatJunctionMatchMap repeat_junction_match_map;  // Map of junction_id to read_name to MatchedJunction

  // stores all junction ids that we have encountered
  map<string,uint32_t> all_junction_ids;
  
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
  
  if (verbose)
  {
    cout << "Total junction ids: " << all_junction_ids.size() << endl;
  }
  
	///
	// Determine which junctions are real, prefer ones with most matches
	///
  
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
                   junction_info_list
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

    if (test_info.total_non_overlap_reads > 0) {
      junction_test_info_list.push_back(test_info);
    }
    else {
      resolve_junction(
                       settings,
                       summary,
                       ref_seq_info,
                       trims_list,
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
  
  // @JEB Better algorithm:
  //  Score every junction and sort
  //  Accept top match...
  //  (Optional: rescore every junction, for lost degenerate matches, and sort)
  //  Accept next top match...
  // To implement this, we need to add a separate score function.

	//sort junction ids based on size of vector
  junction_test_info_list.sort();
  
  PosHashProbabilityTable pos_hash_p_value_calculator(summary);
  
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
                   junction_info_list
                   );
    
    
    JunctionInfo junction_info(junction_id);
    
    // Test the best-scoring junction.
    bool failed = false;  

    int32_t possible_overlap_positions = avg_read_length - 1 - abs(junction_info.alignment_overlap);
    ASSERT(possible_overlap_positions > 0, "Possible overlap positions <= 0");
    
    /* @JEB Deprecate...
    uint32_t pos_hash_score_cutoff_1, pos_hash_score_cutoff_2;
    int32_t distance_score_cutoff_1, distance_score_cutoff_2;
    
      pos_hash_score_cutoff_1 = pos_hash_cutoffs[junction_info.sides[0].seq_id][possible_overlap_positions];
      pos_hash_score_cutoff_2 = pos_hash_cutoffs[junction_info.sides[1].seq_id][possible_overlap_positions];
      
      distance_score_cutoff_1 = possible_overlap_positions - distance_cutoffs[junction_info.sides[0].seq_id];
      distance_score_cutoff_2 = possible_overlap_positions - distance_cutoffs[junction_info.sides[1].seq_id];
    */
    
    // table is by overlap, then pos hash score
    double neg_log10_p_value_1 = 999999999.99;
    double neg_log10_p_value_2 = 999999999.99;    
    
    double neg_log10_p_value_1_2 = pos_hash_p_value_calculator.probability(junction_info.sides[0].seq_id, junction_test_info.pos_hash_score, abs(junction_info.alignment_overlap));
    double neg_log10_p_value_2_2 = pos_hash_p_value_calculator.probability(junction_info.sides[1].seq_id, junction_test_info.pos_hash_score, abs(junction_info.alignment_overlap));

    //@JEB Eventually deprecate this. It's just a check.
    if (settings.use_r_junction_p_value_table_check) {      
      if (p_value_table_map[junction_info.sides[0].seq_id].size())
          neg_log10_p_value_1 = p_value_table_map[junction_info.sides[0].seq_id][abs(junction_info.alignment_overlap)][junction_test_info.pos_hash_score];
      if (p_value_table_map[junction_info.sides[1].seq_id].size())
        neg_log10_p_value_2 = p_value_table_map[junction_info.sides[1].seq_id][abs(junction_info.alignment_overlap)][junction_test_info.pos_hash_score];
      ASSERT(neg_log10_p_value_1 - neg_log10_p_value_1_2 < 0.01, "Variance in calculated p-values: " + to_string(neg_log10_p_value_1) + " " + to_string(neg_log10_p_value_1_2) );
      ASSERT(neg_log10_p_value_2 - neg_log10_p_value_2_2 < 0.01, "Variance in calculated p-values: " + to_string(neg_log10_p_value_2) + " " + to_string(neg_log10_p_value_2_2) );
    }
    
    neg_log10_p_value_1 = neg_log10_p_value_1_2;
    neg_log10_p_value_2 = neg_log10_p_value_2_2;
    
    // take the most significantly below pos_hash cutoff
    junction_test_info.neg_log10_pos_hash_p_value = max(neg_log10_p_value_1, neg_log10_p_value_2);
    
    // Both pos_hash score cutoffs might be zero - indicating these are missing contigs 
    // that are basically deleted. Always fail in this case.
    //failed = failed || (pos_hash_score_cutoff_1 == 0); 
    //failed = failed || (pos_hash_score_cutoff_2 == 0);
    
    failed = failed || (junction_test_info.neg_log10_pos_hash_p_value > settings.junction_pos_hash_neg_log10_p_value_cutoff);
    
    //failed = failed || ( junction_test_info.pos_hash_score < pos_hash_score_cutoff_1 );
    //failed = failed || ( junction_test_info.pos_hash_score < pos_hash_score_cutoff_2 );
    
    ///failed = failed || ( junction_test_info.max_left < distance_score_cutoff_1 );
    //failed = failed || ( junction_test_info.max_right < distance_score_cutoff_2 );
        
    if (verbose) 
    {
      cout << "Testing Junction: " << junction_id << endl;
      cout <<  (failed ? "FAILED" : "SUCCESS") << endl;
      cout << "  Neg log10 pos hash p-value: " << junction_test_info.neg_log10_pos_hash_p_value
           << " [ " << settings.junction_pos_hash_neg_log10_p_value_cutoff << " ]" << endl;
      //cout << "  Pos hash score: " << junction_test_info.pos_hash_score << " [ " << pos_hash_score_cutoff_1 << " / " << pos_hash_score_cutoff_2 << " ]" << endl;
      //cout << "  Distance score: " << junction_test_info.max_left << " | " << junction_test_info.max_right 
      //     << " [ " << distance_score_cutoff_1 << " / " << distance_score_cutoff_1 << " ]" << endl;
      cout << "  Number of unique matches: " << unique_junction_match_map[junction_id].size() << endl;
      size_t num_degenerate_matches = repeat_junction_match_map.count(junction_id) ? repeat_junction_match_map[junction_id].size() : 0;
      cout << "  Number of degenerate matches: " << num_degenerate_matches << endl;
      cout << "  Number of total_non_overlap reads: " << junction_test_info.total_non_overlap_reads  << endl;
    }
    
    //@JEB: still need to resolve ones with no overlap reads to give their repeat matches to the reference.
    //ASSERT(failed || (junction_test_info.total_non_overlap_reads > 0), "Junction passed with no non-overlap reads.");

    
    resolve_junction(
                     settings,
                     summary,
                     ref_seq_info,
                     trims_list,
                     junction_test_info.junction_id,
                     unique_junction_match_map,
                     repeat_junction_match_map,
                     resolved_reference_tam,
                     resolved_junction_tam,
                     failed,
                     junction_test_info.total_non_overlap_reads > 0
                     ); 
    
    if (!failed) 
      passed_junction_test_info_list.push_back(junction_test_info);
    else
      rejected_junction_test_info_list.push_back(junction_test_info);
    
    junction_test_info_list.pop_back();
    
    // @JEB TODO: Re-score ones that might have changed due to removing repeat matches and re-sort
    // to do this efficiently, we need a list of their junction id's to be passed back by resolve_junction
    //junction_test_info_list.sort();
  }
    
  map<int32_t, int32_t> accepted_pos_hash_score_distribution;
	map<int32_t, int32_t> observed_pos_hash_score_distribution;
  for(list<JunctionTestInfo>::iterator it = passed_junction_test_info_list.begin(); it != passed_junction_test_info_list.end(); it++)
  {
    JunctionTestInfo& junction_test_info = *it;
		string key = junction_test_info.junction_id;
		if (verbose) cout << key << endl;
		diff_entry item = junction_to_diff_entry(key, ref_seq_info, junction_test_info);
		gd.add(item);

		// save the score in the distribution
		add_score_to_distribution(accepted_pos_hash_score_distribution, junction_test_info.pos_hash_score);

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

				// Do not count for coverage if it is redundant!!
				if (from_string<int32_t>(item[side_key + "_redundant"])) continue;
        
				// Write out match corresponding to this part to SAM file
				// By trimming in the candidate junctions sequence, rather than on each half,
				// this is done properly.
				Trims trims = _trim_ambiguous_ends(a, junction_trims_list);
        
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
					&trims
				);
			}
		}
	}

  for(list<JunctionTestInfo>::iterator it = rejected_junction_test_info_list.begin(); it != rejected_junction_test_info_list.end(); it++)
  {
    JunctionTestInfo& junction_test_info = *it;
		string key = junction_test_info.junction_id;
		diff_entry item = junction_to_diff_entry(key, ref_seq_info, junction_test_info);
		add_reject_reason(item, "NJ");
		gd.add(item);
	}

  // Save summary statistics
	summary.alignment_resolution.observed_pos_hash_score_distribution = observed_pos_hash_score_distribution;
	summary.alignment_resolution.accepted_pos_hash_score_distribution = accepted_pos_hash_score_distribution;

  // Write the genome diff file
	gd.write(settings.jc_genome_diff_file_name);
}
  

  
/*! Passes back calculated values as part of summary
 */
/*
void calculate_cutoffs(const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info)
{
  bool verbose = false;
  
  // distance from read end that must be achieved to pass threshold for different overlap values
  // int32_t max_read_length = summary.sequence_conversion.max_read_length;
  int32_t avg_read_length = round(summary.sequence_conversion.avg_read_length);
  if (verbose) cout << "Average read length: " << avg_read_length << endl;
  
  map<string,int32_t> distance_cutoffs;
  storable_map<string, storable_vector<int32_t> > pos_hash_cutoffs;
  
  for(vector<cAnnotatedSequence>::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); it++) {
    const string& seq_id = it->m_seq_id;
    
    uint32_t sequence_length = ref_seq_info[seq_id].m_length;
    
    double distance_pr_cutoff = sqrt(settings.junction_accept_pr / static_cast<double>(sequence_length));
    double pos_hash_pr_cutoff = settings.junction_accept_pr;
    //double pos_hash_pr_cutoff = sqrt(settings.junction_accept_pr / static_cast<double>(sequence_length));
    
    double pr_no_coverage_position_strand = summary.error_count[seq_id].no_pos_hash_per_position_pr;    
    if (verbose) cout << pr_no_coverage_position_strand << endl;
    if (verbose) cout << "Probability of read starting at position: " << seq_id << " " << (static_cast<double>(1.0)-pr_no_coverage_position_strand) << endl;
      
    if (verbose) cout << "deletion_coverage_propagation_cutoff " << summary.preprocess_coverage[seq_id].deletion_coverage_propagation_cutoff << endl;
    if (verbose) cout << "nbinom_mean_parameter " << summary.preprocess_coverage[seq_id].nbinom_mean_parameter << endl;
    
    double pr_coverage_position_strand = (static_cast<double>(1.0)-pr_no_coverage_position_strand);
    
    
    double junction_coverage_cutoff = summary.preprocess_coverage[seq_id].junction_coverage_cutoff;
    if (verbose) cout << "Junction coverage cutoff: " << seq_id << " " << junction_coverage_cutoff << endl;
    
    double average = summary.preprocess_coverage[seq_id].average;
    if (verbose) cout << "Average coverage: " << seq_id << " " << average << endl;
    
    double adjusted_pr_no_coverage_position_strand = -1;
    if (pr_no_coverage_position_strand > 0) {
      adjusted_pr_no_coverage_position_strand = exp((junction_coverage_cutoff / average) * log(pr_no_coverage_position_strand));
    }
    pr_no_coverage_position_strand = adjusted_pr_no_coverage_position_strand;
    if (verbose) cout << "Adjusted coverage cutoff: " << seq_id << " " << adjusted_pr_no_coverage_position_strand << endl;
    
    for(int32_t i = 0; i <= summary.sequence_conversion.avg_read_length-1; i++) {
      Coverage& cov = summary.preprocess_coverage[seq_id];
      
      // pr_no_coverage_position_strand can be -1 if the fit failed, assign a value of zero for this case
      pos_hash_cutoffs[seq_id].push_back( 
                                         (pr_no_coverage_position_strand > 0) &&  (pr_no_coverage_position_strand < 1)
                                         ? qbinomial(
                                                     pos_hash_pr_cutoff,  
                                                     (2*i), 
                                                     1-pr_no_coverage_position_strand 
                                                     )
                                         : 0
                                         );
      
      if (verbose) cout << "Pos hash cutoffs: " << seq_id << " " << i << " " << pos_hash_cutoffs[seq_id][i] << endl;
      
    }
    
    distance_cutoffs[seq_id] = (pr_no_coverage_position_strand > 0) &&  (pr_no_coverage_position_strand < 1)
    ? qmissing(distance_pr_cutoff, pr_no_coverage_position_strand)
    : 0;
    if (verbose) cout << "Distance cutoff: " << seq_id << " " << distance_cutoffs[seq_id] << endl;
    
  }

  // return values
  summary.alignment_resolution.pos_hash_cutoffs = pos_hash_cutoffs;
  summary.alignment_resolution.distance_cutoffs = distance_cutoffs;
}
*/
    
void load_junction_alignments(
                              const Settings& settings, 
                              Summary& summary, 
                              cReadFiles& read_files, 
                              cReferenceSequences& ref_seq_info,
                              cReferenceSequences& junction_ref_seq_info,
                              SequenceTrimsList& trims_list,
                              map<string,uint32_t>& all_junction_ids,
                              bool junction_prediction,
                              const vector<JunctionInfo>& junction_info_list,
                              UniqueJunctionMatchMap& unique_junction_match_map,
                              RepeatJunctionMatchMap& repeat_junction_match_map,
                              tam_file& resolved_reference_tam
                              )
{
  bool verbose = false;
  uint32_t reads_processed = 0;
  
  tam_file* reference_tam = NULL;
  tam_file* junction_tam = NULL;
  
  for (uint32_t fastq_file_index = 0; fastq_file_index < read_files.size(); fastq_file_index++)
  {    
    const cReadFile& rf = read_files[fastq_file_index];
    string fastq_file_name = read_files.base_name_to_read_file_name(rf.m_base_name);
    
    cerr << "  READ FILE:" << rf.m_base_name << endl;
    
    Summary::AlignmentResolution::ReadFile summary_info;
    
    // Traverse the original fastq files to keep track of order
    // b/c some matches may exist in only one or the other file
    
    cFastqFile in_fastq(fastq_file_name, ios::in);
    
    string this_unmatched_file_name = settings.data_path + "/unmatched."
    + rf.m_base_name + ".fastq";
    cFastqFile out_unmatched_fastq(this_unmatched_file_name, ios::out);
    assert(!out_unmatched_fastq.fail());
    
    string reference_sam_file_name = settings.file_name(settings.reference_sam_file_name, "#", rf.m_base_name);
    string reference_fasta = settings.reference_fasta_file_name;
    
    reference_tam = new tam_file(reference_sam_file_name, settings.reference_fasta_file_name, ios::in); 
    
    if (junction_prediction)
    {
      string junction_sam_file_name = settings.file_name(settings.candidate_junction_sam_file_name, "#", rf.m_base_name);
      junction_tam = new tam_file(junction_sam_file_name, settings.candidate_junction_fasta_file_name, ios::in); 
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
    while (in_fastq.read_sequence(seq)) // READ
    {
      if ((settings.resolve_alignment_read_limit) && (reads_processed >= settings.resolve_alignment_read_limit))
        break; // to next file
      
      reads_processed++;

      
      if (reads_processed % 10000 == 0)
        cerr << "    READS:" << reads_processed << endl;
      
      if (verbose)
        cerr << "===> Read: " << seq.m_name << endl;
      
      uint32_t best_junction_score = 0;
      uint32_t best_reference_score = 0;
      
      // Does this read have eligible candidate junction matches?
      alignment_list this_junction_alignments;
      
      if (verbose)
      {
        cerr << " Before Overlap Reference alignments = " << reference_alignments.size() << endl;
        cerr << " Before Overlap Junction alignments = " << junction_alignments.size() << endl;
      }
      
      if (verbose && (junction_alignments.size() > 0))
      {
        cerr << " Junction SAM read name: " << junction_alignments.front()->read_name() <<endl;
      }
      
      if ((junction_alignments.size() > 0) && (seq.m_name == junction_alignments.front()->read_name()))
      {
        
        this_junction_alignments = junction_alignments;
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
        
        for (alignment_list::iterator it = this_junction_alignments.begin(); it != this_junction_alignments.end(); )
        {
          if (!alignment_overlaps_junction(junction_info_list, it->get()))
            it = this_junction_alignments.erase(it);
          else
            it++; 
        }
        
        best_junction_score = eligible_read_alignments(settings, junction_ref_seq_info, this_junction_alignments);
        
        if (verbose)
          cerr << " Best junction score: " << best_junction_score << endl;
      }
      
      if (verbose && (reference_alignments.size() > 0))
      {
        cerr << " Reference SAM read name: " << reference_alignments.front()->read_name() <<endl;
      }
      
      // Does this read have eligible reference sequence matches?
      alignment_list this_reference_alignments;
      if ((reference_alignments.size() > 0) && (seq.m_name == reference_alignments.front()->read_name()))
      {
        
        this_reference_alignments = reference_alignments;
        reference_tam->read_alignments(reference_alignments, false);
        
        best_reference_score = eligible_read_alignments(settings, ref_seq_info, this_reference_alignments);
        
        if (verbose)
          cerr << " Best reference score: " << best_reference_score << endl;
      }
      
      // Nothing to be done if there were no eligible matches to either
      // Record in the unmatched FASTQ data file
      if ((this_junction_alignments.size() == 0) && (this_reference_alignments.size() == 0))
      {
        summary_info.num_unmatched_reads++;
        out_unmatched_fastq.write_sequence(seq);
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
        cerr << " Mapping quality difference: " << best_reference_score << endl;
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
      if (mapping_quality_difference < 0)
      {
        if (verbose)
          cout << "Best alignment to reference. MQD: " << mapping_quality_difference << endl;
        
        _write_reference_matches(settings, ref_seq_info, trims_list, this_reference_alignments, resolved_reference_tam, fastq_file_index);
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
        // Multiple equivalent matches to junctions and reference, ones with most hits later will win theserepeat matches
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
        
    // save statistics
    summary.alignment_resolution.read_file[read_files[fastq_file_index].m_base_name] = summary_info;
    
    // safe only because we know they are always or never used
    if (junction_tam != NULL) delete junction_tam;
    if (reference_tam != NULL) delete reference_tam;
    
  } // End of Read File loop
  
  // for jump to end of alignments
}
  

/*! Tests whether a read alignment to a candidate junction extends across
 *  the entire junction (past overlapping or unique sequence) and thus
 *  can be used as real evidence for the junction.
 */
bool alignment_overlaps_junction(const vector<JunctionInfo>& junction_info_list, const alignment_wrapper& a)
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


void _write_reference_matches(const Settings& settings, cReferenceSequences& ref_seq_info, const SequenceTrimsList& trims_list, alignment_list& reference_alignments, tam_file& reference_tam, uint32_t fastq_file_index)
{
  (void)settings; //TODO: unused?
	// Nice try, no alignments
	if (reference_alignments.size() == 0) return;

	vector<Trims> trims;

  for(alignment_list::iterator it=reference_alignments.begin(); it!=reference_alignments.end(); it++)
  {
    Trims t = _trim_ambiguous_ends(**it, trims_list);
		trims.push_back(t);
  }
  
	reference_tam.write_alignments((int32_t)fastq_file_index, reference_alignments, &trims, &ref_seq_info, true);
}
  
/*! returns whether it has non overlap alignment
 */
void score_junction(
                    const Settings& settings, 
                    Summary& summary, 
                    const string& junction_id, 
                    UniqueJunctionMatchMap& unique_junction_match_map, 
                    RepeatJunctionMatchMap& repeat_junction_match_map, 
                    tam_file& resolved_junction_tam, 
                    JunctionTestInfo& junction_test_info, 
                    vector<JunctionInfo>& junction_info_list
                  )
{
  // may not need these
  (void) settings;
  (void) summary;
  
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
  
	if (verbose) {
		cout << "Testing Junction Candidate: " << junction_id << endl;
    size_t unique_matches_size = (unique_matches) ? unique_matches->size() : 0;
    size_t repeat_matches_size = (repeat_matches) ? repeat_matches->size() : 0;
		cout << "Unique Matches: " << unique_matches_size << " Degenerate Matches: " << repeat_matches_size << endl;
	}
  
	//// TEST 1: Reads that go a certain number of bp into the nonoverlap sequence on each side of the junction on each strand
	map<bool,int32_t> max_left_per_strand = make_map<bool,int32_t>(true,0)(false,0);
	map<bool,int32_t> max_right_per_strand = make_map<bool,int32_t>(true,0)(false,0);
	map<bool,int32_t> max_min_left_per_strand = make_map<bool,int32_t>(true,0)(false,0);
	map<bool,int32_t> max_min_right_per_strand = make_map<bool,int32_t>(true,0)(false,0);
	map<bool,int32_t> count_per_strand = make_map<bool,int32_t>(true,0)(false,0);
	uint32_t total_non_overlap_reads = 0;
	map<int32_t,bool> pos_hash[2];
  uint32_t pos_hash_count(0);
  
	// basic information about the junction
	JunctionInfo scj(junction_id);
	int32_t alignment_overlap = scj.alignment_overlap;
	int32_t flanking_left = scj.flanking_left;
    
	// We also need to count degenerate matches b/c sometimes ambiguity unfairly penalizes real reads...
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
    
		//! Do not count reads that map the reference equally well toward the score.
		if (item->mapping_quality_difference == 0) {
      if (verbose) cout << "  Degenerate:" << item->junction_alignments.front()->read_name() << endl;
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
    
    // doesn't span 
    if ((this_right <= 0) || (this_left <= 0)) continue;
    
    ////
    // COUNT reads that overlap both sides toward the pos hash score and other statistics
    ////
    
    total_non_overlap_reads++;
    
		bool rev_key = a->reversed();
		count_per_strand[rev_key]++;
    
    // Note that reference here is the junction's sequence, not the reference genome sequence!
    int32_t begin_read_coord = a->reference_start_1();
    
    if (verbose)
			cout << "  " << item->junction_alignments.front()->read_name() << ' ' << static_cast<int32_t>(rev_key) << ' ' << begin_read_coord << endl;
    
    if (!pos_hash[rev_key].count(begin_read_coord))
    {
      pos_hash[rev_key][begin_read_coord] = true;
      pos_hash_count++;
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
  
  // @JEB: TODO broken... needs to construct a unique list of all junction_ids supported by all degenerate matches
  
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
    JunctionInfo main = junction_info_list[*it];
    for (it++; it != repeat_junction_tid_list.end(); it++)
    {
      JunctionInfo test = junction_info_list[*it]; // the junction key
      
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
          if (verbose) cout << "Marking side " << best_side_index << " as redundant." << endl;
        }
      }
    }
  }
  
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
    redundant[0],
    redundant[1],
    junction_id,
    999999999.99
	};
  
  junction_test_info = this_junction_test_info;
}
  
/*! deals with the reads corresponding to a successful or failed junction
 */
void resolve_junction(
                      const Settings& settings,
                      Summary& summary,
                      cReferenceSequences& ref_seq_info,
                      const SequenceTrimsList& trims_list,
                      const string& junction_id,
                      UniqueJunctionMatchMap& unique_junction_match_map,
                      RepeatJunctionMatchMap& repeat_junction_match_map,
                      tam_file& resolved_reference_tam,
                      tam_file& resolved_junction_tam,
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
        if (verbose) cout << "Before size: " << (unique_matches ? unique_matches->size() : 0) << endl;
        
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
        
        if (verbose) cout << "After size: " << (unique_matches ? unique_matches->size() : 0) << endl;
        
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
					_write_reference_matches(settings, ref_seq_info, trims_list, this_reference_al, resolved_reference_tam, fastq_file_index);
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
          alignments.push_back(matched_alignment);
          resolved_junction_tam.write_alignments(fastq_file_index, alignments, NULL, &ref_seq_info, true);
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
        _write_reference_matches(settings, ref_seq_info, trims_list, this_reference_al, resolved_reference_tam, fastq_file_index);
      }
      
      // REGARDLESS of success: write matches to the candidate junction SAM file
      resolved_junction_tam.write_alignments(fastq_file_index, item.junction_alignments, NULL, &ref_seq_info, true);
    }
  }

}

diff_entry junction_to_diff_entry(
                                         const string& key, 
                                         cReferenceSequences& ref_seq_info, 
                                         JunctionTestInfo& test_info
                                         )
{
	// split the key to an item with information about the junction
	JunctionInfo jc(key);

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

	cSequenceFeaturePtr repeat_ptr(NULL);
	for (int32_t i = 0; i <= 1; i++)
	{
		// Determine IS elements
		// Is it within an IS or near the boundary of an IS in the direction leading up to the junction?
		repeat_ptr = cReferenceSequences::find_closest_repeat_region_boundary(jc.sides[i].position, ref_seq_info[jc.sides[i].seq_id].m_repeats, 20, jc.sides[i].strand);
		if (repeat_ptr.get() != NULL)
		{
			jc.sides[i].is = repeat_ptr;
			jc.sides[i].is_interval = (repeat_ptr->m_strand == 1) 
        ? to_string(repeat_ptr->m_start) + "-" + to_string(repeat_ptr->m_end) 
        : to_string(repeat_ptr->m_end) + "-" + to_string(repeat_ptr->m_start);
		}
	}

	// Determine which side of the junction is the IS and which is unique
	// these point to the correct initial interval...
	jc.is_side = UNDEFINED_UINT32;
	if (jc.sides[0].is.get() && !jc.sides[1].is.get())
	{
		if (abs(static_cast<int32_t>(jc.sides[0].is->m_start) - static_cast<int32_t>(jc.sides[0].position)) <= 20)
		{
			jc.is_side = 0;
			jc.sides[jc.is_side].is_side_key = "start";
		}
		else if (abs(static_cast<int32_t>(jc.sides[0].is->m_end) - static_cast<int32_t>(jc.sides[0].position)) <= 20 )
		{
			jc.is_side = 0;
			jc.sides[jc.is_side].is_side_key = "end";
		}
		jc.unique_side = 1;
	}

	else if (!jc.sides[0].is.get() && jc.sides[1].is.get())
	{
		if (abs(static_cast<int32_t>(jc.sides[1].is->m_start) - static_cast<int32_t>(jc.sides[1].position)) <= 20)
		{
			jc.is_side = 1;
			jc.sides[jc.is_side].is_side_key = "start";
		}
		else if (abs(static_cast<int32_t>(jc.sides[1].is->m_end) - static_cast<int32_t>(jc.sides[1].position)) <= 20 )
		{
			jc.is_side = 1;
			jc.sides[jc.is_side].is_side_key = "end";
		}
		jc.unique_side = 0;
	}

// @JEB TODO: Testing removal
// both were IS! -- define as redundant here
//	else if (jc.sides[0].is.name.size() > 0)
//		jc.sides[0].redundant = true;
//	else if (jc.sides[1].is.name.size() > 0)
//		jc.sides[1].redundant = true;

//  @JEB TODO: Testing removal
// Carry over redundancy from degenerate matches
//  if (test_info.redundant_1) jc.sides[0].redundant = true;
//  if (test_info.redundant_2) jc.sides[1].redundant = true;
  
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
          ? jc.sides[jc.is_side].is->m_start : jc.sides[jc.is_side].is->m_end)) - jc.sides[jc.is_side].position);

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
			uint32_t strand_direction = (jc.sides[1].strand > 0 ? 1 : -1);
			jc.sides[1].position += jc.overlap * strand_direction;
			jc.sides[1].overlap = 0;
			jc.overlap = 0;
		}
		else  // side_1 was redundant, give overlap to side_2
		{
			uint32_t strand_direction = (jc.sides[0].strand > 0 ? -1 : 1);
			jc.sides[0].position += jc.overlap * strand_direction;
			jc.sides[0].overlap = 0;
			jc.overlap = 0;
		}

		// If both sides were redundant, no adjustment because we are not going to count coverage
	}

	// flatten things to only what we want to keep
  diff_entry item(JC);
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
  
//	## may want to take only selected of these fields..
  
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
  ("neg_log10_pos_hash_p_value", to_string(test_info.neg_log10_pos_hash_p_value))
  ;

	/// Note: Other adjustments to overlap can happen at the later annotation stage
	/// and they will not affect coverage for calling deletions or mutations
	/// because they will be in REDUNDANTLY matched sides of junctions
	return item;
}

  
Trims _trim_ambiguous_ends(const alignment_wrapper& a, const SequenceTrimsList& trims)
{
	bool verbose = false;
  
	// which reference sequence?
	uint32_t tid = a.reference_target_id();
  
  Trims t;
  t.L= trims[tid].left_trim_0(a.reference_start_0());
  t.R = trims[tid].right_trim_0(a.reference_end_0());

  t.L += a.query_start_0();
  t.R += a.read_length() - a.query_end_1();

//  cerr << a.read_name() << endl;
//  cerr << "start: " << a.reference_start_1() << " end: " << a.reference_end_1() << endl;
//  cerr << "left: " << t.L << " right: " << t.R << endl;

  return t;
}
  
void read_trims(SequenceTrimsList& trims, const cReferenceSequences& ref_seqs, const string &in_trims_file_name ) 
{
  trims.resize(ref_seqs.size());
  for(uint32_t i = 0; i < ref_seqs.size(); i++) {
    string this_file_name = Settings::file_name(in_trims_file_name, "@", ref_seqs[i].m_seq_id);
    trims[i].ReadFile(this_file_name, ref_seqs[i].m_length);
  }
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
  
} // namespace breseq

