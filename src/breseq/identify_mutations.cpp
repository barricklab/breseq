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

#include "common.h"
#include "pileup.h"
#include "identify_mutations.h"
#include "error_count.h"
#include "reference_sequence.h"

using namespace std;

namespace breseq {

/*! Convenience wrapper around the identify_mutations_pileup class.
 */
  
void identify_mutations(
                const Settings& settings,
                const Summary& summary,
								const string& bam,
								const string& fasta,
								const string& gd_file,
                const cReferenceSequences& ref_seq_info,
                const vector<double>& deletion_propagation_cutoff,
                const vector<double>& deletion_seed_cutoffs,
								double mutation_cutoff,
								double polymorphism_cutoff,
								double polymorphism_frequency_cutoff,
                double polymorphism_precision_decimal,
                uint32_t polymorphism_precision_places,
								bool print_per_position_file
 ) {

	// do the mutation identification:
	identify_mutations_pileup imp(
                settings,
                summary,
								bam,
								fasta,
                deletion_propagation_cutoff,
                deletion_seed_cutoffs,
								mutation_cutoff,
								polymorphism_cutoff,
								polymorphism_frequency_cutoff,
                polymorphism_precision_decimal,
                polymorphism_precision_places,
								print_per_position_file
							);
	imp.do_pileup(settings.call_mutations_seq_id_set());
  if (settings.predict_soft_clipping) imp.add_sc_evidence(summary, ref_seq_info);
  // Write discordant-pair candidate regions accumulated during the pileup (single operation).
  imp.write_dp_candidate_regions(settings.dp_candidate_regions_file_name);
  imp.write_gd(gd_file);
}

/*! Chooses which RA predictions to reject and for what reasons
 */
enum ePredictionType {unknown, consensus, polymorphism};


bool rejected_RA_polymorphism_bias(cDiffEntry& ra,
                                       cReferenceSequences& ref_seq_info,
                                       const Settings& settings
                                       )
{
  (void)ref_seq_info;
  bool rejected = false;


  ////// Strand and base quality score biases
  
  // Test quality score bias
  if (settings.polymorphism_ks_quality_p_value_cutoff) {
    if (ra.entry_exists("ks_quality_p_value")) {
      if (from_string<double>(ra["ks_quality_p_value"]) < settings.polymorphism_ks_quality_p_value_cutoff) {
        rejected = true;
        ra.add_reject_reason("KS_BASE_QUALITY");
      }
    }
  }
  
  if (settings.polymorphism_fisher_strand_p_value_cutoff) {
    if (ra.entry_exists("fisher_strand_p_value")) {
      if (from_string<double>(ra["fisher_strand_p_value"]) < settings.polymorphism_fisher_strand_p_value_cutoff) {
        rejected = true;
        ra.add_reject_reason("FISHER_STRAND");
      }
    }
  }
  
  return rejected;
}
  
// ---------------------------------------------------------------------------------------------
// Exact (Clopper-Pearson) confidence bounds on the variant frequency at an RA position, so the
// frequency cutoffs are a confidence statement rather than a point-estimate comparison: a call is
// never rejected merely for having shallow coverage, only for being confidently on the wrong side.
//
// The interval is CENTRED on the reported POLYMORPHISM_FREQUENCY -- a maximum-likelihood fit that
// weights by base quality -- rather than on the raw new_cov/total_cov ratio, and given the width
// implied by the total read depth. The two centres agree to a median of 4e-5 across the test
// goldens but can differ by up to 0.11 where quality is heterogeneous, and the bound has to
// describe the same number the report shows; otherwise "rejected at frequency 0.43" carrying a
// bound derived from 0.32 is unreadable. Depth still sets the width, so the interval widens
// correctly where coverage is thin, which is the entire point.
//
// Returns false if the position has no usable depth, in which case callers fall back to the
// point estimate.
static bool RA_frequency_bounds(const cDiffEntry& ra, double& lower, double& upper)
{
  if (!ra.entry_exists(TOTAL_COV) || !ra.entry_exists(POLYMORPHISM_FREQUENCY)) return false;
  vector<string> top_bot = split(ra.get(TOTAL_COV), "/");
  if (top_bot.size() < 2) return false;
  double n = from_string<double>(top_bot[0]) + from_string<double>(top_bot[1]);
  if (!(n > 0.0)) return false;
  double k = from_string<double>(ra.get(POLYMORPHISM_FREQUENCY)) * n;
  lower = binomial_frequency_lower_bound(k, n);
  upper = binomial_frequency_upper_bound(k, n);
  return true;
}



// The bounds are deliberately NOT written onto the entry. They are recoverable at any time from
// total_cov + polymorphism_frequency (see ra_frequency_confidence_bounds in output.cpp, which
// recomputes them for the report), and storing them would add two fields to all ~840 RA entries in
// the test goldens, burying the handful of genuine classification changes in cosmetic churn.

bool rejected_RA_polymorphism_coverage(cDiffEntry& ra,
                               cReferenceSequences& ref_seq_info,
                               const Settings& settings
                               )
{
  (void)ref_seq_info;
  bool rejected = false;

  // Minimum coverage on both strands for the NON-REFERENCE allele only.
  //
  // This guard exists to catch a VARIANT that appears on only one strand -- the signature of a
  // chemistry artifact. A REFERENCE allele that is thin on one strand says nothing about that; it
  // says the position is closer to fixed, which is a statement about frequency and is already
  // answered by the confidence bounds. Checking it here rejected well-supported calls for the
  // wrong reason: a position 82% variant, with the variant itself at 14/18, failed because the 7
  // residual reference reads happened to fall 1/6.
  if (settings.polymorphism_minimum_variant_coverage_each_strand > 0) {
    bool passed = true;
    string ref_base = ra.entry_exists(REF_BASE) ? ra[REF_BASE] : "";
    bool major_is_variant = ra.entry_exists(MAJOR_BASE) && (ra[MAJOR_BASE] != ref_base);
    bool minor_is_variant = ra.entry_exists(MINOR_BASE) && (ra[MINOR_BASE] != ref_base);
    // Degenerate case -- neither allele differs from the reference. Check both, as before.
    if (!major_is_variant && !minor_is_variant) { major_is_variant = true; minor_is_variant = true; }

    if (major_is_variant && ra.entry_exists(MAJOR_COV)) {
      vector<string> top_bot = split(ra[MAJOR_COV], "/");
      passed = passed && (from_string<double>(top_bot[0]) >= settings.polymorphism_minimum_variant_coverage_each_strand);
      passed = passed && (from_string<double>(top_bot[1]) >= settings.polymorphism_minimum_variant_coverage_each_strand);
    }
    if (minor_is_variant && ra.entry_exists(MINOR_COV)) {
      vector<string> top_bot = split(ra[MINOR_COV], "/");
      passed = passed && (from_string<double>(top_bot[0]) >= settings.polymorphism_minimum_variant_coverage_each_strand);
      passed = passed && (from_string<double>(top_bot[1]) >= settings.polymorphism_minimum_variant_coverage_each_strand);
    }

    if (!passed) {
      rejected = true;
      ra.add_reject_reason("VARIANT_STRAND_COVERAGE");
    }
  }
  
  if (settings.polymorphism_minimum_total_coverage_each_strand > 0) {
    bool passed = true;
    vector<string> top_bot;
    double top;
    double bot;
    
    top_bot = split(ra[TOTAL_COV], "/");
    top = from_string<double>(top_bot[0]);
    bot = from_string<double>(top_bot[1]);
    passed = passed && (top >= settings.polymorphism_minimum_total_coverage_each_strand);
    passed = passed && (bot >= settings.polymorphism_minimum_total_coverage_each_strand);
    
    if (!passed) {
      rejected = true;
      ra.add_reject_reason("TOTAL_STRAND_COVERAGE");
    }
  }
  
  if (settings.polymorphism_minimum_variant_coverage > 0) {
    
    bool passed = true;
    vector<string> top_bot;
    double top;
    double bot;
    
    top_bot = split(ra[MAJOR_COV], "/");
    top = from_string<double>(top_bot[0]);
    bot = from_string<double>(top_bot[1]);
    passed = passed && (top+bot >= settings.polymorphism_minimum_variant_coverage);
    
    top_bot = split(ra[MINOR_COV], "/");
    top = from_string<double>(top_bot[0]);
    bot = from_string<double>(top_bot[1]);
    passed = passed && (top+bot >= settings.polymorphism_minimum_variant_coverage);
    
    if (!passed) {
      rejected = true;
      ra.add_reject_reason("VARIANT_COVERAGE");
    }
  }
  
  if (settings.polymorphism_minimum_total_coverage > 0) {
    vector<string> top_bot = split(ra[TOTAL_COV], "/");
    if ( from_string<double>(top_bot[0]) + from_string<double>(top_bot[1]) < settings.polymorphism_minimum_total_coverage ) {
      ra.add_reject_reason("TOTAL_COVERAGE");
      rejected = true;
    }
  }
  
  return rejected;
}
  
bool rejected_RA_consensus_coverage(cDiffEntry& ra,
                                    cReferenceSequences& ref_seq_info,
                                    const Settings& settings
                                    )
{
  (void)ref_seq_info;
  bool rejected = false;
  
  if (settings.consensus_minimum_variant_coverage_each_strand > 0) {
    vector<string> top_bot = split(ra[MAJOR_COV], "/");
    if ( (from_string<double>(top_bot[0]) < settings.consensus_minimum_variant_coverage_each_strand) ||
        (from_string<double>(top_bot[1]) < settings.consensus_minimum_variant_coverage_each_strand) ) {
      ra.add_reject_reason("VARIANT_STRAND_COVERAGE");
      rejected = true;
    }
  }
  
  if (settings.consensus_minimum_total_coverage_each_strand > 0) {
    vector<string> top_bot = split(ra[TOTAL_COV], "/");
    if ( (from_string<double>(top_bot[0]) < settings.consensus_minimum_variant_coverage_each_strand) ||
        (from_string<double>(top_bot[1]) < settings.consensus_minimum_variant_coverage_each_strand) ) {
      ra.add_reject_reason("TOTAL_STRAND_COVERAGE");
      rejected = true;
    }
  }
  
  if (settings.consensus_minimum_variant_coverage > 0) {
    vector<string> top_bot = split(ra[MAJOR_COV], "/");
    if ( from_string<double>(top_bot[0]) + from_string<double>(top_bot[1]) < settings.consensus_minimum_total_coverage ) {
      ra.add_reject_reason("VARIANT_COVERAGE");
      rejected = true;
    }
  }
  
  if (settings.consensus_minimum_total_coverage > 0) {
    vector<string> top_bot = split(ra[TOTAL_COV], "/");
    if ( from_string<double>(top_bot[0]) + from_string<double>(top_bot[1]) < settings.consensus_minimum_total_coverage ) {
      ra.add_reject_reason("TOTAL_COVERAGE");
      rejected = true;
    }
  }
  
  return rejected;
}
  
bool rejected_RA_indel_homopolymer(cDiffEntry& ra,
                              cReferenceSequences& ref_seq_info,
                              uint32_t reject_indel_homopolymer_length,
                              uint32_t reject_surrounding_homopolymer_length,
                              bool no_indel_polymorphisms
                               )
{
  bool rejected = false;
  
  ////// Optionally, ignore if in a homopolymer stretch longer than this
  // There are several cases here that may be of interest for removing false-positives
  // 1) indel in what was originally a homopolymer adding a base that is the same
  // 2) A mutation in the middle of a stretch of one base to what is in the rest of that stretch
  //    TTTTTTATTTTTT
  //    TTTTTTTTTTTTT
  //
  // First case,
  // 1) indel in what was originally a homopolymer adding a base that is the same
  if (reject_indel_homopolymer_length && ((ra[REF_BASE] == ".") || (ra[NEW_BASE] == ".")))
  {
    // Code needs to be robust to the mutation being at the beginning
    // OR end of the homopolymer tract
    
    string seq_id = ra["seq_id"];
    int32_t mut_pos = from_string<uint32_t>(ra["position"]);
    bool is_insertion = (from_string<int32_t>(ra[INSERT_POSITION]) > 0);
    string mut_base = (ra[REF_BASE] == ".") ? ra[NEW_BASE] : ra[REF_BASE];
    int32_t homopolymer_length = 0;
    
    // If we are an insertion, we need to start on the base before or the base after,
    // depending on which one we match (if we match neither than our length is zero
    bool no_match = false;
    if (is_insertion) {
      
      // This is the base before
      if (ref_seq_info.get_sequence_1(seq_id, mut_pos, mut_pos) != mut_base) {
        
        // No match? -- Now check the base after
        mut_pos++;
        if (mut_pos > static_cast<int32_t>(ref_seq_info.get_sequence_length(seq_id))) {
          no_match = true;
        } else if (ref_seq_info.get_sequence_1(seq_id, mut_pos, mut_pos) != mut_base) {
          no_match = true;
        }
      }
    }
    
    // extend match
    if (!no_match) {
      
      // check before
      int32_t start_pos = mut_pos - 1;
      while ( (ref_seq_info.get_sequence_1(seq_id, start_pos, start_pos) == mut_base) && (start_pos > 0) ) {
        start_pos--;
      }
      start_pos++;
      
      int32_t end_pos = mut_pos + 1;
      while ( (ref_seq_info.get_sequence_1(seq_id, end_pos, end_pos) == mut_base) && (end_pos <= static_cast<int32_t>(ref_seq_info.get_sequence_length(seq_id)) ) ) {
        end_pos++;
      }
      end_pos--;
      homopolymer_length = end_pos - start_pos + 1;
    }
    
    if (homopolymer_length >= static_cast<int32_t>(reject_indel_homopolymer_length))
    {
      rejected = true;
      ra.add_reject_reason("INDEL_HOMOPOLYMER");
    }
  }
  
  // Second case
  // 2) A mutation in the middle of a stretch of one base to what is in the rest of that stretch
  //    TTTTTTATTTTTT
  //    TTTTTTTTTTTTT
  if (reject_surrounding_homopolymer_length && (ra[REF_BASE] != ".") && (ra[NEW_BASE] != "."))    {
    
    string seq_id = ra["seq_id"];
    int32_t mut_pos = from_string<int32_t>(ra["position"]);
    string ref_base = ra[REF_BASE];
    string mut_base = ra[NEW_BASE];
    
    // determine the start and end of homopolymer created by substitution
    int32_t start_pos = mut_pos - 1;
    while (start_pos >= 1)
    {
      if (ref_seq_info.get_sequence_1(seq_id, start_pos, start_pos) != mut_base) {
        break;
      }
      start_pos--;
    }
    start_pos++;
    
    int32_t end_pos = mut_pos + 1;
    while (end_pos <= static_cast<int32_t>(ref_seq_info.get_sequence_length(seq_id)))
    {
      if (ref_seq_info.get_sequence_1(seq_id, end_pos, end_pos) != mut_base) {
        break;
      }
      end_pos++;
    }
    end_pos--;
    
    // check bounds
    if ( (start_pos < mut_pos) && (end_pos > mut_pos) && (end_pos - start_pos + 1 >= static_cast<int32_t>(reject_surrounding_homopolymer_length)) ) {
      rejected = true;
      ra.add_reject_reason("SURROUNDING_HOMOPOLYMER");
    }
  }
  ////// END checking on homopolymers
  
  if (no_indel_polymorphisms && ((ra[REF_BASE] == ".") || (ra[NEW_BASE] == ".")))
  {
    rejected = true;
    ra.add_reject_reason("POLYMORPHIC_INDEL");
  }
  return rejected;
}
  
// Returns whether it should be deleted
//
// CONSENSUS MODE asks one question of every position: is there evidence that the variant is the
// MAJORITY allele. That gives two flat attempts, in order, each gated on a score and on the exact
// 95% lower confidence bound of the variant frequency:
//
//   1. consensus     -- consensus_score clears its cutoff AND L >= consensus_frequency_cutoff
//                       (0.50 by default: confidently more than half the reads)
//   2. polymorphism  -- otherwise, polymorphism_score clears its cutoff AND
//                       L >= polymorphism_frequency_cutoff (0.10 by default: confidently present)
//   3. otherwise the position is dropped.
//
// Both tests read the LOWER bound, so a position is never accepted on thin coverage and never
// rejected merely for having thin coverage -- only for failing to establish the claim being made.
// The strand-minimum, Fisher strand-bias, K-S quality-bias and homopolymer filters are all OFF by
// default here and run only when explicitly enabled; they are applied within whichever attempt
// they belong to so those command-line options keep working.
bool test_RA_evidence_CONSENSUS_mode(
                                     cDiffEntry& ra,
                                     cReferenceSequences& ref_seq_info,
                                     const Settings& settings
                                     )
{
  double consensus_score = double_from_string(ra[CONSENSUS_SCORE]);
  double polymorphism_score = double_from_string(ra[POLYMORPHISM_SCORE]);

  // Exact 95% lower bound on the variant frequency; falls back to the point estimate if the
  // position has no usable depth.
  double lower, upper;
  if (!RA_frequency_bounds(ra, lower, upper))
    lower = from_string<double>(ra[POLYMORPHISM_FREQUENCY]);

  /////////////////////////////////
  // 1. Is it the majority allele?
  /////////////////////////////////

  if (consensus_score < settings.mutation_log10_e_value_cutoff) {
    ra.add_reject_reason("SCORE_CUTOFF");
  }
  if ((settings.consensus_frequency_cutoff > 0.0) && (lower < settings.consensus_frequency_cutoff)) {
    ra.add_reject_reason("FREQUENCY_CUTOFF");
  }
  rejected_RA_consensus_coverage(ra, ref_seq_info, settings);
  rejected_RA_indel_homopolymer(ra,
                                ref_seq_info,
                                settings.consensus_reject_indel_homopolymer_length,
                                settings.consensus_reject_surrounding_homopolymer_length,
                                false
                                );

  if (!ra.entry_exists(REJECT)) {
    ra[PREDICTION] = "consensus";
    ra[FREQUENCY] = "1";
    // Delete if we are just the reference base!
    return (ra[REF_BASE] == ra[MAJOR_BASE]);
  }
  ra[CONSENSUS_REJECT] = ra[REJECT];
  ra.clear_reject_reasons();

  /////////////////////////////////
  // 2. Is it present at all?
  /////////////////////////////////

  if (polymorphism_score < settings.polymorphism_log10_e_value_cutoff) {
    ra.add_reject_reason("SCORE_CUTOFF");
  }
  if ((settings.polymorphism_frequency_cutoff > 0.0) && (lower < settings.polymorphism_frequency_cutoff)) {
    ra.add_reject_reason("FREQUENCY_CUTOFF");
  }
  rejected_RA_polymorphism_bias(ra, ref_seq_info, settings);
  rejected_RA_polymorphism_coverage(ra, ref_seq_info, settings);
  rejected_RA_indel_homopolymer(ra,
                                ref_seq_info,
                                settings.polymorphism_reject_indel_homopolymer_length,
                                settings.polymorphism_reject_surrounding_homopolymer_length,
                                settings.polymorphism_no_indels
                                );

  if (!ra.entry_exists(REJECT)) {
    ra[PREDICTION] = "polymorphism";
    ra[FREQUENCY] = ra[POLYMORPHISM_FREQUENCY];
    return false;
  }
  ra[POLYMORPHISM_REJECT] = ra[REJECT];
  ra.clear_reject_reasons();

  /////////////////////////////////
  // 3. Neither claim holds.
  /////////////////////////////////

  return true;
}

// Returns whether it should be deleted
//
// POLYMORPHISM MODE runs the same two attempts as consensus mode, but the null hypothesis is
// reversed: here a position is assumed to be POLYMORPHIC unless the data say otherwise, so the
// bounds swap roles.
//
//   1. consensus     -- consensus_score clears its cutoff AND U >= consensus_frequency_cutoff
//                       (0.99 by default). The upper bound, not the lower: a call snaps to a fixed
//                       100% difference only when we cannot rule out that it IS fixed. Under
//                       consensus mode the same slot asks "is it confidently the majority" with the
//                       lower bound, because there the fixed interpretation is the default.
//   2. polymorphism  -- otherwise, polymorphism_score clears its cutoff AND
//                       L >= polymorphism_frequency_cutoff (0.01 by default: confidently present)
//   3. otherwise the position is KEPT and marked rejected -- unlike consensus mode, which drops it.
//
// The optional guards sit in whichever attempt they belong to and are all OFF by default.
bool test_RA_evidence_POLYMORPHISM_mode(
                                     cDiffEntry& ra,
                                     cReferenceSequences& ref_seq_info,
                                     const Settings& settings
                                     )
{
  double consensus_score = double_from_string(ra[CONSENSUS_SCORE]);
  double polymorphism_score = double_from_string(ra[POLYMORPHISM_SCORE]);

  double lower, upper;
  if (!RA_frequency_bounds(ra, lower, upper))
    lower = upper = from_string<double>(ra[POLYMORPHISM_FREQUENCY]);

  /////////////////////////////////
  // 1. Can we still call it fixed?
  /////////////////////////////////

  if (consensus_score < settings.mutation_log10_e_value_cutoff) {
    ra.add_reject_reason("SCORE_CUTOFF");
  }
  if ((settings.consensus_frequency_cutoff > 0.0) && (upper < settings.consensus_frequency_cutoff)) {
    ra.add_reject_reason("FREQUENCY_CUTOFF");
  }
  rejected_RA_consensus_coverage(ra, ref_seq_info, settings);
  rejected_RA_indel_homopolymer(ra,
                                ref_seq_info,
                                settings.consensus_reject_indel_homopolymer_length,
                                settings.consensus_reject_surrounding_homopolymer_length,
                                false
                                );

  if (!ra.entry_exists(REJECT)) {
    ra[PREDICTION] = "consensus";
    ra[FREQUENCY] = "1";
    // Delete if we are just the reference base!
    return (ra[REF_BASE] == ra[MAJOR_BASE]);
  }
  ra[CONSENSUS_REJECT] = ra[REJECT];
  ra.clear_reject_reasons();

  /////////////////////////////////
  // 2. Is it present at all?
  /////////////////////////////////

  if (polymorphism_score < settings.polymorphism_log10_e_value_cutoff) {
    ra.add_reject_reason("SCORE_CUTOFF");
  }
  if ((settings.polymorphism_frequency_cutoff > 0.0) && (lower < settings.polymorphism_frequency_cutoff)) {
    ra.add_reject_reason("FREQUENCY_CUTOFF");
  }
  rejected_RA_polymorphism_bias(ra, ref_seq_info, settings);
  rejected_RA_polymorphism_coverage(ra, ref_seq_info, settings);
  rejected_RA_indel_homopolymer(ra,
                                ref_seq_info,
                                settings.polymorphism_reject_indel_homopolymer_length,
                                settings.polymorphism_reject_surrounding_homopolymer_length,
                                settings.polymorphism_no_indels
                                );

  if (!ra.entry_exists(REJECT)) {
    ra[PREDICTION] = "polymorphism";
    ra[FREQUENCY] = ra[POLYMORPHISM_FREQUENCY];
    return false;
  }

  /////////////////////////////////
  // 3. Neither claim holds -- keep it, marked rejected.
  /////////////////////////////////

  ra[POLYMORPHISM_REJECT] = ra[REJECT];
  ra.clear_reject_reasons();
  ra[PREDICTION] = "polymorphism";
  ra[FREQUENCY] = ra[POLYMORPHISM_FREQUENCY];
  ra[REJECT] = ra[POLYMORPHISM_REJECT];
  ra.erase(POLYMORPHISM_REJECT);

  return false;
}

void test_RA_evidence(
                        cGenomeDiff& gd,
                        cReferenceSequences& ref_seq_info,
                        const Settings& settings
                        )
{
  
  // Assumes these entries are present initially
  //  * REF_BASE
  //  * CONSENSUS_SCORE is present for all entries        = maximumum likelihood frequency for 100% mutated model
  //  * POLYMORPHISM_FREQUENCY is present for all entries = maximumum likelihood frequency for mixed model
  //  * QUALITY is present for all entries                = consensus model
  //  * POLYMORPHISM_SCORE may or may not be present      = mixed model
  //  * NEW_COV, REF_COV present for all entries
  
  diff_entry_list_t list = gd.get_list();
  // Note nonstandard non-increment, since we remove items
  for (diff_entry_list_t::iterator it = list.begin(); it != list.end(); ) {
    
    cDiffEntry& ra = **it;
    if (ra._type != RA) {
      it++;
      continue;
    }
    
    // At this point, both consensus and polymorphism quality may be > 0,
    // we see if either is below the significance + frequency cutoff first.

    if (!ra.entry_exists(REF_BASE)) {
      ERROR("Expected field 'ref_base' in evidence item\n" + ra.as_string());
    }
    if (!ra.entry_exists(CONSENSUS_SCORE)) {
      ERROR("Expected field 'consensus_score' in evidence item\n" + ra.as_string());
    }
    if (!ra.entry_exists(POLYMORPHISM_SCORE)) {
      ERROR("Expected field 'polymorphism_score' in evidence item\n" + ra.as_string());
    }
    if (!ra.entry_exists(POLYMORPHISM_FREQUENCY)) {
      ERROR("Expected field '" + cString(POLYMORPHISM_FREQUENCY) + "' in evidence item\n" + ra.as_string());
    }
    
    bool delete_entry = false;

    if (settings.polymorphism_prediction) {
      delete_entry = test_RA_evidence_POLYMORPHISM_mode(ra, ref_seq_info, settings);
    } else {
      delete_entry = test_RA_evidence_CONSENSUS_mode(ra, ref_seq_info, settings);
    }
  
    // User evidence items should be treated normally (to define polymorphism vs. consensus predictions)
    // but they should not be rejected, and they should not be deleted.
    if (ra.entry_exists(USER_DEFINED)) {
      delete_entry=false;
    }

    if (delete_entry) {
      it = list.erase(it);
    } else {
      it++;
    }
  
  } // end loop over ra evidence items

  // Resave the list (from which we may have removed items)
  gd.set_list(list);
}


/*! Constructor.
 */
identify_mutations_pileup::identify_mutations_pileup(
                              const Settings& settings,
                              const Summary& summary,
															const string& bam,
															const string& fasta,
                              const vector<double>& deletion_propagation_cutoffs,
                              const vector<double>& deletion_seed_cutoffs,
															double consensus_score_cutoff,
															double polymorphism_score_cutoff,
															double polymorphism_frequency_cutoff,
                              double polymorphism_precision_decimal,
                              uint32_t polymorphism_precision_places,
															bool print_per_position_file
                                                            )
: pileup_base(bam, fasta)
, _settings(settings)
, _gd()
, _deletion_seed_cutoffs(deletion_seed_cutoffs)
, _deletion_propagation_cutoffs(deletion_propagation_cutoffs)
, _consensus_score_cutoff(consensus_score_cutoff)
, _polymorphism_score_cutoff(polymorphism_score_cutoff)
, _polymorphism_frequency_cutoff(polymorphism_frequency_cutoff)
, _polymorphism_precision_decimal(polymorphism_precision_decimal)
, _polymorphism_precision_places(polymorphism_precision_places)
, _log10_ref_length(0)
, _total_reference_length(summary.sequence_conversion.total_reference_sequence_length)
, _snp_caller("haploid", summary.sequence_conversion.total_reference_sequence_length)
, _this_deletion_reaches_seed_value(false)
, _this_deletion_redundant_reached_zero(false)
, _last_position_coverage_printed(0)
, _print_per_position_file(print_per_position_file)
, _dp_seed(settings.discordant_pair_seed)
{

  // remove once used
  (void)settings;

  // Initialize per-bin discordant-pair region state (bin = strand*3 + orient; see kDPnBins).
  for (int bin = 0; bin < kDPnBins; bin++) {
    _dp_metric[bin] = 0;
    _dp_region_start[bin] = UNDEFINED_UINT32;
    _dp_region_max_count[bin] = 0;
    _dp_region_redundant_count[bin] = 0;
  }

  set_print_progress(true);
  
  ASSERT(m_bam_header->n_targets == (int32_t)_deletion_propagation_cutoffs.size(), 
         "Number of targets in BAM file [" + to_string(m_bam_header->n_targets) + "] " +
         "does not match + number in cutoff table [" + to_string(_deletion_propagation_cutoffs.size()) + "]."
         );
  ASSERT(m_bam_header->n_targets == (int32_t)_deletion_seed_cutoffs.size(),
         "Number of targets in BAM file [" + to_string(m_bam_header->n_targets) + "] " +
         "does not match + number in cutoff table [" + to_string(_deletion_propagation_cutoffs.size()) + "]."
         );
    
	// reserve enough space for the sequence info:
	_seq_info.resize(m_bam_header->n_targets);
	
  // tally up the reference lengths from the bam file
	for(int i=0; i<m_bam_header->n_targets; ++i) {
		_log10_ref_length += static_cast<double>(m_bam_header->target_len[i]);
	}
	assert(_log10_ref_length != 0);
	_log10_ref_length = log10(_log10_ref_length);
  
  // are we printing detailed coverage information? (only needed when --predict_copy_number reads it back)
  _print_coverage_data = _settings.predict_copy_number;
  
  // load the error table file and convert back to probabilities
  _error_table.read_log10_prob_table(settings.error_rates_file_name);
  _error_table.log10_prob_to_prob();
  
  if (_print_per_position_file) {
    _per_position_file.open(_settings.mutation_identification_per_position_file_name.c_str());
  }
  
  // load the list of user RA evidence
  if (_settings.user_evidence_genome_diff_file_name != "") {
    load_user_ra_evidence_from_gd();
  }

  // Set up discordant-pair (DP) candidate-region detection.
  // For each PAIRED read group with a valid paired-mapping-distance fit, build a sliding-window
  // group with width W = median + 2.42 * MAD, and map each member file's BAM index (cReadFile::m_id)
  // to that group so we can classify a read's group by its fastq_file_index() during the pileup.
  const PairedMappingDistanceDistributionSummaries& pmdd = summary.preliminary_paired_mapping_distance_distribution;
  for (vector<cReadFileSet>::const_iterator rfs = settings.read_file_sets.begin(); rfs != settings.read_file_sets.end(); rfs++) {
    if (!rfs->is_paired()) continue;
    PairedMappingDistanceDistributionSummaries::const_iterator it = pmdd.find(rfs->m_base_name);
    if (it == pmdd.end()) continue;
    if (it->second.median <= 0.0) continue;

    dp_group g;
    g.window_width = it->second.median + 2.42 * it->second.mad;
    g.r1_m_id = rfs->m_files[0].m_id;  // R1 fastq file index
    g.r2_m_id = rfs->m_files[1].m_id;  // R2 fastq file index
    int group_index = static_cast<int>(_dp_groups.size());
    _dp_groups.push_back(g);
    for (vector<cReadFile>::const_iterator rf = rfs->m_files.begin(); rf != rfs->m_files.end(); rf++) {
      _fastq_index_to_dp_group[rf->m_id] = group_index;
    }
  }

}

void identify_mutations_pileup::load_user_ra_evidence_from_gd()
{
  cGenomeDiff gd(_settings.user_evidence_genome_diff_file_name);
  gd.sort(); // very important to be sorted!!
  
  // get rid of any previous annotations
  gd.strip_to_spec();
  
  // Return just RA entries
  _user_evidence_ra_list = gd.get_list(make_vector<gd_entry_type>(RA));
}

void identify_mutations_pileup::add_sc_evidence(const Summary& summary, const cReferenceSequences& ref_seq_info)
{
  if (!file_exists(_settings.soft_clipping_counts_file_name.c_str())) return;

  // The null model is estimated at tabulation time (soft_clipping.cpp), where the
  // zero-clip positions it needs are still enumerable. The counts file only lists
  // positions with at least one clip event, so it cannot be re-derived here.
  double p0 = summary.soft_clipping.soft_clipping_null_rate;
  double rho = summary.soft_clipping.soft_clipping_dispersion;

  // Beta-binomial shape parameters. rho is the intra-class correlation, so the
  // variance is inflated by (1 + (n-1)*rho) relative to a binomial; rho == 0 falls
  // back to the plain binomial via bdtrc below.
  double alpha = 0.0, beta = 0.0;
  bool use_beta_binomial = (rho > 0.0) && (rho < 1.0);
  if (use_beta_binomial) {
    alpha = p0 * (1.0 - rho) / rho;
    beta = (1.0 - p0) * (1.0 - rho) / rho;
    if (!(alpha > 0.0) || !(beta > 0.0)) use_beta_binomial = false;
  }

  // Bonferroni: 2 directions × total genome length positions
  uint64_t total_genome_length = 0;
  for (uint32_t i = 0; i < ref_seq_info.size(); i++) {
    total_genome_length += ref_seq_info[i].m_length;
  }
  double log10_n_tests = log10(2.0 * static_cast<double>(total_genome_length));

  ifstream in(_settings.soft_clipping_counts_file_name.c_str());
  if (!in.is_open()) return;

  // Candidates are collected first and only turned into evidence at the end, because the
  // local non-maximum suppression below needs to compare each one against its neighbors.
  struct sc_candidate {
    string   seq_id;
    uint32_t position;
    int32_t  direction;
    uint32_t clipped_count;
    uint32_t total_count;
    uint32_t agree_count;
    string   consensus_tail;
    double   score;
    double   frequency;
    double   consensus_fraction;
    bool     suppressed;

    sc_candidate() : position(0), direction(0), clipped_count(0), total_count(0),
                     agree_count(0), score(0.0), frequency(0.0), consensus_fraction(0.0),
                     suppressed(false) {}

    static bool by_seq_direction_position(const sc_candidate& a, const sc_candidate& b) {
      if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
      if (a.direction != b.direction) return a.direction < b.direction;
      return a.position < b.position;
    }
  };
  vector<sc_candidate> candidates;

  string line;

  // Parameter line. Guards against a counts file left over from a run with different
  // tabulation settings, which the error-count done-file would otherwise let through.
  getline(in, line);
  ASSERT(!line.empty() && (line[0] == '#'),
         "Soft-clipping counts file is missing its parameter header and is probably from an older run:\n  "
         + _settings.soft_clipping_counts_file_name
         + "\nDelete 07_error_calibration/error_counts.done and re-run to regenerate it.");
  {
    string expected = "#sc_format=2\tsoft_clipping_minimum_bases=" + to_string(_settings.soft_clipping_minimum_bases);
    ASSERT(line.substr(0, expected.size()) == expected,
           "Soft-clipping counts file was tabulated with different settings than the current run:\n  "
           + line + "\nDelete 07_error_calibration/error_counts.done and re-run to regenerate it.");
  }

  getline(in, line); // column header

  // Checked only after the header, so that a stale counts file reports the actionable
  // error above rather than silently producing no evidence here.
  ASSERT(p0 > 0.0,
         "Soft-clipping null rate is zero. The summary in 07_error_calibration is probably from an "
         "older run:\n  " + _settings.soft_clipping_summary_file_name
         + "\nDelete 07_error_calibration/error_counts.done and re-run to regenerate it.");

  while (getline(in, line)) {
    vector<string> fields = split(line, "\t");
    ASSERT(fields.size() >= 7,
           "Soft-clipping counts file has too few columns and is probably from an older run:\n  "
           + _settings.soft_clipping_counts_file_name
           + "\nDelete 07_error_calibration/error_counts.done and re-run to regenerate it.");

    string seq_id = fields[0];
    uint32_t position = from_string<uint32_t>(fields[1]);
    int32_t direction = from_string<int32_t>(fields[2]);
    uint32_t clipped_count = from_string<uint32_t>(fields[3]);
    uint32_t total_count = from_string<uint32_t>(fields[4]);
    uint32_t agree_count = from_string<uint32_t>(fields[5]);
    string consensus_tail = fields[6];

    if (clipped_count == 0 || total_count == 0) continue;

    // The quantity tested is the number of reads clipped here *with the same tail*.
    // Reads clipped for uninteresting reasons disagree with each other, so this both
    // sharpens the test and keeps such positions from inflating the null.
    uint32_t test_count = agree_count;

    //// TIER 1: hard discard. These entries never enter the genome diff at all.

    if ((_settings.soft_clipping_minimum_read_count > 0) &&
        (test_count < _settings.soft_clipping_minimum_read_count)) continue;
    if (test_count == 0) continue;

    double p_value;
    if (use_beta_binomial) {
      p_value = beta_binomial_sf(static_cast<double>(test_count),
                                 static_cast<double>(total_count),
                                 alpha, beta);
    } else {
      p_value = bdtrc(static_cast<double>(test_count - 1),
                      static_cast<double>(total_count),
                      p0);
    }
    double log10_p_value = (p_value > 0.0) ? log10(p_value) : -300.0;
    double score = -(log10_p_value + log10_n_tests);

    // Expected to occur by chance somewhere in the genome: hopeless, and keeping these
    // would swamp the marginal evidence table. Mirrors the early-out in
    // test_RA_evidence_CONSENSUS_mode().
    if (score < 0.0) continue;

    sc_candidate c;
    c.seq_id             = seq_id;
    c.position           = position;
    c.direction          = direction;
    c.clipped_count      = clipped_count;
    c.total_count        = total_count;
    c.agree_count        = agree_count;
    c.consensus_tail     = consensus_tail;
    c.score              = score;
    c.frequency          = static_cast<double>(test_count) / static_cast<double>(total_count);
    c.consensus_fraction = static_cast<double>(agree_count) / static_cast<double>(clipped_count);
    candidates.push_back(c);
  }

  /*
   * Local non-maximum suppression.
   *
   * The exact base at which an aligner places a clip is ambiguous when the reference and the
   * donor share sequence at the junction, so the reads supporting one breakpoint smear over
   * several adjacent positions. Measured on ADP1: position 2151230 carries 290 clipped reads
   * and its neighbor 2151231 carries 12 whose consensus tail (TCAATACTCCTT) is one of
   * 2151230's donor groups (CTCAATACTCCT) shifted by a single base -- the same reads, clipped
   * one over. 2151232/33/34 hold a further 2-3 reads each. Only the strongest position in such
   * a run is a real event; the rest are shadows of it.
   *
   * Suppression is per (seq_id, direction). It must NOT reach across directions: a mobile
   * element insertion produces a -1 and a +1 junction a target-site duplication apart (ADP1
   * 288908/288909 one base apart, 2311060/2311061, 600249/600253), and both are genuine.
   *
   * The ranking key is the score, not the raw clipped count. At 2151229-31 the highest clipped
   * count belongs to 2151230, whose reads carry seven different donors and which is already
   * discarded above for having no consensus; ranking by clipped count would let it suppress
   * the genuine call at 2151229 (71 reads, 93% consensus, score 150). Ranking by score keeps
   * 2151229 and removes only the shadow at 2151231.
   */
  {
    const int32_t window = static_cast<int32_t>(_settings.soft_clipping_minimum_bases);

    // Sort by (seq_id, direction, position) so each run of nearby same-direction candidates
    // is contiguous and the sweep below only has to look forward.
    sort(candidates.begin(), candidates.end(), sc_candidate::by_seq_direction_position);

    for (size_t i = 0; i < candidates.size(); i++) {
      for (size_t j = i + 1; j < candidates.size(); j++) {
        if (candidates[j].seq_id != candidates[i].seq_id) break;
        if (candidates[j].direction != candidates[i].direction) break;
        if (static_cast<int32_t>(candidates[j].position) -
            static_cast<int32_t>(candidates[i].position) > window) break;

        // Strictly better score wins; ties go to the lower position so the result does not
        // depend on input order.
        if (candidates[j].score > candidates[i].score) candidates[i].suppressed = true;
        else                                           candidates[j].suppressed = true;
      }
    }
  }

  for (vector<sc_candidate>::const_iterator it = candidates.begin(); it != candidates.end(); it++) {
    const sc_candidate& c = *it;

    cDiffEntry sc_entry(SC);
    sc_entry[SEQ_ID]           = c.seq_id;
    sc_entry[POSITION]         = to_string(c.position);
    sc_entry[STRAND]           = to_string(c.direction);
    sc_entry[SC_READ_COUNT]    = to_string(c.clipped_count);
    sc_entry[SC_AGREE_COUNT]   = to_string(c.agree_count);
    sc_entry[SC_TOTAL_COUNT]   = to_string(c.total_count);
    // Frequency of the event actually tested: consensus-supporting clipped reads over
    // all reads reaching the position (total_count = clipped + read-through).
    sc_entry[FREQUENCY]        = formatted_double(c.frequency, 4).to_string();
    sc_entry[SC_CONSENSUS_FRACTION] = formatted_double(c.consensus_fraction, 4).to_string();
    if (c.consensus_tail != ".") sc_entry[SC_CONSENSUS_TAIL] = c.consensus_tail;
    sc_entry[SC_LOG10_E_VALUE] = formatted_double(c.score, kMutationScorePrecision).to_string();

    //// TIER 2: soft reject. The entry is kept and shown as marginal evidence.

    if (c.suppressed) {
      sc_entry.add_reject_reason("NEARBY_BETTER_SOFT_CLIPPING");
    }
    if ((_settings.soft_clipping_log10_e_value_cutoff > 0.0) &&
        (c.score < _settings.soft_clipping_log10_e_value_cutoff)) {
      sc_entry.add_reject_reason("SCORE_CUTOFF");
    }
    // Only a frequency that is too LOW is a reason to reject: an unusually high clipped
    // fraction is the signal itself, not a problem. (Unlike RA, which rejects at both ends.)
    if ((_settings.soft_clipping_frequency_cutoff > 0.0) &&
        (c.frequency < _settings.soft_clipping_frequency_cutoff - _settings.polymorphism_precision_decimal)) {
      sc_entry.add_reject_reason("FREQUENCY_BELOW_CUTOFF");
    }
    if ((_settings.soft_clipping_consensus_base_fraction > 0.0) &&
        (_settings.soft_clipping_consensus_fraction_cutoff > 0.0) &&
        (c.consensus_fraction < _settings.soft_clipping_consensus_fraction_cutoff - _settings.polymorphism_precision_decimal)) {
      sc_entry.add_reject_reason("CLIPPED_TAIL_CONSENSUS");
    }

    _gd.add(sc_entry);
  }
}

/*! Destructor.
 */
identify_mutations_pileup::~identify_mutations_pileup()
{
}

/*! Called for each reference genome position.
 */
void identify_mutations_pileup::pileup_callback(const pileup& p) {
  
  bool verbose = false;
  ASSERT(p.target() < _seq_info.size(), "Unknown target id: " + to_string<uint32_t>(p.target()));
  if (verbose) cout << "Target id: " << p.target() << " position: " << p.position_1() << endl;

  _this_deletion_propagation_cutoff = _deletion_propagation_cutoffs[p.target()];
  _this_deletion_seed_cutoff = _deletion_seed_cutoffs[p.target()];
  
  // if the propagation cutoff is zero then the coverage distribution failed
  if (_this_deletion_propagation_cutoff < 0.0) return;
  
  uint32_t position = p.position_1();
  
	int32_t insert_count=-1;
	bool next_insert_count_exists=true;
	
  // User RA evidence
  // We need to continue looping in some cases where there is user evidence RA
  // but no reads actually have a base at that insert position. Here's how:
  int32_t force_insert_count_max = 0;
  if (_user_evidence_ra_list.size()) {
    
    diff_entry_list_t::iterator user_ra = _user_evidence_ra_list.begin();
    
    while ( (user_ra != _user_evidence_ra_list.end()) && (from_string<uint32_t>((**user_ra)[POSITION]) == position)) {
      force_insert_count_max = from_string<int32_t>((**user_ra)[INSERT_POSITION]);
      user_ra++;
    }
  }
  
  // BIG outer loop to traverse this position AND any positions after this one
  // that are inserted relative to the reference genome
	while(next_insert_count_exists || (insert_count < force_insert_count_max)) {
		++insert_count;
    next_insert_count_exists = false;
		
		base_char ref_base_char = '.';
		if(!insert_count) {
			ref_base_char = p.reference_base_char_1(position);
		}
		
		//## zero out the info about this position
		position_base_info pos_info;
    position_base_info redundant_pos_info;
		
		//## keep track of coverage for deletion prediction
		position_coverage this_position_coverage;
		bool this_position_unique_only_coverage=true;
		
    //## reset SNP caller
    _snp_caller.reset(basechar2index(ref_base_char));
        
		//## polymorphism prediction data
		vector<polymorphism_data> pdata;
    
		//## for each alignment within this pileup:
		for(pileup::const_iterator i=p.begin(); i!=p.end(); ++i) {

      //## Discordant-pair (DP) region detection: incremental "enter" step.
      //## Add each discordant read to its paired read group's sliding window exactly ONCE, at its
      //## leftmost reference column (so it is counted a single time as it enters the window). Done
      //## here (piggybacking on this per-read loop) to avoid a second per-column scan; guarded to
      //## the first insert slot and to the read's start column.
      if ((insert_count == 0) && !_dp_groups.empty() && (i->reference_start_1() == position)) {
        if (i->is_paired() && !i->unmapped() && !i->proper_pair() &&
            !(i->flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FMUNMAP))) {
          map<uint32_t, int>::const_iterator gi = _fastq_index_to_dp_group.find(i->fastq_file_index());
          if (gi != _fastq_index_to_dp_group.end()) {
            const dp_group& g = _dp_groups[gi->second];

            // Bin by (focal-read strand) x (pair orientation). Strand: 0 = forward, 1 = reverse.
            int s = i->reversed() ? 1 : 0;
            // Orientation from the XP tag written by mark_pair_info: FR / RF / FF (RR folded to FF).
            string orientation;
            int o;
            if (!i->aux_get_Z("XP", orientation)) continue;  // no orientation -> not a countable DP read
            if      (orientation == "FR") o = 0;
            else if (orientation == "RF") o = 1;
            else if (orientation == "FF") o = 2;
            else continue;                                    // unexpected orientation -> skip
            int bin = s * 3 + o;

            // Build the condensed read-pair key: <read1_name>__<read2_name>__<r1_insert_size>.
            // breseq names mates <file_index>:<read_num> sharing read_num (fastq.cpp), and never sets
            // BAM_FREAD1/2, so R1/R2 role comes from fastq_file_index() vs the group's R1/R2 file ids.
            bool focal_is_r1 = (i->fastq_file_index() == g.r1_m_id);
            string fname = i->read_name();
            size_t colon = fname.find(':');
            string read1_name, read2_name;
            if (colon != string::npos) {
              string counter = fname.substr(colon + 1);
              read1_name = to_string(g.r1_m_id + 1) + ":" + counter;
              read2_name = to_string(g.r2_m_id + 1) + ":" + counter;
            } else {
              // Fallback for unexpected naming: use the focal name for its own slot only.
              read1_name = focal_is_r1 ? fname : "";
              read2_name = focal_is_r1 ? "" : fname;
            }
            // Insert size from R1's perspective (identical for both mates): negate when focal is R2.
            int32_t r1_insert_size = focal_is_r1 ? i->insert_size() : -i->insert_size();
            string key = read1_name + "__" + read2_name + "__" + to_string(r1_insert_size);

            dp_read dr;
            dr.start_pos = position;
            dr.key = key;
            // A tie-broken redundant discordant pair keeps X1>1 on its marked alignment (see
            // resolve_alignments.cpp); flag its DP side redundant so it propagates to DP evidence.
            dr.redundant = (i->redundancy() > 1);
            _dp_groups[gi->second].reads[bin].push_back(dr);
            ++_dp_metric[bin];
            // If a region in this bin is already open (from a previous column), this newly entering
            // read belongs to it, so append its key now. Reads entering on the column that OPENS a
            // region are instead captured by the open-time snapshot in check_discordant_completion
            // (which runs after this loop), so they are never double-counted.
            if (_dp_region_start[bin] != UNDEFINED_UINT32) {
              _dp_region_descriptors[bin].push_back(key);
              if (dr.redundant) ++_dp_region_redundant_count[bin];
            }
          }
        }
      }

      //## After these substitutions...
      //## Indel is -1 if the ref base is deleted in the read,
      //## Zero if the read base is aligned to a ref base, and
      //## Positive if the read base is an insertion relative to the ref base
      int indel=i->indel();
      if(indel < 0) {
        indel = 0;
      }
      if(i->is_del()) {
        indel = -1;
      }
      
      base_bam read_base_bam='.';
      bool on_insert_position_past_base = true; 
      if(indel >= insert_count) {
        read_base_bam = i->read_base_bam_0(i->query_position_0() + insert_count);
        on_insert_position_past_base = false; 
      }
            
      // Don't count bases without qualities!! -- not safe to even count them for coverage
      if(_base_bam_is_N(read_base_bam)) continue;
      
      //## gather information about the aligned base
      int32_t redundancy = i->redundancy();
      int32_t fastq_file_index = i->fastq_file_index();
      int strand = i->strand();
      bool trimmed = i->is_trimmed(on_insert_position_past_base);
      
      //##### update coverage if this is not a deletion in read relative to reference
      //### note that we count trimmed reads here, but not when looking for short indel mutations...	
      
      if(redundancy == 1) {
        //## keep track of unique coverage	
        ++this_position_coverage.unique[1+strand];
        
        //## we don't continue to consider further insertions relative
        //## to the reference unless uniquely aligned reads have them 
        if(indel > insert_count)
          next_insert_count_exists = true;
        
      } else {		
        //## mark that this position has some non-unique coverage
        this_position_unique_only_coverage = false;
        //## keep track of redundant coverage
        this_position_coverage.redundant[1+strand] += 1.0/redundancy;
        ++this_position_coverage.raw_redundant[1+strand];
        
        if (!trimmed) 
          ++redundant_pos_info[basebam2char(read_base_bam)][1+strand];
      }
      
      
      // When predicting base substitutions and short indels...

      // 1. Don't use information from redundantly mapped reads
      if (redundancy > 1) continue;
      
      // 2. Don't use information from trimmed bases in mapped reads
      if (trimmed) continue;

      // Right now this happens during resolution of read mappings
      // Later: Re-implement this here, along with other guards that use saved alignment values that can be easily looked up
      // 3. Don't use reads with low mapping quality (if requested in settings)
      //if (i->mapping_quality() < _settings.minimum_mapping_quality) continue;

			//##### deal with base calls
      //cerr << "POSITION:" << position << endl;
      
      covariate_values_t cv; 

      bool is_ok = _error_table.alignment_position_to_covariates(*i, insert_count, cv);
      //cv.obs_base is still not a char here...
      
      if (is_ok)  {
      
        if (cv.quality() < _settings.base_quality_cutoff) {
          //cout << cv.quality()  << endl;
          continue;
        }

        //## this is the coverage for SNP counts, tabulate AFTER skipping trimmed reads
        ++pos_info[baseindex2char(cv.obs_base())][1+strand];
        
        //##### this is for polymorphism prediction and making strings
        pdata.push_back(polymorphism_data(baseindex2char(cv.obs_base()),cv.quality(),i->strand(),cv.read_set(), cv));
        
        //cerr << " " << cv.obs_base() << " " << (char)ref_base << endl;

        _snp_caller.update(cv, strand == 1, i->mapping_quality(), _error_table);
      }
		} // end for-each read
		
    
    //#############################
		//## PER POSITION/INSERT COUNT
    //#############################
    
		//#sum up coverage observations	
		this_position_coverage.sum();
				
		//#we are trying to find the base with the most support
    cSNPCall snp_call = _snp_caller.get_prediction();
    
    base_char best_base_char('N');
    double consensus_bonferroni_score(numeric_limits<double>::quiet_NaN());
    double polymorphism_bonferroni_score(numeric_limits<double>::quiet_NaN());

    // SNP caller returns one genotype
    best_base_char = snp_call.genotype[0];
    // Here is where we convert to a Bonferroni type correction for phred score from SNP caller
    consensus_bonferroni_score = snp_call.score - (_log10_ref_length);
    
    //Do we predict a base at this position?
    bool base_predicted=false;
    if(consensus_bonferroni_score >= _consensus_score_cutoff) {
      base_predicted = true;
    }
    
    //cout << position << " " << best_base_char << " " << e_value_call << " " << (base_predicted ? "T" : "F") << endl;
    
		int total_cov[3]={0,0,0}; // triple, same as above
    
    //// BEGIN Print per-position output file
		ostringstream line;
		if (_print_per_position_file) {
      line << position << " " << insert_count << " " << ref_base_char << " " << consensus_bonferroni_score;
		}
    
    for(size_t j=0; j<base_list_size; ++j) {
			double top_cov = pos_info[base_char_list[j]][2];
			double bot_cov = pos_info[base_char_list[j]][0];
			total_cov[2] += (int)round(top_cov);
			total_cov[0] += (int)round(bot_cov);
    }
    
    //// Summing coverage here should be moved to where coverage is updated?
    
    // Debug: print additional information to file.
    if (_print_per_position_file) {
       
      // Print unique bases
      for(size_t j=0; j<base_list_size; ++j) {
        double top_cov = pos_info[base_char_list[j]][2];
        double bot_cov = pos_info[base_char_list[j]][0];
        if (_print_per_position_file) {
          line << " " << base_char_list[j] << " (" << bot_cov << "/" << top_cov << ")";
        }
      }
    
      // Print redundant bases
      for(size_t j=0; j<base_list_size; ++j) {
        double top_cov = redundant_pos_info[base_char_list[j]][2];
        double bot_cov = redundant_pos_info[base_char_list[j]][0];
        if (_print_per_position_file) {
          line << " r" << base_char_list[j] << " (" << bot_cov << "/" << top_cov << ")";
        }
      }
    }
    
    // Finally print line (kept separate from above because the line is
    // added to at various points in the code).
    if (_print_per_position_file) {
      _per_position_file << line.str() << endl;
    }
    //// END Per-position output file
		
		//###
		//## DELETION DELETION DELETION
		//###
		
    if(!_settings.skip_missing_coverage_prediction && (insert_count == 0))
      check_deletion_completion(p.target(), position, ref_base_char, this_position_coverage, consensus_bonferroni_score);

		//###
		//## DISCORDANT PAIR DISCORDANT PAIR DISCORDANT PAIR
		//###

    // Runs once per reference column (including empty columns) so the sliding window advances and
    // regions close correctly even where there is no coverage.
    if(!_dp_groups.empty() && (insert_count == 0))
      check_discordant_completion(p.target(), position);

		//###
		//## POLYMORPHISM POLYMORPHISM POLYMORPHISM
		//###								

		bool passed_as_polymorphism_prediction(false);
    polymorphism_prediction ppred;
    base_char second_best_base_char('N');

    cDiffEntry mut(RA);
    
    //## evaluate whether to call an actual mutation!
    // -- note that we accept > 0 and only reject later
    // so that these can potentially make it into the marginal data
    bool passed_as_consensus_prediction = (best_base_char != ref_base_char) && (!std::isnan(consensus_bonferroni_score) && (consensus_bonferroni_score > 0));
    
    // Find the bases with the highest and second highest coverage
    // We only predict polymorphisms involving these 'major' and 'minor' alleles
    base_index best_base_index(base_list_N_index);
    int best_base_coverage(0);
    base_index second_best_base_index(base_list_N_index);
    int second_best_base_coverage(0);

    vector<double> snp_probs = _snp_caller.get_genotype_log10_probabilities();
            
    for (uint8_t i=0; i<base_list_size; i++) {
      base_char this_base_char = base_char_list[i];
      base_index this_base_index = i;
      int this_base_coverage = pos_info[this_base_char][0] + pos_info[this_base_char][2];
      
      if (this_base_coverage==0) continue;
      
      // if better coverage or tied in coverage and better probability
      if ((this_base_coverage > best_base_coverage) ||
        ((this_base_coverage == best_base_coverage) && ((best_base_index==base_list_N_index) || (snp_probs[this_base_index] > snp_probs[best_base_index])))) {
        second_best_base_index = best_base_index;
        second_best_base_coverage = best_base_coverage;
        best_base_index = this_base_index;
        best_base_coverage = this_base_coverage;
      }
      else if ((this_base_coverage > second_best_base_coverage) 
        || ((this_base_coverage == second_best_base_coverage) && (((second_best_base_index==base_list_N_index) ||snp_probs[this_base_index] > snp_probs[second_best_base_index])))) {
        second_best_base_index = this_base_index;
        second_best_base_coverage = this_base_coverage;
      }
    }
    
    int this_base_coverage = min(pos_info[best_base_char][0], pos_info[best_base_char][2]);
    
    // Only try mixed SNP model if there is coverage for more than one base!
    ppred.frequency = 1;
    if (second_best_base_index != base_list_N_index) {
      best_base_char = base_char_list[best_base_index];
      second_best_base_char = base_char_list[second_best_base_index];

      // tries all frequencies of the best two
      if (_settings.polymorphism_prediction) {
        ppred = predict_polymorphism(best_base_char, second_best_base_char, pdata);
      }
      
      // tries only the raw ML frequency of the best two
      else /* if (_settings.mixed_base_prediction) */ {
        ppred = predict_mixed_base(best_base_char, second_best_base_char, pdata);
      }
      
      // Calculate E-value for polymorphism score (E theta)
      polymorphism_bonferroni_score = (-(log(ppred.likelihood_ratio_test_p_value)/log(10)) - _log10_ref_length);
      
      // Do we accept this as a polymorphism?
      if (polymorphism_bonferroni_score >= _polymorphism_score_cutoff)
        passed_as_polymorphism_prediction = true;
      
      //cerr << ppred.frequency << " " << ppred.log10_base_likelihood << " " << ppred.p_value << endl;
    }
		
		//###
		//## UNKNOWN UNKNOWN UNKNOWN
		//###
		if(insert_count == 0) {
			update_unknown_intervals(position, p.target(), base_predicted, this_position_unique_only_coverage);
		}
    
    //###
		//## Does any RA evidence pass tests?
		//###
				
    // This is a pointer to the added mutation, so that we can change it later for user-RA
    diff_entry_ptr_t added_mut_p(NULL);
    if (passed_as_consensus_prediction || passed_as_polymorphism_prediction) {
      
      //###
      //## Create new RA evidence mutation for genome diff
      //###		
      
      //## Fields common to consensus mutations and polymorphisms
      mut[SEQ_ID] = p.target_name();
      mut[POSITION] = to_string<uint32_t>(position);
      mut[INSERT_POSITION] = to_string<uint32_t>(insert_count);
      // Genotype quality is for the top called genotype
      mut[CONSENSUS_SCORE] = formatted_double(consensus_bonferroni_score, kMutationScorePrecision).to_string();
      mut[POLYMORPHISM_SCORE] = formatted_double(polymorphism_bonferroni_score, kMutationScorePrecision).to_string();
      
      //## Specific initializations for polymorphisms. Must take precedence.

      //# the frequency returned is the probability of the FIRST base
      //# we want to quote the probability of the second base (the minor allele from the reference).
      mut[REF_BASE] = ref_base_char;
      mut[NEW_BASE] = (best_base_char == ref_base_char) ? second_best_base_char : best_base_char;
      
      mut[MAJOR_BASE] = best_base_char;
      mut[MINOR_BASE] = second_best_base_char;
      mut[MAJOR_FREQUENCY] = formatted_double(ppred.frequency, _polymorphism_precision_places, true).to_string();
      
      double variant_frequency = ppred.frequency;
      if (mut[REF_BASE] == mut[MAJOR_BASE]) {
        variant_frequency = 1.0 - variant_frequency;
      }
      mut[POLYMORPHISM_FREQUENCY] = formatted_double(variant_frequency, _polymorphism_precision_places, true).to_string();


      // Add line to the polymorphism statistics input file if we are only a polymorphism

      if (ppred.frequency != 1) {
        annotate_polymorphism_statistics(mut, best_base_char, second_best_base_char, pos_info, pdata);
      }

      //## More fields common to consensus mutations and polymorphisms
      //## ...now that ref_base and new_base are defined
      
      vector<uint32_t>& ref_cov = pos_info[from_string<base_char>(mut[REF_BASE])];
      mut[REF_COV] = to_string(make_pair(static_cast<int32_t>(ref_cov[2]), static_cast<int32_t>(ref_cov[0])));
      
      vector<uint32_t>& new_cov = pos_info[from_string<base_char>(mut[NEW_BASE])];
      mut[NEW_COV] = to_string(make_pair(static_cast<int32_t>(new_cov[2]), static_cast<int32_t>(new_cov[0])));
      
      vector<uint32_t>& major_cov = pos_info[from_string<base_char>(mut[MAJOR_BASE])];
      mut[MAJOR_COV] = to_string(make_pair(static_cast<int32_t>(major_cov[2]), static_cast<int32_t>(major_cov[0])));
      
      vector<uint32_t>& minor_cov = pos_info[from_string<base_char>(mut[MINOR_BASE])];
      mut[MINOR_COV] = to_string(make_pair(static_cast<int32_t>(minor_cov[2]), static_cast<int32_t>(minor_cov[0])));
      
      mut[TOTAL_COV] = to_string(make_pair(total_cov[2], total_cov[0]));
      
      _gd.add(mut);
      
      //cout << "Added:" << _gd.evidence_list().back()->as_string() << endl;
      
      // what's going on here? we may need to change a value latter,
      // and add added a copy not the current one
      added_mut_p = _gd.get_list().back();
    } // END ra_output
    
    // Now we print additional RA items as user= if they have not already been printed.
    // for this position and insert_count
    while (_user_evidence_ra_list.size()
           && (((*(_user_evidence_ra_list.front()))[SEQ_ID]) == p.target_name())
           && (from_string<uint32_t>((*(_user_evidence_ra_list.front()))[POSITION]) == position)
           && (from_string<int32_t>((*(_user_evidence_ra_list.front()))[INSERT_POSITION]) == insert_count)) {
      
      // Note! mut only exists as a valid record for comparison
      // Must check everything here so multiple user_evidence items added at same position can pass
      if ( (passed_as_consensus_prediction || passed_as_polymorphism_prediction) && (*added_mut_p == *(_user_evidence_ra_list.front())))  {
        
        //cout << "Comparing:" << added_mut_p->as_string() << endl;
        //cout << "       to:" << (_user_evidence_ra_list.front())->as_string() << endl;
        
        //  cout << "FOUND MATCH in user evidence. User evidence not added" << endl;
          
        // in this case we have already created a valid RA item
        (*added_mut_p)[USER_DEFINED] = "1";
          
        // do not mark it as user_defined or there will be problems:
        // (it will not show up in HTML and it will mess up polymorphism stats)
        
      } else {
        
        //cout << "User evidence for position:" << (_user_evidence_ra_list.front().get())->as_string() << endl;
        //cout << "NO MATCH to new mutation prediction. Added user evidence entry" << endl;
        
        // copy the input entry and fill in more information
        cDiffEntry mut = *(_user_evidence_ra_list.front().get());
        mut.to_spec(); //remove additional fields that might be left over!!
        mut[USER_DEFINED] = "1";
        
        // These are already assigned correctly by copy of RA
        //mut[SEQ_ID] = p.target_name();
        //mut[POSITION] = to_string<uint32_t>(position);
        //mut[INSERT_POSITION] = to_string<uint32_t>(insert_count);

        // tries all frequencies of the best two
        
        base_char best_base_char = from_string<base_char>(mut[REF_BASE]);
        base_char second_best_base_char = from_string<base_char>(mut[NEW_BASE]);
        
        if (_settings.polymorphism_prediction) {
          ppred = predict_polymorphism(best_base_char, second_best_base_char, pdata);
        }
        // tries only the raw ML frequency of the best two
        else {
          ppred = predict_mixed_base(best_base_char, second_best_base_char, pdata);
        }
        
        double polymorphism_bonferroni_score = 10 * (-(log(ppred.likelihood_ratio_test_p_value)/log(10)) - _log10_ref_length);
        
        // These defs of major and minor base are not quite always accurate.
        // because they only take into account ref_base and new_base of the RA line
        // (and there could be a 3rd base that matters)
        // ---> in the future it may be better to KEEP the major and minor alleles and work with
        //      these instead, but that leads to its own issues
        mut[MAJOR_BASE] = (ppred.frequency > 0.5) ? best_base_char : second_best_base_char;
        mut[MINOR_BASE] = (ppred.frequency > 0.5) ? second_best_base_char : best_base_char;
        mut[MAJOR_FREQUENCY] = formatted_double( ((ppred.frequency > 0.5) ? ppred.frequency : 1 - ppred.frequency), _polymorphism_precision_places, true).to_string();
        
        // The frequency of the variant is always this, due to the way the ref_base is set as best_base_char
        mut[POLYMORPHISM_FREQUENCY] = formatted_double(1 - ppred.frequency, _polymorphism_precision_places, true).to_string();
        
        // Consensus mode
        if (!_settings.polymorphism_prediction) {
          if (from_string<double>(mut[POLYMORPHISM_FREQUENCY]) > 0.5 ) {
            mut[FREQUENCY] = "1";
          } else {
            mut[FREQUENCY] = "0";
          }
        } else { // Polymorphism mode
          mut[FREQUENCY] =  mut[POLYMORPHISM_FREQUENCY];
        }
        
        // Genotype quality is for the top called genotype
        mut[CONSENSUS_SCORE] = formatted_double(consensus_bonferroni_score, kMutationScorePrecision).to_string();
        mut[POLYMORPHISM_SCORE] = formatted_double(polymorphism_bonferroni_score, kMutationScorePrecision).to_string();
        
        vector<uint32_t>& ref_cov = pos_info[from_string<base_char>(mut[REF_BASE])];
        mut[REF_COV] = to_string(make_pair(static_cast<int32_t>(ref_cov[2]), static_cast<int32_t>(ref_cov[0])));
        
        vector<uint32_t>& new_cov = pos_info[from_string<base_char>(mut[NEW_BASE])];
        mut[NEW_COV] = to_string(make_pair(static_cast<int32_t>(new_cov[2]), static_cast<int32_t>(new_cov[0])));
        
        vector<uint32_t>& major_cov = pos_info[from_string<base_char>(mut[MAJOR_BASE])];
        mut[MAJOR_COV] = to_string(make_pair(static_cast<int>(major_cov[2]), static_cast<int>(major_cov[0])));
        
        vector<uint32_t>& minor_cov = pos_info[from_string<base_char>(mut[MINOR_BASE])];
        mut[MINOR_COV] = to_string(make_pair(static_cast<int>(minor_cov[2]), static_cast<int>(minor_cov[0])));
        
        mut[TOTAL_COV] = to_string(make_pair(total_cov[2], total_cov[0]));
        
        _gd.add(mut);
        
        //cout << "Added:" << _gd.evidence_list().back()->as_string() << endl;
      }
      
      _user_evidence_ra_list.pop_front();
    }
    
	}
}

/*! Called at the beginning of a reference sequence fragment
    Open per-reference files
*/

void identify_mutations_pileup::at_target_start(const uint32_t tid)
{
    
  // Open per-reference coverage file:
	if(_print_coverage_data) {
		string filename = _settings.file_name(_settings.complete_coverage_text_file_name, "@", target_name(tid));
		_coverage_data.open(filename.c_str());
    ASSERT(!_coverage_data.fail(), "Could not open output file:" + filename);
    _coverage_data << "position" << "\t";
    _coverage_data << "ref_base" << "\t";
    _coverage_data << "unique_top_cov" << "\t";
    _coverage_data << "unique_bot_cov" << "\t";
    _coverage_data << "redundant_top_cov" << "\t";
    _coverage_data << "redundant_bot_cov" << "\t";
    _coverage_data << "raw_redundant_top_cov" << "\t";
    _coverage_data << "raw_redundant_bot_cov" << endl;
	}
  
  // Reset the Missing Coverage evidence variables
  _last_deletion_start_position = UNDEFINED_UINT32;
	_last_deletion_end_position = UNDEFINED_UINT32;
	_last_deletion_redundant_start_position = UNDEFINED_UINT32;
	_last_deletion_redundant_end_position = UNDEFINED_UINT32;
	_last_start_unknown_interval = UNDEFINED_UINT32;

  // Reset the Discordant Pair (DP) evidence variables (regions do not span reference boundaries)
  for (int bin = 0; bin < kDPnBins; bin++) {
    for (vector<dp_group>::iterator g = _dp_groups.begin(); g != _dp_groups.end(); g++) g->reads[bin].clear();
    _dp_metric[bin] = 0;
    _dp_region_start[bin] = UNDEFINED_UINT32;
    _dp_region_max_count[bin] = 0;
    _dp_region_descriptors[bin].clear();
  }

}
  
/*! Called at the end of a reference sequence fragment
    Close per-reference files
 */
void identify_mutations_pileup::at_target_end(const uint32_t tid) {

  // end "open" Missing Coverahge and Unknown intervals
  if (!_settings.skip_missing_coverage_prediction) {
    check_deletion_completion(tid, target_length(tid)+1, '.', position_coverage(numeric_limits<double>::quiet_NaN()), numeric_limits<double>::quiet_NaN());
  }
  update_unknown_intervals(target_length(tid)+1, tid, true, false);

  // Flush any open discordant-pair (DP) region at the end of this reference sequence.
  if (!_dp_groups.empty()) {
    check_discordant_completion(tid, target_length(tid)+1);
  }

  // if this target failed to have its coverage fit, mark the entire thing as a deletion
  double _this_deletion_propagation_cutoff = _deletion_propagation_cutoffs[tid];
  // if the propagation cutoff is -1 then the coverage distribution failed
  if (!_settings.skip_missing_coverage_prediction && (_this_deletion_propagation_cutoff < 0.0))
  {
    cDiffEntry del(MC);
    del[SEQ_ID] = target_name(tid);
    
    del[START] = to_string<uint32_t>(1);
    del[END] = to_string<uint32_t>(target_length(tid));
    del[START_RANGE] = to_string<uint32_t>(0);
    del[END_RANGE] = to_string<uint32_t>(0);
    
    del[LEFT_OUTSIDE_COV] = "NA";
    del[LEFT_INSIDE_COV] = formatted_double(0.0, 0).to_string();
    del[RIGHT_INSIDE_COV] = formatted_double(0.0, 0).to_string();
    del[RIGHT_OUTSIDE_COV] = "NA";

    _gd.add(del);
  }
  
  // Close per-reference coverage file:
	if(_print_coverage_data)
		_coverage_data.close();
}


/*! Helper method to track information about putative deleted regions.
 
 Used at each pileup iteration and at the end.
 //## when called at the end of a fragment, the position is fragment_length+1
 //## and this_position_coverage is undefined
 
 @JEB This function expects 1-indexed positions!!!
 
 */
void identify_mutations_pileup::check_deletion_completion(uint32_t seq_id, uint32_t position, char ref_base_char, const position_coverage& this_position_coverage, double e_value_call) {

	//cerr << position << " " << e_value_call << endl;
	
  // special case = beginning of new seq_id
  if (position == 1) 
    _last_position_coverage = position_coverage(numeric_limits<double>::quiet_NaN());
  
  // print to optional output file
  if (!std::isnan(this_position_coverage.unique[1]) && _coverage_data.is_open()) {
    _coverage_data << position << "\t";
    _coverage_data << ref_base_char << "\t";
    _coverage_data << this_position_coverage.unique[2] << "\t"; // top
    _coverage_data << this_position_coverage.unique[0] << "\t"; // bottom
    _coverage_data << this_position_coverage.redundant[2] << "\t"; // top
    _coverage_data << this_position_coverage.redundant[0] << "\t"; // bottom
    _coverage_data << this_position_coverage.raw_redundant[2] << "\t"; // top
    _coverage_data  << this_position_coverage.raw_redundant[0] << endl; // bottom
  }
  
  
  //## UNIQUE COVERAGE
  //#start a new possible deletion if we fall below the propagation cutoff
  if(this_position_coverage.unique[1] <= _this_deletion_propagation_cutoff) {
    if(_last_deletion_start_position == UNDEFINED_UINT32) {
      _last_deletion_start_position = position;
      _left_outside_coverage_item = _last_position_coverage;
      _left_inside_coverage_item = this_position_coverage;
    }
  }
		
  //##keep track of whether we've encountered the seed value
  if(!std::isnan(this_position_coverage.unique[1]) && (this_position_coverage.total <= _this_deletion_seed_cutoff)) {
    _this_deletion_reaches_seed_value = true;
  }
	
	//## If we are in a deletion and rise back above the propagation cutoff OR we are at the end of this fragment (NAN),
  //## then record the current deletion.
	if( (_last_deletion_start_position != UNDEFINED_UINT32) && 
     ( std::isnan(this_position_coverage.unique[1]) || (this_position_coverage.unique[1] > _this_deletion_propagation_cutoff) ) )
  {
		
		if(_this_deletion_reaches_seed_value) {

      _last_deletion_end_position = position-1;
      if (_last_deletion_redundant_end_position == UNDEFINED_UINT32) 
        _last_deletion_redundant_end_position = _last_deletion_end_position;
      if (_last_deletion_redundant_start_position == UNDEFINED_UINT32) 
        _last_deletion_redundant_start_position = _last_deletion_start_position;

      cDiffEntry del(MC);
			del[SEQ_ID] = target_name(seq_id);
			del[START] = to_string<uint32_t>(_last_deletion_start_position);
			del[END] = to_string<uint32_t>(_last_deletion_end_position);
			del[START_RANGE] = to_string<uint32_t>(_last_deletion_redundant_start_position - _last_deletion_start_position);
			del[END_RANGE] = to_string<uint32_t>(_last_deletion_end_position - _last_deletion_redundant_end_position);
			
      del[LEFT_OUTSIDE_COV] = formatted_double(_left_outside_coverage_item.unique[1], 0).to_string();
      del[LEFT_INSIDE_COV] = formatted_double(_left_inside_coverage_item.unique[1], 0).to_string();
      del[RIGHT_INSIDE_COV] = formatted_double(_last_position_coverage.unique[1], 0).to_string();
      del[RIGHT_OUTSIDE_COV] = formatted_double(this_position_coverage.unique[1], 0).to_string();
      
			_gd.add(del);
		}
		
		//#reset the search
		_this_deletion_reaches_seed_value = false;
		_this_deletion_redundant_reached_zero = false;
		_last_deletion_start_position = UNDEFINED_UINT32;
		_last_deletion_end_position = UNDEFINED_UINT32;
		_last_deletion_redundant_start_position = UNDEFINED_UINT32;
		_last_deletion_redundant_end_position = UNDEFINED_UINT32;
	}
  
  //## REDUNDANT COVERAGE
  //## updated only if we are still within a deletion
  if (_last_deletion_start_position != UNDEFINED_UINT32) {
    
    if (this_position_coverage.redundant[1] == 0) {
      _this_deletion_redundant_reached_zero = true;
      _last_deletion_redundant_end_position = UNDEFINED_UINT32;
    }
    else if (this_position_coverage.redundant[1] > 0) {
      //## if there is any redundant coverage remember the start (until we find zero redundant coverage)
      if (!_this_deletion_redundant_reached_zero) {
        _last_deletion_redundant_start_position = position;
      }
      else if (_last_deletion_redundant_end_position == UNDEFINED_UINT32) {
        _last_deletion_redundant_end_position = position;
      }
    }
  }
	
  
  _last_position_coverage = this_position_coverage;
}


/*! Helper method to track discordant-pair (DP) candidate regions.

 Called once per reference column (including empty columns), and once past the end of each reference
 sequence (position = target_length+1) to flush an open region. Maintains, per paired read group, a
 sliding window of width (median + 2.42*MAD) over the reference start positions of discordant reads.
 Reads are added once at their start column ("enter", in pileup_callback) and removed once here when
 they fall out of the window ("exit"). A candidate region is a maximal run of columns where the total
 in-window discordant count reaches --discordant-pair-seed.

 @JEB expects 1-indexed positions.
 */
void identify_mutations_pileup::check_discordant_completion(uint32_t seq_id, uint32_t position)
{
  static const char* kOrientName[3] = { "FR", "RF", "FF" };

  // Run an independent detector for each (focal strand x orientation) bin, so a breakpoint's
  // forward/reverse shoulders and its distinct orientations each become separate regions.
  for (int bin = 0; bin < kDPnBins; bin++) {

    // "Exit" step: drop discordant reads that have fallen out of each group's bin window.
    for (vector<dp_group>::iterator g = _dp_groups.begin(); g != _dp_groups.end(); g++) {
      while (!g->reads[bin].empty() &&
             (static_cast<double>(position) - static_cast<double>(g->reads[bin].front().start_pos) >= g->window_width)) {
        g->reads[bin].pop_front();
        --_dp_metric[bin];
      }
    }

    // Never "above" once past the end of the sequence: closes an open region at the flush sentinel
    // and prevents opening a spurious region that would never close.
    bool above = (_dp_metric[bin] >= _dp_seed) && (position <= target_length(seq_id));

    if (above) {
      if (_dp_region_start[bin] == UNDEFINED_UINT32) {
        // Open a new region; snapshot the bin's reads currently in-window as its initial keys.
        _dp_region_start[bin] = position;
        _dp_region_max_count[bin] = static_cast<uint32_t>(_dp_metric[bin]);
        _dp_region_descriptors[bin].clear();
        _dp_region_redundant_count[bin] = 0;
        for (vector<dp_group>::iterator g = _dp_groups.begin(); g != _dp_groups.end(); g++) {
          for (deque<dp_read>::iterator r = g->reads[bin].begin(); r != g->reads[bin].end(); r++) {
            _dp_region_descriptors[bin].push_back(r->key);
            if (r->redundant) ++_dp_region_redundant_count[bin];
          }
        }
      } else if (static_cast<uint32_t>(_dp_metric[bin]) > _dp_region_max_count[bin]) {
        _dp_region_max_count[bin] = static_cast<uint32_t>(_dp_metric[bin]);
      }
    } else if (_dp_region_start[bin] != UNDEFINED_UINT32) {
      // Close the open region (its last in-region column was position-1).
      dp_candidate_region region;
      region.seq_id = target_name(seq_id);
      region.start = _dp_region_start[bin];
      region.end = position - 1;
      region.strand = (bin / 3 == 0) ? 'F' : 'R';
      region.orientation = kOrientName[bin % 3];
      region.max_count = _dp_region_max_count[bin];
      // Redundant side: a majority of the region's discordant reads mapped redundantly (i.e. the
      // reads pile up here as a tie-broken multicopy element side). Propagated to DP evidence's
      // side_N_redundant flag, mirroring how JC marks a redundant junction side.
      region.redundant = (2 * _dp_region_redundant_count[bin] > _dp_region_descriptors[bin].size());
      string joined;
      for (size_t k = 0; k < _dp_region_descriptors[bin].size(); k++) {
        if (k) joined += ";";
        joined += _dp_region_descriptors[bin][k];
      }
      region.discordant_pairs = joined;
      _dp_candidate_regions.push_back(region);

      _dp_region_start[bin] = UNDEFINED_UINT32;
      _dp_region_max_count[bin] = 0;
      _dp_region_descriptors[bin].clear();
      _dp_region_redundant_count[bin] = 0;
    }
  }
}


/*! Write all accumulated DP candidate regions to a CSV in one pass (after the pileup completes). */
void identify_mutations_pileup::write_dp_candidate_regions(const string& filename)
{
  // Sort by (seq_id, start, strand, orientation) so subregions of one breakpoint appear adjacent.
  sort(_dp_candidate_regions.begin(), _dp_candidate_regions.end(),
       [](const dp_candidate_region& a, const dp_candidate_region& b) {
         if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
         if (a.start != b.start) return a.start < b.start;
         if (a.strand != b.strand) return a.strand < b.strand;
         return a.orientation < b.orientation;
       });

  ofstream out(filename.c_str());
  ASSERT(!out.fail(), "Could not open output file: " + filename);
  // 'redundant' is a new column BEFORE 'discordant_pairs' (which must stay last: it is ';'-joined
  // and parsed as the final field by dp_evidence). 1 = a majority of this region's reads were
  // redundant (tie-broken multicopy side); 0 otherwise.
  out << "seq_id,start,end,strand,orientation,length,max_discordant_count,redundant,discordant_pairs" << endl;
  for (vector<dp_candidate_region>::const_iterator r = _dp_candidate_regions.begin(); r != _dp_candidate_regions.end(); r++) {
    out << r->seq_id << "," << r->start << "," << r->end << "," << r->strand << "," << r->orientation << ","
        << (r->end - r->start + 1) << "," << r->max_count << "," << (r->redundant ? 1 : 0) << "," << r->discordant_pairs << endl;
  }
  out.close();
}


/*! Helper method to track unknowns.
 */
void identify_mutations_pileup::update_unknown_intervals(uint32_t position, uint32_t seq_id, bool base_predicted, bool this_position_unique_only_coverage)
{
  //debug
  /*
	cerr << position << " " << base_predicted << " " << this_position_unique_only_coverage << endl;
	if(_last_start_unknown_interval != UNDEFINED_UINT32) {
		cerr << _last_start_unknown_interval << endl;
	} else {
		cerr << "undef" << endl;
	}
  */
	
	if(!base_predicted) {
		if(this_position_unique_only_coverage) {
			++s.coverage_unique_uncalled;
		}
		if(_last_start_unknown_interval == UNDEFINED_UINT32) {
			_last_start_unknown_interval = position;
		}
	}	else {
		if(this_position_unique_only_coverage) {
			++s.coverage_unique_called;
		}
			
		//#end interval where we were unable to call mutations
		if(_last_start_unknown_interval != UNDEFINED_UINT32) {
      cDiffEntry new_interval(UN);
			new_interval[SEQ_ID] = target_name(seq_id);
			new_interval[START] = to_string<uint32_t>(_last_start_unknown_interval);
			new_interval[END] = to_string<uint32_t>(position - 1);
			_gd.add(new_interval);
			
			_last_start_unknown_interval = UNDEFINED_UINT32;
		}
	}
}
  
void identify_mutations_pileup::annotate_polymorphism_statistics(cDiffEntry& mut, char best_base_char, char second_best_base_char, position_base_info& pos_info, const vector<polymorphism_data>& pdata)
{
  uint32_t major_top_strand = pos_info[best_base_char][2];
  uint32_t major_bot_strand = pos_info[best_base_char][0];
  uint32_t minor_top_strand = pos_info[second_best_base_char][2];
  uint32_t minor_bot_strand = pos_info[second_best_base_char][0];

  vector<double> major_quals;
  vector<double> minor_quals;
  for (vector<polymorphism_data>::const_iterator it = pdata.begin(); it != pdata.end(); ++it) {
    if (it->_base_char == best_base_char)
      major_quals.push_back(static_cast<double>(it->_quality));
    if (it->_base_char == second_best_base_char)
      minor_quals.push_back(static_cast<double>(it->_quality));
  }

  double ks_quality_p_value = 1.0;
  if (!major_quals.empty() && !minor_quals.empty())
    ks_quality_p_value = ks_test_two_sample_less(minor_quals, major_quals);

  double fisher_strand_p_value = fisher_exact_test_2x2(minor_top_strand, minor_bot_strand, major_top_strand, major_bot_strand);

  double combined_log = -2.0 * (log(ks_quality_p_value) + log(fisher_strand_p_value));
  double bias_p_value = isinf(combined_log) ? 0.0 : incompletegamma(2.0, combined_log / 2.0, true);
  double bias_e_value = bias_p_value * static_cast<double>(_total_reference_length);

  mut["ks_quality_p_value"]    = formatted_double(ks_quality_p_value,    5, true).to_string();
  mut["fisher_strand_p_value"] = formatted_double(fisher_strand_p_value, 5, true).to_string();
  mut["bias_p_value"]          = formatted_double(bias_p_value,          5, true).to_string();
  mut["bias_e_value"]          = formatted_double(bias_e_value,          5, true).to_string();
}


/*! Predict the significance of putative polymorphisms.
 */
polymorphism_prediction identify_mutations_pileup::predict_polymorphism(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata ) {
  
  //#calculate the likelihood of observed reads given this position is 100% the best base
	double log10_likelihood_of_one_base_model = 0;
  for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
  
    double log10_correct_pr;

    covariate_values_t this_cv = it->_cv;
      
    if(it->_strand == 1) {
      this_cv.ref_base() = basechar2index(best_base_char);
    } else {
      this_cv.ref_base() = basechar2index(complement_base_char(best_base_char)); 
      this_cv.obs_base() = complement_base_index(this_cv.obs_base());
    }
    log10_correct_pr = _error_table.get_log10_prob(this_cv);
       
    log10_likelihood_of_one_base_model  += log10_correct_pr;
  }

	vector<uint8_t> best_base_qualities;
	vector<uint8_t> second_best_base_qualities;
	uint32_t best_base_strand_hash[] = {0, 0};
	uint32_t second_best_base_strand_hash[] = {0, 0};
  
  for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
  
    int8_t zp_strand = (it->_strand == +1) ? 1 : 0;
    if (it->_base_char == best_base_char) {
      best_base_qualities.push_back(it->_quality);
      best_base_strand_hash[zp_strand]++;
    }
    else if (it->_base_char == second_best_base_char) {
      second_best_base_qualities.push_back(it->_quality);
      second_best_base_strand_hash[zp_strand]++;
    }
  }
  
	//## Maximum likelihood of observing alignment if sequenced bases were a mixture of the top two bases  
  pair<double,double> best_two_base_model = best_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata);
  double max_likelihood_fr_first_base = best_two_base_model.first;
  double log10_likelihood_of_two_base_model = best_two_base_model.second;
    
  //## Likelihood ratio test
  double log10_likelihood_difference = log10_likelihood_of_one_base_model - log10_likelihood_of_two_base_model;

  //debug output 
  /*
  cerr  << "ML Best Base Fraction: " << max_likelihood_fr_first_base << endl;
  cerr  << " Log10 Likelihood (one base model): " << log10_likelihood_of_one_base_model << endl;
  cerr  << " Log10 Likelihood (two base model): " << log10_likelihood_of_two_base_model << endl;
  cerr  << " Log10 Likelihood (different): " << log10_likelihood_difference << endl;
  */
    
  long double p_value = 1;
  if (max_likelihood_fr_first_base != 1.0) {
    double likelihood_ratio_test_value = -2*log(10)*log10_likelihood_difference;
    
    p_value = pchisq(1.0L, likelihood_ratio_test_value);
    //cerr << "likelihood_ratio_test_value: " << likelihood_ratio_test_value << " p-value: " << p_value << endl;
  }

  //debug output 
  /*
  cerr
    << " Log10 Likelihood Difference (one vs two base model): " << (log10_likelihood_of_one_base_model - log10_likelihood_of_two_base_model) 
    << " P-value: " << p_value 
    << endl;
  */
   
  polymorphism_prediction p(max_likelihood_fr_first_base, log10_likelihood_of_one_base_model - log10_likelihood_of_two_base_model, p_value);
		
	return p;
}
  
/*! Predict the significance of putative polymorphisms.
 */
polymorphism_prediction identify_mutations_pileup::predict_mixed_base(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata ) {
  
  //#calculate the likelihood of observed reads given this position is 100% the best base  
  double log10_likelihood_of_one_base_model = 0;
  for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
    
    double log10_correct_pr;
    
    covariate_values_t this_cv = it->_cv;
    
    if(it->_strand == 1) {
      this_cv.ref_base() = basechar2index(best_base_char);
    } else {
      this_cv.ref_base() = basechar2index(complement_base_char(best_base_char)); 
      this_cv.obs_base() = complement_base_index(this_cv.obs_base());
    }
    log10_correct_pr = _error_table.get_log10_prob(this_cv);
    
    log10_likelihood_of_one_base_model  += log10_correct_pr;
  }
  
  vector<uint8_t> best_base_qualities;
  vector<uint8_t> second_best_base_qualities;
  uint32_t best_base_strand_hash[] = {0, 0};
  uint32_t second_best_base_strand_hash[] = {0, 0};
  
  for(vector<polymorphism_data>::iterator it=pdata.begin(); it<pdata.end(); ++it) {
    
    int8_t zp_strand = (it->_strand == +1) ? 1 : 0;

    if (it->_base_char == best_base_char) {
      best_base_qualities.push_back(it->_quality);
      best_base_strand_hash[zp_strand]++;
    }
    else if (it->_base_char == second_best_base_char) {
      second_best_base_qualities.push_back(it->_quality);
      second_best_base_strand_hash[zp_strand]++;
    }
  }
  
  
  
  // Unlike full polymorphism prediction, we test just the raw frequency of the two bases
  // and do not check for bias later.
  
  // This would calculate the true best frequency -- slower, but more accurate
  // pair<double,double> best_two_base_model = best_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata);
  
  double max_likelihood_fr_first_base = static_cast<double>(best_base_strand_hash[0] + best_base_strand_hash[1]) 
    / static_cast<double>(best_base_strand_hash[0] + best_base_strand_hash[1] + second_best_base_strand_hash[0] + second_best_base_strand_hash[1]);
  
  double log10_likelihood_of_two_base_model = calculate_two_base_model_log10_likelihood(
                                                                                        best_base_char, 
                                                                                        second_best_base_char, 
                                                                                        pdata, 
                                                                                        max_likelihood_fr_first_base
                                                                                        );
  
  //## Likelihood ratio test
  double log10_likelihood_difference = log10_likelihood_of_one_base_model - log10_likelihood_of_two_base_model;
  
  //debug output 
  /*
   cerr  << "ML Best Base Fraction: " << max_likelihood_fr_first_base << endl;
   cerr  << " Log10 Likelihood (one base model): " << log10_likelihood_of_one_base_model << endl;
   cerr  << " Log10 Likelihood (two base model): " << log10_likelihood_of_two_base_model << endl;
   cerr  << " Log10 Likelihood (different): " << log10_likelihood_difference << endl;
   */
  
  long double p_value = 1;
  if (max_likelihood_fr_first_base != 1.0) {
    double likelihood_ratio_test_value = -2*log(10)*log10_likelihood_difference;
    
    p_value = pchisq(1.0L, likelihood_ratio_test_value);
    //cerr << "likelihood_ratio_test_value: " << likelihood_ratio_test_value << " p-value: " << p_value << endl;
  }
  
  polymorphism_prediction p(max_likelihood_fr_first_base, log10_likelihood_of_one_base_model - log10_likelihood_of_two_base_model, p_value);
  
  return p;
}
 
  
double identify_mutations_pileup::slope_at_percentage_best_base(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata, const double guess, const double precision, double& middle_point_log10_likelihood)  
{
  // precision is a fraction of the value
  double point_1 = max(guess * (1 - precision), 0.0);  
  double point_2 = min(guess * (1 + precision), 1.0);
  
  double point_1_log10_likelihood = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, point_1);
  double point_2_log10_likelihood = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, point_2);
  
  middle_point_log10_likelihood = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, guess);

  // Check for local minimum!
  if ( (middle_point_log10_likelihood < point_1_log10_likelihood) && (middle_point_log10_likelihood < point_2_log10_likelihood) ) {
    return 0;
  }
  
  return (point_2_log10_likelihood - point_1_log10_likelihood) / (point_2 - point_1);
}

/*! Find the best fraction for the best base at a polymorphic site.
 */
  
pair<double,double> identify_mutations_pileup::best_two_base_model_log10_likelihood(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata)
{	
  uint32_t iterations = 0;
  double current_upper_pr_first_base = 1.0;
  double current_lower_pr_first_base = 0.0;
  
  double current_upper_pr_first_base_log10_likelihood = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, current_upper_pr_first_base);
  double current_lower_pr_first_base_log10_likelihood = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, current_lower_pr_first_base);
  
  // precision is a fraction of the value...
  while (current_upper_pr_first_base - current_lower_pr_first_base > (current_upper_pr_first_base + current_lower_pr_first_base) / 2 * _polymorphism_precision_decimal) {
    
    iterations++;
    //cout << iterations << " "  << current_lower_pr_first_base << " " << current_upper_pr_first_base << endl;
    
    double current_middle_pr_first_base = (current_upper_pr_first_base + current_lower_pr_first_base) / 2;
    
    double current_middle_pr_first_base_log10_likelihood;
    double middle_slope = slope_at_percentage_best_base(best_base_char, second_best_base_char, pdata, current_middle_pr_first_base, _polymorphism_precision_decimal, current_middle_pr_first_base_log10_likelihood);
    
    // Slope is set to zero if the tested point is better than the ones 
    // on either side, when calculating the slope.
    if ( middle_slope == 0) {
      return make_pair(current_middle_pr_first_base, current_middle_pr_first_base_log10_likelihood);
    
    } else if ( middle_slope < 0) {
      
      current_upper_pr_first_base = current_middle_pr_first_base;
      current_upper_pr_first_base_log10_likelihood = current_middle_pr_first_base_log10_likelihood;
      
    } else {
      
      current_lower_pr_first_base = current_middle_pr_first_base;
      current_lower_pr_first_base_log10_likelihood = current_middle_pr_first_base_log10_likelihood;
    }
  }  
    
  if (current_lower_pr_first_base_log10_likelihood > current_upper_pr_first_base_log10_likelihood) {

    return make_pair(current_lower_pr_first_base, current_lower_pr_first_base_log10_likelihood);
  }
  
  return make_pair(current_upper_pr_first_base, current_upper_pr_first_base_log10_likelihood);
}
  
/*
pair<double,double> identify_mutations_pileup::best_two_base_model_log10_likelihood(base_char best_base_char, base_char second_best_base_char, vector<polymorphism_data>& pdata)
{	
  
	double cur_pr_first_base = 1;
	double cur_log_pr = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, cur_pr_first_base);

	double last_pr_first_base = 1;
	double last_log_pr = cur_log_pr;

	//print "$cur_pr_first_base $cur_log_pr\n" if ($verbose);

	while (cur_log_pr >= last_log_pr)
	{
		last_log_pr = cur_log_pr;
		last_pr_first_base = cur_pr_first_base;
    if (cur_pr_first_base < 0) break;

		cur_pr_first_base -= 0.001;
		cur_log_pr = calculate_two_base_model_log10_likelihood(best_base_char, second_best_base_char, pdata, cur_pr_first_base);
	}
	
	return make_pair(last_pr_first_base, last_log_pr);
}
 */

/*! Calculate the likelihood of a mixture model of two bases leading to the observed read bases.
 */
double identify_mutations_pileup::calculate_two_base_model_log10_likelihood(
                                                                            base_char best_base_char, 
                                                                            base_char second_best_base_char, 
                                                                            const vector<polymorphism_data>& pdata, 
                                                                            double best_base_freq
                                                                            )
{
	double log10_likelihood = 0;	
	
  for(vector<polymorphism_data>::const_iterator it=pdata.begin(); it<pdata.end(); ++it) {
  
    //## the first value is pr_base, second is pr_not_base
    double best_base_log10pr;
    double second_best_base_log10pr;
    
    covariate_values_t this_cv = it->_cv;
    
    if (it->_strand == -1) {
      this_cv.obs_base() = complement_base_index(this_cv.obs_base());
    }
    
    if(it->_strand == 1) {
      this_cv.ref_base() = basechar2index(best_base_char);
    } else {
      this_cv.ref_base() = basechar2index(complement_base_char(best_base_char));
    }
    best_base_log10pr = _error_table.get_log10_prob(this_cv);
    
    if(it->_strand == 1) {
      this_cv.ref_base() = basechar2index(second_best_base_char);
    } else {
      this_cv.ref_base() = basechar2index(complement_base_char(second_best_base_char));
    }
    second_best_base_log10pr = _error_table.get_log10_prob(this_cv);

    //debug output
    //cerr << "Base in Read: " << it->base << " Read Strand: " << it->strand << endl;
    //cerr << "Best Base: " << best_base << " Key: " << best_base_key << " Chance of Observing: " << pow(10,best_base_log10pr) << endl;
    //cerr << "Second Best Base: " << second_best_base << " Key: " << second_best_base_key << " Chance of Observing: " << pow(10,second_best_base_log10pr) << endl;

    double pr_ref_base_given_obs = best_base_freq * pow(10, best_base_log10pr) + (1-best_base_freq) * pow(10, second_best_base_log10pr);

    log10_likelihood += log(pr_ref_base_given_obs);		
  }

	log10_likelihood /= log(10);
  
  //debug output
  /*
  cerr << "Best Base: " << best_base << " Second Best Base: " << second_best_base << " Fraction Best Base: " << best_base_freq << " Log10 Likelihood " << log10_likelihood << endl;
  */
	return log10_likelihood;
}

cDiscreteSNPCaller::cDiscreteSNPCaller(
                                       const string& type,
                                       uint32_t reference_length
                                       ) 
: _type(type)
{
  

  // Uniform priors across all bases.
  if (_type == "haploid") {
    
    double uniform_probability = 1.0 / static_cast<double>(base_list_size);
    add_genotype("A", uniform_probability);
    add_genotype("C", uniform_probability);
    add_genotype("G", uniform_probability);
    add_genotype("T", uniform_probability);
    add_genotype(".", uniform_probability);
  }
  
  /*
  // Prior that we expect one change from reference
  else if (_type == "haploid-change") {
    
    // recall that first one counts as reference
    double uniform_probability = 1.0 / reference_length;    
    add_genotype("A", 1.0 - 4.0 * uniform_probability);
    add_genotype("C", uniform_probability);
    add_genotype("G", uniform_probability);
    add_genotype("T", uniform_probability);
    add_genotype(".", uniform_probability);
  }
  */
  
  // Extra states and priors for unexpected mixed states... experimental
  else if (_type == "haploid-cnv") {
    
    double mixed_probability = 1.0 / reference_length;    
    double uniform_probability = (1.0 - 10 * mixed_probability) / base_list_size;    
    
    add_genotype("A", uniform_probability);
    add_genotype("C", uniform_probability);
    add_genotype("G", uniform_probability);
    add_genotype("T", uniform_probability);
    add_genotype(".", uniform_probability);
    
    add_genotype("AC", mixed_probability);
    add_genotype("AG", mixed_probability);
    add_genotype("AT", mixed_probability);
    add_genotype("A.", mixed_probability);
    
    add_genotype("CG", mixed_probability);
    add_genotype("CT", mixed_probability);
    add_genotype("C.", mixed_probability);
    
    add_genotype("GT", mixed_probability);
    add_genotype("G.", mixed_probability);
    
    add_genotype("T.", mixed_probability);
  }  
  
  else
  {
    ERROR("Unknown SNP Caller type:" + type);
  }
  
  // Check priors
  double total_probability = 0;
  for(size_t i=0; i<_log10_genotype_prior_probabilities.size(); i++) {
    total_probability += pow(10, _log10_genotype_prior_probabilities[i]);
  }
  ostringstream ss;
  ss << setprecision(5) << total_probability;
  
  ASSERT( from_string<double>(ss.str()) == 1.0, "Prior probabilities do not sum to 1. (" + to_string(total_probability) + ").")
  
  reset(0);
}
  
void cDiscreteSNPCaller::add_genotype(const string& genotype, double probability) {
  
  _log10_genotype_prior_probabilities.push_back(log10(probability));
  
  vector<base_index> gv;
  for(size_t i=0; i<genotype.length(); i++) {
    gv.push_back(basechar2index(genotype[i]));
  }
  _genotype_vector.push_back(gv);
  
}
  
  
void cDiscreteSNPCaller::reset(uint8_t ref_base_index) {
  _best_genotype_index = 0;
  _observations = 0;
  _normalized_observations = 0;
  _log10_genotype_probabilities = _log10_genotype_prior_probabilities;
  
  (void) ref_base_index;
  //this is for where there are unbalanced priors -- haploid-change
  //they do not behave properly when the reference is 'N'
  //swap(_genotype_probability[0], _genotype_probability[ref_base_index]);
}
  
void cDiscreteSNPCaller::update(const covariate_values_t& cv, bool obs_top_strand, int32_t mapping_quality, cErrorTable& et) {

  covariate_values_t this_cv = cv;
  //update probabilities give observation using Bayes rule
  
  if (!obs_top_strand) {
    this_cv.obs_base() = complement_base_index(this_cv.obs_base()); 
  }
  
  double incorrect_mapping_prob = pow(10, -static_cast<double>(mapping_quality) / 10);
  double correct_mapping_prob = 1 - incorrect_mapping_prob;
  this->_normalized_observations += correct_mapping_prob;
  
  double total_prob = 0.0;
  for (uint32_t i=0; i<_genotype_vector.size(); i++) {
  
    vector<base_index>& gv = this->_genotype_vector[i];
    double this_pr = 0.0;
    
    for (uint32_t j=0; j < gv.size(); j++) {
      
      this_cv.ref_base() = gv[j];
      
      if (!obs_top_strand) {
        this_cv.ref_base() = complement_base_index(this_cv.ref_base()); 
      }
      this_pr += (correct_mapping_prob * et.get_prob(this_cv) + incorrect_mapping_prob * 1.0 / _genotype_vector.size()) * pow(10, this->_log10_genotype_probabilities[i]) / gv.size();
      
      // Floating point error can make this a very  negative number
      if (this_pr < 0.0) this_pr = 0.0;
    }
    
    total_prob += this_pr;
  }

  double highest_pr = -numeric_limits<double>::max();
  for (uint32_t i=0; i<this->_genotype_vector.size(); i++) {
    
    vector<base_index>& gv = this->_genotype_vector[i];
    double this_pr = 0.0;
    
    for (uint32_t j=0; j < gv.size(); j++) {
      
      this_cv.ref_base() = gv[j];
      
      if (!obs_top_strand) {
        this_cv.ref_base() = complement_base_index(this_cv.ref_base()); 
      }
      this_pr += (correct_mapping_prob * et.get_prob(this_cv) + incorrect_mapping_prob * 1 / _genotype_vector.size()) / gv.size();
    }
    
    this->_log10_genotype_probabilities[i] += log10(this_pr) - log10(total_prob);
    
    if (this->_log10_genotype_probabilities[i] > highest_pr) {
      this->_best_genotype_index = i;
      highest_pr = this->_log10_genotype_probabilities[i];
    }
  }
  
  //print();
  
  _observations++;
}
  
void cDiscreteSNPCaller::print() {
  
  for (uint32_t i=0; i<_genotype_vector.size(); i++) {
    
    vector<base_index>& gv = _genotype_vector[i];
    
    cout << "Genotype: " ;
    for (uint32_t j=0; j < gv.size(); j++) {
      
      cout << baseindex2char(gv[j]);
    }
    
    cout << " " << _log10_genotype_probabilities[i] << endl;
  }
}


cSNPCall cDiscreteSNPCaller::get_prediction()
{
  cSNPCall snp_call;
  
  if (_observations == 0) {
    //Best base is 'N' and E-value is NAN
    return snp_call;
  }
  
  // need to go through and find the most probable
  snp_call.genotype = "";
  for(size_t i=0; i<_genotype_vector[_best_genotype_index].size(); i++)
    snp_call.genotype += baseindex2char(_genotype_vector[_best_genotype_index][i]);
    
  // we want to normalize the probabilities, but avoid floating point errors    
  double log10_offset_probability = -numeric_limits<double>::max();
  for(size_t i=0; i < _log10_genotype_probabilities.size(); i++) {
    if (i != _best_genotype_index)
      log10_offset_probability = max(log10_offset_probability, _log10_genotype_probabilities[i]);
  }
  
  double total_error_probability = 0;
  for (uint32_t i=0; i<_genotype_vector.size()-1; i++) {
    if (i != _best_genotype_index)
      total_error_probability += pow(10, _log10_genotype_probabilities[i] - log10_offset_probability);
  }
  double log10_total_error_probability = log10(total_error_probability);
  log10_total_error_probability += log10_offset_probability;
  
  snp_call.score = (_log10_genotype_probabilities[_best_genotype_index] - log10_total_error_probability);
  
  return snp_call;
}

} // namespace breseq

