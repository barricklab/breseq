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
#include "stats.h"
#include "pd_evidence.h"
#include "coverage_output.h"

using namespace std;

namespace breseq {

//! Coverage columns of the table written for CNery, in emission order.
/*! Shared by the header, the per-read-group repeats and the rows, so the three cannot drift apart.
    These are the names 'bam2cov --total-only' writes (coverage_output.cpp), which is the schema
    CNery expects; keeping them identical is what lets it read either table with one code path.
 */
static const vector<string>& coverage_column_names()
{
  static const vector<string> names = make_vector<string>("unique_cov")("redundant_cov")("total_cov");
  return names;
}

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
                polymorphism_precision_decimal,
                polymorphism_precision_places,
								print_per_position_file
							);
	imp.do_pileup(settings.call_mutations_seq_id_set());
  if (settings.predict_soft_clipping) imp.add_sc_evidence(summary, ref_seq_info);
  // Write discordant-pair candidate regions accumulated during the pileup (single operation).
  imp.write_dp_candidate_regions(settings.dp_candidate_regions_file_name);
  // Likewise for missing-pair regions -- but only when they were actually collected, so a run
  // without --predict-missing-pairs leaves no empty CSV behind.
  if (settings.predict_missing_pairs) imp.write_mp_candidate_regions(settings.mp_candidate_regions_file_name);
  // Likewise for pair-distance regions.
  if (settings.predict_pair_distance) imp.write_pd_candidate_regions(settings.pd_candidate_regions_file_name);
  imp.write_gd(gd_file);
}

/*! Build the null distribution the pair-distance (PD) seed tests against, for one read group.

 Turns the shared length-bias-corrected histogram (see pd_load_weighted_histogram, which is where the
 correction and the deliberate absence of an upper truncation are explained) into a mid-CDF in
 kPDuScale fixed point, indexed by distance. Sharing the loader with pd_evidence is what makes the
 null a region was found under and the null it is then judged under the same by construction.

 The MID-CDF (G(d-1) + G(d))/2 rather than the CDF G(d): the distance is discrete, and E[G(D)] for a
 discrete D is 1/2 + E[P(D=d)]/2, strictly greater than 1/2. Using G directly would give every
 column a small positive bias and seed the LONG bin genome-wide.

 Also returns the mean gap of a covering pair, which is the mean of (d - trunc) under this same
 length-biased distribution. That is the distance over which the set of pairs covering a column turns
 over completely, so it is the correlation length of the seed statistic and hence what the effective
 number of independent tests is the reference length divided by. It is computed here because the
 weighted histogram is exactly what it is a moment of, and because the obvious alternative --
 (median distance - 2 * read_length) -- is not the same quantity and goes NEGATIVE on any library
 whose fragments are shorter than two reads, which silently disables the correction entirely.

 Returns false (leaving cdf empty, which disables PD for this group's reads) if the histogram is
 missing or carries no usable weight.
 */
static bool pd_build_cdf(const string& hist_file_name, int32_t trunc, int32_t max_distance,
                         vector<int64_t>& cdf, double& mean_covering_gap)
{
  cdf.clear();
  mean_covering_gap = 0.0;

  vector<double> weighted;
  double total = 0.0;
  if (!pd_load_weighted_histogram(hist_file_name, trunc, max_distance, weighted, total)) return false;

  cdf.resize(weighted.size(), 0);
  double below = 0.0;
  double gap_sum = 0.0;
  for (size_t d = 0; d < weighted.size(); d++) {
    double mid = (below + 0.5 * weighted[d]) / total;
    below += weighted[d];
    gap_sum += weighted[d] * (static_cast<double>(d) - static_cast<double>(trunc));
    int64_t scaled = static_cast<int64_t>(mid * static_cast<double>(identify_mutations_pileup::kPDuScale) + 0.5);
    if (scaled < 0) scaled = 0;
    if (scaled > identify_mutations_pileup::kPDuScale) scaled = identify_mutations_pileup::kPDuScale;
    cdf[d] = scaled;
  }
  mean_covering_gap = gap_sum / total;
  return true;
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
  
// The single RA score: log10 evidence that the position carries a non-reference allele at all.
//
// Evidence files written before the two scores were merged carry consensus_score and
// polymorphism_score instead. Those are still readable -- users re-run the Output step over an old
// evidence .gd routinely -- and the better-supported of the two stands in, since each was the
// admitting score for one of the two claims this one score now gates.
static double RA_score(const cDiffEntry& ra)
{
  if (ra.entry_exists(SCORE)) return double_from_string(ra.get(SCORE));

  double legacy = numeric_limits<double>::quiet_NaN();
  if (ra.entry_exists(CONSENSUS_SCORE)) legacy = double_from_string(ra.get(CONSENSUS_SCORE));
  if (ra.entry_exists(POLYMORPHISM_SCORE)) {
    double p = double_from_string(ra.get(POLYMORPHISM_SCORE));
    if (std::isnan(legacy) || (p > legacy)) legacy = p;
  }
  return legacy;
}

// ---------------------------------------------------------------------------------------------
// Exact (Clopper-Pearson) confidence bounds on the variant frequency at an RA position, so the
// frequency cutoffs are a confidence statement rather than a point-estimate comparison: a call is
// never rejected merely for having shallow coverage, only for being confidently on the wrong side.
//
// The bounds are computed where the model lives (write_RA_frequency_bounds, below) and read back
// from the entry here, so the interval that decides the call is the one the .gd and the report
// show. It is centred on the reported frequency and widened by the error-model effective depth,
// which is what a binomial "n" has to be once reads differ in how much they actually say: a
// quality-2 base and a quality-40 base are not one trial each.
//
// The fallback path reconstructs the old interval from total_cov -- raw read count as n -- for
// evidence written before the bounds were recorded. That is exactly the calibration-blind
// construction this replaced, so it is a compatibility shim and not a second opinion.
//
// Returns false if the position has no usable depth, in which case callers fall back to the
// point estimate.
static bool RA_frequency_bounds(const cDiffEntry& ra, double& lower, double& upper)
{
  if (ra.entry_exists(FREQUENCY_LOWER) && ra.entry_exists(FREQUENCY_UPPER)) {
    lower = from_string<double>(ra.get(FREQUENCY_LOWER));
    upper = from_string<double>(ra.get(FREQUENCY_UPPER));
    return true;
  }

  if (!ra.entry_exists(TOTAL_COV) || !ra.entry_exists(FREQUENCY)) return false;
  vector<string> top_bot = split(ra.get(TOTAL_COV), "/");
  if (top_bot.size() < 2) return false;
  double n = from_string<double>(top_bot[0]) + from_string<double>(top_bot[1]);
  if (!(n > 0.0)) return false;
  double k = from_string<double>(ra.get(FREQUENCY)) * n;
  lower = binomial_frequency_lower_bound(k, n);
  upper = binomial_frequency_upper_bound(k, n);
  return true;
}

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
    if ( (from_string<double>(top_bot[0]) < settings.consensus_minimum_total_coverage_each_strand) ||
        (from_string<double>(top_bot[1]) < settings.consensus_minimum_total_coverage_each_strand) ) {
      ra.add_reject_reason("TOTAL_STRAND_COVERAGE");
      rejected = true;
    }
  }
  
  if (settings.consensus_minimum_variant_coverage > 0) {
    vector<string> top_bot = split(ra[MAJOR_COV], "/");
    if ( from_string<double>(top_bot[0]) + from_string<double>(top_bot[1]) < settings.consensus_minimum_variant_coverage ) {
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
// MAJORITY allele. That gives two flat attempts, in order, each gated on the score and on the 95%
// lower confidence bound of the variant frequency:
//
//   1. consensus     -- score clears mutation_log10_e_value_cutoff AND
//                       L >= consensus_frequency_cutoff (0.50 by default: confidently a majority)
//   2. polymorphism  -- otherwise, score clears polymorphism_log10_e_value_cutoff AND
//                       L >= polymorphism_frequency_cutoff (0.10 by default: confidently present)
//   3. otherwise the position is dropped.
//
// Both tests read the LOWER bound, so a position is never accepted on thin coverage and never
// rejected merely for having thin coverage -- only for failing to establish the claim being made.
//
// DO NOT "unify" this with polymorphism mode, which reads the UPPER bound in the same slot. The
// asymmetry is the whole point: the two modes carry opposite defaults. Here the fixed
// interpretation is the default, so the question is "is the variant confidently the majority" and
// the lower bound answers it. There the polymorphic interpretation is the default, so the question
// is "can we still rule out that it is fixed" and the upper bound answers that.
//
// The temptation to unify comes from an argument that no longer applies. Under the exact
// Clopper-Pearson bounds this used to use, the identity L(k of k) = alpha^(1/k) meant a lower-bound
// test against 0.50 demanded ln(0.05)/ln(0.50) = 4.32 variant reads even at 100% frequency, so
// genuine fixed variants were dropped for want of depth. Profile-likelihood bounds do not have that
// floor -- at two decisive reads the lower bound is already 0.505 -- and measured over the test
// goldens, 0 of 218 consensus-mode fixed variants now fail this test. Meanwhile switching this slot
// to the upper bound would relabel two well-measured ~48% variants (at 115x and 292x coverage) as
// fixed 100% mutations, which is exactly what the lower bound is here to prevent.
//
// The Fisher strand-bias and variant-strand-minimum filters are ON by default; the K-S quality-bias
// and homopolymer filters are off. They are applied within whichever attempt they belong to so
// those command-line options keep working.
bool test_RA_evidence_CONSENSUS_mode(
                                     cDiffEntry& ra,
                                     cReferenceSequences& ref_seq_info,
                                     const Settings& settings
                                     )
{
  double score = RA_score(ra);

  // Falls back to the point estimate if the entry carries neither recorded bounds nor usable depth.
  double lower, upper;
  if (!RA_frequency_bounds(ra, lower, upper))
    lower = upper = from_string<double>(ra[FREQUENCY]);

  /////////////////////////////////
  // 1. Is it the majority allele?
  /////////////////////////////////

  if (score < settings.mutation_log10_e_value_cutoff) {
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
    // Delete if we are just the reference base!
    return (ra[REF_BASE] == ra[MAJOR_BASE]);
  }
  ra[CONSENSUS_REJECT] = ra[REJECT];
  ra.clear_reject_reasons();

  /////////////////////////////////
  // 2. Is it present at all?
  /////////////////////////////////

  if (score < settings.polymorphism_log10_e_value_cutoff) {
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
//   1. consensus     -- score clears mutation_log10_e_value_cutoff AND
//                       U >= consensus_frequency_cutoff (0.95 by default). The upper bound, not the
//                       lower: a call snaps to a fixed 100% difference only when we cannot rule out
//                       that it IS fixed. Under consensus mode the same slot asks "is it confidently
//                       the majority" with the lower bound, because there the fixed interpretation
//                       is the default. See the note there before trying to make these agree.
//   2. polymorphism  -- otherwise, score clears polymorphism_log10_e_value_cutoff AND
//                       L >= polymorphism_frequency_cutoff (0.05 by default: confidently present)
//   3. otherwise the position is KEPT and marked rejected -- unlike consensus mode, which drops it.
//
// The Fisher strand-bias and variant-strand-minimum guards are ON by default here; the K-S
// quality-bias and homopolymer guards are off. Each sits in whichever attempt it belongs to.
bool test_RA_evidence_POLYMORPHISM_mode(
                                     cDiffEntry& ra,
                                     cReferenceSequences& ref_seq_info,
                                     const Settings& settings
                                     )
{
  double score = RA_score(ra);

  double lower, upper;
  if (!RA_frequency_bounds(ra, lower, upper))
    lower = upper = from_string<double>(ra[FREQUENCY]);

  /////////////////////////////////
  // 1. Can we still call it fixed?
  /////////////////////////////////

  if (score < settings.mutation_log10_e_value_cutoff) {
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
    // Delete if we are just the reference base!
    return (ra[REF_BASE] == ra[MAJOR_BASE]);
  }
  ra[CONSENSUS_REJECT] = ra[REJECT];
  ra.clear_reject_reasons();

  /////////////////////////////////
  // 2. Is it present at all?
  /////////////////////////////////

  if (score < settings.polymorphism_log10_e_value_cutoff) {
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
    return false;
  }

  /////////////////////////////////
  // 3. Neither claim holds -- keep it, marked rejected.
  /////////////////////////////////

  ra[POLYMORPHISM_REJECT] = ra[REJECT];
  ra.clear_reject_reasons();
  ra[PREDICTION] = "polymorphism";
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
  //  * SCORE                  log10 evidence that a non-reference allele is present at all
  //                           (or, on evidence written before the scores were merged, at least one
  //                            of CONSENSUS_SCORE / POLYMORPHISM_SCORE -- see RA_score())
  //  * FREQUENCY              the variant allele's fitted frequency, as a fraction of total depth
  //  * NEW_COV, REF_COV       present for all entries
  
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
    if (!ra.entry_exists(SCORE) && !ra.entry_exists(CONSENSUS_SCORE) && !ra.entry_exists(POLYMORPHISM_SCORE)) {
      ERROR("Expected field 'score' in evidence item\n" + ra.as_string());
    }
    if (!ra.entry_exists(FREQUENCY)) {
      ERROR("Expected field '" + cString(FREQUENCY) + "' in evidence item\n" + ra.as_string());
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
, _polymorphism_precision_decimal(polymorphism_precision_decimal)
, _polymorphism_precision_places(polymorphism_precision_places)
, _log10_ref_length(0)
, _total_reference_length(summary.sequence_conversion.total_reference_sequence_length)
, _this_deletion_reaches_seed_value(false)
, _this_deletion_redundant_reached_zero(false)
, _last_position_coverage_printed(0)
, _print_per_position_file(print_per_position_file)
, _dp_seed(settings.discordant_pair_seed)
, _mp_enabled(false)
, _mp_seed(settings.missing_pair_seed)
, _mp_seed_fraction(settings.missing_pair_seed_fraction)
, _mp_S_k(0.0)
, _mp_S_k2_over_n(0.0)
, _mp_S_n(0.0)
, _mp_S_n_trimmed(0.0)
, _mp_total_k(0)
, _mp_columns(0)
, _mp_columns_trimmed(0)
, _mp_window_width(0.0)
, _pd_enabled(false)
, _pd_seed(settings.pair_distance_seed)
, _pd_pairs_considered(0)
, _pd_pairs_ambiguous_placement(0)
, _pd_z_seed(0.0)
, _pd_max_span(0)
, _pd_ring_w(0)
, _pd_n(0)
, _pd_long(0)
, _pd_short(0)
, _pd_u(0)
, _pd_last_b(0)
, _pd_mean_covering_gap(0.0)
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

  // Initialize per-bin pair-distance region state (bin = shift direction; see kPDnBins).
  for (int bin = 0; bin < kPDnBins; bin++) {
    _pd_region_start[bin] = UNDEFINED_UINT32;
    _pd_region_peak_position[bin] = 0;
    _pd_region_peak_covering[bin] = 0;
    _pd_region_peak_tail[bin] = 0;
    _pd_region_peak_z[bin] = 0.0;
    for (int step = 0; step < kPDLadderSteps; step++) {
      _pd_ladder_open[bin][step] = false;
      _pd_ladder_count[bin][step] = 0;
    }
  }

  // Initialize per-bin missing-pair region state (bin = focal strand; see kMPnBins).
  for (int bin = 0; bin < kMPnBins; bin++) {
    _mp_metric[bin] = 0;
    _mp_total_metric[bin] = 0;
    _mp_region_start[bin] = UNDEFINED_UINT32;
    _mp_region_max_count[bin] = 0;
    _mp_region_max_total[bin] = 0;
    _mp_region_redundant_count[bin] = 0;
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

  // ...and does splitting that coverage by read group tell CNery anything? Unlike
  // 'bam2cov --per-read-group', which always writes at least an RG-0 set, this table gets the extra
  // columns only when the BAM declares two or more read groups. With one group they would be exact
  // duplicates of the aggregate columns -- pure width, no information -- so the whole per-group
  // accumulator is skipped in the (usual) single-group case.
  _print_read_group_coverage = _print_coverage_data && (read_groups().size() > 1);

  // load the error table file and convert back to probabilities
  _error_table.read_log10_prob_table(settings.error_rates_file_name);
  _error_table.log10_prob_to_prob();

  // The k_read_set covariate indexes a read FILE, but a record names its read GROUP -- which for a
  // paired set spans two files. Give the table this BAM's header groups so it can recombine the
  // group with the mate flag and land on the same read file the rates were fit for.
  _error_table.set_read_file_partition(&read_groups(),
                                       make_read_file_partition(read_groups(), settings.read_file_sets));
  
  if (_print_per_position_file) {
    _per_position_file.open(_settings.mutation_identification_per_position_file_name.c_str());
  }
  
  // load the list of user RA evidence
  if (_settings.user_evidence_genome_diff_file_name != "") {
    load_user_ra_evidence_from_gd();
  }

  // Set up discordant-pair (DP) candidate-region detection.
  // For each PAIRED read group with a valid paired-mapping-distance fit, build a sliding-window
  // group with width W = median + 2.42 * MAD, and map that read group's BAM index to it so we can
  // classify a read's group by its read_group_index() during the pileup.
  //
  // Keyed on the BAM's own read group index -- which is what the RG:Z: tag resolves to -- rather
  // than on cReadFile::m_id. m_id is assigned in command-line order and is never renumbered when
  // read files are dropped by conversion or --limit-fold-coverage, so it was a different index
  // space from the one the lookup used.
  const PairedMappingDistanceDistributionSummaries& pmdd = summary.preliminary_paired_mapping_distance_distribution;
  _read_group_to_dp_group.assign(settings.read_file_sets.size(), -1);
  for (size_t set_index = 0; set_index < settings.read_file_sets.size(); set_index++) {
    const cReadFileSet& rfs = settings.read_file_sets[set_index];
    if (!rfs.is_paired()) continue;
    PairedMappingDistanceDistributionSummaries::const_iterator it = pmdd.find(rfs.m_base_name);
    if (it == pmdd.end()) continue;
    if (it->second.median <= 0.0) continue;

    dp_group g;
    g.window_width = it->second.median + 2.42 * it->second.mad;
    g.mp_window_width = it->second.median;   // see dp_group::mp_window_width
    g.r1_read_name_prefix = rfs.m_files[0].m_id + 1;  // see dp_group for why this is m_id-derived
    g.r2_read_name_prefix = rfs.m_files[1].m_id + 1;
    int group_index = static_cast<int>(_dp_groups.size());
    _dp_groups.push_back(g);
    _read_group_to_dp_group[set_index] = group_index;
  }

  // Missing-pair (MP) detection rides on the same groups, so it can only run once they exist. Unlike
  // DP -- whose seeding is always on -- MP seeding is gated by its option, so a run without
  // --predict-missing-pairs pays exactly one predictable branch per read in the pileup loop.
  _mp_enabled = _settings.predict_missing_pairs && !_dp_groups.empty();

  // The width the MP null is tabulated at, and therefore the width mp_evidence must count over.
  // Taking the max across groups matches paired_library_params' convention for pair_median.
  for (vector<dp_group>::const_iterator g = _dp_groups.begin(); g != _dp_groups.end(); g++)
    if (g->mp_window_width > _mp_window_width) _mp_window_width = g->mp_window_width;

  // Set up pair-distance (PD) candidate-region detection. Same read groups again, and the same
  // option gating as MP: without --predict-pair-distance the pileup pays one predictable branch per
  // read. That gating matters more here than for DP or MP, because PD's qualifying pairs are the
  // ORDINARY concordant ones -- it cannot short-circuit on !proper_pair() the way DP does, so when
  // it is on it runs on essentially every read.
  if (_settings.predict_pair_distance && !_dp_groups.empty()) {

    int32_t read_length = static_cast<int32_t>(summary.sequence_conversion.read_length_avg + 0.5);
    int32_t trunc = 2 * read_length;
    double derived_span = 0.0;

    // Two passes, because the span bounds the null and so has to be known before any CDF is built --
    // see pd_load_weighted_histogram for what an unbounded null does to the seed test. pd_evidence
    // derives the span exactly this way, which is what keeps the null a region is FOUND under and the
    // null it is later JUDGED under identical.
    for (size_t set_index = 0; set_index < settings.read_file_sets.size(); set_index++) {
      const cReadFileSet& rfs = settings.read_file_sets[set_index];
      if (!rfs.is_paired()) continue;
      if (_read_group_to_dp_group[set_index] < 0) continue;
      PairedMappingDistanceDistributionSummaries::const_iterator it = pmdd.find(rfs.m_base_name);
      if (it == pmdd.end()) continue;

      derived_span = max(derived_span, 2.0 * it->second.distance_cutoff);
    }

    _pd_max_span = _settings.pair_distance_maximum_span;
    if (_pd_max_span <= 0) _pd_max_span = static_cast<int32_t>(derived_span + 0.5);

    // Matches the fallback in pd_evidence: a span no larger than the reads leaves nothing to bound.
    int32_t null_max_distance = (_pd_max_span > trunc) ? _pd_max_span : 0;

    _pd_z_hist.assign(kPDzHistBins, 0);

    int usable_groups = 0;
    for (size_t set_index = 0; set_index < settings.read_file_sets.size(); set_index++) {
      const cReadFileSet& rfs = settings.read_file_sets[set_index];
      if (!rfs.is_paired()) continue;
      int group_index = _read_group_to_dp_group[set_index];
      if (group_index < 0) continue;
      PairedMappingDistanceDistributionSummaries::const_iterator it = pmdd.find(rfs.m_base_name);
      if (it == pmdd.end()) continue;

      string hist_file_name = Settings::file_name(settings.paired_mapping_distance_histogram_file_name, "#", rfs.m_base_name);
      double group_gap = 0.0;
      if (!pd_build_cdf(hist_file_name, trunc, null_max_distance, _dp_groups[group_index].pd_cdf, group_gap)) continue;
      // The orientation the null was built from. Pairs of any other orientation are not scored
      // against it; see dp_group::pd_orientation.
      _dp_groups[group_index].pd_orientation = it->second.majority_orientation;
      _pd_mean_covering_gap = max(_pd_mean_covering_gap, group_gap);
      usable_groups++;
    }

    // The seed is now only a sensitivity filter -- what a PD call has to clear is the genome-wide
    // e-value computed in pd_evidence, which is calibrated from the region counts collected below.
    // So this threshold is set LOW deliberately, low enough that the regions it produces are
    // overwhelmingly noise, because that noise is the sample the calibration is fitted to. Setting it
    // where a call would be believable would leave nothing to measure the null with.
    //
    // The analytic derivation is kept as the upper bound: whatever else, there is no point seeding
    // where a normal null with the correct effective test count already expects many hits per genome.
    // The effective number of tests is (reference length / mean covering gap) -- one per block over
    // which the covering set turns over. Note this uses the measured covering gap and not
    // (median - 2*read_length), which is negative on short-insert libraries and used to make this
    // whole correction a no-op.
    _pd_z_seed = _settings.pair_distance_seed_z;
    if (_pd_z_seed <= 0.0) {
      double alpha_half = 0.005;
      if ((_pd_mean_covering_gap > 0.0) && (_total_reference_length > 0))
        alpha_half *= _pd_mean_covering_gap / static_cast<double>(_total_reference_length);
      if (alpha_half > 0.5) alpha_half = 0.5;
      if (alpha_half < 1e-12) alpha_half = 1e-12;
      double analytic_z = ndtri(1.0 - alpha_half);
      // Loose enough that the regions are overwhelmingly noise and there are several hundred of them
      // to fit a tail to, tight enough that the per-region BAM rescan stays affordable. On the 3.6 Mb
      // library this defect was found on it yields roughly 500 regions where the old threshold of 3.0
      // yielded 113. Never seed ABOVE the analytic threshold: there is nothing to learn out there.
      const double kPDSeedSensitivityFloor = 2.5;
      _pd_z_seed = min(analytic_z, kPDSeedSensitivityFloor);
      if (_pd_z_seed < 1.0) _pd_z_seed = 1.0;
    }

    _pd_enabled = (usable_groups > 0) && (_pd_max_span > 0);
    if (_pd_enabled) {
      _pd_ring_w = _pd_max_span + 2;
      _pd_ring_n.assign(_pd_ring_w, 0);
      _pd_ring_long.assign(_pd_ring_w, 0);
      _pd_ring_short.assign(_pd_ring_w, 0);
      _pd_ring_u.assign(_pd_ring_w, 0);
    } else {
      cerr << "WARNING: --predict-pair-distance was given, but no paired read group has a usable" << endl;
      cerr << "         mapping-distance histogram. No pair-distance (PD) evidence will be predicted." << endl;
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

// Size the genome-wide spanning strand split is rescaled to when it stands in for a position's
// own (absent) read-through population in the SC strand test. Large enough that the resulting
// hypergeometric is numerically a binomial against that ratio, small enough to stay well inside
// the range where fisher_exact_test_2x2's lgamma arithmetic is exact.
const double kSCStrandFallbackScale = 10000.0;

/*! Is a consensus clipped tail low-complexity enough to be an end-of-read artifact rather than
    donor sequence?

    Two ways to fail, because the same artifact shows up in both shapes. A dark-cycle poly-G tail
    is one unbroken run (GGGGGGGGGGGG, or CCCCCCCCCCCC once stored reference-forward for a leading
    clip); a tail that straddles the start of the dark cycles is broken but still nearly all one
    base (AAAAAACCCCCC, GGGGGGGGTTTT). Measured over 29 LTEE clones, 929 of 1040 accepted SC calls
    had a run of >= 8 in 12 bases, and no plausible real breakpoint did.

    Both thresholds are fractions of the compared length rather than absolute counts, so they hold
    when --soft-clipping-minimum-bases changes. Either fraction at 0 turns that half off.
 */
static bool sc_tail_is_low_complexity(const string& tail,
                                      double maximum_homopolymer_fraction,
                                      double maximum_base_fraction)
{
  if (tail.empty() || (tail == ".")) return false;

  // Non-ACGT columns are not evidence of anything either way, so they are excluded from both
  // the numerators and the denominator.
  map<char, uint32_t> counts;
  uint32_t informative = 0;
  uint32_t longest_run = 0, run = 0;
  char run_base = '\0';
  for (size_t i = 0; i < tail.size(); i++) {
    char b = static_cast<char>(toupper(static_cast<unsigned char>(tail[i])));
    if ((b != 'A') && (b != 'C') && (b != 'G') && (b != 'T')) { run = 0; run_base = '\0'; continue; }
    informative++;
    counts[b]++;
    if (b == run_base) run++; else { run_base = b; run = 1; }
    if (run > longest_run) longest_run = run;
  }
  if (informative == 0) return false;

  // The epsilon is not cosmetic: the intended reading of the 0.66 default at 12 compared bases is
  // "8 of 12", and a fraction that multiplies out to exactly the integer count must include it
  // rather than fall on the wrong side of a rounding step.
  const double kEpsilon = 1e-9;
  double n = static_cast<double>(informative);
  if ((maximum_homopolymer_fraction > 0.0) &&
      (static_cast<double>(longest_run) >= maximum_homopolymer_fraction * n - kEpsilon)) return true;

  uint32_t most_common = 0;
  for (map<char, uint32_t>::const_iterator it = counts.begin(); it != counts.end(); it++) {
    if (it->second > most_common) most_common = it->second;
  }
  if ((maximum_base_fraction > 0.0) &&
      (static_cast<double>(most_common) >= maximum_base_fraction * n - kEpsilon)) return true;

  return false;
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
    uint32_t agree_count_fw;   // agree_count split by the strand of the clipped read
    uint32_t agree_count_rv;
    uint32_t spanning_fw;      // read-through reads here, split the same way
    uint32_t spanning_rv;
    string   consensus_tail;
    double   score;
    double   frequency;
    double   consensus_fraction;
    double   fisher_strand_p_value;
    bool     suppressed;

    sc_candidate() : position(0), direction(0), clipped_count(0), total_count(0),
                     agree_count(0), agree_count_fw(0), agree_count_rv(0),
                     spanning_fw(0), spanning_rv(0), score(0.0), frequency(0.0),
                     consensus_fraction(0.0), fisher_strand_p_value(1.0),
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
    string expected = "#sc_format=3\tsoft_clipping_minimum_bases=" + to_string(_settings.soft_clipping_minimum_bases);
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
    ASSERT(fields.size() >= 11,
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
    uint32_t agree_count_fw = from_string<uint32_t>(fields[7]);
    uint32_t agree_count_rv = from_string<uint32_t>(fields[8]);
    uint32_t spanning_fw = from_string<uint32_t>(fields[9]);
    uint32_t spanning_rv = from_string<uint32_t>(fields[10]);

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
    c.agree_count_fw     = agree_count_fw;
    c.agree_count_rv     = agree_count_rv;
    c.spanning_fw        = spanning_fw;
    c.spanning_rv        = spanning_rv;
    c.consensus_tail     = consensus_tail;
    c.score              = score;
    c.frequency          = static_cast<double>(test_count) / static_cast<double>(total_count);
    c.consensus_fraction = static_cast<double>(agree_count) / static_cast<double>(clipped_count);

    /*
     * Strand test. Same 2x2 as the RA polymorphism one (see _add_polymorphism_bias_statistics):
     * the agreeing clipped reads are the "minor allele" and the reads that read through the
     * position are the "major allele", so a clip population drawn from a different strand mix
     * than the local coverage is what gets rejected.
     *
     * This is the discriminator that matters on real Illumina data. A dark-cycle poly-G tail --
     * the dominant SC false positive -- is always the 3' end of the read, so for a given clip
     * direction it can only come from one strand: direction +1 clips from forward reads,
     * direction -1 clips from reverse reads. Over 29 LTEE clones, 993 of 1040 accepted SC calls
     * had every clipped read on a single strand, while every real-looking breakpoint was
     * balanced. The consensus test cannot see this at all, because poly-G tails agree with each
     * other perfectly.
     *
     * Comparing against the LOCAL spanning strand split rather than against 50/50 is what keeps
     * a genuinely strand-skewed pileup from being read as a strand-skewed clip. Where there is
     * no local read-through at all -- exactly the frequency == 1.000 positions, which carry the
     * highest clip counts -- the contingency row is empty and Fisher returns 1 for anything;
     * falling back to the genome-wide spanning split restores the test there (it becomes, in
     * effect, a binomial test against the run's overall strand ratio).
     */
    uint32_t major_fw = spanning_fw;
    uint32_t major_rv = spanning_rv;
    if (major_fw + major_rv == 0) {
      // Rescaled to kSCStrandFallbackScale rather than passed at full size: the genome-wide
      // totals run into the hundreds of millions, where the lgamma differences inside
      // fisher_exact_test_2x2 lose precision and row1 + row2 can overflow its uint32_t. At this
      // size the hypergeometric is already indistinguishable from the binomial the fallback is
      // meant to be, so only the ratio matters.
      const double scale = kSCStrandFallbackScale;
      double gw_fw = static_cast<double>(summary.soft_clipping.total_spanning_read_bases_forward);
      double gw_rv = static_cast<double>(summary.soft_clipping.total_spanning_read_bases_reverse);
      double gw_total = gw_fw + gw_rv;
      if (gw_total > 0.0) {
        major_fw = static_cast<uint32_t>(floor(scale * gw_fw / gw_total + 0.5));
        major_rv = static_cast<uint32_t>(scale) - major_fw;
      }
    }
    // With no reference population at all -- neither local nor genome-wide -- there is nothing to
    // compare against, so leave the p-value at 1 rather than inventing a 50/50 expectation.
    c.fisher_strand_p_value = (major_fw + major_rv > 0)
      ? fisher_exact_test_2x2(c.agree_count_fw, c.agree_count_rv, major_fw, major_rv)
      : 1.0;

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
    // Exact (Clopper-Pearson) bounds on that same ratio -- k = agreeing clipped reads, n = all
    // reads reaching the position -- so the cutoff below is a confidence statement rather than a
    // point-estimate comparison, matching RA, JC and DP. A shallow position is then never rejected
    // merely for being shallow, only for being confidently below the cutoff.
    double sc_freq_lower = binomial_frequency_lower_bound(c.agree_count, c.total_count);
    double sc_freq_upper = binomial_frequency_upper_bound(c.agree_count, c.total_count);
    sc_entry[FREQUENCY_LOWER]  = formatted_double(sc_freq_lower, 4).to_string();
    sc_entry[FREQUENCY_UPPER]  = formatted_double(sc_freq_upper, 4).to_string();
    sc_entry[SC_CONSENSUS_FRACTION] = formatted_double(c.consensus_fraction, 4).to_string();
    if (c.consensus_tail != ".") sc_entry[SC_CONSENSUS_TAIL] = c.consensus_tail;
    sc_entry[SC_LOG10_E_VALUE] = formatted_double(c.score, kMutationScorePrecision).to_string();
    // Strand split of the reads the score was computed from, and of the read-through population
    // they were compared against. Reported whether or not the test rejected, because "all on one
    // strand" is the first thing to look at when judging an SC call by eye.
    sc_entry[SC_AGREE_COUNT_FORWARD] = to_string(c.agree_count_fw);
    sc_entry[SC_AGREE_COUNT_REVERSE] = to_string(c.agree_count_rv);
    sc_entry[SC_SPANNING_COUNT_FORWARD] = to_string(c.spanning_fw);
    sc_entry[SC_SPANNING_COUNT_REVERSE] = to_string(c.spanning_rv);
    // Spelled out rather than using output.h's FISHER_STRAND_P_VALUE, matching how the RA
    // polymorphism code in this file writes the same key.
    sc_entry["fisher_strand_p_value"] = formatted_double(c.fisher_strand_p_value, 5, true).to_string();

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
    //
    // Tested on the LOWER confidence bound, so this asks "is it confidently below the cutoff"
    // rather than "is the point estimate below the cutoff", as RA, JC and DP all do. No epsilon:
    // that guarded a rational agree/total landing exactly on a round cutoff, and an incomplete-beta
    // inverse never does. The bound sits below the point estimate, so this is strictly stricter --
    // see the cutoff note in Settings for what has and has not been recalibrated for it.
    if ((_settings.soft_clipping_frequency_cutoff > 0.0) &&
        (sc_freq_lower < _settings.soft_clipping_frequency_cutoff)) {
      sc_entry.add_reject_reason("FREQUENCY_BELOW_CUTOFF");
    }
    if ((_settings.soft_clipping_consensus_base_fraction > 0.0) &&
        (_settings.soft_clipping_consensus_fraction_cutoff > 0.0) &&
        (c.consensus_fraction < _settings.soft_clipping_consensus_fraction_cutoff - _settings.polymorphism_precision_decimal)) {
      sc_entry.add_reject_reason("CLIPPED_TAIL_CONSENSUS");
    }
    // Clipped reads drawn from a different strand mix than the reads that read through the
    // position; see the computation above for why this is the strongest SC filter there is.
    if ((_settings.soft_clipping_fisher_strand_p_value_cutoff > 0.0) &&
        (c.fisher_strand_p_value < _settings.soft_clipping_fisher_strand_p_value_cutoff)) {
      sc_entry.add_reject_reason("FISHER_STRAND");
    }
    // The clipped tail is a homopolymer or near-homopolymer: a dark-cycle or adapter artifact,
    // not donor sequence. Independent of the strand test, and it reaches the low-count positions
    // where the strand test has no power.
    if (sc_tail_is_low_complexity(c.consensus_tail,
                                  _settings.soft_clipping_maximum_tail_homopolymer_fraction,
                                  _settings.soft_clipping_maximum_tail_base_fraction)) {
      sc_entry.add_reject_reason("LOW_COMPLEXITY_TAIL");
    }

    _gd.add(sc_entry);
  }
}

/*! Destructor.
 */
identify_mutations_pileup::~identify_mutations_pileup()
{
}

/*! The DP/MP/PD group a read belongs to, or -1 if its read group has none.
 */
int identify_mutations_pileup::dp_group_for(const pileup& p, const alignment_wrapper& a) const
{
  uint32_t rg = p.read_group_index(a);
  if (rg >= _read_group_to_dp_group.size()) return -1;
  return _read_group_to_dp_group[rg];
}

/*! One entry of a candidate region's discordant_pairs list: the read-pair key followed by the
    ALIGNED extent of this region's read of that pair, "<key>__<read_start>__<read_end>".

    dp_evidence seeds each candidate junction side's boundary search from these extents (the reads
    that actually support the edge), rather than from the region's span -- which is a product of the
    sliding window (its end is a read's start plus the window width) and is shared by every partner
    of a hub region. The key is itself '__'-separated (<read1>__<read2>__<insert_size>), so the
    identity key stays recoverable as the first three fields.
 */
static string dp_descriptor(const identify_mutations_pileup::dp_read& r)
{
  return r.key + "__" + to_string(r.start_pos) + "__" + to_string(r.end_pos);
}

/*! One entry of an MP candidate region's unpaired_reads list: the ALIGNED extent of one
    mate-unmapped read, "<read_start>__<read_end>".

    mp_evidence seeds each region's boundary search from these extents, for the same reason DP does
    (the region span's far bound is an artifact of the sliding window). There is no identity key
    here: with the mate unmapped there is no second region to join this read to, so nothing needs to
    be matched up across regions -- the extents alone are what the refinement consumes.
 */
static string mp_descriptor(const identify_mutations_pileup::mp_read& r)
{
  return to_string(r.start_pos) + "__" + to_string(r.end_pos);
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
  if (_this_deletion_propagation_cutoff < 0.0) {

    // ...which is exactly what a reference sequence with no reads mapped to it looks like, and that
    // reference still needs its coverage table. CNery is handed the table paths and raises on a
    // missing one, and its own "this reference had nothing to measure" reporting needs a file it can
    // open; a header row with nothing under it used to be all it got. Tallied separately here
    // because the ordinary tally is built further down, past the DP/MP/PD accumulators that this
    // return is what keeps out of a reference breseq is not calling mutations on.
    if (_print_coverage_data && _coverage_data.is_open()) {
      position_coverage pc;
      vector<position_coverage> pc_by_rg;
      if (_print_read_group_coverage) pc_by_rg.resize(read_groups().size());
      tabulate_position_coverage(p, pc, pc_by_rg);
      write_coverage_row(p.position_1(), p.reference_base_char_1(p.position_1()), pc, pc_by_rg);
    }

    return;
  }

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

    //## the same, split by read group, for the coverage table CNery reads. Kept separate from the
    //## aggregate above rather than summed out of it, so the numbers every other consumer sees are
    //## unchanged. Empty (and never touched) unless there is more than one read group.
    vector<position_coverage> this_position_coverage_by_rg;
    if (_print_read_group_coverage) this_position_coverage_by_rg.resize(read_groups().size());

		//## polymorphism prediction data
		vector<polymorphism_data> pdata;

    //## Summed log10 P(observed bases | this position is 100% base b), one entry per candidate
    //## base. Every pure-genotype quantity at this position is a function of just these five
    //## numbers; see pure_genotype_call().
    double log10_pr_sum[base_list_size];
    for (uint8_t j=0; j<base_list_size; j++) log10_pr_sum[j] = 0.0;
    
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
          int gi = dp_group_for(p, *i);
          if (gi >= 0) {
            const dp_group& g = _dp_groups[gi];

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
            // breseq names mates <file_index>:<read_num> sharing read_num (fastq.cpp); both mates
            // are in one read group, so R1/R2 role comes from the mate flag.
            bool focal_is_r1 = !(i->flag() & BAM_FREAD2);
            string fname = i->read_name();
            size_t colon = fname.find(':');
            string read1_name, read2_name;
            if (colon != string::npos) {
              string counter = fname.substr(colon + 1);
              read1_name = to_string(g.r1_read_name_prefix) + ":" + counter;
              read2_name = to_string(g.r2_read_name_prefix) + ":" + counter;
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
            dr.end_pos = i->reference_end_1();
            dr.key = key;
            // A tie-broken redundant discordant pair keeps X1>1 on its marked alignment (see
            // resolve_alignments.cpp); flag its DP side redundant so it propagates to DP evidence.
            dr.redundant = (i->redundancy() > 1);
            _dp_groups[gi].reads[bin].push_back(dr);
            ++_dp_metric[bin];
            // If a region in this bin is already open (from a previous column), this newly entering
            // read belongs to it, so append its key now. Reads entering on the column that OPENS a
            // region are instead captured by the open-time snapshot in check_discordant_completion
            // (which runs after this loop), so they are never double-counted.
            if (_dp_region_start[bin] != UNDEFINED_UINT32) {
              _dp_region_descriptors[bin].push_back(dp_descriptor(dr));
              if (dr.redundant) ++_dp_region_redundant_count[bin];
            }
          }
        }
      }

      //## Missing-pair (MP) region detection: incremental "enter" step.
      //## Same once-per-read-at-its-start-column bookkeeping as DP above, but the qualifying read is
      //## a SINGLETON: mapped and flagged paired with its mate unmapped anywhere (see
      //## mark_mate_unmapped in resolve_alignments.cpp). Such a read has no XP tag and no meaningful
      //## insert size, so the only bin dimension is the focal read's strand. A pile of these facing
      //## one point is the signature of a novel sequence inserted there: the mates landed in
      //## sequence that is in neither the reference nor any candidate junction.
      if (_mp_enabled && (insert_count == 0) && (i->reference_start_1() == position)) {
        if (!i->unmapped() && !(i->flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) {
          int gi = dp_group_for(p, *i);
          if (gi >= 0) {
            int bin = i->reversed() ? 1 : 0;

            // Every primary mapped read of this strand enters the denominator window, whether or not
            // its mate is missing. The seed test below is an ENRICHMENT test, and it needs both terms.
            _dp_groups[gi].mp_all_reads[bin].push_back(position);
            ++_mp_total_metric[bin];

            // BAM_FMUNMAP alone is not enough: it is also set for a mate the aligner placed but
            // breseq rejected on --require-match-fraction / mapping quality / mismatches. Such a
            // mate is a partially-aligning read (SC's business), not sequence missing from the
            // reference. Requiring the narrower tag is what keeps this count from tracking
            // alignment stringency -- see mate_never_aligned() and kBreseqMateNeverAlignedBAMTag.
            if (i->is_paired() && (i->flag() & BAM_FMUNMAP) && mate_never_aligned(*i)) {
              mp_read mr;
              mr.start_pos = position;
              mr.end_pos = i->reference_end_1();
              mr.redundant = (i->redundancy() > 1);
              _dp_groups[gi].mp_reads[bin].push_back(mr);
              ++_mp_metric[bin];
              // As for DP: reads entering on the column that OPENS a region are captured by the
              // open-time snapshot in check_missing_pair_completion (which runs after this loop),
              // so appending here only for an already-open region cannot double-count.
              if (_mp_region_start[bin] != UNDEFINED_UINT32) {
                _mp_region_descriptors[bin].push_back(mp_descriptor(mr));
                if (mr.redundant) ++_mp_region_redundant_count[bin];
              }
            }
          }
        }
      }

      //## Pair-distance (PD) region detection: incremental "enter" step.
      //## Unlike DP and MP, the qualifying read here is an ORDINARY one: a normally-oriented pair
      //## whose two mates both mapped. Nothing about it is individually anomalous, and that is the
      //## point -- an event of a few hundred bases shifts every spanning pair by the same amount
      //## without pushing any single one past the discordant cutoff.
      //##
      //## A pair constrains a breakpoint only if the breakpoint lies in its unsequenced middle gap:
      //## anywhere else and one of the mates would have crossed it and been clipped into a junction
      //## read instead. In between-base coordinates (b = "between reference base b and b+1") that
      //## gap is [E1, S2-1] for a leftmost mate ending at E1 whose mate starts at S2. Both bounds
      //## come from this alignment alone, so no mate lookup is needed. Note the gap is EMPTY when
      //## the mates overlap (S2 <= E1) -- such a pair carries no positional information at all, and
      //## this is why PD's reach for insertions stops around (median distance - 2 * read length).
      //##
      //## The requirement S2 > E1 also self-selects the leftmost mate of the pair (the rightmost
      //## mate always has mate_start <= its own end), so each pair is counted exactly once without
      //## an explicit leftmost test -- including the equal-starts tie where mark_pair_info marks
      //## both mates leftmost.
      if (_pd_enabled && (insert_count == 0) && (i->reference_start_1() == position)) {
        // Qualification is deliberately all integer flag tests, with no aux-tag lookup: this block
        // sees ~100% of reads rather than DP's ~1%, so the string built by aux_get_Z("XP") would be
        // a per-read cost. Orientation comes from the flags instead -- for the pair whose leftmost
        // mate we are looking at, FR means this read is forward and its mate reverse, RF the other
        // way around -- and is then required to MATCH the group's majority orientation, the one the
        // null was built from. Only FF/RR are excluded here; the FR/RF discrimination needs the
        // group, so it happens below.
        if (i->is_paired() && !i->unmapped()
            && !(i->flag() & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FMUNMAP))
            && (i->mate_reference_target_id() == i->reference_target_id())
            && (i->insert_size() > 0)
            && (i->reversed() != ((i->flag() & BAM_FMREVERSE) != 0))
            // A tie-broken multi-mapping placement has an arbitrary distance, so it would add noise
            // to the distribution rather than signal. This is the opposite of DP, which USES
            // redundancy as evidence that a side sits on a multicopy element.
            && (i->redundancy() <= 1)) {

          // ...and the same for the MATE, which redundancy() cannot see. Only the leftmost mate of a
          // pair reaches here (insert_size() > 0), so the test above covers whichever mate that
          // happens to be and leaves the other unexamined -- half the affected pairs, chosen by a
          // coin flip. Worse, X1 is written from the alignments that SURVIVED downselection, so a
          // mate whose copy was picked to make the pair look concordant reads back as unique.
          // kBreseqAmbiguousPairPlacementBAMTag records the choice at the point it was made.
          //
          // Tested here rather than in the chain above so the two outcomes can be tallied: the claim
          // that refusing these pairs is free has to be checkable on every run.
          _pd_pairs_considered++;
          bool chosen_placement = ambiguous_pair_placement(*i);
          if (chosen_placement) _pd_pairs_ambiguous_placement++;

          int gi = chosen_placement ? -1 : dp_group_for(p, *i);
          if (gi >= 0 && pd_orientation_matches(_dp_groups[gi].pd_orientation, *i)) {
            const vector<int64_t>& cdf = _dp_groups[gi].pd_cdf;
            int32_t d = i->insert_size();

            int32_t lo = static_cast<int32_t>(i->reference_end_1());
            int32_t hi = i->mate_start_1() - 1;

            // Empty gap (overlapping mates), no null for this group, or a pair longer than the ring
            // can hold. The span cap is not an optimization: a pair writes into ring slots up to
            // position+d, so one longer than the ring would alias a slot the drain has not reached.
            if ((hi >= lo) && !cdf.empty() && (d <= _pd_max_span)) {
              int64_t u = (d < static_cast<int32_t>(cdf.size())) ? cdf[d] : kPDuScale;

              int bin = -1;
              int64_t lower_tail = static_cast<int64_t>(_settings.pair_distance_tail_quantile * static_cast<double>(kPDuScale));
              if (u >= kPDuScale - lower_tail)   bin = kPDbinLong;
              else if (u <= lower_tail)          bin = kPDbinShort;

              pd_add_interval(lo, hi, bin, u);
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
      int strand = i->strand();
      bool trimmed = i->is_trimmed(on_insert_position_past_base);
      
      //##### update coverage if this is not a deletion in read relative to reference
      //### note that we count trimmed reads here, but not when looking for short indel mutations...	
      
      position_coverage* this_position_coverage_rg = NULL;
      if (!this_position_coverage_by_rg.empty())
        this_position_coverage_rg = &this_position_coverage_by_rg[p.read_group_index(*i)];

      if(redundancy == 1) {
        //## keep track of unique coverage
        ++this_position_coverage.unique[1+strand];
        if (this_position_coverage_rg) ++this_position_coverage_rg->unique[1+strand];

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
        if (this_position_coverage_rg) {
          this_position_coverage_rg->redundant[1+strand] += 1.0/redundancy;
          ++this_position_coverage_rg->raw_redundant[1+strand];
        }

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
        pdata.push_back(polymorphism_data(baseindex2char(cv.obs_base()),cv.quality(),i->strand(), i->mapping_quality(), cv));

        //cerr << " " << cv.obs_base() << " " << (char)ref_base << endl;

        // One error-table query per hypothesis per read base, done once here. Everything
        // downstream reads the cache.
        fill_read_base_likelihoods(pdata.back());
        for (uint8_t j=0; j<base_list_size; j++) log10_pr_sum[j] += pdata.back()._log10_pr[j];
      }
		} // end for-each read
		
    
    //#############################
		//## PER POSITION/INSERT COUNT
    //#############################
    
		//#sum up coverage observations
		this_position_coverage.sum();
    for (size_t g = 0; g < this_position_coverage_by_rg.size(); g++)
      this_position_coverage_by_rg[g].sum();

		//#we are trying to find the base with the most support
    cSNPCall snp_call = pure_genotype_call(log10_pr_sum, static_cast<uint32_t>(pdata.size()));
    
    base_char best_base_char('N');
    double consensus_bonferroni_score(numeric_limits<double>::quiet_NaN());
    double variant_score(numeric_limits<double>::quiet_NaN());

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
		
    if (insert_count == 0) {

      // The coverage table is written whatever happens to missing-coverage prediction. It used to be
      // written from inside check_deletion_completion(), so --skip-MC-prediction and
      // --targeted-sequencing (which forces it on) left CNery a header row and no data for every
      // reference sequence, however well covered.
      write_coverage_row(position, ref_base_char, this_position_coverage, this_position_coverage_by_rg);

      if(!_settings.skip_missing_coverage_prediction)
        check_deletion_completion(p.target(), position, this_position_coverage, consensus_bonferroni_score);
    }

		//###
		//## DISCORDANT PAIR DISCORDANT PAIR DISCORDANT PAIR
		//###

    // Runs once per reference column (including empty columns) so the sliding window advances and
    // regions close correctly even where there is no coverage.
    if(!_dp_groups.empty() && (insert_count == 0))
      check_discordant_completion(p.target(), position);

		//###
		//## MISSING PAIR MISSING PAIR MISSING PAIR
		//###

    // Same per-column cadence as DP, for the same reason.
    if(_mp_enabled && (insert_count == 0))
      check_missing_pair_completion(p.target(), position);

		//###
		//## PAIR DISTANCE PAIR DISTANCE PAIR DISTANCE
		//###

    // Same per-column cadence again. For PD this is not merely so regions close correctly: the
    // running totals are a prefix sum over a difference array, so every column must be drained.
    if(_pd_enabled && (insert_count == 0))
      check_pair_distance_completion(p.target(), position);

		//###
		//## POLYMORPHISM POLYMORPHISM POLYMORPHISM
		//###								

		bool passed_as_polymorphism_prediction(false);

    cDiffEntry mut(RA);

    //## evaluate whether to call an actual mutation!
    // -- note that we accept > 0 and only reject later
    // so that these can potentially make it into the marginal data
    bool passed_as_consensus_prediction = (best_base_char != ref_base_char) && (!std::isnan(consensus_bonferroni_score) && (consensus_bonferroni_score > 0));

    //## Fit all five alleles at once, so every frequency below is a fraction of TOTAL depth and
    //## the reference allele is in the model whether or not it is well supported. The two-base
    //## fit this replaces chose its pair by raw coverage and renormalized within it, which made
    //## the reported frequency a share of that pair rather than of the position.
    bool all_alleles[base_list_size];
    for (uint8_t i=0; i<base_list_size; i++) all_alleles[i] = true;
    allele_model amodel = fit_allele_frequencies(pdata, all_alleles);

    //## Three questions, three bases, none of them reassigned into another:
    //##   best_base_char (above)  the most probable single genotype -- consensus score, UN calling
    //##   major_base_char         the highest-frequency allele in the fit
    //##   variant_base_char       the highest-frequency allele that is not the reference
    //## They agree except at genuinely mixed sites, where the disagreement is information. The
    //## old code overwrote the first with the second partway through and reported a mixture of
    //## the two.
    const base_index ref_base_index = (ref_base_char == 'N') ? base_list_N_index : basechar2index(ref_base_char);
    const base_index major_base_index = amodel.major_index();
    const base_index minor_base_index = amodel.next_index(major_base_index);
    const base_index variant_base_index = amodel.next_index(ref_base_index);

    base_char major_base_char   = (major_base_index   == base_list_N_index) ? 'N' : base_char_list[major_base_index];
    base_char minor_base_char   = (minor_base_index   == base_list_N_index) ? 'N' : base_char_list[minor_base_index];
    base_char variant_base_char = (variant_base_index == base_list_N_index) ? 'N' : base_char_list[variant_base_index];

    if (variant_base_index != base_list_N_index) {
      variant_score = variant_presence_score(pdata, amodel, variant_base_index);

      // Do we accept this as a polymorphism?
      if (variant_score >= _polymorphism_score_cutoff)
        passed_as_polymorphism_prediction = true;
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
      //## One score: the log10 evidence that this position carries a non-reference allele at all.
      //## Whether that allele is then called fixed or polymorphic is decided by the frequency
      //## cutoffs, not by a second score. The pure-genotype log-odds is still computed -- it gates
      //## UN calling and the consensus emission test above -- but it answers a different question
      //## ("which single base is this") and reporting both invited them to be compared.
      mut[SCORE] = formatted_double(variant_score, kMutationScorePrecision).to_string();
      
      //## Specific initializations for polymorphisms. Must take precedence.

      mut[REF_BASE] = ref_base_char;
      mut[NEW_BASE] = variant_base_char;

      mut[MAJOR_BASE] = major_base_char;
      mut[MINOR_BASE] = minor_base_char;

      //## Both of these are fractions of TOTAL depth now, so unlike the old pair-relative
      //## frequencies they do NOT sum to 1 -- the whole fitted spectrum does.
      mut[MAJOR_FREQUENCY] = formatted_double(
          amodel.reported_frequency(major_base_index), _polymorphism_precision_places, true).to_string();
      mut[FREQUENCY] = formatted_double(
          amodel.reported_frequency(variant_base_index), _polymorphism_precision_places, true).to_string();

      //## The whole fitted spectrum, so a position with more than one credible non-reference
      //## allele is inspectable rather than silently reduced to its strongest one.
      mut[ALLELE_FREQUENCIES] = amodel.spectrum_string(_polymorphism_precision_places);

      write_RA_frequency_bounds(mut, pdata, amodel, variant_base_index);

      //## Strand and quality bias statistics, computed for every entry.
      //##
      //## These used to be skipped whenever the position had only one allele with coverage, which
      //## sounds harmless -- with nothing to compare against, the tests are uninformative -- but it
      //## silently exempted those entries from a filter that is ON by default:
      //## polymorphism_fisher_strand_p_value_cutoff is 0.05 in BOTH modes, and
      //## rejected_RA_polymorphism_bias() only applies it when fisher_strand_p_value exists. So the
      //## absence of the field was doing the work of a passing test. Compute it always and let the
      //## filter decide; annotate_polymorphism_statistics() already returns 1.0 for an empty
      //## comparison, which is a pass on the merits rather than by omission.
      annotate_polymorphism_statistics(mut, major_base_char, minor_base_char, pos_info, pdata);

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

        //## Read the user's allele straight off the position's own fit, rather than re-fitting a
        //## two-base model of (ref_base, new_base). The old private fit is what forced the major
        //## and minor alleles here to be defined pair-relatively -- with the acknowledged caveat
        //## that "there could be a 3rd base that matters" -- and it took the reference base as an
        //## error-table lookup key, which asserts outright when the reference is 'N'. Neither
        //## problem survives sharing the five-allele fit: a third allele is already in the model,
        //## and the reference base is never used as a key.
        const base_char user_ref_base_char = from_string<base_char>(mut[REF_BASE]);
        const base_char user_new_base_char = from_string<base_char>(mut[NEW_BASE]);
        const base_index user_ref_index = basechar2index(user_ref_base_char);
        const base_index user_variant_index = basechar2index(user_new_base_char);

        double user_polymorphism_score = numeric_limits<double>::quiet_NaN();
        double user_variant_frequency = 0.0;
        if (user_variant_index < base_list_size) {
          user_polymorphism_score = variant_presence_score(pdata, amodel, user_variant_index);
          user_variant_frequency = amodel.reported_frequency(user_variant_index);
        }
        // A user entry may name a reference base of 'N', which has no allele in the model.
        double user_ref_frequency = amodel.reported_frequency(user_ref_index);

        //## Major and minor stay defined WITHIN the pair the user named, even though the
        //## frequencies now come from the position's full fit. Reporting the position's own major
        //## and minor here instead would silently discard the user's question:
        //## mutation_predictor recovers the allele as (major == ref) ? minor : major, so at a
        //## position carrying no real variant every user entry would come back as 'N'.
        bool user_variant_is_major = (user_variant_frequency > user_ref_frequency);
        mut[MAJOR_BASE] = user_variant_is_major ? user_new_base_char : user_ref_base_char;
        mut[MINOR_BASE] = user_variant_is_major ? user_ref_base_char : user_new_base_char;
        mut[MAJOR_FREQUENCY] = formatted_double(
            user_variant_is_major ? user_variant_frequency : user_ref_frequency,
            _polymorphism_precision_places, true).to_string();

        mut[FREQUENCY] = formatted_double(user_variant_frequency, _polymorphism_precision_places, true).to_string();
        mut[ALLELE_FREQUENCIES] = amodel.spectrum_string(_polymorphism_precision_places);
        write_RA_frequency_bounds(mut, pdata, amodel, user_variant_index);

        // [frequency] keeps the fitted estimate in both modes; the call itself goes in [prediction],
        // the same place test_RA_evidence() records it for non-user entries. Consensus mode asks
        // only whether the user's variant is the majority allele. A "polymorphism" verdict here is
        // what drops a sub-majority user entry in consensus mode -- mutation_predictor keeps only
        // consensus evidence when not predicting polymorphisms.
        if (!_settings.polymorphism_prediction) {
          mut[PREDICTION] = (user_variant_frequency > 0.5) ? "consensus" : "polymorphism";
        } else {
          mut[PREDICTION] = "polymorphism";
        }

        mut[SCORE] = formatted_double(user_polymorphism_score, kMutationScorePrecision).to_string();

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

    // The same schema 'bam2cov --total-only' writes, so CNery reads this table with exactly the
    // code path it uses for one it generated itself. Strand-split counts would only be summed
    // again on the way in, so they are not written; nothing else reads this file.
    _coverage_data << "position" << "\t";
    _coverage_data << "ref_base";
    for (size_t c = 0; c < coverage_column_names().size(); c++)
      _coverage_data << "\t" << coverage_column_names()[c];

    // Group-major, after the aggregate columns, as in coverage_output::table().
    if (_print_read_group_coverage) {
      for (size_t g = 0; g < read_groups().size(); g++)
        for (size_t c = 0; c < coverage_column_names().size(); c++)
          _coverage_data << "\t" << coverage_output::read_group_column_prefix(g) << coverage_column_names()[c];
    }
    _coverage_data << endl;
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

  // Reset the Missing Pair (MP) evidence variables (regions do not span reference boundaries)
  for (int bin = 0; bin < kMPnBins; bin++) {
    for (vector<dp_group>::iterator g = _dp_groups.begin(); g != _dp_groups.end(); g++) {
      g->mp_reads[bin].clear();
      g->mp_all_reads[bin].clear();
    }
    _mp_metric[bin] = 0;
    _mp_total_metric[bin] = 0;
    _mp_region_start[bin] = UNDEFINED_UINT32;
    _mp_region_max_count[bin] = 0;
    _mp_region_max_total[bin] = 0;
    _mp_region_descriptors[bin].clear();
    _mp_region_redundant_count[bin] = 0;
  }

  // Reset the Pair Distance (PD) evidence variables (regions do not span reference boundaries).
  // The rings must be zeroed wholesale, not just drained: pairs near the end of the previous
  // reference wrote into slots ahead of where its last column drained, and those pending deltas
  // are meaningless here. The running totals go with them.
  if (_pd_enabled) {
    _pd_ring_n.assign(_pd_ring_w, 0);
    _pd_ring_long.assign(_pd_ring_w, 0);
    _pd_ring_short.assign(_pd_ring_w, 0);
    _pd_ring_u.assign(_pd_ring_w, 0);
    _pd_n = 0;
    _pd_long = 0;
    _pd_short = 0;
    _pd_u = 0;
    _pd_last_b = 0;
    for (int bin = 0; bin < kPDnBins; bin++) {
      _pd_region_start[bin] = UNDEFINED_UINT32;
      _pd_region_peak_position[bin] = 0;
      _pd_region_peak_covering[bin] = 0;
      _pd_region_peak_tail[bin] = 0;
      _pd_region_peak_z[bin] = 0.0;
      // Only the open flags: the ladder COUNTS are a whole-reference-set tally, which is the scale
      // the multiple-testing correction is defined at. The flush call at the end of the previous
      // target already closed anything open, so this is belt and braces.
      for (int step = 0; step < kPDLadderSteps; step++) _pd_ladder_open[bin][step] = false;
    }
  }

}

/*! Called at the end of a reference sequence fragment
    Close per-reference files
 */
void identify_mutations_pileup::at_target_end(const uint32_t tid) {

  // end "open" Missing Coverahge and Unknown intervals
  if (!_settings.skip_missing_coverage_prediction) {
    check_deletion_completion(tid, target_length(tid)+1, position_coverage(numeric_limits<double>::quiet_NaN()), numeric_limits<double>::quiet_NaN());
  }
  update_unknown_intervals(target_length(tid)+1, tid, true, false);

  // Flush any open discordant-pair (DP) region at the end of this reference sequence.
  if (!_dp_groups.empty()) {
    check_discordant_completion(tid, target_length(tid)+1);
  }

  // Flush any open missing-pair (MP) region at the end of this reference sequence.
  if (_mp_enabled) {
    check_missing_pair_completion(tid, target_length(tid)+1);
  }

  // Flush any open pair-distance (PD) region at the end of this reference sequence.
  if (_pd_enabled) {
    check_pair_distance_completion(tid, target_length(tid)+1);
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


/*! Write the unique/redundant/total triple of one coverage tally, in coverage_column_names() order.

    'total_cov' is the plain sum rather than position_coverage::total, which sum() rounds each strand
    tally before adding. bam2cov's column of the same name is the unrounded sum, and this table has
    to mean what that one means.
 */
void identify_mutations_pileup::write_coverage_columns(const position_coverage& pc)
{
  _coverage_data << "\t" << pc.unique[1];
  _coverage_data << "\t" << pc.redundant[1];
  _coverage_data << "\t" << (pc.unique[1] + pc.redundant[1]);
}


/*! Write one position's row of the coverage table CNery reads.

    Moved out of check_deletion_completion(), where it used to live, because that call is gated on
    --skip-MC-prediction being off AND on this reference's coverage-distribution fit having
    succeeded. Neither has anything to do with the table, and both silently emptied it: a reference
    with no coverage gets propagation cutoff -1, every position returns from pileup_callback() before
    reaching here, and CNery was handed a file holding nothing but its header row.

    The isnan() guard is what used to keep the once-per-fragment end call (position = length+1, with
    an undefined position_coverage) from writing a row past the end of the sequence. Nothing reaches
    here from there any more, but the guard is kept: it costs one comparison, and it is the only
    thing standing between an uninitialized tally and a row of NaNs in CNery's input.
 */
void identify_mutations_pileup::write_coverage_row(uint32_t position, char ref_base_char, const position_coverage& pc, const vector<position_coverage>& pc_by_rg)
{
  if (std::isnan(pc.unique[1]) || !_coverage_data.is_open()) return;

  _coverage_data << position << "\t";
  _coverage_data << ref_base_char;
  write_coverage_columns(pc);
  for (size_t g = 0; g < pc_by_rg.size(); g++)
    write_coverage_columns(pc_by_rg[g]);
  _coverage_data << endl;
}


/*! Tally coverage at one position for the coverage table alone.

    Only reached for a reference whose coverage-distribution fit failed, where pileup_callback()
    returns before building its own tally. Counting here rather than reusing that tally costs an
    extra pass over the pileup, so it is deliberately confined to those references: on every other
    one this is never called and the ordinary tally is used.

    Mirrors what the main loop counts at insert_count 0 -- every aligned read that is not a deletion
    at this position and whose base carries a quality, split by strand, redundant reads weighted
    1/redundancy -- so a low-coverage reference reports its real depth rather than zeros.
 */
void identify_mutations_pileup::tabulate_position_coverage(const pileup& p, position_coverage& pc, vector<position_coverage>& pc_by_rg)
{
  for (pileup::const_iterator i = p.begin(); i != p.end(); ++i) {

    // insert_count is 0 here, so the only reads without a base at this position are deletions.
    if (i->is_del()) continue;

    // Don't count bases without qualities -- not safe to even count them for coverage.
    base_bam read_base_bam = i->read_base_bam_0(i->query_position_0());
    if (_base_bam_is_N(read_base_bam)) continue;

    int32_t redundancy = i->redundancy();
    int strand = i->strand();

    position_coverage* pc_rg = NULL;
    if (!pc_by_rg.empty()) pc_rg = &pc_by_rg[p.read_group_index(*i)];

    if (redundancy == 1) {
      ++pc.unique[1+strand];
      if (pc_rg) ++pc_rg->unique[1+strand];
    } else {
      pc.redundant[1+strand] += 1.0/redundancy;
      ++pc.raw_redundant[1+strand];
      if (pc_rg) {
        pc_rg->redundant[1+strand] += 1.0/redundancy;
        ++pc_rg->raw_redundant[1+strand];
      }
    }
  }

  pc.sum();
  for (size_t g = 0; g < pc_by_rg.size(); g++) pc_by_rg[g].sum();
}


/*! Helper method to track information about putative deleted regions.
 
 Used at each pileup iteration and at the end.
 //## when called at the end of a fragment, the position is fragment_length+1
 //## and this_position_coverage is undefined
 
 @JEB This function expects 1-indexed positions!!!
 
 */
void identify_mutations_pileup::check_deletion_completion(uint32_t seq_id, uint32_t position, const position_coverage& this_position_coverage, double e_value_call) {

	//cerr << position << " " << e_value_call << endl;
	
  // special case = beginning of new seq_id
  if (position == 1) 
    _last_position_coverage = position_coverage(numeric_limits<double>::quiet_NaN());
  
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
            _dp_region_descriptors[bin].push_back(dp_descriptor(*r));
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


/*! Helper method to track missing-pair (MP) candidate regions.

 Called once per reference column (including empty columns), and once past the end of each reference
 sequence (position = target_length+1) to flush an open region. Maintains, per paired read group, a
 sliding window of width (median + 2.42*MAD) over the reference start positions of reads whose mate
 did not map anywhere. Reads are added once at their start column ("enter", in pileup_callback) and
 removed once here when they fall out of the window ("exit"). A candidate region is a maximal run of
 columns where the total in-window count reaches --missing-pair-seed.

 The only difference from check_discordant_completion is the bin space: MP has no pair-orientation
 dimension, so there are two bins (focal strand) rather than six.

 @JEB expects 1-indexed positions.
 */
void identify_mutations_pileup::check_missing_pair_completion(uint32_t seq_id, uint32_t position)
{
  if (!_mp_enabled) return;

  // Run an independent detector for each focal strand, so a breakpoint's forward shoulder (mates
  // reaching right into the insert) and its reverse shoulder become separate regions.
  for (int bin = 0; bin < kMPnBins; bin++) {

    // "Exit" step: drop reads that have fallen out of each group's bin window -- both the
    // mate-unmapped reads and the same-strand denominator.
    for (vector<dp_group>::iterator g = _dp_groups.begin(); g != _dp_groups.end(); g++) {
      while (!g->mp_reads[bin].empty() &&
             (static_cast<double>(position) - static_cast<double>(g->mp_reads[bin].front().start_pos) >= g->mp_window_width)) {
        g->mp_reads[bin].pop_front();
        --_mp_metric[bin];
      }
      while (!g->mp_all_reads[bin].empty() &&
             (static_cast<double>(position) - static_cast<double>(g->mp_all_reads[bin].front()) >= g->mp_window_width)) {
        g->mp_all_reads[bin].pop_front();
        --_mp_total_metric[bin];
      }
    }

    // Accumulate the null BEFORE the seed test, and over every column -- seeded or not. The quiet
    // columns are the whole sample: they are what says how often a read's mate fails to align here
    // and how unevenly that is spread, and they never appear in the candidate-region CSV.
    //
    // Each read is counted once per column it is in-window for, about W times over. That
    // multiplicity is uniform, so it cancels in p0 = S_k / S_n. It does NOT cancel in the column
    // count, so _mp_columns is a redundant sample size, not an independent one -- see fill_mp_summary.
    if (position <= target_length(seq_id)) {
      int32_t k = _mp_metric[bin];
      int32_t n = _mp_total_metric[bin];
      if (n > 0) {
        _mp_columns++;
        _mp_S_n += static_cast<double>(n);
        _mp_total_k += static_cast<uint64_t>(k);
        double f = static_cast<double>(k) / static_cast<double>(n);
        // A real insertion must not define the background it is judged against. Measured on real
        // data this matters a lot -- the same point soft_clipping.cpp makes about its own trim.
        if ((_settings.missing_pair_dispersion_trim_frequency > 0.0) &&
            (f >= _settings.missing_pair_dispersion_trim_frequency)) {
          _mp_columns_trimmed++;
          _mp_S_n_trimmed += static_cast<double>(n);
        } else {
          _mp_S_k += static_cast<double>(k);
          _mp_S_k2_over_n += static_cast<double>(k) * static_cast<double>(k) / static_cast<double>(n);
        }
      }
    }

    // A region opens only where mate-unmapped reads are ENRICHED, not merely present: an absolute
    // count alone clears any fixed threshold at every position once coverage is decent, and the
    // entire covered region becomes one candidate.
    //
    // This test is a SENSITIVITY FILTER and nothing more. What decides an MP call is the genome-wide
    // score computed in mp_evidence against the null accumulated just above.
    //
    // This comment used to claim the mate-unmapped fraction is "roughly constant along a genome --
    // a property of the library and the aligner, not of the locus", and the whole model rested on
    // it. It is not true, and it is not even the same number from run to run. Two 2x150 E. coli
    // libraries off one flowcell measure 0.72% and 1.00% by the definition used here (a mate the
    // aligner placed nowhere), with Pearson dispersion phi of 1.8 and 1.6. Under the older, wider
    // definition -- any mate without an ELIGIBLE alignment -- the same libraries read 2.2%, phi 4.8,
    // with 3.9% of all windows above a local fraction of 0.10.
    //
    // So a fixed fraction like the 0.25 here has no fixed relationship to the background, and the
    // fixed ACCEPTANCE cutoff of 0.10 this model used to rely on sat below the 96th percentile of
    // pure noise. Keep the seed loose and cheap; let the score, which knows what the background
    // actually is, do the deciding.
    //
    // Never "above" once past the end of the sequence: closes an open region at the flush sentinel
    // and prevents opening a spurious region that would never close.
    bool enriched = (_mp_total_metric[bin] > 0) &&
                    (static_cast<double>(_mp_metric[bin]) >= _mp_seed_fraction * static_cast<double>(_mp_total_metric[bin]));
    bool above = (_mp_metric[bin] >= _mp_seed) && enriched && (position <= target_length(seq_id));

    if (above) {
      if (_mp_region_start[bin] == UNDEFINED_UINT32) {
        // Open a new region; snapshot the bin's reads currently in-window as its initial extents.
        _mp_region_start[bin] = position;
        _mp_region_max_count[bin] = static_cast<uint32_t>(_mp_metric[bin]);
        _mp_region_max_total[bin] = static_cast<uint32_t>(_mp_total_metric[bin]);
        _mp_region_descriptors[bin].clear();
        _mp_region_redundant_count[bin] = 0;
        for (vector<dp_group>::iterator g = _dp_groups.begin(); g != _dp_groups.end(); g++) {
          for (deque<mp_read>::iterator r = g->mp_reads[bin].begin(); r != g->mp_reads[bin].end(); r++) {
            _mp_region_descriptors[bin].push_back(mp_descriptor(*r));
            if (r->redundant) ++_mp_region_redundant_count[bin];
          }
        }
      } else if (static_cast<uint32_t>(_mp_metric[bin]) > _mp_region_max_count[bin]) {
        _mp_region_max_count[bin] = static_cast<uint32_t>(_mp_metric[bin]);
        _mp_region_max_total[bin] = static_cast<uint32_t>(_mp_total_metric[bin]);
      }
    } else if (_mp_region_start[bin] != UNDEFINED_UINT32) {
      // Close the open region (its last in-region column was position-1).
      mp_candidate_region region;
      region.seq_id = target_name(seq_id);
      region.start = _mp_region_start[bin];
      region.end = position - 1;
      region.strand = (bin == 0) ? 'F' : 'R';
      region.max_count = _mp_region_max_count[bin];
      region.max_total = _mp_region_max_total[bin];
      region.redundant = (2 * _mp_region_redundant_count[bin] > _mp_region_descriptors[bin].size());
      string joined;
      for (size_t k = 0; k < _mp_region_descriptors[bin].size(); k++) {
        if (k) joined += ";";
        joined += _mp_region_descriptors[bin][k];
      }
      region.unpaired_reads = joined;
      _mp_candidate_regions.push_back(region);

      _mp_region_start[bin] = UNDEFINED_UINT32;
      _mp_region_max_count[bin] = 0;
      _mp_region_max_total[bin] = 0;
      _mp_region_descriptors[bin].clear();
      _mp_region_redundant_count[bin] = 0;
    }
  }
}


/*! Turn the accumulated sliding-window sufficient statistics into the MP null.

 Same estimator as soft_clipping.cpp's, for the same reason: the quantity being tested is a count out
 of a local total, its genome-wide rate is not zero, and that rate is spread far too unevenly along a
 real genome for a plain binomial to describe. p0 is the pooled rate; rho is the intra-class
 correlation implied by a Pearson chi-square, which inflates the variance by (1 + (n-1)*rho).

 One difference from SC worth knowing. SC's fit has one observation per reference position and they
 are very nearly independent. MP's columns overlap almost completely -- adjacent columns share all
 but a couple of reads -- so tested_columns counts roughly W-fold redundant observations, and the
 effective sample size behind rho is nearer reference length / W. phi is still a consistent estimator
 (correlation inflates its variance, not its expectation), but it is noisier than the column count
 suggests, which is why --missing-pair-maximum-dispersion exists and must be honoured. Do not read
 tested_columns as a sample size.
*/
void identify_mutations_pileup::fill_mp_summary(MissingPairSummary& s) const
{
  s.window_width = _mp_window_width;
  s.seed_fraction_used = _mp_seed_fraction;
  s.regions_seeded = _mp_candidate_regions.size();
  s.score_cutoff = _settings.missing_pair_log10_e_value_cutoff;

  double p0_raw = (_mp_S_n > 0.0) ? (static_cast<double>(_mp_total_k) / _mp_S_n) : 0.0;
  double p0 = p0_raw;
  if ((_settings.missing_pair_minimum_rate > 0.0) && (p0 < _settings.missing_pair_minimum_rate))
    p0 = _settings.missing_pair_minimum_rate;

  uint64_t N_prime = (_mp_columns > _mp_columns_trimmed) ? (_mp_columns - _mp_columns_trimmed) : 0;
  double S_n_prime = _mp_S_n - _mp_S_n_trimmed;
  double mean_n_prime = (N_prime > 0) ? (S_n_prime / static_cast<double>(N_prime)) : 0.0;

  double phi = 0.0;
  double rho = 0.0;
  if ((N_prime > 1) && (mean_n_prime > 1.0) && (p0 > 0.0) && (p0 < 1.0)) {
    double chi2 = (_mp_S_k2_over_n - 2.0 * p0 * _mp_S_k + p0 * p0 * S_n_prime) / (p0 * (1.0 - p0));
    phi = chi2 / static_cast<double>(N_prime - 1);
    rho = (phi - 1.0) / (mean_n_prime - 1.0);
  }
  double rho_raw = rho;
  if (!std::isfinite(rho) || (rho < 0.0)) rho = 0.0;          // degenerate -> plain binomial
  if ((_settings.missing_pair_maximum_dispersion > 0.0) &&
      (rho > _settings.missing_pair_maximum_dispersion))
    rho = _settings.missing_pair_maximum_dispersion;
  if (rho < 1e-6) rho = 0.0;                                  // numerically identical to binomial

  s.null_rate = p0;
  s.null_rate_raw = p0_raw;
  s.dispersion = rho;
  s.dispersion_raw = rho_raw;
  s.pearson_phi = phi;
  s.tested_columns = N_prime;
  s.trimmed_columns = _mp_columns_trimmed;
  s.mean_tested_reads = mean_n_prime;

  // The MP statistic is a sliding window of width W, so adjacent columns share nearly all of their
  // reads and the window only turns over completely every W bases: the reference offers about
  // length/W independent chances per strand, not one per base. Using 2 * genome length -- SC's
  // count, which is right for SC because its statistic really is per-base -- would over-correct by
  // log10(W), enough on a kilobase-scale library to reject the real calls along with the noise.
  // Same reasoning as PD's n_effective_tests = reference length / mean covering gap.
  uint64_t total_length = 0;
  for (uint32_t i = 0; i < num_targets(); i++) total_length += target_length(i);
  s.n_effective_tests = (_mp_window_width > 0.0)
    ? (2.0 * static_cast<double>(total_length) / _mp_window_width) : 1.0;
  if (s.n_effective_tests < 1.0) s.n_effective_tests = 1.0;
}


/*! Write all accumulated MP candidate regions to a CSV in one pass (after the pileup completes). */
void identify_mutations_pileup::write_mp_candidate_regions(const string& filename)
{
  MissingPairSummary calibration;
  fill_mp_summary(calibration);

  // Sort by (seq_id, start, strand) so the two shoulders of one breakpoint appear adjacent.
  sort(_mp_candidate_regions.begin(), _mp_candidate_regions.end(),
       [](const mp_candidate_region& a, const mp_candidate_region& b) {
         if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
         if (a.start != b.start) return a.start < b.start;
         return a.strand < b.strand;
       });

  ofstream out(filename.c_str());
  ASSERT(!out.fail(), "Could not open output file: " + filename);

  // Calibration first, as `#key=value` lines, exactly as write_pd_candidate_regions does. This CSV
  // is the only channel from the pileup to the MP prediction step, and the null is not recoverable
  // from the regions alone -- it is measured over the columns where nothing was seeded.
  //
  // The leading format token exists so a CSV written by an older binary fails loudly rather than
  // being silently reinterpreted. Bump it whenever the meaning of a column or of the null changes.
  out << std::fixed << std::setprecision(8);
  out << "#mp_format=1" << endl;
  out << "#window_width=" << calibration.window_width << endl;
  out << "#n_effective_tests=" << calibration.n_effective_tests << endl;
  out << "#null_rate=" << calibration.null_rate << endl;
  out << "#null_rate_raw=" << calibration.null_rate_raw << endl;
  out << "#dispersion=" << calibration.dispersion << endl;
  out << "#dispersion_raw=" << calibration.dispersion_raw << endl;
  out << "#pearson_phi=" << calibration.pearson_phi << endl;
  out << "#tested_columns=" << calibration.tested_columns << endl;
  out << "#trimmed_columns=" << calibration.trimmed_columns << endl;
  out << "#mean_tested_reads=" << calibration.mean_tested_reads << endl;
  out << "#seed_fraction=" << calibration.seed_fraction_used << endl;
  out.unsetf(std::ios::fixed);

  // 'unpaired_reads' must stay last: it is ';'-joined and parsed as the final field by mp_evidence.
  // Each entry is <read_start>__<read_end> (see mp_descriptor). 'redundant' is 1 when a majority of
  // this region's reads mapped redundantly.
  out << "seq_id,start,end,strand,length,max_unpaired_count,max_window_total_count,redundant,unpaired_reads" << endl;
  for (vector<mp_candidate_region>::const_iterator r = _mp_candidate_regions.begin(); r != _mp_candidate_regions.end(); r++) {
    out << r->seq_id << "," << r->start << "," << r->end << "," << r->strand << ","
        << (r->end - r->start + 1) << "," << r->max_count << "," << r->max_total << ","
        << (r->redundant ? 1 : 0) << "," << r->unpaired_reads << endl;
  }
  out.close();

  cerr << "  MP seeded " << _mp_candidate_regions.size() << " candidate regions"
       << " (mate-unmapped rate " << std::scientific << std::setprecision(3) << calibration.null_rate;
  if (calibration.null_rate != calibration.null_rate_raw)
    cerr << ", raised from " << calibration.null_rate_raw << " by --missing-pair-minimum-rate";
  cerr << std::fixed << std::setprecision(5) << ", dispersion " << calibration.dispersion;
  if (calibration.dispersion != calibration.dispersion_raw)
    cerr << " (clamped from " << calibration.dispersion_raw << ")";
  cerr << std::setprecision(0) << ", " << calibration.n_effective_tests << " effective tests)." << endl;
  cerr.unsetf(std::ios::fixed | std::ios::scientific);
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
  // redundant (tie-broken multicopy side); 0 otherwise. Each 'discordant_pairs' entry is
  // <read1>__<read2>__<insert_size>__<read_start>__<read_end> (see dp_descriptor).
  out << "seq_id,start,end,strand,orientation,length,max_discordant_count,redundant,discordant_pairs" << endl;
  for (vector<dp_candidate_region>::const_iterator r = _dp_candidate_regions.begin(); r != _dp_candidate_regions.end(); r++) {
    out << r->seq_id << "," << r->start << "," << r->end << "," << r->strand << "," << r->orientation << ","
        << (r->end - r->start + 1) << "," << r->max_count << "," << (r->redundant ? 1 : 0) << "," << r->discordant_pairs << endl;
  }
  out.close();
}


/*! Add one read pair's gap interval [lo,hi] (between-base coordinates) to the PD ring.

 Four parallel difference arrays share one circular buffer: the covering count, the two tail counts,
 and the sum of distance quantiles. `bin` is kPDbinLong / kPDbinShort, or -1 for a pair in neither
 tail -- which still contributes to the covering count and the quantile sum, because those are the
 denominator and the statistic the seed test is built from.

 The ASSERT is the load-bearing check on the whole scheme. Writes must land in slots the per-column
 drain has not yet reached and will reach before the ring wraps: lo is at or after the next column to
 be drained, and hi+1 is no more than _pd_max_span past it. --pair-distance-maximum-span is what
 guarantees the upper bound; without it a long-range pair would silently corrupt a live slot.
 */
void identify_mutations_pileup::pd_add_interval(int32_t lo, int32_t hi, int bin, int64_t u_scaled)
{
  int32_t next_b = static_cast<int32_t>(_pd_last_b) + 1;
  ASSERT((lo >= next_b) && (hi + 1 <= next_b + _pd_max_span),
         "Pair-distance ring write out of range: [" + to_string(lo) + "," + to_string(hi) +
         "] with next column " + to_string(next_b) + " and span " + to_string(_pd_max_span));

  int32_t i_lo = lo % _pd_ring_w;
  int32_t i_hi = (hi + 1) % _pd_ring_w;

  _pd_ring_n[i_lo] += 1;
  _pd_ring_n[i_hi] -= 1;
  _pd_ring_u[i_lo] += u_scaled;
  _pd_ring_u[i_hi] -= u_scaled;

  if (bin == kPDbinLong) {
    _pd_ring_long[i_lo] += 1;
    _pd_ring_long[i_hi] -= 1;
  } else if (bin == kPDbinShort) {
    _pd_ring_short[i_lo] += 1;
    _pd_ring_short[i_hi] -= 1;
  }
}


/*! Helper method to track pair-distance (PD) candidate regions.

 Called once per reference column (including empty columns), and once past the end of each reference
 sequence (position = target_length+1) to flush an open region.

 Where DP and MP keep a deque of individually-anomalous reads and threshold its size, PD has to test
 a DISTRIBUTION: the pairs whose gap covers this column are collectively shifted, or they are not.
 So the per-column quantities are maintained as a prefix sum over the difference arrays that
 pd_add_interval writes -- an interval-stabbing count in O(1) per pair, which is what makes it
 affordable to run this over every read in the pileup.

 A deque could not do this job. DP's works only because its key (a read's start) is monotone in visit
 order and its expiry is a fixed width; PD's interval bounds are a read's END and its MATE's start,
 which are neither.

 Two gates must both pass for a column to seed, and they fail in opposite directions:

   - the rank-sum z over ALL covering pairs. Under the null each pair's distance quantile is uniform,
     so the mean quantile has mean 1/2 and variance 1/(12n); a coherent shift of even a fraction of a
     standard deviation moves it decisively once n is in the dozens. This is the gate that catches
     what DP cannot -- the shift that makes no individual pair an outlier. Without it, three pairs
     that happen to land in a tail would seed anywhere.

   - a raw count of pairs in the matching tail. A mild systematic bias (coverage, GC, mapping) can
     drag the mean quantile off 1/2 across a whole region without any pair being unusual; requiring
     real tail mass keeps that from seeding.

 @JEB expects 1-indexed positions.
 */
void identify_mutations_pileup::check_pair_distance_completion(uint32_t seq_id, uint32_t position)
{
  if (!_pd_enabled) return;

  // Drain every column up to and including this one. This is written as an explicit loop rather
  // than a single step because a prefix sum is only correct if no column is ever skipped, and
  // pileup_base::handle_position can in principle skip columns (clipping, downsampling). Clearing
  // each slot as it is consumed is what lets the buffer wrap safely.
  for (uint32_t b = _pd_last_b + 1; b <= position; b++) {
    int32_t i = static_cast<int32_t>(b % _pd_ring_w);
    _pd_n     += _pd_ring_n[i];      _pd_ring_n[i] = 0;
    _pd_long  += _pd_ring_long[i];   _pd_ring_long[i] = 0;
    _pd_short += _pd_ring_short[i];  _pd_ring_short[i] = 0;
    _pd_u     += _pd_ring_u[i];      _pd_ring_u[i] = 0;
  }
  _pd_last_b = position;

  ASSERT((_pd_n >= 0) && (_pd_long >= 0) && (_pd_short >= 0),
         "Pair-distance running count went negative at position " + to_string(position) +
         " (n=" + to_string(_pd_n) + ", long=" + to_string(_pd_long) + ", short=" + to_string(_pd_short) + ")");

  // Rank-sum z of the covering pairs' quantiles. The accumulation is integer so it is bit-identical
  // across platforms; only this comparison is floating point, and its inputs are exact integers.
  double z = 0.0;
  if (_pd_n > 0) {
    double mean_offset = (static_cast<double>(_pd_u) - static_cast<double>(_pd_n) * static_cast<double>(kPDuScale) / 2.0)
                         / static_cast<double>(kPDuScale);
    z = mean_offset * sqrt(12.0 / static_cast<double>(_pd_n));

    // Record it, so the spread of z over the whole genome can be measured and the threshold
    // corrected for it once the pileup is done. See pd_z_inflation(). Every scored column counts,
    // including the ones that go on to seed a region: real events are far too rare to move an
    // interquartile range, and excluding them would be the circular step of using the threshold to
    // decide what may inform the threshold.
    if (_pd_n >= kPDzHistMinCovering) {
      int bin = static_cast<int>(floor(z * kPDzHistPerUnit)) + kPDzHistLimit * kPDzHistPerUnit;
      if (bin < 0) bin = 0;
      if (bin >= kPDzHistBins) bin = kPDzHistBins - 1;
      _pd_z_hist[bin]++;
    }
  }

  for (int bin = 0; bin < kPDnBins; bin++) {
    int32_t tail = (bin == kPDbinLong) ? _pd_long : _pd_short;
    bool right_way = (bin == kPDbinLong) ? (z >= _pd_z_seed) : (z <= -_pd_z_seed);

    // Never "above" once past the end of the sequence: closes an open region at the flush sentinel
    // and prevents opening a spurious region that would never close.
    bool above = (tail >= _pd_seed) && right_way && (position <= target_length(seq_id));

    // The same seed test at a ladder of thresholds, counting each excursion once as it closes. This
    // measures how many INDEPENDENT chances to seed the reference actually offers -- the multiple
    // testing correction PD needs -- without assuming either a correlation length or a tail shape.
    // The gate is identical to the real one apart from the threshold, so the objects counted are the
    // same kind of object. See pd_fit_region_tail().
    {
      bool in_sequence = (position <= target_length(seq_id));
      bool tail_ok = (tail >= _pd_seed) && in_sequence;
      double directed_z = (bin == kPDbinLong) ? z : -z;
      for (int step = 0; step < kPDLadderSteps; step++) {
        bool rung_above = tail_ok && (directed_z >= pd_ladder_z(step));
        if (rung_above) {
          _pd_ladder_open[bin][step] = true;
        } else if (_pd_ladder_open[bin][step]) {
          _pd_ladder_count[bin][step]++;
          _pd_ladder_open[bin][step] = false;
        }
      }
    }

    if (above) {
      if (_pd_region_start[bin] == UNDEFINED_UINT32) {
        _pd_region_start[bin] = position;
        _pd_region_peak_position[bin] = position;
        _pd_region_peak_covering[bin] = _pd_n;
        _pd_region_peak_tail[bin] = tail;
        _pd_region_peak_z[bin] = z;
      } else if (fabs(z) > fabs(_pd_region_peak_z[bin])) {
        // Track the strongest column rather than the widest: prediction starts its BAM rescan from
        // this one point, and the region's own span is a by-product of where the seed test happens
        // to lapse, which can sit well away from the breakpoint.
        _pd_region_peak_position[bin] = position;
        _pd_region_peak_covering[bin] = _pd_n;
        _pd_region_peak_tail[bin] = tail;
        _pd_region_peak_z[bin] = z;
      }
    } else if (_pd_region_start[bin] != UNDEFINED_UINT32) {
      // Close the open region (its last in-region column was position-1).
      pd_candidate_region region;
      region.seq_id = target_name(seq_id);
      region.start = _pd_region_start[bin];
      region.end = position - 1;
      region.direction = (bin == kPDbinLong) ? 'L' : 'S';
      region.peak_position = _pd_region_peak_position[bin];
      region.peak_covering = _pd_region_peak_covering[bin];
      region.peak_tail = _pd_region_peak_tail[bin];
      region.peak_z = _pd_region_peak_z[bin];
      _pd_candidate_regions.push_back(region);

      _pd_region_start[bin] = UNDEFINED_UINT32;
      _pd_region_peak_position[bin] = 0;
      _pd_region_peak_covering[bin] = 0;
      _pd_region_peak_tail[bin] = 0;
      _pd_region_peak_z[bin] = 0.0;
    }
  }
}


/*! Write all accumulated PD candidate regions to a CSV in one pass (after the pileup completes).

 Unlike the DP and MP CSVs this carries no per-pair descriptor list. Those exist because a DP or MP
 rescan window can miss reads the seed saw; PD's cannot, since every pair whose gap covers the peak
 has its leftmost mate within one maximum span of it, so the rescan is guaranteed to re-derive the
 same set. A difference array has no membership list to write out in any case.
 */
double identify_mutations_pileup::pd_z_inflation() const
{
  if (_pd_z_hist.empty()) return 1.0;

  uint64_t total = 0;
  for (size_t i = 0; i < _pd_z_hist.size(); i++) total += _pd_z_hist[i];
  // Below this there is no useful quantile to read, and a whole-genome pileup is never this small
  // unless PD was effectively inactive anyway.
  if (total < 1000) return 1.0;

  // Interquartile range, read off the histogram. Bin centres, so a half-bin offset either side.
  uint64_t want_lo = total / 4, want_hi = (total * 3) / 4;
  double q_lo = 0.0, q_hi = 0.0;
  uint64_t run = 0;
  bool have_lo = false, have_hi = false;
  for (size_t i = 0; i < _pd_z_hist.size(); i++) {
    run += _pd_z_hist[i];
    double centre = (static_cast<double>(i) + 0.5) / kPDzHistPerUnit - kPDzHistLimit;
    if (!have_lo && (run >= want_lo)) { q_lo = centre; have_lo = true; }
    if (!have_hi && (run >= want_hi)) { q_hi = centre; have_hi = true; break; }
  }
  if (!have_lo || !have_hi) return 1.0;

  // 1.34898 = the IQR of a standard normal, so this is the sd the null WOULD have if it were normal
  // with this spread. Under the model the seed assumes, it is 1.
  double sd = (q_hi - q_lo) / 1.34898038;
  if (!(sd > 1.0)) return 1.0;
  return sd;
}

bool identify_mutations_pileup::pd_fit_region_tail(double& out_log_intercept, double& out_sigma,
                                                   double& out_fit_lo, double& out_fit_hi) const
{
  out_log_intercept = 0.0;
  out_sigma = 1.0;
  out_fit_lo = 0.0;
  out_fit_hi = 0.0;
  if (!_pd_enabled) return false;

  // Below this a rung's count is too small to distinguish Poisson wobble, and too easily a
  // meaningful fraction real events rather than noise. Both would bias the fitted tail.
  const uint64_t kPDFitMinRegions = 20;
  const int      kPDFitMinRungs   = 3;
  // ... and above this the ladder is not yet a tail. Its lowest rungs are governed by the seed's
  // tail-COUNT gate rather than by z, so they sit on a plateau -- 1061, 981, 851, 744 regions on the
  // library this was written against, against 104 at |z| = 3. Fitting a Gaussian exponent through
  // that plateau mixes bulk with tail and reports a spread that is neither. Start the window a factor
  // of kPDFitTailFactor below the busiest rung instead, which is where the count has actually begun
  // to fall like a tail.
  const uint64_t kPDFitTailFactor = 3;

  uint64_t max_count = 0;
  for (int step = 0; step < kPDLadderSteps; step++) {
    uint64_t count = 0;
    for (int bin = 0; bin < kPDnBins; bin++) count += _pd_ladder_count[bin][step];
    max_count = max(max_count, count);
  }

  // x = t^2 and y = log R(t): the model log R = a - t^2 / (2 sigma^2) is linear in x. Weight by the
  // count, since var(log R) is about 1/R for a Poisson count.
  double sw = 0.0, swx = 0.0, swy = 0.0, swxx = 0.0, swxy = 0.0;
  int rungs = 0;
  for (int step = 0; step < kPDLadderSteps; step++) {
    uint64_t count = 0;
    for (int bin = 0; bin < kPDnBins; bin++) count += _pd_ladder_count[bin][step];
    if (count < kPDFitMinRegions) continue;
    if (count * kPDFitTailFactor > max_count) continue;
    double t = pd_ladder_z(step);
    double x = t * t;
    double y = log(static_cast<double>(count));
    double w = static_cast<double>(count);
    sw += w; swx += w * x; swy += w * y; swxx += w * x * x; swxy += w * x * y;
    if (rungs == 0) out_fit_lo = t;
    out_fit_hi = t;
    rungs++;
  }
  if (rungs < kPDFitMinRungs) return false;

  double denom = sw * swxx - swx * swx;
  if (!(fabs(denom) > 0.0)) return false;
  double slope = (sw * swxy - swx * swy) / denom;
  double intercept = (swy - slope * swx) / sw;

  // slope = -1 / (2 sigma^2), so a non-negative slope means R(t) does not fall with the threshold at
  // all -- not a tail, and nothing to extrapolate from.
  if (!(slope < 0.0)) return false;
  double sigma = sqrt(-1.0 / (2.0 * slope));
  if (!(sigma > 0.0) || !isfinite(sigma)) return false;

  out_log_intercept = intercept;
  out_sigma = sigma;
  return true;
}

void identify_mutations_pileup::fill_pd_summary(PairDistanceSummary& s) const
{
  s.mean_covering_gap = _pd_mean_covering_gap;
  s.n_effective_tests = (_pd_mean_covering_gap > 0.0)
    ? static_cast<double>(_total_reference_length) / _pd_mean_covering_gap : 0.0;
  s.seed_z_used = _pd_z_seed;
  s.z_iqr_inflation = pd_z_inflation();
  s.regions_seeded = _pd_candidate_regions.size();
  s.pairs_considered = _pd_pairs_considered;
  s.pairs_ambiguous_placement = _pd_pairs_ambiguous_placement;
  s.region_tail_fit_ok = pd_fit_region_tail(s.region_tail_log_intercept, s.region_tail_sigma,
                                            s.region_tail_fit_lo, s.region_tail_fit_hi);
}

void identify_mutations_pileup::write_pd_candidate_regions(const string& filename)
{
  // The regions are NOT filtered by pd_z_inflation() any more, and the seed z they cleared is
  // deliberately low. Both changes serve the same end: what decides a PD call is now a genome-wide
  // e-value computed in pd_evidence, and the calibration that e-value needs is fitted to these very
  // regions. Throwing the weak ones away here would be discarding the null sample. It would also be
  // inconsistent with the ladder counts written below, which are measured on the uncorrected z.
  //
  // pd_z_inflation() survives as a reported diagnostic only. It is an interquartile range, so it
  // describes the bulk of the z distribution and is blind to the tail -- on the library that exposed
  // this it read 1.04x while the tail was about 2.5x too heavy for a normal. pd_fit_region_tail()
  // measures the tail itself.
  PairDistanceSummary calibration;
  fill_pd_summary(calibration);

  // Seeding loosely costs a BAM rescan per region, and on an over-dispersed library -- a mate-pair
  // run, typically -- the low rungs can hold tens of thousands of them. The ladder already says how
  // many regions each threshold yields, so raise the threshold to the lowest rung that stays inside
  // the budget rather than guessing one up front. Several hundred regions is ample to fit a tail to;
  // the rest only cost time. Whatever threshold this settles on is the one the calibration is
  // conditioned on, since the retained regions are both the null sample and the tested population.
  {
    const uint64_t kPDMaxCandidateRegions = 3000;
    double capped_z = calibration.seed_z_used;
    for (int step = 0; step < kPDLadderSteps; step++) {
      double t = pd_ladder_z(step);
      if (t < calibration.seed_z_used) continue;
      uint64_t count = 0;
      for (int bin = 0; bin < kPDnBins; bin++) count += _pd_ladder_count[bin][step];
      capped_z = t;
      if (count <= kPDMaxCandidateRegions) break;
    }
    if (capped_z > calibration.seed_z_used) {
      size_t before = _pd_candidate_regions.size();
      vector<pd_candidate_region> kept;
      kept.reserve(_pd_candidate_regions.size());
      for (vector<pd_candidate_region>::const_iterator r = _pd_candidate_regions.begin(); r != _pd_candidate_regions.end(); r++)
        if (fabs(r->peak_z) >= capped_z) kept.push_back(*r);
      _pd_candidate_regions.swap(kept);
      cerr << "  PD seeded " << before << " candidate regions at |z| >= "
           << std::fixed << std::setprecision(2) << calibration.seed_z_used
           << ", more than the rescan budget; raising the seed to " << capped_z
           << " leaves " << _pd_candidate_regions.size() << "." << endl;
      cerr.unsetf(std::ios::fixed);
      calibration.seed_z_used = capped_z;
      calibration.regions_seeded = _pd_candidate_regions.size();
      _pd_z_seed = capped_z;
    }
  }

  // Sort by (seq_id, start, direction) so a nearby deletion and insertion appear adjacent.
  sort(_pd_candidate_regions.begin(), _pd_candidate_regions.end(),
       [](const pd_candidate_region& a, const pd_candidate_region& b) {
         if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
         if (a.start != b.start) return a.start < b.start;
         return a.direction < b.direction;
       });

  ofstream out(filename.c_str());
  ASSERT(!out.fail(), "Could not open output file: " + filename);

  // Calibration first, as `#key=value` lines. This CSV is the only channel from the pileup to the PD
  // prediction step, and none of this is recoverable from the regions alone.
  out << std::fixed << std::setprecision(6);
  out << "#mean_covering_gap=" << calibration.mean_covering_gap << endl;
  out << "#n_effective_tests=" << calibration.n_effective_tests << endl;
  out << "#seed_z=" << calibration.seed_z_used << endl;
  out << "#pairs_considered=" << calibration.pairs_considered << endl;
  out << "#pairs_ambiguous_placement=" << calibration.pairs_ambiguous_placement << endl;
  out << "#z_iqr_inflation=" << calibration.z_iqr_inflation << endl;
  out << "#region_tail_fit_ok=" << (calibration.region_tail_fit_ok ? 1 : 0) << endl;
  out << "#region_tail_log_intercept=" << calibration.region_tail_log_intercept << endl;
  out << "#region_tail_sigma=" << calibration.region_tail_sigma << endl;
  out << "#region_tail_fit_lo=" << calibration.region_tail_fit_lo << endl;
  out << "#region_tail_fit_hi=" << calibration.region_tail_fit_hi << endl;
  // The raw ladder, so the fit can be second-guessed from the output alone.
  out << "#region_ladder=";
  for (int step = 0; step < kPDLadderSteps; step++) {
    uint64_t count = 0;
    for (int bin = 0; bin < kPDnBins; bin++) count += _pd_ladder_count[bin][step];
    if (step) out << ";";
    out << std::setprecision(2) << pd_ladder_z(step) << ":" << count;
  }
  out << endl;

  // 'direction' is L (pairs mapping farther apart than expected -- deletion-like) or S (closer --
  // insertion-like). The peak_* columns describe the single strongest column in the region, which
  // is where pd_evidence starts its rescan.
  out << "seq_id,start,end,direction,length,peak_position,peak_covering_count,peak_tail_count,peak_z" << endl;
  out << std::setprecision(4);
  for (vector<pd_candidate_region>::const_iterator r = _pd_candidate_regions.begin(); r != _pd_candidate_regions.end(); r++) {
    out << r->seq_id << "," << r->start << "," << r->end << "," << r->direction << ","
        << (r->end - r->start + 1) << "," << r->peak_position << "," << r->peak_covering << ","
        << r->peak_tail << "," << r->peak_z << endl;
  }
  out.close();

  cerr << "  PD seeded " << _pd_candidate_regions.size() << " candidate regions at |z| >= "
       << std::fixed << std::setprecision(2) << calibration.seed_z_used
       << " (mean covering gap " << std::setprecision(0) << calibration.mean_covering_gap
       << " bases, " << calibration.n_effective_tests << " effective tests";
  if (calibration.region_tail_fit_ok)
    cerr << ", null tail sigma " << std::setprecision(2) << calibration.region_tail_sigma;
  cerr << ")." << endl;
  cerr.unsetf(std::ios::fixed);
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

  mut["ks_quality_p_value"]    = formatted_double(ks_quality_p_value,    5, true).to_string();
  mut["fisher_strand_p_value"] = formatted_double(fisher_strand_p_value, 5, true).to_string();
}


/*! Report a fitted frequency, flooring an allele the fit places below the half-read level to zero.

  EM approaches an absent component asymptotically instead of reaching it, so an allele with no
  supporting reads at all settles around 1e-13 rather than 0. Printing that is just noise, and it
  is the same "below half a read is not a called allele" rule next_index() applies.
*/
double allele_model::reported_frequency(base_index b) const
{
  if ((n == 0) || (b >= base_list_size)) return 0.0;
  return (f[b] < 0.5 / static_cast<double>(n)) ? 0.0 : f[b];
}

string allele_model::spectrum_string(uint32_t precision_places) const
{
  string s;
  for (uint8_t b=0; b<base_list_size; b++) {
    double fr = reported_frequency(b);
    if (fr <= 0.0) continue;
    if (!s.empty()) s += ",";
    s += string(1, baseindex2char(b)) + ":" + formatted_double(fr, precision_places, true).to_string();
  }
  return s;
}

base_index allele_model::major_index() const
{
  if (n == 0) return base_list_N_index;
  base_index best = 0;
  for (uint8_t b=1; b<base_list_size; b++) {
    if (f[b] > f[best]) best = b;
  }
  return (f[best] > 0.0) ? best : base_list_N_index;
}

base_index allele_model::next_index(base_index exclude_index) const
{
  if (n == 0) return base_list_N_index;

  // An allele the fit puts below the half-read level is not a called allele. EM approaches an
  // absent component asymptotically rather than reaching zero, so without a floor every position
  // would report all five alleles as present at frequencies like 1e-9.
  const double present_threshold = 0.5 / static_cast<double>(n);

  base_index best = base_list_N_index;
  for (uint8_t b=0; b<base_list_size; b++) {
    if (b == exclude_index) continue;
    if (f[b] < present_threshold) continue;
    if ((best == base_list_N_index) || (f[b] > f[best])) best = b;
  }
  return best;
}

/*! Maximum log10 likelihood with one allele's frequency held fixed.

  The same EM as the free fit, except that f[variant_index] is pinned and the remaining alleles are
  renormalized to share what is left. Warm-started from the free fit, so it usually converges in a
  handful of iterations.
*/
double identify_mutations_pileup::profile_log10_likelihood(const vector<polymorphism_data>& pdata, const allele_model& full, base_index variant_index, double f_fixed) const
{
  if ((full.n == 0) || (variant_index >= base_list_size)) return 0.0;

  double f[base_list_size];
  double other_total = 0.0;
  for (uint8_t b=0; b<base_list_size; b++) { if (b != variant_index) other_total += full.f[b]; }
  for (uint8_t b=0; b<base_list_size; b++) {
    if (b == variant_index) f[b] = f_fixed;
    else f[b] = (other_total > 0.0) ? (1.0 - f_fixed) * full.f[b] / other_total
                                    : (1.0 - f_fixed) / static_cast<double>(base_list_size - 1);
  }

  const uint32_t k_max_iterations = 50;
  double log10_likelihood = 0.0;

  for (uint32_t iter = 0; iter < k_max_iterations; iter++) {

    double sum_w[base_list_size];
    for (uint8_t b=0; b<base_list_size; b++) sum_w[b] = 0.0;
    log10_likelihood = 0.0;

    for (vector<polymorphism_data>::const_iterator it=pdata.begin(); it!=pdata.end(); ++it) {
      double s = 0.0;
      for (uint8_t b=0; b<base_list_size; b++) s += f[b] * it->_r[b];
      if (s > 0.0) {
        log10_likelihood += log10(s) + it->_log10_pr_max;
        for (uint8_t b=0; b<base_list_size; b++) sum_w[b] += f[b] * it->_r[b] / s;
      } else {
        for (uint8_t b=0; b<base_list_size; b++) sum_w[b] += f[b];
      }
    }

    double others = 0.0;
    for (uint8_t b=0; b<base_list_size; b++) { if (b != variant_index) others += sum_w[b]; }

    double max_delta = 0.0;
    for (uint8_t b=0; b<base_list_size; b++) {
      if (b == variant_index) continue;
      double f_new = (others > 0.0) ? (1.0 - f_fixed) * sum_w[b] / others
                                    : (1.0 - f_fixed) / static_cast<double>(base_list_size - 1);
      max_delta = max(max_delta, fabs(f_new - f[b]));
      f[b] = f_new;
    }

    if (max_delta < _polymorphism_precision_decimal) break;
  }

  return log10_likelihood;
}

/*! Write the variant frequency's confidence bounds onto an RA entry.

  A profile-likelihood interval on the variant allele's frequency: the set of f where holding the
  allele at f and re-fitting everything else costs less than a fixed amount of log likelihood.

    { f : log10 L_max - log10 L_profile(f) <= kProfileLikelihoodLog10Drop }

  This is the interval the error model itself implies. Every read enters through the calibrated
  per-base likelihoods, so a quality-2 base widens the interval and a quality-40 base narrows it
  because of what they say about the data, not because of a weight assigned to them beforehand.

  It replaces a binomial bound fed an "effective depth". That construction needed a scalar count of
  how many reads "really" informed the call, and there is no non-arbitrary way to produce one: the
  derivation from likelihood curvature is degenerate at f -> 1, which is where fixed variants sit,
  and the discrimination heuristic that worked at both ends was not derived from anything. Since a
  likelihood interval needs no such count, the whole question disappears -- and it handles the
  boundary natively, becoming one-sided when the fit is pinned at 0 or 1.

  The bounds are recorded rather than recomputed downstream from text, so the interval that decides
  the call is the one the .gd and the report show. Previously they were derived during
  classification and discarded, which meant a rejection reading "FREQUENCY_CUTOFF at frequency
  0.43" could not be explained from the output at all.
*/

// Each endpoint is a ONE-SIDED 95% bound, matching the Clopper-Pearson bounds this replaces (and
// the ones JC still uses). The interval { f : log10 L_max - log10 L(f) <= D } has two-sided
// coverage P(chi2_1 <= 2*ln(10)*D), so one-sided 95% per side means 90% two-sided:
//   2*ln(10)*D = chi2_1(0.90) = 2.705543  ->  D = 0.587566
static const double kProfileLikelihoodLog10Drop = 0.587566;

void identify_mutations_pileup::write_RA_frequency_bounds(cDiffEntry& mut, const vector<polymorphism_data>& pdata, const allele_model& amodel, base_index variant_index) const
{
  double lower = 0.0, upper = 1.0;

  if ((amodel.n > 0) && (variant_index < base_list_size)) {

    const double f_hat = amodel.f[variant_index];

    // Evaluate the peak through the same routine as every other point, so that the profile is
    // self-consistent and pl(f_hat) is exactly the maximum the threshold is measured from.
    const double pl_max = profile_log10_likelihood(pdata, amodel, variant_index, f_hat);
    const double target = pl_max - kProfileLikelihoodLog10Drop;

    // Lower endpoint: the smallest f whose profile still clears the threshold.
    if (profile_log10_likelihood(pdata, amodel, variant_index, 0.0) >= target) {
      lower = 0.0;
    } else {
      double lo = 0.0, hi = f_hat;
      for (uint32_t i = 0; (i < 40) && ((hi - lo) > _polymorphism_precision_decimal); i++) {
        double mid = 0.5 * (lo + hi);
        if (profile_log10_likelihood(pdata, amodel, variant_index, mid) >= target) hi = mid;
        else                                                                       lo = mid;
      }
      lower = hi;
    }

    // Upper endpoint: the largest such f.
    if (profile_log10_likelihood(pdata, amodel, variant_index, 1.0) >= target) {
      upper = 1.0;
    } else {
      double lo = f_hat, hi = 1.0;
      for (uint32_t i = 0; (i < 40) && ((hi - lo) > _polymorphism_precision_decimal); i++) {
        double mid = 0.5 * (lo + hi);
        if (profile_log10_likelihood(pdata, amodel, variant_index, mid) >= target) lo = mid;
        else                                                                       hi = mid;
      }
      upper = lo;
    }
  }

  mut[FREQUENCY_LOWER] = formatted_double(lower, _polymorphism_precision_places, true).to_string();
  mut[FREQUENCY_UPPER] = formatted_double(upper, _polymorphism_precision_places, true).to_string();
}

/*! Fit the allele frequency mixture at one position by EM.

  The position carries frequencies f over the five candidate alleles. Read i draws a true base b
  with probability f[b] and is then observed with probability P_i(b), which the error table already
  gave us (cached in polymorphism_data by fill_read_base_likelihoods):

    E-step:  w_i(b) = f[b] r_i(b) / s_i,   s_i = Sum_b' f[b'] r_i(b')
    M-step:  f[b]   = (1/n) Sum_i w_i(b)

  Fitting all five alleles at once, rather than a mixture of the top two by coverage, is the whole
  point: it puts every frequency on the same denominator (total depth) and keeps the reference
  allele in the model even when it is not one of the two best-supported bases. That case is not
  exotic -- across the test goldens, 44% of RA entries have neither the major nor the minor allele
  equal to the reference.

  The objective is concave on the simplex, so the fixed point EM reaches is the global maximum.
  That is what lets this replace a bisection that needed an explicit local-minimum guard.

  Everything is evaluated on the max-normalized r_i(b) rather than raw probabilities, so the
  mixture sums cannot underflow; the discarded exponent comes back per read via _log10_pr_max.
*/
allele_model
identify_mutations_pileup::fit_allele_frequencies(const vector<polymorphism_data>& pdata, const bool allowed[base_list_size]) const
{
  allele_model m;
  m.n = static_cast<uint32_t>(pdata.size());
  if (m.n == 0) return m;

  // Start from Laplace-smoothed observed base counts. The half-count keeps every allowed component
  // strictly interior: a component initialized at exactly zero is a fixed point of the E-step and
  // could never be revived, so a hard zero would decide the answer instead of the data.
  double init_total = 0.0;
  uint32_t n_allowed = 0;
  for (uint8_t b=0; b<base_list_size; b++) {
    if (!allowed[b]) continue;
    n_allowed++;
    m.f[b] = 0.5;
  }
  if (n_allowed == 0) return m;

  for (vector<polymorphism_data>::const_iterator it=pdata.begin(); it!=pdata.end(); ++it) {
    base_index obs = basechar2index(it->_base_char);
    if ((obs < base_list_size) && allowed[obs]) m.f[obs] += 1.0;
  }
  for (uint8_t b=0; b<base_list_size; b++) init_total += m.f[b];
  for (uint8_t b=0; b<base_list_size; b++) m.f[b] /= init_total;

  const uint32_t k_max_iterations = 50;
  const double k_tolerance = _polymorphism_precision_decimal;

  for (m.iterations = 1; m.iterations <= k_max_iterations; m.iterations++) {

    double sum_w[base_list_size];
    for (uint8_t b=0; b<base_list_size; b++) sum_w[b] = 0.0;
    double log10_likelihood = 0.0;

    for (vector<polymorphism_data>::const_iterator it=pdata.begin(); it!=pdata.end(); ++it) {

      double s = 0.0;
      for (uint8_t b=0; b<base_list_size; b++) {
        if (allowed[b]) s += m.f[b] * it->_r[b];
      }

      if (s > 0.0) {
        log10_likelihood += log10(s) + it->_log10_pr_max;
        for (uint8_t b=0; b<base_list_size; b++) {
          if (!allowed[b]) continue;
          sum_w[b] += m.f[b] * it->_r[b] / s;
        }
      } else {
        // Unreachable with a calibrated table: every entry is Laplace-smoothed and the
        // mapping-quality term adds a strictly positive floor. Fall back to the current
        // frequencies so that Sum_b Sum_i w_i(b) == n stays exact regardless.
        for (uint8_t b=0; b<base_list_size; b++) {
          if (!allowed[b]) continue;
          sum_w[b] += m.f[b];
        }
      }
    }

    double max_delta = 0.0;
    for (uint8_t b=0; b<base_list_size; b++) {
      if (!allowed[b]) continue;
      double f_new = sum_w[b] / static_cast<double>(m.n);
      max_delta = max(max_delta, fabs(f_new - m.f[b]));
      m.f[b] = f_new;
    }

    // Commit the sums that produced the frequencies we just set. Because the M-step defines
    // f[b] = sum_w[b]/n, the identity sum_w[b] == n * f[b] holds exactly for the committed pair,
    // which is what the effective depth downstream relies on.
    for (uint8_t b=0; b<base_list_size; b++) m.sum_w[b] = sum_w[b];
    m.log10_likelihood = log10_likelihood;

    if (max_delta < k_tolerance) break;
  }
  if (m.iterations > k_max_iterations) m.iterations = k_max_iterations;

  return m;
}

/*! log10 evidence that a given allele is present at this position at all.

  A profile likelihood ratio between the full fit and the best fit that holds this allele out
  entirely, Bonferroni-corrected by the reference length as every score in this file is.

  The null keeps the reference allele AND any third allele, which is what makes the number
  meaningful when the two best-supported bases are both non-reference: the question asked is
  "is THIS allele present", not "are the top two bases in different proportions".
*/
double identify_mutations_pileup::variant_presence_score(const vector<polymorphism_data>& pdata, const allele_model& full, base_index variant_index) const
{
  if ((full.n == 0) || (variant_index >= base_list_size)) return numeric_limits<double>::quiet_NaN();

  bool without[base_list_size];
  uint32_t n_remaining = 0;
  for (uint8_t b=0; b<base_list_size; b++) {
    without[b] = (b != variant_index);
    if (without[b]) n_remaining++;
  }
  if (n_remaining == 0) return numeric_limits<double>::quiet_NaN();

  allele_model null_fit = fit_allele_frequencies(pdata, without);

  return (full.log10_likelihood - null_fit.log10_likelihood) - _log10_ref_length;
}


/*! Per-read-base likelihoods under each of the five candidate true bases.

  The calibrated error table is keyed in READ-strand space: count_alignment_position()
  (error_count.cpp) complements both the reference and the observed base when the read is reversed,
  so a hypothesis about the reference base has to be complemented to match before lookup.

  The mapping-quality term mixes every hypothesis with a uniform base at the read's probability of
  being misplaced, eps = 10^(-MQ/10). It puts a floor of eps/5 under each hypothesis -- about 1.3e-5
  at bowtie2's typical MAPQ 42 -- which caps any single read base's log10 likelihood ratio near 4.9.
  That is deliberate: a read that is probably somewhere else entirely should not get to decide a
  position no matter how confident its base call is.
*/
void identify_mutations_pileup::fill_read_base_likelihoods(polymorphism_data& pd)
{
  const double incorrect_mapping_prob = pow(10, -static_cast<double>(pd._mapping_quality) / 10);
  const double correct_mapping_prob = 1 - incorrect_mapping_prob;
  const double uniform_prob = 1.0 / static_cast<double>(base_list_size);

  covariate_values_t this_cv = pd._cv;
  const bool obs_top_strand = (pd._strand == 1);
  if (!obs_top_strand) {
    this_cv.obs_base() = complement_base_index(this_cv.obs_base());
  }

  pd._log10_pr_max = -numeric_limits<double>::max();
  for (uint8_t b=0; b<base_list_size; b++) {
    this_cv.ref_base() = obs_top_strand ? b : complement_base_index(b);
    double pr = correct_mapping_prob * _error_table.get_prob(this_cv)
              + incorrect_mapping_prob * uniform_prob;
    // Floating point error can make this a very slightly negative number.
    if (pr < 0.0) pr = 0.0;
    pd._log10_pr[b] = log10(pr);
    pd._log10_pr_max = max(pd._log10_pr_max, pd._log10_pr[b]);
  }
  for (uint8_t b=0; b<base_list_size; b++) {
    pd._r[b] = pow(10, pd._log10_pr[b] - pd._log10_pr_max);
  }
}

/*! Call the single most probable pure genotype, and score it against the alternatives.

  With uniform priors over the five pure genotypes, the per-read renormalization that the old
  incremental caller applied is a factor common to every genotype, so it cancels out of both the
  argmax and the reported score. What is left is

    score = log10 P(data | best) - log10 Sum_{b != best} P(data | b)

  evaluated by log-sum-exp against the largest competing term -- which is exactly the number
  cDiscreteSNPCaller::get_prediction() produced, from five accumulated sums instead of a
  renormalized posterior vector carried across every read.
*/
cSNPCall identify_mutations_pileup::pure_genotype_call(const double log10_pr_sum[base_list_size], uint32_t observations) const
{
  cSNPCall snp_call;

  // Best base is 'N' and the score is NaN when there is nothing to call.
  if (observations == 0) return snp_call;

  base_index best = 0;
  for (uint8_t b=1; b<base_list_size; b++) {
    if (log10_pr_sum[b] > log10_pr_sum[best]) best = b;
  }
  snp_call.genotype = string(1, baseindex2char(best));

  // Offset by the largest competing genotype so the exponentials below cannot overflow.
  double log10_offset_probability = -numeric_limits<double>::max();
  for (uint8_t b=0; b<base_list_size; b++) {
    if (b != best) log10_offset_probability = max(log10_offset_probability, log10_pr_sum[b]);
  }

  // '.' (the gap state) is a genotype like any other, so it belongs in the error mass. Leaving it
  // out while the offset above included it made the score claim more confidence than the data
  // supported at exactly the positions where the competing genotype is a deletion: with '.' as the
  // runner-up every remaining term is negligible against the offset, the sum underflows to zero,
  // and the score comes back +inf. tests/bull_1/expected.gd had two such entries, one of them a
  // 52/48 G/'.' mixture reported as an infinitely confident consensus call.
  double total_error_probability = 0;
  for (uint8_t b=0; b<base_list_size; b++) {
    if (b != best) total_error_probability += pow(10, log10_pr_sum[b] - log10_offset_probability);
  }
  double log10_total_error_probability = log10(total_error_probability);
  log10_total_error_probability += log10_offset_probability;

  snp_call.score = log10_pr_sum[best] - log10_total_error_probability;

  return snp_call;
}


} // namespace breseq

