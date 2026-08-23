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
    uint64_t num_unmapped_reads;
    uint64_t num_unmapped_read_bases;
    uint64_t num_total_reads;
    uint64_t num_total_bases;
    
    ReadFileSummary()
    : num_unmapped_reads(0)
    , num_unmapped_read_bases(0)
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
    uint64_t total_unmapped_reads;
    uint64_t total_unmapped_read_bases;
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
    : total_unmapped_reads(0)
    , total_unmapped_read_bases(0)
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
    double nbinom_relative_variance;
    double average;
    double variance;
    double relative_variance;
    
    CoverageSummary()
    : deletion_coverage_propagation_cutoff(0.0)
    , deletion_coverage_seed_cutoff(0.0)
    , nbinom_size_parameter(0.0)
    , nbinom_mean_parameter(0.0)
    , nbinom_prob_parameter(0.0)
    , nbinom_variance(0.0)
    , nbinom_relative_variance(0.0)
    , average(0.0)
    , variance(0.0)
    , relative_variance(0.0)
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
    bool reads_were_split;
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
    , reads_were_split(false)
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
                 bool _reads_were_split,
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
    , reads_were_split(_reads_were_split)
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

  class PairedMappingDistanceDistributionSummary : public JSONStorable<PairedMappingDistanceDistributionSummary>
  {
  public:
    double median;
    double mad;              // upper (one-sided) median absolute deviation
    double distance_cutoff;  // median + z*mad/0.6745 (Iglewicz & Hoaglin modified z-score rule)
    string majority_orientation; // most common of FF/FR/RF by total observation count
    double mapped_pairs;     // number of read pairs (both mates mapped, concordant or discordant)
    double concordant_pairs; // number of those pairs called concordant

    PairedMappingDistanceDistributionSummary()
    : median(0.0)
    , mad(0.0)
    , distance_cutoff(0.0)
    , majority_orientation("")
    , mapped_pairs(0.0)
    , concordant_pairs(0.0)
    {}
  };

  class PairedMappingDistanceDistributionSummaries : public map<string, PairedMappingDistanceDistributionSummary>, public JSONStorable<PairedMappingDistanceDistributionSummaries>
  {
  public:
    PairedMappingDistanceDistributionSummaries() {}
  };

  class SoftClippingSummary : public JSONStorable<SoftClippingSummary>
  {
  public:
    uint64_t total_spanning_read_bases;      // read-through opportunities: reads spanning a position with
                                             // >= min_bases aligned on BOTH sides, counted once per direction
    uint64_t total_spanning_read_bases_forward;  // ...split by the strand of the spanning read. The strand
    uint64_t total_spanning_read_bases_reverse;  // test falls back to this ratio where a position has no
                                                 // read-through of its own (frequency == 1.000).
    uint64_t total_clipped_read_ends;  // total soft-clip events; a read with both ends clipped counts twice
    double   soft_clipping_rate;       // raw clip rate: total_clipped_read_ends / total opportunities

    // Null model actually used by the test, computed at tabulation time. The
    // per-position counts file only lists positions with >= 1 clip event, so the
    // zero-clip positions needed for these estimates are not recoverable later --
    // they must be carried forward here.
    uint64_t total_agreeing_clipped_read_ends;  // clip events whose tail matches the position consensus
    // Diagnostic, not part of the null: how many of those sit at positions where every agreeing
    // clipped read came from one strand. Near 100% means the run's clip population is dominated by
    // an end-of-read artifact (dark-cycle poly-G, adapter read-through), not by breakpoints.
    uint64_t total_strand_pure_agreeing_clipped_read_ends;
    double   soft_clipping_null_rate;           // p0 used (agreeing rate, after the minimum-rate floor)
    double   soft_clipping_dispersion;          // rho used; 0 => plain binomial
    double   soft_clipping_pearson_phi;         // diagnostic: Pearson chi2 / (N-1) over the fitted positions
    uint64_t soft_clipping_tested_positions;    // N' -- (position, direction) pairs entering the fit
    uint64_t soft_clipping_trimmed_positions;   // positions excluded from the fit by the trim
    double   soft_clipping_mean_tested_reads;   // mean n over the fitted positions

    SoftClippingSummary()
    : total_spanning_read_bases(0)
    , total_spanning_read_bases_forward(0)
    , total_spanning_read_bases_reverse(0)
    , total_clipped_read_ends(0)
    , soft_clipping_rate(0.0)
    , total_agreeing_clipped_read_ends(0)
    , total_strand_pure_agreeing_clipped_read_ends(0)
    , soft_clipping_null_rate(0.0)
    , soft_clipping_dispersion(0.0)
    , soft_clipping_pearson_phi(0.0)
    , soft_clipping_tested_positions(0)
    , soft_clipping_trimmed_positions(0)
    , soft_clipping_mean_tested_reads(0.0)
    {}
  };

  // Everything that decided which discordant-pair (DP) junctions were kept: the library shape the
  // tests are derived from, the cutoffs actually in force, the background fitted to this run, and the
  // resulting tally. Reported so a user can see WHY a run predicted few or many DP items -- most of
  // these are derived per run rather than fixed, so the command line alone does not tell you.
  class DiscordantPairSummary : public JSONStorable<DiscordantPairSummary>
  {
  public:
    // Library shape (the paired-mapping distance fit the DP tests are built on).
    double   read_length_avg;
    double   pair_distance_median;
    double   pair_distance_mad;
    double   pair_distance_cutoff;
    string   pair_orientation;                  // FR or RF (FF/RR are unsupported)

    // Concordant-pair crossing model: the null for the skew score, and its effect size.
    string   crossing_reference_seq_id;
    double   crossing_reference_coverage;
    double   expected_concordant_crossing;      // uncensored mean crossing at the reference = lambda
    bool     crossing_use_empirical;            // false => negative-binomial fallback
    uint64_t crossing_normal_positions;         // N inside the censor window

    // Cutoffs in force for this run.
    double   frequency_cutoff;                  // tracks the prediction mode
    double   skew_cutoff;
    double   minimum_crossing;                  // power gate on the skew test; 0 = always apply
    bool     skew_in_force;                     // expected_concordant_crossing cleared minimum_crossing
    double   background_e_value_cutoff;
    int32_t  minimum_pairs_option;              // the user's absolute floor
    int32_t  minimum_pairs_used;                // the floor actually used (raised by the background fit)

    // Spurious-pair background fitted to this run's edge-weight spectrum.
    uint64_t candidate_junctions;               // region pairs sharing at least one read pair
    double   background_mean;
    double   background_size;                   // negative binomial size; 0 => Poisson

    // Outcome.
    uint64_t items_tested;                      // DP items emitted (weight >= minimum_pairs_used)
    uint64_t items_dropped_unsupported;         // cleared the floor, but no pair bridged the placed sides
    uint64_t items_merged_duplicate;            // folded into an item already placed at the same breakpoint
    uint64_t items_accepted;
    uint64_t items_rejected_frequency;
    uint64_t items_rejected_skew;
    uint64_t items_ignored_circular;

    DiscordantPairSummary()
    : read_length_avg(0.0)
    , pair_distance_median(0.0)
    , pair_distance_mad(0.0)
    , pair_distance_cutoff(0.0)
    , crossing_reference_coverage(0.0)
    , expected_concordant_crossing(0.0)
    , crossing_use_empirical(true)
    , crossing_normal_positions(0)
    , frequency_cutoff(0.0)
    , skew_cutoff(0.0)
    , minimum_crossing(0.0)
    , skew_in_force(false)
    , background_e_value_cutoff(0.0)
    , minimum_pairs_option(0)
    , minimum_pairs_used(0)
    , candidate_junctions(0)
    , background_mean(0.0)
    , background_size(0.0)
    , items_tested(0)
    , items_dropped_unsupported(0)
    , items_merged_duplicate(0)
    , items_accepted(0)
    , items_rejected_frequency(0)
    , items_rejected_skew(0)
    , items_ignored_circular(0)
    {}
  };

  // Everything that decided which pair-distance (PD) calls were kept. PD's criterion is a genome-wide
  // e-value calibrated FROM THIS RUN rather than a fixed cutoff, so none of it is recoverable from the
  // command line or from the evidence .gd -- it has to be recorded here to be explicable at all.
  class PairDistanceSummary : public JSONStorable<PairDistanceSummary>
  {
  public:
    // Library shape. mean_covering_gap is the length-biased mean gap of the pairs that can cover a
    // point, E[(d-trunc)^2]/E[(d-trunc)], which is the correlation length of the seed statistic and
    // therefore what sets the effective number of independent tests. It is NOT median - 2*read_length,
    // which goes negative whenever fragments are shorter than two reads and silently disabled the
    // genome-wide correction entirely.
    double   read_length_avg;
    double   pair_distance_median;
    string   pair_orientation;                  // majority orientation; pairs of any other are ignored
    double   mean_covering_gap;
    double   n_effective_tests;                 // reference length / mean_covering_gap

    // Seed stage: the loose sensitivity filter, and the empirical excursion multiplicity measured
    // over the whole reference at a ladder of thresholds.
    double   seed_z_used;
    double   z_iqr_inflation;                   // diagnostic only: the old bulk (IQR) spread estimate
    bool     region_tail_fit_ok;
    double   region_tail_log_intercept;         // log R(t) = a - t^2 / (2 sigma^2)
    double   region_tail_sigma;
    double   region_tail_fit_lo;                // fitting window in t ...
    double   region_tail_fit_hi;                // ... over which real events are negligible
    uint64_t regions_seeded;

    // How many pairs the intake examined, and how many it refused because their placement was
    // CHOSEN rather than measured (kBreseqAmbiguousPairPlacementBAMTag). Reported because the whole
    // claim that refusing them is free rests on this ratio being small on ordinary data -- so it
    // should be visible on every run rather than inferred.
    uint64_t pairs_considered;
    uint64_t pairs_ambiguous_placement;

    // Accept stage.
    double   score_cutoff;                      // --pair-distance-score-cutoff in force

    // Outcome.
    uint64_t items_tested;
    uint64_t items_dropped_score;               // score < 0: expected by chance genome-wide, never emitted
    uint64_t items_accepted;
    uint64_t items_rejected_score;
    uint64_t items_rejected_other;              // count, frequency, duplicates, size, suppression, ...

    PairDistanceSummary()
    : read_length_avg(0.0)
    , pair_distance_median(0.0)
    , mean_covering_gap(0.0)
    , n_effective_tests(0.0)
    , seed_z_used(0.0)
    , z_iqr_inflation(1.0)
    , region_tail_fit_ok(false)
    , region_tail_log_intercept(0.0)
    , region_tail_sigma(1.0)
    , region_tail_fit_lo(0.0)
    , region_tail_fit_hi(0.0)
    , regions_seeded(0)
    , pairs_considered(0)
    , pairs_ambiguous_placement(0)
    , score_cutoff(0.0)
    , items_tested(0)
    , items_dropped_score(0)
    , items_accepted(0)
    , items_rejected_score(0)
    , items_rejected_other(0)
    {}

    //! Expected number of independent noise regions genome-wide at seed threshold t.
    double expected_regions(double t) const
    {
      if (!region_tail_fit_ok || !(region_tail_sigma > 0.0)) return 0.0;
      return exp(region_tail_log_intercept - (t * t) / (2.0 * region_tail_sigma * region_tail_sigma));
    }
  };

  // Everything that decided which missing-pair (MP) calls were kept. Like PD's, MP's criterion is a
  // genome-wide e-value calibrated FROM THIS RUN, so none of it can be read off the command line or
  // the evidence .gd -- without this record a user cannot tell why a run predicted few or many.
  //
  // The null it records is the rate at which a mapped read's mate produced no alignment at all.
  // Nothing about that rate is a constant of nature: it depends on the library, the read length and
  // the aligner, and it is spread very unevenly along a genome. Measuring it per run, rather than
  // assuming it is negligible, is the difference between MP evidence and a list of the noisiest
  // windows in the sample.
  class MissingPairSummary : public JSONStorable<MissingPairSummary>
  {
  public:
    // Library / window geometry. window_width is the width BOTH the seed statistic and the score's
    // denominator are counted over; the null below is only meaningful at that width.
    double   pair_distance_median;
    double   window_width;
    double   n_effective_tests;                 // 2 * reference length / window_width

    // The null, measured over every reference column during the stage-08 pileup.
    double   null_rate;                         // p0 in force (after the minimum-rate floor)
    double   null_rate_raw;                     // before the floor
    double   dispersion;                        // rho in force; 0 => plain binomial
    double   dispersion_raw;                    // before clamping
    double   pearson_phi;                       // chi2 / (N-1) over the fitted columns
    uint64_t tested_columns;                    // columns entering the dispersion fit
    uint64_t trimmed_columns;                   // columns excluded by the trim frequency
    double   mean_tested_reads;                 // mean n over the fitted columns

    // Seed: a sensitivity filter only. What decides a call is score_cutoff.
    double   seed_fraction_used;
    uint64_t regions_seeded;

    // Accept stage.
    double   score_cutoff;                      // --missing-pair-score-cutoff in force

    // Outcome.
    uint64_t items_tested;
    uint64_t items_dropped_score;               // score < 0: expected by chance genome-wide, never emitted
    uint64_t items_dropped_unplaced;            // no supporting read survived the rescan window
    uint64_t items_accepted;
    uint64_t items_rejected_score;
    uint64_t items_rejected_other;              // count, frequency, duplicates, suppression, ...

    MissingPairSummary()
    : pair_distance_median(0.0)
    , window_width(0.0)
    , n_effective_tests(0.0)
    , null_rate(0.0)
    , null_rate_raw(0.0)
    , dispersion(0.0)
    , dispersion_raw(0.0)
    , pearson_phi(0.0)
    , tested_columns(0)
    , trimmed_columns(0)
    , mean_tested_reads(0.0)
    , seed_fraction_used(0.0)
    , regions_seeded(0)
    , score_cutoff(0.0)
    , items_tested(0)
    , items_dropped_score(0)
    , items_dropped_unplaced(0)
    , items_accepted(0)
    , items_rejected_score(0)
    , items_rejected_other(0)
    {}
  };

  // How CNery's copy number analysis of one reference sequence went: the replication bias it fit and
  // divided out, and how much flatter each successive correction actually made the coverage.
  //
  // Kept because none of it is recoverable afterwards -- CNery writes it into
  // 09_copy_number_variation/, which the pipeline deletes when it finishes, and the CN evidence
  // entries record only the calls, not the fit that produced them.
  class CopyNumberSummary : public JSONStorable<CopyNumberSummary>
  {
  public:
    // The ori-ter replication bias. When otr_detected is false CNery found none, and the
    // coordinates below are the placeholders it writes in that case (the first and last position of
    // the sequence) rather than a real origin and terminus.
    bool     otr_detected;
    int32_t  origin;              // 1-based genomic coordinate (CNery's "Origin window" is not an index)
    int32_t  terminus;
    double   origin_coverage;     // the fitted ramp's height at the origin ...
    double   terminus_coverage;   // ... and at the terminus
    double   otr_ratio;           // origin_coverage / terminus_coverage: the peak-to-trough ratio

    // Reported whether or not a ramp was fit, because these are what say WHY one was not.
    //
    // relative_copy_number is this sequence's coverage against the LONGEST sequence of the run,
    // which reads exactly 1.0 -- the number that says a plasmid sits at 2.96x. Non-integral on
    // purpose: it is a measurement, and rounding it would discard what makes it worth reporting.
    // 0.0 => CNery did not report one.
    double   relative_copy_number;
    // How the coordinates above were arrived at, and which arm supplied them ("coverage fit",
    // "GC skew", "not corrected"). Empty => this CNery did not report it.
    string   correction_type;
    string   breakpoint_source;
    // Non-empty => the sequence had nothing to measure at all, and says which way. Distinct from
    // "no ori-ter bias detected": one means the fit was rejected, the other that no reads landed.
    string   no_coverage_reason;

    // Resolution the analysis was run at.
    int32_t  window_size;
    int32_t  window_step;

    // How large the GC correction was: the range of the LOWESS factor divided out. A range close to
    // 1.0-1.0 means there was hardly any GC bias to remove.
    double   gc_correction_min;
    double   gc_correction_max;

    // Spread of the same coverage measurement after each successive correction, as a robust
    // coefficient of variation (1.4826 * MAD / median). Measured over single-copy windows only, so
    // that a real amplification or deletion -- which inflates all three equally -- cannot mask the
    // improvement the corrections made.
    uint64_t spread_windows;      // how many windows the three values below were computed over
    bool     spread_single_copy;  // false => no HMM calls to select on, so all windows were used
    double   cv_uncorrected;
    double   cv_gc_corrected;
    double   cv_otr_gc_corrected;

    CopyNumberSummary()
    : otr_detected(false)
    , origin(0)
    , terminus(0)
    , origin_coverage(0.0)
    , terminus_coverage(0.0)
    , otr_ratio(0.0)
    , relative_copy_number(0.0)
    , window_size(0)
    , window_step(0)
    , gc_correction_min(0.0)
    , gc_correction_max(0.0)
    , spread_windows(0)
    , spread_single_copy(false)
    , cv_uncorrected(0.0)
    , cv_gc_corrected(0.0)
    , cv_otr_gc_corrected(0.0)
    {}
  };

  class CopyNumberSummaries : public map<string, CopyNumberSummary>, public JSONStorable<CopyNumberSummaries>
  {
  public:
    CopyNumberSummaries() {}
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
    SoftClippingSummary soft_clipping;
    DiscordantPairSummary discordant_pair;
    PairDistanceSummary pair_distance;
    MissingPairSummary missing_pair;
    // Per reference sequence, and only under --predict-copy-number.
    CopyNumberSummaries copy_number;
    // Stage-03 preprocessing fit (median/upper-MAD/cutoff/orientation + preprocessing-tabulated counts).
    PairedMappingDistanceDistributionSummaries preliminary_paired_mapping_distance_distribution;
    // Final pass: same fit fields, plus mapped/concordant counts from BAM pair-flag assignment.
    PairedMappingDistanceDistributionSummaries paired_mapping_distance_distribution;

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

  // PairedMappingDistanceDistributionSummary
  void to_json(json& j, const PairedMappingDistanceDistributionSummary& s);
  void from_json(const json& j, PairedMappingDistanceDistributionSummary& s);

  // PairedMappingDistanceDistributionSummaries
  void to_json(json& j, const PairedMappingDistanceDistributionSummaries& s);
  void from_json(const json& j, PairedMappingDistanceDistributionSummaries& s);

  // SoftClippingSummary
  void to_json(json& j, const SoftClippingSummary& s);
  void from_json(const json& j, SoftClippingSummary& s);

  // DiscordantPairSummary
  void to_json(json& j, const DiscordantPairSummary& s);
  void from_json(const json& j, DiscordantPairSummary& s);

  // PairDistanceSummary
  void to_json(json& j, const PairDistanceSummary& s);
  void from_json(const json& j, PairDistanceSummary& s);

  // MissingPairSummary
  void to_json(json& j, const MissingPairSummary& s);
  void from_json(const json& j, MissingPairSummary& s);

  // CopyNumberSummary
  void to_json(json& j, const CopyNumberSummary& s);
  void from_json(const json& j, CopyNumberSummary& s);

  // CopyNumberSummaries
  void to_json(json& j, const CopyNumberSummaries& s);
  void from_json(const json& j, CopyNumberSummaries& s);

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
    bool reads_were_split;
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
    double coverage_nbinom_relative_variance;
    double coverage_average;
    double coverage_variance;
    double coverage_relative_variance;
    
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
    string sample_name;
    string genbank_field_for_seq_id;
    int32_t num_processors;
    bool skip_read_filtering;
    bool skip_new_junction_prediction;
    bool skip_read_alignment_and_missing_coverage_prediction;
    bool skip_missing_coverage_prediction;
    bool no_evidence_html;

    //! Settings: Read File
    bool aligned_sam_mode;
    double  read_file_coverage_fold_limit;
    uint32_t read_file_read_length_min;
    double read_file_max_same_base_fraction;
    double read_file_max_N_fraction;
    uint32_t read_file_long_read_trigger_length;
    uint32_t read_file_long_read_split_length;
    bool read_file_long_read_distribute_remainder;

    
    //! Settings: Read Alignment
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
    int32_t  junction_indel_split_length;
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
    bool junction_allow_suboptimal_matches;
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
    double polymorphism_fisher_strand_p_value_cutoff;
    double polymorphism_ks_quality_p_value_cutoff;
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
    int32_t ignore_within_this_multiple_of_average_read_length_of_contig_end;
    
    //! Settings: Output
    uint32_t max_flanking_bases;
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
    SoftClippingSummary soft_clipping;
    DiscordantPairSummary discordant_pair;
    PairDistanceSummary pair_distance;
    MissingPairSummary missing_pair;
    CopyNumberSummaries copy_number;
    PairedMappingDistanceDistributionSummaries preliminary_paired_mapping_distance_distribution;
    PairedMappingDistanceDistributionSummaries paired_mapping_distance_distribution;

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
