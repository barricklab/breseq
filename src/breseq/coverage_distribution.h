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

// @DTF: VERSION: BINARY EDGE SEARCH 2

#ifndef _BRESEQ_COVERAGE_DISTRIBUTION_H_
#define _BRESEQ_COVERAGE_DISTRIBUTION_H_

#include "common.h"

#include "reference_sequence.h"
#include "settings.h"
#include "pileup_base.h"

using namespace std;

namespace breseq {

  // A unique-only coverage histogram: n[k] is the number of reference positions with unique-only
  // coverage depth k, for k=1..max_coverage. n[0] is unused (the *.unique_only_coverage_distribution.tab
  // file has no row for coverage=0; see error_count_pileup::print_coverage).
  struct coverage_histogram_t {
    vector<double> n;
    uint32_t max_coverage;
  };

  // Read a two-column (coverage, n) *.unique_only_coverage_distribution.tab file into a histogram.
  coverage_histogram_t read_coverage_histogram(const string& distribution_file_name);

  // Decide whether the JC pos_hash "skew" p-value should marginalize coverage over the EMPIRICAL
  // coverage histogram (returns true) or the fitted negative binomial (returns false). Empirical is
  // used only when the histogram has enough non-deletion positions that its resolution ceiling
  // (~log10(N)) clears the reject cutoff with a fixed 1-decade margin:
  //   deletion_floor = max(1, round(0.05*average));  N = sum of n[c] for c > deletion_floor;
  //   empirical iff N > 0 and log10(N) >= neg_log10_p_value_cutoff + 1.0.
  // Sets out_N and out_deletion_floor. (For small references N is tiny, so the empirical p-value would
  // be floored near 1/N and could never fail skew -- hence the fall back to the parametric nbinom.)
  bool use_empirical_pos_hash_coverage(const coverage_histogram_t& hist, double average,
                                       double neg_log10_p_value_cutoff,
                                       double& out_N, int32_t& out_deletion_floor);

  // Fit a negative binomial (size, mu) to a count histogram (n[i] = # observations with value i, i in
  // [0,maxidx]) via the censored least-squares fitter. Returns true and sets out_size/out_mu on a valid
  // fit (both > 0). Public wrapper so e.g. dp_evidence can build a parametric fallback for the
  // concordant-pair crossing distribution when the empirical one is too sparse.
  bool fit_negative_binomial_histogram(const vector<double>& n, uint32_t maxidx,
                                       double& out_size, double& out_mu);

  // Result of fitting a negative binomial distribution to a unique-only
  // coverage histogram (see CoverageDistribution::fit).
  struct CoverageDistributionFitResult
  {
    double nb_fit_size = 0;
    double nb_fit_mu = 0;
    double average = 0;
    double variance = 0;
    double relative_variance = 0;
    double deletion_coverage_propagation_cutoff = -1;
  };

	class CoverageDistribution
	{
	public:

		string path;

		static CoverageDistributionFitResult fit(
                              string distribution_file_name,
                              string plot_file,
                              double deletion_propagation_pr_cutoff
                              );
    
		static void analyze_unique_coverage_distribution(
                                                     Settings& settings, 
                                                     Summary& summary, 
                                                     cReferenceSequences& ref_seq_info, 
                                                     uint32_t coverage_group_id, 
                                                     string plot_file_name, 
                                                     string distribution_file_name,
                                                     string step_key
                                                     );
    
		static void analyze_unique_coverage_distributions(
                                                      Settings& settings, 
                                                      Summary& summary, 
                                                      cReferenceSequences& ref_seq_info, 
                                                      string plot_key, 
                                                      string distribution_file_name,
                                                      string step_key
                                                      );
    
    // Output GC content of all reads (of a theorized length)
    static void analyze_coverage_bias (
                                       const string& _fasta_file_name,
                                       const string& _bam_file_name,
                                       const string& _output_file_name,
                                       const int32_t _read_length
                                      );

	}; // class CoverageDistribution

  // HEAVY-TAILED LIBRARIES (see PairedMappingDistanceDistribution::fit).
  //
  // The distance cutoff takes the spread of the upper half of the distribution to describe its
  // whole tail, which is true of most libraries and badly false of some. When more than
  // kPairDistanceMaximumTailFraction of pairs lie beyond the cutoff, it is moved out to the
  // distance that leaves that fraction beyond it -- but only on evidence that those pairs are the
  // library's own long fragments and not pairs spanning a rearrangement, which the histogram of
  // distances cannot tell apart. Where they are tells them apart: a library's tail is spread
  // evenly over the genome, and a rearrangement's pairs are at its breakpoints. So the pairs are
  // also tallied by REGION of the reference (by the leftmost coordinate of the pair), regions
  // with an excess of far pairs are set aside, and the tail is measured in the rest. A reference
  // too small or too thinly covered to give kPairDistanceMinimumRegions usable regions keeps the
  // rule's cutoff: there, one deletion in much of the sample can be most of the far pairs.
  const double   kPairDistanceMaximumTailFraction = 0.001;
  const uint32_t kPairDistanceRegionsPerReference = 100;   // regions the reference is divided into...
  const uint32_t kPairDistanceMinimumRegionSize = 1000;    // ...but none smaller than this
  const uint32_t kPairDistanceMinimumRegionPairs = 200;    // pairs a region needs to be used
  const uint32_t kPairDistanceMinimumRegions = 20;         // usable regions needed to move the cutoff
  const double   kPairDistanceChimericMultiple = 10;       // pairs beyond this many times the rule's cutoff are chimeric, not tail
  //! Size of the regions pairs are tallied by, for a reference of this total length.
  inline uint32_t pair_distance_region_size(const uint64_t total_reference_length)
  { return static_cast<uint32_t>(max<uint64_t>(kPairDistanceMinimumRegionSize, total_reference_length / kPairDistanceRegionsPerReference)); }
  //! Name of the by-region tally written next to a pair distance distribution file.
  inline string pair_distance_region_file_name(const string& distribution_file_name)
  { return distribution_file_name + ".regions.csv"; }

  // Result of computing a robust outlier cutoff for a paired-mapping distance histogram (see
  // PairedMappingDistanceDistribution::fit) via the median absolute deviation (MAD) and the
  // Iglewicz & Hoaglin modified z-score rule.
  struct PairedMappingDistanceDistributionFitResult
  {
    double median = 0;
    double mad = 0;             // upper (one-sided) median absolute deviation
    double distance_cutoff = 0; // median + z*mad/0.6745
    string majority_orientation; // most common of FF/FR/RF by total observation count
    double mapped_pairs = 0;     // total same-reference mapped pairs (all orientations) in the preprocessing tabulation
    double concordant_pairs = 0; // majority-orientation pairs within the cutoff
  };

  class PairedMappingDistanceDistribution
  {
  public:

    // Reads the condensed (orientation,distance,count) CSV written by
    // PreprocessAlignments::merge_two_sets_of_paired_sam_files/preprocess_one_set_of_paired_sam_files,
    // determines the majority orientation by total observation count, computes the median,
    // MAD, and a modified-z-score outlier cutoff for the majority-orientation distance
    // histogram (ignoring rows for other orientations), and draws a diagnostic plot.
    // min_nonoverlapping_distance is the smallest distance a pair of non-overlapping reads can
    // span (2 x the read length); it is shaded grey on the plot. Pass 0 to omit the band.
    static PairedMappingDistanceDistributionFitResult fit(
                                                          string distribution_file_name,
                                                          string plot_file,
                                                          string histogram_file_name = "",
                                                          uint32_t min_nonoverlapping_distance = 0
                                                          );

    // Fits and plots the paired-mapping distance distribution for one paired read file set,
    // storing the result into summary.paired_mapping_distance_distribution[read_file_set.m_base_name].
    // The CSV/plot are not registered as deletable intermediates -- like the final coverage
    // distribution plot, they must survive to the end of the run for the HTML report.
    static void fit_paired_mapping_distance_distribution(
                                                         Settings& settings,
                                                         Summary& summary,
                                                         const cReadFileSet& read_file_set
                                                         );

    // Entry point called from the pipeline. Loops over settings.read_file_sets, calling
    // fit_paired_mapping_distance_distribution for every paired set (unpaired sets have no CSV).
    static void fit_paired_mapping_distance_distributions(
                                                          Settings& settings,
                                                          Summary& summary
                                                          );

  }; // class PairedMappingDistanceDistribution

} // namespace breseq

#endif
