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

  // Result of computing a robust outlier cutoff for a paired-mapping distance histogram (see
  // PairedMappingDistanceDistribution::fit) via the median absolute deviation (MAD) and the
  // Iglewicz & Hoaglin modified z-score rule.
  struct PairedMappingDistanceDistributionFitResult
  {
    double median = 0;
    double mad = 0;             // median absolute deviation
    double distance_cutoff = 0; // median + 3.5*mad/0.6745
  };

  class PairedMappingDistanceDistribution
  {
  public:

    // Reads the condensed (orientation,distance,count) CSV written by
    // PreprocessAlignments::merge_two_sets_of_paired_sam_files/preprocess_one_set_of_paired_sam_files,
    // determines the majority orientation by total observation count, computes the median,
    // MAD, and a modified-z-score outlier cutoff for the majority-orientation distance
    // histogram (ignoring rows for other orientations), and draws a diagnostic plot.
    static PairedMappingDistanceDistributionFitResult fit(
                                                          string distribution_file_name,
                                                          string plot_file
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
