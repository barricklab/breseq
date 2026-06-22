/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2022 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the
  terms the GNU General Public License as published by the Free Software
  Foundation; either version 1, or (at your option) any later version.

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

} // namespace breseq

#endif
