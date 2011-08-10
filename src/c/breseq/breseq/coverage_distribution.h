/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the
  terms the GNU General Public License as published by the Free Software
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_COVERAGE_DISTRIBUTION_H_
#define _BRESEQ_COVERAGE_DISTRIBUTION_H_

#include "common.h"

#include "annotated_sequence.h"
#include "settings.h"

using namespace std;

namespace breseq {

	class CoverageDistribution
	{
	public:

		string path;
		string r_script;

		CoverageDistribution();
		vector<string> fit(string distribution_file, string plot_file, uint32_t deletion_propagation_pr_cutoff, uint32_t junction_coverage_pr_cutoff, uint32_t junction_accept_pr_cutoff, uint32_t junction_keep_pr_cutoff, uint32_t junction_max_score);
		static void analyze_unique_coverage_distribution(Settings& settings, string seq_id, Summary& summary, string plot_key, string distribution_key);
		static void analyze_unique_coverage_distributions(Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, string plot_key, string distribution_key);

	}; // class CoverageDistribution

} // namespace breseq

#endif