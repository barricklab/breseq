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

		static vector<string> fit(Settings& settings, string distribution_file_name, string plot_file, double deletion_propagation_pr_cutoff, double junction_coverage_pr_cutoff, double junction_accept_pr_cutoff, double junction_keep_pr_cutoff, double junction_max_score);
		static void analyze_unique_coverage_distribution(Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, string seq_id, string plot_file_name, string distribution_file_name);
		static void analyze_unique_coverage_distributions(Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, string plot_key, string distribution_file_name);

	}; // class CoverageDistribution

} // namespace breseq

#endif
