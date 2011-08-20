
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

#include "breseq/coverage_distribution.h"

using namespace std;

namespace breseq {

	CoverageDistribution::CoverageDistribution()
	{
		if (path.size() == 0) path = ".";
		r_script = /*$FindBin::Bin .*/ "/../lib/perl5/Breseq/coverage_distribution.r";
	}

	vector<string> CoverageDistribution::fit(string distribution_file, string plot_file, uint32_t deletion_propagation_pr_cutoff, uint32_t junction_coverage_pr_cutoff, uint32_t junction_accept_pr_cutoff, uint32_t junction_keep_pr_cutoff, uint32_t junction_max_score)
	{
		string log_file_name = path + "/$$.r.log";
		string command = "R --vanilla < " + r_script + " > " + log_file_name;
		command += " distribution_file=" + distribution_file;
		command += " plot_file=" + plot_file;
		command += " deletion_propagation_pr_cutoff=" + deletion_propagation_pr_cutoff;
		command += " junction_coverage_pr_cutoff=" + junction_coverage_pr_cutoff;
		command += " junction_accept_pr_cutoff=" + junction_accept_pr_cutoff;
		command += " junction_keep_pr_cutoff=" + junction_keep_pr_cutoff;
		command += " junction_max_score=" + junction_max_score;

		int retval = system(command.c_str());

		ifstream ROUT(log_file_name.c_str());
		string line;
		vector<string> lines;
		while (getline(ROUT, line))
			if (line.find("[1]") == 0)
				lines.push_back(line);
		ROUT.close();
		remove(log_file_name.c_str());

		return(lines);
	}

	// helper functions

	void CoverageDistribution::analyze_unique_coverage_distribution(Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, string seq_id, string plot_key, string distribution_key)
	{
		//initialize summary information
		summary.unique_coverage[seq_id].nbinom_size_parameter = "ND";
		summary.unique_coverage[seq_id].nbinom_mean_parameter = "ND";
		summary.unique_coverage[seq_id].nbinom_prob_parameter = "ND";
		summary.unique_coverage[seq_id].average = 1;
		summary.unique_coverage[seq_id].variance = "ND";
		summary.unique_coverage[seq_id].dispersion = "ND";
		summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff = 5;

		string unique_only_coverage_plot_file_name = settings.file_name(plot_key, "@", seq_id);
		string unique_only_coverage_distribution_file_name = settings.file_name(distribution_key, "@", seq_id);

		// Define various coverage thresholds...
		uint32_t sequence_length = ref_seq_info[ref_seq_info.seq_id_to_index(seq_id)].m_length;

		/// DELETION PROPAGATION CUTOFF
		// One-tailed test p=0.05, Bonferroni correction
		//# my del_propagation_pr_cutoff = 0.05 / sequence_length;

		// One-tailed test p=0.01, no Bonferroni correction
		//#my del_propagation_pr_cutoff = 0.01;

		// We really want somewhere between these two, try this...
		float deletion_propagation_pr_cutoff = 0.05 / sqrt(sequence_length);

		/// NEW JUNCTION COVERAGE CUTOFFS
		// Arbitrary value that seems to work....
		float junction_coverage_pr_cutoff = 1/sequence_length; //# *0.05

		// We really want somewhere between these two, try this...
		float junction_accept_pr_cutoff = 0.01;
		float junction_keep_pr_cutoff = 0.01 / sqrt(sequence_length);
		int32_t junction_max_score = int(2 * summary.sequence_conversion.avg_read_length);

		CoverageDistribution dist;
		vector<string> lines = dist.fit(unique_only_coverage_distribution_file_name, unique_only_coverage_plot_file_name,
				deletion_propagation_pr_cutoff, junction_coverage_pr_cutoff, junction_accept_pr_cutoff, junction_keep_pr_cutoff, junction_max_score);

		// First two lines are negative binomial parameters.
		// Next three lines are average, standard deviation, and index of overdispersion

		// Put these into summary
		summary.unique_coverage[seq_id].nbinom_size_parameter = lines[0];
		summary.unique_coverage[seq_id].nbinom_mean_parameter = lines[1];
		// Calculated by formula, prob = size/(size + mu)
		summary.unique_coverage[seq_id].nbinom_prob_parameter = from_string<float>(lines[0]) / (from_string<float>(lines[0]) + from_string<float>(lines[1]));
		summary.unique_coverage[seq_id].average = from_string<uint32_t>(lines[2]);
		summary.unique_coverage[seq_id].variance = lines[3];
		summary.unique_coverage[seq_id].dispersion = lines[4];

		summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff = from_string<uint32_t>(lines[5]);
		summary.unique_coverage[seq_id].junction_coverage_cutoff = from_string<uint32_t>(lines[6]);
		summary.unique_coverage[seq_id].junction_accept_score_cutoff = from_string<uint32_t>(lines[7]);
		summary.unique_coverage[seq_id].junction_keep_score_cutoff = from_string<uint32_t>(lines[8]);
	}

	void CoverageDistribution::analyze_unique_coverage_distributions(Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info, string plot_key, string distribution_key)
	{
		for (uint32_t i = 0; i < ref_seq_info.seq_ids.size(); i++)
			analyze_unique_coverage_distribution(settings, summary, ref_seq_info, ref_seq_info.seq_ids[i], plot_key, distribution_key);
	}

} // namespace breseq