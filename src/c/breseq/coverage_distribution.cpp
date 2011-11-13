
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

#include "libbreseq/coverage_distribution.h"

using namespace std;

namespace breseq {
  /*! fit
    @abstract Function for the interaction with R using the script coverage_distribution.r
    @param settings Used for file_paths
    @param distribution_file_name Input file for R, created by error_count() and saved as *.unique_only_coverage_distribution.tab
    @param plot_file Output by R
    @return vector<string> Each line contains a parameter set by R
    !*/

	vector<string> CoverageDistribution::fit(
                                           Settings& settings, 
                                           string distribution_file_name, 
                                           string plot_file, 
                                           double deletion_propagation_pr_cutoff, 
                                           double junction_coverage_pr_cutoff, 
                                           double junction_accept_pr_cutoff, 
                                           double junction_keep_pr_cutoff, 
                                           double junction_max_score
                                           )
	{
    pid_t pid = getpid();
		string log_file_name = distribution_file_name + ".r.log";
		string command = "R --vanilla < " + settings.program_data_path + "/coverage_distribution.r" + " > " + log_file_name;
		command += " distribution_file=" + distribution_file_name;
		command += " plot_file=" + plot_file;
		command += " deletion_propagation_pr_cutoff=" + to_string<double>(deletion_propagation_pr_cutoff);
		command += " junction_coverage_pr_cutoff=" + to_string<double>(junction_coverage_pr_cutoff);
		command += " junction_accept_pr_cutoff=" + to_string<double>(junction_accept_pr_cutoff);
		command += " junction_keep_pr_cutoff=" + to_string<double>(junction_keep_pr_cutoff);
		command += " junction_max_score=" + to_string<double>(junction_max_score);
    
		SYSTEM(command);

		ifstream ROUT(log_file_name.c_str());
		string line;
		vector<string> lines;
		while (getline(ROUT, line))
    {
      size_t pos = line.find("[1]");
			if (pos == 0)
      {
				lines.push_back(line.substr(pos+3));
      }
    }
    ROUT.close();
//		remove(log_file_name.c_str());

		return(lines);
	}

  // helper functions
  /*! analyze_unique_coverage_distribution
    @abstract Assigns variables to be sent off to the R script coverage_distribution.r


    @param settings
    @param summary
    @param ref_seq_info
    !*/

	void CoverageDistribution::analyze_unique_coverage_distribution(
                                                                  Settings& settings, 
                                                                  Summary& summary, 
                                                                  cReferenceSequences& ref_seq_info, 
                                                                  string seq_id, 
                                                                  string plot_key, 
                                                                  string distribution_file_name
                                                                  )
	{
		//initialize summary information
		summary.unique_coverage[seq_id].nbinom_size_parameter = NAN;
		summary.unique_coverage[seq_id].nbinom_mean_parameter = NAN;
		summary.unique_coverage[seq_id].nbinom_prob_parameter = NAN;
		summary.unique_coverage[seq_id].average = 1.0;
		summary.unique_coverage[seq_id].variance = NAN;
		summary.unique_coverage[seq_id].dispersion = NAN;
    summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff = 5.0;
    summary.unique_coverage[seq_id].deletion_coverage_seed_cutoff = 0;

		string unique_only_coverage_plot_file_name = settings.file_name(plot_key, "@", seq_id);
		string unique_only_coverage_distribution_file_name = settings.file_name(distribution_file_name, "@", seq_id);

		// Define various coverage thresholds...
		uint32_t sequence_length = ref_seq_info[ref_seq_info.seq_id_to_index(seq_id)].m_length;

		/// DELETION PROPAGATION CUTOFF
		// One-tailed test p=0.05, Bonferroni correction
		//# my del_propagation_pr_cutoff = 0.05 / sequence_length;

		// One-tailed test p=0.01, no Bonferroni correction
		//#my del_propagation_pr_cutoff = 0.01;

		// We really want somewhere between these two, try this...
		double deletion_propagation_pr_cutoff = 0.05 / sqrt(sequence_length);

		/// NEW JUNCTION COVERAGE CUTOFFS
		// Arbitrary value that seems to work....
    //double junction_coverage_pr_cutoff =  sqrt(settings.junction_accept_pr / static_cast<double>(sequence_length));
    double junction_coverage_pr_cutoff = 0.01;
    
		// We really want somewhere between these two, try this...
    double junction_accept_pr_cutoff = 0.01;
		double junction_keep_pr_cutoff = 0.01 / sqrt(sequence_length);
		int32_t junction_max_score = int(2 * summary.sequence_conversion.avg_read_length);
    
		CoverageDistribution dist;
		vector<string> lines = dist.fit(settings, 
                                    unique_only_coverage_distribution_file_name, 
                                    unique_only_coverage_plot_file_name,
                                    deletion_propagation_pr_cutoff, 
                                    junction_coverage_pr_cutoff, 
                                    junction_accept_pr_cutoff, 
                                    junction_keep_pr_cutoff, 
                                    junction_max_score
                                    );

		// First two lines are negative binomial parameters.
		// Next three lines are average, standard deviation, and index of overdispersion

		// Put these into summary
		summary.unique_coverage[seq_id].nbinom_size_parameter = from_string<double>(lines[0]);
		summary.unique_coverage[seq_id].nbinom_mean_parameter = from_string<double>(lines[1]);
		// Calculated by formula, prob = size/(size + mu)
		summary.unique_coverage[seq_id].nbinom_prob_parameter = summary.unique_coverage[seq_id].nbinom_size_parameter / (summary.unique_coverage[seq_id].nbinom_mean_parameter + summary.unique_coverage[seq_id].nbinom_size_parameter);
		summary.unique_coverage[seq_id].average = from_string<double>(lines[2]);
		summary.unique_coverage[seq_id].variance = from_string<double>(lines[3]);
		summary.unique_coverage[seq_id].dispersion = from_string<double>(lines[4]);

    summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff = from_string<double>(lines[5]);
    summary.unique_coverage[seq_id].junction_coverage_cutoff = from_string<double>(lines[6]);
		
    // deprecated statistics
    //summary.unique_coverage[seq_id].junction_accept_score_cutoff = from_string<double>(lines[7]);
		//summary.unique_coverage[seq_id].junction_keep_score_cutoff = from_string<double>(lines[8]);
        
    bool verbose = false;
    if (verbose)
    {
      cout << seq_id << endl;
      cout << "nbinom_size_parameter " << summary.unique_coverage[seq_id].nbinom_size_parameter << endl;
      cout << "nbinom_mean_parameter " << summary.unique_coverage[seq_id].nbinom_mean_parameter << endl;
      cout << "nbinom_prob_parameter " << summary.unique_coverage[seq_id].nbinom_prob_parameter << endl;
      cout << "average " << summary.unique_coverage[seq_id].average << endl;
      cout << "variance " << summary.unique_coverage[seq_id].variance << endl;
      cout << "dispersion " << summary.unique_coverage[seq_id].dispersion << endl;
      cout << "deletion_coverage_propagation_cutoff " << summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff << endl;
      cout << "junction_coverage_cutoff " << summary.unique_coverage[seq_id].junction_coverage_cutoff << endl;
      //cout << "junction_accept_score_cutoff " << summary.unique_coverage[seq_id].junction_accept_score_cutoff << endl;
      //cout << "junction_keep_score_cutoff " << summary.unique_coverage[seq_id].junction_keep_score_cutoff << endl;
      //cout << "pr_no_coverage_position_strand " << summary.unique_coverage[seq_id].pr_no_coverage_position_strand << endl;

    }
	}

	void CoverageDistribution::analyze_unique_coverage_distributions(
                                                                   Settings& settings, 
                                                                   Summary& summary, 
                                                                   cReferenceSequences& ref_seq_info, 
                                                                   string plot_file_name, 
                                                                   string distribution_file_name
                                                                   )
  {
		for (uint32_t i = 0; i < ref_seq_info.size(); i++) {
          
			analyze_unique_coverage_distribution(
                                           settings, 
                                           summary, 
                                           ref_seq_info, 
                                           ref_seq_info[i].m_seq_id, 
                                           plot_file_name, 
                                           distribution_file_name
                                           );
    }
	}

} // namespace breseq
