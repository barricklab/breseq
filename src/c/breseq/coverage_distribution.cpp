
/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2008-2010 Michigan State University
 Copyright (c) 2011-2022 The University of Texas at Austin
 
 breseq is free software; you can redistribute it and/or modify it under the
 terms the GNU General Public License as published by the Free Software
 Foundation; either version 1, or (at your option) any later version.
 
 *****************************************************************************/

// @DTF: VERSION: BINARY EDGE SEARCH 9

#include "libbreseq/coverage_distribution.h"

#include "libbreseq/genome_diff.h"

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
                                         double deletion_propagation_pr_cutoff
                                         )
{
  pid_t pid = getpid();
  string log_file_name = distribution_file_name + ".r.log";
  
  string command = "R --vanilla < " + double_quote(settings.program_data_path +
  "/coverage_distribution.r") + " > " + double_quote(log_file_name);
  command += " --args";
  command += " distribution_file=" + double_quote(distribution_file_name);
  command += " plot_file=" + double_quote(plot_file);
  command += " deletion_propagation_pr_cutoff=" + to_string<double>(deletion_propagation_pr_cutoff);
  
  SYSTEM(command, false, false, false); //NOTE: Not escaping shell characters here.
  
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
  //    remove(log_file_name.c_str());
  
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
                                                                uint32_t coverage_group_id,
                                                                string plot_key,
                                                                string distribution_file_name,
                                                                string step_key
                                                                )
{
  
  vector<string> seq_ids = settings.refseq_settings.m_seq_ids_by_coverage_group[coverage_group_id];
  
  // We need to know the total sequence length of just this group
  uint32_t sequence_length = 0;
  for (vector<string>::iterator it=seq_ids.begin(); it!=seq_ids.end(); it++) {
    string& seq_id = *it;
    sequence_length += ref_seq_info.get_sequence_length(seq_id);
  }
  
  // Perform no fitting if we are in targeted_sequencing mode.
  if (settings.targeted_sequencing) return;
  
  string unique_only_coverage_plot_file_name = settings.file_name(plot_key, "@", to_string<uint32_t>(coverage_group_id));
  string unique_only_coverage_distribution_file_name = settings.file_name(distribution_file_name, "@", to_string<uint32_t>(coverage_group_id));
  
  // Define various coverage thresholds...
  
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
  int32_t junction_max_score = int(2 * summary.sequence_conversion.read_length_avg);
  
  CoverageDistribution dist;
  vector<string> lines = dist.fit(settings,
                                  unique_only_coverage_distribution_file_name,
                                  unique_only_coverage_plot_file_name,
                                  deletion_propagation_pr_cutoff
                                  );
  settings.track_intermediate_file(step_key, unique_only_coverage_plot_file_name);
  settings.track_intermediate_file(step_key, unique_only_coverage_distribution_file_name);
  settings.track_intermediate_file(step_key, unique_only_coverage_distribution_file_name + ".r.log");
  
  // First two lines are negative binomial parameters.
  // Next three lines are average, standard deviation, and index of overdispersion
  
  // Put these into summary
  
  for (vector<string>::iterator it=seq_ids.begin(); it!=seq_ids.end(); it++) {
    string seq_id = *it;
    summary.unique_coverage[seq_id].nbinom_size_parameter = from_string<double>(lines[0]);
    summary.unique_coverage[seq_id].nbinom_mean_parameter = from_string<double>(lines[1]);
    // Calculated by formula, prob = size/(size + mu)
    
    // These remain at their defaults of zero if this is zero = fit failed
    if (summary.unique_coverage[seq_id].nbinom_mean_parameter != 0) {
      summary.unique_coverage[seq_id].nbinom_prob_parameter = summary.unique_coverage[seq_id].nbinom_size_parameter
      / (summary.unique_coverage[seq_id].nbinom_mean_parameter + summary.unique_coverage[seq_id].nbinom_size_parameter);
      // Calculated by formula variance = mu + mu ^ 2 / size
      summary.unique_coverage[seq_id].nbinom_variance = summary.unique_coverage[seq_id].nbinom_mean_parameter
      + pow(summary.unique_coverage[seq_id].nbinom_mean_parameter, 2) / summary.unique_coverage[seq_id].nbinom_size_parameter;
      // Calculated by formula relative_variance = variance / mu
      summary.unique_coverage[seq_id].nbinom_relative_variance = summary.unique_coverage[seq_id].nbinom_variance / summary.unique_coverage[seq_id].nbinom_mean_parameter;
    }
    
    summary.unique_coverage[seq_id].average = from_string<double>(lines[2]);
    summary.unique_coverage[seq_id].variance = from_string<double>(lines[3]);
    summary.unique_coverage[seq_id].relative_variance = from_string<double>(lines[4]);
    
    summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff = from_string<double>(lines[5]);
    
    bool verbose = false;
    if (verbose)
    {
      cout << seq_id << endl;
      cout << "nbinom_size_parameter " << summary.unique_coverage[seq_id].nbinom_size_parameter << endl;
      cout << "nbinom_mean_parameter " << summary.unique_coverage[seq_id].nbinom_mean_parameter << endl;
      cout << "nbinom_prob_parameter " << summary.unique_coverage[seq_id].nbinom_prob_parameter << endl;
      cout << "average " << summary.unique_coverage[seq_id].average << endl;
      cout << "variance " << summary.unique_coverage[seq_id].variance << endl;
      cout << "relative_variance " << summary.unique_coverage[seq_id].relative_variance << endl;
      cout << "deletion_coverage_propagation_cutoff " << summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff << endl;
    }
  }
}

void CoverageDistribution::analyze_unique_coverage_distributions(
                                                                 Settings& settings,
                                                                 Summary& summary,
                                                                 cReferenceSequences& ref_seq_info,
                                                                 string plot_file_name,
                                                                 string distribution_file_name,
                                                                 string step_key
                                                                 )
{
  vector<vector<string> > seq_ids_by_coverage_group = settings.seq_ids_by_coverage_group();
  
  for (uint32_t i=0; i< settings.seq_ids_by_coverage_group().size(); i++) {
    
    
    analyze_unique_coverage_distribution(
                                         settings,
                                         summary,
                                         ref_seq_info,
                                         i,
                                         plot_file_name,
                                         distribution_file_name,
                                         step_key
                                         );
    
    //Warning
  }
}




/*
 * Function: calculate_periodicity
 * --------------------------------
 * SUMMARY GOES HERE
 */
void CoverageDistribution::analyze_coverage_bias (
                                                   const string& _fasta_file_name,
                                                   const string& _bam_file_name,
                                                   const string& _output_file_name_prefix,
                                                   const int32_t _read_length
                                                   )
  
{
  
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFile(_fasta_file_name);
  
  bam_file bam;
  bam.open_read(_bam_file_name, _fasta_file_name);
  
  
  string read_output_file_name = _output_file_name_prefix + ".read.csv";
  ofstream read_out(read_output_file_name.c_str());
  
  // First line of output is read_length
  read_out << _read_length << endl;
  
  // Second line is lengths of reference sequences
  for (cReferenceSequences::iterator it=ref_seq_info.begin(); it != ref_seq_info.end(); it++) {
    if (it != ref_seq_info.begin()) read_out << ",";
    read_out << it->get_sequence_length() << endl;
  }
  
  // Third line is GC of reference sequences
  for (cReferenceSequences::iterator it=ref_seq_info.begin(); it != ref_seq_info.end(); it++) {
    if (it != ref_seq_info.begin()) read_out << ",";
    read_out << gc_percentage_string(it->m_fasta_sequence.get_sequence()) << endl;
  }
  
  alignment_list alignments;
  while (bam.read_alignments(alignments, false)) {
    double gc(0);

    for(alignment_list::iterator it=alignments.begin(); it!=alignments.end(); it++)
    {
      bam_alignment& a = *(it->get());
      
      int64_t start = (a.strand()==-1) ? a.reference_end_1() - _read_length + 1 : a.reference_start_1();
      string read_string = ref_seq_info[a.reference_target_id()].get_sequence_1_start_size(start, _read_length);
      // Don't need to reverse complement for GC content
      
      gc += gc_percentage_string(read_string);
    }

    gc /= alignments.size();
    
    read_out << gc << endl;
  }
  
  // Now the reference file
  string ref_output_file_name = _output_file_name_prefix + ".ref.csv";
  ofstream ref_out(ref_output_file_name.c_str());
  for (cReferenceSequences::iterator it=ref_seq_info.begin(); it != ref_seq_info.end(); it++) {

    for (size_t i=1; i<= it->get_sequence_length(); i++) {
      string read_string = it->get_sequence_1_start_size(i, _read_length);
      ref_out << gc_percentage_string(read_string) << endl;
    }
  }
}
  

} // namespace breseq
