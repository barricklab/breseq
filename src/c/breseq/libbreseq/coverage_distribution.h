/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2012 The University of Texas at Austin

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

using namespace std;

namespace breseq {

	class CoverageDistribution
	{
	public:

		string path;
		string r_script;

		static vector<string> fit(
                              Settings& settings, 
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
    
    //Tiling is taking a coverage file and grouping entries
    //together in sets of a specified constant in order to produce a condensed
    //file. A constant of 500 will take all entries with positions:
    //1-500, 501-1000, 1001-1500...
    //and average the values in each section to produce a file like so:
    //0       0
    //500     a
    //1000    b
    //1500    c
    //Where a, b, and c are the average of all values in the second column
    //in their corresponding section.
    //Written by Aaron Reba
    static void tile (
                      bool ignore_redundant_coverage,
                      string in_file_name,
                      string out_file_name,
                      string tiled_for_edging_file_name,
                      uint32_t tile_size
                      );

    //Given a coverage file, this will perform the circular binary
    //segmentation algorithm on the file. This separates the positions of the
    //file into distinct ranges.
    //The produced file might look like so:
    //1 10
    //12 30
    //11 11
    //This function also accepts files created by coverage_distribution::tile()
    //Written by Aaron Reba
    static void find_segments (const Settings& settings,
                               double summary_average,
                               string in_file_name,
                               string tiled_for_edging_file_name,
                               string out_file_name,
                               string history_file_name
                               );

    //Given a coverage file and a range file produced from
    //coverage_distribution::find_segments, this will smooth the coverage inside
    //the ranges to the mean of the range rounded to the nearest integer.
    //Written by Aaron Reba
    static void smooth_segments (const Settings& settings,
                                const string& seq_id,
                                double summary_average,
                                string tile_file_name,
                                string segment_file_name,
                                string out_file_name,
                                string final_file_name,
                                string gd_file_name
                                );
    
    //Given a coverage file, this will find the sum of all differences squared
    //of the specified coverage values element-by-element.
    
    //If the method is "ut,ub" the data taken will be the unique_top coverage
    //and the unique_bot coverage.
    //If the method is "ut+rt,ut+rb" the data taken will be the unique_top
    //coverage + the redundant_top coverage and the unique_top coverage +
    //the redundant_bot coverage.
    
    //The offset will shift the second set of data (either the unique_bot or
    //the unique_top coverage + the redundant_bot coverage), to the left and
    //append whatever is misaligned to the end of the data.
    
    //Written by Aaron Reba
    static void calculate_periodicity (
                                      string coverage_file_name,
                                      string period_file_name,
                                      uint32_t method,
                                      uint32_t offset_start,
                                      uint32_t offset_end,
                                      uint32_t offset_step
                                      );

	}; // class CoverageDistribution

  

struct search_pair_t {
  
  search_pair_t(): first(0), second(0), t_score_exists(-1.0), p_value_exists(-1.0)
  {}
  
  int32_t first;
  int32_t second;
  double t_score_exists;
  double p_value_exists;
  
};
  
  
/*
 * Data structure to store a range and its associated copy number,
 * for range merging within ::smooth_segments function.
 * Written by Tyler Fields.
 */
struct range_t {
  
  range_t(): copy_number(-1), copy_number_float(-1)
           , t_score_within(-1), p_value_within(-1)
           , t_score_exists(-1), p_value_exists(-1)
    {}
  
  // A range pair is a terminal "range", that is, a contiguous region in a
  // genome with the same copy number. Each range is thus either
  // an outlier, or the remaining area between or beside outliers. So each range
  // demarcates where an outlier either starts or ends.
  // The unit is position (in base pairs), not tile
  pair<int32_t, int32_t> range_pair;
  
  // Copy number
  double copy_number;
  double copy_number_float;

  // values for whether there was something within this
  double t_score_within;
  double p_value_within;

  // values for whether this tile "exists" (was different from those outside it)
  double t_score_exists;
  double p_value_exists;

  // Overloading "<" operator so the c++ standard sort() fn will work properly  
  bool operator < (const range_t& r) const
  {
    return range_pair.first < r.range_pair.first;
  }
    
};
  
} // namespace breseq

#endif
