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

		static vector<string> fit(
                              Settings& settings, 
                              string distribution_file_name, 
                              string plot_file, 
                              double deletion_propagation_pr_cutoff, 
                              double junction_coverage_pr_cutoff, 
                              double junction_accept_pr_cutoff, 
                              double junction_keep_pr_cutoff, 
                              double junction_max_score
                              );
    
		static void analyze_unique_coverage_distribution(
                                                     Settings& settings, 
                                                     Summary& summary, 
                                                     cReferenceSequences& ref_seq_info, 
                                                     string seq_id, 
                                                     string plot_file_name, 
                                                     string distribution_file_name
                                                     );
    
		static void analyze_unique_coverage_distributions(
                                                      Settings& settings, 
                                                      Summary& summary, 
                                                      cReferenceSequences& ref_seq_info, 
                                                      string plot_key, 
                                                      string distribution_file_name
                                                      );
    
    //Tiling is taking a coverage file and grouping entries
    //together in sets of a specified constant in order to produce a condensed
    //file. A constant of 500 will take all entries with positions:
    //0-499, 500-999, 1000-1499...
    //and average the values in each section to produce a file like so:
    //0       a
    //500     b
    //1000    c
    //Where a, b, and c are the average of all values in the second column
    //in their corresponding section.
    //Written by Aaron Reba
    static void tile (
                      string in_file_name,
                      string out_file_name,
                      uint32_t tile_size
                      );
    
    //Given a coverage file, this will perform the circular binary
    //segmentation algorithm on the file. This separates the positions of the
    //file into distinct ranges.
    //The produced file might look like so:
    //0 10
    //12 30
    //11 11
    //This function also accepts files created by coverage_distribution::tile()
    //Written by Aaron Reba
    static void find_segments (
                               string in_file_name,
                               string out_file_name
                               );
    
    //Given a coverage file and a range file produced from
    //coverage_distribution::find_segments, this will smooth the coverage inside
    //the ranges to the mean of the range rounded to the nearest integer.
    //Written by Aaron Reba
    //static void smooth_segments (
    //                            string in_file_name,
    //                            string range_file_name,
    //                            string out_file_name
    //                            );

	}; // class CoverageDistribution

} // namespace breseq

#endif
