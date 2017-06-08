
/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2008-2010 Michigan State University
 Copyright (c) 2011-2012 The University of Texas at Austin
 
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
  
  string command = "R --vanilla < " + cString(settings.program_data_path).escape_shell_chars() +
  "/coverage_distribution.r" + " > " + cString(log_file_name).escape_shell_chars();
  command += " distribution_file=" + cString(distribution_file_name).escape_shell_chars();
  command += " plot_file=" + cString(plot_file).escape_shell_chars();
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
  
  //initialize summary information
  uint32_t sequence_length = 0;
  for (vector<string>::iterator it=seq_ids.begin(); it!=seq_ids.end(); it++) {
    string seq_id = *it;
    summary.unique_coverage[seq_id].nbinom_size_parameter = NAN;
    summary.unique_coverage[seq_id].nbinom_mean_parameter = NAN;
    summary.unique_coverage[seq_id].nbinom_prob_parameter = NAN;
    summary.unique_coverage[seq_id].average = 1.0;
    summary.unique_coverage[seq_id].variance = NAN;
    summary.unique_coverage[seq_id].dispersion = NAN;
    summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff = 0.0;
    summary.unique_coverage[seq_id].deletion_coverage_seed_cutoff = 0;
    
    // Keep running total of sequence length
    sequence_length += ref_seq_info[ref_seq_info.seq_id_to_index(seq_id)].m_length;
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
  int32_t junction_max_score = int(2 * summary.sequence_conversion.avg_read_length);
  
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
    summary.unique_coverage[seq_id].nbinom_prob_parameter = summary.unique_coverage[seq_id].nbinom_size_parameter
    / (summary.unique_coverage[seq_id].nbinom_mean_parameter + summary.unique_coverage[seq_id].nbinom_size_parameter);
    // Calculated by formula variance = mu + mu ^ 2 / size
    summary.unique_coverage[seq_id].nbinom_variance = summary.unique_coverage[seq_id].nbinom_mean_parameter
    + pow(summary.unique_coverage[seq_id].nbinom_mean_parameter, 2) / summary.unique_coverage[seq_id].nbinom_size_parameter;
    // Calculated by formula dispersion = variance / mu
    summary.unique_coverage[seq_id].nbinom_dispersion = summary.unique_coverage[seq_id].nbinom_variance / summary.unique_coverage[seq_id].nbinom_mean_parameter;
    
    summary.unique_coverage[seq_id].average = from_string<double>(lines[2]);
    summary.unique_coverage[seq_id].variance = from_string<double>(lines[3]);
    summary.unique_coverage[seq_id].dispersion = from_string<double>(lines[4]);
    
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
      cout << "dispersion " << summary.unique_coverage[seq_id].dispersion << endl;
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
  }
}


/*
 * Function: tile
 * --------------------------------
 * Reads in [genome].coverage.tab and generates [genome].tiled.tab.
 
 * Summary: Breaks the [genome].coverage.tab file into sections
 * of equal length (ex: 500 bp) called tiles, calculates the
 * average coverage within each tile, and then writes this data
 * to [genome].tiled.tab
 *
 * THE DOCUMENTATION FOR THIS FUNCTION IS A WORK IN PROGESS
 */

/* Notes & Explanations (below, collapse until needed) */
/* --------------------------------
 * Conceptual Overview: The functions within the CoverageDistribution class keep track of location within a genome or segment of a genome (hereby referred to as just "segment") in two ways: by "position" and by "tile". POSITION is referring to base pairs, if every base pair within a segment were numbered sequentially from beginning to end, starting at 1 (why 1? because that's where biologists like to start their lists). Each position has a corresponding coverage value. You can see this by looking at the coverage.tab file (in the 08_mutation_identification directory). However, for our purposes, analyzing every single position in an entire segment would take too long. So for efficiency we use this function to break the the segment into sections of equal length (e.g. 500 base pairs long) called TILES. We then calculate the average coverage within within each tile and write this data to a new file called tiled.tab (in the 09_copy_number_variation directory). This file will then be analyzed by the CoverageDistribution::find_segments function to find "outliers" (see find_segments for more about outliers).
 
 */
void CoverageDistribution::tile(
                                bool ignore_redundant_coverage,
                                string in_file_name,
                                string out_file_name,
                                string tiled_for_edging_file_name,
                                uint32_t tile_size
                                )
{
  // @DTF: make this a default setting later on to pass as parameter
  // to this function:
  int32_t tile_size_for_edging = 10;
  
  cout << "Performing tiling..." << endl;
  cout << "  1: Identification tile size    " << setw(6) << right <<  tile_size << " bp"<< endl;
  cout << "  2: Edge optimization tile size " << setw(6) << right << tile_size_for_edging << " bp"<< endl;

  
  //The index is the marker of the current tile.
  //It is always a multiple of tile_size
  uint32_t tile_index;
  uint32_t e_tile_index; // Same process repeated for edging file
  
  //these are for taking running sums and averaging.
  double tile_sum;
  double e_tile_sum; // Same process repeated for edging file
  double average;
  double e_average; // Same process repeated for edging file
  
  //position and coverage are taken from the file.
  uint32_t position;
  double coverage;
  double sum_coverage;
  
  //used for checking if there is any redundant coverage
  double redundant1;
  double redundant2;
  double redundant3;
  double redundant4;
  
  //if ignoring redundant coverage is enabled, this will disable writing the
  //tile
  bool do_write;
  
  //used to make sure the very last position is not written twice.
  uint32_t last_position;
  
  //used to ignore values in the input file.
  string skip;
  
  ifstream in_file;
  ofstream out_file;
  ofstream out_file_for_edging;
  
  in_file.open ( in_file_name.c_str() );
  out_file.open ( out_file_name.c_str() );
  out_file_for_edging.open( tiled_for_edging_file_name.c_str() );
  
  //skipping the first line in the input file.
  getline(in_file, skip);
  
  out_file << "position\tcoverage\t" << tile_size << "\n";
  out_file_for_edging << "position\tcoverage\t" << tile_size_for_edging << "\n";
  
  out_file << "0\t0\n";
  out_file_for_edging << "0\t0\n";
  
  tile_index = 0;
  e_tile_index = 0;
  
  tile_sum = 0;
  e_tile_sum = 0;
  
  do_write = true;
  
  while ( in_file >> coverage )
  {
    sum_coverage = coverage;
    
    in_file >> coverage;
    sum_coverage += coverage;
    
    //next 4 are values that represent redundancy
    in_file >> redundant1 >> redundant2 >> redundant3 >> redundant4;
    
    //the next value in the line isn't used.
    in_file >> skip;
    
    //last is position (@DTF: begins at 1, not 0)
    in_file >> position;
    
    coverage = sum_coverage;
    
    
    // Check if the current position is in the current tile.
    // If it is, add the coverage to the sum and increment the count.
    // Else, since we just read in coverage data from the first position of a new tile
    // (i.e. 1, 501, 1001) calculate the avg without this new coverage value, print the
    // avg to file, then update all variables, including clearing the sum of coverage
    // and resetting to the last coverage value that was read in.
    if ( position <= tile_index + tile_size )
    {
      tile_sum += coverage;
    }
    else
    {
      if ( do_write )
      {
        average = tile_sum / tile_size;
        out_file << tile_index + tile_size << "\t" << average << "\n";
      }
      
      //Set tile_index to position rounded down to the nearest multiple
      //of tile_size
      tile_index = position / tile_size * tile_size;
      tile_sum = coverage;
      last_position = position;
      do_write = true;
    }
    
    // Do the same as above, only for the edging file
    if ( position <= e_tile_index + tile_size_for_edging )
    {
      e_tile_sum += coverage;
    }
    else
    {
      if ( do_write )
      {
        e_average = e_tile_sum / tile_size_for_edging;
        out_file_for_edging << e_tile_index + tile_size_for_edging << "\t" << e_average << "\n";
      }
      
      //Set tile_index to position rounded down to the nearest multiple
      //of tile_size
      e_tile_index = position / tile_size_for_edging * tile_size_for_edging;
      e_tile_sum = coverage;
      last_position = position;
      do_write = true;
    }
    
    
    if ( ignore_redundant_coverage )
    {
      
      if ( redundant1 != 0 ||
          redundant2 != 0 ||
          redundant3 != 0 ||
          redundant4 != 0 )
      {
        //skip this tile.
        do_write = false;
      }
    }
    
  }
  
  
  // Now we've reached the end of the file, having just read in the last line...
  
  // When the last position in a tile is read in, the avg coverage of that tile
  // is calculated and printed to the file during the next run of the loop.
  // But if the last position in a tile is also the last position in the file, the
  // loop will break without having printed the avg coverage of the last tile.
  // So, if we've reached a point where avg coverage should be printed to file,
  // (i.e. if position is a multiple of tile size), and if it hasn't had an
  // opportunity to print yet (we know this when the above line "last_position = position"
  // hasn't had an opportunity to run), then calculate avg coverage for the
  // last tile and print to file.
  if ( position % tile_size == 0 && position != last_position )
  {
    average = tile_sum / tile_size;
    
    out_file << position << "\t" << average << "\n";
  }
  
  // Do the same as above, only for the edging file
  if ( position % tile_size_for_edging == 0 && position != last_position )
  {
    e_average = e_tile_sum / tile_size_for_edging;
    
    out_file_for_edging << position << "\t" << e_average << "\n";
  }
  
  in_file.close();
  out_file.close();
  out_file_for_edging.close();
}


/*
 * Function: find_segments
 * --------------------------------
 * Reads in [genome].tiled.tab and generates [genome].ranges.tab and
 * [genome].history.tab.
 *
 * Summary: Applies the "circular binary segmentation" algorithm (or CBS) to
 * the coverage data within [genome].tiled.tab in order to find "outliers",
 * or regions where the read coverage depth is significantly greater or
 * lesser than what is typical across the entire genome. Outputs the outliers
 * it finds to [genome].ranges.tab, and logs the progression of the algorithm
 * to [genome].history.tab for debugging purposes.
 */

/* Notes & Explanations (below, collapse until needed) */
/* --------------------------------
 * Overview, Circular Binary Segmentation: Although the theory behind the "circular binary segmentation" algorithm (or CBS) is backed by some serious math, in practice this algorithm is simple to follow and understand. Imagine the following oversimplification: a white picket fence that extends as far as you can see. Each white picket (vertical board/panel) represents one tile of the genome (see "Tiles vs Position" below). You notice the pickets are wildly uneven in height, ranging a few feet off the ground to tens of feet tall. This represents coverage depth. Although the arrangement of the pickets seems mostly random, occasionally there are contiguous regions where all the pickets are at least a certain height. When these regions are large enough compared to the rest of the fence (as measured by something called the "test statistic" or "t_score"), they're considered outliers. The CBS algorithm searches for these outliers, and when it finds one it searches through it again, treating it as its own separate genome, looking for outliers within the outlier. It continues this process, finding bumps within bumps within bumps, until no more are found.
 
 @ @   @ @
 @ @   @ @ @ @ @               @ @
 @ @   @ @ @ @ @   @ @     @ @ @ @ @ @
 @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @
 ... @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ ...
 
 * Details, Circular Binary Segmentation: The algorithm is a loop that consists of two main steps: First, find the region in the current search segment with the largest "test statistic", or "t_score". (For now, all you need to know about the t_score is that it is a single number that takes into account all the different ways you can measure the size of a potential outlier; the larger the potential outlier, the larger the t_score.) The second step is to determine if the aforementioned largest region is statisticially significant, and thus a true outlier. The first time the loop runs, the "search segment" is the entire genome, and if no outlier is found, the loop stops, having only run one time. However, if an outlier IS FOUND, the loop runs again to try to find other outliers that may have been overshadowed before. The loop keeps running until no more outliers are found. We'll now explore these steps in greater detail:
 
 1) Find the region in the current search segment with the largest t_score: To do this, we must calculate the t_score for every possible region within the search segment. Every possible region can be represented by "i" and "j", where i is the start of the region, and j is the end of the region. We use a double for-loop to iterate through every possible combination of i and j, calculating the t_score each time. NOTE: There is not just one pair of i and j variables in the code. Location comes in the form of both tile and position (see "Tiles vs Position" below), and so there are pairs of i and j variables to keep track of both.
 
 2) Determine if that largest region is statisticially significant, and thus a true outlier. This is important because random "noise" naturally occurs in the data we are processing, like static on a phone call, and we don't want to mistake static for an outlier. So step 2 asks, "Is this potential outlier significantly larger than the random static around it?" To do this, we run the following test 100 times: We clone the search segment and then randomly shuffle the tiles around, creating a new random segment. Then, just as in step 1, we find the largest t_score within this random segment. We then compare this largest random t_score to the original largest t_score, and keep track of the results. If the original largest t_score was greater 95 times out of 100, we accept it (i.e. the region it represents) as an outlier. Otherwise, we reject it.
 
 If the answer to step 2 is YES/ACCEPT, you break the search segment into three new subsegments or ranges: the middle subsegment (which is the newly found outlier), the left subsegment (everything to the left of the middle), and the right subsegment (everything to the right of the middle). You then repeat steps 1 and 2 for each of these new subsegments, treating each as its own entity with its own potential outliers to be found. We're looking for "bumbs within bumps".
 
 If the answer to step 2 is NO/REJECT, you output the search segment to the ranges.tab file. If no outliers were found within the entire genome (i.e. after the first run of the loop), of course the "search segment" outputted will span the entire genome. In every other case, the search segment outputted to ranges.tab is terminal, meaning that it doesn't have any other outliers within it (although it may itself be outlier).
 
 * T_Score Formula: Please read "Details, Circular Binary Segmentation" above before going any further. The formula for calculating the t_score, or "test statistic", comes from this academic paper: <Biostat (2004) 5 (4): 557-572. doi: 10.1093/biostatistics/kxh008>. Fully understanding the formula isn't necessary to work on the code. If you wish to do so anyway, you should refer the paper, and we borrow some of the same variable names for ease of identification in that respect (e.g. "y", "z", "d1", "d2"). All you need to know is that the simplest steps of the formula determine the avg coverage of the tiles between i and j, and the avg coverage of the tiles outside of i and j, then subtract one from the other to find the difference. Further steps, for example, factor in the width of the potential outlier (between i and j) compared to the width of the entire search segment. The end result is a single number that takes into account all relevent measures of size. Here is the formula:
 
 y = (sum_of_coverage_to_j_pos - sum_of_coverage_to_i_pos) / (j_position - i_position);
 z = (sum_of_coverage_in_segment - sum_of_coverage_to_j_pos + sum_of_coverage_to_i_pos) / (num_positions_in_segment - j_position + i_position);
 d1 = (1.0 / (j_position - i_position));
 d2 = (1.0 / (num_positions_in_segment - j_position + i_position));
 t = (y - z) / sqrt( d1 + d2 );
 
 * Tiles vs Position: Please read documentation for the CoverageDistribution::tile function before going any further. This function keeps track of location within a segment in two ways: by "position" and by "tile". Each entry in our input file, tiled.tab, represents the average coverage of every 500 positions across the entire genome (or whatever tile_size is). So of course each POSITION within tiled.tab will be a multiple of tile_size (i.e. 0, 500, 1000, 1500, etc). Say we're looking in tiled.tab at position 30500 which shows coverage of 60.2. This means that the avg coverage is 60.2 for the previous 500 positions (inclusive), so from 30001 to 30500. Each of these 500 position-long sections is considered a TILE, so there are as many tiles as entries in tiled.tab. Although not explicitly labeled in the file, in the code we refer to each tile with a number in sequence, where the first entry in tiled.tab represents tile 0, the next is tile 1, then tile 2, tile 3, and so on. We do this to make it easier to manipulate the entries later on, like loading them into data structures and iterating through them with loops. NOTE: The first entry in every tiled.tab is a pair of zeros, meaning that at position 0, avg coverage is 0. Although there technically is no position 0, because the coverage data in coverage.tab starts at position 1, we make one up because it allows us to eliminate special cases during for-loops throughout the code. (Coverage.tab is the file used to generate tiled.tab)
 
 * Single-Tile Outliers: The algorithm will have trouble finding outliers that are only one tile wide if they are isolated and not connected to other outlying segments. It is much better at finding these single-tile outliers if they border the edge of another larger outlier that is two or more tiles wide (for an example, look at the HTML output for genome REL10951). This trade-off is intentional because there isn't a good way to find a single-tile outlier: the noise we are trying to screen-out often occurs in spikes that are only one tile wide. The chances of finding something meaningful increase significantly when you require the outlier to be at least two tiles wide because, statistically, it's a much rarer occurrence.
 
 * Floating Point Comparisons (and the inherent problems within): Using comparison operators (i.e. "<", ">", "==", "!=", etc) with floating point numbers like doubles can cause problems because floating point math is not exact. Here is a good url to read about it: <bit.ly/xIW2Mf>. To quote: "Many values cannot be precisely represented using binary floating point numbers, and the limited precision of floating point numbers means that slight changes in the order of operations or the precision of intermediates can change the result." For example, storing the same number, even an integer, within two different floats and comparing to see if they are equal will often return false. This function uses less-than and greater-than floating point comparisons, but the difference between the compared values is usually large enough that it yields the expected results anyway. However, at some point it might be a good idea to follow some of the suggestions in the referenced url and replace the current method of comparing exact values with one that compares ranges, like ULP or FLT_EPSILON.
 */
void CoverageDistribution::find_segments(const Settings& settings,
                                         double summary_average,       // Avg coverage depth
                                         // across entire genome.
                                         string in_file_name,          // tiled.tab
                                         string tiled_for_edging_file_name,// tiled_for_edging.tab
                                         string out_file_name,         // ranges.tab
                                         string history_file_name)     // history.tab

{
  cout << "Finding significant regions..." << endl;
  
  (void)settings;
  
  // A search segment is represented by a pair of tiles: the first tile of the
  // segment and the last tile of the segment (inclusive)
  search_pair_t current_search_pair_by_tile;
  
  // Vectors keep track of every search segment that will be examined for outliers,
  // both at normal tile size and during improved edge detection
  vector< search_pair_t > all_search_pairs_vector;
  vector< search_pair_t > all_search_pairs_vector_for_edging;
  
  // REMOVE FOR BINARY EDGE SEARCH
  // Vectors keep track of every segment that is found to be a terminal "range",
  // that is, a contiguous region with the same copy number, both at normal
  // tile size and during improved edge detection. Each range is thus either
  // an outlier, or the remaining area between or beside outliers. So each range
  // demarcates where an outlier either starts or ends.
  vector< search_pair_t > ranged_pairs_vector;
  vector< search_pair_t > ranged_pairs_vector_copy;
  vector< search_pair_t > ranged_pairs_vector_for_edging;
  
  // For binary edge search:
  // Vectors keep track of every segment that is found to be a terminal "range",
  // that is, a contiguous region with the same copy number, both at normal
  // tile size and during improved edge detection. Each range is thus either
  // an outlier, or the remaining area between or beside outliers. So each range
  // demarcates where an outlier either starts or ends.
  vector< range_t > range_t_vector;
  vector< range_t > range_t_vector_copy;
  
  // Used to add elems to range_t_vector
  // Members:
  //   pair<int32_t, int32_t> range_pair
  //   double copy_number
  //   double copy_number_float
  //   double t_score
  //   double p_value
  range_t range;
  
  // Used to populate data structures needed by the CBS algorithm
  int32_t current_position;
  int32_t previous_position;
  
  // Relative positions of i and j within the search segment.
  // Can be reused for both high and low coveage calculations.
  int32_t i_position;
  int32_t j_position;
  
  // Used to determine absolute position within the entire genome
  int32_t position_hash;
  
  // Defines the beginning and end tiles of potential high-coverage
  // and low-coverage outliers
  int32_t best_i_tile_hi;
  int32_t best_j_tile_hi;
  int32_t best_i_tile_lo;
  int32_t best_j_tile_lo;
  
  // Defines the beginning and end tiles of the highest and lowest t_score within
  // a random segment. Used only for stepping thru the code in debug mode.
  int32_t r_best_i_tile_hi;
  int32_t r_best_j_tile_hi;
  int32_t r_best_i_tile_lo;
  int32_t r_best_j_tile_lo;
  
  // Keeps track of the total number of positions in the search segment
  // and in the entire genome, respectively.
  // i.e. the width of the search segment and genome in bp.
  int32_t num_positions_in_segment;
  int32_t num_positions_in_genome;
  
  // Maps the relative sum of coverage up to every position within the search segment
  map<int32_t, double> sum_of_coverage_to_position;
  double running_sum = 0; // Used to populate above
  
  // Relative sums of coverage up to i_position and j_position within the search segment.
  // Can be reused for both high and low coveage calculations.
  double sum_of_coverage_to_i_pos;
  double sum_of_coverage_to_j_pos;
  
  // Relative sum of coverage within the entire search segment.
  // Can be reused for both high and low coveage calculations.
  double sum_of_coverage_in_segment;
  
  // For keeping track of the highest and lowest t_score within the search segment
  double best_t_score_hi;
  double best_t_score_lo;
  int32_t t_length;
  int32_t r_t_length;
  int32_t best_t_length_hi;
  int32_t best_t_length_lo;
  
  // For use in the t_score formula (see below):
  double y; // Avg coverage inside of i and j
  double z; // Avg coverage outside of i and j
  double d1; // To factor in the width of i and j
  double d2; // To factor in the width of i and j
  double t_score;
  
  /*
   t_score formula (for reference):
   
   y = (sum_of_coverage_to_j_pos - sum_of_coverage_to_i_pos) / (j_position - i_position);
   z = (sum_of_coverage_in_segment - sum_of_coverage_to_j_pos + sum_of_coverage_to_i_pos) / (num_positions_in_segment - j_position + i_position);
   d1 = (1.0 / (j_position - i_position));
   d2 = (1.0 / (num_positions_in_segment - j_position + i_position));
   t_score = (y - z) / sqrt( d1 + d2 );
   */
  
  // For keeping track of the highest and lowest t_score among the random segments
  double r_t_score;
  double r_best_t_score_hi;
  double r_best_t_score_lo;
  int32_t r_best_t_length_hi;
  int32_t r_best_t_length_lo;
  
  // Extra data written to the history.tab file, for debugging purposes
  double mean_inside;
  double mean_left;
  double mean_right;
  int32_t copy_inside;
  int32_t copy_left;
  int32_t copy_right;
  
  // Legacy code used for finding mean error
  //double inside_mean;
  //double outside_mean;
  //double mean_error;
  
  // Used for randomizing a segment
  vector<double> vector_of_randomized_coverage;
  list<int32_t>::iterator list_iter;
  list<int32_t> list_of_tiles_in_segment;
  int32_t random_num;
  int32_t best_t_was_higher_counter;
  int32_t best_t_was_lower_counter;
  int32_t max_error_threshold;
  
  // Used for randomized testing for statistical significance
  int32_t randomizations = 100; // How many times the random test loop will run
  max_error_threshold = int(randomizations * .05) + 1; // Ex: 6
  srand(0);
  
  // Define filestreams
  ifstream tiles_file;
  ifstream tiled_for_edging_file;
  ofstream ranges_file;
  ofstream history_file;
  
  // Open filestreams
  tiles_file.open ( in_file_name.c_str() );
  tiled_for_edging_file.open ( tiled_for_edging_file_name.c_str() );
  ranges_file.open ( out_file_name.c_str() );
  history_file.open ( history_file_name.c_str() );
  
  // Print labels for all the columns in the output files
  ranges_file << "Start_Position\tEnd_Position\tT_Score_Within\tP_Value_Within\tT_Score_Exists\tP_Value_Exists\n";
  
  history_file << "Start_Search\tEnd_Search\tStart_Segment\tEnd_Segment\t"
  << "Mean_Coverage_Inside\tMean_Coverage_Left\tMean_Coverage_Right\t"
  << "CN_Inside\tCN_Left\tCN_Right\t"
  << "T_Score_Within\tP_Value_Within\n";
  
  /*
   * Get data from tiled.tab file
   */
  
  // For reading input from tiled.tab
  int32_t num_tiles = 0;
  int32_t num_tiles_in_genome = 0;
  int32_t position;
  double coverage;
  map<int32_t, double> coverage_map;
  vector<int32_t> ordered_positions_in_segment;
  
  // Skip header of tiled.tab and get tile_size
  string skip;
  int32_t tile_size;
  tiles_file >> skip >> skip >> tile_size;
  
  // Read each line of data (position & coverage) from tiled.tab file
  while ( tiles_file >> position >> coverage )
  {
    // Each position and it's corresponding coverage is populated into this map
    coverage_map[position] = coverage;
    
    // Each position that exists within tiled.tab is populated into this vector.
    // Since the positions are populated in order, each can be accessed by its
    // corresponding tile index.
    ordered_positions_in_segment.push_back(position);
    
    // Count the total number of tiles (lines of data)
    num_tiles++;
  }
  
  
  /*
   * Populate the first search segment
   */
  
  // The first search segment is the entire genome, so the total number of positions
  // in the first seach segment is highest position value in tiled.tab
  num_positions_in_genome = position;
  num_positions_in_segment = position; // NOTE: This may be unnecessary
  num_tiles_in_genome = num_tiles;
  
  // Add the first search segment (the entire genome) to the end of all_search_pairs_vector
  current_search_pair_by_tile.first = 1;
  current_search_pair_by_tile.second = num_tiles - 1;
  current_search_pair_by_tile.t_score_exists = -1.0;
  current_search_pair_by_tile.p_value_exists = -1.0;

  all_search_pairs_vector.push_back(current_search_pair_by_tile);
  
  /*
   * Beginning of CBS algorithm loop
   */
  
  // Until there are no more search segments in which to look for outliers...
  while ( all_search_pairs_vector.size() > 0 ) {
    
    
    /* Set the current search segment: */
    
    // Remove the last search pair from the all_search_pairs_vector and
    // set it as the current search pair. NOTE: using pop_back() is why
    // results are printed to the ranges and history files in seemingly
    // reverse order from which they are added to the queue (later down the line).
    current_search_pair_by_tile = all_search_pairs_vector.back();
    all_search_pairs_vector.pop_back();

    int32_t start_coord = ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size + 1;
    int32_t end_coord = ordered_positions_in_segment[current_search_pair_by_tile.second];

    // Print to stdout
    cout << "  REGION: " << setw(8) << right << start_coord << " - " << setw(8) << right << end_coord << endl;
    
    // If the current search segment consists of 1 tile, it cannot be examined for
    // outliers because no outliers can exist within it.
    if (current_search_pair_by_tile.first == current_search_pair_by_tile.second)
    {
      // Print to ranges file
      ranges_file << start_coord << "\t"
      << end_coord << "\t"
      << 0.0 << "\t"
      << 1.0 << "\t"
      << 0.0 << "\t"
      << 1.0 << endl;

      // Print to history file
      history_file << start_coord << "\t"
      << end_coord << "\t"
      << start_coord << "\t"
      << end_coord << "\t"
      << -1 << "\t"
      << -1 << "\t"
      << -1 << "\t"
      << -1 << "\t"
      << -1 << "\t"
      << -1 << "\t"
      << 0.0 << "\t"
      << 1.0 << endl;
      continue;
    }
    
    /*
     * Prep work to populate data structures needed by CBS algorithm
     */
    
    // VERBOSE, set to TRUE for debugging output
    const bool VERBOSE = false;
    
    if (VERBOSE) {
      cout << "\n############## VERBOSE OUTPUT ###############" << endl;
      cout << "\ncurrent_search_pair_by_tile_index: [" << current_search_pair_by_tile.first
      << "," << current_search_pair_by_tile.second << "]" << endl;
    }
    
    // Set to determine absolute position within the entire genome.
    // Ex: Assuming tile_size=500, if the first tile of the search segment is 61,
    // and the position value of that tile is 30500, then the position_hash is the
    // position value of the tile whose index is 60, which is 30000. */
    position_hash = ordered_positions_in_segment[current_search_pair_by_tile.first - 1];
    
    if (VERBOSE) {
      cout << "position_hash: " << position_hash << endl;
      cout << "Press ENTER to continue..." << endl;
      cin.get();
    }
    
    // Clear map before populating it
    sum_of_coverage_to_position.clear();
    
    // Keep track of our relative position within the search segment. This always sets
    // current_position to 0.
    current_position = ordered_positions_in_segment[current_search_pair_by_tile.first - 1] - position_hash;
    
    // The sum of coverage at position 0 is always 0
    sum_of_coverage_to_position[current_position] = 0;
    num_tiles = 0;
    
    // This loop is the prep work. Before you can search for outliers in this
    // segment, you must figure out how many positions it has (i.e. base pairs),
    // and you must populate a map with the sum of coverage up to each position.
    // NOTE: Don't forget we are treating this current search segment as its own
    // entity. Thus, no matter where this segment is within the actual genome,
    // we are treating the first position of this segment as zero.
    for (int32_t tile_index = current_search_pair_by_tile.first;
         tile_index <= current_search_pair_by_tile.second; tile_index++)
    {
      current_position = ordered_positions_in_segment[tile_index] - position_hash;
      previous_position = ordered_positions_in_segment[tile_index - 1] - position_hash;
      
      sum_of_coverage_to_position[current_position] = coverage_map[current_position + position_hash] + sum_of_coverage_to_position[previous_position];
      num_tiles++;
      
      if (VERBOSE) {
        cout << "**********INSIDE PREP-LOOP*************" << endl;
        cout << "current_position:" << current_position << " "
        << "= ordered_positions_in_segment[tile_index]:" << ordered_positions_in_segment[tile_index] << " "
        << "- position_hash:" << position_hash << endl;
        cout << "previous_position:" << previous_position << " "
        << "= ordered_positions_in_segment[tile_index-1]:" << ordered_positions_in_segment[tile_index - 1] << " "
        << "- position_hash:" << position_hash << endl;
        cout << "sum_of_coverage_to_position[current_position]:" << sum_of_coverage_to_position[current_position] << " "
        << "= coverage_map[current_position+position_hash]:"
        << coverage_map[current_position + position_hash] << " "
        << "+ sum_of_coverage_to_position[previous_position]:"
        << sum_of_coverage_to_position[previous_position] << endl;
        cout << "num_tiles: " << num_tiles << endl;
      }
    }
    
    
    
    
    
    /*
     * Step 1 of algorithm: Find the region in the search segment with the largest t_score
     */
    
    // Set the current_position to the last position in this search segment.
    // NOTE: This may be redundant with the first line in the above prep work loop
    current_position = ordered_positions_in_segment[current_search_pair_by_tile.second] - position_hash;
    
    // Set values for upcoming step 1 of algorithm
    sum_of_coverage_in_segment = sum_of_coverage_to_position[current_position];
    num_positions_in_segment = current_position;
    
    // Set to largest possible negative number so the first t_score found
    // is always greater than the initialized value
    best_t_score_hi = -1.7e-308;
    
    // Zero out values
    best_t_was_higher_counter = 0;
    best_t_length_hi = 0;
    
    
    
    // For every possible value of i within the current search segment...
    for (int32_t i_tile = current_search_pair_by_tile.first - 1;
         i_tile < current_search_pair_by_tile.second; i_tile++)
    {
      /* NOTE: Setting i_tile to [current_search_pair.first - 1] works above because
       every search segment always starts at zero in both position and coverage. Thus,
       say you're analyzing tiles 61 thru 97, you don't have to worry that the
       coverage data from tile 60 is being used because that data is at i_position 0
       (relative to this segment), which always has a coverage value of 0. */
      
      // These always start at zero and increase with each iteration of the for-loop
      i_position = ordered_positions_in_segment[i_tile] - position_hash;
      sum_of_coverage_to_i_pos = sum_of_coverage_to_position[ i_position ];
      
      
      
      // For every possible value of j within the current search segment...
      for (int32_t j_tile = i_tile + 1;
           j_tile <= current_search_pair_by_tile.second; j_tile++)
      {
        // These always start at zero and increase with each
        // iteration of the for-loop
        j_position = ordered_positions_in_segment[j_tile] - position_hash;
        sum_of_coverage_to_j_pos = sum_of_coverage_to_position[ j_position ];
        
        // Calculate the t_score of the subsegment with the given i and j
        y = (sum_of_coverage_to_j_pos - sum_of_coverage_to_i_pos) / (j_position - i_position);
        z = (sum_of_coverage_in_segment - sum_of_coverage_to_j_pos + sum_of_coverage_to_i_pos) / (num_positions_in_segment - j_position + i_position);
        d1 = (1.0 / (j_position - i_position));
        d2 = (1.0 / (num_positions_in_segment - j_position + i_position));
        t_score = (y - z) / sqrt( d1 + d2 );
        
        // Determine number of tiles between i and j
        t_length = j_tile - i_tile;
        
        // Keep track of subsegment with the largest t_score.
        // Note: Floats are being compared here, and thus the equality
        // comparison ("==") will likely never evaluate to true.
        if ( abs(t_score) > best_t_score_hi ||
            ( abs(t_score) == best_t_score_hi && t_length > best_t_length_hi) )
        {
          best_t_score_hi = abs(t_score);
          best_i_tile_hi = i_tile;
          best_j_tile_hi = j_tile;
          best_t_length_hi = t_length;
        }
        
      } // End of j_tile for-loop
    } // End of i_tile for-loop
    
    // Save info about the potential outlier that has just been found
    i_position = ordered_positions_in_segment[best_i_tile_hi] - position_hash;
    j_position = ordered_positions_in_segment[best_j_tile_hi] - position_hash;
    sum_of_coverage_to_i_pos = sum_of_coverage_to_position[ i_position ];
    sum_of_coverage_to_j_pos = sum_of_coverage_to_position[ j_position ];
    
    /*
     * Make debugging info for the history file.
     * Note: Doing calculations for the left and right subsegments
     * may no longer be applicable. This code should be refactored
     * if this is the case.
     */
    
    // Calculate the mean coverage between i and j
    mean_inside = (sum_of_coverage_to_j_pos - sum_of_coverage_to_i_pos) / double(best_j_tile_hi - best_i_tile_hi);
    
    // Calculate the copy number between i and j by dividing mean_inside by the avg
    // coverage across the entire genome, and then rounding to nearest whole value
    copy_inside = static_cast<int32_t>(floor((mean_inside / summary_average) + .5));
    
    // If the left edge of the potential outlier equals the left edge of the
    // search segment, there's no left subsegment in which to search for further outlers
    if (best_i_tile_hi == current_search_pair_by_tile.first - 1) {
      mean_left = -1;
      copy_left = -1;
    }
    // Otherwise, the subsegment to the left will be added to the queue
    else {
      mean_left = sum_of_coverage_to_i_pos / double(best_i_tile_hi);
      copy_left = static_cast<int32_t>(floor((mean_left / summary_average) + .5));
    }
    
    // If the right edge of the potential outlier equals the right edge of the
    // search segment, there's no right subsegment in which to search for further outlers
    if ( best_j_tile_hi == current_search_pair_by_tile.second ) {
      mean_right = -1;
      copy_right = -1;
    }
    // Otherwise, the subsegment to the right will be added to the queue
    else {
      mean_right = (sum_of_coverage_in_segment - sum_of_coverage_to_j_pos) / double(current_search_pair_by_tile.second - best_j_tile_hi);
      copy_right = static_cast<int32_t>(floor((mean_right / summary_average) + .5));
    }
    
    
    /*
     * Step 2 of algorithm: Randomized testing of potential outlier
     */
    
    /* Conceptually, this for-loop uses a randomized testing procedure to generate
     data that will later be used to determine if the best t_score from the
     search segment is statistically significant (i.e. is an outlier).
     
     The loop runs 100 times, each time comparing the best t_score of the search
     segment to the best t_score of a new segment randomly generated by shuffling the
     tiles of the search segment. It counts how many times the search segment's best t_score
     is larger than the random best t_score.
     
     NOTE: The actual check for statistical significance does not happen in this for-loop.
     Later on it will determine if the current best t_score was larger than the random
     best t_score 95 times out of 100; if yes, the current best t_score is "accepted"
     as an outlier. */
    
    for ( int32_t examinations = 0; examinations < randomizations; examinations++ )
    {
      
      /*
       * Generate random segment
       */
      
      /* Below is an unusual yet efficient way to randomly select tiles within the
       current search segment and then populate them into a new random segment to be
       used for testing purposes: */
      
      list_of_tiles_in_segment.clear();
      vector_of_randomized_coverage.clear();
      
      // First, create a reference list of all the tile indexes in the
      // current search segment that will need to be shuffled:
      // For each tile index value in the current segment (ex: 61 thru 97)...
      for (int32_t tile_index = current_search_pair_by_tile.first;
           tile_index <= current_search_pair_by_tile.second; tile_index++)
      {
        // Populate the tile index into a linked list
        list_of_tiles_in_segment.push_back(tile_index);
      }
      
      // Next, we need to randomly select a tile in the current search segment and
      // populate its corresponding coverage value into a vector for later use:
      // Run the following for-loop the same number of times as the number of
      // tiles in the current search segment...
      for (int32_t tile_count = current_search_pair_by_tile.first;
           tile_count <= current_search_pair_by_tile.second; tile_count++)
      {
        // Generate a random number that is less than the total number of
        // tiles in the current segment (but greater than or equal to zero)
        random_num = rand() % list_of_tiles_in_segment.size();
        
        /* The next for-loop is how the random tile is actually chosen. It is not
         designed to run from beginning to end. Instead, it is meant to run a random
         number of times as specified by the random_num above. It iterates thru the
         list of tile indexes, decrementing the random_num by 1 each time, and
         whatever tile the iterator is pointing to when random_num = 0 becomes the
         next randomly chosen tile. */
        for (list_iter = list_of_tiles_in_segment.begin();
             list_iter != list_of_tiles_in_segment.end(); list_iter++)
        {
          if ( random_num == 0 )
          {
            /* Note that list_of_tiles_in_segment represents the same tiles of the
             search segment that exist within ordered_positions_in_segment (ex: 61
             thru 97). Thus, the list_iter is just a pointer to a value (a tile
             index, ex: 63), and that value is accessed by *list_iter. Using this
             tile to access a specific element within ordered_positions_in_segment
             returns the corresponding position (ex: 31500), and then this position
             is used to access the corresponding coverage at this position within
             coverage_map. Finally, this average coverage value is pushed onto the
             end of the vector_of_randomized_coverage. */
            vector_of_randomized_coverage.push_back(coverage_map[ordered_positions_in_segment[*list_iter]]);
            
            // To prevent this same coverage value from being used again, its
            // corresponding tile is removed from the list_of_tiles_in_segment list.
            list_of_tiles_in_segment.erase(list_iter);
            break;
          }
          
          // Decrement random_num by 1
          random_num--;
          
        }
      } // End of tile_count for-loop

      
      /*
       * More prep work (still in randomized testing)
       */
      
      // This loop repopulates the sum_of_coverage_to_position vector with the randomly
      // chosen values from vector_of_randomized_coverage.
      // Note: sum_of_coverage_to_position does not need to be zeroed-out before this
      // loop is run; see "sum_of_coverage_to_position" below...
      for (size_t iter = 0; iter < vector_of_randomized_coverage.size(); iter++)
      {
        // current_position always starts at tile_size and then increases
        // by tile_size with every loop (ex: 500, 1000, 1500, ...)
        current_position = ordered_positions_in_segment[current_search_pair_by_tile.first + iter] - position_hash;
        
        // ?2? previous_position always starts at zero and then increases
        // by tile_size with every loop (ex: 0, 500, 1000, ...)
        previous_position = ordered_positions_in_segment[current_search_pair_by_tile.first + iter - 1] - position_hash;
        
        // The first time through this loop (i.e. when iter=0 and previous_position=0),
        // the sum_of_coverage_to_position[previous_position] will always
        // equal zero, because sum_of_coverage_to_position[0] will always equal zero.
        sum_of_coverage_to_position[current_position] = sum_of_coverage_to_position[previous_position] + vector_of_randomized_coverage[iter];
      }
      
      /*
       * Find the best t_score within the random segment
       */
      
      // Set to largest possible negative number so the first t_score found
      // is always greater than the initialized value
      r_best_t_score_hi = -1.7e-308;
      
      // For every possible value of i within the current search segment...
      for ( int32_t i_tile = current_search_pair_by_tile.first - 1;
           i_tile < current_search_pair_by_tile.second; i_tile++ )
      {
        // These always start at zero and increase with each iteration of the for-loop
        i_position = ordered_positions_in_segment[i_tile] - position_hash;
        sum_of_coverage_to_i_pos = sum_of_coverage_to_position[ i_position ];

        // For every possible value of j within the current search segment...
        for ( int32_t j_tile = i_tile + 1;
             j_tile <= current_search_pair_by_tile.second; j_tile++ )
        {
          // These always start at zero and increase with each
          // iteration of the for-loop
          j_position = ordered_positions_in_segment[j_tile] - position_hash;
          sum_of_coverage_to_j_pos = sum_of_coverage_to_position[ j_position ];
          
          // Calculate the t_score of the subsegment with the given i and j
          y = (sum_of_coverage_to_j_pos - sum_of_coverage_to_i_pos) / (j_position - i_position);
          z = (sum_of_coverage_in_segment - sum_of_coverage_to_j_pos + sum_of_coverage_to_i_pos) / (num_positions_in_segment - j_position + i_position);
          d1 = 1.0 / (j_position - i_position);
          d2 = 1.0 / (num_positions_in_segment - j_position + i_position);
          r_t_score = (y - z) / sqrt( d1 + d2 );
          r_t_score = abs(r_t_score);
          
          
          // Keep track of subsegment with the largest t_score
          // Note: Floats are being compared here; also, an "equality" comparison
          // (that actually compares within a range) might need to be added so that
          // it mirrors the methodology of the similar code block near line 740.
          if ( r_t_score > r_best_t_score_hi )
          {
            r_best_t_score_hi = abs(r_t_score);
            
            // Used for stepping thru the code in debug mode
            r_best_i_tile_hi = i_tile;
            r_best_j_tile_hi = j_tile;
          }
          
        } // End of for-loop j_tile
      } // End of for-loop i_tile
      
      //  Now we know the best t_score of the random segment
      if (VERBOSE) {
        cout << "\n\n#####################################";
        cout << "\n best_t_score: " << fixed << setprecision(40) << best_t_score_hi;
        cout << "\r_best_t_score: " << fixed << setprecision(40) << r_best_t_score_hi
        << endl;
      }
      
      /*
       * Count the times best_t_score is larger than r_best_t_score
       */
      
      if ( best_t_score_hi > r_best_t_score_hi )
      {
        best_t_was_higher_counter++;
      }
      
      // This heuristic breaks out of the randomized testing early if it becomes
      // impossible to reach the number of counts required to be significant.
      // Ex: if ([(96 + 1) - 90] > 6) { ...
      // NOTE: This could probably be simplified to the following:
      /*  if (examinations - best_t_was_bigger_counter >= max_error_threshold){ ... */
      
      /*
      if ((examinations + 1) - best_t_was_higher_counter > max_error_threshold){
        
        // Break the entire randomized testing loop
        break;
      }
      */
      
    } // End of randomized testing (aka "examinations") for-loop
    
    
    /*
     * If answer to step 2 is YES/ACCEPT (the potential outlier is significant)
     */
    
    double t_score_within = best_t_score_hi;
    double p_value_within = 1 - (double)(best_t_was_higher_counter) / randomizations;
    
    // If the potential outlier is significant, and if the potential outlier is NOT
    // just the entire genome...
    if ( (p_value_within <= 0.05) &&
        !((current_search_pair_by_tile.first == best_i_tile_hi + 1) && (current_search_pair_by_tile.second == best_j_tile_hi))
        )
    {
      
      /*
       * Break search segment into three subsegments and add each to search queue
       */
      
      // Record this within the history file
      history_file << start_coord << "\t" <<
      end_coord << "\t" <<
      ordered_positions_in_segment[best_i_tile_hi + 1] - tile_size + 1 << "\t" <<
      ordered_positions_in_segment[best_j_tile_hi] << "\t" <<
      mean_inside << "\t" <<
      mean_left << "\t" <<
      mean_right << "\t" <<
      copy_inside << "\t" <<
      copy_left << "\t" <<
      copy_right << "\t" <<
      t_score_within << "\t" <<
      p_value_within << endl;
      
      // Add the middle subsegment (which is the newly found outlier between i and j)
      // to the search queue (i.e. all_search_pairs_vector)
      search_pair_t middle_search_pair_by_tile;
      middle_search_pair_by_tile.first = best_i_tile_hi + 1;
      middle_search_pair_by_tile.second = best_j_tile_hi;
      middle_search_pair_by_tile.t_score_exists = t_score_within;
      middle_search_pair_by_tile.p_value_exists = p_value_within;
      all_search_pairs_vector.push_back( middle_search_pair_by_tile );
      
      // Add the left subsegment (if there is one) to the search queue.
      // i.e. If the middle subsegment doesn't start at the same tile
      // as the search segment.
      if ( best_i_tile_hi + 1 > current_search_pair_by_tile.first )
      {
        search_pair_t left_search_pair_by_tile;
        left_search_pair_by_tile.first = current_search_pair_by_tile.first;
        left_search_pair_by_tile.second = best_i_tile_hi;
        left_search_pair_by_tile.t_score_exists = t_score_within;
        left_search_pair_by_tile.p_value_exists = p_value_within;
        all_search_pairs_vector.push_back( left_search_pair_by_tile );
      }
      
      // Add the right subsegment (if there is one) to the search queue.
      // i.e. If the middle subsegment doesn't end at the same tile
      // as the search segment.
      if ( best_j_tile_hi < current_search_pair_by_tile.second )
      {
        search_pair_t right_search_pair_by_tile;
        right_search_pair_by_tile.first = best_j_tile_hi + 1;
        right_search_pair_by_tile.second = current_search_pair_by_tile.second;
        right_search_pair_by_tile.t_score_exists = t_score_within;
        right_search_pair_by_tile.p_value_exists = p_value_within;
        all_search_pairs_vector.push_back( right_search_pair_by_tile );
      }
      
    } // End of YES/ACCEPT if-statement
    
    
    /*
     * If answer to step 2 is NO/REJECT (the potential outlier is NOT significant)
     */
    
    /* Output the search segment to the ranges.tab file. If no outliers were found
     within the entire genome (i.e. after the first run of the loop), of course the
     "search segment" outputted will span the entire genome. In every other case,
     the search segment outputted to ranges.tab is terminal, meaning that it doesn't
     have any other outliers within it (although it may itself be outlier). */
    else
    {
      
      /* Print to ranges file */
      
      // If the only "range" found spans the entire genome because no outliers
      // were found, don't perform improved edge detection. Instead, just print
      // the range to the ranges file, and end the function.
      if (current_search_pair_by_tile.first == 1 &&
          current_search_pair_by_tile.second == num_tiles_in_genome - 1)
      {
        ranges_file << ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size + 1 << "\t"
        << ordered_positions_in_segment[current_search_pair_by_tile.second] << "\t"
        << t_score_within << "\t"
        << p_value_within << endl;
        
        tiles_file.close();
        tiled_for_edging_file.close();
        ranges_file.close();
        history_file.close();
        
        return;
      }
      
      
      // Save current search pair as a range for later binary edge search.
      // Note that we convert the location in genome from tile to position
      range.range_pair.first = ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size;
      range.range_pair.second = ordered_positions_in_segment[current_search_pair_by_tile.second];
      range.t_score_within = t_score_within;
      range.p_value_within = p_value_within;
      range.t_score_exists = current_search_pair_by_tile.t_score_exists;
      range.p_value_exists = current_search_pair_by_tile.p_value_exists;
      range_t_vector.push_back(range);
      
      
      // REMOVE FOR BINARY EDGE SEARCH
      /* Because the edges of the formerly found parent ranges will probably shift
       due to the improved detection process, we should discard these two edge tiles
       from their respective ranges (i.e. cut each range short by one tile on each end),
       and then later make up for it (within ::smooth_segments function) by merging the
       ranges found at the smaller tile size with their bordering parent ranges based on
       having the same copy number. */
      
      /*
       // If the first tile in the range is the first tile of the genome,
       // chop off ONLY the right side of the parent range and print to ranges file.
       if (current_search_pair_by_tile.first == 1)
       {
       int32_t right_side_without_edge = current_search_pair_by_tile.second - 1;
       
       ranges_file << ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size + 1 << "\t"
       << ordered_positions_in_segment[right_side_without_edge] << "\t"
       << best_t_score_hi << "\t"
       << 1 - (double)(best_t_was_higher_counter) / randomizations << endl;
       }
       
       // Else, if the last tile in the range is the last tile of the genome,
       // chop off ONLY the left side of the parent range and print to ranges file.
       else if (current_search_pair_by_tile.second == num_tiles_in_genome - 1)
       {
       int32_t left_side_without_edge = current_search_pair_by_tile.first + 1;
       
       ranges_file << ordered_positions_in_segment[left_side_without_edge] - tile_size + 1 << "\t"
       << ordered_positions_in_segment[current_search_pair_by_tile.second] << "\t"
       << best_t_score_hi << "\t"
       << 1 - (double)(best_t_was_higher_counter) / randomizations << endl;
       }
       
       // Else, chop off BOTH the left and right sides of the parent range
       // and print to ranges file.
       else
       {
       int32_t left_side_without_edge = current_search_pair_by_tile.first + 1;
       int32_t right_side_without_edge = current_search_pair_by_tile.second - 1;
       
       ranges_file << ordered_positions_in_segment[left_side_without_edge] - tile_size + 1 << "\t"
       << ordered_positions_in_segment[right_side_without_edge] << "\t"
       << best_t_score_hi << "\t"
       << 1 - (double)(best_t_was_higher_counter) / randomizations << endl;
       }
       */
      
      // Print ORIGINAL, unmodified segments to history file for debugging
      history_file << ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size + 1 << "\t" <<
      ordered_positions_in_segment[current_search_pair_by_tile.second] << "\t" <<
      ordered_positions_in_segment[best_i_tile_hi + 1] - tile_size + 1 << "\t" <<
      ordered_positions_in_segment[best_j_tile_hi] << "\t" <<
      mean_inside << "\t" <<
      mean_left << "\t" <<
      mean_right << "\t" <<
      copy_inside << "\t" <<
      copy_left << "\t" <<
      copy_right << "\t" <<
      t_score_within << "\t" <<
      p_value_within << endl;
      
      
    } // End of NO/REJECT else-statement
    
    
  } // End of CBS algorithm while-loop
  
  
  // Note: range_t_vector should have a size of at least 2 or 3 at
  // this point: one outlier, and one or two ranges that represent the rest
  // of the genome (depending on whether the outlier starts or ends on
  // an edge of the genome).
  
  
  // Copy range_t_vector for safekeeping
  range_t_vector_copy = range_t_vector;
  
  /**
   ** Sort the ranges by position (base pair) in genome, lowest to highest.
   ** (This works because the "<" operator is overloaded in the range_t struct
   ** in coverage_distribution.h)
   **/
  sort(range_t_vector.begin(), range_t_vector.end());
  
  /**
   **
   **
   ** Time to apply binary edge search to previously found outliers
   **
   **
   **/
  
  /*
   * Get data from tiled_for_edging.tab file
   */
  
  // For reading input from file
  num_tiles = 0;
  uint32_t num_large_tiles_in_genome = num_tiles_in_genome;
  num_tiles_in_genome = 0;
  position = 0;
  coverage = 0;
  coverage_map.clear();
  ordered_positions_in_segment.clear();
  
  // Skip header of file and get tile_size
  skip = "";
  int32_t old_tile_size = tile_size;
  tile_size = 0;
  tiled_for_edging_file >> skip >> skip >> tile_size;
  
  // Read each line of data (position & coverage) from file
  while ( tiled_for_edging_file >> position >> coverage )
  {
    // Each position and it's corresponding coverage is populated into this map
    coverage_map[position] = coverage;
    
    // Each position that exists within file is populated into this vector.
    // Since the positions are populated in order, each can be accessed by its
    // corresponding tile index.
    ordered_positions_in_segment.push_back(position);
    
    // Count the total number of tiles (lines of data)
    num_tiles++;
  }
  
  cout << "Optimizing edges of significant regions..." << endl;
  
  // Reset these values based on tiled_for_edging file
  num_tiles_in_genome = num_tiles;
  num_positions_in_genome = position; // However this value should not change
  
  
  /*
   * Analyze the 1000 bp range that surrounds that edge of each outlier (aka the "edging range") in order to narrow in on where the true edge should be.
   *
   * Note: In loops that calculate avg coverage, we count in increments using tile_size
   */
  
  // Initial search radius in base pairs (around the old edge from which the new edge will be searched for)
  int32_t search_radius = 200;
  
  int32_t edge_range_start = -1; // initialize to impossible value
  int32_t edge_range_end = -1;
  int32_t edge_range_middle = -1;
  int32_t edge = -1;
  string high_side = "";
  
  // For each edge of a range in range_t_vector (except the first, bc it is just the left side of the genome, so we start "i" at 1 instead of 0)
  for (uint32_t i = 1; i < range_t_vector.size(); i++)
  {
    
    // Store the position of this edge using the range_t_vector[i].range_pair.first
    edge = range_t_vector[i].range_pair.first;
    
    // Calculate the start and end of the edging range by subtracting and adding 500 bp to the edge, respectively
    edge_range_start = edge - search_radius;
    edge_range_end =  edge + search_radius;
    
    // But, ignore this edge (and skip to next one) if it's within 500 bp of the first or last position in the genome, because data near the beginning and end of the genome is often inaccurate, and it makes the code cleaner to always have a 1000 bp edging range to analyze
    if (edge_range_start < 0 || edge_range_end > num_positions_in_genome)
    {
      continue;
    }
    
    
    
    
    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    // EXPERIMENTAL:
    // Alternate way of finding edge, by taking a running avg of both five tiles to the left and five tiles to the right in a position for-loop, and declaring the new edge as the position where the difference (abs val) between the sums is greatest.
    
    double diff_best = -1.0; // initialize to impossible value
    int32_t position_best = -1; // initialize to impossible value
    
    for (int32_t position = edge_range_start + (tile_size * 2);
         position <= edge_range_end - (tile_size * 3);
         position += tile_size)
    {
      double avg_cov_left = 0.0;
      avg_cov_left = (
                      coverage_map[position] +
                      coverage_map[position - (tile_size * 1)] +
                      coverage_map[position - (tile_size * 2)] //+
                      //coverage_map[position - (tile_size * 3)] +
                      //coverage_map[position - (tile_size * 4)]
                      ) / 5.0;
      
      double avg_cov_right = 0.0;
      avg_cov_right = (
                       coverage_map[position + (tile_size * 1)] +
                       coverage_map[position + (tile_size * 2)] +
                       coverage_map[position + (tile_size * 3)] //+
                       //coverage_map[position + (tile_size * 4)] +
                       //coverage_map[position + (tile_size * 5)]
                       ) / 5.0;
      
      double diff = abs(avg_cov_left - avg_cov_right);
      if (diff > diff_best)
      {
        diff_best = diff;
        position_best = position;
      }
    }
    
    edge_range_middle = position_best;
    
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    
    
    
    
    // Determine which side of the edging range is the high side (higher coverage) and which is the low (the inside and outside of an outlier's edge, respectively). To do this, we need to find the avg coverage of both the left and right halves (each 500 bp long), and the side with the higher number will be the high side.
    /*
     double sum = 0.0;
     double avg_coverage_left = 0.0;
     double avg_coverage_right = 0.0;
     for (int32_t position = edge_range_start; position <= edge; position += tile_size)
     {
     sum += coverage_map[position];
     }
     avg_coverage_left = sum / (search_radius / tile_size);
     
     sum = 0.0;
     for (int32_t position = edge; position <= edge_range_end; position += tile_size)
     {
     sum += coverage_map[position];
     }
     avg_coverage_right = sum / (search_radius / tile_size);
     
     if (avg_coverage_left > avg_coverage_right)
     high_side = "left";
     else
     high_side = "right";
     
     
     // Find both the highest and lowest coverage among all the 10bp tiles. Then find the avg between these two numbers. This is what the ideal avg SHOULD BE if we've found the new edge of the range (aka the true edge). Now calculate the avg coverage for the entire edging range. If it's greater than the ideal avg, it's an overhang, so narrow search to the outer half of the range. If it's less than the ideal avg, it's an underhang, so narrow search to the inner half of the range. But if it's very close to the ideal, narrow search to the middle. Repeat this todo step for each narrowed range (don't forget to recalculate the ideal avg each time we narrow) until there are only two 10bp tiles left (middle is new edge)(should always be two tiles remaining since starting from an even 1000 and halving each time).
     edge_range_middle = edge;
     int32_t edge_range_start_prev = -1; // initialized to non-possible base pair value
     int32_t edge_range_end_prev = -1; // initialized to non-possible base pair value
     //bool is_initial_center_zoom = true;
     while ((edge_range_end - edge_range_start > tile_size)
     && !(edge_range_start == edge_range_start_prev && edge_range_end == edge_range_end_prev))
     {
     
     edge_range_start_prev = edge_range_start;
     edge_range_end_prev = edge_range_end;
     
     double highest_coverage = 0.0;
     double lowest_coverage = 1E+37; // largest possible double
     double sum = 0.0;
     for (int32_t position = edge_range_start; position <= edge_range_end; position += tile_size)
     {
     if (coverage_map[position] > highest_coverage)
     highest_coverage = coverage_map[position];
     if (coverage_map[position] < lowest_coverage)
     lowest_coverage = coverage_map[position];
     
     sum += coverage_map[position];
     }
     
     double ideal_avg_coverage = (highest_coverage + lowest_coverage) / 2.0;
     double avg_coverage = sum / ((edge_range_end - edge_range_start) / tile_size);
     
     
     // If the avg coverage within the edging range is very close to the avg between the highest and lowest points of coverage, zoom in on middle of range
     double margin_of_error = 0.05;
     
     //if (is_initial_center_zoom)
     //    margin_of_error = 0.2;
     
     int32_t edge_range_width = edge_range_end - edge_range_start;
     int32_t edge_range_quarter = edge_range_width / 4;
     if (abs(ideal_avg_coverage - avg_coverage) < (margin_of_error * ideal_avg_coverage))
     {
     // Find 1/4 of edge range width, rounded to nearest multiple of tile_size
     
     int32_t remainder = edge_range_quarter % tile_size;
     if (remainder != 0)
     {
     if (remainder < tile_size / 2)
     edge_range_quarter = edge_range_quarter - tile_size + remainder;
     else
     edge_range_quarter = edge_range_quarter + tile_size - remainder;
     }
     
     // Zoom in on the middle of the edging range
     edge_range_start = edge_range_start + edge_range_quarter;
     edge_range_end = edge_range_end - edge_range_quarter;
     }
     
     // If it's an overhang, narrow search to outer half of edging range
     else if (avg_coverage > ideal_avg_coverage)
     {
     //is_initial_center_zoom = false;
     
     // If high side is the left
     if (high_side == "left")
     {
     edge_range_start = edge_range_middle;
     //edge_range_start = edge_range_start + edge_range_quarter;
     }
     // Else, since high side is the right
     else
     {
     edge_range_end = edge_range_middle;
     //edge_range_end = edge_range_end - edge_range_quarter;
     }
     
     }
     // Else it's an underhang, so narrow search to inner half of edging range
     else
     {
     //is_initial_center_zoom = false;
     
     // If high side is the left
     if (high_side == "left")
     {
     edge_range_end = edge_range_middle;
     //edge_range_end = edge_range_end - edge_range_quarter;
     }
     // Else, since high side is the right
     else
     {
     edge_range_start = edge_range_middle;
     //edge_range_start = edge_range_start + edge_range_quarter;
     }
     }
     
     // Recalculate the middle of the range, rounding to the nearest multiple of tile_size
     edge_range_middle = (edge_range_start + edge_range_end) / 2;
     int32_t remainder = edge_range_middle % tile_size;
     if (remainder != 0)
     {
     if (remainder < tile_size / 2)
     edge_range_middle = edge_range_middle - tile_size + remainder;
     else
     edge_range_middle = edge_range_middle + tile_size - remainder;
     }
     
     //cout << "\nedge_range_start: " << edge_range_start << endl;
     //cout << "edge_range_middle: " << edge_range_middle << endl;
     //cout << "edge_range_end: " << edge_range_end << endl << endl;
     }
     */
    
    //cout << "\n\n************************" << endl << endl;
    
    
    // Shift the edge of this outlier by rewriting the start and end points based on the binary edge search
    range_t_vector[i - 1].range_pair.second = edge_range_middle;
    range_t_vector[i].range_pair.first = edge_range_middle + 1;
    
  } // Repeat above FOR LOOP for each additional edge within the range_t_vector
  
  
  // Now write the newly shifted original outliers to the ranges file
  for (uint32_t i = 0; i < range_t_vector.size(); i++)
  {
    ranges_file << range_t_vector[i].range_pair.first << "\t"
    << range_t_vector[i].range_pair.second << "\t"
    << range_t_vector[i].t_score_within << "\t"
    << range_t_vector[i].p_value_within << "\t"
    << range_t_vector[i].t_score_exists << "\t"
    << range_t_vector[i].p_value_exists << endl;
  }
  
  // Close the filestreams
  tiles_file.close();
  tiled_for_edging_file.close();
  ranges_file.close();
  history_file.close();
  
}



/*
 * Function: smooth_segments
 * --------------------------------
 * Reads in [genome].ranges.tab and generates [genome].smoothed_ranges.tab,
 * [genome].cnv_final.tab, and [genome].cn_evidence.gd.
 
 * Summary: Takes the segments (aka ranges) within ranges.tab and calculates
 * the copy number within each, outputting this to smoothed_ranges.tab. This is done
 * by taking the avg coverage within each range and dividing it by the avg coverage
 * across the entire genome (see note), then rounding to the nearest whole number.
 * NOTE: This algorithm could be improved by giving a range of possible copy numbers
 * when appropriate, i.e. 150/60 = 2.5, so coverage could be either 2 or 3.
 *
 * THE DOCUMENTATION FOR THIS FUNCTION IS A WORK IN PROGESS
 */

void CoverageDistribution::smooth_segments(
                                           const Settings& settings,
                                           const string& seq_id,
                                           double summary_average,   // @DTF: avg coverage over entire genome
                                           string tile_file_name,
                                           string segment_file_name,
                                           string out_file_name,
                                           string final_file_name,
                                           string gd_file_name
                                           )
{
  cout << "Smoothing significant regions..." << endl;

  
  cGenomeDiff gd; /* cGenomeDiff is a data type defined by breseq. it holds
                   entries like SNP, INS, DEL, CNV, and has more information about
                   those types like their sizes. CNV (copy number variation)
                   is the type that will be added here. CNV is the outliers
                   you're finding. */
  
  //used for finding the mean of a segment in the tile file.
  double segment_sum;
  double segment_mean;
  uint32_t num_tiles;
  
  //keeps track of a segment's bounds.
  uint32_t position_start;
  uint32_t position_end;
  
  //keeps track of an entry in the tile file.
  uint32_t tile_position;
  double tile_coverage;
  double t_score_exists;
  double p_value_exists;
  double t_score_within;
  double p_value_within;
  
  //holds the value a segment's coverage will be changed to.
  double copy_number;
  double copy_number_float;
  
  //used for easily accessing the tile file.
  pair<uint32_t, double> tile_entry;
  vector< pair<uint32_t, double> > tile_data;
  
  //initial parameter.
  uint32_t tile_size;
  
  //used for skipping lines.
  string skip;
  
  bool written_final;
  
  // Used to merge ranges with the same copy number.
  // Members:
  //   pair<int32_t, int32_t> range_pair
  //   double copy_number
  //   double copy_number_float
  //   double t_score
  //   double p_value
  range_t range;
  
  // Vectors used to store range_t variables,
  // used to merge ranges with the same copy number
  vector< range_t > range_t_vector;
  vector< range_t > range_t_vector_merged;
  bool add_to_merged_vector;
  
  ifstream tile_file;
  ifstream segment_file;
  ofstream out_file;
  ofstream final_file;
  
  tile_file.open ( tile_file_name.c_str() );
  segment_file.open ( segment_file_name.c_str() );
  
//@JEB: this output file is now redundant with GD output
//#define _FINAL_FILE_
#ifdef _FINAL_FILE_
  final_file.open( final_file_name.c_str() );
  final_file << "Start_Position" << "\t" <<
  "End_Position" << "\t" <<
  "Copy_Number" << endl;
#else
  (void) final_file_name;
#endif
  
//@JEB: this optional file is across the entire genome,
//      it could be used to plot the modeled coverage
//#define _SMOOTH_RANGES_FILE_
#ifdef _SMOOTH_RANGES_FILE_
  out_file.open ( out_file_name.c_str() );
  out_file << "Position\tSmooth_Coverage\n";
#else
  (void) out_file_name;
#endif
  
  /*
   * Populate the tile_data vector
   */
  tile_file >> skip >> skip >> tile_size;
  while ( tile_file >> tile_position >> tile_coverage )
  {
    tile_entry.first = tile_position;
    tile_entry.second = tile_coverage;
    tile_data.push_back( tile_entry );
  }
  
  
  /*
   * Populate data from ranges file into range_t_vector, calculating
   * copy number for each range along the way
   */
  getline(segment_file, skip);
  while (segment_file >> position_start >> position_end >> t_score_within >> p_value_within >> t_score_exists >> p_value_exists)
  {
    segment_sum = 0;
    num_tiles = 0;
    //get the mean of the segment in the tile file.
    //search through tile_data and find the segment.
    for ( uint32_t i = 0; i < tile_data.size(); i++ )
    {
      tile_entry = tile_data[i];
      
      /* ??? What and why: finding a tile that exists between the current
       segment.*/
      if ( position_start <= tile_entry.first &&
          position_end >= tile_entry.first )
      {
        segment_sum += tile_entry.second;
        num_tiles++;
      }
      /* ??? What and why: the for loop starts at the beginning of the genome.
       If the positions you're now looking at are past the segment you're
       filling in, you've gone past the segment. break.*/
      else if ( tile_entry.first > position_end )
      {
        break;
      }
    }
    segment_mean = segment_sum / num_tiles;
    
    //calculate the new value of the segment by taking the mean of the
    //segment divided by the mean of the tile data.
    copy_number_float = segment_mean / summary_average;
    
    //this rounding method doesn't handle negative values
    //fortunately, there's no such thing as negative coverage
    copy_number = floor(copy_number_float + .5);
    
    // Set values for range
    range.range_pair.first = position_start;
    range.range_pair.second = position_end;
    range.copy_number_float = copy_number_float;
    range.copy_number = copy_number;
    range.t_score_exists = t_score_exists;
    range.p_value_exists = p_value_exists;
    range.t_score_within = t_score_within;
    range.p_value_within = p_value_within;
    
    // Add range to end of vector
    range_t_vector.push_back(range);
  }
  
  // Debug output
  /*
   for (uint32_t i = 0; i < range_t_vector.size(); i++) {
   cout << "\nrange_t_vector[" << i << "]:" << endl;
   cout << "range_pair.first = " << range_t_vector[i].range_pair.first << endl;
   cout << "range_pair.second = " << range_t_vector[i].range_pair.second << endl;
   cout << "copy_number_float = " << range_t_vector[i].copy_number_float << endl;
   cout << "copy_number = " << range_t_vector[i].copy_number << endl;
   }
   cout << endl;
   */
  
  
  
  
  /**
   **
   **
   ** Print ranges to various output files
   **
   **
   **/
  
  // For each range within range_t_vector
  for (uint32_t i = 0; i < range_t_vector.size(); i++)
  {
    uint32_t position_start = range_t_vector[i].range_pair.first;
    uint32_t position_end = range_t_vector[i].range_pair.second;
    uint32_t copy_number = range_t_vector[i].copy_number;
    formatted_double copy_number_float(range_t_vector[i].copy_number_float, 2);
    double p_value = range_t_vector[i].p_value_exists;
    
    // @JEB create the genome diff evidence entry if copy_number is not one
    
    if (copy_number > 1) {
      cDiffEntry item(CN);
      item[SEQ_ID] = seq_id;
      item[START] = to_string<uint32_t>(position_start);
      item[END] = to_string<uint32_t>(position_end);
      item["tile_size"] = to_string<double>(settings.copy_number_variation_tile_size);
      item["copy_number"] = to_string<double>(copy_number);
      item["relative_coverage"] = copy_number_float.to_string();
      item["p-value"] = to_string<double>(p_value);
      gd.add(item);
    }
    

    
#ifdef _FINAL_FILE_
    
    //This writes a file for each final region with lines of the form:
    // start \t end \t rounded coverage
    
    written_final = false;
    for ( uint32_t i = 0; i < tile_data.size(); i++ )
    {
      tile_entry = tile_data[i];
      
      if ( position_start <= tile_entry.first &&
          position_end >= tile_entry.first )
      {
        tile_data[i].second = copy_number;
        
        if (! written_final)
        {
          final_file << position_start << "\t" <<
          position_end << "\t" <<
          tile_data[i].second << endl;
          written_final = true;
        }
      }
      
      else if ( tile_entry.first > position_end )
      {
        break;
      }
    }
#endif
  } // End while-loop

  
#ifdef _SMOOTH_RANGES_FILE_
  // write data to smoothed_ranges file
  for ( uint32_t i = 0; i < tile_data.size(); i++ )
  {
    tile_entry = tile_data[i];
    out_file << tile_entry.first << "\t" << tile_entry.second << endl;
  }
#endif
  
  gd.write(gd_file_name);
  
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
    read_out << gc_percentage_string(it->m_fasta_sequence.m_sequence) << endl;
  }
  
  alignment_list alignments;
  while (bam.read_alignments(alignments, false)) {
    double gc(0);

    for(alignment_list::iterator it=alignments.begin(); it!=alignments.end(); it++)
    {
      bam_alignment& a = *(it->get());
      
      int64_t start = (a.strand()==-1) ? a.reference_end_1() - _read_length + 1 : a.reference_start_1();
      string read_string = ref_seq_info[a.reference_target_id()].get_circular_sequence_1(start, _read_length);
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
      string read_string = it->get_circular_sequence_1(i, _read_length);
      ref_out << gc_percentage_string(read_string) << endl;
    }
  }
}
  

} // namespace breseq
