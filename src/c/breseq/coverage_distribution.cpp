
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
    summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff = 0.0;
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
                                    deletion_propagation_pr_cutoff
                                    );

    // First two lines are negative binomial parameters.
    // Next three lines are average, standard deviation, and index of overdispersion

    // Put these into summary
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
                              uint32_t tile_size
                              )
{
    //The index is the marker of the current tile.
    //It is always a multiple of tile_size
    uint32_t tile_index;

    //these are for taking running sums and averaging.
    double tile_sum;
    double average;

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

    in_file.open ( in_file_name.c_str() );
    out_file.open ( out_file_name.c_str() );

    //skipping the first line in the input file.
    getline(in_file, skip);

    out_file << "position\tcoverage\t" << tile_size << "\n";

    out_file << "0\t0\n";

    tile_index = 0;
    tile_sum = 0;

    do_write = true;

    while ( in_file >> coverage )
    {
        sum_coverage = coverage;

        in_file >> coverage;
        sum_coverage += coverage;

        //next 4 are redundant values;
        in_file >> redundant1 >> redundant2 >> redundant3 >> redundant4;

        //the next value in the line isn't used.
        in_file >> skip;

        //last is position (@DTF: begins at 1, not 0)
        in_file >> position;

        coverage = sum_coverage;

        //Check if the current position is in the current tile.
        //If it is, add the coverage to the sum and increment the count.
        //If it isn't, write the last tile and set up the variables for the new 
        //tile
        if ( position <= tile_index + tile_size ) {
            tile_sum += coverage;
        }
        
        else {
            
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

        if ( ignore_redundant_coverage ) {
            
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

    // @DTF: If the last position in the file also happens to be the last
    // position in the current tile (i.e. a multiple of tile size), and
    // if this is the last tile in the file
    if ( position % tile_size == 0 && position != last_position )
    {
        average = tile_sum / tile_size;
      
        out_file << position << "\t" << average << "\n";
    }

    in_file.close();
    out_file.close();
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
   
        If the answer to step 2 is YES/ACCEPT, you break the search segment in to three new subsegments or ranges: the middle subsegment (which is the newly found outlier), the left subsegment (everything to the left of the middle), and the right subsegment (everything to the right of the middle). You then repeat steps 1 and 2 for each of these new subsegments, treating each as its own entity with its own potential outliers to be found. We're looking for "bumbs within bumps".
 
        If the answer to step 2 is NO/REJECT, you output the search segment to the ranges.tab file. If no outliers were found within the entire genome (i.e. after the first run of the loop), of course the "search segment" outputted will span the entire genome. In every other case, the search segment outputted to ranges.tab is terminal, meaning that it doesn't have any other outliers within it (although it may itself be outlier).
 
 * T_Score Formula: Please read "Details, Circular Binary Segmentation" above before going any further. The formula for calculating the t_score, or "test statistic", comes this academic paper: <Biostat (2004) 5 (4): 557-572. doi: 10.1093/biostatistics/kxh008>. Fully understanding the formula isn't necessary to work on the code. If you wish to do so anyway, you should refer the paper, and we borrow some of the same variable names for ease of identification in that respect (e.g. "y", "z", "d1", "d2"). All you need to know is that the simplest steps of the formula determine the avg coverage of the tiles between i and j, and the avg coverage of the tiles outside of i and j, then subtract one from the other to find the difference. Further steps, for example, factor in the width of the potential outlier (between i and j) compared to the width of the entire search segment. The end result is a single number that takes into account all relevent measures of size. Here is the formula:
 
         y = (sum_of_coverage_to_j_pos - sum_of_coverage_to_i_pos) / (j_position - i_position);
         z = (sum_of_coverage_in_segment - sum_of_coverage_to_j_pos + sum_of_coverage_to_i_pos) / (num_positions_in_segment - j_position + i_position);
         d1 = (1.0 / (j_position - i_position));
         d2 = (1.0 / (num_positions_in_segment - j_position + i_position));
         t = (y - z) / sqrt( d1 + d2 );
 
 * Tiles vs Position: Please read the top of the CoverageDistribution::tile function before going any further. This function keeps track of location within a segment in two ways: by "position" and by "tile". Each entry in our input file, tiled.tab, represents the average coverage of every 500 positions across the entire genome (or whatever tile_size is). So of course each POSITION within tiled.tab will be a multiple of tile_size (i.e. 0, 500, 1000, 1500, etc). Say we're looking in tiled.tab at position 30500 which shows coverage of 60.2. This means that the avg coverage is 60.2 for the previous 500 positions (inclusive), so from 30001 to 30500. Each of these 500 position sections is considered a TILE, so there are as many tiles as entries in tiled.tab. Although not explicitly labeled in the file, in the code we refer to each tile with a number in sequence, where the first entry in tiled.tab represents tile 0, the next is tile 1, then tile 2, tile 3, and so on. We do this to make it easier to manipulate the entries later on, like loading them into data structures and iterating through them with loops. NOTE: The first entry in every tiled.tab is a pair of zeros, meaning that at position 0, avg coverage is 0. Although there technically is no position 0, because the coverage data in coverage.tab starts at position 1, we make one up because it allows us to eliminate special cases during for-loops throughout the code. (Coverage.tab is the file from which tiled.tab is generated.)
 
 * Single-Tile Outliers: The algorithm will have trouble finding outliers that are only one tile wide if they are isolated and not connected to other outlying segments. It is much better at finding these single-tile outliers if they border the edge of another larger outlier that is two or more tiles wide (for an example, look at the HTML output for genome REL10951). This trade-off is intentional because there isn't a good way to find a single-tile outlier: the noise we are trying to screen-out often occurs in spikes that are only one tile wide. The chances of finding something meaningful increase significantly when you require the outlier to be at least two tiles wide because, statistically, it's a much rarer occurrence.
 
 * Floating Point Comparisons (and the inherent problems within): Using comparison operators (i.e. "<", ">", "==", "!=", etc) with floating point numbers like doubles can cause problems because floating point math is not exact. Here is a good url to read about it: <bit.ly/xIW2Mf>. To quote: "Many values cannot be precisely represented using binary floating point numbers, and the limited precision of floating point numbers means that slight changes in the order of operations or the precision of intermediates can change the result." For example, storing the same number, even an integer, within two different floats and comparing to see if they are equal will often return false. This function uses less-than and greater-than floating point comparisons, but the difference between the compared values is usually large enough that it yields the expected results anyway. However, at some point it might be a good idea to follow some of the suggestions in the referenced url and replace the current method of comparing exact values with one that compares ranges, like ULP or FLT_EPSILON.
 */
void CoverageDistribution::find_segments(const Settings& settings,
                                         double summary_average,       // Avg coverage depth
                                                                       // across entire genome.
                                         string in_file_name,          // tiled.tab
                                         string out_file_name,         // ranges.tab
                                         string history_file_name)     // history.tab
                                         
{

    (void)settings;

    // A search segment is represented by a pair of tiles: the first tile of the
    // segment and the last tile of the segment (inclusive)
    pair<int32_t, int32_t> current_search_pair_by_tile;
    pair<int32_t, int32_t> next_search_pair_by_tile;
    
    // A vector keeps track of every search segment that will be examined for outliers
    vector< pair<uint32_t, uint32_t> > all_search_pairs_vector;
    
    // Used to populate data structures needed by the CBS algorithm
    int32_t current_position;
    int32_t previous_position;
    
    // Relative positions of i and j within the search segment
    int32_t i_position;
    int32_t j_position;
    
    // Used to determine absolute position within the entire genome
    int32_t position_hash;
    
    // Defines the beginning and end tiles of a potential outlier
    int32_t best_i_tile;
    int32_t best_j_tile;
    
    // Defines the beginning and end tiles of the largest t_score within
    // a random segment. Used only for stepping thru the code in debug mode.
    int32_t r_best_i_tile;
    int32_t r_best_j_tile;

    // Keeps track of the total number of positions in the search segment,
    // i.e. the width of the search segment in bp.
    int32_t num_positions_in_segment;
    
    // Maps the relative sum of coverage up to every position within the search segment
    map<int32_t, double> sum_of_coverage_to_position;
    double running_sum = 0; // Used to populate above
    
    // Relative sums of coverage up to i_position and j_position within the search segment
    double sum_of_coverage_to_i_pos;
    double sum_of_coverage_to_j_pos;
    
    // Relative sum of coverage within the entire search segment
    double sum_of_coverage_in_segment;

    // For keeping track of the largest t_score within the search segment
    double best_t_score;
    int32_t t_length;
    int32_t best_t_length;
    
    // For use in the t_score formula (see below):
    double y; // Avg coverage inside of i and j 
    double z; // Avg coverage outside of i and j
    double d1; // To factor in the width of i and j
    double d2; // To factor in the width of i and j
    double t_score;
    
    /*
     t_score formula:
     
         y = (sum_of_coverage_to_j_pos - sum_of_coverage_to_i_pos) / (j_position - i_position);
         z = (sum_of_coverage_in_segment - sum_of_coverage_to_j_pos + sum_of_coverage_to_i_pos) / (num_positions_in_segment - j_position + i_position);
         d1 = (1.0 / (j_position - i_position));
         d2 = (1.0 / (num_positions_in_segment - j_position + i_position));
         t_score = (y - z) / sqrt( d1 + d2 );
    */

    // For keeping track of the best t_score among the random segments
    double r_t_score;
    double r_best_t_score;
    
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
    int32_t best_t_was_bigger_counter;
    int32_t max_error_threshold;

    // Used for randomized testing for statistical significance
    int32_t randomizations = 100; // How many times the random test loop will run
    max_error_threshold = int(randomizations * .05) + 1; // Ex: 6
    srand(0);
    
    // Define filestreams
    ifstream tiles_file;
    ofstream ranges_file;
    ofstream history_file;
    
    // Open filestreams
    tiles_file.open ( in_file_name.c_str() );
    ranges_file.open ( out_file_name.c_str() );
    history_file.open ( history_file_name.c_str() );

    // Print labels for all the columns in the output files
    ranges_file << "Start_Position\tEnd_Position\tT_Score\tP_Value\n";

    history_file << "Start_Search\tEnd_Search\tStart_Position\tEnd_Position\t"
                 << "Start_Segment\tEnd_Segment\t"
                 << "T_Score\tMean_Coverage_Inside\tMean_Coverage_Left\tMean_Coverage_Right\t"
                 << "CN_Inside\tCN_Left\tCN_Right\t"
                 << "P_Value_Inside\n";

    
    

    
    /*
     * Get data from tiled.tab file
     */
    
    // For reading input from tiled.tab
    uint32_t num_tiles = 0;
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
    num_positions_in_segment = position;
    
    // Add the first search segment (the entire genome) to the end of all_search_pairs_vector
    current_search_pair_by_tile.first = 1;
    current_search_pair_by_tile.second = num_tiles - 1;
    all_search_pairs_vector.push_back(current_search_pair_by_tile);
    
    
    
    
    
    /*
     * Beginning of CBT algorithm loop
     */
    
    // Until there are no more search segments in which to look for outliers...
    while ( all_search_pairs_vector.size() > 0 ) {
        
        /* Set the current search segment: */
        
        // Remove the last search pair from the all_search_pairs_vector and
        // set it as the current search pair. NOTE: using pop_back() is why
        // results are printed to the ranges and history files in seemingly
        // reverse order from which they are added to the queue (later down the line).
        current_search_pair_by_tile.first = all_search_pairs_vector.back().first;
        current_search_pair_by_tile.second = all_search_pairs_vector.back().second;
        all_search_pairs_vector.pop_back();
        
        // If the current search segment has 0 or 1 tiles, it cannot be examined for
        // outliers because no outliers can exist within it.
        if (current_search_pair_by_tile.first == current_search_pair_by_tile.second ||
            current_search_pair_by_tile.first + 1 == current_search_pair_by_tile.second)
        {
            ranges_file << ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size + 1 << "\t"
                     << ordered_positions_in_segment[current_search_pair_by_tile.second] << "\t"
                     << 0.0 << "\t"
                     << 1.0 << endl;

            history_file << ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size + 1 << "\t"
                         << ordered_positions_in_segment[current_search_pair_by_tile.second] << "\t"
                         << ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size + 1 << "\t"
                         << ordered_positions_in_segment[current_search_pair_by_tile.second] << "\t"
                         << -1 << "\t"
                         << -1 << "\t"
                         << 0.0 << "\t"
                         << -1 << "\t"
                         << -1 << "\t"
                         << -1 << "\t"
                         << -1 << "\t"
                         << -1 << "\t"
                         << -1 << "\t"
                         << 1.0 << endl;
            continue;
        }
      
        
        
        
        
        /*  
         * Prep work to populate data structures needed by CBS algorithm
         */
        
        // VERBOSE, set to TRUE for debugging output
        const bool VERBOSE = false;
        
        if (VERBOSE) {
            cout << "\n############## PREP_LOOP_VERBOSE OUTPUT ###############" << endl;
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
        
        // Set to -1 so the first t_score found is always greater than the initialized value
        best_t_score = -1;
        
        // Zero out values
        best_t_was_bigger_counter = 0;
        best_t_length = 0;
        
        
        
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
                if ( abs(t_score) > best_t_score ||
                    ( abs(t_score) == best_t_score && t_length > best_t_length) )
                {
                    best_t_score = abs(t_score);
                    best_i_tile = i_tile;
                    best_j_tile = j_tile;
                    best_t_length = t_length;
                }
            
            } // End of j_tile for-loop
        } // End of i_tile for-loop
        
        // Save info about the potential outlier that has just been found
        i_position = ordered_positions_in_segment[best_i_tile] - position_hash;
        j_position = ordered_positions_in_segment[best_j_tile] - position_hash;
        sum_of_coverage_to_i_pos = sum_of_coverage_to_position[ i_position ];
        sum_of_coverage_to_j_pos = sum_of_coverage_to_position[ j_position ];

        

        
        
        /*
         * Make debugging info for the history file
         */
        
        // Calculate the mean coverage between i and j
        mean_inside = (sum_of_coverage_to_j_pos - sum_of_coverage_to_i_pos) / double(best_j_tile - best_i_tile);
        
        // Calculate the copy number between i and j by dividing mean_inside by the avg
        // coverage across the entire genome, and then rounding to nearest whole value
        copy_inside = static_cast<int32_t>(floor((mean_inside / summary_average) + .5));
        
        // If the left edge of the potential outlier equals the left edge of the
        // search segment, there's no left subsegment in which to search for further outlers
        if (best_i_tile == current_search_pair_by_tile.first - 1) {
            mean_left = -1;
            copy_left = -1;
        }
        // Otherwise, the subsegment to the left will be added to the queue
        else {
            mean_left = sum_of_coverage_to_i_pos / double(best_i_tile);
            copy_left = static_cast<int32_t>(floor((mean_left / summary_average) + .5));
        }
        
        // If the right edge of the potential outlier equals the right edge of the
        // search segment, there's no right subsegment in which to search for further outlers
        if ( best_j_tile == current_search_pair_by_tile.second ) {
            mean_right = -1;
            copy_right = -1;
        }
        // Otherwise, the subsegment to the right will be added to the queue
        else {
            mean_right = (sum_of_coverage_in_segment - sum_of_coverage_to_j_pos) / double(current_search_pair_by_tile.second - best_j_tile);
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
            
            /* Below is an unsusual yet efficient way to randomly select tiles within the
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
            
            // Note: This line may again be redundant
            current_position = ordered_positions_in_segment
            [current_search_pair_by_tile.first  - 1] - position_hash;
            
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
            
            // Set to -1 so the first t_score found is always greater than the initialized value
            r_best_t_score = -1.0;

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
                    // Legacy code used for finding mean error
                    /*
                    inside_mean = (sum_of_coverage_to_j_pos - sum_of_coverage_to_i_pos) / (j_tile - i_tile);
                    outside_mean = (sum_of_coverage_in_segment - sum_of_coverage_to_j_pos + sum_of_coverage_to_i_pos) / (num_tiles - j_tile + i_tile);
                    mean_error = 0;
                    for ( int32_t k = 0; k < ordered_positions_in_segment.size(); k++ )
                    {
                        if ( i_tile < k && k <= j_tile )
                        {
                            //it's usually faster to use x * x than x ** 2.
                            mean_error += ((coverage_map[ordered_positions_in_segment[k]] - inside_mean) * (coverage_map[ordered_positions_in_segment[k]] - inside_mean));
                        }
                        else
                        {
                            mean_error += ((coverage_map[ordered_positions_in_segment[k]] - outside_mean) * (coverage_map[ordered_positions_in_segment[k]] - outside_mean));
                        }
                    }
                    mean_error = sqrt( mean_error );
                    */
                    
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
                    if ( r_t_score > r_best_t_score )
                    {
                        r_best_t_score = abs(r_t_score);
                        
                        // Used for stepping thru the code in debug mode
                        r_best_i_tile = i_tile;
                        r_best_j_tile = j_tile;
                    }

                } // End of for-loop j_tile
            } // End of for-loop i_tile
            
            //  Now we know the best t_score of the random segment
            if (VERBOSE) {
                cout << "\n\n#####################################";
                cout << "\n best_t_score: " << fixed << setprecision(40) << best_t_score;
                cout << "\r_best_t_score: " << fixed << setprecision(40) << r_best_t_score
                     << endl;
            }
            
            
            
            
            
            /*
             * Count the times best_t_score is larger than r_best_t_score
             */
            
            if ( best_t_score > r_best_t_score )
            {
                best_t_was_bigger_counter++;
            }
            
            // This heuristic breaks out of the randomized testing early if it becomes
            // impossible to reach the number of counts required to be significant. 
            // Ex: if ([(96 + 1) - 90] > 6) { ...
            // NOTE: This could probably be simplified to the following:
            /*  if (examinations - best_t_was_bigger_counter >= max_error_threshold){ ... */
            if ((examinations + 1) - best_t_was_bigger_counter > max_error_threshold){

                // Break the entire randomized testing loop
                break;
            }
       
        } // End of randomized testing (aka "examinations") for-loop
        
        
      
        
        
        /*
         * If answer to step 2 is YES/ACCEPT (the potential outlier is significant)
         */
        
        // If the potential outlier is significant, and if the potential outlier is NOT
        // just the entire genome...
        if ( (double)(best_t_was_bigger_counter) / randomizations >= .95 &&
           !((current_search_pair_by_tile.first == best_i_tile + 1) && (current_search_pair_by_tile.second == best_j_tile))
         )
        {
            
            
            
            /*
             * Break search segment into three subsegments and add each to search queue
             */
            
            // Add the middle subsegment (which is the newly found outlier between i and j)
            // to the search queue (i.e. all_search_pairs_vector)
            next_search_pair_by_tile.first = best_i_tile + 1;
            next_search_pair_by_tile.second = best_j_tile;
            all_search_pairs_vector.push_back( next_search_pair_by_tile );

            // Record this within the history file
            history_file << ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size + 1 << "\t" <<
                          ordered_positions_in_segment[current_search_pair_by_tile.second] << "\t" <<
                          ordered_positions_in_segment[best_i_tile + 1] - tile_size + 1 << "\t" <<
                          ordered_positions_in_segment[best_j_tile] << "\t" <<
                          ordered_positions_in_segment[next_search_pair_by_tile.first] - tile_size + 1 << "\t" <<
                          ordered_positions_in_segment[next_search_pair_by_tile.second] << "\t" <<
                          best_t_score << "\t" <<
                          mean_inside << "\t" <<
                          mean_left << "\t" <<
                          mean_right << "\t" <<
                          copy_inside << "\t" <<
                          copy_left << "\t" <<
                          copy_right << "\t" <<
                          1 - (double)(best_t_was_bigger_counter) / randomizations << endl;

            //if best_i_tile doesn't include the first element, add the segment before
            //best_i_tile
            /* ??? checking to make sure the segment doesn't
             border the edge. if it borders the edge,
             there's nothing to add.*/
            
            // Add the left subsegment (if there is one, i.e. if the middle subsegment
            // doesn't border the start of the search segment) to the search queue
            if ( best_i_tile + 1 > current_search_pair_by_tile.first )
            {
                next_search_pair_by_tile.first = current_search_pair_by_tile.first;
                next_search_pair_by_tile.second = best_i_tile;
                all_search_pairs_vector.push_back( next_search_pair_by_tile );
            }

            // Add the right subsegment (if there is one, i.e. if the middle subsegment
            // doesn't border the end of the search segment) to the search queue
            if ( best_j_tile < current_search_pair_by_tile.second )
            {
                next_search_pair_by_tile.first = best_j_tile + 1;
                next_search_pair_by_tile.second = current_search_pair_by_tile.second;
                all_search_pairs_vector.push_back( next_search_pair_by_tile );
            }
          
        } // End of if-conditional YES/ACCEPT
        
        
        
        
        
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
            ranges_file << ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size + 1 << "\t"
                     << ordered_positions_in_segment[current_search_pair_by_tile.second] << "\t"
                     << best_t_score << "\t"
                     << 1 - (double)(best_t_was_bigger_counter) / randomizations << endl;
            
            history_file << ordered_positions_in_segment[current_search_pair_by_tile.first] - tile_size + 1 << "\t" <<
                            ordered_positions_in_segment[current_search_pair_by_tile.second] << "\t" <<
                            ordered_positions_in_segment[best_i_tile + 1] - tile_size + 1 << "\t" <<
                            ordered_positions_in_segment[best_j_tile] << "\t" <<
                            -1 << "\t" <<
                            -1 << "\t" <<
                            best_t_score << "\t" <<
                            mean_inside << "\t" <<
                            mean_left << "\t" <<
                            mean_right << "\t" <<
                            copy_inside << "\t" <<
                            copy_left << "\t" <<
                            copy_right << "\t" <<
                            1 - (double)(best_t_was_bigger_counter) / randomizations << endl;
        }
      
      
    } // End of CBT algorithm while-loop

    
    // If the search queue is empty, close the filestreams
    tiles_file.close();
    ranges_file.close();
}
  
    
    
/*
 * Function: smooth_segments
 * --------------------------------
 * Reads in [genome].ranges.tab and generates [genome].smoothed_ranges.tab,
 * [genome].cnv_final.tab, and [genome].cn_evidence.tab.
 
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
    cGenomeDiff gd; /* ???: cGenomeDiff is a data type defined by breseq. it holds
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
    double z_score; /* ???: synonymous with t_score. */
    double greater_than; /* ???: reads in the fourth value from segment_file.
                          The p-value of a segment or (1 - greater_than / randomizations). */


    //holds the value a segment's coverage will be changed to.
    double new_segment_mean;

    //used for easily accessing the tile file.
    pair<uint32_t, double> tile_entry;
    vector< pair<uint32_t, double> > tile_data;

    //initial parameter.
    uint32_t tile_size;

    //used for skipping lines.
    string skip;

    bool written_final;

    ifstream tile_file;
    ifstream segment_file;
    ofstream out_file;
    ofstream final_file;

    tile_file.open ( tile_file_name.c_str() );
    segment_file.open ( segment_file_name.c_str() );
    out_file.open ( out_file_name.c_str() );
    final_file.open( final_file_name.c_str() );

    out_file << "Position\tSmooth_Coverage\n";

    final_file << "Start_Position" << "\t" <<
                  "End_Position" << "\t" <<
                  "T_Score" << "\t" <<          // @DTF: possible typo here, was Z_Score
                  "Greater_Than" << "\t" <<
                  "Copy_Number" << endl;

    tile_file >> skip >> skip >> tile_size;

    //populate the tile_data vector

    while ( tile_file >> tile_position >> tile_coverage )
    {
      tile_entry.first = tile_position;
      tile_entry.second = tile_coverage;
      tile_data.push_back( tile_entry );
    }


    //go through each entry in the segment file.
    getline(segment_file, skip);
    while ( segment_file >> position_start >> position_end >> z_score >> greater_than )
    {
      segment_sum = 0;
      num_tiles = 0;
      written_final = false;
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
      new_segment_mean = segment_mean / summary_average;
      
      //this rounding method doesn't handle negative values
      //fortunately, there's no such thing as negative coverage
      new_segment_mean = floor(new_segment_mean + .5);
      
        // @JEB create the genome diff evidence entry if mean is not one
        /* ??? What does this block do: JEB are Barrick's initials. He's doing
         all of the cGenomeDiff coding so if you're really interested, ask him.
         Otherwise, it's not important to the algorithm.*/
      
      if (new_segment_mean > 1) {
        cDiffEntry item(CN);
        item[SEQ_ID] = seq_id;
        item[START] = to_string<uint32_t>(position_start);
        item[END] = to_string<uint32_t>(position_end);
        item["tile_size"] = to_string<double>(settings.copy_number_variation_tile_size);
        item["copy_number"] = to_string<double>(new_segment_mean);
        stringstream num;
        num << fixed << setprecision(2) << (segment_mean / summary_average);
        item["relative_coverage"] = num.str();
        gd.add(item);
      }
      
      
        //the following part of the code could probably be done much more elegantly.
        //i'm not entirely sure but, what i think it's doing is actually writing
        //a rounded coverage value for the single tile currently being examined.
        
        //set the new value of the coverage for this segment to this value.
        /* ??? Clarify: new value being the rounded coverage value, as in 2. */
      for ( uint32_t i = 0; i < tile_data.size(); i++ )
      {
        tile_entry = tile_data[i];
        
        if ( position_start <= tile_entry.first &&
             position_end >= tile_entry.first )
        {
          tile_data[i].second = new_segment_mean;
          
          if (! written_final)
          {
            final_file << position_start << "\t" <<
                          position_end << "\t" <<
                          z_score << "\t" <<
                          greater_than << "\t" <<
                          tile_data[i].second << endl;
            written_final = true;
          }
        }
        
        else if ( tile_entry.first > position_end )
        {
          break;
        }
      }
      
      
    }

    //write data to smoothed_ranges file
    for ( uint32_t i = 0; i < tile_data.size(); i++ )
    {
      tile_entry = tile_data[i];
      out_file << tile_entry.first << "\t" << tile_entry.second << endl;
    }

    tile_file.close();
    segment_file.close();
    out_file.close();
    final_file.close();

    // merge intervals
    /* ??? Clarify: more cGenomeDiff code. again ask Barrick if you really
     want to know. but unnecessary for the actual algorithm. */
    //that applies to everything below this line.
    gd.sort();
    diff_entry_list_t muts = gd.list();
    diff_entry_ptr_t last_de(NULL);
    cGenomeDiff gd_merged;

    for(diff_entry_list_t::iterator it = muts.begin(); it != muts.end(); it++)
    {
      //cout << **it << endl;
      
      if (!last_de.get())
        last_de = *it;
      else {
        diff_entry_ptr_t this_de = *it;
        
        if (  ( from_string<uint32_t>((*this_de)[START]) - 1 == from_string<uint32_t>((*last_de)[END]))
          &&  ( (*this_de)["copy_number"] == (*last_de)["copy_number"] )
          &&  ( (*this_de)[SEQ_ID] == (*last_de)[SEQ_ID] ) ) {
          
          (*last_de)[END] = (*this_de)[END];
          
          double length_last = from_string<uint32_t>((*last_de)[END]) - from_string<uint32_t>((*last_de)[START]) + 1;
          double length_this = from_string<uint32_t>((*this_de)[END]) - from_string<uint32_t>((*this_de)[START]) + 1;
          double rel_cov_last = from_string<double>((*last_de)["relative_coverage"]);
          double rel_cov_this = from_string<double>((*this_de)["relative_coverage"]);
          double new_rel_cov = (rel_cov_last * length_last + rel_cov_this * length_this) / (length_last + length_this);
          
          stringstream num;
          num << fixed << setprecision(2) << (new_rel_cov);
          (*last_de)["relative_coverage"] = num.str();
        }
        else {
          gd_merged.add(*last_de);
          last_de = this_de;
        }
      }
    }

    if (last_de.get())
      gd_merged.add(*last_de);


    gd_merged.write(gd_file_name);
}


/*
 * Function: calculate_periodicity
 * --------------------------------
 * SUMMARY GOES HERE
 */
void CoverageDistribution::calculate_periodicity (
                                string coverage_file_name,
                                string period_file_name,
                                uint32_t method,
                                uint32_t start_range,
                                uint32_t end_range,
                                uint32_t step
                                ){

//used for indexing the coverage_file's values.
stringstream line_stream;
string line;

string column_title;
uint8_t current_column;
uint8_t column_count;

uint8_t unique_top_column;
uint8_t unique_bot_column;
uint8_t redundant_top_column;
uint8_t redundant_bot_column;

int32_t unique_top;
int32_t unique_bot;
int32_t redundant_top;
int32_t redundant_bot;

//used to ignore values.
string skip;

ifstream coverage_file;
ofstream period_file;

//used for holding file data
vector<int32_t> unique_top_vector;
vector<int32_t> unique_bot_vector;
vector<int32_t> redundant_top_vector;
vector<int32_t> redundant_bot_vector;

//used to hold values for operations
vector<int32_t> offset_range;
vector<int32_t> top_vector;
vector<int32_t> bot_vector;
vector<int32_t> offset_bot_vector;
int32_t bot_value;
int32_t offset_sum;

//reserve vector space
unique_top_vector.reserve(4700000);
unique_bot_vector.reserve(4700000);
redundant_top_vector.reserve(4700000);
redundant_bot_vector.reserve(4700000);

coverage_file.open( coverage_file_name.c_str() );
period_file.open( period_file_name.c_str() );

//populate offset_range
for (uint32_t i = start_range; i <= end_range; i += step){
  offset_range.push_back(i);
}

//parse header line so that each line can be correctly indexed.
getline(coverage_file, line);
line_stream.str(line);

column_count = 0;
while (getline(line_stream, column_title, '\t')){
  if (column_title == "unique_top_cov"){
    unique_top_column = column_count;
  }
  else if (column_title == "unique_bot_cov"){
    unique_bot_column = column_count;
  }
  else if (column_title == "redundant_top_cov"){
    redundant_top_column = column_count;
  }
  else if (column_title == "redundant_bot_cov"){
    redundant_bot_column = column_count;
  }
  
  column_count++;
}

//read the entire file and store it in the 4 coverage vectors
line_stream.clear();
while ( getline(coverage_file, line) ){
  
  if (line == ""){
    break;
  }
  
  line_stream.clear();
  line_stream.str(line);
  
  for (current_column = 0; current_column < column_count; current_column++){
    if (current_column == unique_top_column){
      line_stream >> unique_top;
    }
    else if (current_column == unique_bot_column){
      line_stream >> unique_bot;
    }
    else if (current_column == redundant_top_column){
      line_stream >> redundant_top;
    }
    else if (current_column == redundant_bot_column){
      line_stream >> redundant_bot;
    }
    else {
      line_stream >> skip;
    }
  }
  
  unique_top_vector.push_back(unique_top);
  unique_bot_vector.push_back(unique_bot);
  redundant_top_vector.push_back(redundant_top);
  redundant_bot_vector.push_back(redundant_bot);
  
}

//set up template top and bot vectors for calculations
if ( method == 1 ){
  top_vector = unique_top_vector;
  bot_vector = unique_bot_vector;
}
else if ( method == 2 ){
  top_vector = unique_top_vector;
  bot_vector = unique_top_vector;
  
  for ( uint32_t i = 0; i < unique_top_vector.size(); i++ ){
    top_vector[i] += redundant_top_vector[i];
    bot_vector[i] += redundant_bot_vector[i];
  }
}

//go through each offset and find the offset sum
for ( uint32_t i = 0; i < offset_range.size(); i++ ){
  //cout << offset_range[i] << "\r";
  //cout.flush();
  
  //set up bottom vector for this offset
  offset_bot_vector = bot_vector;
  
  for ( int32_t j = 0; j < offset_range[i]; j++ ){
    bot_value = offset_bot_vector[0];
    offset_bot_vector.erase(offset_bot_vector.begin());
    offset_bot_vector.push_back( bot_value );
  }
  
  //finding sum
  
  offset_sum = 0;
  
  for (uint32_t j = 0; j < top_vector.size(); j++ ){
    //ignoring spots where there is 0 coverage in either the top or bottom
    if (top_vector[j] == 0 || offset_bot_vector[j] == 0 ) continue;
    offset_sum += ((top_vector[j] - offset_bot_vector[j]) * (top_vector[j] - offset_bot_vector[j]));
  }
  
  period_file << offset_range[i] << "\t" << offset_sum << endl;
}
//cout << endl;
coverage_file.close();
period_file.close();
}
  
} // namespace breseq
