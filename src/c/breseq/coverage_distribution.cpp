
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
      
      //last is position
      in_file >> position;
      
      coverage = sum_coverage;
      
      //Check if the current position is in the current tile.
      //If it is, add the coverage to the sum and increment the count.
      //If it isn't, write the last tile and set up the variables for the new 
      //tile
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
    
    if ( position % tile_size == 0 && position != last_position )
    {
      average = tile_sum / tile_size;
      
      out_file << position << "\t" << average << "\n";
    }
    
    in_file.close();
    out_file.close();
  }
  
  void CoverageDistribution::find_segments(const Settings& settings,
                                           double summary_average,
                                          string in_file_name,
                                          string out_file_name,
                                          string history_file_name
                                          )
  {
    (void)settings;
    //a search entry represents a segment of positions in the file that will be
    //examined. start and end are inclusive.
    //search entries are referencing the index in the file, not position
    pair<int32_t, int32_t> next_search;
    pair<int32_t, int32_t> current_search;
    
    //this holds all search entries
    vector< pair<uint32_t, uint32_t> > searching;
    
    //used to populate saved_sums
    double running_sum;
    
    //a saved_sum entry represents the sum of all coverage up to a position.
    map<int32_t, double> saved_sums;
    
    //ordered_sums holds the order that saved_sums was populated.
    vector<int32_t> ordered_sums;
    
    //position and coverage are taken from the file.
    int32_t position;
    double coverage;
    map<int32_t, double> all_coverage;
    
    //p_i and p_j represent the position of the indeces of i and j in the later
    //nested for-loops
    int32_t p_i;
    int32_t p_j;
    int32_t position_downset;
    int32_t check_position;
    int32_t previous_position;
    
    //best_i and best_j define the segments that scored the best in an iteration
    //these represent indeces, not positions.
    int32_t best_i;
    int32_t best_j;
    int32_t r_best_i;
    int32_t r_best_j;
    
    //number of entries in the input file
    uint32_t count;
    
    //highest position value in the input file
    int32_t m;
    
    double Si;
    double Sj;
    double Sm;
    
    //the t value as per the CBS algorithm is the score of a given segment.
    double t;
    double t_left;
    double t_right;
    double best_t;
    double best_t_left;
    double best_t_right;
    int32_t t_length;
    int32_t best_t_length;
    
    //to easily see the t calculation these are used:
    double y;
    double z;
    double d1;
    double d2;
    
    int32_t randomizations;
    
    double r_t_left;
    double r_t_right;
    double r_t;
    double r_best_left;
    double r_best_right;
    double r_best_t;
    
    //extra useful data written
    double mean_inside;
    double mean_left;
    double mean_right;
    int32_t copy_inside;
    int32_t copy_left;
    int32_t copy_right;
    
    //used for finding error
    //double inside_mean;
    //double outside_mean;
    //double mean_error;
    
    //used for randomizing a segment
    vector<double> randomized_coverage;
    list<int32_t>::iterator available_coverage_it;
    list<int32_t> available_coverage;
    int32_t random_index;
    uint32_t best_greater_than;
    uint32_t left_greater_than;
    uint32_t right_greater_than;
    
    //to skip the header
    string skip;
    int32_t tile_size;
    
    ifstream in_file;
    ofstream out_file;
    ofstream history_file;
    
    randomizations = 100;
    
    srand(0);
    
    in_file.open ( in_file_name.c_str() );
    out_file.open ( out_file_name.c_str() );
    history_file.open ( history_file_name.c_str() );
    
    in_file >> skip >> skip >> tile_size;
    
    out_file << "Start_Position\tEnd_Position\tT_Score\tP_Value\n";
    history_file << "Start_Search\tEnd_Search\tStart_Position\tEnd_Position\t"
                 << "Start_Segment\tEnd_Segment\t"
                 << "T_Score\tMean_Coverage_Inside\tMean_Coverage_Left\tMean_Coverage_Right\t"
                 << "CN_Inside\tCN_Left\tCN_Right\t"
                 << "P_Value_Inside\tP_Value_Left\tP_Value_Right\n";
    
    //First, populate saved_sums and ordered_positions.
    //and get length
    count = 0;
    running_sum = 0;
    
    while ( in_file >> position >> coverage )
    {
      all_coverage[position] = coverage;
      ordered_sums.push_back( position );
      count++;
    }
    
    //the value of position is now the highest position value.
    m = position;
    
    //The first search entry goes from 1 to count - 1
    //0 is only used as a place holder.
    current_search.first = 1;
    current_search.second = count - 1;
    searching.push_back ( current_search );
    
    //now calculate z for each segment inside of current_search
    //repeat until there are no more search entries.
    while ( searching.size() > 0 )
    {
      current_search.first = searching.back().first;
      current_search.second = searching.back().second;
      
      searching.pop_back();
      
      if (current_search.first == current_search.second ||
          current_search.first + 1 == current_search.second)
      {
        out_file << ordered_sums[current_search.first] - tile_size + 1 << "\t"
                 << ordered_sums[current_search.second] << "\t"
                 << 0.0 << "\t"
                 << 1.0 << endl;
        history_file << ordered_sums[current_search.first] - tile_size + 1 << "\t"
                     << ordered_sums[current_search.second] << "\t"
                     << ordered_sums[current_search.first] - tile_size + 1 << "\t"
                     << ordered_sums[current_search.second] << "\t"
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
      
      //cout << "search" << current_search.first << " " << current_search.second << endl;
      //cin.get();
      
      position_downset = ordered_sums[current_search.first - 1];
      
      saved_sums.clear();
      check_position = ordered_sums[current_search.first - 1] - position_downset;
      saved_sums[check_position] = 0;
      //cout << "initial";
      //cout << saved_sums[check_position] << " ";
      count = 0;
      
      for (int32_t i = current_search.first; i <= current_search.second; i ++)
      {
        check_position = ordered_sums[i] - position_downset;
        previous_position = ordered_sums[i - 1] - position_downset;
        
        saved_sums[check_position] = all_coverage[check_position + position_downset] + saved_sums[previous_position];
        //cout << saved_sums[check_position] << " ";
        count ++;
      }
      check_position = ordered_sums[current_search.second] - position_downset;
      Sm = saved_sums[check_position];
      m = check_position;
      //cout << endl;
      
      
      //cout << ordered_sums[current_search.first] << " " << ordered_sums[current_search.second] << " " << endl;
      
      
      best_t = -1;
      best_t_left = -1;
      best_t_right = -1;
      best_greater_than = 0;
      left_greater_than = 0;
      right_greater_than = 0;
      best_t_length = 0;
      
      
      //i and j are inclusive boundaries so <= is used in the comparison.
      
      //find the best_t
      for ( int32_t i = current_search.first - 1; i < current_search.second; i++ )
      {
        p_i = ordered_sums[i] - position_downset;
        Si = saved_sums[ p_i ];
        for ( int32_t j = i + 1; j <= current_search.second; j++ )
        {
          //get means of i to j and outside i to j.
          
          /*inside_mean = (Sj - Si) / (j - i);
          outside_mean = (Sm - Sj + Si) / (count - j + i);
          
          mean_error = 0;
          for ( int32_t k = 0; k < ordered_sums.size(); k++ )
          {
            if ( i < k && k <= j )
            {
              //it's usually faster to use x * x than x ** 2.
              mean_error += ((all_coverage[ordered_sums[k]] - inside_mean) *
                             (all_coverage[ordered_sums[k]] - inside_mean));
            }
            else
            {
              mean_error += ((all_coverage[ordered_sums[k]] - outside_mean) *
                             (all_coverage[ordered_sums[k]] - outside_mean));
            }
          }
          
          mean_error = sqrt( mean_error );*/
          
          
          p_j = ordered_sums[j] - position_downset;
          Sj = saved_sums[ p_j ];
          
          y = (Sj - Si) / (p_j - p_i);
          z = (Sm - Sj + Si) / (m - p_j + p_i);
          d1 = (1.0 / (p_j - p_i));
          d2 = (1.0 / (m - p_j + p_i));
          t = (y - z) / sqrt( d1 + d2 );
          
          t_length = j - i;
          
          if ( abs(t) > best_t ||
               abs(t) == best_t && t_length > best_t_length )
          {
            best_t = abs(t);
            best_i = i;
            best_j = j;
            best_t_length = t_length;
            
            y = Si / p_i;
            z = (Sm - Si) / ( m - p_i );
            d1 = 1.0 / p_i;
            d2 = 1.0 / (m - p_i);
            t_left = (y - z) / sqrt ( d1 + d2 );
            best_t_left = abs(t_left);
            
            y = Sj / p_j;
            z = (Sm - Sj) / ( m - p_j );
            d1 = 1.0 / p_j;
            d2 = 1.0 / (m - p_j);
            t_right = (y - z) / sqrt ( d1 + d2 );
            best_t_right = abs(t_right);
          }
        }
      }
      p_i = ordered_sums[best_i] - position_downset;
      Si = saved_sums[ p_i ];
      p_j = ordered_sums[best_j] - position_downset;
      Sj = saved_sums[ p_j ];
      
      //making useful information that will be written to the history file
      if (best_i == current_search.first - 1){
        mean_left = -1;
        copy_left = -1;
      }
      else{
        mean_left = Si / double(best_i);
        copy_left = floor((mean_left / summary_average) + .5);
      }
      
      mean_inside = (Sj - Si) / double(best_j - best_i);
      copy_inside = floor((mean_inside / summary_average) + .5);
      
      if ( best_j == current_search.second ){
        mean_right = -1;
        copy_right = -1;
      }
      else{
        mean_right = (Sm - Sj) / double(current_search.second - best_j);
        copy_right = floor((mean_right / summary_average) + .5);
      }
      
      //cout << "population" << Sm << " " << m << endl;
      //cout << "best_initial " << best_i << " " << best_j << endl;
      //cout << "searched through " << current_search.first << " " << current_search.second << endl;
      //cout << "best_pos" << ordered_sums[best_i] << " " << ordered_sums[best_j] << endl;
      //cout << "best_cov" << saved_sums[ordered_sums[best_i]] << " " << saved_sums[ordered_sums[best_j]] << endl;
      //cout << "t" << best_t << endl;
      //cout << "calc  " << best_y << " " << best_z << " " << best_d1 << " " << best_d2 << endl;
      //cin.get();
      //examine the segment i + 1 to j further by randomizing its distribution and scoring it.
      //randomized_segment_scores.empty();
      for ( int32_t examinations = 0; examinations < randomizations; examinations++ )
      {
        available_coverage.clear();
        randomized_coverage.clear();
        
        for ( int32_t i = current_search.first; i <= current_search.second; i++ )
        {
          available_coverage.push_back(i);
        }
        
        for ( int32_t i = current_search.first; i <= current_search.second; i++ )
        {
          random_index = rand() % available_coverage.size();
          
          //remove whatever random index was used from the available_coverage vector.
          for (available_coverage_it = available_coverage.begin(); available_coverage_it != available_coverage.end(); available_coverage_it++)
          {
            if ( random_index == 0 )
            {
              randomized_coverage.push_back(all_coverage[ordered_sums[*available_coverage_it]]);
              available_coverage.erase(available_coverage_it);
              break;
            }
            random_index--;
          }
          
        }
        
        //create the new saved sums
        //cout << "rand";
        check_position = ordered_sums[current_search.first  - 1] - position_downset;
        //cout << saved_sums[check_position] << " ";
        for (size_t j = 0; j < randomized_coverage.size(); j++)
        {
          check_position = ordered_sums[current_search.first + j] - position_downset;
          previous_position = ordered_sums[current_search.first + j - 1] - position_downset;
          saved_sums[check_position] = saved_sums[previous_position] + randomized_coverage[j];
          //cout << saved_sums[check_position] << " ";
        }
        //cout << endl;
        //randomized_segment now has a random segment of coverage keeping
        //the positions in place.
        
        r_best_t = -1;
        r_best_left = -1;
        r_best_right = -1;
        
        //find the t score of this new segment of coverage.
        for ( int32_t i = current_search.first - 1; i < current_search.second; i++ )
        {
          p_i = ordered_sums[i] - position_downset;
          Si = saved_sums[ p_i ];
          for ( int32_t j = i + 1; j <= current_search.second; j++ )
          {
            /*inside_mean = (Sj - Si) / (j - i);
            outside_mean = (Sm - Sj + Si) / (count - j + i);
            
            mean_error = 0;
            for ( int32_t k = 0; k < ordered_sums.size(); k++ )
            {
              if ( i < k && k <= j )
              {
                //it's usually faster to use x * x than x ** 2.
                mean_error += ((all_coverage[ordered_sums[k]] - inside_mean) *
                               (all_coverage[ordered_sums[k]] - inside_mean));
              }
              else
              {
                mean_error += ((all_coverage[ordered_sums[k]] - outside_mean) *
                               (all_coverage[ordered_sums[k]] - outside_mean));
              }
            }
            
            mean_error = sqrt( mean_error );*/
            
            p_j = ordered_sums[j] - position_downset;
            Sj = saved_sums[ p_j ];
            
            //calculating t score of middle segment.
            
            y = (Sj - Si) / (p_j - p_i);
            z = (Sm - Sj + Si) / (m - p_j + p_i);
            d1 = 1.0 / (p_j - p_i);
            d2 = 1.0 / (m - p_j + p_i);
            r_t = (y - z) / sqrt( d1 + d2 );
            r_t = abs(r_t);
            
            y = Si / p_i;
            z = (Sm - Si) / ( m - p_i );
            d1 = 1.0 / p_i;
            d2 = 1.0 / (m - p_i);
            t_left = (y - z) / sqrt ( d1 + d2 );
            r_t_left = abs(t_left);
            
            y = Sj / p_j;
            z = (Sm - Sj) / ( m - p_j );
            d1 = 1.0 / p_j;
            d2 = 1.0 / (m - p_j);
            t_right = (y - z) / sqrt ( d1 + d2 );
            r_t_right = abs(t_right);
            
            if ( r_t > r_best_t )
            {
              r_best_t = abs(r_t);
            }
            
            if ( r_t_left > r_best_left )
            {
              r_best_left = r_t_left;
              r_best_i = i;
            }
            
            if ( r_t_right > r_best_right )
            {
              r_best_right = r_t_right;
              r_best_j = j;
            }
            
          }
        }
        //cout << "b " << r_best_t << " " << best_t << " " << best_greater_than << endl;
        //cout << "rcacout << "lol formatting" << 1 - (double)(best_greater_than) / randomizations << endl;lc " << y << " " << z << " " << d1 << " " << d2 << endl;
        //cin.get();
        //r_best_* variables have the best of this random segment.
        
        if ( best_t > r_best_t )
        {
          best_greater_than++;
        }
        
        if ( best_t_left > r_best_left ){
          left_greater_than++;
        }
        
        if ( best_t_right > r_best_right ){
          right_greater_than++;
        }
      }
      //cout << "best" << best_greater_than << endl;
      //if 95% of the t scores of the randomized coverage are less than
      //the initial best_t, add a new search for this segment
      
      //if the range to be added is identical to the current search, don't add it
      //only write it.
      //cout << best_greater_than << endl;
      //cout << best_greater_than / 10.0 << endl;
      //cout << current_search.first << " " << best_i + 1 << endl;
      //cout << current_search.second << " " << best_j << endl;
      //cin.get();
      if ( (double)(best_greater_than) / randomizations >= .95 &&
           !((current_search.first == best_i + 1) && (current_search.second == best_j))
         )
      {
        if ( (double) (left_greater_than) / randomizations >= .95 &&
             (double) (right_greater_than) / randomizations >= .95){
          //cout << "add_segment\n";
          next_search.first = best_i + 1;
          next_search.second = best_j;
          searching.push_back( next_search );
          
          history_file << ordered_sums[current_search.first] - tile_size + 1 << "\t" <<
                          ordered_sums[current_search.second] << "\t" <<
                          ordered_sums[best_i + 1] - tile_size + 1 << "\t" <<
                          ordered_sums[best_j] << "\t" <<
                          ordered_sums[next_search.first] - tile_size + 1 << "\t" <<
                          ordered_sums[next_search.second] << "\t" <<
                          best_t << "\t" <<
                          mean_inside << "\t" <<
                          mean_left << "\t" <<
                          mean_right << "\t" <<
                          copy_inside << "\t" <<
                          copy_left << "\t" <<
                          copy_right << "\t" <<
                          1 - (double)(best_greater_than) / randomizations <<  "\t" <<
                          1 - (double)(left_greater_than) / randomizations <<  "\t" <<
                          1 - (double)(right_greater_than) / randomizations << endl;
          
          //if best_i doesn't include the first element, add the segment before
          //best_i
          if ( best_i + 1 > current_search.first )
          {
            
            next_search.first = current_search.first;
            next_search.second = best_i;
            searching.push_back( next_search );
            
          }
          
          //if best_j doesn't include the last element, add the segment after
          //best_j
          if ( best_j < current_search.second )
          {
            
            next_search.first = best_j + 1;
            next_search.second = current_search.second;
            searching.push_back( next_search );
            
          }
        }
        
        else if ((double) (left_greater_than) / randomizations >= .95){
          next_search.first = best_i + 1;
          next_search.second = current_search.second;
          searching.push_back( next_search );
          
          history_file << ordered_sums[current_search.first] - tile_size + 1 << "\t" <<
                          ordered_sums[current_search.second] << "\t" <<
                          ordered_sums[best_i + 1] - tile_size + 1 << "\t" <<
                          ordered_sums[best_j] << "\t" <<
                          ordered_sums[next_search.first] - tile_size + 1 << "\t" <<
                          ordered_sums[next_search.second] << "\t" <<
                          best_t << "\t" <<
                          mean_inside << "\t" <<
                          mean_left << "\t" <<
                          mean_right << "\t" <<
                          copy_inside << "\t" <<
                          copy_left << "\t" <<
                          copy_right << "\t" <<
                          1 - (double)(best_greater_than) / randomizations <<  "\t" <<
                          1 - (double)(left_greater_than) / randomizations <<  "\t" <<
                          1 - (double)(right_greater_than) / randomizations << endl;
          
          if ( best_i > current_search.first )
          {
            next_search.first = current_search.first;
            next_search.second = best_i;
            searching.push_back( next_search );
          }
        }
        
        else if ((double) (right_greater_than) / randomizations >= .95){
          next_search.first = current_search.first;
          next_search.second = best_j;
          searching.push_back( next_search );
          
          history_file << ordered_sums[current_search.first] - tile_size + 1 << "\t" <<
                          ordered_sums[current_search.second] << "\t" <<
                          ordered_sums[best_i + 1] - tile_size + 1 << "\t" <<
                          ordered_sums[best_j] << "\t" <<
                          ordered_sums[next_search.first] - tile_size + 1 << "\t" <<
                          ordered_sums[next_search.second] << "\t" <<
                          best_t << "\t" <<
                          mean_inside << "\t" <<
                          mean_left << "\t" <<
                          mean_right << "\t" <<
                          copy_inside << "\t" <<
                          copy_left << "\t" <<
                          copy_right << "\t" <<
                          1 - (double)(best_greater_than) / randomizations <<  "\t" <<
                          1 - (double)(left_greater_than) / randomizations <<  "\t" <<
                          1 - (double)(right_greater_than) / randomizations << endl;
          
          if ( best_j + 1 < current_search.second )
          {
            next_search.first = best_j + 1;
            next_search.second = current_search.second;
            searching.push_back( next_search );
          }
        }
        
        else{
          out_file << ordered_sums[current_search.first] - tile_size + 1 << "\t"
                   << ordered_sums[best_i] << "\t"
                   << best_t << "\t"
                   << 1 - (double)(best_greater_than) / randomizations << endl;
          out_file << ordered_sums[best_i + 1] - tile_size + 1 << "\t"
                   << ordered_sums[best_j] << "\t"
                   << best_t << "\t"
                   << 1 - (double)(best_greater_than) / randomizations << endl;
          out_file << ordered_sums[best_j + 1] - tile_size + 1 << "\t"
                   << ordered_sums[current_search.second] << "\t"
                   << best_t << "\t"
                   << 1 - (double)(best_greater_than) / randomizations << endl;
                 
          history_file << ordered_sums[current_search.first] - tile_size + 1 << "\t" <<
                          ordered_sums[current_search.second] << "\t" <<
                          ordered_sums[best_i + 1] - tile_size + 1 << "\t" <<
                          ordered_sums[best_j] << "\t" <<
                          -1 << "\t" <<
                          -1 << "\t" <<
                          best_t << "\t" <<
                          mean_inside << "\t" <<
                          mean_left << "\t" <<
                          mean_right << "\t" <<
                          copy_inside << "\t" <<
                          copy_left << "\t" <<
                          copy_right << "\t" <<
                          1 - (double)(best_greater_than) / randomizations <<  "\t" <<
                          1 - (double)(left_greater_than) / randomizations <<  "\t" <<
                          1 - (double)(right_greater_than) / randomizations << endl;
        }
      }
      //otherwise, just write the segment.
      else
      {
        out_file << ordered_sums[current_search.first] - tile_size + 1 << "\t"
                 << ordered_sums[current_search.second] << "\t"
                 << best_t << "\t"
                 << 1 - (double)(best_greater_than) / randomizations << endl;
        history_file << ordered_sums[current_search.first] - tile_size + 1 << "\t" <<
                        ordered_sums[current_search.second] << "\t" <<
                        ordered_sums[best_i + 1] - tile_size + 1 << "\t" <<
                        ordered_sums[best_j] << "\t" <<
                        -1 << "\t" <<
                        -1 << "\t" <<
                        best_t << "\t" <<
                        mean_inside << "\t" <<
                        mean_left << "\t" <<
                        mean_right << "\t" <<
                        copy_inside << "\t" <<
                        copy_left << "\t" <<
                        copy_right << "\t" <<
                        1 - (double)(best_greater_than) / randomizations <<  "\t" <<
                        1 - (double)(left_greater_than) / randomizations <<  "\t" <<
                        1 - (double)(right_greater_than) / randomizations << endl;
      }
      
      
    }
    
    in_file.close();
    out_file.close();
  }
  
  void CoverageDistribution::smooth_segments(
                                              const Settings& settings,
                                              const string& seq_id,
                                              double summary_average,
                                              string tile_file_name,
                                              string segment_file_name,
                                              string out_file_name,
                                              string final_file_name,
                                              string gd_file_name
                                              )
  {
    //cout << "sum_avg" << summary_average << endl;
    cGenomeDiff gd;
    
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
    double z_score;
    double greater_than;
    
    
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
                  "Z_Score" << "\t" <<
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
        
        if ( position_start <= tile_entry.first &&
             position_end >= tile_entry.first )
        {
          segment_sum += tile_entry.second;
          num_tiles++;
        }
        
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
      
      if (new_segment_mean != 1.0) {
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
      
      
      //set the new value of the coverage for this segment to this value.
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
