/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2010 Michigan State University

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "breseq/candidate_junctions.cpp"

using namespace breseq;
using namespace std;

namespace po = boost::program_options;

/*! 
 
 Working to condense all internal breseq steps that can be called separately from perl into one binary...
 Initially work to
 
 */

/* Perform an analysis of aligned reads to predict the best possible new candidate junctions.
 */

int do_candidate_junctions(int argc, char* argv[]) {
    
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("fasta,f", po::value<string>(), "reference sequences in FASTA format")
  ("sam,s", po::value<string>(), "text SAM file of input alignments")
  ("output,o", po::value<string>(), "output FASTA file of candidate junctions")

// additional options  
//  $required_both_unique_length_per_side = $settings->{required_both_unique_length_per_side};
//	$required_one_unique_length_per_side = $settings->{required_one_unique_length_per_side};
//	$maximum_inserted_junction_sequence_length = $settings->{maximum_inserted_junction_sequence_length};
// 	$required_match_length = $settings->{required_match_length};
//	$required_extra_pair_total_length = $settings->{required_extra_pair_total_length};
  
//  ("output,o", po::value<string>(),"out to file") // outputs to STDOUT for now
  ;

  po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("fasta")
		 || !options.count("sam")
     || !options.count("output")
		 ) {
		cout << "Usage: breseq_utils CANDIDATE_JUNCTIONS --fasta=reference.fasta " 
         << " --sam=reference.sam --output=output.fasta" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
    
    // plain function
    candidate_junctions(options["fasta"].as<string>(),
                        options["sam"].as<string>(),
                        options["output"].as<string>() );
    
  } catch(...) {
		// failed; 
		return -1;
	}
  
  return 0;
}

/*! fastq_utils command line tool
 
    First argument is a sub-command that should be removed from argv.
 
 */
int main(int argc, char* argv[]) {
	
  // Extract the sub-command argument
  std::string subcommand;
  if (argc > 1) {
    subcommand = argv[1];
  } else {
    cout << "No [command] provided." << endl << endl;
    cout << "Usage: breseq_utils [command] options ..." << endl;
		return -1;
  }
  
  boost::to_upper(subcommand);
  
  if (subcommand == "CANDIDATE_JUNCTIONS") {
    return do_candidate_junctions(argc, argv);
  }
  
  cout << "Unrecognized command" << endl;

	return -1;
}
