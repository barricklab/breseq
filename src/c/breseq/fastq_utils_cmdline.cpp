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

#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "breseq/fastq.h"

using namespace breseq;
using namespace std;

namespace po = boost::program_options;

/*! analyze_fastq
 
 Extract information about 
 
 */

int analyze_fastq(int argc, char* argv[]) {
    
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("input,i", po::value<string>(), "input FASTQ File")
  ("output,o", po::value<string>(),"out to file")
	("summary,s", po::value<string>(), "output summary file (optional)")
  ;

  po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("input")
		 ) {
		cout << "Usage: fastq_utils --input input.fastq [--summary <summary.tab> --qscore <qscore.tab>]" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
    
    
    cFastqSequence sequence;
    cFastqFile fastqparse(options["input"].as<std::string>(), options["output"].as<std::string>(), std::fstream::in);
    
    fastqparse.check_if_file_opened();
    fastqparse.read_sequence(sequence);
    fastqparse.write_summary_file();
    
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
    cout << "Usage: fastq_utils [command] options ..." << endl;
		return -1;
  }
  
  if (subcommand == "analyze") {
    return analyze_fastq(argc, argv);
  }
  
  cout << "Unrecognized command" << endl;

	return -1;
}
