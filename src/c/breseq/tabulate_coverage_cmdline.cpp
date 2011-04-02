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
#include <fstream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "breseq/tabulate_coverage.h"

bool file_exists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile;
}

/*! Tabulate coverage.
 */
int main(int argc, char* argv[]) {
	using namespace std;
	namespace po = boost::program_options;
	
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("bam,b", po::value<string>(), "bam file containing sequences to be aligned")
	("fasta,f", po::value<string>(), "FASTA file of reference sequence")
	("output,o", po::value<string>(), "name of output file")
	("region,r", po::value<string>()->default_value(""), "region to print (accession:start-end)")
	("downsample,d", po::value<int>()->default_value(1), "Only print information every this many positions")
	("read_start_output,d", po::value<string>()->default_value(""), "Create additional data file binned by read start bases")
	("gc_output,d", po::value<string>()->default_value(""), "Create additional data file binned by GC content of reads")
  ;
  
	po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("bam")
		 || !options.count("fasta")
		 || !options.count("output")     
  ) {
		cout << "Usage: tabulate_coverage --bam <sequences.bam> --fasta <reference.fasta> --output <coverage.tab> [--downsample <float> --detailed]" << endl;
		cout << cmdline_options << endl;
		return -1;
	}  
                     
  
	// attempt to calculate error calibrations:
	try {
		breseq::tabulate_coverage(
                                 options["bam"].as<string>(),
                                 options["fasta"].as<string>(),
                                 options["output"].as<string>(),
                                 options["region"].as<string>(),
                                 options["downsample"].as<int>(),
                                 options["read_start_output"].as<string>(),
                                 options["gc_output"].as<string>()

                              );
	} catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}
