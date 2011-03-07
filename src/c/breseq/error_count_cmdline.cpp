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
#include "breseq/error_count.h"

/*! Calculate error calibrations from FASTA and BAM reference files.
 
 This file only does the command-line parsing bit; the real work is over in
 error_count.cpp.
 */
int main(int argc, char* argv[]) {
	namespace po = boost::program_options;
	
  std::cout << "temp" << std::endl;
  
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("bam,b", po::value<std::string>(), "bam file containing sequences to be aligned")
	("fasta,f", po::value<std::string>(), "FASTA file of reference sequence")
	("output,o", po::value<std::string>(), "output directory")
	("readfile,r", po::value<std::vector<std::string> >(), "names of readfiles (no extension)")
	("coverage", "generate unique coverage distribution output")
  ("errors", "generate unique error count output")
  ("minimum-quality-score", po::value<int>()->default_value(0), "ignore base quality scores lower than this")
  ("covariates", po::value<std::string>()->default_value(""), "covariates for error model")
  ;

	po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("bam")
		 || !options.count("fasta")
		 || !options.count("output")
		 || !options.count("readfile")
		 || (!options.count("coverage") && !options.count("errors")) ) {
		std::cout << "Usage: error_count --bam <sequences.bam> --fasta <reference.fasta> --output <path> --readfile <filename> [--coverage] [--errors] [--minimum-quality-score 3]" << std::endl;
		std::cout << cmdline_options << std::endl;
		return -1;
	}
	
	// attempt to calculate error calibrations:
	try {
		breseq::error_count(options["bam"].as<std::string>(),
												options["fasta"].as<std::string>(),
												options["output"].as<std::string>(),
												options["readfile"].as<std::vector<std::string> >(),
												options.count("coverage"),
                        options.count("errors"),
                        options["minimum-quality-score"].as<int>(),
                        options["covariates"].as<std::string>()
                      );
	} catch(...) {
		// failed; 
    std::cout << "<<<Failed>>>" << std::endl;
		return -1;
	}
  
	return 0;
}
