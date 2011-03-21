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

#include "breseq/annotated_sequence.h"

using namespace breseq;
using namespace std;
  
void convert_genbank(std::string in, std::string fasta, std::string ft) {

  cAnnotatedSequence seq;
  
  // Load the GenBank file
  LoadGenBankFile(in, seq);

  // Output sequence
  seq.WriteFASTA(fasta);
  
  // Output feature table
  seq.WriteFeatureTable(ft);
}

/*! Identify mutations.
 
 This file only does the command-line parsing bit; the real work is over in
 identify_mutations.cpp.
 */
int main(int argc, char* argv[]) {
	namespace po = boost::program_options;
	
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("input,i", po::value<string>(), "input Genbank Flat File")
	("feature_table,i", po::value<string>(), "output feature table (optional)")
	("fasta,f", po::value<string>(), "output FASTA sequence file (optional)")
  ;
  
	po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("input")
		 || (!options.count("feature_table") && !options.count("fasta"))  
		 ) {
		cout << "Usage: convert_genbank --input <sequence.gbk> [--fasta <output.fasta> --feature_table <output.tab>]" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
		convert_genbank(  options["input"].as<string>(),
                      options["fasta"].as<string>(),
                      options["feature_table"].as<string>()
                  );
    } catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}
