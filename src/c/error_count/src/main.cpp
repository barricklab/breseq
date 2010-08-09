#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "error_count.h"

/*! Calculate error calibrations from FASTA and BAM reference files.
 
 Usage:
	 error_count -fasta reference.fasta -bam reference.bam
 
 This file only does the command-line parsing bit; the real work is over in
 error_count.cpp.
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
	("output,o", po::value<string>(), "output directory")
	("readfile,r", po::value<vector<string> >(), "names of readfiles (no extension)")
	("coverage", "generate unique coverage distribution output")
  ("errors", "generate unique error count output");
  
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
		cout << "Usage: error_count --bam <sequences.bam> --fasta <reference.fasta> --output <path> --readfile <filename> [--coverage] [--errors]" << endl;
		cout << cmdline_options << endl;
		exit(0);
	}
	
	// attempt to calculate error calibrations:
	try {
		breseq::error_count(options["bam"].as<string>(),
												options["fasta"].as<string>(),
												options["output"].as<string>(),
												options["readfile"].as<vector<string> >(),
                        options.count("coverage"),
                        options.count("errors"));
	} catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}
