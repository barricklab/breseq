#include <boost/program_options.hpp>
#include <iostream>
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
	("fasta,f", po::value<string>(), "indexed fasta file of reference sequence");
	
	po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("bam")
		 || !options.count("fasta")) {
		cout << "Usage: error_count -bam <sequences.bam> -fasta <reference.fai>" << endl;
		cout << cmdline_options << endl;
		return -1;
	}
	
	// attempt to calculate error calibrations:
	try {
		breseq::error_count(options["bam"].as<string>(), options["fasta"].as<string>());
	} catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}
