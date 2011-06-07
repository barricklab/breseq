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

#include <boost/program_options.hpp>

#include "breseq/alignment_output.h"

using namespace std;
using namespace breseq;
namespace po = boost::program_options;

int main(int argc, char* argv[]) {
//	cerr << "OUTPUT:    " << argv[0] << endl;
  // Options should eventually be the same as Perl bam2aln.
  // Concentrate first on HTML output.

	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("bam,b", po::value<string>(), "bam file containing sequences to be aligned")
	("fasta,f", po::value<string>(), "FASTA file of reference sequence")
  ("region,r", po::value<string>()->default_value(""), "region to print (accession:start-end)")
  ("output,o", po::value<string>(), "name of output file")
	("max-reads,n", po::value<uint32_t>()->default_value(1000), "maximum number of reads to show in alignment")
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
		cout << "Usage: bam2aln --bam=<reference.bam> --fasta=<reference.fasta> --region=<accession:start-end> "
         << "--output=<output.html> [--max-reads=1000]" << endl;
		cout << cmdline_options << endl;
		
		return -1;
	}

  // generate alignment!
	try {
		alignment_output ao(options["bam"].as<string>(),
                        options["fasta"].as<string>(),
                        options["max-reads"].as<uint32_t>());

    string aln = ao.html_alignment(options["region"].as<string>());

    // add: write alignment to output file....

  } catch(...) {
		// failed;
		return -1;
	}

	return 0;
}
