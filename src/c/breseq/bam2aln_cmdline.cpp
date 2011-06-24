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

#include "breseq/anyoption.h"
#include "breseq/alignment_output.h"

using namespace std;
using namespace breseq;

int main(int argc, char* argv[]) {
//	cerr << "OUTPUT:    " << argv[0] << endl;
  // Options should eventually be the same as Perl bam2aln.
  // Concentrate first on HTML output.

  // setup and parse configuration options:
	AnyOption options("Usage: bam2aln --bam=<reference.bam> --fasta=<reference.fasta> --region=<accession:start-end> --output=<output.html> [--max-reads=1000]");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  ("bam,b", "bam file containing sequences to be aligned")
	("fasta,f", "FASTA file of reference sequence")
  ("output,o", "name of output file")
  ("region,r", "region to print (accession:start-end)", "")
  ("max-reads,n", "maximum number of reads to show in alignment", 1000)
  .processCommandArgs(argc, argv);

  //("output,o", "out to files") // outputs to STDOUT for now

  
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("bam")
		 || !options.count("fasta")
		 || !options.count("output")) {
		options.printUsage();
		return -1;
	}

  // generate alignment!
	try {
		alignment_output ao(
                        options["bam"],
                        options["fasta"],
                        from_string<uint32_t>(options["max-reads"])
                        );
    
		ao.create_alignment(
                        options["bam"],
                        options["fasta"],
                        options["region"]
                        );
   string html_output = ao.html_alignment(options["region"]);
   cout << html_output << endl;


  } catch(...) {
		// failed;
		return -1;
	}	return 0;
}
