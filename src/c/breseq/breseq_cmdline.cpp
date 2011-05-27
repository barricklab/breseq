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

#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "breseq/fastq.h"
#include "breseq/alignment_output.h"
#include "breseq/annotated_sequence.h"
#include "breseq/calculate_trims.h"
#include "breseq/error_count.h"
#include "breseq/identify_mutations.h"
#include "breseq/candidate_junctions.h"
#include "breseq/tabulate_coverage.h"


using namespace breseq;
using namespace std;
namespace po = boost::program_options;

// Utility function
bool file_exists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile;
}

/*! Analyze FASTQ
 
 Extract information about reads in a FASTQ file.
 
 */

int do_analyze_fastq(int argc, char* argv[]) {
  
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("input,i", po::value<string>(), "input FASTQ File")
  //  ("output,o", po::value<string>(),"out to file") // outputs to STDOUT for now
  ;
  
  po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("input")
		 ) {
		cout << "Usage: breseq ANALYZE_FASTQ --input input.fastq" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	try {
    cFastqSequence sequence;
    cFastqFile fastqparse(options["input"].as<std::string>(), std::fstream::in);
    
    fastqparse.check_if_file_opened();
    fastqparse.read_sequence(sequence);
    fastqparse.write_summary_file();
    
  } catch(...) {
		// failed; 
		return -1;
	}
  
  return 0;
}

/*! Convert Genbank
 
 Create a tab-delimited file of information about genes from a 
 */

// Helper function
void convert_genbank(std::string in, std::string fasta, std::string ft) {
  
  cAnnotatedSequence seq;
  
  // Load the GenBank file
  LoadGenBankFile(in, seq);
  
  // Output sequence
  seq.WriteFASTA(fasta);
  
  // Output feature table
  seq.WriteFeatureTable(ft);
}


int do_convert_genbank(int argc, char* argv[]) {
	
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
		cout << "Usage: breseq CONVERT_GENBANK --input <sequence.gbk> [--fasta <output.fasta> --feature_table <output.tab>]" << endl;
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


/*! Calculate Trims
 
 Calculate how much to ignore on the end of reads due to ambiguous alignments
 of those bases.
 
 */
int do_calculate_trims(int argc, char* argv[]) {
	using namespace std;
	namespace po = boost::program_options;
	
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("fasta,f", po::value<string>(), "FASTA file of reference sequence")
	("output,o", po::value<string>(), "output directory")
  ;
  
	po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("fasta")
		 || !options.count("output")
		 ) {
		cout << "Usage: breseq CALCULATE_TRIMS --bam <sequences.bam> --fasta <reference.fasta> --error_dir <path> --genome_diff <path> --output <path> --readfiles <filename> --coverage_dir <dirname> [--minimum-quality-score 3]" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
		calculate_trims(options["fasta"].as<string>(),options["output"].as<string>());
  } catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}

/*! Candidate Junctions
 
 Calculate how much to ignore on the end of reads due to ambiguous alignments
 of those bases.
 
 */
int do_candidate_junctions(int argc, char* argv[]) {
	using namespace std;
	namespace po = boost::program_options;
	
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("fasta,f", po::value<string>(), "FASTA file of reference sequence")
	("output,o", po::value<string>(), "output directory")
  
  //These options are almost always default values
  ("required-both-unique-length-per-side,1", po::value<uint32> >(), 
   "Only count reads where both matches extend this many bases outside of the overlap.")
  ("required-one-unique-length-per-side,2", po::value<uint32> >(), 
   "Only count reads where at least one match extends this many bases outside of the overlap.")
  ("maximum-inserted-junction-sequence-length,3", po::value<uint32> >(), 
   "Maximum number of bases allowed in the overlapping part of a candidate junction.")
  ("required-match-length,4", po::value<uint32> >(), 
   "At least this many bases in the read must match the reference genome for it to count.")
  ("required-extra-pair-total-length,5", po::value<uint32> >(), 
   "Each match pair must have at least this many bases not overlapping for it to count.")
  ;
  
	po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("fasta")
		 || !options.count("output")
		 ) {
		cout << "Usage: breseq CALCULATE_TRIMS --bam <sequences.bam> --fasta <reference.fasta> --error_dir <path> --genome_diff <path> --output <path> --readfiles <filename> --coverage_dir <dirname> [--minimum-quality-score 3]" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
		candidate_junctions cj(
      options["fasta"].as<string>(),
      options["output"].as<string>());
    cj.identify_candidate_junctions();
    
  } catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}


/*! Error Count
 
Calculate error calibrations from FASTA and BAM reference files.
 
 */
int do_error_count(int argc, char* argv[]) {
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
		std::cout << "Usage: breseq ERROR_COUNT --bam <sequences.bam> --fasta <reference.fasta> --output <path> --readfile <filename> [--coverage] [--errors] [--minimum-quality-score 3]" << std::endl;
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


/*! Identify mutations
 
 Identify read-alignment and missing coverage mutations by going through
 the pileup of reads at each position in the reference genome.
 
 */
int do_identify_mutations(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("bam,b", po::value<string>(), "bam file containing sequences to be aligned")
	("fasta,f", po::value<string>(), "FASTA file of reference sequence")
	("readfiles,r", po::value<vector<string> >(), "names of readfiles (no extension)")
	("error_dir,e", po::value<string>(), "Directory containing error rates files")
	("error_table", po::value<string>()->default_value(""), "Error rates files")  
	("genome_diff,g", po::value<string>(), "Genome diff file")
	("output,o", po::value<string>(), "output directory")
	("coverage_dir", po::value<string>()->default_value(string("")), "directory for coverage files")
	("mutation_cutoff,c", po::value<double>()->default_value(2.0), "mutation cutoff (log10 e-value)")
	("deletion_propagation_cutoff,u", po::value<vector<double> >(), "number after which to cutoff deletions")	
  ("minimum_quality_score", po::value<int>()->default_value(0), "ignore base quality scores lower than this")
	("predict_deletions,d", po::value<bool>()->default_value(true), "whether to predict deletions")
	("predict_polymorphisms,p", po::value<bool>()->default_value(false), "whether to predict polymorphisms")
  ("polymorphism_cutoff", po::value<double>()->default_value(2.0), "polymorphism cutoff (log10 e-value)") 
  ("polymorphism_frequency_cutoff", po::value<double>()->default_value(0.0), "ignore polymorphism predictions below this frequency")
  ("per_position_file", po::value<bool>()->default_value(false), "print out verbose per position file")
  ;
  
	po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
  
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("bam")
		 || !options.count("fasta")
		 || !options.count("error_dir")
		 || !options.count("genome_diff")
		 || !options.count("output")
		 || !options.count("readfiles")
		 || !options.count("coverage_dir")
		 || !options.count("deletion_propagation_cutoff") 
     
		 ) {
		cout << "Usage: breseq IDENTIFY_MUTATIONS --bam <sequences.bam> --fasta <reference.fasta> --error_dir <path> --genome_diff <path> --output <path> --readfiles <filename> --coverage_dir <dirname>" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
		identify_mutations(
                         options["bam"].as<string>(),
                         options["fasta"].as<string>(),
                         options["error_dir"].as<string>(),
                         options["genome_diff"].as<string>(),
                         options["output"].as<string>(),
                         options["readfiles"].as<vector<string> >(),
                         options["coverage_dir"].as<string>(),
                         options["deletion_propagation_cutoff"].as<vector<double> >(),
                         options["mutation_cutoff"].as<double>(),
                         options["predict_deletions"].as<bool>(),
                         options["predict_polymorphisms"].as<bool>(),
                         options["minimum_quality_score"].as<int>(),
                         options["polymorphism_cutoff"].as<double>(),
                         options["polymorphism_frequency_cutoff"].as<double>(),
                         options["error_table"].as<string>(),
                         options["per_position_file"].as<bool>()
                         );
	} catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}


/* Candidate Junctions
 
 Perform an analysis of aligned reads to predict the best possible new candidate junctions.
 
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

/*! Tabulate coverage.
 */
int do_tabulate_coverage(int argc, char* argv[]) {
	
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


/*! breseq commands
 
    First argument is a command that should be removed from argv.
    Down the road -- if there is no command, pass to actual breseq pipeline
 
 */
int main(int argc, char* argv[]) {
	
  // Extract the sub-command argument
  string command;
  if (argc > 1) {
    command = argv[1];
  } else {
    cout << "No [command] provided." << endl;
    cout << "Usage: breseq [command] options ..." << endl;
		return -1;
  }
  
  // Pass the command to the proper handler
  boost::to_upper(command);
  if (command == "ANALYZE_FASTQ") {
    return do_analyze_fastq(argc, argv);
  } else if (command == "CONVERT_GENBANK") {
      return do_convert_genbank(argc, argv);
  } else if (command == "CALCULATE_TRIMS") {
    return do_calculate_trims(argc, argv);    
  } else if (command == "ERROR_COUNT") {
    return do_error_count(argc, argv); 
  } else if (command == "IDENTIFY_MUTATIONS") {
    return do_identify_mutations(argc, argv);
  } else if (command == "CANDIDATE_JUNCTIONS") {
  return do_candidate_junctions(argc, argv); 
  } else if (command == "TABULATE_COVERAGE") {
    return do_tabulate_coverage(argc, argv); 
  }
  
  // Command was not recognized. Should output valid commands.
  cout << "Unrecognized command" << endl;
  return -1;
}
