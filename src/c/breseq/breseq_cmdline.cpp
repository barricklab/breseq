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
#include "breseq/resolve_alignments.h"
#include "breseq/tabulate_coverage.h"
#include "breseq/settings.h"


using namespace breseq;
using namespace std;
namespace po = boost::program_options;

/*! Analyze FASTQ
 
 Extract information about reads in a FASTQ file.
 
 */

int do_analyze_fastq(int argc, char* argv[]) {
  
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("input,i", po::value<string>(), "input FASTQ file")
  ("convert,c", po::value<string>(), "converted FASTQ file (created only if necessary)")
  //  ("output,o", po::value<string>(),"out to file") // outputs to STDOUT for now
  ;
  
  po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("input")
     || !options.count("convert")
		 ) {
		cout << "Usage: breseq ANALYZE_FASTQ --input input.fastq --convert converted.fastq" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	try {
    
    analyze_fastq(options["input"].as<std::string>(), options["convert"].as<std::string>());
        
  } catch(...) {
		// failed; 
		return -1;
	}
  
  return 0;
}

/*! Convert Genbank
 
 Create a tab-delimited file of information about genes and a
 FASTA sequence file from an input GenBank file.
 */

// Helper function
void convert_genbank(const vector<string>& in, const string& fasta, const string& ft) {
  
  cReferenceSequences refseqs;
  
  // Load the GenBank file
  LoadGenBankFile(refseqs, in);
  
  // Output sequence
  if (fasta != "") refseqs.WriteFASTA(fasta);
  
  // Output feature table
  if (ft != "") refseqs.WriteFeatureTable(ft);
}


int do_convert_genbank(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("input,i", po::value<vector<string> >(), "input GenBank flatfile (multiple allowed)")
	("features,g", po::value<string>()->default_value(""), "output feature table (optional)")
	("fasta,f", po::value<string>()->default_value(""), "output FASTA sequence file (optional)")
  ;
  
	po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("input")
		 || (!options.count("features") && !options.count("fasta"))  
		 ) {
		cout << "Usage: breseq CONVERT_GENBANK --input <sequence.gbk> [--fasta <output.fasta> --features <output.tab>]" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
		convert_genbank(  
                    options["input"].as<vector<string> >(),
                    options["fasta"].as<string>(),
                    options["features"].as<string>()
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
	namespace po = boost::program_options;
	
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
	("fasta,f", po::value<string>(), "FASTA file of reference sequence")
	("output,o", po::value<string>(), "output directory")
  
  //These options are almost always default values
  ("required-both-unique-length-per-side,1", po::value<uint32_t>(), 
   "Only count reads where both matches extend this many bases outside of the overlap.")
  ("required-one-unique-length-per-side,2", po::value<uint32_t>(), 
   "Only count reads where at least one match extends this many bases outside of the overlap.")
  ("maximum-inserted-junction-sequence-length,3", po::value<uint32_t>(), 
   "Maximum number of bases allowed in the overlapping part of a candidate junction.")
  ("required-match-length,4", po::value<uint32_t>(), 
   "At least this many bases in the read must match the reference genome for it to count.")
  ("required-extra-pair-total-length,5", po::value<uint32_t>(), 
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
    
  } catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}

/*!  Resolve alignments
 
 Compare matches to candidate junctions and matches to reference sequence
 to determine what junctions have support. Also splits mosaic matches.
 
 */
int do_resolve_alignments(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
	("help,h", "produce this help message")
  ("junction-prediction,p", po::value<bool>(), "whether to predict new junctions")
	("reference-fasta,f", po::value<string>(), "FASTA file of reference sequences")
  ("junction-fasta,j", po::value<string>(), "FASTA file of candidate junction sequences")
  ("reference-sam-path", po::value<string>(), "path to SAM files of read alignments to reference sequences")
  ("junction-sam-path", po::value<string>(), "path to SAM files of read alignments to candidate junction sequences")
  ("resolved-path,o", po::value<string>(), "output path for resolved sam files")
  ("data-path,o", po::value<string>(), "data path")
  ("features,g", po::value<string>(), "feature table file for reference sequences")
	("read-file,r", po::value<vector<string> >(), "FASTQ read files (multiple allowed) ") 
  ("max-read-length,m", po::value<uint32_t>(), "number of flanking bases in candidate junctions")
  ;
  
	po::variables_map options;
	po::store(po::parse_command_line(argc, argv, cmdline_options), options);
	po::notify(options);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("junction-prediction")
		 || !options.count("reference-fasta")
     || !options.count("junction-fasta")
		 || !options.count("reference-sam-path")
		 || !options.count("junction-sam-path")
		 || !options.count("resolved-path")
     || !options.count("data-path")
		 || !options.count("features")
     || !options.count("read-file")
		 || !options.count("max-read-length")

		 ) {
		cout << "Usage: breseq RESOLVE_ALIGNMENTS ... " << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
    cReadFiles rf(options["read-file"].as<vector<string> >());
    
    resolve_alignments(
      options["junction-prediction"].as<bool>(),
      options["reference-fasta"].as<string>(),
      options["junction-fasta"].as<string>(),
      options["reference-sam-path"].as<string>(),
      options["junction-sam-path"].as<string>(),
      options["resolved-path"].as<string>(),
      options["data-path"].as<string>(),
      options["features"].as<string>(),
      rf,
      options["max-read-length"].as<uint32_t>()
    );
    
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

int do_preprocess_alignments(int argc, char* argv[]) {

	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
		("help,h", "produce this help message")
		("fasta,f", po::value<string>(), "reference sequences in FASTA format")
		("sam,s", po::value<string>(), "text SAM file of input alignments")
		("output,o", po::value<string>(), "output FASTA file of candidate junctions")
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
		cout << "Usage: breseq_utils PREPROCESS_ALIGNMENTS --fasta=reference.fasta "
         << " --sam=reference.sam --output=output.fasta" << endl;
		cout << cmdline_options << endl;
		return -1;
	}

	// attempt to calculate error calibrations:
	try {

    // plain function

	CandidateJunction::Settings settings;
	CandidateJunction::Summary summary;
	/*struct CandidateJunction::RefSeqInfo refSeqInfo(
		options["sam"].as<string>(),
		options["fasta"].as<string>(),
		options["output"].as<string>()
	);*/

	settings.candidate_junction_fasta_file_name = options["candidate_junction_fasta_file_name"].as<string>();
	settings.candidate_junction_faidx_file_name = options["candidate_junction_faidx_file_name"].as<string>();
	settings.candidate_junction_sam_file_name = options["candidate_junction_sam_file_name"].as<string>();
	settings.candidate_junction_score_method = options["candidate_junction_score_method"].as<string>();

	settings.jc_genome_diff_file_name = options["jc_genome_diff_file_name"].as<string>();
	settings.preprocess_junction_split_sam_file_name = options["preprocess_junction_split_sam_file_name"].as<string>();
	settings.preprocess_junction_best_sam_file_name = options["preprocess_junction_best_sam_file_name"].as<string>();
	settings.reference_fasta_file_name = options["reference_fasta_file_name"].as<string>();
	settings.reference_faidx_file_name = options["reference_faidx_file_name"].as<string>();
	settings.reference_sam_file_name = options["reference_sam_file_name"].as<string>();
	settings.resolved_reference_sam_file_name = options["resolved_reference_sam_file_name"].as<string>();
	settings.resolved_junction_sam_file_name = options["resolved_junction_sam_file_name"].as<string>();
	settings.unmatched_read_file_name = options["unmatched_read_file_name"].as<string>();

	settings.no_junction_prediction = options["no_junction_prediction"].as<bool>();
	settings.unmatched_reads = options["unmatched_reads"].as<bool>();
	settings.add_split_junction_sides = options["add_split_junction_sides"].as<bool>();
	settings.require_complete_match = options["require_complete_match"].as<bool>();

	settings.alignment_read_limit = options["alignment_read_limit"].as<int>();
	settings.candidate_junction_read_limit = options["candidate_junction_read_limit"].as<int>();
	settings.max_read_length = options["max_read_length"].as<int>();
	settings.maximum_read_mismatches = options["maximum_read_mismatches"].as<int>();
	settings.required_match_length = options["required_match_length"].as<int>();

	if (options.count("preprocess_junction_min_indel_split_length"))
		settings.preprocess_junction_min_indel_split_length = options["preprocess_junction_min_indel_split_length"].as<int>();

	cReferenceSequences ref_seq_info;
	breseq::LoadFeatureIndexedFastaFile(ref_seq_info, "", options["fasta"].as<string>());
	//ref_seq_info[0].m_fasta_sequence.m_sequence

	CandidateJunction::preprocess_alignments(settings, summary, ref_seq_info);

  } catch(...) {
		// failed;
		return -1;
	}

  return 0;
}

int do_identify_candidate_junctions(int argc, char* argv[]) {

	// setup and parse configuration options:
	po::options_description cmdline_options("Allowed options");
	cmdline_options.add_options()
		("help,h", "produce this help message")
		("fasta,f", po::value<string>(), "reference sequences in FASTA format")
		("sam,s", po::value<string>(), "text SAM file of input alignments")
		("output,o", po::value<string>(), "output FASTA file of candidate junctions")
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
		cout << "Usage: breseq_utils IDENTIFY_CANDIDATE_JUNCTIONS --fasta=reference.fasta "
         << " --sam=reference.sam --output=output.fasta" << endl;
		cout << cmdline_options << endl;
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
    
    // plain function

	CandidateJunction::Settings settings;
	CandidateJunction::Summary summary;

	cReferenceSequences ref_seq_info;
	breseq::LoadFeatureIndexedFastaFile(ref_seq_info, "", options["fasta"].as<string>());
	//ref_seq_info[0].m_fasta_sequence.m_sequence

	CandidateJunction::identify_candidate_junctions(settings, summary, ref_seq_info);
    
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
  } else if (command == "PREPROCESS_ALIGNMENTS") {
    return do_preprocess_alignments(argc, argv);
  } else if (command == "IDENTIFY_CANDIDATE_JUNCTIONS") {
    return do_identify_candidate_junctions(argc, argv);
  } else if (command == "TABULATE_COVERAGE") {
    return do_tabulate_coverage(argc, argv); 
  } else if (command == "RESOLVE_ALIGNMENTS") {
  return do_resolve_alignments(argc, argv); 
}
  
  // Command was not recognized. Should output valid commands.
  cout << "Unrecognized command" << endl;
  return -1;
}
