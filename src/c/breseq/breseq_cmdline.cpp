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

#include "breseq/anyoption.h"
#include "breseq/alignment_output.h"
#include "breseq/annotated_sequence.h"
#include "breseq/calculate_trims.h"
#include "breseq/candidate_junctions.h"
#include "breseq/contingency_loci.h"
#include "breseq/error_count.h"
#include "breseq/fastq.h"
#include "breseq/genome_diff.h"
#include "breseq/identify_mutations.h"
#include "breseq/resolve_alignments.h"
#include "breseq/settings.h"
#include "breseq/tabulate_coverage.h"


using namespace breseq;
using namespace std;


/*! Analyze FASTQ
 
 Extract information about reads in a FASTQ file.
 
 */

int do_analyze_fastq(int argc, char* argv[]) {
  
	// setup and parse configuration options:
	AnyOption options("Usage: breseq ANALYZE_FASTQ --input input.fastq --convert converted.fastq");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		("input,i", "input FASTQ file")
		("convert,c", "converted FASTQ file (created only if necessary)")
		//("output,o", "out to files") // outputs to STDOUT for now
	.processCommandArgs(argc, argv);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("input")
     || !options.count("convert")
		 ) {
		options.printUsage();
		return -1;
	}                       
  
	try {
    
    analyze_fastq(options["input"], options["convert"]);
        
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
void convert_genbank(const vector<string>& in, const string& fasta, const string& ft, const string& gff3 ) {
  
  cReferenceSequences refseqs;
  
  // Load the GenBank file
  LoadGenBankFile(refseqs, in);
  
  // Output sequence
  if (fasta != "") refseqs.WriteFASTA(fasta);
  
  // Output feature table
  if (ft != "") refseqs.WriteFeatureTable(ft);
     
  if (gff3 != "" ) refseqs.WriteGFF( gff3 );
}


int do_convert_genbank(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq CONVERT_GENBANK --input <sequence.gbk> [--fasta <output.fasta> --features <output.tab>]");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		("input,i", "input GenBank flatfile (multiple allowed, comma-separated)")
		("features,g", "output feature table", "")
		("fasta,f", "FASTA file of reference sequences", "")
    ("gff3,v", "GFF file of features", "" )
	.processCommandArgs(argc, argv);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("input")
		 || (!options.count("features") && !options.count("fasta"))  
		 ) {
		options.printUsage();
		return -1;
	}
  
	// attempt to calculate error calibrations:
	try {
        
		convert_genbank(  
                    from_string<vector<string> >(options["input"]),
                    options["fasta"],
                    options["features"],
                    options["gff3"]
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
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq CALCULATE_TRIMS --bam <sequences.bam> --fasta <reference.fasta> --error_dir <path> --genome_diff <path> --output <path> --readfile <filename> --coverage_dir <dirname> [--minimum-quality-score 3]");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		("fasta,f", "FASTA file of reference sequence")
		("output,o", "output directory")
	.processCommandArgs(argc, argv);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("fasta")
		 || !options.count("output")
		 ) {
		options.printUsage();
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
		calculate_trims(options["fasta"],options["output"]);
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
	AnyOption options("Usage: breseq RESOLVE_ALIGNMENTS ... ");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		("junction-prediction,p", "whether to predict new junctions")
		("junction-fasta,j", "FASTA file of candidate junction sequences")
// convert to basing everything off the main output path, so we don't have to set so many options
    ("path", "path to breseq output")
		("reference-sam-path", "path to SAM files of read alignments to reference sequences")
		("junction-sam-path", "path to SAM files of read alignments to candidate junction sequences")
		("resolved-path", "output path for resolved sam files")
		("data-path", "data path")
		("read-file,r", "FASTQ read files (multiple allowed, comma-separated) ")
		("max-read-length,m", "number of flanking bases in candidate junctions")
		("alignment-read-limit", "maximum number of alignments to process. DEFAULT = 0 (OFF).", 0)
	.processCommandArgs(argc, argv);

	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("junction-prediction")
     || !options.count("junction-fasta")
     || !options.count("path")
		 || !options.count("reference-sam-path")
		 || !options.count("junction-sam-path")
		 || !options.count("resolved-path")
     || !options.count("data-path")
     || !options.count("read-file")
		 || !options.count("max-read-length")

		 ) {
		options.printUsage();
		return -1;
	}                       
  
	try {
    
    Settings settings(options["path"]);

    
    Summary summary;

    cReadFiles rf(from_string<vector<string> >(options["read-file"]));
    
    // Load the reference sequence info
    cReferenceSequences ref_seq_info;
    LoadFeatureIndexedFastaFile(ref_seq_info, settings.reference_features_file_name, settings.reference_fasta_file_name);
    
    resolve_alignments(
      settings,
      summary,
      ref_seq_info,
      from_string<bool>(options["junction-prediction"]),
      options["reference-sam-path"],
      options["junction-sam-path"],
      options["data-path"],
      rf,
      from_string<uint32_t>(options["max-read-length"]),
      from_string<uint32_t>(options["alignment-read-limit"])
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
	  
	// setup and parse configuration options:
	AnyOption options("Usage: breseq ERROR_COUNT --bam <sequences.bam> --fasta <reference.fasta> --output <path> --readfile <filename> [--coverage] [--errors] [--minimum-quality-score 3]");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		("bam,b", "bam file containing sequences to be aligned")
		("fasta,f", "FASTA file of reference sequence")
		("output,o", "output directory")
		("readfile,r", "name of readfile (no extension). may occur multiple times")
		("coverage", "generate unique coverage distribution output", TAKES_NO_ARGUMENT)
		("errors", "generate unique error count output", TAKES_NO_ARGUMENT)
    ("covariates", "covariates for error model", "")
    ("minimum-quality-score", "ignore base quality scores lower than this", 0)
	.processCommandArgs(argc, argv);
  
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("bam")
		 || !options.count("fasta")
		 || !options.count("output")
		 || !options.count("readfile")
		 || (!options.count("coverage") && !options.count("errors")) ) {
		options.printUsage();
		return -1;
	}
	
	// attempt to calculate error calibrations:
	try {
		breseq::error_count(options["bam"],
												options["fasta"],
												options["output"],
												split(options["readfile"], "\n"),
												options.count("coverage"),
                        options.count("errors"),
                        from_string<uint32_t>(options["minimum-quality-score"]),
                        options["covariates"]
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
	AnyOption options("Usage: breseq IDENTIFY_MUTATIONS --bam <sequences.bam> --fasta <reference.fasta> --error_dir <path> --genome_diff <path> --output <path> --readfile <filename> --coverage_dir <dirname>");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		("bam,b", "bam file containing sequences to be aligned")
		("fasta,f", "FASTA file of reference sequence")
		("readfile,r", "names of readfile (no extension)")
		("error_dir,e", "Directory containing error rates files")
		("error_table", "Error rates files", "")
		("genome_diff,g", "Genome diff file")
		("output,o", "output directory")
		("coverage_dir", "directory for coverage files", "")
		("mutation_cutoff,c", "mutation cutoff (log10 e-value)", 2.0L)
		("deletion_propagation_cutoff,u", "number after which to cutoff deletions")
		("minimum_quality_score", "ignore base quality scores lower than this", 0)
		("predict_deletions,d", "whether to predict deletions", TAKES_NO_ARGUMENT)
		("predict_polymorphisms,p", "whether to predict polymorphisms", TAKES_NO_ARGUMENT)
		("polymorphism_cutoff", "polymorphism cutoff (log10 e-value)", 2.0L)
		("polymorphism_frequency_cutoff", "ignore polymorphism predictions below this frequency", 0.0L)
		("per_position_file", "print out verbose per position file", TAKES_NO_ARGUMENT)
	.processCommandArgs(argc, argv);

	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("bam")
		 || !options.count("fasta")
		 || !options.count("error_dir")
		 || !options.count("genome_diff")
		 || !options.count("output")
		 || !options.count("readfile")
		 || !options.count("coverage_dir")
		 || !options.count("deletion_propagation_cutoff") 
     
		 ) {
		options.printUsage();
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
		identify_mutations(
                         options["bam"],
                         options["fasta"],
                         options["error_dir"],
                         options["genome_diff"],
                         options["output"],
                         split(options["readfile"], "\n"),
                         options["coverage_dir"],
                         from_string<vector<double> >(options["deletion_propagation_cutoff"]),
                         from_string<double>(options["mutation_cutoff"]),
                         options.count("predict_deletions"),
                         options.count("predict_polymorphisms"),
                         from_string<int>(options["minimum_quality_score"]),
                         from_string<double>(options["polymorphism_cutoff"]),
                         from_string<double>(options["polymorphism_frequency_cutoff"]),
                         options["error_table"],
                         options.count("per_position_file")
                         );
	} catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}

/*! Contingency Loci
 
 Analyze lengths of homopolymer repeats in mixed samples.
 
 */
int do_contingency_loci(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq CONTINGENCY_LOCI --bam <sequences.bam> --fasta <reference.fasta> --output <path>");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  ("bam,b", "bam file containing sequences to be aligned")
  ("fasta,f", "FASTA file of reference sequence")
  ("output,o", "output file")
	.processCommandArgs(argc, argv);
  
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("bam")
		 || !options.count("fasta")
		 || !options.count("output")
     
		 ) {
		options.printUsage();
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
		analyze_contingency_loci(
                       options["bam"],
                       options["fasta"],
                       options["output"]
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
	AnyOption options("Usage: breseq_utils PREPROCESS_ALIGNMENTS --fasta=reference.fasta --sam=reference.sam --output=output.fasta");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		("data-path", "path of data")
		("reference-alignment-path", "path where input alignment")
		("candidate-junction-path", "path where candidate junction files will be created")
		("read-file,r", "FASTQ read files (multiple allowed, comma-separated) ")

		("candidate-junction-score-method", "scoring method", "POS_HASH")
		("min-indel-split-length", "split indels this long in matches", 3)
		("max-read-mismatches", "ignore reads with more than this number of mismatches", uint32_t(0))
		("require-complete-match", "require the complete read to match (both end bases", TAKES_NO_ARGUMENT)
		("required-match-length", "require this length of sequence -- on the read -- to match", static_cast<uint32_t>(28))
		("candidate-junction-read-limit", "limit handled reads to this many", -1)
    .processCommandArgs(argc, argv);

	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("data-path")
     || !options.count("reference-alignment-path")
     || !options.count("candidate-junction-path")
     || !options.count("read-file")
		 ) {
		options.printUsage();
		return -1;
	}

	try {

  Summary summary;
    
  // Set the things we need...
	Settings settings;
    
  settings.read_structures.Init(from_string<vector<string> >(options["read-file"]));
 
	settings.candidate_junction_fasta_file_name = options["candidate-junction-path"];
  settings.candidate_junction_fasta_file_name += "/candidate_junctions.fasta";
	settings.candidate_junction_faidx_file_name = settings.candidate_junction_fasta_file_name + ".fai";
    
	settings.candidate_junction_sam_file_name = options["candidate-junction-path"];
  settings.candidate_junction_sam_file_name += "/#.candidate_junction.sam";
    
	settings.candidate_junction_score_method = options["candidate-junction-score-method"];
    
	settings.preprocess_junction_split_sam_file_name = options["candidate-junction-path"];
  settings.preprocess_junction_split_sam_file_name += "/#.split.sam";
    
	settings.preprocess_junction_best_sam_file_name = options["candidate-junction-path"];
  settings.preprocess_junction_best_sam_file_name += "/best.sam";
    
	settings.reference_fasta_file_name = options["data-path"] + "/reference.fasta";
  settings.reference_faidx_file_name = settings.reference_fasta_file_name + ".fai";
    
	settings.reference_sam_file_name = options["reference-alignment-path"];
  settings.reference_sam_file_name += "/#.reference.sam";
    
  settings.max_read_mismatches = from_string<int32_t>(options["max-read-mismatches"]);
  settings.require_complete_match = options.count("require-complete-match");
	settings.candidate_junction_read_limit = from_string<int32_t>(options["candidate-junction-read-limit"]);
	settings.required_match_length = from_string<int32_t>(options["required-match-length"]);
  settings.preprocess_junction_min_indel_split_length = from_string<int32_t>(options["min-indel-split-length"]);
 
	cReferenceSequences ref_seqs;
	breseq::LoadFeatureIndexedFastaFile(ref_seqs, "", settings.reference_fasta_file_name);

	CandidateJunctions::preprocess_alignments(settings, summary, ref_seqs);

  } catch(...) {
		// failed;
		return -1;
	}

  return 0;
}

int do_identify_candidate_junctions(int argc, char* argv[]) {

	// setup and parse configuration options:
	AnyOption options("Usage: breseq_utils IDENTIFY_CANDIDATE_JUNCTIONS --fasta=reference.fasta --sam=reference.sam --output=output.fasta");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("candidate-junction-path", "path where candidate junction files will be created")
    ("data-path", "path of data")
    ("read-file,r", "FASTQ read files (multiple allowed, comma-separated)")

    ("candidate-junction-read-limit", "limit handled reads to this many", static_cast<unsigned long>(0))
    ("required-both-unique-length-per-side,1",
     "Only count reads where both matches extend this many bases outside of the overlap.", static_cast<unsigned long>(5))
    ("required-one-unique-length-per-side,2",
     "Only count reads where at least one match extends this many bases outside of the overlap.", static_cast<unsigned long>(10))
    ("maximum-inserted-junction-sequence-length,3",
     "Maximum number of bases allowed in the overlapping part of a candidate junction.", static_cast<unsigned long>(20))
    ("required-match-length,4",
     "At least this many bases in the read must match the reference genome for it to count.", static_cast<unsigned long>(28))
    ("required-extra-pair-total-length,5",
     "Each match pair must have at least this many bases not overlapping for it to count.", static_cast<unsigned long>(2))
    ("maximum-read-length", "Length of the longest read.")

    // Defaults should be moved to Settings constructor
    ("maximum-candidate-junctions",
     "Maximum number of candidate junction to create.", static_cast<unsigned long>(5000))
    ("minimum-candidate-junctions",
     "Minimum number of candidate junctions to create.", static_cast<unsigned long>(200))
    // This should be in the summary...
    ("reference-sequence-length",
     "Total length of reference sequences.")  
    ("maximum-candidate-junction-length-factor",
     "Total length of candidate junction sequences must be no more than this times reference sequence length.", static_cast<double>(0.1))    
	.processCommandArgs(argc, argv);
  
  //These options are almost always default values
  //TODO: Supply default values?

	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("data-path")
     || !options.count("candidate-junction-path")
     || !options.count("read-file")
     || !options.count("maximum-read-length")
     || !options.count("reference-sequence-length")
		 ) {
		options.printUsage();
		return -1;
	}                       
  
	try {
    
    // plain function

    Settings settings;
      
    // File settings
    settings.read_structures.Init(from_string<vector<string> >(options["read-file"]));
    settings.preprocess_junction_split_sam_file_name = options["candidate-junction-path"] + "/#.split.sam";
    settings.reference_fasta_file_name = options["data-path"] + "/reference.fasta";     
    settings.candidate_junction_fasta_file_name = options["candidate-junction-path"] + "/candidate_junction.fasta";

    // Other settings
    settings.candidate_junction_read_limit = from_string<int32_t>(options["candidate-junction-read-limit"]);
    settings.required_extra_pair_total_length = from_string<int32_t>(options["required-extra-pair-total-length"]);
    settings.required_match_length  = from_string<int32_t>(options["required-match-length"]);
    settings.maximum_inserted_junction_sequence_length = from_string<int32_t>(options["maximum-inserted-junction-sequence-length"]);
    settings.required_one_unique_length_per_side = from_string<int32_t>(options["required-one-unique-length-per-side"]);
    settings.required_both_unique_length_per_side = from_string<int32_t>(options["required-both-unique-length-per-side"]);
    settings.maximum_read_length = from_string<int32_t>(options["maximum-read-length"]);
    
    settings.maximum_candidate_junctions = from_string<int32_t>(options["maximum-candidate-junctions"]);
    settings.minimum_candidate_junctions = from_string<int32_t>(options["minimum-candidate-junctions"]);
    settings.maximum_candidate_junction_length_factor = from_string<double>(options["maximum-candidate-junction-length-factor"]);

    // We should inherit the summary object from earlier steps
    Summary summary;
    summary.sequence_conversion.total_reference_sequence_length = from_string<int32_t>(options["reference-sequence-length"]);
    
    cReferenceSequences ref_seq_info;
    breseq::LoadFeatureIndexedFastaFile(ref_seq_info, "", options["data-path"] + "/reference.fasta");
        
    CandidateJunctions::identify_candidate_junctions(settings, summary, ref_seq_info);
    
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
	AnyOption options("Usage: tabulate_coverage --bam <sequences.bam> --fasta <reference.fasta> --output <coverage.tab> [--downsample <float> --detailed]");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		("bam,b", "bam file containing sequences to be aligned")
		("fasta,f", "FASTA file of reference sequence")
		("output,o", "name of output file")
		("region,r", "region to print (accession:start-end)", "")
		("downsample,d", "Only print information every this many positions", 1)
		("read_start_output,r", "Create additional data file binned by read start bases", "")
		("gc_output,g", "Create additional data file binned by GC content of reads", "")
	.processCommandArgs(argc, argv);

	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("bam")
		 || !options.count("fasta")
		 || !options.count("output")
     ) {
		options.printUsage();
		return -1;
	}  
  
  
	// attempt to calculate error calibrations:
	try {
		breseq::tabulate_coverage(
                              options["bam"],
                              options["fasta"],
                              options["output"],
                              options["region"],
                              from_string<int>(options["downsample"]),
                              options["read_start_output"],
                              options["gc_output"]
						  );
	} catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}


int do_convert_gvf( int argc, char* argv[]){
    AnyOption options("Usage: GD2GVF --gd <genomediff.gd> --output <gvf.gvf>"); options( "gd,i","gd file to convert") ("output,o","name of output file").processCommandArgs( argc,argv);
    
    if( !options.count("gd") || !options.count("output") ){
        options.printUsage(); return -1;
    }
    
    try{
        GDtoGVF( options["gd"], options["output"] );
    } 
    catch(...){ 
        return -1; // failed 
    }
    
    return 0;
}

int do_convert_gd( int argc, char* argv[]){
    AnyOption options("Usage: VCF2GD --vcf <vcf.vcf> --output <gd.gd>"); options( "vcf,i","gd file to convert") ("output,o","name of output file").processCommandArgs( argc,argv);
    
    if( !options.count("vcf") || !options.count("output") ){
        options.printUsage(); return -1;
    }
    
    try{
        VCFtoGD( options["vcf"], options["output"] );
    } 
    catch(...){ 
        return -1; // failed 
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
	char* argv_new[argc];
	int argc_new = argc - 1;

	if (argc > 1) {

		command = argv[1];
		argv_new[0] = argv[0];
		for (int32_t i = 1; i < argc; i++)
			argv_new[i] = argv[i + 1];

	} else {
		cout << "No [command] provided." << endl;
		cout << "Usage: cbreseq [command] options ..." << endl;
		return -1;
	}

	// Pass the command to the proper handler
	command = to_upper(command);
	if (command == "ANALYZE_FASTQ") {
		return do_analyze_fastq(argc_new, argv_new);
	} else if (command == "CALCULATE_TRIMS") {
		return do_calculate_trims(argc_new, argv_new);
	} else if (command == "CONVERT_GENBANK") {
		return do_convert_genbank(argc_new, argv_new);
	} else if (command == "CONTINGENCY_LOCI") {
		return do_contingency_loci(argc_new, argv_new);
	} else if (command == "ERROR_COUNT") {
		return do_error_count(argc_new, argv_new);
	} else if (command == "IDENTIFY_MUTATIONS") {
		return do_identify_mutations(argc_new, argv_new);
	} else if (command == "PREPROCESS_ALIGNMENTS") {
		return do_preprocess_alignments(argc_new, argv_new);
	} else if (command == "IDENTIFY_CANDIDATE_JUNCTIONS") {
		return do_identify_candidate_junctions(argc_new, argv_new);
	} else if (command == "TABULATE_COVERAGE") {
		return do_tabulate_coverage(argc_new, argv_new);
	} else if (command == "RESOLVE_ALIGNMENTS") {
		return do_resolve_alignments(argc_new, argv_new);
	} else if( command == "GD2GVF"){
        return do_convert_gvf(argc_new, argv_new);
    } else if( command == "VCF2GD" ){
        return do_convert_gd( argc_new, argv_new);
    }

	// Command was not recognized. Should output valid commands.
	cout << "Unrecognized command" << endl;
	return -1;
}
