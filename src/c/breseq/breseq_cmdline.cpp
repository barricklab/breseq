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
#include "breseq/mutation_predictor.h"


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

void convert_bull_form(const vector<string>& in, const string& fasta, const string& ft, const string& gff3 ) {
  cReferenceSequences refseqs;
  
  // Load the GenBank file
  LoadBullFile(refseqs, in);
  
  // Output sequence
  //if (fasta != "") refseqs.WriteFASTA(fasta);
  
  // Output feature table
  if (ft != "") refseqs.WriteFeatureTable(ft);
  
  //if (gff3 != "" ) refseqs.WriteGFF( gff3 );
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

int do_convert_bull_form(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq CONVERT_BULL_FORM --input <bull_form.txt> [--features <output.tab>]");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  ("input,i", "input bull form flatfile (multiple allowed, comma-separated)")
  ("features,g", "output feature table", "")
  ("fasta,f", "FASTA file of reference sequences", "")
  ("gff3,v", "GFF file of features", "" )
	.processCommandArgs(argc, argv);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("input")
		 || (!options.count("features") /*&& !options.count("fasta")*/)  
		 ) {
		options.printUsage();
		return -1;
	}
  
	// attempt to calculate error calibrations:
	try {
    
		convert_bull_form(  
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
		("no-junction-prediction,p", "whether to predict new junctions", TAKES_NO_ARGUMENT)
// convert to basing everything off the main output path, so we don't have to set so many options
    ("path", "path to breseq output")
		("readfile,r", "FASTQ read files (multiple allowed, comma-separated) ")
		("maximum-read-length,m", "number of flanking bases in candidate junctions")
		("alignment-read-limit", "maximum number of alignments to process. DEFAULT = 0 (OFF).", 0)
    ("junction-cutoff", "coverage cutoffs for different reference sequences")

	.processCommandArgs(argc, argv);

	// make sure that the config options are good:
	if(options.count("help")
     || !options.count("path")
     || !options.count("readfile")
		 || !options.count("maximum-read-length")
		 || !options.count("junction-cutoff")
		 ) {
		options.printUsage();
		return -1;
	}                       
  
	try {
    
    Settings settings(options["path"]);

    Summary summary;
    

    cReadFiles rf(from_string<vector<string> >(options["readfile"]));
    
    // Load the reference sequence info
    cReferenceSequences ref_seq_info;
    LoadFeatureIndexedFastaFile(ref_seq_info, settings.reference_features_file_name, settings.reference_fasta_file_name);
    
    // should be one coverage cutoff value for each reference sequence
    vector<double> coverage_cutoffs = from_string<vector<double> >(options["junction-cutoff"]);
    assert(coverage_cutoffs.size() == ref_seq_info.size());
    
    for (uint32_t i=0; i<ref_seq_info.size(); i++)
    {
      summary.preprocess_coverage[ref_seq_info[i].m_seq_id].junction_accept_score_cutoff = coverage_cutoffs[i];
    }
        
    resolve_alignments(
      settings,
      summary,
      ref_seq_info,
      !options.count("no-junction-prediction"),
      rf,
      from_string<uint32_t>(options["maximum-read-length"]),
      from_string<uint32_t>(options["alignment-read-limit"])
    );
    
  } catch(...) {
		// failed; 
    
		return -1;
	}
	
	return 0;
}


/*!  Predict Mutations
 
 Predict mutations from evidence in a genome diff file.
 
 */

int do_predict_mutations(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq PREDICT_MUTATIONS ... ");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  // convert to basing everything off the main output path, so we don't have to set so many options
  ("path", "path to breseq output")
  ("maximum-read-length,m", "number of flanking bases in candidate junctions")
	.processCommandArgs(argc, argv);
  
	// make sure that the config options are good:
	if(options.count("help")
     || !options.count("path")
		 || !options.count("maximum-read-length")
		 ) {
		options.printUsage();
		return -1;
	}                       
  
	try {
    
    Settings settings(options["path"]);
            
    // Load the reference sequence info
    cReferenceSequences ref_seq_info;
    LoadFeatureIndexedFastaFile(ref_seq_info, settings.reference_features_file_name, settings.reference_fasta_file_name);
        
    MutationPredictor mp(ref_seq_info);
    
    genome_diff gd( settings.evidence_genome_diff_file_name );
    gd.read();
    
    mp.predict(
               settings,
               gd,
               from_string<uint32_t>(options["maximum-read-length"])
               );
    
    gd.write(settings.final_genome_diff_file_name);
    
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
     "Maximum number of candidate junction to create.", static_cast<uint32_t>(5000))
    ("minimum-candidate-junctions",
     "Minimum number of candidate junctions to create.", static_cast<uint32_t>(200))
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

int do_read_gd( int argc, char* argv[]){
  AnyOption options("Usage: --gd <gd.gd>"); options( "gd,i","gd file to read").processCommandArgs( argc,argv);
  if( !options.count("gd")){
        options.printUsage(); return -1;
    }
    
    try{
        genome_diff( string(options["gd"])); ///TODO Why string()?
    } 
    catch(...){ 
        return -1; // failed 
    }
    
    return 0;
}

int breseq_default_action(int argc, char* argv[])
{
	// setup and parse configuration options:
	AnyOption options("Usage: breseq -r reference.gbk reads1.fastq [reads2.fastq, reads3.fastq...]");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		// convert to basing everything off the main output path, so we don't have to set so many options
		("path", "path to breseq output")
		("reference,r", "reference GenBank flatfile")
	.processCommandArgs(argc, argv);

	// make sure that the config options are good:
	if(options.count("help")
		|| !options.count("path")
		|| !options.count("reference")
		|| options.getArgc() < 1
	) {
		options.printUsage();
		return -1;
	}

	cReferenceSequences ref_seq_info;

	vector<string> read_file_names;
	for (int32_t i = 0; i < options.getArgc(); i++)
	{
		string read_file_name = options.getArgv(i);
		read_file_names.push_back(read_file_name);
		ref_seq_info.ReadFASTA(read_file_name);
	}

	vector<string> reference_file_names;
	reference_file_names.push_back(options["reference"]);
	LoadGenBankFile(ref_seq_info, reference_file_names);

	//---- Configuration Options ----

	// Keep a summary of certain statistics.
	Summary summary;

	///
	/// Get options from the command line
	///    handles all GetOpt and filling in many other settings
	///
	Settings settings(options["path"]);
	settings.check_installed();

	//
	// Convert the input reference GenBank into FASTA for alignment
	// sub sequence_conversion {}
	//

	string sequence_converson_summary_file_name = settings.file_name("sequence_conversion_summary_file_name");
	settings.create_path("sequence_conversion_path");
	settings.create_path("data_path");
/*
	if (settings.do_step("sequence_conversion_done_file_name", "Read and reference sequence file input"))
	{
		my $s;

		//Check the FASTQ format and collect some information about the input read files at the same time
		print STDERR "  Analyzing fastq read files...\n";
		my $overall_max_read_length;
		my $overall_max_qual = 0;

		$s->{num_reads} = 0;
		foreach my $read_file ($settings->read_files)
		{
			print STDERR "    READ FILE::$read_file\n";
			my $max_read_length;
			my $num_bases = 0;
			my $num_reads = 0;
			my $fastq_file_name = $settings->read_file_to_fastq_file_name($read_file);
			my $cbreseq = $settings->ctool("cbreseq");
			my $convert_file_name =  $settings->file_name("converted_fastq_file_name", {"//"=>$read_file});
			my $command = "$cbreseq ANALYZE_FASTQ --input=$fastq_file_name --convert=$convert_file_name";
			my $output = Breseq::Shared::capture_system($command);
			//print "$output\n";

			// Parse output
			my $s_rf;
			foreach my $line ( split(/\n/, $output) ) {
				if ($line =~ m/^(\S+)\s+(\S+)$/) {
					$s_rf->{$1} = $2;
				}
			}

			die if (!defined $s_rf->{max_read_length});
			die if (!defined $s_rf->{max_quality_score});
			die if (!defined $s_rf->{num_reads});
			die if (!defined $s_rf->{num_bases});
			die if (!defined $s_rf->{qual_format});
			die if (!defined $s_rf->{original_qual_format});
			die if (!defined $s_rf->{converted_fastq_name});

			// Save the converted file name -- have to save it in summary because only that
			// is reloaded if we skip this step.
			$s->{converted_fastq_name}->{$read_file} = $s_rf->{converted_fastq_name};

			// Record statistics
			$overall_max_read_length = $s_rf->{max_read_length} if ((!defined $overall_max_read_length) || ($s_rf->{max_read_length} > $overall_max_read_length));
			$overall_max_qual = $s_rf->{max_quality_score} if ((!defined $overall_max_qual) || ($s_rf->{max_quality_score} > $overall_max_qual));
			$s->{num_reads} += $s_rf->{num_reads};
			$s->{num_bases} += $s_rf->{num_bases};

			$s->{reads}->{$read_file} = $s_rf;
		}
		$s->{avg_read_length} = $s->{num_bases} / $s->{num_reads};
		$s->{max_read_length} = $overall_max_read_length;
		$s->{max_qual} = $overall_max_qual;
		$summary->{sequence_conversion} = $s;

		my @genbank_file_names = $settings->file_name("reference_genbank_file_names");
		my @junction_only_genbank_file_names = $settings->file_name("junction_only_reference_genbank_file_names");
		my $reference_fasta_file_name = $settings->file_name("reference_fasta_file_name");

		// C++ version of reference sequence processing
		my $reference_features_file_name = $settings->file_name("reference_features_file_name");
		my $reference_gff3_file_name = $settings->file_name("reference_gff3_file_name");
		my $cbreseq = $settings->ctool("cbreseq");
		my $command = "$cbreseq CONVERT_GENBANK --fasta $reference_fasta_file_name --features $reference_features_file_name --gff3 $reference_gff3_file_name";
		foreach my $input (@genbank_file_names) {
			$command .= " --input " . $input;
		}
		Breseq::Shared::system($command);

		// create SAM faidx
		my $samtools = $settings->ctool("samtools");
		Breseq::Shared::system("$samtools faidx $reference_fasta_file_name", 1);

		// calculate trim files
		my $output_path = $settings->file_name("sequence_conversion_path");
		$command = "$cbreseq CALCULATE_TRIMS -f $reference_fasta_file_name -o $output_path";
		Breseq::Shared::system($command);

		// store summary information
		Storable::store($summary->{sequence_conversion}, $sequence_converson_summary_file_name) or die "Can"t store data in file $sequence_converson_summary_file_name!\n";
		$settings->done_step("sequence_conversion_done_file_name");
	}

	$summary->{sequence_conversion} = Storable::retrieve($sequence_converson_summary_file_name);
	die "Can"t retrieve data from file $sequence_converson_summary_file_name!\n" if (!$summary->{sequence_conversion});
	(defined $summary->{sequence_conversion}->{max_read_length}) or die "Can"t retrieve max read length from file $sequence_converson_summary_file_name\n";

	//load C++ info
	my $reference_features_file_name = $settings->file_name("reference_features_file_name");
	my $reference_fasta_file_name = $settings->file_name("reference_fasta_file_name");
	$ref_seq_info = Breseq::ReferenceSequence::bridge_load_ref_seq_info($summary, $reference_features_file_name, $reference_fasta_file_name);
	die "Can"t retrieve reference sequence data from files $reference_features_file_name and $reference_fasta_file_name!\n" if (!$ref_seq_info);

	// reload certain information into $settings from $summary
	foreach my $read_file (keys %{$summary->{sequence_conversion}->{reads}})
	{
		if (defined $summary->{sequence_conversion}->{reads}->{$read_file}->{converted_fastq_name}) {
			$settings->{read_file_to_converted_fastq_file}->{$read_file} = $summary->{sequence_conversion}->{reads}->{$read_file}->{converted_fastq_name};
		}
	}
	$settings->{max_read_length} = $summary->{sequence_conversion}->{max_read_length};
	$settings->{max_smalt_diff} = floor($settings->{max_read_length}/2);
	$settings->{total_reference_sequence_length} = $summary->{sequence_conversion}->{total_reference_sequence_length};

	//load trim data per refseq
	foreach my $seq_id (sort keys %{$ref_seq_info->{ref_strings}})
	{
		my $trim_file_name = $settings->file_name("reference_trim_file_name", {"@"=>$seq_id});
		open TRIM, "<$trim_file_name" or die "Could not open trims file: $trim_file_name\n";
		$ref_seq_info->{trims}->{$seq_id} = "";
		my $block;
		my $block_size = 512;
		while (read TRIM, $block, $block_size)
		{
			$ref_seq_info->{trims}->{$seq_id} .= $block;
		}
		(length($ref_seq_info->{trims}->{$seq_id}) > 0) or die "Error reading trims file: $trim_file_name\n";
		close TRIM;

		//print "$seq_id " . length($ref_seq_info->{trims}->{$seq_id}) . "\n";
		//print Dumper($ref_seq_info->{trims}->{$seq_id});
	}

	//
	// Match all reads against the reference genome
	sub alignment_to_reference {}
	//

	if ($settings->do_step("reference_alignment_done_file_name", "Read alignment to reference genome"))
	{
		$settings->create_path("reference_alignment_path");

	//
	//nohup ssaha2 -save REL606_kmer1_skip1 -kmer 10 -skip 1 -output sam_soft -outfile JEB559_REL606_kmer1_skip1.sam JEB559.fastq &> JEB559_REL606_kmer1_skip1.out
	//nohup ssaha2 -save REL606_solexa -rtype solexa -output sam_soft -outfile JEB559_REL606_solexa.sam JEB559.fastq &> JEB559_REL606_solexa.out
	//
	//the lower seed value is important for finding split matches
	//nohup ssaha2 -save REL606_solexa -kmer 13 -skip 2 -seeds 1 -score 12 -cmatch 9 -ckmer 6 -output sam_soft -outfile JEB559_REL606_solexa_3.sam JEB559.fastq >&JEB559_REL606_solexa_3.out &

		/// create ssaha2 hash
		my $reference_hash_file_name = $settings->file_name("reference_hash_file_name");
		my $reference_fasta_file_name = $settings->file_name("reference_fasta_file_name");

		if (!$settings->{smalt})
		{
			Breseq::Shared::system("ssaha2Build -rtype solexa -skip 1 -save $reference_hash_file_name $reference_fasta_file_name");
		}
		else
		{
			my $smalt = $settings->ctool("smalt");
			Breseq::Shared::system("$smalt index -k 13 -s 1 $reference_hash_file_name $reference_fasta_file_name");
		}
		/// ssaha2 align reads to reference sequences
		foreach my $read_struct ($settings->read_structures)
		{
			//reads are paired
			if (defined $read_struct->{min_pair_dist} && defined $read_struct->{max_pair_dist})
			{
				// JEB this is not working currently
				die "Paired end mapping is broken.";
				die if (scalar @{$read_struct->{read_fastq_list}} != 2);

				my $fastq_1 = $read_struct->{read_fastq_list}->[0];
				my $fastq_2 = $read_struct->{read_fastq_list}->[1];
				my $min = $read_struct->{min_pair_dist};
				my $max = $read_struct->{max_pair_dist};

				my $reference_sam_file_name = $settings->file_name("reference_sam_file_name", {"//"=>$read_struct->{base_name}});
				Breseq::Shared::system("ssaha2 -save $reference_hash_file_name -kmer 13 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 1 -output sam_soft -outfile $reference_sam_file_name -multi 1 -mthresh 9 -pair $min,$max $fastq_1 $fastq_2");
			}

			//reads are not paired
			else
			{
				die if (scalar @{$read_struct->{base_names}} != 1);
				my $read_name = $read_struct->{base_names}->[0];
				my $read_fastq_file = $settings->read_file_to_fastq_file_name($read_name);
				my $reference_sam_file_name = $settings->file_name("reference_sam_file_name", {"//"=>$read_name});

				if (!$settings->{smalt})
				{
					Breseq::Shared::system("ssaha2 -save $reference_hash_file_name -kmer 13 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 1 -output sam_soft -outfile $reference_sam_file_name $read_fastq_file");
				}
				else
				{
					my $smalt = $settings->ctool("smalt");
					Breseq::Shared::system("$smalt map -n 2 -d $settings->{max_smalt_diff} -f samsoft -o $reference_sam_file_name $reference_hash_file_name $read_fastq_file");
					// -m 12
				}
			}
		}

		/// Delete the hash files immediately
		if (!$settings->{keep_all_intermediates})
		{
			unlink "$reference_hash_file_name.base";
			unlink "$reference_hash_file_name.body";
			unlink "$reference_hash_file_name.head";
			unlink "$reference_hash_file_name.name";
			unlink "$reference_hash_file_name.size";
		}

		$settings->done_step("reference_alignment_done_file_name");
	}

	//
	// Identify candidate junctions from split read alignments
	sub identify_candidate_junctions {}
	//


	if (!$settings->{no_junction_prediction})
	{
		$settings->create_path("candidate_junction_path");

		if ( (defined $settings->{preprocess_junction_min_indel_split_length}) || ($settings->{candidate_junction_score_method} eq "POS_HASH"))
		{
			my $preprocess_junction_done_file_name = $settings->file_name("preprocess_junction_done_file_name");

			if ($settings->do_step("preprocess_junction_done_file_name", "Preprocessing alignments for candidate junction identification"))
			{
				my $cbreseq = $settings->ctool("cbreseq");
				my $readfiles = join(" --read-file ", $settings->read_files);
				my $data_path = $settings->file_name("data_path");
				my $reference_alignment_path = $settings->file_name("reference_alignment_path");
				my $candidate_junction_path = $settings->file_name("candidate_junction_path");

				my $cmdline = "$cbreseq PREPROCESS_ALIGNMENTS --data-path $data_path";
				$cmdline .= " --reference-alignment-path $reference_alignment_path";
				$cmdline .= " --candidate-junction-path $candidate_junction_path";
				$cmdline .= " --read-file $readfiles";

				Breseq::Shared::system($cmdline);

				$settings->done_step("preprocess_junction_done_file_name");
			}
		}

		if ($settings->{candidate_junction_score_method} eq "POS_HASH")
		{
			my $coverage_junction_summary_file_name = $settings->file_name("coverage_junction_summary_file_name");

			if ($settings->do_step("coverage_junction_done_file_name", "Preliminary analysis of coverage distribution"))
			{
				my $reference_faidx_file_name = $settings->file_name("reference_faidx_file_name");
				my $preprocess_junction_best_sam_file_name = $settings->file_name("preprocess_junction_best_sam_file_name");
				my $coverage_junction_best_bam_file_name = $settings->file_name("coverage_junction_best_bam_file_name");
				my $coverage_junction_best_bam_prefix = $settings->file_name("coverage_junction_best_bam_prefix");
				my $coverage_junction_best_bam_unsorted_file_name = $settings->file_name("coverage_junction_best_bam_unsorted_file_name");

				my $samtools = $settings->ctool("samtools");

				Breseq::Shared::system("$samtools import $reference_faidx_file_name $preprocess_junction_best_sam_file_name $coverage_junction_best_bam_unsorted_file_name");
				Breseq::Shared::system("$samtools sort $coverage_junction_best_bam_unsorted_file_name $coverage_junction_best_bam_prefix");
				unlink $coverage_junction_best_bam_unsorted_file_name if (!$settings->{keep_all_intermediates});
				Breseq::Shared::system("$samtools index $coverage_junction_best_bam_file_name");

				// Count errors
				my $reference_fasta_file_name = $settings->file_name("reference_fasta_file_name");
				my $reference_bam_file_name = $settings->file_name("coverage_junction_best_bam_file_name");

				my $cbreseq = $settings->ctool("cbreseq");
				// deal with distribution or error count keys being undefined...
				my $coverage_fn = $settings->file_name("coverage_junction_distribution_file_name", {"@"=>""});
				my $outputdir = `dirname $coverage_fn`;
				chomp $outputdir; $outputdir .= "/";
				my $readfiles = join(" --readfile ", $settings->read_files);
				my $cmdline = "$cbreseq ERROR_COUNT --bam $reference_bam_file_name --fasta $reference_fasta_file_name --output $outputdir --readfile $readfiles --coverage";
				// IMPORTANT: WE DON"T COUNT ERRORS HERE, ONLY COVERAGE @JEB
				//$cmdline .= " --errors";
				$cmdline .= " --minimum-quality-score $settings->{base_quality_cutoff}" if ($settings->{base_quality_cutoff});
				Breseq::Shared::system($cmdline);

				my $error_rates_summary_file_name = $settings->file_name("error_rates_summary_file_name");
				Breseq::CoverageDistribution::analyze_unique_coverage_distributions($settings, $summary, $ref_seq_info,
					"coverage_junction_plot_file_name", "coverage_junction_distribution_file_name");

				Storable::store($summary->{unique_coverage}, $coverage_junction_summary_file_name) or die "Can"t store data in file $coverage_junction_summary_file_name!\n";
				$settings->done_step("coverage_junction_done_file_name");
			}

			$summary->{preprocess_coverage} = Storable::retrieve($coverage_junction_summary_file_name);
			die "Can"t retrieve data from file $coverage_junction_summary_file_name!\n" if (!$summary->{preprocess_coverage});
		}

		my $candidate_junction_summary_file_name = $settings->file_name("candidate_junction_summary_file_name");
		if ($settings->do_step("candidate_junction_done_file_name", "Identifying candidate junctions"))
		{
			print STDERR "Identifying candidate junctions...\n";

			my $cbreseq = $settings->ctool("cbreseq");
			my $readfiles = join(" --read-file ", $settings->read_files);
			my $data_path = $settings->file_name("data_path");
			my $candidate_junction_path = $settings->file_name("candidate_junction_path");

			my $cmdline = "$cbreseq IDENTIFY_CANDIDATE_JUNCTIONS --data-path $data_path";
			$cmdline .= " --candidate-junction-path $candidate_junction_path";
			$cmdline .= " --read-file $readfiles";
			my $maximum_read_length = $summary->{sequence_conversion}->{max_read_length};
			$cmdline .= " --maximum-read-length $maximum_read_length";
			my $reference_sequence_length = $settings->{total_reference_sequence_length};
			$cmdline .= " --reference-sequence-length $reference_sequence_length";

			Breseq::Shared::system($cmdline);

			my $samtools = $settings->ctool("samtools");
			my $faidx_command = "$samtools faidx $candidate_junction_path/candidate_junction.fasta";
			if (-s "$candidate_junction_path/candidate_junction.fasta" > 0)
			{
				Breseq::Shared::system($faidx_command);
			}
			// @JEB Fix this -- no summary currently...
			$summary->{candidate_junction} = {};

			Storable::store($summary->{candidate_junction}, $candidate_junction_summary_file_name)
				or die "Can"t store data in file $candidate_junction_summary_file_name!\n";

			Breseq::Output::record_time("Candidate junction identification");

			$settings->done_step("candidate_junction_done_file_name");
		}

		//load this info
		$summary->{candidate_junction} = Storable::retrieve($candidate_junction_summary_file_name);
		die "Can"t retrieve data from file $candidate_junction_summary_file_name!\n" if (!$summary->{candidate_junction});


	//
	// Find matches to new junction candidates
	sub candidate_junction_alignment {}
	//

		if ($settings->do_step("candidate_junction_alignment_done_file_name", "Candidate junction alignment"))
		{
			$settings->create_path("candidate_junction_alignment_path");

			/// create ssaha2 hash
			my $candidate_junction_hash_file_name = $settings->file_name("candidate_junction_hash_file_name");
			my $candidate_junction_fasta_file_name = $settings->file_name("candidate_junction_fasta_file_name");

			if (-s $candidate_junction_fasta_file_name > 0)
			{
				if (!$settings->{smalt})
				{
					Breseq::Shared::system("ssaha2Build -rtype solexa -skip 1 -save $candidate_junction_hash_file_name $candidate_junction_fasta_file_name");
				}
				else
				{
					my $smalt = $settings->ctool("smalt");
					Breseq::Shared::system("$smalt index -k 13 -s 1 $candidate_junction_hash_file_name $candidate_junction_fasta_file_name");
				}
			}


			/// ssaha2 align reads to candidate junction sequences

			foreach my $read_name ($settings->read_files)
			{
				my $candidate_junction_sam_file_name = $settings->file_name("candidate_junction_sam_file_name", {"//"=>$read_name});

				my $read_fastq_file = $settings->read_file_to_fastq_file_name($read_name);
				print "$read_fastq_file\n";

				if (!$settings->{smalt} && (-e "$candidate_junction_hash_file_name.base"))
				{
					Breseq::Shared::system("ssaha2 -save $candidate_junction_hash_file_name -best 1 -rtype solexa -skip 1 -seeds 1 -output sam_soft -outfile $candidate_junction_sam_file_name $read_fastq_file");
					// Note: Added -best parameter to try to avoid too many matches to redundant junctions!
				}
				elsif (-e "$candidate_junction_hash_file_name.sma")
				{
					my $smalt = $settings->ctool("smalt");
					Breseq::Shared::system("$smalt map -c 0.8 -x -n 2 -d 1 -f samsoft -o $candidate_junction_sam_file_name $candidate_junction_hash_file_name $read_fastq_file");
					//-m 12
				}
			}

			/// Delete the hash files immediately
			if (!$settings->{keep_all_intermediates})
			{
				unlink "$candidate_junction_hash_file_name.base";
				unlink "$candidate_junction_hash_file_name.body";
				unlink "$candidate_junction_hash_file_name.head";
				unlink "$candidate_junction_hash_file_name.name";
				unlink "$candidate_junction_hash_file_name.size";
			}

			$settings->done_step("candidate_junction_alignment_done_file_name");
		}
	}

	//
	// Resolve matches to new junction candidates
	sub alignment_correction {}
	//

	my @hybrids;
	my $alignment_correction_summary_file_name = $settings->file_name("alignment_correction_summary_file_name");
	if ($settings->do_step("alignment_correction_done_file_name", "Resolving alignments with candidate junctions"))
	{
		$settings->create_path("alignment_correction_path");

		my $cbreseq = $settings->ctool("cbreseq");

		my $cmdline = "$cbreseq RESOLVE_ALIGNMENTS";
		my $base_output_path = $settings->file_name("base_output_path");
		$cmdline .= " --path $base_output_path";

		//print Dumper($settings);
		// complicated because we need full paths to original or converted files

		foreach my $read_file ($settings->read_files)
		{
			my $fastq_file = $settings->read_file_to_fastq_file_name($read_file);
			$cmdline .= " --readfile $fastq_file";
		}

		my $maximum_read_length = $summary->{sequence_conversion}->{max_read_length};
		$cmdline .= " --maximum-read-length $maximum_read_length";

		foreach my $seq_id (@{$ref_seq_info->{seq_ids}})
		{
			my $junction_cutoff = $summary->{preprocess_coverage}->{$seq_id}->{junction_accept_score_cutoff};
			$cmdline .= " --junction-cutoff $junction_cutoff";
		}

		Breseq::Shared::system($cmdline);

		// need to fill this in
		$summary->{alignment_correction} = {};

		my $alignment_correction_summary_file_name = $settings->file_name("alignment_correction_summary_file_name");
		Storable::store($summary->{alignment_correction}, $alignment_correction_summary_file_name)
			or die "Can"t store data in file $alignment_correction_summary_file_name!\n";
		Breseq::Output::record_time("Resolve candidate junctions");
		$settings->done_step("alignment_correction_done_file_name");
	}
	$summary->{alignment_correction} = Storable::retrieve($alignment_correction_summary_file_name) if (-e $alignment_correction_summary_file_name);


	//
	// Create BAM files
	sub bam_creation {}
	//

	if ($settings->do_step("bam_done_file_name", "Creating BAM files"))
	{
		$settings->create_path("bam_path");

		my $reference_faidx_file_name = $settings->file_name("reference_faidx_file_name");
		my $candidate_junction_faidx_file_name = $settings->file_name("candidate_junction_faidx_file_name");

		my $resolved_junction_sam_file_name = $settings->file_name("resolved_junction_sam_file_name");
		my $junction_bam_unsorted_file_name = $settings->file_name("junction_bam_unsorted_file_name");
		my $junction_bam_prefix = $settings->file_name("junction_bam_prefix");
		my $junction_bam_file_name = $settings->file_name("junction_bam_file_name");

		my $samtools = $settings->ctool("samtools");

		if (!$settings->{no_junction_prediction})
		{
			Breseq::Shared::system("$samtools import $candidate_junction_faidx_file_name $resolved_junction_sam_file_name $junction_bam_unsorted_file_name");
			Breseq::Shared::system("$samtools sort $junction_bam_unsorted_file_name $junction_bam_prefix");
			unlink $junction_bam_unsorted_file_name if (!$settings->{keep_all_intermediates});
			Breseq::Shared::system("$samtools index $junction_bam_file_name");
		}

		my $resolved_reference_sam_file_name = $settings->file_name("resolved_reference_sam_file_name");
		my $reference_bam_unsorted_file_name = $settings->file_name("reference_bam_unsorted_file_name");
		my $reference_bam_prefix = $settings->file_name("reference_bam_prefix");
		my $reference_bam_file_name = $settings->file_name("reference_bam_file_name");

		Breseq::Shared::system("$samtools import $reference_faidx_file_name $resolved_reference_sam_file_name $reference_bam_unsorted_file_name");
		Breseq::Shared::system("$samtools sort $reference_bam_unsorted_file_name $reference_bam_prefix");
		unlink $reference_bam_unsorted_file_name if (!$settings->{keep_all_intermediates});
		Breseq::Shared::system("$samtools index $reference_bam_file_name");

		// delete unneeded files
		unlink $reference_bam_unsorted_file_name;

		$settings->done_step("bam_done_file_name");
	}*/

	//
	//  Graph paired read outliers (experimental)
	// sub paired_read_distances {}
	//
	// {
	// 	my @rs = $settings->read_structures;
	//
	// 	my @min_pair_dist;
	// 	my @max_pair_dist;
	//
	// 	my $paired = 0;
	//
	// 	my $i=0;
	// 	foreach my $rfi (@{$settings->{read_file_index_to_struct_index}})
	// 	{
	// 		$min_pair_dist[$i] = 0;
	// 		$max_pair_dist[$i] = 0;
	//
	// 		if ($rs[$rfi]->{paired})
	// 		{
	// 			$paired = 1;
	// 			$min_pair_dist[$i] = $rs[$rfi]->{min_pair_dist};
	// 			$max_pair_dist[$i] = $rs[$rfi]->{max_pair_dist};
	// 		}
	// 		$i++;
	// 	}
	//
	// 	my $long_pairs_file_name = $settings->file_name("long_pairs_file_name");
	//
	// 	if ($paired && (!-e $long_pairs_file_name))
	// 	{
	//
	// 		my $reference_sam_file_name = $settings->file_name("resolved_reference_sam_file_name");
	// 		my $reference_tam = Bio::DB::Tam->open($reference_sam_file_name) or die "Could not open $reference_sam_file_name";
	//
	// 		my $reference_faidx_file_name = $settings->file_name("reference_faidx_file_name");
	// 		my $reference_header = $reference_tam->header_read2($reference_faidx_file_name) or throw("Error reading reference fasta index file: $reference_faidx_file_name");
	// 		my $target_names = $reference_header->target_name;
	//
	// 		my $save;
	// 		my $on_alignment = 0;
	// 		my $last;
	//
	// 		while (1)
	// 		{
	// 			$a = Bio::DB::Bam::Alignment->new();
	// 			my $bytes = $reference_tam->read1($reference_header, $a);
	// 			last if ($bytes <= 0);
	//
	//
	// 			my $start       = $a->start;
	// 		    my $end         = $a->end;
	// 		    my $seqid       = $target_names->[$a->tid];
	//
	// 			$on_alignment++;
	// 			print "$on_alignment\n" if ($on_alignment % 10000 == 0);
	//
	// 			//last if ($on_alignment > 100000);
	//
	// 			//print $a->qname . "\n";
	//
	// 			if (!$a->unmapped)
	// 			{
	// 				my $mate_insert_size = abs($a->isize);
	// 				my $mate_end = $a->mate_end;
	// 				my $mate_start = $a->mate_start;
	// 				my $mate_reversed = 2*$a->mreversed + $a->reversed;
	// 		 		my $mreversed = $a->mreversed;
	// 		 		my $reversed = $a->reversed;
	//
	// 				my $fastq_file_index = $a->aux_get("X2");
	// 				//print "$mate_insert_size $min_pair_dist[$fastq_file_index] $max_pair_dist[$fastq_file_index]\n";
	// 				//if (($mate_insert_size < $min_pair_dist[$fastq_file_index]) || ($mate_insert_size > $max_pair_dist[$fastq_file_index]))
	// 				if ((($mate_insert_size >= 400) && ($mate_insert_size < $min_pair_dist[$fastq_file_index])) || ($mate_insert_size > $max_pair_dist[$fastq_file_index]))
	// 				{
	// 					//correct pair
	//
	// 					if ($last && ($last->{start} == $mate_start))
	// 					{
	// 						$save->{int($start/100)}->{int($mate_start/100)}->{$mate_reversed}++;
	// 						$save->{int($last->{start}/100)}->{int($last->{mate_start}/100)}->{$last->{mate_reversed}}++;
	// 						undef $last;
	// 					}
	// 					else
	// 					{
	// 						($last->{mate_reversed}, $last->{start}, $last->{mate_start}) = ($mate_reversed, $start, $mate_start);
	// 					}
	//
	// 					//$save->{$mate_reversed}->{int($start/100)}->{int($mate_start/100)}++;
	// 				    //print $a->qname," aligns to $seqid:$start..$end, $mate_start $mate_reversed ($mreversed $reversed) $mate_insert_size\n";
	// 				}
	//
	// 			}
	// 		}
	//
	// 		open LP, ">$long_pairs_file_name" or die;
	//
	// 		foreach my $key_1 (sort {$a <=> $b} keys %$save)
	// 		{
	// 			foreach my $key_2 (sort {$a <=> $b} keys %{$save->{$key_1}})
	// 			{
	// 				foreach my $key_reversed (sort {$a <=> $b} keys %{$save->{$key_1}->{$key_2}})
	// 				{
	// 					print LP "$key_1\t$key_2\t$key_reversed\t$save->{$key_1}->{$key_2}->{$key_reversed}\n";
	// 				}
	// 			}
	// 		}
	// 		close LP;
	// 	}
	//
	// 	if ($paired)
	// 	{
	// 		open LP, "$long_pairs_file_name" or die;
	// 		while ($_ = <LP>)
	// 		{
	// 			chomp $_;
	// 			my ($start, $end, $key_reversed);
	// 		}
	// 	}
	// }
	//

	//
	// Tabulate error counts and coverage distribution at unique only sites
	//sub error_count {}
	//

	/*if ($settings->do_step("error_counts_done_file_name", "Tabulating error counts"))
	{
		$settings->create_path("error_calibration_path");

		my $reference_fasta_file_name = $settings->file_name("reference_fasta_file_name");
		my $reference_bam_file_name = $settings->file_name("reference_bam_file_name");

		my $cbreseq = $settings->ctool("cbreseq");

		// deal with distribution or error count keys being undefined...
		my $coverage_fn = $settings->file_name("unique_only_coverage_distribution_file_name", {"@"=>""});
		my $outputdir = `dirname $coverage_fn`;
		chomp $outputdir; $outputdir .= "/";
		my $readfiles = join(" --readfile ", $settings->read_files);
		my $cmdline = "$cbreseq ERROR_COUNT --bam $reference_bam_file_name --fasta $reference_fasta_file_name --output $outputdir --readfile $readfiles --coverage";
		$cmdline .= " --errors";
		$cmdline .= " --minimum-quality-score $settings->{base_quality_cutoff}" if ($settings->{base_quality_cutoff});

	//----->//this line should be calculated
		my $num_read_files = scalar(keys %{$summary->{sequence_conversion}->{reads}});
		my $num_qual = $summary->{sequence_conversion}->{max_qual} + 1;
		$cmdline .= " --covariates=read_set=$num_read_files,obs_base,ref_base,quality=$num_qual";
	//		$cmdline .= " --covariates=read_set=$num_read_files,obs_base,ref_base,quality=$summary->{sequence_conversion}->{max_qual},base_repeat=5";
		Breseq::Shared::system($cmdline);

		$settings->done_step("error_counts_done_file_name");
	}


	//
	// Calculate error rates
	sub error_rates {}
	//

	$settings->create_path("output_path"); //need output for plots
	my $error_rates_summary_file_name = $settings->file_name("error_rates_summary_file_name");

	if ($settings->do_step("error_rates_done_file_name", "Re-calibrating base error rates"))
	{
		$summary->{unique_coverage} = {};
		if (!$settings->{no_deletion_prediction}) {
			Breseq::CoverageDistribution::analyze_unique_coverage_distributions($settings, $summary, $ref_seq_info,
				"unique_only_coverage_plot_file_name", "unique_only_coverage_distribution_file_name");
		}

		foreach my $read_file ($settings->read_files) {
			my $error_rates_base_qual_error_prob_file_name = $settings->file_name("error_rates_base_qual_error_prob_file_name", {"//" => $read_file});
			my $plot_error_rates_r_script_file_name = $settings->file_name("plot_error_rates_r_script_file_name");
			my $plot_error_rates_r_script_log_file_name = $settings->file_name("plot_error_rates_r_script_log_file_name", {"//" => $read_file});
			my $error_rates_plot_file_name = $settings->file_name("error_rates_plot_file_name", {"//" => $read_file});
			Breseq::Shared::system("R --vanilla in_file=$error_rates_base_qual_error_prob_file_name out_file=$error_rates_plot_file_name < $plot_error_rates_r_script_file_name > $plot_error_rates_r_script_log_file_name");
		}

		Storable::store($summary->{unique_coverage}, $error_rates_summary_file_name) or die "Can"t store data in file $error_rates_summary_file_name!\n";
		$settings->done_step("error_rates_done_file_name");
	}
	$summary->{unique_coverage} = Storable::retrieve($error_rates_summary_file_name);
	die "Can"t retrieve data from file $error_rates_summary_file_name!\n" if (!$summary->{unique_coverage});
	//these are determined by the loaded summary information
	$settings->{unique_coverage} = $summary->{unique_coverage};

	//
	// Make predictions of point mutations, small indels, and large deletions
	sub mutation_prediction {}
	//

	my @mutations;
	my @deletions;
	my @unknowns;

	if (!$settings->{no_mutation_prediction})
	{
		$settings->create_path("mutation_identification_path");

		if ($settings->do_step("mutation_identification_done_file_name", "Read alignment mutations"))
		{
			my $error_rates;

			my $reference_fasta_file_name = $settings->file_name("reference_fasta_file_name");
			my $reference_bam_file_name = $settings->file_name("reference_bam_file_name");

			my $cbreseq = $settings->ctool("cbreseq");

			my $coverage_fn = $settings->file_name("unique_only_coverage_distribution_file_name", {"@"=>""});
			my $error_dir = `dirname $coverage_fn`;
			chomp $error_dir; $error_dir .= "/";
			my $this_predicted_mutation_file_name = $settings->file_name("predicted_mutation_file_name", {"@"=>""});
			my $output_dir = `dirname $this_predicted_mutation_file_name`;
			chomp $output_dir; $output_dir .= "/";
			my $readfiles = join(" --readfile ", $settings->read_files);
			my $cmdline = "$cbreseq IDENTIFY_MUTATIONS --bam $reference_bam_file_name --fasta $reference_fasta_file_name --readfile $readfiles";
			$cmdline .= " --error_dir $error_dir";
			my $ra_mc_genome_diff_file_name = $settings->file_name("ra_mc_genome_diff_file_name");
			$cmdline .= " --genome_diff $ra_mc_genome_diff_file_name";
			$cmdline .= " --output $output_dir";
			if(defined $settings->{mutation_log10_e_value_cutoff}) {
				$cmdline .= " --mutation_cutoff $settings->{mutation_log10_e_value_cutoff}"; // defaults to 2.0.
			}
			my $coverage_tab_file_name = $settings->file_name("complete_coverage_text_file_name", {"@"=>""});
			my $coverage_dir = `dirname $coverage_tab_file_name`;
			chomp $coverage_dir; $coverage_dir .= "/";
			$cmdline .= " --coverage_dir $coverage_dir";

			// It is important that these are in consistent order with the fasta file!!
			foreach my $seq_id ( @{$ref_seq_info->{seq_ids}})
			{
				my $deletion_propagation_cutoff = $settings->{unique_coverage}->{$seq_id}->{deletion_coverage_propagation_cutoff};
				$cmdline .= " --deletion_propagation_cutoff $deletion_propagation_cutoff";
			}

			if((defined $settings->{no_deletion_prediction}) && ($settings->{no_deletion_prediction})) {
				$cmdline .= "";
			} else {
				$cmdline .= " --predict_deletions"; // defaults TO predicting deletions.
			}

			if (defined $settings->{base_quality_cutoff})
			{
				$cmdline .= " --minimum_quality_score $settings->{base_quality_cutoff}";
			}

			if ($settings->{polymorphism_prediction})
			{
				$cmdline .= " --predict_polymorphisms";
				$cmdline .= " --polymorphism_cutoff $settings->{polymorphism_log10_e_value_cutoff}";
				$cmdline .= " --polymorphism_frequency_cutoff $settings->{polymorphism_frequency_cutoff}";
			}

			$cmdline .= " --error_table $error_dir/error_rates.tab";

			Breseq::Shared::system($cmdline);

			$settings->done_step("mutation_identification_done_file_name");
		}

		my $polymorphism_statistics_done_file_name = $settings->file_name("polymorphism_statistics_done_file_name");
		if ($settings->{polymorphism_prediction} && $settings->do_step("polymorphism_statistics_done_file_name", "Polymorphism statistics"))
		{
			Breseq::Shared::polymorphism_statistics($settings, $summary, $ref_seq_info);
			$settings->done_step("polymorphism_statistics_done_file_name");
		}
	}

	//rewire which GenomeDiff we get data from if we have the elaborated polymorphism_statistics version
	$settings->{ra_mc_genome_diff_file_name} = $settings->{polymorphism_statistics_ra_mc_genome_diff_file_name} if ($settings->{polymorphism_prediction});

	if ($settings->do_step("output_done_file_name", "Output"))
	{
		///
		// Output Genome Diff File
		sub genome_diff_output {}
		///
		print STDERR "Creating merged genome diff evidence file...\n";

		// merge all of the evidence GenomeDiff files into one...
		$settings->create_path("evidence_path");
		my $jc_genome_diff_file_name = $settings->file_name("jc_genome_diff_file_name");
		my $jc_gd = GenomeDiff->new( {in => $jc_genome_diff_file_name} );
		my $ra_mc_genome_diff_file_name = $settings->file_name("ra_mc_genome_diff_file_name");
		my $ra_mc_gd = GenomeDiff->new( {in => $ra_mc_genome_diff_file_name} );
		my $evidence_genome_diff_file_name = $settings->file_name("evidence_genome_diff_file_name");
		my $evidence_gd = GenomeDiff::merge($jc_gd, $ra_mc_gd);
		$evidence_gd->write($evidence_genome_diff_file_name);


		// predict mutations from evidence in the GenomeDiff
		print STDERR "Predicting mutations from evidence...\n";

		my $gd;
		my $final_genome_diff_file_name = $settings->file_name("final_genome_diff_file_name");

		my $cbreseq = $settings->ctool("cbreseq");
		my $maximum_read_length = $summary->{sequence_conversion}->{max_read_length};
		my $base_output_path = $settings->file_name("base_output_path");
		my $cmdline = "$cbreseq PREDICT_MUTATIONS --maximum-read-length $maximum_read_length --path $base_output_path";
		Breseq::Shared::system($cmdline);
		$gd = GenomeDiff->new({file_name => $final_genome_diff_file_name});
		//unlink $evidence_genome_diff_file_name;

		sub mutation_annotation {}

		my @genbank_file_names = $settings->file_name("reference_genbank_file_names");
		my @junction_only_genbank_file_names = $settings->file_name("junction_only_reference_genbank_file_names");


		//
		// Annotate mutations
		//
		print STDERR "Annotating mutations...\n";
		Breseq::ReferenceSequence::annotate_mutations($ref_seq_info, $gd);

		//
		// Plot coverage of genome and large deletions
		//
		print STDERR "Drawing coverage plots...\n";
		Breseq::Output::draw_coverage($settings, $ref_seq_info, $gd);

		//
		// Mark lowest RA evidence items as no-show, or we may be drawing way too many alignments
		//

		my @ra = $gd->filter_used_as_evidence($gd->list("RA"));

		@ra = grep { ($_->{frequency} != 0) && ($_->{frequency} != 1) && (!$_->{no_show})  } @ra;
		@ra = sort { -($a->{quality} <=> $b->{quality}) } @ra;

		for (my $i = $settings->{max_rejected_polymorphisms_to_show}; $i< scalar @ra; $i++)
		{
			$ra[$i]->{no_show} = 1;
		}

		// require a certain amount of coverage
		foreach my $item ($gd->filter_used_as_evidence($gd->list("RA")))
		{
			my ($top, $bot) = split /\//, $item->{tot_cov};
			$item->{no_show} = 1 if ($top + $bot <= 2);
		}
		@ra = grep { !$_->{coverage} && !$_->{no_show} } @ra;

		//
		// Mark lowest scoring reject junctions as no-show
		//
		my @jc = $gd->filter_used_as_evidence($gd->list("JC"));
		@jc = grep { $_->{reject} } @jc;

		@jc = sort { -($a->{pos_hash_score} <=> $b->{pos_hash_score}) || -($a->{min_overlap_score} <=> $b->{min_overlap_score})  || ($a->{total_reads} <=> $a->{total_reads}) } @jc;
		for (my $i = $settings->{max_rejected_junctions_to_show}; $i< scalar @jc; $i++)
		{
			$jc[$i]->{no_show} = 1;
		}

		//
		// Create evidence files containing alignments and coverage plots
		//
		if (!$settings->{no_alignment_generation})
		{
			Breseq::Output::create_evidence_files($settings, $gd);
		}

		///
		// HTML output
		///

		print STDERR "Creating index HTML table...\n";

		my $index_html_file_name = $settings->file_name("index_html_file_name");
		Breseq::Output::html_index($index_html_file_name, $settings, $summary, $ref_seq_info, $gd);

		my $marginal_html_file_name = $settings->file_name("marginal_html_file_name");
		Breseq::Output::html_marginal_predictions($marginal_html_file_name, $settings, $summary, $ref_seq_info, $gd);

		///
		// Temporary debug output using Data::Dumper
		///

		my $summary_text_file_name = $settings->file_name("summary_text_file_name");
		open SUM, ">$summary_text_file_name";
		print SUM Dumper($summary);
		close SUM;

		my $settings_text_file_name = $settings->file_name("settings_text_file_name");
		open SETTINGS, ">$settings_text_file_name";
		print SETTINGS Dumper($settings);
		close SETTINGS;

		// record the final time and print summary table
		$settings->record_end_time("Output");

		Breseq::Output::html_statistics($settings->{summary_html_file_name}, $settings, $summary, $ref_seq_info);

		$settings->done_step("output_done_file_name");
	}*/
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
  } else if (command == "CONVERT_BULL_FORM") {
    return do_convert_bull_form(argc_new, argv_new);
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
  } else if (command == "PREDICT_MUTATIONS") {
		return do_predict_mutations(argc_new, argv_new);
	} else if (command == "GD2GVF") {
    return do_convert_gvf(argc_new, argv_new);
  } else if (command == "VCF2GD") {
    return do_convert_gd( argc_new, argv_new);
  } else if (command == "READGD") {
    return do_read_gd(argc_new, argv_new);
  }

	// Command was not recognized. Should output valid commands.
	cout << "Unrecognized command" << endl;
	return -1;
}
