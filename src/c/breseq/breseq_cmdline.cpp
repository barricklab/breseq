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

#include "libbreseq/anyoption.h"
#include "libbreseq/alignment_output.h"
#include "libbreseq/annotated_sequence.h"
#include "libbreseq/calculate_trims.h"
#include "libbreseq/candidate_junctions.h"
#include "libbreseq/contingency_loci.h"
#include "libbreseq/coverage_distribution.h"
#include "libbreseq/coverage_output.h"
#include "libbreseq/error_count.h"
#include "libbreseq/fastq.h"
#include "libbreseq/genome_diff.h"
#include "libbreseq/identify_mutations.h"
#include "libbreseq/resolve_alignments.h"
#include "libbreseq/settings.h"
#include "libbreseq/contingency_loci.h"
#include "libbreseq/mutation_predictor.h"
#include "libbreseq/output.h"

using namespace breseq;
using namespace std;


/*! bam2aln
 Draw HTML alignment from BAM
 */
int do_bam2aln(int argc, char* argv[]) {
    // setup and parse configuration options:
	AnyOption options("Usage: bam2aln --bam=<reference.bam> --fasta=<reference.fasta> --region=<accession:start-end> --output=<output.html> [--max-reads=1000]");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  ("bam,b", "bam file containing sequences to be aligned", "data/reference.bam")
	("fasta,f", "FASTA file of reference sequence", "data/reference.fasta")
  ("output,o", "output to file [region.html]")
  ("region,r", "region to print (accession:start-end)", "")
  ("max-reads,n", "maximum number of reads to show in alignment", 200)
  ("quality-score-cutoff,c", "quality score cutoff", 0)
  ("stdout", "write output to stdout", TAKES_NO_ARGUMENT)
  .processCommandArgs(argc, argv);
  
  
  
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("region")
     || !file_exists(options["fasta"].c_str())
     || !file_exists(options["bam"].c_str()) )
  {
		options.printUsage();
		return -1;
	}
  
  // generate alignment!
  alignment_output ao(
                      options["bam"],
                      options["fasta"],
                      from_string<uint32_t>(options["max-reads"]),
                      from_string<uint32_t>(options["quality-score-cutoff"])
                      );
  
  string html_output = ao.html_alignment(options["region"]);
  //cout << html_output << endl;
  
  if (options.count("stdout"))
  {
    cout << html_output << endl;
  }
  else
  {
    ///Write to html file
    string file_name = options["region"] + ".html";
    if (options.count("output")) {
      file_name = options["output"];
    }
    
    ofstream myfile (file_name.c_str());
    if (myfile.is_open())
    {
      myfile << html_output;
      myfile.close();
    }
    else cerr << "Unable to open file";
  }
    
  return 0;
}

/*! bam2aln
 Draw HTML alignment from BAM
 */
int do_bam2cov(int argc, char* argv[]) {
  // setup and parse configuration options:
	AnyOption options("Usage: bam2aln --bam=<reference.bam> --fasta=<reference.fasta> --region=<accession:start-end> --output=<output.html> [--max-reads=1000]");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  // required options
  ("bam,b", "bam file containing sequences to be aligned", "data/reference.bam")
	("fasta,f", "FASTA file of reference sequence", "data/reference.fasta")
  ("output,o", "base name of output file. region appended if multiple output files")
  // which regions to create files for
  ("region,r", "region to print (accession:start-end)", "")
  ("tile-size", "size of tiles")
  ("tile-overlap", "overlap between tiles (1/2 on each side)")
  // options controlling what files are output
  ("plot,p", "create graphical plot of coverage", TAKES_NO_ARGUMENT)
  ("plot-format", "plot format: PNG or PDF", "PNG")
  ("table,t", "create text table of coverage", TAKES_NO_ARGUMENT)
//  ("read_start_output,r", "file name for table file binned by read start bases (DEFAULT: OFF)")
//  ("gc_output,g", "create additional table file binned by GC content of reads (DEFAULT: OFF)")
  // options controlling information that is output
  ("total-only,1", "only plot/tabulate total coverage, not per strand", TAKES_NO_ARGUMENT)
  ("resolution", "rough number of positions to output coverage information for in interval (0=ALL)", 600)
  .processCommandArgs(argc, argv);
  
  // make sure that the required config options are good:
	if(options.count("help")
     || !file_exists(options["fasta"].c_str())
     || !file_exists(options["bam"].c_str()) )
  {
		options.printUsage();
		return -1;
	}
  
  ASSERT(options.count("plot") || options.count("table"), "Must specify either --plot or --table.");

  vector<string> region_list;
  if (options.count("region"))
    region_list= from_string<vector<string> >(options["region"]);
  
  // also take regions off the command line
  for (int32_t i = 0; i < options.getArgc(); i++)
  {
    string region = options.getArgv(i);
    region_list.push_back(region);
  }
  
  bool tiling_mode = options.count("tile-size") && options.count("tile-overlap");
  ASSERT(tiling_mode || (!options.count("tile-size") && !options.count("tile-overlap")),
          "--tile-size and --tile-overlap args must both be provided to activate tile mode");
  
  if (tiling_mode && (region_list.size() > 0))
  {
    WARN("Tiling mode activated. Ignoring " + to_string(region_list.size()) + "regions that were specified.");
    region_list.clear();
  }
  
  ASSERT(tiling_mode || (region_list.size() > 0), "No regions specified.");
  
  // create empty settings object to have R script name
  Settings settings;
    
  // generate coverage table/plot!
  coverage_output co(
                      options["bam"],
                      options["fasta"],
                      settings.coverage_plot_r_script_file_name
                      );
  
  // Set options
  co.total_only(options.count("total-only"));
  co.output_format( options["plot-format"] );
  
  // create regions that tile the genome
  if (tiling_mode)
  {
    int32_t tile_size = from_string<int32_t>(options["tile-size"]);
    int32_t tile_overlap = from_string<int32_t>(options["tile-overlap"]);
    tile_overlap = floor( static_cast<double>(tile_overlap) / 2.0);

    for(uint32_t target_id=0; target_id < co.num_targets(); target_id++)
    {
      const char* target_name = co.target_name(target_id);
      int32_t target_length = co.target_length(target_id);
      
      int32_t start = 1;
      while (start < target_length)
      {
        int32_t end = start + tile_size - 1;
        
        int32_t offset_start = start - tile_overlap;
        if (offset_start < 1) offset_start = 1;
        
        int32_t offset_end = end + tile_overlap;
        if (offset_end > target_length) offset_end = target_length;
        
        string region = target_name;
        region += ":" + to_string(offset_start) + "-" + to_string(offset_end);
        
        region_list.push_back(region);
        
        start += tile_size;
      }
    }
  }
  
  for(vector<string>::iterator it = region_list.begin(); it!= region_list.end(); it++)
  {
// these are experimental... additional table files
//    if (options.count("read_start_output"))
//      co.read_begin_output_file_name(options["read_start_output"]);
//    if (options.count("gc_output"))
//      co.gc_output_file_name(options["gc_output"]);
    
    string file_name = options["output"];
    if ((region_list.size() > 0) || (file_name == ""))  file_name += *it;
    cout << "Coverage for region: " << *it << endl;
    
    if (options.count("table"))
      co.table(*it, file_name + ".tab", from_string<uint32_t>(options["resolution"]));
    
    if (options.count("plot"))
      co.plot(*it, file_name + "." + to_lower(options["plot-format"]) , from_string<uint32_t>(options["resolution"]));
  }
  
  return 0;
}

int do_convert_fastq(int argc, char* argv[])
{
  
	// setup and parse configuration options:
	AnyOption options("Usage: breseq ANALYZE_FASTQ --input input.fastq --convert converted.fastq");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  ("input,i", "input FASTQ file")
  ("output,o", "output FASTQ file")
  ("in-format,1", "format to convert from")
  ("out-format,2", "format to convert to")
  
  //("output,o", "out to files") // outputs to STDOUT for now
	.processCommandArgs(argc, argv);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("input")
     || !options.count("output")
     || !options.count("in-format")
     || !options.count("out-format")
		 ) {
		options.printUsage();
		return -1;
	}                       
  
	try {
    
    convert_fastq(options["input"], options["output"], options["in-format"], options["out-format"]);
    
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
  
  
  for (vector<string>::const_iterator it = in.begin(); it < in.end(); it++) 
  {
    refseqs.ReadGenBank(*it);
  }
  
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
	AnyOption options("Usage: breseq CALCULATE_TRIMS --bam <sequences.bam> --fasta <reference.fasta> --error_dir <path> --cGenomeDiff <path> --output <path> --readfile <filename> --coverage_dir <dirname> [--minimum-quality-score 3]");
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
    ref_seq_info.LoadFile(settings.reference_gff3_file_name);
        
    MutationPredictor mp(ref_seq_info);
    
    cGenomeDiff gd( settings.evidence_cGenomeDiff_file_name );
    
    mp.predict(
               settings,
               gd,
               from_string<uint32_t>(options["maximum-read-length"])
               );
    
    gd.write(settings.final_cGenomeDiff_file_name);
    
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
	AnyOption options("Usage: breseq ERROR_COUNT --bam reference.bam --fasta reference.fasta --output test --readfile reads.fastq --covariates ref_base,obs_base,quality=40 [--coverage] [--errors] [--minimum-quality-score 3]");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		("bam,b", "bam file containing sequences to be aligned", "data/reference.bam")
		("fasta,f", "FASTA file of reference sequence", "data/reference.fasta")
		("output,o", "output directory", "./")
		("readfile,r", "name of readfile (no extension). may occur multiple times")
		("coverage", "generate unique coverage distribution output", TAKES_NO_ARGUMENT)
		("errors", "generate unique error count output", TAKES_NO_ARGUMENT)
    ("covariates", "covariates for error model. a comma separated list (no spaces) of these choices: ref_base, obs_base, prev_base, quality, read_set, ref_pos, read_pos, base_repeat. For quality, read_pos, and base_repeat you must specify the maximum value possible, e.g. quality=40")
    ("minimum-quality-score", "ignore base quality scores lower than this", 0)
	.processCommandArgs(argc, argv);
  
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("readfile")
		 || (!options.count("coverage") && !options.count("errors")) ) {
		options.printUsage();
		return -1;
	}
  
  if(options.count("errors") && !options.count("covariates") ) {
    WARN("Must provide --covariates when --errors specified.");
		options.printUsage();
		return -1;
	}
	
	// attempt to calculate error calibrations:
	try {
    Summary summary;
		breseq::error_count(
                        summary,
                        options["bam"],
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
	AnyOption options("Usage: breseq IDENTIFY_MUTATIONS --bam <sequences.bam> --fasta <reference.fasta> --error_dir <path> --cGenomeDiff <path> --output <path> --readfile <filename> --coverage_dir <dirname>");
	options
		("help,h", "produce this help message", TAKES_NO_ARGUMENT)
		("bam,b", "bam file containing sequences to be aligned")
		("fasta,f", "FASTA file of reference sequence")
		("readfile,r", "names of readfile (no extension)")
		("error_dir,e", "Directory containing error rates files")
		("error_table", "Error rates files", "")
		("cGenomeDiff,g", "Genome diff file")
		("output,o", "output directory")
		("coverage_dir", "directory for coverage files", "")
		("mutation_cutoff,c", "mutation cutoff (log10 e-value)", 2.0L)
    ("deletion_propagation_cutoff,u", "number after which to cutoff deletions")
    ("deletion_seed_cutoff", "value below which to cutoff deletions", 0)
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
		 || !options.count("cGenomeDiff")
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
                         options["cGenomeDiff"],
                         options["output"],
                         split(options["readfile"], "\n"),
                         options["coverage_dir"],
                         from_string<vector<double> >(options["deletion_propagation_cutoff"]),
                         from_string<vector<double> >(options["deletion_seed_cutoff"]),
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
int do_tabulate_contingency_loci(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq TABULATE_CL --bam <sequences.bam> --fasta <reference.fasta> --output <path> --loci <loci.txt> --strict");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  ("bam,b", "bam file containing sequences to be aligned")
  ("fasta,f", "FASTA file of reference sequence")
  ("output,o", "output file")
  ("loci,l", "Contingency loci coordinates" )
  ("strict,s", "exclude non-perfect matches in surrounding 5 bases", TAKES_NO_ARGUMENT)
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
    
    string loci_file = "";
    if (options.count("loci")) loci_file = options["loci"];
    
    analyze_contingency_loci(
                             options["bam"],
                             options["fasta"],
                             options["output"],
                             loci_file,
                             options.count("strict")
                             );
	} catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}


/*! Graph Contingency Loci
 
 Combines output of tabulate contingency loci across difference genomes, 
 and uses R to produce bar plots from this combined file
 
 */
int do_graph_contingency_loci(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq GRAPH_CL --bam <sequences.bam> --fasta <reference.fasta> --output <path> --loci <loci.txt> --strict");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  ("output,o", "output file")
  ("loci,l", "Contingency loci coordinates" )
  ("strict,s", "exclude non-perfect matches in surrounding 5 bases", TAKES_NO_ARGUMENT)
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
    
    vector<string> v = options.getRemainingArgs();
    
    string loci_file = "";
    if (options.count("loci")) loci_file = options["loci"];
    
    analyze_contingency_loci(
                             options["bam"],
                             options["fasta"],
                             options["output"],
                             loci_file,
                             options.count("strict")
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

		("min-indel-split-length", "split indels this long in matches", 3)
		("max-read-mismatches", "ignore reads with more than this number of mismatches", int32_t(-1))
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
     
  settings.read_files.Init(from_string<vector<string> >(options["read-file"]));
 
	settings.candidate_junction_fasta_file_name = options["candidate-junction-path"];
  settings.candidate_junction_fasta_file_name += "/candidate_junctions.fasta";
	settings.candidate_junction_faidx_file_name = settings.candidate_junction_fasta_file_name + ".fai";
    
	settings.candidate_junction_sam_file_name = options["candidate-junction-path"];
  settings.candidate_junction_sam_file_name += "/#.candidate_junction.sam";
        
	settings.preprocess_junction_split_sam_file_name = options["candidate-junction-path"];
  settings.preprocess_junction_split_sam_file_name += "/#.split.sam";
    
	settings.preprocess_junction_best_sam_file_name = options["candidate-junction-path"];
  settings.preprocess_junction_best_sam_file_name += "/best.sam";
    
	settings.reference_fasta_file_name = options["data-path"] + "/reference.fasta";
  settings.reference_faidx_file_name = settings.reference_fasta_file_name + ".fai";
    
	settings.reference_sam_file_name = options["reference-alignment-path"];
  settings.reference_sam_file_name += "/#.reference.sam";
    
  settings.max_read_mismatches = from_string<int32_t>(options["max-read-mismatches"]);
	settings.candidate_junction_read_limit = from_string<int32_t>(options["candidate-junction-read-limit"]);

  // settings controlling which alignments are used
  settings.require_complete_match = options.count("require-complete-match");
  settings.require_match_length = from_string<int32_t>(options["require-match-length"]);
  settings.require_match_fraction = from_string<int32_t>(options["require-match-fraction"]);
    
  settings.preprocess_junction_min_indel_split_length = from_string<int32_t>(options["min-indel-split-length"]);
 
	cReferenceSequences ref_seqs;
	ref_seqs.ReadFASTA(settings.reference_fasta_file_name);

	PreprocessAlignments::preprocess_alignments(settings, summary, ref_seqs);

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
    ("require-match-length,4",
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
    settings.read_files.Init(from_string<vector<string> >(options["read-file"]));
    settings.preprocess_junction_split_sam_file_name = options["candidate-junction-path"] + "/#.split.sam";
    settings.reference_fasta_file_name = options["data-path"] + "/reference.fasta";     
    settings.candidate_junction_fasta_file_name = options["candidate-junction-path"] + "/candidate_junction.fasta";

    // Other settings
    settings.candidate_junction_read_limit = from_string<int32_t>(options["candidate-junction-read-limit"]);
    settings.require_match_length  = from_string<int32_t>(options["require-match-length"]);
    settings.required_one_unique_length_per_side = from_string<int32_t>(options["required-one-unique-length-per-side"]);
    settings.required_both_unique_length_per_side = from_string<int32_t>(options["required-both-unique-length-per-side"]);
    
    settings.maximum_candidate_junctions = from_string<int32_t>(options["maximum-candidate-junctions"]);
    settings.minimum_candidate_junctions = from_string<int32_t>(options["minimum-candidate-junctions"]);
    settings.maximum_candidate_junction_length_factor = from_string<double>(options["maximum-candidate-junction-length-factor"]);

    // We should inherit the summary object from earlier steps
    Summary summary;
    summary.sequence_conversion.total_reference_sequence_length = from_string<uint32_t>(options["reference-sequence-length"]);
    summary.sequence_conversion.max_read_length = from_string<int32_t>(options["maximum-read-length"]);

    cReferenceSequences ref_seq_info;
    ref_seq_info.ReadFASTA(options["data-path"] + "/reference.fasta");
        
    CandidateJunctions::identify_candidate_junctions(settings, summary, ref_seq_info);
    
  } catch(...) {
		// failed; 
		return -1;
	}
  
  return 0;
}



int do_convert_gvf( int argc, char* argv[]){
  AnyOption options("Usage: breseq GD2GVF --gd <genomediff.gd> --output <gvf.gvf>"); 

  options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("input,i","gd file to convert") 
    ("output,o","name of output file")
    ("snv-only", "only output SNV entries", TAKES_NO_ARGUMENT)
    ;
  options.processCommandArgs( argc,argv);
  
  if( !options.count("input") || !options.count("output") ){
      options.printUsage(); return -1;
  }
  
  try{
      GDtoGVF( options["input"], options["output"], options.count("snv-only") );
  } 
  catch(...){ 
      return -1; // failed 
  }
  
  return 0;
}

int do_convert_gd( int argc, char* argv[])
{
    AnyOption options("Usage: VCF2GD --vcf <vcf.vcf> --output <gd.gd>");
    options("vcf,i","gd file to convert");
    options("output,o","name of output file");
    options.processCommandArgs( argc,argv);
    
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

int do_mutate(int argc, char *argv[])
{
  AnyOption options("Usage: breseq APPLY -g <file.gd> -r <reference>");
  options("genomediff,g", "genome diff file");
  options("reference,r",".gbk/.gff3/.fasta/.bull reference sequence file");
  options("fasta,f","output FASTA file");
  options("gff3,3","output GFF3 file");
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Input a single GenomeDiff, and as many reference files");
  options.addUsage("as you like.  Using the GenomeDiff, this will apply all");
  options.addUsage("the mutations to the reference sequences, output is to");
  options.addUsage("a single file that includes all the references.");
  
  if (!options.count("genomediff") ||
      !options.count("reference")) {
    options.addUsage("");
    options.addUsage("You must supply BOTH the --genomediff and --reference");
    options.addUsage("options for input.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("gff3") && !options.count("fasta")) {
    options.addUsage("");
    options.addUsage("You must supply at least one of the --fasta or --gff3 options for output.");
    options.printUsage();
    return -1;
  }  
  
  cGenomeDiff gd(options["genomediff"]);
  cReferenceSequences ref_seq_info;  
  ref_seq_info.LoadFiles(from_string<vector<string> >(options["reference"]));
  
  //Check to see if every item in the loaded .gd is
  // applicable to the reference file.
  ASSERT(gd.is_valid(ref_seq_info, options.count("verbose")), "Reference file and GenomeDiff file don't match.");
  
  cReferenceSequences new_ref_seq_info = gd.apply_to_sequences(ref_seq_info, options.count("verbose"));
  
  if (options.count("fasta"))
    new_ref_seq_info.WriteFASTA(options["fasta"], options.count("verbose"));
  if (options.count("gff3"))
    new_ref_seq_info.WriteGFF(options["gff3"], options.count("verbose"));
  
  return 0;
}

int do_subtract(int argc, char *argv[])
{
  AnyOption options("Usage: breseq SUBTRACT -1 <file.gd> -2 <file.gd> -o <output.gd>");
  options("input1,1","input GD file 1");
  options("input2,2","input GD file 2");
  options("output,o","output GD file");
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Input two GenomeDiff files and compare mutations.");
  options.addUsage("Mutations that appear in BOTH will be subtracted from");
  options.addUsage("--input1 and results will output to a new file specified.");
  
  if (!options.count("input1") ||
      !options.count("input2")) {
    options.addUsage("");
    options.addUsage("You must supply BOTH the --input1 and --input2 options");
    options.addUsage("for input.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("You must supply the --output option for output.");
    options.printUsage();
    return -1;
  }
  
  cGenomeDiff gd1(options["input1"]);
  cGenomeDiff gd2(options["input2"]);
  
  gd1.subtract(gd2, options.count("verbose"));
  
  gd1.write(options["output"]);
  
  return 0;
}

int do_merge(int argc, char *argv[])
{
  AnyOption options("Usage: breseq MERGE -g <file1.gd file2.gd file3.gd ...> -o <output.gd>");
  options("genomediff,g","input GD files");
  options("output,o","output GD file");
  options("unique,u","Unique Entries Only (Flag)", TAKES_NO_ARGUMENT);  
  options("id,i","Reorder IDs (Flag)", TAKES_NO_ARGUMENT);
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Input as many GenomeDiff files as you want, and have them");
  options.addUsage("merged together into a single GenomeDiff file specified.");
  options.addUsage("Unique IDs will remain unique across files, any IDs that");
  options.addUsage("aren't unique will get new ones.  Mutations that use those");
  options.addUsage("IDs will be properly updated to acknowledge the change.");
  
  if (!options.count("genomediff")) {
    options.addUsage("");
    options.addUsage("You must supply the --genomediff option for input.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("You must supply the --output option for output.");
    options.printUsage();
    return -1;
  }
  
  vector<string> gd_list = from_string<vector<string> >(options["genomediff"]);  
  cGenomeDiff gd1(gd_list[0]);
  
  //Load all the GD files that were input with -g.
  for(uint32_t i = 1; i < gd_list.size(); i++)
  {
    cGenomeDiff gd2(gd_list[i]);
    gd1.merge(gd2, options.count("unique"), options.count("id"), options.count("verbose"));
  }
  
  //Treat EVERY extra argument passed as another GD file.
  for(int32_t i = 0; i < options.getArgc() ; i++)
  {
    cGenomeDiff gd2(options.getArgv(i));
    gd1.merge(gd2, options.count("unique"), options.count("id"), options.count("verbose"));
  }
  
  gd1.write(options["output"]);
  
  return 0;
}

int do_not_evidence(int argc, char *argv[])
{
  AnyOption options("Usage: breseq NOT_EVIDENCE -g <genomediff.gd> -o <output.gd>");
  options("genomediff,g","input GD files");
  options("output,o","output GD file");
  options("id,i","Reorder IDs (Flag)", TAKES_NO_ARGUMENT);
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Takes a GenomeDiff file and removes all of the entries that");
  options.addUsage("are NOT used as evidence by a mutation.  Outputs to a new");
  options.addUsage("GenomeDiff file if specified.  If no output is specified,");
  options.addUsage("verbose will still inform what evidence isn't being used.");
  
  if (!options.count("genomediff")) {
    options.addUsage("");
    options.addUsage("You must supply the --genomediff option for input.");
    options.printUsage();
    return -1;
  }
  
  cGenomeDiff gd1(options["genomediff"]);
  
  gd1.filter_not_used_as_evidence(options.count("verbose"));
  
  if(options.count("output"))
  {
    if(options.count("id"))
    {
      cGenomeDiff gd2;
      gd2.merge(gd1, true, true);
      gd2.write(options["output"]);
    }
    else  {
      gd1.write(options["output"]);  }
  }
  
  return 0;
}

int do_compare(int argc, char *argv[])
{
  AnyOption options("Usage: breseq COMPARE -c <control.gd> -t <test.gd> -o <output.gd>");
  options("control,c", "control genome diff file, mutations within are assumed to be accurate");
  options("test,t",    "test genome diff file, mutations not been checked for accuracy");
  options("output,o",  "output compared genome diff file name");
  options.processCommandArgs(argc, argv);

  options.addUsage("Usage: breseq COMPARE -c <control.gd> -t <test.gd> -o <output.gd>");

  if (!options.count("control") ||
      !options.count("test")    ||
      !options.count("output")) {
    options.printUsage();
    return -1;
  }

  cGenomeDiff control_gd(options["control"]);
  cGenomeDiff test_gd(options["test"]);

  cGenomeDiff compare_gd = cGenomeDiff::compare_cGenomeDiff_files(control_gd, test_gd);

  compare_gd.write(options["output"]);


  return 0;
}

int do_intersection(int argc, char *argv[])
{
  AnyOption options("Usage: breseq INTERSECTION -o <output.gd> <file1.gd file2.gd file3.gd ...>");
  options("output,o", "output intersection mutations to this file");
  options.processCommandArgs(argc, argv);

  typedef vector<string> string_vector_t;
  string_vector_t gd_file_names;
  //Get GD Files
  for (int32_t i = 0; i < options.getArgc() ; i++) {
    const string &file_name = options.getArgv(i);
    gd_file_names.push_back(file_name);
  }

  if (!options.count("output") ||
      gd_file_names.size() < 2) {
    printf("\n");
    printf("ERROR!\n");
    printf("Ouput:%s\n", options["output"].c_str());
    printf("GenomeDiff file count: %i\n", static_cast<int>(gd_file_names.size()));
    printf("\n");
    options.printUsage();
    return -1;
  }

  cGenomeDiff first_gd(*gd_file_names.begin());
  const diff_entry_list_t &first_mutations = first_gd.mutation_list();

  //Strip counted_ptr
  typedef set<cDiffEntry> cDiffEntry_set_t;
  cDiffEntry_set_t first_mutations_set;
  for (diff_entry_list_t::const_iterator it = first_mutations.begin();
       it != first_mutations.end(); it++) {
    first_mutations_set.insert(**it);
  }

  /*! Step: Compare to other genome diffs and reassign first_mutations_set when applicable.
    Stop when there is no intersection or we run out of file_names
  */
  string_vector_t::const_iterator it_file_name = gd_file_names.begin();
  advance(it_file_name,1);
  while(it_file_name != gd_file_names.end()) {

    cGenomeDiff second_gd(*it_file_name);
    const diff_entry_list_t &second_mutations = second_gd.mutation_list();

    //Strip counted_ptr
    cDiffEntry_set_t second_mutations_set;
    for (diff_entry_list_t::const_iterator it = second_mutations.begin();
         it != second_mutations.end(); it++) {
      second_mutations_set.insert(**it);
    }

    //Find intersection
    cDiffEntry_set_t intersecting_mutations_set;
    set_intersection(first_mutations_set.begin(), first_mutations_set.end(),
                     second_mutations_set.begin(), second_mutations_set.end(),
                     inserter(intersecting_mutations_set, intersecting_mutations_set.begin()));

    if (intersecting_mutations_set.empty()) {
      printf("No intersecting mutations were found across the current file: %s\n",
             it_file_name->c_str());
      return -1;
    }


    first_mutations_set = intersecting_mutations_set;

    it_file_name++;
  }

  printf("Found %i intersecting mutations across files: %s.\n",
         static_cast<int>(first_mutations_set.size()), join(gd_file_names, ", ").c_str());

  cGenomeDiff intersecting_gd;

  for (cDiffEntry_set_t::iterator it = first_mutations_set.begin();
       it != first_mutations_set.end(); it++) {
    //cDiffEntry de = *it;
    intersecting_gd.add(*it);
  }

  intersecting_gd.write(options["output"]);

	return 0;
}

int do_annotate(int argc, char* argv[])
{
  AnyOption options("Usage: breseq ANNOTATE -r reference.gbk  -i input.gd -o annotated.gd");
  
  options
  ("help,h", "produce advanced help message", TAKES_NO_ARGUMENT)
  // convert to basing everything off the main output path, so we don't have to set so many options
  ("input,i", "path to input genome diff (REQUIRED)")
  ("output,o", "path to output genome diff with added mutation data (REQUIRED)")
  ("reference,r", "reference sequence in GenBank flatfile format (REQUIRED)")
  ("ignore-pseudogenes", "treats pseudogenes as normal genes for callings AA changes", TAKES_NO_ARGUMENT)
  ;
  options.processCommandArgs(argc, argv);
  
  if ( !options.count("input")
    || !options.count("output")
    || !options.count("reference")
      ) {
    options.printUsage();
    return -1;
  }
  
  cGenomeDiff gd(options["input"]);
  
  vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(reference_file_names);
  ref_seq_info.annotate_mutations(gd, false, options.count("ignore-pseudogenes"));
  
  gd.write(options["output"]);
  
  return 0;
}

int do_simulate_read(int argc, char *argv[])
{
  AnyOption options("Usage: breseq SIMULATE-READ -g <genome diff> -r <reference file> -c <average coverage> -o <output file>");

  options
  ("cGenomeDiff,g", "Genome diff file.")
  ("reference,r", "Reference file for input.")
  ("coverage,c", "Average coverage value to simulate.", static_cast<uint32_t>(10))
  ("length,l", "Read length to simulate.", static_cast<uint32_t>(36))
  ("output,o", "Output fastq file name.")  
  ("gff3,3", "Output Applied GFF3 File. (Flag)", TAKES_NO_ARGUMENT)
  ("verbose,v", "Verbose Mode (Flag)", TAKES_NO_ARGUMENT)
  ;
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Takes a GenomeDiff file and applies the mutations to the reference.");
  options.addUsage("Then using the applied reference, it simulates reads based on it.");
  options.addUsage("Output is a .fastq that if run normally against the reference");
  options.addUsage("should produce the same GenomeDiff file you entered.");
  
  if (!options.count("cGenomeDiff")) {
    options.addUsage("");
    options.addUsage("You must supply the --cGenomeDiff option for input.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("reference")) {
    options.addUsage("");
    options.addUsage("You must supply the --reference option for input.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("You must supply the --output option for output.");
    options.printUsage();
    return -1;
  }

  //! Step: Load reference sequence file.
  cReferenceSequences ref_seq_info;
  const string &ref_file_name = options["reference"];
  ref_seq_info.LoadFile(ref_file_name);


  //! Step: Apply genome diff mutations to reference sequence.
  const string &gd_file_name = options["cGenomeDiff"];
  bool verbose = options.count("verbose");
  cGenomeDiff gd(gd_file_name);

  cReferenceSequences new_ref_seq_info = gd.apply_to_sequences(ref_seq_info, verbose);  
  
  //! Write applied GFF3 file if requested.
  if(options.count("gff3"))new_ref_seq_info.WriteGFF(options["output"] + ".gff3", options.count("verbose"));

  const cFastaSequence &fasta_sequence = new_ref_seq_info.front().m_fasta_sequence;

  //! Step: Create simulated read sequence.
  //Parameters to create simulated read.
  const cSequence &sequence = fasta_sequence.m_sequence;
  const uint32_t average_coverage = from_string<uint32_t>(options["coverage"]);
  const uint32_t read_length = from_string<uint32_t>(options["length"]);

  const cFastqSequenceVector &fsv =
      cFastqSequenceVector::simulateFromSequence(sequence, average_coverage, read_length, verbose);

  //! Step: Write simulated reads to file.
  if(verbose){cout << "Writing FASTQ\n" << "\t" << options["output"] << endl;}
  cFastqFile ff(options["output"], ios_base::out);

  const size_t &fsv_size = fsv.size();
  for (size_t i = 0; i < fsv_size; i++) {
    ff.write_sequence(fsv[i]);
    if(verbose && !(i % 10000) && i){cout << "\tREAD: " << i << endl;}
  }
  if(verbose){cout << "\t**FASTQ Complete**" << endl;}
  return 0;
}


int breseq_default_action(int argc, char* argv[])
{  
  
	///
	/// Get options from the command line
	///
  Summary summary;
	Settings settings(argc, argv);
	settings.check_installed();
  
	//
	// 01_sequence_conversion 
  // * Convert the input reference into FASTA for alignment and GFF3 for reloading features
  // * Rename reads in the input FASTQ and change quality scores to Sanger
	//
	create_path(settings.sequence_conversion_path);
  create_path(settings.data_path);

  cReferenceSequences ref_seq_info;

	if (settings.do_step(settings.sequence_conversion_done_file_name, "Read and reference sequence file input"))
	{
		Summary::SequenceConversion s;
    cReferenceSequences conv_ref_seq_info;
    
    // Load all of the reference sequences and convert to FASTA and GFF3
    conv_ref_seq_info.LoadFiles(settings.reference_file_names);
    conv_ref_seq_info.WriteFASTA(settings.reference_fasta_file_name);
    conv_ref_seq_info.WriteGFF(settings.reference_gff3_file_name);

    //Check the FASTQ format and collect some information about the input read files at the same time
		cerr << "  Analyzing fastq read files..." << endl;
		uint32_t overall_max_read_length = UNDEFINED_UINT32;
		uint32_t overall_max_qual = 0;

		s.num_reads = 0;
    s.num_bases = 0;
		for (uint32_t i = 0; i < settings.read_files.size(); i++)
		{
			string base_name = settings.read_files[i].m_base_name;
			cerr << "    READ FILE::" << base_name << endl;
			string fastq_file_name = settings.base_name_to_read_file_name(base_name);
			string convert_file_name =  settings.file_name(settings.converted_fastq_file_name, "#", base_name);

			// Parse output
			AnalyzeFastq s_rf = normalize_fastq(fastq_file_name, convert_file_name, i+1);
      
			// Save the converted file name -- have to save it in summary because only that
			// is reloaded if we skip this step.
			s.converted_fastq_name[base_name] = s_rf.converted_fastq_name;

			// Record statistics
			if ((overall_max_read_length == UNDEFINED_UINT32) || (s_rf.max_read_length > overall_max_read_length))
				overall_max_read_length = s_rf.max_read_length;
			if ((overall_max_qual == UNDEFINED_UINT32) || (s_rf.max_quality_score > overall_max_qual))
				overall_max_qual = s_rf.max_quality_score;
			s.num_reads += s_rf.num_reads;
			s.num_bases += s_rf.num_bases;

			s.reads[base_name] = s_rf;
		}
		s.avg_read_length = s.num_bases / s.num_reads;
		s.max_read_length = overall_max_read_length;
		s.max_qual = overall_max_qual;
		summary.sequence_conversion = s;

		// create SAM faidx
		string samtools = settings.ctool("samtools");
		string command = samtools + " faidx " + settings.reference_fasta_file_name;
		SYSTEM(command.c_str());
    
		// calculate trim files
		calculate_trims(settings.reference_fasta_file_name, settings.sequence_conversion_path);

		// store summary information
		summary.sequence_conversion.store(settings.sequence_conversion_summary_file_name);
		settings.done_step(settings.sequence_conversion_done_file_name);
	}

	summary.sequence_conversion.retrieve(settings.sequence_conversion_summary_file_name);
	ASSERT(summary.sequence_conversion.max_read_length != UNDEFINED_UINT32, "Can't retrieve max read length from file: " + settings.sequence_conversion_summary_file_name);

  //(re)load the reference sequences from our converted files
  ref_seq_info.ReadGFF(settings.reference_gff3_file_name);
  
  // Calculate the total reference sequence length
  summary.sequence_conversion.total_reference_sequence_length = ref_seq_info.total_length();
  
  // @JEB -- This is a bit of an ugly wart from when converting the input file was optional.
	// reload certain information into $settings from $summary  
	for (map<string, AnalyzeFastq>::iterator it = summary.sequence_conversion.reads.begin(); it != summary.sequence_conversion.reads.end(); it++)
	{
		string read_file = it->first;
		if (it->second.converted_fastq_name.size() > 0)
			settings.read_files.read_file_to_converted_fastq_file_name_map[read_file] = it->second.converted_fastq_name;
	}

	//
  // 02_reference_alignment
	// * Match all reads against the reference genome
	//

	if (settings.do_step(settings.reference_alignment_done_file_name, "Read alignment to reference genome"))
	{
		create_path(settings.reference_alignment_path);

		/// create ssaha2 hash
		string reference_hash_file_name = settings.reference_hash_file_name;
		string reference_fasta_file_name = settings.reference_fasta_file_name;

		if (!settings.smalt)
		{
			string command = "ssaha2Build -rtype solexa -skip " + to_string(settings.ssaha2_skip_length) + " -save " + reference_hash_file_name + " " + reference_fasta_file_name;
			SYSTEM(command.c_str());
		}
		else
		{
			string smalt = settings.ctool("smalt");
			string command = smalt + " index -k " + to_string(settings.ssaha2_seed_length) + " -s " + to_string(settings.ssaha2_skip_length) + " " + reference_hash_file_name + " " + reference_fasta_file_name;
      SYSTEM(command.c_str());
		}
		/// ssaha2 align reads to reference sequences
		for (uint32_t i = 0; i < settings.read_files.size(); i++)
		{
			cReadFile read_file = settings.read_files[i];

			//reads are paired -- Needs to be re-implemented with major changes elsewhere. @JEB
			//if (is_defined(read_struct.min_pair_dist) && is_defined(read_struct.max_pair_dist))
			//if (is_defined(read_file.m_paired_end_group))
			//{
				// JEB this is not working currently
			//	breseq_assert(false);
				/*
				die "Paired end mapping is broken.";
				die if (scalar @{$read_struct->{read_fastq_list}} != 2);

				my $fastq_1 = $read_struct->{read_fastq_list}->[0];
				my $fastq_2 = $read_struct->{read_fastq_list}->[1];
				my $min = $read_struct->{min_pair_dist};
				my $max = $read_struct->{max_pair_dist};

				my $reference_sam_file_name = $settings->file_name("reference_sam_file_name", {"//"=>$read_struct->{base_name}});
				Breseq::Shared::system("ssaha2 -disk 2 -save $reference_hash_file_name -kmer 13 -skip 1 -seeds 1 -score 12 -cmatch 9 -ckmer 1 -output sam_soft -outfile $reference_sam_file_name -multi 1 -mthresh 9 -pair $min,$max $fastq_1 $fastq_2");
				*/
			//}
			//else //reads are not paired
			{
				string base_read_file_name = read_file.base_name();
				string read_fastq_file = settings.base_name_to_read_file_name(base_read_file_name);
				string reference_sam_file_name = settings.file_name(settings.reference_sam_file_name, "#", base_read_file_name);

				if (!settings.smalt)
				{
					string command = "ssaha2 -disk 2 -save " + reference_hash_file_name + " -kmer " + to_string(settings.ssaha2_seed_length) + " -skip " + to_string(settings.ssaha2_skip_length) + " -seeds 1 -score 12 -cmatch 9 -ckmer 1 -output sam_soft -outfile " + reference_sam_file_name + " " + read_fastq_file;
					SYSTEM(command.c_str());
				}
				else
				{
					string smalt = settings.ctool("smalt");
					string command = smalt + " map -n 2 -d " + to_string(settings.max_smalt_diff) + " -f samsoft -o " + reference_sam_file_name + " " + reference_hash_file_name + " " + read_fastq_file; // -m 12
					SYSTEM(command.c_str());
				}
			}
		}

		/// Delete the hash files immediately
		if (!settings.keep_all_intermediates)
		{
			remove( (reference_hash_file_name + ".base").c_str() );
			remove( (reference_hash_file_name + ".body").c_str() );
			remove( (reference_hash_file_name + ".head").c_str() );
			remove( (reference_hash_file_name + ".name").c_str() );
			remove( (reference_hash_file_name + ".size").c_str() );
		}

		settings.done_step(settings.reference_alignment_done_file_name);
	}

  //
	// 03_candidate_junctions
	// * Identify candidate junctions from split read alignments
	//
	if (settings.junction_prediction)
	{
		create_path(settings.candidate_junction_path);

    string preprocess_junction_done_file_name = settings.preprocess_junction_done_file_name;

    if (settings.do_step(settings.preprocess_junction_done_file_name, "Preprocessing alignments for candidate junction identification"))
    {
      PreprocessAlignments::preprocess_alignments(settings, summary, ref_seq_info);
      settings.done_step(settings.preprocess_junction_done_file_name);
    }

    
    if (settings.do_step(settings.coverage_junction_done_file_name, "Preliminary analysis of coverage distribution"))
    {
      string reference_faidx_file_name = settings.reference_faidx_file_name;
      string preprocess_junction_best_sam_file_name = settings.preprocess_junction_best_sam_file_name;
      string coverage_junction_best_bam_file_name = settings.coverage_junction_best_bam_file_name;
      string coverage_junction_best_bam_prefix = settings.coverage_junction_best_bam_prefix;
      string coverage_junction_best_bam_unsorted_file_name = settings.coverage_junction_best_bam_unsorted_file_name;

      string samtools = settings.ctool("samtools");

      string command = samtools + " import " + reference_faidx_file_name + " " + preprocess_junction_best_sam_file_name + " " + coverage_junction_best_bam_unsorted_file_name;
      SYSTEM(command.c_str());
      command = samtools + " sort " + coverage_junction_best_bam_unsorted_file_name + " " + coverage_junction_best_bam_prefix;
      SYSTEM(command.c_str());
      if (!settings.keep_all_intermediates)
        remove(coverage_junction_best_bam_unsorted_file_name.c_str());
      command = samtools + " index " + coverage_junction_best_bam_file_name;
      SYSTEM(command.c_str());

      // Count errors
      string reference_fasta_file_name = settings.reference_fasta_file_name;
      string reference_bam_file_name = settings.coverage_junction_best_bam_file_name;

      error_count(
        summary,
        reference_bam_file_name,
        reference_fasta_file_name,
        settings.candidate_junction_path,
        settings.read_file_names,
        true, // coverage
        false, // errors
        settings.base_quality_cutoff,
        "" //covariates
      );

      CoverageDistribution::analyze_unique_coverage_distributions(settings, 
                                                                  summary, 
                                                                  ref_seq_info,
                                                                  settings.coverage_junction_plot_file_name, 
                                                                  settings.coverage_junction_distribution_file_name
                                                                  );

      // Note that storing from unique_coverage and reloading in preprocess_coverage is by design
      summary.unique_coverage.store(settings.coverage_junction_summary_file_name);
      summary.preprocess_error_count.store(settings.coverage_junction_error_count_summary_file_name);
      settings.done_step(settings.coverage_junction_done_file_name);
		}
    summary.preprocess_coverage.retrieve(settings.coverage_junction_summary_file_name);
    summary.preprocess_error_count.retrieve(settings.coverage_junction_error_count_summary_file_name);
    
		string candidate_junction_summary_file_name = settings.candidate_junction_summary_file_name;
		if (settings.do_step(settings.candidate_junction_done_file_name, "Identifying candidate junctions"))
		{
			cerr << "Identifying candidate junctions..." << endl;
      CandidateJunctions::identify_candidate_junctions(settings, summary, ref_seq_info);

			string samtools = settings.ctool("samtools");
			string faidx_command = samtools + " faidx " + settings.candidate_junction_fasta_file_name;
			if (!file_empty(settings.candidate_junction_fasta_file_name.c_str()))
				SYSTEM(faidx_command.c_str());

			summary.candidate_junction.store(candidate_junction_summary_file_name);
			settings.done_step(settings.candidate_junction_done_file_name);
		}
		summary.candidate_junction.retrieve(candidate_junction_summary_file_name);

    
    //
    // 04 candidate_junction_alignment
    // * Align reads to new junction candidates
    //
		if (settings.do_step(settings.candidate_junction_alignment_done_file_name, "Candidate junction alignment"))
		{
			create_path(settings.candidate_junction_alignment_path);

			/// create ssaha2 hash
			string candidate_junction_hash_file_name = settings.candidate_junction_hash_file_name;
			string candidate_junction_fasta_file_name = settings.candidate_junction_fasta_file_name;

			if (!file_empty(candidate_junction_fasta_file_name.c_str()))
			{
				if (!settings.smalt)
				{
					string command = "ssaha2Build -rtype solexa -skip " + to_string(settings.ssaha2_skip_length) + " -save " + candidate_junction_hash_file_name + " " + candidate_junction_fasta_file_name;
					SYSTEM(command.c_str());
				}
				else
				{
					string smalt = settings.ctool("smalt");
					string command = smalt + " index -k " + to_string(settings.ssaha2_seed_length) + " -s " + to_string(settings.ssaha2_skip_length) + " " + candidate_junction_hash_file_name + " "+ candidate_junction_fasta_file_name;
					SYSTEM(command.c_str());
				}
			}

			/// ssaha2 align reads to candidate junction sequences
			for (uint32_t i = 0; i < settings.read_files.size(); i++)
			{
				string base_read_file_name = settings.read_files[i].m_base_name;
				string candidate_junction_sam_file_name = settings.file_name(settings.candidate_junction_sam_file_name, "#", base_read_file_name);

				string read_fastq_file = settings.base_name_to_read_file_name(base_read_file_name);
				string filename = candidate_junction_hash_file_name + ".base";
				if (!settings.smalt && file_exists(filename.c_str()))
				{
					string command = "ssaha2 -disk 2 -save " + candidate_junction_hash_file_name + " -best 1 -rtype solexa -kmer " + to_string(settings.ssaha2_seed_length) + " -skip " + to_string(settings.ssaha2_skip_length) + " -seeds 1 -output sam_soft -outfile " + candidate_junction_sam_file_name + " " + read_fastq_file;
					SYSTEM(command.c_str());
					// Note: Added -best parameter to try to avoid too many matches to redundant junctions!
				}
				else
				{
					filename = candidate_junction_hash_file_name + ".sma";
					if (file_exists(filename.c_str()))
					{
						string smalt = settings.ctool("smalt");
						string command = smalt + " map -c 0.8 -x -n 2 -d 1 -f samsoft -o " + candidate_junction_sam_file_name + " " + candidate_junction_hash_file_name + " " + read_fastq_file;
						SYSTEM(command.c_str());
					}
				}
			}

			/// Delete the hash files immediately
			if (!settings.keep_all_intermediates)
			{
				remove((candidate_junction_hash_file_name + ".base").c_str());
				remove((candidate_junction_hash_file_name + ".body").c_str());
				remove((candidate_junction_hash_file_name + ".head").c_str());
				remove((candidate_junction_hash_file_name + ".name").c_str());
				remove((candidate_junction_hash_file_name + ".size").c_str());
			}

			settings.done_step(settings.candidate_junction_alignment_done_file_name);
		}
  }

	
	//
  // 05 alignment_correction
	// * Resolve matches to new junction candidates
	//
	if (settings.do_step(settings.alignment_correction_done_file_name, "Resolving alignments with candidate junctions"))
	{
		create_path(settings.alignment_correction_path);
    
		// should be one coverage cutoff value for each reference sequence
		//vector<double> coverage_cutoffs;
		//for (uint32_t i = 0; i < ref_seq_info.size(); i++)
    //{
    //  Coverage f = summary.preprocess_coverage[ref_seq_info[i].m_seq_id];
	//		coverage_cutoffs.push_back(summary.preprocess_coverage[ref_seq_info[i].m_seq_id].junction_accept_score_cutoff);
		//}
    //assert(coverage_cutoffs.size() == ref_seq_info.size());

    bool junction_prediction = settings.junction_prediction;
    if (junction_prediction && file_empty(settings.candidate_junction_fasta_file_name.c_str())) junction_prediction = false;
    
		resolve_alignments(
			settings,
			summary,
			ref_seq_info,
      junction_prediction,
			settings.read_files
		);

		summary.alignment_resolution.store(settings.alignment_resolution_summary_file_name);
		settings.done_step(settings.alignment_correction_done_file_name);
	}
  
	if (file_exists(settings.alignment_resolution_summary_file_name.c_str()))
    summary.alignment_resolution.retrieve(settings.alignment_resolution_summary_file_name);
  
	//
  // 05 bam
	// * Create BAM read alignment database files
	//

	if (settings.do_step(settings.bam_done_file_name, "Creating BAM files"))
	{
		create_path(settings.bam_path);

		string reference_faidx_file_name = settings.reference_faidx_file_name;
		string candidate_junction_faidx_file_name = settings.candidate_junction_faidx_file_name;

		string resolved_junction_sam_file_name = settings.resolved_junction_sam_file_name;
		string junction_bam_unsorted_file_name = settings.junction_bam_unsorted_file_name;
		string junction_bam_prefix = settings.junction_bam_prefix;
		string junction_bam_file_name = settings.junction_bam_file_name;

		string samtools = settings.ctool("samtools");
		string command;

    // only run samtools if we are predicting junctions and there were results in the sam file
    // first part of conditional really not necessary @JEB
		if (settings.junction_prediction && !file_empty(resolved_junction_sam_file_name.c_str()))
		{
			command = samtools + " import " + candidate_junction_faidx_file_name + " " + resolved_junction_sam_file_name + " " + junction_bam_unsorted_file_name;
			SYSTEM(command);
			command = samtools + " sort " + junction_bam_unsorted_file_name + " " + junction_bam_prefix;
      SYSTEM(command);
			if (!settings.keep_all_intermediates)
				remove(junction_bam_unsorted_file_name.c_str());
			command = samtools + " index " + junction_bam_file_name;
			SYSTEM(command);
		}

		string resolved_reference_sam_file_name = settings.resolved_reference_sam_file_name;
		string reference_bam_unsorted_file_name = settings.reference_bam_unsorted_file_name;
		string reference_bam_prefix = settings.reference_bam_prefix;
		string reference_bam_file_name = settings.reference_bam_file_name;

		command = samtools + " import " + reference_faidx_file_name + " " + resolved_reference_sam_file_name + " " + reference_bam_unsorted_file_name;
    SYSTEM(command);
		command = samtools + " sort " + reference_bam_unsorted_file_name + " " + reference_bam_prefix;
    SYSTEM(command);
		if (!settings.keep_all_intermediates)
			remove(reference_bam_unsorted_file_name.c_str());
		command = samtools + " index " + reference_bam_file_name;
    SYSTEM(command);

		settings.done_step(settings.bam_done_file_name);
	}

	//
	//#  Graph paired read outliers (experimental)
	//# sub paired_read_distances {}
	//#
	//# {
	//# 	my @rs = settings.read_structures;
	//#
	//# 	my @min_pair_dist;
	//# 	my @max_pair_dist;
	//#
	//# 	my $paired = 0;
	//#
	//# 	my $i=0;
	//# 	foreach my $rfi (@{settings.{read_file_index_to_struct_index}})
	//# 	{
	//# 		$min_pair_dist[$i] = 0;
	//# 		$max_pair_dist[$i] = 0;
	//#
	//# 		if ($rs[$rfi]->{paired})
	//# 		{
	//# 			$paired = 1;
	//# 			$min_pair_dist[$i] = $rs[$rfi]->{min_pair_dist};
	//# 			$max_pair_dist[$i] = $rs[$rfi]->{max_pair_dist};
	//# 		}
	//# 		$i++;
	//# 	}
	//#
	//# 	my $long_pairs_file_name = settings.file_name("long_pairs_file_name");
	//#
	//# 	if ($paired && (!-e $long_pairs_file_name))
	//# 	{
	//#
	//# 		my $reference_sam_file_name = settings.file_name("resolved_reference_sam_file_name");
	//# 		my $reference_tam = Bio::DB::Tam->open($reference_sam_file_name) or die "Could not open $reference_sam_file_name";
	//#
	//# 		my $reference_faidx_file_name = settings.file_name("reference_faidx_file_name");
	//# 		my $reference_header = $reference_tam->header_read2($reference_faidx_file_name) or throw("Error reading reference fasta index file: $reference_faidx_file_name");
	//# 		my $target_names = $reference_header->target_name;
	//#
	//# 		my $save;
	//# 		my $on_alignment = 0;
	//# 		my $last;
	//#
	//# 		while (1)
	//# 		{
	//# 			$a = Bio::DB::Bam::Alignment->new();
	//# 			my $bytes = $reference_tam->read1($reference_header, $a);
	//# 			last if ($bytes <= 0);
	//#
	//#
	//# 			my $start       = $a->start;
	//# 		    my $end         = $a->end;
	//# 		    my $seqid       = $target_names->[$a->tid];
	//#
	//# 			$on_alignment++;
	//# 			print "$on_alignment\n" if ($on_alignment % 10000 == 0);
	//#
	//# 			//#last if ($on_alignment > 100000);
	//#
	//# 			//#print $a->qname . "" << endl;
	//#
	//# 			if (!$a->unmapped)
	//# 			{
	//# 				my $mate_insert_size = abs($a->isize);
	//# 				my $mate_end = $a->mate_end;
	//# 				my $mate_start = $a->mate_start;
	//# 				my $mate_reversed = 2*$a->mreversed + $a->reversed;
	//# 		 		my $mreversed = $a->mreversed;
	//# 		 		my $reversed = $a->reversed;
	//#
	//# 				my $fastq_file_index = $a->aux_get("X2");
	//# 				//#print "$mate_insert_size $min_pair_dist[$fastq_file_index] $max_pair_dist[$fastq_file_index]" << endl;
	//# 				//#if (($mate_insert_size < $min_pair_dist[$fastq_file_index]) || ($mate_insert_size > $max_pair_dist[$fastq_file_index]))
	//# 				if ((($mate_insert_size >= 400) && ($mate_insert_size < $min_pair_dist[$fastq_file_index])) || ($mate_insert_size > $max_pair_dist[$fastq_file_index]))
	//# 				{
	//# 					//#correct pair
	//#
	//# 					if ($last && ($last->{start} == $mate_start))
	//# 					{
	//# 						$save->{int($start/100)}->{int($mate_start/100)}->{$mate_reversed}++;
	//# 						$save->{int($last->{start}/100)}->{int($last->{mate_start}/100)}->{$last->{mate_reversed}}++;
	//# 						undef $last;
	//# 					}
	//# 					else
	//# 					{
	//# 						($last->{mate_reversed}, $last->{start}, $last->{mate_start}) = ($mate_reversed, $start, $mate_start);
	//# 					}
	//#
	//# 					//#$save->{$mate_reversed}->{int($start/100)}->{int($mate_start/100)}++;
	//# 				    //print $a->qname," aligns to $seqid:$start..$end, $mate_start $mate_reversed ($mreversed $reversed) $mate_insert_size" << endl;
	//# 				}
	//#
	//# 			}
	//# 		}
	//#
	//# 		open LP, ">$long_pairs_file_name" or die;
	//#
	//# 		foreach my $key_1 (sort {$a <=> $b} keys %$save)
	//# 		{
	//# 			foreach my $key_2 (sort {$a <=> $b} keys %{$save->{$key_1}})
	//# 			{
	//# 				foreach my $key_reversed (sort {$a <=> $b} keys %{$save->{$key_1}->{$key_2}})
	//# 				{
	//# 					print LP "$key_1\t$key_2\t$key_reversed\t$save->{$key_1}->{$key_2}->{$key_reversed}" << endl;
	//# 				}
	//# 			}
	//# 		}
	//# 		close LP;
	//# 	}
	//#
	//# 	if ($paired)
	//# 	{
	//# 		open LP, "$long_pairs_file_name" or die;
	//# 		while ($_ = <LP>)
	//# 		{
	//# 			chomp $_;
	//# 			my ($start, $end, $key_reversed);
	//# 		}
	//# 	}
	//# }
	//#

	//
	// Tabulate error counts and coverage distribution at unique only sites
	//

	if (settings.do_step(settings.error_counts_done_file_name, "Tabulating error counts"))
	{
		create_path(settings.error_calibration_path);

		string reference_fasta_file_name = settings.reference_fasta_file_name;
		string reference_bam_file_name = settings.reference_bam_file_name;

		//my $cbreseq = settings.ctool("cbreseq");

		// deal with distribution or error count keys being undefined...
		string coverage_fn = settings.file_name(settings.unique_only_coverage_distribution_file_name, "@", "");
		string outputdir = dirname(coverage_fn) + "/";

		uint32_t num_read_files = summary.sequence_conversion.reads.size();
		uint32_t num_qual = summary.sequence_conversion.max_qual + 1;

		error_count(
      summary,
			reference_bam_file_name, // bam
			reference_fasta_file_name, // fasta
			settings.error_calibration_path, // output
			settings.read_files.base_names(), // readfile
			true, // coverage
			true, // errors
			settings.base_quality_cutoff, // minimum quality score
			"read_set=" + to_string(num_read_files) + ",obs_base,ref_base,quality=" + to_string(num_qual) // covariates
		);

		settings.done_step(settings.error_counts_done_file_name);
	}


	//
	// Calculate error rates
	//sub error_rates {}
	//

	create_path(settings.output_path); //need output for plots
  create_path(settings.output_calibration_path);

	if (settings.do_step(settings.error_rates_done_file_name, "Re-calibrating base error rates"))
	{
		if (!settings.no_deletion_prediction)
			CoverageDistribution::analyze_unique_coverage_distributions(
                                                                  settings, 
                                                                  summary, 
                                                                  ref_seq_info,
                                                                  settings.unique_only_coverage_plot_file_name, 
                                                                  settings.unique_only_coverage_distribution_file_name
                                                                  );

    //Coverage distribution user option --deletion-coverage-propagation-cutoff
    if (settings.deletion_coverage_propagation_cutoff) {
      if (settings.deletion_coverage_propagation_cutoff < 1) {
        for (uint32_t i = 0; i < ref_seq_info.size(); i++) {
          string seq_id = ref_seq_info[i].m_seq_id;
          double average = summary.unique_coverage[seq_id].average;
          double &deletion_coverage_propagation_cutoff = summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff;

          deletion_coverage_propagation_cutoff = average * settings.deletion_coverage_propagation_cutoff;
        }
      } else if (settings.deletion_coverage_propagation_cutoff >= 1) {
        for (uint32_t i = 0; i < ref_seq_info.size(); i++) {
          string seq_id = ref_seq_info[i].m_seq_id;
          double &deletion_coverage_propagation_cutoff = summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff;

          deletion_coverage_propagation_cutoff = settings.deletion_coverage_propagation_cutoff;
        }
      }
    }

    //Coverage distribution user option --deletion-coverage-seed-cutoff
    if (settings.deletion_coverage_seed_cutoff) {
      if (settings.deletion_coverage_seed_cutoff < 1) {
        for (uint32_t i = 0; i < ref_seq_info.size(); i++) {
          string seq_id = ref_seq_info[i].m_seq_id;
          double average = summary.unique_coverage[seq_id].average;
          double &deletion_coverage_seed_cutoff = summary.unique_coverage[seq_id].deletion_coverage_seed_cutoff;

          deletion_coverage_seed_cutoff = average * settings.deletion_coverage_seed_cutoff;
        }
    } else if (settings.deletion_coverage_seed_cutoff >= 1) {
        for (uint32_t i = 0; i < ref_seq_info.size(); i++) {
          string seq_id = ref_seq_info[i].m_seq_id;
          double &deletion_coverage_seed_cutoff = summary.unique_coverage[seq_id].deletion_coverage_seed_cutoff;

          deletion_coverage_seed_cutoff = settings.deletion_coverage_seed_cutoff;
        }
      }
    }
		string command;
		for (uint32_t i = 0; i<settings.read_files.size(); i++) {
			string base_name = settings.read_files[i].base_name();
			string error_rates_base_qual_error_prob_file_name = settings.file_name(settings.error_rates_base_qual_error_prob_file_name, "#", base_name);
			string plot_error_rates_r_script_file_name = settings.plot_error_rates_r_script_file_name;
			string plot_error_rates_r_script_log_file_name = settings.file_name(settings.plot_error_rates_r_script_log_file_name, "#", base_name);
			string error_rates_plot_file_name = settings.file_name(settings.error_rates_plot_file_name, "#", base_name);
			command = "R --vanilla in_file=" + error_rates_base_qual_error_prob_file_name + " out_file=" + error_rates_plot_file_name + " < " + plot_error_rates_r_script_file_name + " > " + plot_error_rates_r_script_log_file_name;
			SYSTEM(command);
		}

		summary.unique_coverage.store(settings.error_rates_summary_file_name);
		settings.done_step(settings.error_rates_done_file_name);
	}
	summary.unique_coverage.retrieve(settings.error_rates_summary_file_name);

	//
	// Make predictions of point mutations, small indels, and large deletions
	//sub mutation_prediction {}
	//

	if (!settings.no_mutation_prediction)
	{
		create_path(settings.mutation_identification_path);

		if (settings.do_step(settings.mutation_identification_done_file_name, "Read alignment mutations"))
		{
			string reference_fasta_file_name = settings.reference_fasta_file_name;
			string reference_bam_file_name = settings.reference_bam_file_name;

			string coverage_fn = settings.file_name(settings.unique_only_coverage_distribution_file_name, "@", "");
			string error_dir = dirname(coverage_fn) + "/";
			string this_predicted_mutation_file_name = settings.file_name(settings.predicted_mutation_file_name, "@", "");
			string output_dir = dirname(this_predicted_mutation_file_name) + "/";
			string ra_mc_cGenomeDiff_file_name = settings.ra_mc_cGenomeDiff_file_name;
			string coverage_tab_file_name = settings.file_name(settings.complete_coverage_text_file_name, "@", "");
			string coverage_dir = dirname(coverage_tab_file_name) + "/";

			// It is important that these are in consistent order with the fasta file!!
      vector<double> deletion_propagation_cutoffs;
      vector<double> deletion_seed_cutoffs;
      for (uint32_t i = 0; i < ref_seq_info.size(); i++) {
        deletion_propagation_cutoffs.push_back(summary.unique_coverage[ref_seq_info[i].m_seq_id].deletion_coverage_propagation_cutoff);
        deletion_seed_cutoffs.push_back(summary.unique_coverage[ref_seq_info[i].m_seq_id].deletion_coverage_seed_cutoff);
        //cout << ref_seq_info[i].m_seq_id << " " << to_string(deletion_propagation_cutoffs.back()) << " " << to_string(deletion_seed_cutoffs.back()) << endl;
      }

			identify_mutations(
				reference_bam_file_name,
				reference_fasta_file_name,
				error_dir,
				ra_mc_cGenomeDiff_file_name,
				output_dir,
				settings.read_file_names,
				coverage_dir,
        deletion_propagation_cutoffs,
        deletion_seed_cutoffs,
				settings.mutation_log10_e_value_cutoff, // mutation_cutoff
				!settings.no_deletion_prediction,
				settings.polymorphism_prediction,
				settings.base_quality_cutoff, // minimum_quality_score
				settings.polymorphism_log10_e_value_cutoff, // polymorphism_cutoff
				settings.polymorphism_frequency_cutoff, //polymorphism_frequency_cutoff
				error_dir + "/error_rates.tab",
				false //per_position_file
			);

			settings.done_step(settings.mutation_identification_done_file_name);
		}

    // extra processing for polymorphisms
		if (settings.polymorphism_prediction && settings.do_step(settings.polymorphism_statistics_done_file_name, "Polymorphism statistics"))
		{
			ref_seq_info.polymorphism_statistics(settings, summary);
			settings.done_step(settings.polymorphism_statistics_done_file_name);
		}
	}

	//rewire which GenomeDiff we get data from if we have the elaborated polymorphism_statistics version
	 if (settings.polymorphism_prediction)
		settings.ra_mc_cGenomeDiff_file_name = settings.polymorphism_statistics_ra_mc_cGenomeDiff_file_name;

  create_path(settings.evidence_path); //need output for plots

	if (settings.do_step(settings.output_done_file_name, "Output"))
	{
		///
		// Output Genome Diff File
		//sub cGenomeDiff_output {}
		///
		cerr << "Creating merged genome diff evidence file..." << endl;

		// merge all of the evidence GenomeDiff files into one...
		create_path(settings.evidence_path);
    cGenomeDiff jc_gd(settings.jc_cGenomeDiff_file_name);
    cGenomeDiff ra_mc_gd(settings.ra_mc_cGenomeDiff_file_name);
        
    cGenomeDiff merged_cGenomeDiff;
    cGenomeDiff evidence_gd = cGenomeDiff::fast_merge(jc_gd, ra_mc_gd);
		evidence_gd.write(settings.evidence_cGenomeDiff_file_name);

		// predict mutations from evidence in the GenomeDiff
		cerr << "Predicting mutations from evidence..." << endl;

		MutationPredictor mp(ref_seq_info);
    cGenomeDiff mpgd(settings.evidence_cGenomeDiff_file_name);
    mp.predict(settings, mpgd, summary.sequence_conversion.max_read_length, summary.sequence_conversion.avg_read_length);

    // Add addition values from breseq pipleline to genome diff here.
    if (settings.values_to_gd){
       for (storable_map<string, Coverage>::iterator it = summary.unique_coverage.begin();
            it != summary.unique_coverage.end(); it ++) {
         const string& seq_id = it->first;
//         mpgd.add_breseq_data("junction_accept_score_cutoff:" + seq_id,
//                              to_string(summary.unique_coverage[seq_id].junction_accept_score_cutoff));
       }

    }

    mpgd.write(settings.final_cGenomeDiff_file_name);


    cGenomeDiff gd(settings.final_cGenomeDiff_file_name);
		//#unlink $evidence_cGenomeDiff_file_name;

		//
		// Mark lowest RA evidence items as no-show, or we may be drawing way too many alignments
		//

    vector<gd_entry_type> ra_types = make_list<gd_entry_type>(RA);
    list<counted_ptr<cDiffEntry> > ra = gd.filter_used_as_evidence(gd.list(ra_types));

    ra.remove_if(cDiffEntry::frequency_less_than_two_or_no_show());
    ra.sort(cDiffEntry::by_scores(make_list<string>("quality")));

    uint32_t i=0;
    for(diff_entry_list_t::iterator item = ra.begin(); item != ra.end(); item++)
		{
      cDiffEntry& ra_item = **item;
      if (!ra_item.entry_exists(FREQUENCY)) continue;
      
      if ( from_string<double>(ra_item[FREQUENCY]) != 0.0 ) continue;
      if ( from_string<double>(ra_item[FREQUENCY]) != 1.0 ) continue;
        
      if (i++ > settings.max_rejected_polymorphisms_to_show)
        ra_item[NO_SHOW] = "1";
		}

		// require a certain amount of coverage
    diff_entry_list_t new_ra = gd.filter_used_as_evidence(gd.list(ra_types));
    for (diff_entry_list_t::iterator item = new_ra.begin(); item != new_ra.end(); item++)
		{
			vector<string> top_bot = split((**item)["tot_cov"], "/");
			uint32_t top = from_string<int32_t>(top_bot[0]);
			uint32_t bot = from_string<int32_t>(top_bot[1]);
			if (top + bot <= 2) (**item)["no_show"] = "1";
		}
		
		//
		// Mark lowest scoring reject junctions as no-show
		//
    vector<gd_entry_type> jc_types = make_list<gd_entry_type>(JC);
    
    diff_entry_list_t jc = gd.filter_used_as_evidence(gd.list(jc_types));
	  
    // This is the correct way to erase from a list! @MDS
    for (diff_entry_list_t::iterator it = jc.begin(); it != jc.end(); )
    {
      if (((*it)->number_reject_reasons() == 0))
        it = jc.erase(it);
      else
        it++;
    }
    
    jc.sort(cDiffEntry::by_scores(make_list<diff_entry_key_t>("pos_hash_score")("total_non_overlap_reads")));

    i = 0;
    for (diff_entry_list_t::iterator it = jc.begin(); it != jc.end(); it++ )
		{
			if (i++ >= settings.max_rejected_junctions_to_show)
				(**it)["no_show"] = "1";
		}


    //
		// Annotate mutations
		//
		cerr << "Annotating mutations..." << endl;
		ref_seq_info.annotate_mutations(gd);
    
		//
		// Plot coverage of genome and large deletions
		//
		cerr << "Drawing coverage plots..." << endl;
		output::draw_coverage(settings, ref_seq_info, gd);
    
		//
		// Create evidence files containing alignments and coverage plots
		//
		if (!settings.no_alignment_generation)
			output::Evidence_Files(settings, gd);

		///
		// HTML output
		///

		cerr << "Creating index HTML table..." << endl;

		output::html_index(settings.index_html_file_name, settings, summary, ref_seq_info, gd);
		output::html_marginal_predictions(settings.marginal_html_file_name, settings, summary, ref_seq_info, gd);
        
		// record the final time and print summary table
		settings.record_end_time("Output");

		output::html_statistics(settings.summary_html_file_name, settings, summary, ref_seq_info);

		settings.done_step(settings.output_done_file_name);
	}
  
  return 0;
}

/*! breseq commands
 
    First argument is a command that should be removed from argv.
 
 */
int main(int argc, char* argv[]) {

	// Extract the sub-command argument
	string command;
	char* argv_new[argc];
	int argc_new = argc - 1;

  // Print out our generic header
  Settings::command_line_run_header();
  
	if (argc > 1) {

		command = argv[1];
		argv_new[0] = argv[0];
		for (int32_t i = 1; i < argc; i++)
			argv_new[i] = argv[i + 1];

	} else {
    breseq_default_action(argc, argv); // Gives default usage in this case.
    return -1; 
	}

	// Pass the command to the proper handler
	command = to_upper(command);
  if (command == "CONVERT_FASTQ") {
		return do_convert_fastq(argc_new, argv_new);
	} else if (command == "CALCULATE_TRIMS") {
		return do_calculate_trims(argc_new, argv_new);
	} else if (command == "CONVERT_GENBANK") {
		return do_convert_genbank(argc_new, argv_new);
	} else if (command == "TABULATE_CL") {
		return do_tabulate_contingency_loci(argc_new, argv_new);
	} else if (command == "GRAPH_CL") {
		return do_graph_contingency_loci(argc_new, argv_new);
	} else if (command == "ERROR_COUNT") {
		return do_error_count(argc_new, argv_new);
	} else if (command == "IDENTIFY_MUTATIONS") {
		return do_identify_mutations(argc_new, argv_new);
	} else if (command == "PREPROCESS_ALIGNMENTS") {
		return do_preprocess_alignments(argc_new, argv_new);
	} else if (command == "IDENTIFY_CANDIDATE_JUNCTIONS") {
		return do_identify_candidate_junctions(argc_new, argv_new);
  } else if (command == "PREDICT_MUTATIONS") {
		return do_predict_mutations(argc_new, argv_new);
	} else if (command == "GD2GVF") {
    return do_convert_gvf(argc_new, argv_new);
  } else if (command == "VCF2GD") {
    return do_convert_gd( argc_new, argv_new);
  } else if (command == "BAM2ALN") {
    return do_bam2aln( argc_new, argv_new);  
  } else if (command == "BAM2COV") {
    return do_bam2cov( argc_new, argv_new);    
  } else if ((command == "APPLY") || (command == "MUTATE")) {
    return do_mutate(argc_new, argv_new);    
  } else if (command == "SUBTRACT") {
    return do_subtract(argc_new, argv_new);    
  } else if (command == "UNION" || command == "MERGE") {
    return do_merge(argc_new, argv_new);    
  } else if (command == "NOT_EVIDENCE") {
    return do_not_evidence(argc_new, argv_new);
  } else if (command == "COMPARE") {
    return do_compare(argc_new, argv_new);
  } else if (command == "INTERSECTION") {
    return do_intersection(argc_new, argv_new);
  } else if (command == "ANNOTATE") {
    return do_annotate(argc_new, argv_new);
  } else if (command == "SIMULATE-READ") {
    return do_simulate_read(argc_new, argv_new);
  } else {
    // Not a sub-command. Use original argument list.
    return breseq_default_action(argc, argv);
  }
	return -1;
}
