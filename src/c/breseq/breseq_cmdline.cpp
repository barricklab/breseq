 /*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2012 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#include <iostream>
#include <string>
#include <vector>

#include "config.h"

#include "libbreseq/anyoption.h"
#include "libbreseq/alignment_output.h"
#include "libbreseq/reference_sequence.h"
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
#include "libbreseq/summary.h"
#include "libbreseq/contingency_loci.h"
#include "libbreseq/mutation_predictor.h"
#include "libbreseq/rna_seq.h"
#include "libbreseq/output.h"

using namespace breseq;
using namespace std;


/*! bam2aln
 Draw HTML alignment from BAM
 */
int do_bam2aln(int argc, char* argv[]) {
  
  // setup and parse configuration options:
	AnyOption options("Usage: breseq BAM2ALN [--bam=reference.bam --fasta=reference.fasta --region=accession:start-end --output=alignment.html --max-reads=200] accession:start-end");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  ("bam,b", "bam file containing sequences to be aligned", "data/reference.bam")
	("fasta,f", "FASTA file of reference sequence", "data/reference.fasta")
  ("output,o", "output to file [region.html]", "alignment.html")
  ("format", "Output format: txt or html", "html")
  ("region,r", "region to print (accession:start-end), may also be provided as unnamed arguments at the end of the command line.")
  ("max-reads,n", "maximum number of reads to show in alignment", 200)
  ("repeat", "show reads with multiple best matches in reference", TAKES_NO_ARGUMENT)
  ("quality-score-cutoff,c", "quality score cutoff", 0)
  ("stdout", "write output to stdout", TAKES_NO_ARGUMENT)
  .processCommandArgs(argc, argv);  
  
	// make sure that the config options are good:
	if(options.count("help")
     || !file_exists(options["fasta"].c_str())
     || !file_exists(options["bam"].c_str()) )
  {
		options.printUsage();
		return -1;
	}
  
  vector<string> region_list;
  if (options.count("region"))
    region_list= from_string<vector<string> >(options["region"]);
  
  // Also take regions off the command line
  for (int32_t i = 0; i < options.getArgc(); i++)
  {
    string region = options.getArgv(i);
    region_list.push_back(region);
  }  
  
  if (!region_list.size()) {
    options.addUsage("");
    options.addUsage("You must supply the --region option for input.");
    options.printUsage();
    return -1;
  }  
  
  for(uint32_t j = 0; j < region_list.size(); j++)
  {
    // clean commas
    region_list[j] = substitute(region_list[j], ",", "");
    
    cerr << "Creating alignment for region: " << region_list[j] << endl;
    
    // Generate Alignment!
    alignment_output ao(
                        options["bam"],
                        options["fasta"],
                        from_string<uint32_t>(options["max-reads"]),
                        from_string<uint32_t>(options["quality-score-cutoff"]),
                        1,
                        options.count("repeat")
                        );
    
    string default_file_name = region_list[j];
    
    string output_string;
    if (options["format"] == "html") {
      output_string = ao.html_alignment(region_list[j]);
      default_file_name += ".html";
    }
    else if (options["format"] == "txt") {
      output_string = ao.text_alignment(region_list[j]);
      default_file_name += ".txt";
    }
    else {
      options.addUsage("");
      options.addUsage("Unknown format requested: " + options["format"]);
      options.printUsage();
      return -1;
    }
      
    if (options.count("stdout"))
    {
      cout << output_string << endl;
    }
    else
    {
      ///Write to file
      string file_name = default_file_name;
            
      if (options.count("output"))
      {
        file_name = options["output"];
        if(region_list.size() > 1)  {
          file_name = (split(options["output"], "."))[0] + "_" + region_list[j] + ".html";  }
      }
      
      ofstream myfile (file_name.c_str());
      if (myfile.is_open())
      {
        myfile << output_string;
        myfile.close();
      }
      else cerr << "Unable to open file";
    }
  }
    
  return 0;
}

/*! bam2cov
 Draw HTML coverage from BAM
 */
int do_bam2cov(int argc, char* argv[]) {
  // setup and parse configuration options:
	AnyOption options("Usage: bam2aln -b reference.bam -f reference.fasta [--format PNG -o output.png] seq_id1:start1-end1 [seq_id2:start2-end2 ...]");
	options
  ("help,h", "Produce advanced help message", TAKES_NO_ARGUMENT)
  // required options
  ("bam,b", "BAM file containing sequences to be aligned", "data/reference.bam")
	("fasta,f", "FASTA file of reference sequence", "data/reference.fasta")
  // options controlling what files are output
  ("output,o", "If there is only one output file, then this option specifies the name of the output file. Otherwise, the name of a path that will contain all output files. By default, the base name of the output file for a region is seq_id:start-end.")
  ("format", "Format of output plot(s): PNG or PDF", "PNG")
  ("table,t", "Create tab delimited file of coverage instead of a plot", TAKES_NO_ARGUMENT)
  ("total-only,1", "Only plot/tabulate total coverage, not per strand coverage", TAKES_NO_ARGUMENT)
  ("show-average,a", "Show the average coverage across the reference sequence as a horizontal line. Only possible if used in the main output directory of breseq output.", TAKES_NO_ARGUMENT)
  ("resolution,p", "Number of positions to output coverage information for in interval (0=ALL)", 600)
  // which regions to create files for
  ("tile-size", "In tiling mode, the size of each tile", 0)
  ("tile-overlap", "In tiling mode, overlap between adjacent tiles (1/2 of this is added to each side of every tile)", 0)
//  ("read_start_output,r", "file name for table file binned by read start bases (DEFAULT: OFF)")
//  ("gc_output,g", "create additional table file binned by GC content of reads (DEFAULT: OFF)")
  // options controlling information that is output
  .processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Create a coverage plot or table for the specified region or regions.");

  // Print advanced help
  if (options.count("help"))
  {
    options.printAdvanced();
    exit(-1);
  }
  
  // make sure that the required config options are good:
	if(   !file_exists(options["fasta"].c_str() )
     || !file_exists(options["bam"].c_str()   ) )
  {
		options.printUsage();
		return -1;
	}
  
  vector<string> region_list;
  
  // Take regions off the command line
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
    WARN("Tiling mode activated. Ignoring " + to_string(region_list.size()) + " regions that were specified.");
    region_list.clear();
  }
  
  // Did they specify some regions?
  if(!tiling_mode && (region_list.size() == 0)) {
    options.addUsage("");
    options.addUsage("No regions specified.");
    options.printUsage();
    return -1;
  }
  
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
  co.output_format( options["format"] );
  if (options.count("show-average")) {
    co.show_average( options.count("show-average"), settings.error_rates_summary_file_name );
  }
  
  // create regions that tile the genome
  if (tiling_mode)
  {
    int32_t tile_size = from_string<int32_t>(options["tile-size"]);
    int32_t tile_overlap = from_string<int32_t>(options["tile-overlap"]);
    tile_overlap = static_cast<int32_t>(floor( static_cast<double>(tile_overlap) / 2.0));

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
  
  // Create output path if necessary
  if ((region_list.size() > 1) && (!options["output"].empty())) {
    create_path(options["output"]);
  }
  
  for(vector<string>::iterator it = region_list.begin(); it!= region_list.end(); it++)
  {
// these are experimental... additional table files
//    if (options.count("read_start_output"))
//      co.read_begin_output_file_name(options["read_start_output"]);
//    if (options.count("gc_output"))
//      co.gc_output_file_name(options["gc_output"]);
    
    // clean commas
    *it = substitute(*it, ",", "");
    
    string file_name;
    if (region_list.size() == 1) {
      if (options["output"].empty()) {
        file_name = *it;
      } else {
        file_name = options["output"];
      }
    } else {
      if (options["output"].empty()) {
        file_name = *it;
      } else {
        file_name = options["output"] + "/" + *it;
      }
    }
    
    if (options.count("table")) {
      cout << "Tabulating coverage for region: " << *it << endl;
      co.table(*it, file_name + ".tab", from_string<uint32_t>(options["resolution"]));
    } else {
      cout << "Plotting coverage for region: " << *it << endl;
      co.plot(*it, file_name + "." + to_lower(options["format"]) , from_string<uint32_t>(options["resolution"]));
    }
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
  ("reverse-complement,r", "reverse complement all reads and add _RC to their names", TAKES_NO_ARGUMENT)
   .processCommandArgs(argc, argv);
	
	// make sure that the config options are good:
	if(options.count("help")
		 || !options.count("input")
     || !options.count("output")
     || !options.count("in-format")
     || !options.count("out-format")
		 ) {
    options.addUsage("Valid formats are 'SANGER', 'SOLEXA', 'ILLUMINA_1.3+'");
		options.printUsage();
		return -1;
	}                       
    
  convert_fastq(options["input"], options["output"], options["in-format"], options["out-format"], options.count("reverse-complement"));

  return 0;
}

int do_convert_genbank(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq CONVERT-REFERENCE --input <reference> [--fasta <output.fasta>] OR [--gff <output.gff>]");
	options
		("input,i", "Input reference file(s). (Format: Fasta, GFF, GenBank)")
		("fasta",   "FASTA format output path.")
    ("gff",     "GFF format output path.")
	.processCommandArgs(argc, argv);
	
	// make sure that the config options are good:
  if (!options.count("input")) {
    options.addUsage("No input reference file given.");
		options.printUsage();
		return -1;
	}

  if (!options.count("fasta") && !options.count("gff")) {
    options.addUsage("No output path and format given.");
		options.printUsage();
		return -1;
  }

  
  cReferenceSequences refs;
  refs.LoadFiles(from_string<vector<string> >(options["input"]));

  if (options.count("fasta")) refs.WriteFASTA(options["fasta"]);

  if (options.count("gff")) refs.WriteGFF(options["gff"]);

	
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
    Settings settings;
		error_count(
                settings,
                summary,
                options["bam"],
                options["fasta"],
                options["output"],
                split(options["readfile"], "\n"),
                options.count("coverage"),
                options.count("errors"),
                false,
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



/*! Contingency Loci
 
 Analyze lengths of homopolymer repeats in mixed samples.
 
 */
int do_tabulate_contingency_loci(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq TABULATE_CL --bam <sequences.bam> --fasta <reference.fasta> --reference <reference.gff> --output <path> [--loci <loci.txt> --strict]");
	options
  ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
  ("bam,b", "bam file containing sequences to be aligned", "data/reference.bam")
  ("fasta,f", "FASTA file of reference sequence", "data/reference.fasta")
  ("reference,r","reference sequence file with annotation", "data/reference.gff3")
  ("output,o", "output file", "contingency_loci.tab")
  ("loci,l", "Contingency loci coordinates", "")
  ("strict,s", "exclude non-perfect matches in surrounding 5 bases", TAKES_NO_ARGUMENT)
	.processCommandArgs(argc, argv);
  
	// make sure that the config options are good:
	if(options.count("help")
		 ) {
		options.printUsage();
		return -1;
	}                       
  
	// attempt to calculate error calibrations:
	try {
    
    analyze_contingency_loci(
                             options["bam"],
                             options["fasta"],
                             from_string<vector<string> >(options["reference"]),
                             options["output"],
                             options["loci"],
                             options.count("strict")
                             );
	} catch(...) {
		// failed; 
		return -1;
	}
	
	return 0;
}

int do_analyze_contingency_loci_significance( int argc, char* argv[]){
    // setup and parse configuration options:
	AnyOption options("Usage: breseq CL_SIGNIFICANCE --output <path> --loci <loci.txt> ");
	options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("output,o", "output file", "contingency_loci.tab")
    ("loci,l", "Contingency loci files", "")
	.processCommandArgs(argc, argv);
    
  vector<string> strain_files = from_string<vector<string> >(options["loci"]);
  analyze_contingency_loci_significance(
                                       options["output"],
                                       strain_files
                                       );
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
    settings.read_files.Init(from_string<vector<string> >(options["read-file"]), false);
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








int do_simulate_read(int argc, char *argv[])
{
  AnyOption options("Usage: breseq SIMULATE-READ -g <genome diff> -r <reference file> -c <average coverage> -o <output file>");

  options
  ("genome_diff,g", "Genome diff file.")
  ("reference,r",   "Reference file for input.")
  ("coverage,c",    "Average coverage value to simulate.", static_cast<uint32_t>(80))
  ("length,l",      "Read length to simulate.", static_cast<uint32_t>(50))
  ("output,o",      "Output fastq file name.")
  ("gff3,3",        "Output Applied GFF3 File. (Flag)", TAKES_NO_ARGUMENT)
  ("verbose,v",     "Verbose Mode (Flag)", TAKES_NO_ARGUMENT)
  ("paired-ends",   "Will create two read files *_1.fastq, *_2.fastq that are mate-paired.(Flag)", TAKES_NO_ARGUMENT)
  ("mean",          "Mean fragment size to use for pair ended reads.", static_cast<uint32_t>(200))
  ("stdev",         "Standard deviation of fragment size to use for pair ended reads.", static_cast<uint32_t>(20))
  ("seed",          "Seed value to use for random number generation.")
  ;
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Takes a genome diff file and applies the mutations to the reference.");
  options.addUsage("Then using the applied reference, it simulates reads based on it.");
  options.addUsage("Output is a .fastq that if run normally against the reference");
  options.addUsage("should produce the same GenomeDiff file you entered.");
  
  if(argc == 1)  {
    options.printUsage();
    return -1;  }
  
  if (!options.count("genome_diff")) {
    options.addUsage("");
    options.addUsage("You must supply the --genome_diff option for input.");
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

  if (options.count("seed")) {
    cSimFastqSequence::SEED_VALUE = from_string<int32_t>(options["seed"]);
  }

  const bool verbose = options.count("verbose");

//! Step: Load reference sequence file.
  cReferenceSequences ref_seq_info;
  cReferenceSequences new_ref_seq_info;
  const string &ref_file_name = options["reference"];
  ref_seq_info.LoadFile(ref_file_name);
  new_ref_seq_info.LoadFile(ref_file_name);


  //! Step: Apply genome diff mutations to reference sequence.
  const string &gd_file_name = options["genome_diff"];
  cGenomeDiff gd(gd_file_name);

  gd.apply_to_sequences(ref_seq_info, new_ref_seq_info, verbose);
  
  //! Write applied GFF3 file if requested.
  if(options.count("gff3"))new_ref_seq_info.WriteGFF(options["output"] + ".gff3", options.count("verbose"));


  const cAnnotatedSequence &sequence = new_ref_seq_info[0];
  uint32_t coverage = from_string<uint32_t>(options["coverage"]);
  uint32_t read_size = from_string<uint32_t>(options["length"]);

  uint32_t n_reads = (sequence.get_sequence_size() / read_size) * coverage;

  if (options.count("paired-ends")) {
    cString output = options["output"];

    string extension = output.get_file_extension();
    output.remove_ending(extension);

    string pair_1_file_name = output + "_1" + extension;
    string pair_2_file_name = output + "_2" + extension;

    uint32_t mean  = from_string<uint32_t>(options["mean"]);
    uint32_t stdev = from_string<uint32_t>(options["stdev"]);

    cSimFastqSequence::simulate_paired_ends(sequence,
                                            n_reads,
                                            read_size,
                                            mean,
                                            stdev,
                                            pair_1_file_name,
                                            pair_2_file_name,
                                            verbose);

  } else {
    cSimFastqSequence::simulate_single_ends(sequence,
                                            n_reads,
                                            read_size,
                                            options["output"],
                                            verbose);
  }

  return 0;
}


/*
 * Function: do_copy_number_variation
 * --------------------------------
 * Called if "CNV" is passed as argument when invoking breseq on command line.
 * Ex: breseq CNV -r lambda.gbk lambda.fastq
 * NOTE: "breseq CNV" just runs the copy number variation part of the pipeline by itself (for testing), 
 * whereas "breseq --cnv" runs the whole pipeline including copy number variation.
 * BUG: Running "CNV" might not produce a .done file (i.e. copy_number_variation.done) in
 * the relevant output directory (i.e. .../09_copy_number_variation/).
 * BUG: Running "CNV" might not update the html output in the relevant directory,
 * (i.e. .../output/).
 */
int do_copy_number_variation(int argc, char *argv[])
{
    Settings settings(argc, argv);

    //(re)load the reference sequences from our converted files
    cReferenceSequences ref_seq_info;
    ref_seq_info.ReadGFF(settings.reference_gff3_file_name);

    // Where error rate summary data will be output
    Summary summary;
    summary.unique_coverage.retrieve(settings.error_rates_summary_file_name);

    // Create copy_number_variation directory
    create_path( settings.copy_number_variation_path );

    for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it) {
        
        // Sequence iterator, one sequence at a time from .fastq file
        cAnnotatedSequence& seq = *it;

        // Create filename: [genome].coverage.txt, this is LONG file where the coverage data is coming from
        string this_complete_coverage_text_file_name = settings.file_name(settings.complete_coverage_text_file_name, "@", seq.m_seq_id);

        // Create filename: [genome].tiled.tab, this is LONG file where the tiled coverage data will be output
        string this_tiled_complete_coverage_text_file_name = settings.file_name(settings.tiled_complete_coverage_text_file_name, "@", seq.m_seq_id);

        // Create filename: [genome].tiled_for_edging.tab, this is where the tiled
        // coverage data at a smaller tile_size will be output for use in improved
        // edge detection (aka "edging")
        string this_tiled_for_edging_text_file_name = settings.file_name(settings.tiled_for_edging_text_file_name, "@", seq.m_seq_id);
        
        // Generates [genome].tiled.tab, one line at a time with each for-loop,
        // where the avg coverage of each tile is calculated.
        CoverageDistribution::tile(settings.ignore_redundant_coverage,
                                   this_complete_coverage_text_file_name, // in_file_name
                                   this_tiled_complete_coverage_text_file_name, // out_file_name
                                   this_tiled_for_edging_text_file_name,
                                   settings.copy_number_variation_tile_size);

        // Create filename: [genome].ranges.tab, this is SHORT file used for ???
        // (contains: Start_Position, End_Position, T_Score, P_Value)
        string this_ranges_text_file_name = settings.file_name(settings.ranges_text_file_name, "@", seq.m_seq_id);

        // Get filename: [genome].history.tab, this is SHORT file used for ???
        // (contains: Start_Search, End_Search, Start_Position, End_Position, Start_Segment, End_Segment... 16 values total)
        string this_cnv_history_text_file_name = settings.file_name(settings.cnv_history_text_file_name, "@", seq.m_seq_id);

        // Generates [genome].ranges.tab & [genome].history.tab, one line at a time with each for-loop,
        CoverageDistribution::find_segments(settings,
                                            summary.unique_coverage[seq.m_seq_id].nbinom_mean_parameter,
                                            this_tiled_complete_coverage_text_file_name,
                                            this_tiled_for_edging_text_file_name, // tiled_for_edging_file_name (tiled_for_edging.tab)
                                            this_ranges_text_file_name,
                                            this_cnv_history_text_file_name
                                            );

        // Create filename: [genome].smoothed_ranges.tab, this is LONG file used for ???
        // (contains: Position, Smooth_Coverage)
        string this_smoothed_ranges_text_file_name = settings.file_name(settings.smoothed_ranges_text_file_name, "@", seq.m_seq_id);

        // Create filename: [genome].cnv_final.tab, this is SHORT file used for ???
        // (contains: Start_Position, End_Position, Z_Score, Greater_Than, Copy_Number)
        string this_final_cnv_file_name = settings.file_name(settings.final_cnv_text_file_name, "@", seq.m_seq_id);

        // Create filename: [genome].cn_evidence.gd, this is SHORT file used for ???
        // (contains no labels)
        string this_copy_number_variation_cn_genome_diff_file_name = settings.file_name(settings.copy_number_variation_cn_genome_diff_file_name, "@", seq.m_seq_id);

        // Generates [genome].smoothed_ranges.tab & [genome].cnv_final.tab & [genome].cn_evidence.gd,
        // one line at a time with each for-loop,
        CoverageDistribution::smooth_segments(settings,
                                              seq.m_seq_id,
                                              summary.unique_coverage[seq.m_seq_id].nbinom_mean_parameter, 
                                              this_tiled_for_edging_text_file_name, 
                                              this_ranges_text_file_name, 
                                              this_smoothed_ranges_text_file_name,
                                              this_final_cnv_file_name,
                                              this_copy_number_variation_cn_genome_diff_file_name
                                              );

    } // End of for loop

    return 0;
}

int do_periodicity(int argc, char *argv[])
{
	Settings settings(argc, argv);
  
  //(re)load the reference sequences from our converted files
  cReferenceSequences ref_seq_info;
  ref_seq_info.ReadGFF(settings.reference_gff3_file_name);
  
  Summary summary;
  summary.unique_coverage.retrieve(settings.error_rates_summary_file_name);
  
  CandidateJunctions::identify_candidate_junctions(settings, summary, ref_seq_info);
  
  //create directory
  create_path( settings.copy_number_variation_path );
  
  for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it)
  {
    cAnnotatedSequence& seq = *it;
    string this_complete_coverage_text_file_name = settings.file_name(settings.complete_coverage_text_file_name, "@", seq.m_seq_id);
    string this_periodicity_complete_coverage_text_file_name = settings.file_name(settings.periodicity_table_file_name, "@", seq.m_seq_id);
    CoverageDistribution::calculate_periodicity(
          this_complete_coverage_text_file_name, 
          this_periodicity_complete_coverage_text_file_name,
          settings.periodicity_method,
          settings.periodicity_start,
          settings.periodicity_end,
          settings.periodicity_step
          );
    
   }
  
  return 0;
}



int do_get_sequence(int argc, char *argv[])
{
  AnyOption options("Usage: breseq GET-SEQUENCE -r <reference> -o <output.fasta> -p <REL606:50-100>");
  options("reference,r",".gbk/.gff3/.fasta reference sequence file", "data/reference.fasta");
  options("output,o","output FASTA file");  
  options("position,p","Sequence ID:Start-End");
  options("reverse-complement,c","Reverse Complement (Flag)", TAKES_NO_ARGUMENT);
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Takes a reference sequence and outputs to FASTA a subset");
  options.addUsage("of the sequence.  Using '0' as the End position will set");  
  options.addUsage("it to the length of the relevant sequence.");
  
  if(argc == 1)  {
    options.printUsage();
    return -1;  }
  
  if(!file_exists(options["reference"].c_str()))  {
    options.addUsage("");
    options.addUsage("You must supply a valid --reference option for input.");
    options.addUsage("Could not find:");
    options.addUsage(options["reference"].c_str());
    options.printUsage();
    return -1;    
  }
   
  vector<string> region_list;  
  
  bool do_reverse_complement = options.count("reverse-complement");
  bool verbose = options.count("verbose") || !options.count("output");
  
  if (options.count("position"))  {
    region_list = from_string<vector<string> >(options["position"]);  }
  
  // Also take positions off the command line
  for (int32_t i = 0; i < options.getArgc(); i++)
  {
    string position = options.getArgv(i);
    region_list.push_back(position);
  }
  
  if (!region_list.size()) {
    options.addUsage("");
    options.addUsage("You must supply the --position option for input.");
    options.printUsage();
    return -1;
  }
  
  cReferenceSequences ref_seq_info, new_seq_info;
  ref_seq_info.LoadFiles(from_string<vector<string> >(options["reference"]));
  
  for(uint32_t j = 0; j < region_list.size(); j++)
  {    
    uint32_t replace_target_id, replace_start, replace_end;
    string seq_name = "";
    
    ref_seq_info.parse_region(region_list[j], replace_target_id, replace_start, replace_end);
    
    if(!replace_end)  {
      replace_end = ref_seq_info[replace_target_id].m_length;  
    }
    
    if (replace_start > replace_end) {
      cout << "START greater than END." << " Creating reverse complement" << endl;
      swap(replace_start, replace_end);
      do_reverse_complement = true;
    }
    
    if(verbose)
    {
      cout << "Sequence ID:      " << ref_seq_info[replace_target_id].m_seq_id << endl;
      cout << "Start Position:   " << replace_start << endl;
      cout << "End Position:     " << replace_end << endl;
      cout << "Genome Strand:    " << (do_reverse_complement ? "Bottom" : "Top") << endl;
    }
    
    ASSERT((uint32_t)ref_seq_info[replace_target_id].m_length >= replace_start && (uint32_t)ref_seq_info[replace_target_id].m_length >= replace_end,
           "START:\t" + to_string(replace_start) + "\n" +
           "END:\t" + to_string(replace_end) + "\n" +
           "SIZE:\t" + to_string(ref_seq_info[replace_target_id].m_length) + "\n" +
           "Neither START or END can be greater than the SIZE of " + ref_seq_info[replace_target_id].m_seq_id + ".");
    
    seq_name = ref_seq_info[replace_target_id].m_seq_id + ":" + to_string(replace_start) + "-" + to_string(replace_end);
    
    new_seq_info.add_new_seq(seq_name, "");
    cAnnotatedSequence& new_seq = new_seq_info[seq_name];
    new_seq.m_fasta_sequence = ref_seq_info[replace_target_id].m_fasta_sequence;    
    new_seq.m_fasta_sequence.m_name = seq_name;
    new_seq.m_fasta_sequence.m_sequence = to_upper(ref_seq_info.get_sequence_1(replace_target_id, replace_start, replace_end));
    if(do_reverse_complement) 
      new_seq.m_fasta_sequence.m_sequence = reverse_complement(new_seq.m_fasta_sequence.m_sequence);
    new_seq.m_seq_id = region_list[j];
    new_seq.m_length = new_seq.m_fasta_sequence.m_sequence.size();
    
    if(verbose)  {
      cout << new_seq.m_fasta_sequence << endl;  
    }
  }
  
  if(options.count("output"))  {
    new_seq_info.WriteFASTA(options["output"], verbose);  
  }  
  
  return 0;
}

int do_junction_polymorphism(int argc, char *argv[])
{
  AnyOption options("Usage: breseq JUNCTION-POLYMORPHISM -g input.gd [ -o output.gd ]");  
  options.addUsage("Counts the number of reads supporting a new junction versus supporting the reference junction.");
  options("genome-diff,g","Input Genome Diff file");  
  options("output,o","Output Genome Diff file", "output.gd");
  options("verbose,v","Verbose output", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  // Check options
  if ( !options.count("genome-diff") ) {
    options.printUsage();
    return -1;
  }

  ASSERT(false, "not implemented")
/* could re-implement with read length
  Settings settings;
  settings.verbose = options.count("verbose");
  cGenomeDiff gd(options["genome-diff"]);
  assign_junction_read_counts(settings, summary, gd);
  gd.write(options["output"]);
*/
  return 0;
}

int do_ref_aln()
{
  SYSTEM("samtools view -bt data/reference.fasta.fai 02_reference_alignment/*.sam > data/ref_aln.bam");
  SYSTEM("samtools sort data/ref_aln.bam data/ref_aln.sorted");
  SYSTEM("samtools index data/ref_aln.sorted.bam");
  
  return 0;
}

int do_assemble_unmatched(int argc, char* argv[])
{
  AnyOption options("Usage: breseq ASSEMBLE-UNMATCHED-PAIRS [-o unmatched_assembly] reads1.fastq [reads2.fastq ...]");  
  options.addUsage("Assembles the unmatched reads from a breseq run given paired-end or mate-paired data.");
  options.addUsage("This command must be run from the main results directory of a breseq run (i.e., it must contain a data directory),");
  options.addUsage("You must provide the exact set of original fastq read files used in the breseq run. It is assumed that each set of two read files, in order, contain the first and second reads from pairs.");
  options.addUsage("Output is in directory: unmatched_assembly");
  options("output,o","Main directory containing output from the breseq run. A directory within this called unmatched_assembly will be created for the output of this command.", ".");
  options("verbose,v","Verbose output", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  // Need empty setting object for trimming files.
  Summary summary;
  Settings settings(options["output"]);
  
  // Read the corresponding unmatched file and extract each pair from the original fastq files
  
  create_path(options["output"] + "/assemble_unmatched");
  string unmatched_paired_read_file_name_1 = options["output"] + "/assemble_unmatched/1.fastq";
  string unmatched_paired_read_file_name_2 = options["output"] + "/assemble_unmatched/2.fastq";
  
  cFastqFile unmatched_paired_read_file_1(unmatched_paired_read_file_name_1, fstream::out);
  cFastqFile unmatched_paired_read_file_2(unmatched_paired_read_file_name_2, fstream::out);
  
  cReadFiles read_files;
  vector<string> read_file_names;
  for (int32_t i = 0; i < options.getArgc(); i++)
  {
    string read_file_name = options.getArgv(i);
    read_file_names.push_back(read_file_name);
  }
  read_files.Init(read_file_names, false);
  
  ASSERT(read_file_names.size() % 2 == 0, "Number of read files provided is not even.");
  
  vector<string> read_file_base_names = read_files.base_names();
  for(vector<string>::iterator it=read_file_base_names.begin(); it!=read_file_base_names.end(); it++)
  {
    string read_file_base_name_1 = *it;
    it++;
    string read_file_base_name_2 = *it;

    string read_file_name_1 = read_files.base_name_to_read_file_name(read_file_base_name_1);
    string read_file_name_2 = read_files.base_name_to_read_file_name(read_file_base_name_2);
    
    string unmatched_read_file_name_1 = settings.file_name(settings.unmatched_read_file_name, "#", read_file_base_name_1);
    string unmatched_read_file_name_2 = settings.file_name(settings.unmatched_read_file_name, "#", read_file_base_name_2);

    
    cFastqFile read_file_1(read_file_name_1, fstream::in);
    cFastqFile read_file_2(read_file_name_2, fstream::in);

    cFastqFile unmatched_read_file_1(unmatched_read_file_name_1, fstream::in);
    cFastqFile unmatched_read_file_2(unmatched_read_file_name_2, fstream::in);
    
    cFastqSequence on_seq_1, on_seq_2, um_on_seq_1, um_on_seq_2;
    
    // Figure out format (but just check one of the files to do this)
    uint64_t original_num_reads;
    uint64_t original_num_bases;
    uint32_t max_read_length;
    uint8_t min_quality_score;
    uint8_t max_quality_score;
    string quality_format = cFastqQualityConverter::predict_fastq_file_format(read_file_name_1, original_num_reads, original_num_bases, max_read_length, min_quality_score, max_quality_score);
    
    cFastqQualityConverter fqc(quality_format, "SANGER");
    cFastqQualityConverter unmatched_fqc("SANGER", "SANGER");

    bool rf_1_ok = read_file_1.read_sequence(on_seq_1, fqc);
    bool rf_2_ok = read_file_2.read_sequence(on_seq_2, fqc);
    bool um_rf_1_ok = unmatched_read_file_1.read_sequence(um_on_seq_1, unmatched_fqc);
    bool um_rf_2_ok = unmatched_read_file_2.read_sequence(um_on_seq_2, unmatched_fqc);

    while (um_rf_1_ok || um_rf_2_ok) {
      
      // if we found the sequence in the original read file
      if (um_rf_1_ok && rf_1_ok && on_seq_1.identical(um_on_seq_1) ) {
        
        // write out both current normal seqs
        unmatched_paired_read_file_1.write_sequence(on_seq_1);
        unmatched_paired_read_file_2.write_sequence(on_seq_2);
        
        // advance 
        um_rf_1_ok = unmatched_read_file_1.read_sequence(um_on_seq_1, unmatched_fqc);
        if (um_rf_2_ok && rf_2_ok && on_seq_2.identical(um_on_seq_2))
          um_rf_2_ok = unmatched_read_file_2.read_sequence(um_on_seq_2, unmatched_fqc);
        rf_1_ok = read_file_1.read_sequence(on_seq_1, fqc);
        rf_2_ok = read_file_2.read_sequence(on_seq_2, fqc);
      }
      else if (um_rf_2_ok && rf_2_ok && on_seq_2.identical(um_on_seq_2) ) {
        // write out both current normal seqs
        unmatched_paired_read_file_1.write_sequence(on_seq_1);
        unmatched_paired_read_file_2.write_sequence(on_seq_2);
        
        // advance (we know 1 didn't match)
        um_rf_2_ok = unmatched_read_file_2.read_sequence(um_on_seq_2, unmatched_fqc);
        rf_1_ok = read_file_1.read_sequence(on_seq_1, fqc);
        rf_2_ok = read_file_2.read_sequence(on_seq_2, fqc);
      }
      else {
        // advance the matched reads
        rf_1_ok = read_file_1.read_sequence(on_seq_1, fqc);
        rf_2_ok = read_file_2.read_sequence(on_seq_2, fqc);
      }
    }
  }
  
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
    
    // Do a quick load of the file to detect formatting errors.
    if (settings.user_junction_genome_diff_file_name != "") {
      cGenomeDiff gd(settings.user_junction_genome_diff_file_name); 
    }
    
    // Load all of the reference sequences and convert to FASTA and GFF3
    conv_ref_seq_info.LoadFiles(settings.all_reference_file_names);
    conv_ref_seq_info.WriteFASTA(settings.reference_fasta_file_name);
    conv_ref_seq_info.WriteGFF(settings.reference_gff3_file_name);

    // No conversion if already is sam mode
    if (!settings.aligned_sam_mode) {
      
      //Check the FASTQ format and collect some information about the input read files at the same time
      uint32_t overall_max_read_length = UNDEFINED_UINT32;
      uint32_t overall_max_qual = 0;

      s.num_reads = 0;
      s.num_bases = 0;
      
      // If limiting reads results in skipping later file, the index of the first one to delete
      uint32_t read_file_start_delete_index = 0; 
      
      uint64_t read_file_base_limit = floor(settings.read_file_coverage_fold_limit * static_cast<double>(conv_ref_seq_info.total_length())); 
      
      for (uint32_t i = 0; i < settings.read_files.size(); i++)
      {
        string base_name = settings.read_files[i].m_base_name;
        cerr << "  READ FILE::" << base_name << endl;
        
        // If we have reached the read limit or within some number of it -- delete further read files 
        // The 1000, is so that we have enough bases counted in a read file to fit the error rates.
        if ((settings.read_file_coverage_fold_limit) && (i!= 0) && (s.num_bases + 1000 >= read_file_base_limit)) {
          cerr << "  ::SKIPPED DUE TO REACHING COVERAGE LIMIT::" << endl;
          if (!read_file_start_delete_index) read_file_start_delete_index = i;
          continue;
        }
        
        string fastq_file_name = settings.base_name_to_read_file_name(base_name);
        string convert_file_name =  settings.file_name(settings.converted_fastq_file_name, "#", base_name);

        // Parse output
        Summary::AnalyzeFastq s_rf = normalize_fastq(fastq_file_name, convert_file_name, i+1, settings.quality_score_trim, !settings.no_read_filtering, s.num_bases, read_file_base_limit);
        settings.track_intermediate_file(settings.alignment_correction_done_file_name, convert_file_name);
        
        // Save the converted file name -- have to save it in summary because only that
        // is reloaded if we skip this step.
        s.converted_fastq_name[base_name] = s_rf.converted_fastq_name;

        // Record statistics
        if ((overall_max_read_length == UNDEFINED_UINT32) || (s_rf.max_read_length > overall_max_read_length))
          overall_max_read_length = s_rf.max_read_length;
        if ((overall_max_qual == UNDEFINED_UINT32) || (s_rf.max_quality_score > overall_max_qual))
          overall_max_qual = s_rf.max_quality_score;
        s.num_reads += s_rf.num_reads;
        s.original_num_reads += s_rf.original_reads;
        s.num_bases += s_rf.num_bases;
        s.original_num_bases += s_rf.original_num_bases;        

        s.reads[base_name] = s_rf;
      }
      if (read_file_start_delete_index) settings.read_files.resize(read_file_start_delete_index);
      
      s.avg_read_length = static_cast<double>(s.num_bases) / static_cast<double>(s.num_reads);
      s.max_read_length = overall_max_read_length;
      s.max_qual = overall_max_qual;
      summary.sequence_conversion = s;
      
      // Print totals
      cerr << "  ::TOTAL::" << endl;
      cerr << "    Original reads: " << s.original_num_reads << " bases: " << s.original_num_bases << endl;
      cerr << "    Analyzed reads: " << s.num_reads << " bases: " << s.num_bases << endl;

    }
    
      
		// create SAM faidx
		string samtools = settings.ctool("samtools");
		string command = samtools + " faidx " + settings.reference_fasta_file_name;
		SYSTEM(command);
    
		// calculate trim files
		calculate_trims(settings.reference_fasta_file_name, settings.sequence_conversion_path);
    settings.track_intermediate_file(settings.output_done_file_name, settings.sequence_conversion_path + "/*.trims");

		// store summary information
		summary.sequence_conversion.store(settings.sequence_conversion_summary_file_name);
		settings.done_step(settings.sequence_conversion_done_file_name);
	}

	summary.sequence_conversion.retrieve(settings.sequence_conversion_summary_file_name);
	ASSERT(summary.sequence_conversion.max_read_length != UNDEFINED_UINT32, "Can't retrieve max read length from file: " + settings.sequence_conversion_summary_file_name);

  //(re)load the reference sequences from our converted files
  ref_seq_info.ReadGFF(settings.reference_gff3_file_name);
  
  // update the normal versus junction-only lists
  settings.init_reference_sets(ref_seq_info);
  
  // Calculate the total reference sequence length
  summary.sequence_conversion.total_reference_sequence_length = ref_seq_info.total_length();
  
  // @JEB -- This is a bit of an ugly wart from when converting the input file was optional.
	// reload certain information into $settings from $summary  
	for (map<string, Summary::AnalyzeFastq>::iterator it = summary.sequence_conversion.reads.begin(); it != summary.sequence_conversion.reads.end(); it++)
	{
		string read_file = it->first;
		if (it->second.converted_fastq_name.size() > 0)
			settings.read_files.read_file_to_converted_fastq_file_name_map[read_file] = it->second.converted_fastq_name;
	}

	//
  // 02_reference_alignment
	// * Match all reads against the reference genome
	//

	if (
      settings.do_step(settings.reference_alignment_done_file_name, "Read alignment to reference genome")
      && !settings.aligned_sam_mode
      )
	{
		create_path(settings.reference_alignment_path);

    
    /////////////////////////
    // BUILD MAPPING INDEX //
    /////////////////////////
    
    string reference_hash_file_name = settings.reference_hash_file_name;
    string reference_fasta_file_name = settings.reference_fasta_file_name;
    string command = "bowtie2-build -q " + settings.reference_fasta_file_name + " " + reference_hash_file_name;
    SYSTEM(command);
    settings.track_intermediate_file(settings.reference_alignment_done_file_name, settings.reference_hash_file_name + "*");
    
    ///////////////////////
    // STAGE 1 ALIGNMENT //
    ///////////////////////
    
    for (uint32_t i = 0; i < settings.read_files.size(); i++)
    {
      cReadFile read_file = settings.read_files[i];
      string base_read_file_name = read_file.base_name();
      string read_fastq_file = settings.base_name_to_read_file_name(base_read_file_name);
              
      //Paths
      string stage1_reference_sam_file_name = settings.file_name(settings.stage1_reference_sam_file_name, "#", base_read_file_name);
      string stage1_unmatched_fastq_file_name = settings.file_name(settings.stage1_unmatched_fastq_file_name, "#", base_read_file_name);
      
      //Split alignment into unmatched and matched files.
      uint32_t bowtie2_seed_substring_size_stringent = trunc(summary.sequence_conversion.reads[settings.read_files[i].base_name()].avg_read_length * 0.5);
      // Check bounds
      bowtie2_seed_substring_size_stringent = max<uint32_t>(9, bowtie2_seed_substring_size_stringent);
      bowtie2_seed_substring_size_stringent = min<uint32_t>(31, bowtie2_seed_substring_size_stringent);
      
      string command = "bowtie2 -t -p " + s(settings.num_processors) + " --local " + " -L " + to_string<uint32_t>(bowtie2_seed_substring_size_stringent) + " "
      + settings.bowtie2_score_parameters + " " + settings.bowtie2_min_score_stringent + " --reorder -x " + reference_hash_file_name + " -U " + read_fastq_file + " -S " + stage1_reference_sam_file_name + " --un " + stage1_unmatched_fastq_file_name; 
      SYSTEM(command);
      
      settings.track_intermediate_file(settings.reference_alignment_done_file_name, stage1_unmatched_fastq_file_name);
      settings.track_intermediate_file(settings.reference_alignment_done_file_name, stage1_reference_sam_file_name);
    }

 
    ///////////////////////
    // STAGE 2 ALIGNMENT //
    ///////////////////////
    
		/// align reads to reference sequences
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
        
        // If we are doing staged alignment -- only align the unmatched reads with SSAHA2 and save to different name initially
        string read_fastq_file = settings.file_name(settings.stage1_unmatched_fastq_file_name, "#", base_read_file_name);
        string reference_sam_file_name = settings.file_name(settings.stage2_reference_sam_file_name, "#", base_read_file_name);
          
        uint32_t bowtie2_seed_substring_size_relaxed = 5 + trunc(summary.sequence_conversion.reads[settings.read_files[i].base_name()].avg_read_length * 0.1);
        // Check bounds
        bowtie2_seed_substring_size_relaxed = max<uint32_t>(9, bowtie2_seed_substring_size_relaxed);
        bowtie2_seed_substring_size_relaxed = min<uint32_t>(31, bowtie2_seed_substring_size_relaxed);
        
        string command = "bowtie2 -t -p " + s(settings.num_processors) + " --local " + " -L " + to_string<uint32_t>(bowtie2_seed_substring_size_relaxed) 
          + " " + settings.bowtie2_score_parameters + " " + settings.bowtie2_min_score_relaxed + " --reorder -x " + reference_hash_file_name + " -U " + read_fastq_file + " -S " + reference_sam_file_name; 
        SYSTEM(command);
        
        settings.track_intermediate_file(settings.reference_alignment_done_file_name, reference_sam_file_name);
        
        /////////////////////
        // MERGE SAM FILES //
        /////////////////////        
        
        // Merge the stage1 and stage2 output files
        {
          string base_read_file_name = read_file.base_name();
          string read_fastq_file = settings.base_name_to_read_file_name(base_read_file_name);
          
          string stage1_reference_sam_file_name = settings.file_name(settings.stage1_reference_sam_file_name, "#", base_read_file_name);
          string stage2_reference_sam_file_name = settings.file_name(settings.stage2_reference_sam_file_name, "#", base_read_file_name);
          string reference_sam_file_name = settings.file_name(settings.reference_sam_file_name, "#", base_read_file_name);
          string reference_fasta_file_name = settings.file_name(settings.reference_fasta_file_name, "#", base_read_file_name);
          
          PreprocessAlignments::merge_sort_sam_files(
                                                     stage1_reference_sam_file_name,
                                                     stage2_reference_sam_file_name,
                                                     reference_sam_file_name
                                                     );
          
          // Not deleted until after resolving alignments
          settings.track_intermediate_file(settings.alignment_correction_done_file_name, reference_sam_file_name);
        }
      }
    }

		settings.done_step(settings.reference_alignment_done_file_name);
	}

  //
	// 03_candidate_junctions
	// * Identify candidate junctions from split read alignments
	//
	if (
      !settings.no_junction_prediction
      && !settings.aligned_sam_mode
      )
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
      SYSTEM(command);
      command = samtools + " sort " + coverage_junction_best_bam_unsorted_file_name + " " + coverage_junction_best_bam_prefix;
      SYSTEM(command);
      command = samtools + " index " + coverage_junction_best_bam_file_name;
      SYSTEM(command);
      
      settings.track_intermediate_file(settings.coverage_junction_done_file_name, coverage_junction_best_bam_unsorted_file_name);
      settings.track_intermediate_file(settings.coverage_junction_done_file_name, preprocess_junction_best_sam_file_name);

      // Count errors
      string reference_fasta_file_name = settings.reference_fasta_file_name;
      string reference_bam_file_name = settings.coverage_junction_best_bam_file_name;

      error_count(
        settings,
        summary,
        reference_bam_file_name,
        reference_fasta_file_name,
        settings.candidate_junction_path,
        settings.read_file_names,
        true, // coverage
        false, // errors
        true, //preprocess
        settings.base_quality_cutoff,
        "" //covariates
      );
      settings.track_intermediate_file(settings.coverage_junction_done_file_name, reference_bam_file_name);
      settings.track_intermediate_file(settings.coverage_junction_done_file_name, reference_bam_file_name + ".bai");



      CoverageDistribution::analyze_unique_coverage_distributions(settings, 
                                                                  summary, 
                                                                  ref_seq_info,
                                                                  settings.coverage_junction_plot_file_name, 
                                                                  settings.coverage_junction_distribution_file_name,
                                                                  settings.coverage_junction_done_file_name
                                                                  );

      // Note that storing from unique_coverage and reloading in preprocess_coverage is by design
      summary.unique_coverage.store(settings.coverage_junction_summary_file_name);
      summary.preprocess_error_count.store(settings.coverage_junction_error_count_summary_file_name);
      settings.done_step(settings.coverage_junction_done_file_name);
		}
    summary.preprocess_coverage.retrieve(settings.coverage_junction_summary_file_name);
    summary.preprocess_error_count.retrieve(settings.coverage_junction_error_count_summary_file_name);
    
		string candidate_junction_summary_file_name = settings.candidate_junction_summary_file_name;
		if (settings.do_step(settings.candidate_junction_done_file_name, "Identifying junction candidates"))
		{
      CandidateJunctions::identify_candidate_junctions(settings, summary, ref_seq_info);
      
			string samtools = settings.ctool("samtools");
			string faidx_command = samtools + " faidx " + settings.candidate_junction_fasta_file_name;
			if (!file_empty(settings.candidate_junction_fasta_file_name.c_str()))
				SYSTEM(faidx_command);
      
      settings.track_intermediate_file(settings.output_done_file_name, settings.candidate_junction_fasta_file_name);
      settings.track_intermediate_file(settings.output_done_file_name, settings.candidate_junction_fasta_file_name + ".fai");



			summary.candidate_junction.store(candidate_junction_summary_file_name);
			settings.done_step(settings.candidate_junction_done_file_name);
		}
		summary.candidate_junction.retrieve(candidate_junction_summary_file_name);

    
    //
    // 04 candidate_junction_alignment
    // * Align reads to new junction candidates
    //
		if (
        settings.do_step(settings.candidate_junction_alignment_done_file_name, "Re-alignment to junction candidates")
        && !settings.aligned_sam_mode
        )
		{
			create_path(settings.candidate_junction_alignment_path);

			/// create ssaha2 hash
			string candidate_junction_hash_file_name = settings.candidate_junction_hash_file_name;
			string candidate_junction_fasta_file_name = settings.candidate_junction_fasta_file_name;

			if (!file_empty(candidate_junction_fasta_file_name.c_str()))
			{
        string command = "bowtie2-build -q " + candidate_junction_fasta_file_name + " " + candidate_junction_hash_file_name;
        SYSTEM(command);
        settings.track_intermediate_file(settings.candidate_junction_alignment_done_file_name, candidate_junction_hash_file_name + "*");
			}

			/// ssaha2 align reads to candidate junction sequences
			for (uint32_t i = 0; i < settings.read_files.size(); i++)
			{
				string base_read_file_name = settings.read_files[i].m_base_name;
				string candidate_junction_sam_file_name = settings.file_name(settings.candidate_junction_sam_file_name, "#", base_read_file_name);

				string read_fastq_file = settings.base_name_to_read_file_name(base_read_file_name);

        string filename = candidate_junction_hash_file_name + ".1.bt2";
        if (!file_exists(filename.c_str())) 
          continue;
    
        uint32_t bowtie2_seed_substring_size_stringent = trunc(summary.sequence_conversion.reads[settings.read_files[i].base_name()].avg_read_length * 0.5);
        // Check bounds
        bowtie2_seed_substring_size_stringent = max<uint32_t>(9, bowtie2_seed_substring_size_stringent);
        bowtie2_seed_substring_size_stringent = min<uint32_t>(31, bowtie2_seed_substring_size_stringent);
        
        string command = "bowtie2 -t -p " + s(settings.num_processors) + " --local " + " -L " + to_string<uint32_t>(bowtie2_seed_substring_size_stringent) + " "
         + settings.bowtie2_score_parameters + " " + settings.bowtie2_min_score_stringent + " --reorder -x " + candidate_junction_hash_file_name + " -U " + read_fastq_file + " -S " + candidate_junction_sam_file_name; 
        SYSTEM(command);
        
        settings.track_intermediate_file(settings.output_done_file_name, candidate_junction_sam_file_name + "*");
      }

			settings.done_step(settings.candidate_junction_alignment_done_file_name);
		}
  }

	
	//
  // 05 alignment_correction
	// * Resolve matches to new junction candidates
	//
	if (settings.do_step(settings.alignment_correction_done_file_name, "Resolving alignments with junction candidates"))
	{
		create_path(settings.alignment_resolution_path);

    bool junction_prediction = !settings.no_junction_prediction;
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
		if (!file_empty(resolved_junction_sam_file_name.c_str()))
		{
			command = samtools + " import " + candidate_junction_faidx_file_name + " " + resolved_junction_sam_file_name + " " + junction_bam_unsorted_file_name;
			SYSTEM(command);
			command = samtools + " sort " + junction_bam_unsorted_file_name + " " + junction_bam_prefix;
      SYSTEM(command);
			if (!settings.keep_all_intermediates)
				remove_file(junction_bam_unsorted_file_name.c_str());
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
			remove_file(reference_bam_unsorted_file_name.c_str());
		command = samtools + " index " + reference_bam_file_name;
    SYSTEM(command);

    settings.track_intermediate_file(settings.output_done_file_name, settings.junction_bam_file_name);
    settings.track_intermediate_file(settings.output_done_file_name, settings.junction_bam_file_name + ".bai");

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

		// deal with distribution or error count keys being undefined...

		uint32_t num_read_files = settings.read_files.size();
    uint32_t num_qual;
    if (!settings.aligned_sam_mode) {
      num_qual = summary.sequence_conversion.max_qual + 1; // only filled in when using FASTQ input
    } else {
      num_qual = summary.alignment_resolution.max_sam_base_quality_score + 1; // only filled in when using aligned_sam_mode
    }

		error_count(
      settings,
      summary,
			reference_bam_file_name, // bam
			reference_fasta_file_name, // fasta
			settings.error_calibration_path, // output
			settings.read_files.base_names(), // readfile
			true, // coverage
			true, // errors
      false, //preprocess
			settings.base_quality_cutoff, // minimum quality score
			"read_set=" + to_string(num_read_files) + ",obs_base,ref_base,quality=" + to_string(num_qual) // covariates
		);

		settings.done_step(settings.error_counts_done_file_name);
	}


	//
	// Calculate error rates
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
                                                                  settings.unique_only_coverage_distribution_file_name,
                                                                  "NULL" // means to never delete intermediates
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
			command = "R --vanilla in_file=" + cString(error_rates_base_qual_error_prob_file_name).escape_shell_chars() +
        " out_file=" + cString(error_rates_plot_file_name).escape_shell_chars() +
        " < "        + cString(plot_error_rates_r_script_file_name).escape_shell_chars() +
        " > "        + cString(plot_error_rates_r_script_log_file_name).escape_shell_chars();
			SYSTEM(command,false, false, false); //NOTE: Not escaping shell characters here.
		}

		summary.unique_coverage.store(settings.error_rates_summary_file_name);
		settings.done_step(settings.error_rates_done_file_name);
	}
	summary.unique_coverage.retrieve(settings.error_rates_summary_file_name);

  //
	// 08 Mutation Identification
	// Make predictions of point mutations, small indels, and large deletions
	//

	if (!settings.no_mutation_prediction)
	{
		create_path(settings.mutation_identification_path);

		if (settings.do_step(settings.mutation_identification_done_file_name, "Examining read alignment evidence"))
		{
			string reference_fasta_file_name = settings.reference_fasta_file_name;
			string reference_bam_file_name = settings.reference_bam_file_name;

			string output_dir = settings.mutation_identification_path;
			string ra_mc_genome_diff_file_name = settings.ra_mc_genome_diff_file_name;

			// It is important that these are in consistent order with the fasta file!!
      vector<double> deletion_propagation_cutoffs;
      vector<double> deletion_seed_cutoffs;
      for (uint32_t i = 0; i < ref_seq_info.size(); i++) {
        deletion_propagation_cutoffs.push_back(summary.unique_coverage[ref_seq_info[i].m_seq_id].deletion_coverage_propagation_cutoff);
        deletion_seed_cutoffs.push_back(summary.unique_coverage[ref_seq_info[i].m_seq_id].deletion_coverage_seed_cutoff);
        //cout << ref_seq_info[i].m_seq_id << " " << to_string(deletion_propagation_cutoffs.back()) << " " << to_string(deletion_seed_cutoffs.back()) << endl;
      }

			identify_mutations(
        settings,
        summary,
				reference_bam_file_name,
				reference_fasta_file_name,
				ra_mc_genome_diff_file_name,
        deletion_propagation_cutoffs,
        deletion_seed_cutoffs,
				settings.mutation_log10_e_value_cutoff, // mutation_cutoff
				settings.polymorphism_log10_e_value_cutoff, // polymorphism_cutoff
				settings.polymorphism_frequency_cutoff, //polymorphism_frequency_cutoff
				settings.print_mutation_identification_per_position_file //per_position_file
			);

			settings.done_step(settings.mutation_identification_done_file_name);
		}
    // extra processing for polymorphisms / mixed-base prediction
		if ( (settings.polymorphism_prediction || settings.mixed_base_prediction) && settings.do_step(settings.polymorphism_statistics_done_file_name, "Polymorphism statistics"))
		{
			ref_seq_info.polymorphism_statistics(settings, summary);
			settings.done_step(settings.polymorphism_statistics_done_file_name);
		}
	}

	//rewire which GenomeDiff we get data from if we have the elaborated polymorphism_statistics version
  if (settings.polymorphism_prediction || settings.mixed_base_prediction)
    settings.ra_mc_genome_diff_file_name = settings.polymorphism_statistics_ra_mc_genome_diff_file_name;

    /*
     * 09 Copy number variation
     * --------------------------------
     * This conditional is run if "--cnv" is passed as argument when invoking breseq on command line.
     * Ex: breseq --cnv -r lambda.gbk lambda.fastq
     * NOTE: "--cnv" runs the whole pipeline including the below copy number variation conditional,
     * whereas "CNV" just runs the copy number variation part of the pipeline by itself, see separate
     * do_copy_number_variation() function.
     */
  if (settings.do_copy_number_variation) {
        
    // Create copy_number_variation directory
    create_path( settings.copy_number_variation_path );

    if (settings.do_step(settings.copy_number_variation_done_file_name, "Predicting copy number variation")) { 

    for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it) {

      // Sequence iterator, one sequence at a time from .fastq file
      cAnnotatedSequence& seq = *it;

      // Create filename: [genome].coverage.txt, this is LONG file where the coverage data is coming from
      string this_complete_coverage_text_file_name = settings.file_name(settings.complete_coverage_text_file_name, "@", seq.m_seq_id);

      // Create filename: [genome].tiled.tab, this is LONG file where the tiled coverage data will be output
      string this_tiled_complete_coverage_text_file_name = settings.file_name(settings.tiled_complete_coverage_text_file_name, "@", seq.m_seq_id);

      // Create filename: [genome].tiled_for_edging.tab, this is where the tiled
      // coverage data at a smaller tile_size will be output for use in improved
      // edge detection (aka "edging")
      string this_tiled_for_edging_text_file_name = settings.file_name(settings.tiled_for_edging_text_file_name, "@", seq.m_seq_id);

      // Generates [genome].tiled.tab, one line at a time with each for-loop,
      // where the avg coverage of each tile is calculated.
      CoverageDistribution::tile(settings.ignore_redundant_coverage,
                                 this_complete_coverage_text_file_name, // in_file_name
                                 this_tiled_complete_coverage_text_file_name, // out_file_name
                                 this_tiled_for_edging_text_file_name,
                                 settings.copy_number_variation_tile_size);

      // Create filename: [genome].ranges.tab, this is SHORT file used for ???
      // (contains: Start_Position, End_Position, T_Score, P_Value)
      string this_ranges_text_file_name = settings.file_name(settings.ranges_text_file_name, "@", seq.m_seq_id);

      // Get filename: [genome].history.tab, this is SHORT file used for ???
      // (contains: Start_Search, End_Search, Start_Position, End_Position, Start_Segment, End_Segment... 16 values total)
      string this_cnv_history_text_file_name = settings.file_name(settings.cnv_history_text_file_name, "@", seq.m_seq_id);

      // Generates [genome].ranges.tab & [genome].history.tab, one line at a time with each for-loop,
      CoverageDistribution::find_segments(settings,
                                          summary.unique_coverage[seq.m_seq_id].nbinom_mean_parameter,
                                          this_tiled_complete_coverage_text_file_name,
                                          this_tiled_for_edging_text_file_name, // tiled_for_edging_file_name (tiled_for_edging.tab)
                                          this_ranges_text_file_name,
                                          this_cnv_history_text_file_name
                                          );

      // Create filename: [genome].smoothed_ranges.tab, this is LONG file used for ???
      // (contains: Position, Smooth_Coverage)
      string this_smoothed_ranges_text_file_name = settings.file_name(settings.smoothed_ranges_text_file_name, "@", seq.m_seq_id);

      // Create filename: [genome].cnv_final.tab, this is SHORT file used for ???
      // (contains: Start_Position, End_Position, Z_Score, Greater_Than, Copy_Number)
      string this_final_cnv_file_name = settings.file_name(settings.final_cnv_text_file_name, "@", seq.m_seq_id);

      // Create filename: [genome].cn_evidence.gd, this is SHORT file used for ???
      // (contains no labels)
      string this_copy_number_variation_cn_genome_diff_file_name = settings.file_name(settings.copy_number_variation_cn_genome_diff_file_name, "@", seq.m_seq_id);

      // Generates [genome].smoothed_ranges.tab & [genome].cnv_final.tab & [genome].cn_evidence.gd,
      // one line at a time with each for-loop,
      CoverageDistribution::smooth_segments(settings,
                                            seq.m_seq_id,
                                            summary.unique_coverage[seq.m_seq_id].nbinom_mean_parameter, 
                                            this_tiled_for_edging_text_file_name, 
                                            this_ranges_text_file_name, 
                                            this_smoothed_ranges_text_file_name,
                                            this_final_cnv_file_name,
                                            this_copy_number_variation_cn_genome_diff_file_name
                                            );

      } // End of foreach reference loop
      
    } // End of if cnv done file
    
    if (settings.do_periodicity) {
      create_path (settings.copy_number_variation_path);
      
      if (settings.do_step(settings.periodicity_done_file_name, "Periodicity")){
        for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it){
          cAnnotatedSequence& seq = *it;
          
          string this_complete_coverage_text_file_name = settings.file_name(settings.complete_coverage_text_file_name, "@", seq.m_seq_id);
          string this_period_complete_coverage_text_file_name = settings.file_name(settings.periodicity_table_file_name, "@", seq.m_seq_id);
          
          CoverageDistribution::calculate_periodicity(
                this_complete_coverage_text_file_name,
                this_period_complete_coverage_text_file_name,
                settings.periodicity_method,
                settings.periodicity_start,
                settings.periodicity_end,
                settings.periodicity_step
                );
        }
        settings.done_step(settings.periodicity_done_file_name);
      }
      
    } // End of if_do_periodicity
  
  } // End of if do_cnv

   
  create_path(settings.evidence_path); //need output for plots

	if (settings.do_step(settings.output_done_file_name, "Output"))
	{
		///
		// Output Genome Diff File
		///
		cerr << "Creating merged genome diff evidence file..." << endl;

		// merge all of the evidence GenomeDiff files into one...
		create_path(settings.evidence_path);
    
    cGenomeDiff jc_gd(settings.jc_genome_diff_file_name);
         
    // Add read count information to the JC entries -- call fails if no predictions, because BAM does not exist
    MutationPredictor mpj(ref_seq_info);
    mpj.prepare_junctions(settings, summary, jc_gd); // this step has to be done twice currently, which is a bit wasteful
    assign_junction_read_counts(settings, summary, jc_gd);
    //assign_junction_read_coverage(settings, summary, jc_gd);

    cGenomeDiff ra_mc_gd(settings.ra_mc_genome_diff_file_name);
    
    cGenomeDiff evidence_gd;
    evidence_gd.fast_merge(jc_gd);
    evidence_gd.fast_merge(ra_mc_gd);
    
    // there is a copy number genome diff for each sequence separately
    if (settings.do_copy_number_variation) {
      for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it) {
        cAnnotatedSequence& seq = *it;
        string this_copy_number_variation_cn_genome_diff_file_name = settings.file_name(settings.copy_number_variation_cn_genome_diff_file_name, "@", seq.m_seq_id);
        cGenomeDiff cn_gd(this_copy_number_variation_cn_genome_diff_file_name);
        evidence_gd.fast_merge(cn_gd);
      }
    }
		evidence_gd.write(settings.evidence_genome_diff_file_name);

		// predict mutations from evidence in the GenomeDiff
		cerr << "Predicting mutations from evidence..." << endl;

    MutationPredictor mp(ref_seq_info);
    cGenomeDiff mpgd(settings.evidence_genome_diff_file_name);
    mp.predict(settings, summary, mpgd);
    mpgd.reassign_unique_ids();

    //#=REFSEQ header lines.
    mpgd.metadata.ref_seqs = settings.all_reference_file_names;
    
    //#=READSEQ header lines.
    mpgd.metadata.read_seqs.resize(settings.read_files.size());
    for (size_t i = 0; i < settings.read_files.size(); i++) {
      mpgd.metadata.read_seqs[i] = settings.read_files[i].file_name();
    }
    
    mpgd.metadata.author = string(PACKAGE_STRING) + " " + string(HG_REVISION);
    mpgd.metadata.breseq_data["CREATED"] = Settings::time2string(time(NULL));

    /* @JEB: Could add information such as average coverage, etc. as metadate here for
       generating more detailed summary files at later steps.
    // Add additional header lines if needed.
    if (settings.add_metadata_to_gd){
      for (storable_map<string, Summary::Coverage>::iterator it = summary.unique_coverage.begin();
            it != summary.unique_coverage.end(); it ++) {
         //Usually needed for gathering breseq data.
       }
    }
    */

    // Write and reload 
    mpgd.write(settings.final_genome_diff_file_name);
    cGenomeDiff gd(settings.final_genome_diff_file_name);
    // Empty metadata... necessary to keep consistency tests
    // when things like time are being added to output.gd
    gd.metadata = cGenomeDiff::Metadata();
    
    // Write VCF conversion
    mpgd.write_vcf(settings.output_vcf_file_name, ref_seq_info);
    
    //
    // Mark marginal items as no_show to prevent further processing
    //
    output::mark_gd_entries_no_show(settings, gd);

    //
		// Annotate mutations
		//
		cerr << "Annotating mutations..." << endl;
		ref_seq_info.annotate_mutations(gd);
    
    gd.write(settings.annotated_genome_diff_file_name);
    
		//
		// Plot coverage of genome and large deletions
		//
		cerr << "Drawing coverage plots..." << endl;
		output::draw_coverage(settings, ref_seq_info, gd);
    
		//
		// Create evidence files containing alignments and coverage plots
		//
		if (!settings.no_alignment_or_plot_generation)
			output::cOutputEvidenceFiles(settings, gd);

		///
		// HTML output
		///

		cerr << "Creating index HTML table..." << endl;

		output::html_index(settings.index_html_file_name, settings, summary, ref_seq_info, gd);
		output::html_marginal_predictions(settings.marginal_html_file_name, settings, summary, ref_seq_info, gd);
        
		// record the final time and print summary table
		settings.record_end_time("Output");

		output::html_summary(settings.summary_html_file_name, settings, summary, ref_seq_info);

		settings.done_step(settings.output_done_file_name);
	}
  cerr << "+++   SUCCESSFULLY COMPLETED" << endl;

  
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
  
  // Sequence Utility Commands:
  if (command == "CONVERT-FASTQ") {
		return do_convert_fastq(argc_new, argv_new);
	} else if (command == "CONVERT-REFERENCE") {
		return do_convert_genbank(argc_new, argv_new);
  } else if ((command == "GET-SEQUENCE") || (command == "SUBSEQUENCE")) {
    return do_get_sequence(argc_new, argv_new);
    
  // Breseq Post-Run Commands:
  } else if (command == "BAM2ALN") {
    return do_bam2aln( argc_new, argv_new);  
  } else if (command == "BAM2COV") {
    return do_bam2cov( argc_new, argv_new);
    
    
  // Experimental and Development Commands:
  //
  // None of these commands are documented for use by others. 
  // They may change without warning.
  } else if ((command == "SIMULATE-READ") ||  (command == "SIMULATE-READS")) {
    return do_simulate_read(argc_new, argv_new); 
  } else if (command == "REFALN") {
    return do_ref_aln();  
  } else if (command == "CNV") {
    return do_copy_number_variation(argc_new, argv_new);
  } else if (command == "PERIODICITY"){
    return do_periodicity(argc_new, argv_new);
  } else if (command == "ERROR_COUNT") {
    return do_error_count(argc_new, argv_new);
  } else if (command == "IDENTIFY_CANDIDATE_JUNCTIONS") {
    return do_identify_candidate_junctions(argc_new, argv_new);
  } else if (command == "JUNCTION-POLYMORPHISM") {
    return do_junction_polymorphism(argc_new, argv_new);
  } else if (command == "CL_TABULATE") {
    return do_tabulate_contingency_loci(argc_new, argv_new);
  } else if (command == "CL_SIGNIFICANCE") {
    return do_analyze_contingency_loci_significance( argc_new, argv_new);
  } else if (command == "ASSEMBLE-UNMATCHED") {
    return do_assemble_unmatched( argc_new, argv_new);
  }
  else {
    // Not a sub-command. Use original argument list.
    return breseq_default_action(argc, argv);
  }
	return -1;
}
