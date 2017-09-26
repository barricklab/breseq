 /*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2017 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#include <iostream>
#include <string>
#include <vector>

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
#include "libbreseq/samtools_commands.h"
#include "libbreseq/settings.h"
#include "libbreseq/summary.h"
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
	AnyOption options("Usage: breseq BAM2ALN [-b reference.bam -f reference.fasta -o alignment.html -n 200] region1 [region2 region3 ...]");
  options.addUsage("");
  options.addUsage("Display reads aligned to the specified region or regions.");
  options.addUsage("");
  options.addUsage("Allowed Options");
  options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("bam,b", "BAM database file of read alignments", "data/reference.bam");
  options("fasta,f", "FASTA file of reference sequences", "data/reference.fasta");
  options("output,o", "Output path. If there is just one region, the name of the output file (DEFAULT=region1.*). If there are multiple regions, this argument must be a directory path, and all output files will be output here with names region1.*, region2.*, ... (DEFAULT=.)");
  options("region,r", "Regions to create alignments for. Must be provided as sequence regions in the format ACCESSION:START-END, where ACCESSION is a valid identifier for one of the sequences in the FASTA file, and START and END are 1-indexed coordinates of the beginning and end positions. Any read overlapping these positions will be shown. A separate output file is created for each region. Regions may be provided at the end of the command line as unnamed arguments");
  options("format", "Format of output alignment(s): HTML or TXT", "HTML");
  options("max-reads,n", "Maximum number of reads to show in alignment", 200);
  options("repeat", "Show reads with multiple best matches in reference", TAKES_NO_ARGUMENT, ADVANCED_OPTION);
  options("quality-score-cutoff,c", "Quality score cutoff below which reads are highlighted as yellow", 0);
  options("stdout", "Write output to stdout", TAKES_NO_ARGUMENT, ADVANCED_OPTION);
  options.processCommandArgs(argc, argv);
  
	// make sure that the config options are good:
  if(options.count("help")) {
		options.printAdvancedUsage();
		return -1;
	}
  
  if (!file_exists(options["fasta"].c_str())) {
    options.addUsage("");
    options.addUsage("Could not open input reference FASTA file (-f):\n  " + options["fasta"]);
    options.printUsage();
    return -1;
  }
                     
  if (!file_exists(options["bam"].c_str())) {
    options.addUsage("");
    options.addUsage("Could not open input BAM file of aligned reads (-b):\n  " + options["bam"]);
    options.printUsage();
    return -1;
  }
  
  vector<string> region_list;
  if (options.count("region")) {
    region_list= from_string<vector<string> >(options["region"]);
  }
  
  // Also take regions off the command line
  for (int32_t i = 0; i < options.getArgc(); i++) {
    string region = options.getArgv(i);
    region_list.push_back(region);
  }
  
  if (region_list.size() == 0) {
    options.addUsage("");
    options.addUsage("No input regions provided (-r or unnamed args).");
    options.printUsage();
    return -1;
  }
  
  string format = to_upper(options["format"]);
  if ((format != "HTML") && (format != "TXT")) {
    options.addUsage("");
    options.addUsage("Unknown format requested: " + format);
    options.printUsage();
    return -1;
  }
  
  cerr << "COMMAND: BAM2ALN" << endl;
  cerr << "+++   Creating alignments..." << endl;

  for(uint32_t j = 0; j < region_list.size(); j++) {
    
    cReferenceSequences::normalize_region(region_list[j]);
    cerr << "  Region : " << region_list[j] << endl;
    
    
    // Generate Alignment!
    alignment_output ao(
                        options["bam"],
                        options["fasta"],
                        from_string<uint32_t>(options["max-reads"]),
                        from_string<uint32_t>(options["quality-score-cutoff"]),
                        1,
                        false,
                        options.count("repeat")
                        );
    
    string default_file_name = region_list[j];
    
    string output_string;
    if (format == "HTML") {
      output_string = ao.html_alignment(region_list[j]);
      default_file_name += ".html";
    }
    else if (format == "TXT") {
      output_string = ao.text_alignment(region_list[j]);
      default_file_name += ".txt";
    }
      
    if (options.count("stdout")) {
      cout << output_string << endl;
    } else {
      ///Write to file
      string file_name = default_file_name;
            
      if (options.count("output")) {
        file_name = options["output"];
        if(region_list.size() > 1)  {
          file_name = (split(options["output"], "."))[0] + "_" + region_list[j] + ".html";  }
      }
      
      ofstream myfile (file_name.c_str());
      cerr << "    File : " << file_name << endl;
      
      if (myfile.is_open()) {
        Settings settings;
        if (to_upper(options["format"]) == "HTML")
          myfile << html_header("BRESEQ :: bam2aln output", settings);
        myfile << output_string;
        if (to_upper(options["format"]) == "HTML")
          myfile << html_footer();
        myfile.close();
      } else {
        cerr << "Unable to open file: " << file_name.c_str();
      }
    }
  }
  
  cerr << "+++   SUCCESSFULLY COMPLETED" << endl;
  return 0;
}

/*! bam2cov
 Draw HTML coverage from BAM
 */
int do_bam2cov(int argc, char* argv[]) {
  // setup and parse configuration options:
	AnyOption options("Usage: breseq BAM2COV [-b reference.bam -f reference.fasta --format PNG -o output.png] region1 [region2 region3 ...]");
  options.addUsage("");
  options.addUsage("Create a coverage plot or table for the specified region or regions.");
  options.addUsage("");
  options.addUsage("Allowed Options");

  options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  // required options
  options("bam,b", "BAM database file of read alignments", "data/reference.bam");
  options("fasta,f", "FASTA file of reference sequences", "data/reference.fasta");
  // options controlling what files are output
  options("output,o", "Output path. If there is just one region, the name of the output file (DEFAULT=region1.*). If there are multiple regions, this argument must be a directory path, and all output files will be output here with names region1.*, region2.*, ... (DEFAULT=.)");
  options("region,r", "Regions to create alignments for. Must be provided as sequence regions in the format ACCESSION:START-END, where ACCESSION is a valid identifier for one of the sequences in the FASTA file, and START and END are 1-indexed coordinates of the beginning and end positions. Any read overlapping these positions will be shown. A separate output file is created for each region. Regions may be provided at the end of the command line as unnamed arguments");
  options("format", "Format of output plot(s): PNG or PDF", "PNG");
  options("table,t", "Create tab-delimited file of coverage instead of a plot", TAKES_NO_ARGUMENT);
  options.addUsage("", ADVANCED_OPTION);
  options.addUsage("Advanced Output Options", ADVANCED_OPTION);
  options("total-only,1", "Only plot/tabulate the total coverage at a position. That is, do not not output the coverage on each genomic strand", TAKES_NO_ARGUMENT, ADVANCED_OPTION);
  options("resolution,p", "Number of positions to output coverage information for in interval (0=ALL)", 600, ADVANCED_OPTION);
  options("show-average,a", "Show the average coverage across the reference sequence as a horizontal line. Only possible if used in the main output directory of breseq output", TAKES_NO_ARGUMENT, ADVANCED_OPTION);
  options("fixed-coverage-scale,s", "Fix the maximum value on plots the coverage scale in plots. If the show-average option is provided, then this is a factor that will be multiplied times the average coverage (e.g., 1.5 x avg). Otherwise, this is a coverage value (e.g., 100-fold coverage)", "", ADVANCED_OPTION);
  
  options.addUsage("", ADVANCED_OPTION);
  options.addUsage("Tiling Mode (produce plots that span reference sequences from end to end)", ADVANCED_OPTION);
  options
  ("tile-size", "In tiling mode, the size of each tile that will be output as a separate file", "", ADVANCED_OPTION)
  ("tile-overlap", "In tiling mode, overlap between adjacent tiles (1/2 of this is added to each side of every tile)", "", ADVANCED_OPTION)
  ;
//  ("read_start_output,r", "file name for table file binned by read start bases (DEFAULT: OFF)")
//  ("gc_output,g", "create additional table file binned by GC content of reads (DEFAULT: OFF)")

  options.processCommandArgs(argc, argv);

  if (options.count("help"))
  {
    options.printAdvancedUsage();
    exit(-1);
  }
  
  if (!file_exists(options["fasta"].c_str())) {
    options.addUsage("");
    options.addUsage("Could not open input reference FASTA file (-f):\n  " + options["fasta"]);
    options.printUsage();
    return -1;
  }
  
  if (!file_exists(options["bam"].c_str())) {
    options.addUsage("");
    options.addUsage("Could not open input BAM file of aligned reads (-b):\n  " + options["bam"]);
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
  
  bool tiling_mode = options.count("tile-size") || options.count("tile-overlap");
  ASSERT(tiling_mode || (!options.count("tile-size") && !options.count("tile-overlap")),
         "--tile-size and --tile-overlap args must both be provided to activate tile mode");
  
  if (!tiling_mode && !region_list.size()) {
    options.addUsage("");
    options.addUsage("You must supply at least one genomic region (-r).");
    options.printUsage();
    return -1;
  }
  
  if (tiling_mode && (region_list.size() > 0))
  {
    options.addUsage("");
    options.addUsage("You cannot both provide specific regions (-r) and use tiling mode.");
    options.printUsage();
    return -1;
  }
  
  cerr << "COMMAND: BAM2COV" << endl;
  
  if (options.count("table")) {
    cerr << "+++   Tabulating coverage..." << endl;
  } else {
    cerr << "Plotting coverage..." << endl;
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
  if (options.count("fixed-coverage-scale")) {
    co.fixed_coverage_scale(from_string<double>(options["fixed-coverage-scale"]));
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
    string& region = *it;
// these are experimental... additional table files
//    if (options.count("read_start_output"))
//      co.read_begin_output_file_name(options["read_start_output"]);
//    if (options.count("gc_output"))
//      co.gc_output_file_name(options["gc_output"]);
    
    cReferenceSequences::normalize_region(region);
    cerr << "  Region : " << region << endl;

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
        file_name = options["output"] + "/" + region;
      }
    }
    
    if (options.count("table")) {
      file_name += ".tab";
      cerr << "    File : " << file_name << endl;
      co.table(region, file_name, from_string<uint32_t>(options["resolution"]));
    } else {
      file_name += "." + to_lower(options["format"]);
      cerr << "    File : " << file_name << endl;
      co.plot(region, file_name , from_string<uint32_t>(options["resolution"]));
    }
  }
  
  cerr << "+++   SUCCESSFULLY COMPLETED" << endl;
  return 0;
}

int do_convert_fastq(int argc, char* argv[])
{
  
	// setup and parse configuration options:
	AnyOption options("Usage: breseq CONVERT-FASTQ [-o output.fastq -1 ILLUMINA_1.3+ -2 SANGER -r] input.fastq");
  options.addUsage("");
  options.addUsage("Convert between FASTQ formats using different base pair quality encodings.");
  options.addUsage("");
  options.addUsage("Valid input/output formats are 'SANGER', 'SOLEXA', 'ILLUMINA_1.3+'. GUESS will attempt to determine the input format by analyzing read quality scores. See http://wikipedia.org/wiki/FASTQ_format for a description of FASTQ formats.");
  options.addUsage("");
  options.addUsage("Allowed Options");
  options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("output,o", "output FASTQ file", "output.fastq");
  options("input-format,1", "format to convert from", "GUESS");
  options("output-format,2", "format to convert to", "SANGER");
  options("reverse-complement,r", "reverse complement all reads and add _RC to their names", TAKES_NO_ARGUMENT);
  
  options.processCommandArgs(argc, argv);
  
	// handle help
  if (options.count("help")) {
		options.printAdvancedUsage();
		return -1;
	}
     
  if (options.getArgc() == 0) {
    options.printUsage();
    return -1;
  }
  
  // make sure that the config options are good:
  if (options.getArgc() != 1) {
    options.addUsage("");
    options.addUsage("Please provide exactly one input FASTQ file.");
    options.printUsage();
    return -1;
  }
  
  string input_file_name = options.getArgv(0);
  string output_file_name = options["output"];
  string input_format = options["input-format"];
  string output_format = options["output-format"];

  cerr << "COMMAND: CONVERT-FASTQ" << endl;
  cerr << " Input file    : " << input_file_name << endl;
  cerr << " Input format  : " << input_format << endl;
  cerr << " Output file   : " << input_file_name << endl;
  cerr << " Output format : " << input_file_name << endl;

  
  bool guessed_format = false;
  if (input_format == "GUESS") {
    cerr << "+++   Predicting input format..." << endl;

    uint64_t original_num_reads;
    uint64_t original_num_bases;
    uint32_t min_read_length;
    uint32_t max_read_length;
    uint8_t min_quality_score;
    uint8_t max_quality_score;
    
    input_format = cFastqQualityConverter::predict_fastq_file_format(input_file_name, original_num_reads, original_num_bases, min_read_length, max_read_length, min_quality_score, max_quality_score);
    
    cerr << " Predicted input format : " << input_format << endl;;
  }

  cerr << "+++   Converting FASTQ..." << endl;
  convert_fastq(input_file_name, output_file_name, input_format, output_format, options.count("reverse-complement"));

  cerr << "+++   SUCCESSFULLY COMPLETED" << endl;
  return 0;
}

int do_convert_reference(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq CONVERT-REFERENCE [-f FASTA -o output.fna] input.gbk");
  options.addUsage("");
  options.addUsage("Convert a reference genome into another format. If multiple sequences");
  options.addUsage("are present in the input files, they will be merged into one output file.");
  options.addUsage("");
  options.addUsage("Allowed Options");
  options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("format,f", "Output format. Valid options: FASTA, GFF3, CSV (Default = FASTA)", "FASTA");
  options("no-sequence,n", "Do not include the nucleotide sequence. The output file will only have features. (Not allowed with FASTA format.)", TAKES_NO_ARGUMENT);
  options("output,o", "Output reference file path (Default = output.*)");

	options.processCommandArgs(argc, argv);
	
  if (options.count("help")) {
    options.printUsage();
    return -1;
  }
  
	// make sure that the config options are good:
  if (options.getArgc() == 0) {
    options.addUsage("");
    options.addUsage("No input reference file(s) provided (-r).");
		options.printUsage();
    return -1;
	}
  
  string output_format = to_upper(options["format"]);
  if (output_format=="GFF") output_format="GFF3";
  if ((output_format != "FASTA") && (output_format != "GFF3") && (output_format != "CSV")) {
    options.addUsage("");
    options.addUsage("Unknown output file format requested: " + to_upper(options["format"]));
    options.printUsage();
    return -1;
  }
  
  if ( (output_format=="FASTA") && (options.count("no-sequence")) ) {
    options.addUsage("");
    options.addUsage("The --no-sequence|n option cannot be used with output format FASTA.");
    options.printUsage();
    return -1;
  }
  
  cerr << "COMMAND: CONVERT-REFERENCE" << endl;
  
  cerr << "+++   Loading reference files..." << endl;
  vector<string> reference_file_names;
  for (int32_t i = 0; i < options.getArgc(); ++i) {
    cout << "  Input : " << options.getArgv(i) << endl;
    reference_file_names.push_back(options.getArgv(i));
  }
  
  cReferenceSequences refs;
  refs.LoadFiles(reference_file_names);

  cerr << "+++   Writing reference file..." << endl;
  
  cerr << "  Output : " << options["output"] << endl;
  cerr << "  Format : " << output_format << endl;
  if (output_format == "FASTA") {
    refs.WriteFASTA(options.count("output") ? options["output"] : "output.fna");
  } else if (output_format == "GFF3") {
    refs.WriteGFF(options.count("output") ? options["output"] : "output.gff", options.count("no-sequence"));
  } else if (output_format == "CSV") {
    refs.WriteCSV(options.count("output") ? options["output"] : "output.csv");
  }
	
  cerr << "+++   SUCCESSFULLY COMPLETED" << endl;
	return 0;
}

int do_get_sequence(int argc, char *argv[])
{
  AnyOption options("Usage: breseq GET-SEQUENCE [-o output.fna -c -r input.gbk] REL606:50-100 [REL606:378-403 ...]");
  options.addUsage("");
  options.addUsage("Takes a reference sequence and outputs the DNA base sequences of specific regions in FASTA format.");
  options.addUsage("");
  options.addUsage("Regions are of the format ACCESSION:START-END. For example, REL606:1234-2345");
  options.addUsage("");
  options.addUsage("If START is greater than END, then the reverse complement sequence will be returned.");
  options.addUsage("");
  options.addUsage("Allowed Options");
  options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("reference,r",  "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (Default=data/reference.fasta)");
  options("output,o","output FASTA file. Will write to STDOUT if not provided.");
  options("reverse-complement,c","reverse complement", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  // Handle help
  if (options.count("help")) {
    options.printAdvancedUsage();
    return -1;
  }
  
  if (options.getArgc() == 0) {
    options.addUsage("");
    options.addUsage("Must provide at least one region to retrieve sequence for.");
    options.printUsage();
    return -1;
  }
  
  cerr << "COMMAND: GET-SEQUENCE" << endl;
  
  // Sets default
  vector<string> reference_file_names;
  if (options.count("reference") == 0) {
    reference_file_names.push_back("data/reference.fasta");
  } else {
    reference_file_names = from_string<vector<string> >(options["reference"]);
  }
  
  if (reference_file_names.size() == 0) {
    options.addUsage("");
    options.addUsage("Must provide reference sequence(s) (-r).");
    options.printUsage();
    return -1;
  }
  
  vector<string> region_list;
  
  bool do_reverse_complement = options.count("reverse-complement");
  bool to_stdout = !options.count("output");
  
  // Also take positions off the command line
  for (int32_t i = 0; i < options.getArgc(); i++)
  {
    string position = options.getArgv(i);
    region_list.push_back(position);
  }
  
  cerr << "+++   Loading reference sequence(s)..." << endl;

  cReferenceSequences ref_seq_info, new_seq_info;
  ref_seq_info.LoadFiles(reference_file_names);
  
  for(uint32_t j = 0; j < region_list.size(); j++)
  {
    // clean commas
    region_list[j] = substitute(region_list[j], ",", "");
    
    cerr << "+++   Retrieving sequence of region..." << endl;
    
    uint32_t replace_target_id, replace_start, replace_end;
    string seq_name = "";
    
    do_reverse_complement = ref_seq_info.normalize_region(region_list[j]);
    if (options.count("reverse-complement")) {
      do_reverse_complement = !do_reverse_complement;
    }
    ref_seq_info.parse_region(region_list[j], replace_target_id, replace_start, replace_end);
    
    cerr << "  ACCESSION : " << ref_seq_info[replace_target_id].m_seq_id << endl;
    cerr << "  START     : " << replace_start << endl;
    cerr << "  END       : " << replace_end << endl;
    cerr << "  STRAND    : " << (do_reverse_complement ? "Reverse (Bottom)" : "Forward (Top)") << endl;
    cerr << endl;
    
    ASSERT((uint32_t)ref_seq_info[replace_target_id].m_length >= replace_start && (uint32_t)ref_seq_info[replace_target_id].m_length >= replace_end,
           "START:\t" + to_string(replace_start) + "\n" +
           "END:\t" + to_string(replace_end) + "\n" +
           "SIZE:\t" + to_string(ref_seq_info[replace_target_id].m_length) + "\n" +
           "Neither START or END can be greater than the SIZE of " + ref_seq_info[replace_target_id].m_seq_id + ".");
    
    seq_name = ref_seq_info[replace_target_id].m_seq_id + ":" + to_string(replace_start) + "-" + to_string(replace_end);
    
    new_seq_info.add_new_seq(seq_name, "");
    cAnnotatedSequence& new_seq = new_seq_info[seq_name];
    new_seq.m_fasta_sequence = ref_seq_info[replace_target_id].m_fasta_sequence;
    new_seq.m_fasta_sequence.set_name(seq_name);
    new_seq.m_fasta_sequence.set_sequence(to_upper(ref_seq_info.get_sequence_1(replace_target_id, replace_start, replace_end)));
    if(do_reverse_complement)
      new_seq.m_fasta_sequence.set_sequence(reverse_complement(new_seq.m_fasta_sequence.get_sequence()));
    new_seq.m_seq_id = region_list[j];
    new_seq.m_length = new_seq.m_fasta_sequence.get_sequence_length();
    
    if(to_stdout)  {
      cout << new_seq.m_fasta_sequence << endl;
    }
  }
  
  if(options.count("output"))  {
    cerr << "+++   Writing sequence(s) to file..." << endl;
    cerr << "  Output : " << options["output"] << endl;

    new_seq_info.WriteFASTA(options["output"]);
  }
  
  cerr << "+++   SUCCESSFULLY COMPLETED" << endl;
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
		("bam,b", "BAM file containing aligned read sequences", "data/reference.bam")
		("fasta,f", "FASTA file of reference sequence", "data/reference.fasta")
		("output,o", "output directory", "./")
		("coverage", "generate unique coverage distribution output", TAKES_NO_ARGUMENT)
		("errors", "generate unique error count output", TAKES_NO_ARGUMENT)
    ("covariates", "covariates for error model. a comma separated list (no spaces) of these choices: ref_base, obs_base, prev_base, quality, read_set, ref_pos, read_pos, base_repeat. For quality, read_pos, and base_repeat you must specify the maximum value possible, e.g. quality=40")
    ("minimum-quality-score", "ignore base quality scores lower than this", 0)
	.processCommandArgs(argc, argv);
  
	// make sure that the config options are good:
	if(options.count("help")
		 || (!options.count("coverage") && !options.count("errors")) ) {
		options.printUsage();
		return -1;
	}
  
  if(options.count("errors") && !options.count("covariates") ) {
    WARN("Must provide --covariates when --errors specified.");
		options.printUsage();
		return -1;
	}

  Summary summary;
  Settings settings;
  
  // This loading required to set up reference
  // sequences that get parsed during error count
  // BETTER: would be to read them from the saved settings
  // of the breseq run?
  cReferenceSequences ref_seq_info;
  vector<string> reference_file_names;
  reference_file_names.push_back(options["fasta"]);
  ref_seq_info.LoadFiles(reference_file_names);
  settings.init_reference_sets(ref_seq_info);
  
  vector<string> no_read_file_names;

  error_count(
              settings,
              summary,
              options["bam"],
              options["fasta"],
              options["output"],
              no_read_file_names,
              options.count("coverage"),
              options.count("errors"),
              false,
              from_string<uint32_t>(options["minimum-quality-score"]),
              options["covariates"]
              );
  
	return 0;
}



/*! Contingency Loci
 
 Analyze lengths of homopolymer repeats in mixed samples.
 
 */
int do_tabulate_contingency_loci(int argc, char* argv[]) {
	
	// setup and parse configuration options:
	AnyOption options("Usage: breseq CL-TABULATE [-s -b data/reference.bam -f data/reference.fasta -r data/reference.gff -o contingency_loci.csv -m 8]");
  options.addUsage("");
  options.addUsage("Tabulate the frequencies of length distributions observed for putative contingency loci (homopolymer tracts).");
  options.addUsage("");
  options.addUsage("Allowed Options");
  options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("bam,b", "BAM file containing aligned read sequences", "data/reference.bam");
  options("fasta,f", "FASTA file of reference sequence", "data/reference.fasta");
  options("reference,r","File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files", "data/reference.gff3");
  options("output,o", "Output CSV file", "contingency_loci.csv");
  options("minimum-length,m", "Minimum length of a homopolymer tract in the reference genome to consider a putative contingency locus", "8");
  //options("loci,l", "Contingency loci coordinates", "");
  options("strict,s", "exclude non-perfect matches in surrounding 5 bases", TAKES_NO_ARGUMENT);
	options.processCommandArgs(argc, argv);
  
	if(options.count("help")
		 ) {
		options.printAdvancedUsage();
		return -1;
	}
  
  if ( !file_exists(options["bam"].c_str()) || !file_exists(options["fasta"].c_str()) || !file_exists(options["reference"].c_str()) ) {
    options.addUsage("");
    options.addUsage("Provide valid bam, fasta, and reference file paths.");
    options.printUsage();
    return -1;
  }
  
	// attempt to calculate error calibrations:
  analyze_contingency_loci(
                           options["bam"],
                           options["fasta"],
                           from_string<vector<string> >(options["reference"]),
                           options["output"],
                           "", //options["loci"], -->leave this blank
                           from_string<int32_t>(options["minimum-length"]),
                           options.count("strict")
                           );
	
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


int do_simulate_read(int argc, char *argv[])
{
  AnyOption options("Usage: breseq SIMULATE-READ -g <genome diff> -r <reference file> -c <average coverage> -o <output file>");

  options
  ("genome_diff,g", "Genome diff file.")
  ("reference,r",   "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)")
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
  ref_seq_info.LoadFiles(make_vector<string>(ref_file_name));
  new_ref_seq_info.LoadFiles(make_vector<string>(ref_file_name));


  //! Step: Apply genome diff mutations to reference sequence.
  const string &gd_file_name = options["genome_diff"];
  cGenomeDiff gd(gd_file_name);

  gd.apply_to_sequences(ref_seq_info, new_ref_seq_info, verbose);
  
  //! Write applied GFF3 file if requested.
  if(options.count("gff3"))new_ref_seq_info.WriteGFF(options["output"] + ".gff3");


  const cAnnotatedSequence &sequence = new_ref_seq_info[0];
  uint32_t coverage = from_string<uint32_t>(options["coverage"]);
  uint32_t read_size = from_string<uint32_t>(options["length"]);

  uint32_t n_reads = (sequence.get_sequence_length() / read_size) * coverage;

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
    ref_seq_info.LoadFiles(make_vector<string>(settings.reference_gff3_file_name));

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


// Analyze biases in the coverage - must be done after a breseq run
int do_coverage_bias(int argc, char *argv[])
{
  // setup and parse configuration options:
  AnyOption options("Usage: breseq COVERAGE-BIAS [-b reference.bam -f reference.fasta -o alignment.html -n 200] region1 [region2 region3 ...]");
  options.addUsage("");
  options.addUsage("Display reads aligned to the specified region or regions.");
  options.addUsage("");
  options.addUsage("Allowed Options");
  options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("bam,b", "BAM database file of read alignments", "data/reference.bam");
  options("fasta,f", "FASTA file of reference sequences", "data/reference.fasta");
  options("read-length,l", "Length of read used to calculate %GC", 200);
  options("output,o", "Output file prefix", "coverage_bias");
  options.processCommandArgs(argc, argv);
  
  if(options.count("help")) {
    options.printAdvancedUsage();
    return -1;
  }
  
  if (!file_exists(options["fasta"].c_str())) {
    options.addUsage("");
    options.addUsage("Could not open input reference FASTA file (-f):\n  " + options["fasta"]);
    options.printUsage();
    return -1;
  }
  
  if (!file_exists(options["bam"].c_str())) {
    options.addUsage("");
    options.addUsage("Could not open input BAM file of aligned reads (-b):\n  " + options["bam"]);
    options.printUsage();
    return -1;
  }
  
  //(re)load the reference sequences from our converted files
  //cReferenceSequences ref_seq_info;
  //ref_seq_info.LoadFiles(make_vector<string>(settings.reference_gff3_file_name));
    
  CoverageDistribution::analyze_coverage_bias(
        options["fasta"].c_str(),
        options["bam"].c_str(),
        options["output"],
        from_string<int32_t>(options["read-length"])
        );
  
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

int do_assemble_unmatched(int argc, char* argv[])
{
  AnyOption options("Usage: breseq ASSEMBLE-UNMATCHED-PAIRS [-o unmatched_assembly] reads1.fastq [reads2.fastq ...]");  
  options.addUsage("Assembles the unmatched reads from a breseq run given paired-end or mate-paired data.");
  options.addUsage("This command must be run from the main results directory of a breseq run (i.e., it must contain a data directory),");
  options.addUsage("You must provide the exact set of original fastq read files used in the breseq run. It is assumed that each set of two read files, in order, contain the first and second reads from pairs.");
  options.addUsage("Output is in directory: unmatched_assembly");
  options("output,o","Main directory containing output from the breseq run. A directory within this called unmatched_assembly will be created for the output of this command.", ".");
  options("verbose,v","Verbose output", TAKES_NO_ARGUMENT);
  options.addUsage("");
  options.addUsage("Example assembly commands:");
  options.addUsage("  velveth assemble_unmatched 51 -fastq -shortPaired assemble_unmatched/1.fastq assemble_unmatched/2.fastq");
  options.addUsage("  velvetg assemble_unmatched -exp_cov auto -cov_cutoff auto");

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
    uint32_t min_read_length;
    uint32_t max_read_length;
    uint8_t min_quality_score;
    uint8_t max_quality_score;
    string quality_format = cFastqQualityConverter::predict_fastq_file_format(read_file_name_1, original_num_reads, original_num_bases, min_read_length, max_read_length, min_quality_score, max_quality_score);
    
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
		SequenceConversionSummary s;
    cReferenceSequences conv_ref_seq_info;
    
    // Load all of the reference sequences and convert to FASTA and GFF3
    conv_ref_seq_info.LoadFiles(settings.all_reference_file_names);
    conv_ref_seq_info.WriteFASTA(settings.reference_fasta_file_name);
    conv_ref_seq_info.WriteGFF(settings.reference_gff3_file_name);
    s.total_reference_sequence_length = conv_ref_seq_info.total_length();
    
    // Do a quick load of the file to detect formatting errors.
    if (settings.user_evidence_genome_diff_file_name != "") {
      cGenomeDiff gd(settings.user_evidence_genome_diff_file_name);
      gd.valid_with_reference_sequences(conv_ref_seq_info, false);
    }
    
    // No conversion if already is sam mode
    if (!settings.aligned_sam_mode) {
      
      //Check the FASTQ format and collect some information about the input read files at the same time
      uint32_t overall_min_read_length = UNDEFINED_UINT32;
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
        AnalyzeFastqSummary s_rf = normalize_fastq(fastq_file_name,
                                                     convert_file_name,
                                                     i+1,
                                                     settings.quality_score_trim,
                                                     !settings.no_read_filtering,
                                                     s.num_bases,
                                                     read_file_base_limit,
                                                     settings.read_file_min_read_length,
                                                     settings.read_file_max_same_base_fraction,
                                                     settings.read_file_max_N_fraction
                                                     );
        settings.track_intermediate_file(settings.alignment_correction_done_file_name, convert_file_name);
        
        // Record statistics
        if ((overall_min_read_length == UNDEFINED_UINT32) || (s_rf.min_read_length > overall_min_read_length))
          overall_min_read_length = s_rf.max_read_length;
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
      
      s.avg_read_length = static_cast<double>(s.num_bases) / static_cast<double>(s.num_reads);
      s.min_read_length = overall_min_read_length;
      s.max_read_length = overall_max_read_length;
      s.max_qual = overall_max_qual;
      summary.sequence_conversion = s;
      
      // Print totals
      cerr << "  ::TOTAL::" << endl;
      cerr << "    Original reads: " << s.original_num_reads << " bases: " << s.original_num_bases << endl;
      cerr << "    Analyzed reads: " << s.num_reads << " bases: " << s.num_bases << endl;
      
    }
    
      
		// create SAM faidx
    /*
		string samtools = settings.ctool("samtools");
		string command = samtools + " faidx " + settings.reference_fasta_file_name;
		SYSTEM(command);
    */
    samtools_faidx(settings.reference_fasta_file_name);
    
		// calculate trim files
		calculate_trims(settings.reference_fasta_file_name, settings.sequence_conversion_path);
    settings.track_intermediate_file(settings.output_done_file_name, settings.sequence_conversion_path + "/*.trims");

		// store summary information
		summary.sequence_conversion.store(settings.sequence_conversion_summary_file_name);
		settings.done_step(settings.sequence_conversion_done_file_name);
	}

	summary.sequence_conversion.retrieve(settings.sequence_conversion_summary_file_name);
	ASSERT(summary.sequence_conversion.max_read_length != UNDEFINED_UINT32, "Can't retrieve max read length from file: " + settings.sequence_conversion_summary_file_name);

  // (re)load the reference sequences from our converted files
  // we must be sure to associate them with their original file names
  // so that contig and junction-only references are correctly flagged
  ref_seq_info.LoadFiles(make_vector<string>(settings.reference_gff3_file_name));
  ref_seq_info.use_original_file_names();
  
  // update the normal versus junction-only lists
  // - must be done after reading from our own GFF3 file
  settings.init_reference_sets(ref_seq_info);
  
  // Calculate the total reference sequence length
  summary.sequence_conversion.total_reference_sequence_length = ref_seq_info.total_length();
  
  
  // Reload certain information into settings from summary to make re-entrant
  // Aligned SAM mode knows nothing of read limits, so unused.
  
  if (!settings.aligned_sam_mode) {
    map<string, bool> read_files_surviving_conversion;
    
    // --> converted names
    for (map<string, AnalyzeFastqSummary>::iterator it = summary.sequence_conversion.reads.begin(); it != summary.sequence_conversion.reads.end(); it++)
    {
      string base_name = it->first;
      read_files_surviving_conversion[base_name] = true;

      if (it->second.converted_fastq_name.size() > 0) {
        settings.read_files.read_file_to_converted_fastq_file_name_map[base_name] = it->second.converted_fastq_name;
      }
    }
    
    // --> remove any read files that were not used due to coverage limits
    for (cReadFiles::iterator it=settings.read_files.begin(); it != settings.read_files.end(); ) {
      if (read_files_surviving_conversion.count(it->base_name())) {
        it++;
      } else {
        settings.read_files.read_file_to_fastq_file_name_map.erase(it->file_name());
        settings.read_files.read_file_to_converted_fastq_file_name_map.erase(it->file_name());
        it = settings.read_files.erase(it);
      }
    }
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
      + settings.bowtie2_score_parameters + " " + settings.bowtie2_genome_alignment_reporting_parameters + " " + settings.bowtie2_min_score_stringent + " --reorder -x " + reference_hash_file_name + " -U " + read_fastq_file + " -S " + stage1_reference_sam_file_name + " --un " + stage1_unmatched_fastq_file_name;
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
        
        // If we are doing staged alignment -- only align the unmatched reads and save to different name initially
        string read_fastq_file = settings.file_name(settings.stage1_unmatched_fastq_file_name, "#", base_read_file_name);
        string reference_sam_file_name = settings.file_name(settings.stage2_reference_sam_file_name, "#", base_read_file_name);
          
        uint32_t bowtie2_seed_substring_size_relaxed = 5 + trunc(summary.sequence_conversion.reads[settings.read_files[i].base_name()].avg_read_length * 0.1);
        // Check bounds
        bowtie2_seed_substring_size_relaxed = max<uint32_t>(9, bowtie2_seed_substring_size_relaxed);
        bowtie2_seed_substring_size_relaxed = min<uint32_t>(31, bowtie2_seed_substring_size_relaxed);
        
        string command = "bowtie2 -t -p " + s(settings.num_processors) + " --local " + " -L " + to_string<uint32_t>(bowtie2_seed_substring_size_relaxed) 
          + " " + settings.bowtie2_score_parameters + " " + settings.bowtie2_genome_alignment_reporting_parameters + " " + settings.bowtie2_min_score_relaxed + " --reorder -x " + reference_hash_file_name + " -U " + read_fastq_file + " -S " + reference_sam_file_name;
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
      string coverage_junction_best_bam_unsorted_file_name = settings.coverage_junction_best_bam_unsorted_file_name;

      samtools_import(reference_faidx_file_name, preprocess_junction_best_sam_file_name, coverage_junction_best_bam_unsorted_file_name);
      
      samtools_sort(coverage_junction_best_bam_unsorted_file_name, coverage_junction_best_bam_file_name, settings.num_processors);
      
      samtools_index(coverage_junction_best_bam_file_name);
      
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
      
      /*
			string samtools = settings.ctool("samtools");
			string faidx_command = samtools + " faidx " + settings.candidate_junction_fasta_file_name;
			if (!file_empty(settings.candidate_junction_fasta_file_name.c_str()))
				SYSTEM(faidx_command);
      */
      if (!file_empty(settings.candidate_junction_fasta_file_name.c_str()))
        samtools_faidx(settings.candidate_junction_fasta_file_name);
      
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

			/// create index
			string candidate_junction_hash_file_name = settings.candidate_junction_hash_file_name;
			string candidate_junction_fasta_file_name = settings.candidate_junction_fasta_file_name;

			if (!file_empty(candidate_junction_fasta_file_name.c_str()))
			{
        string command = "bowtie2-build -q " + candidate_junction_fasta_file_name + " " + candidate_junction_hash_file_name;
        SYSTEM(command);
        settings.track_intermediate_file(settings.candidate_junction_alignment_done_file_name, candidate_junction_hash_file_name + "*");
			}

			/// align reads to candidate junction sequences
			for (uint32_t i = 0; i < settings.read_files.size(); i++)
			{
				string base_read_file_name = settings.read_files[i].m_base_name;
				string candidate_junction_sam_file_name = settings.file_name(settings.candidate_junction_sam_file_name, "#", base_read_file_name);

				string read_fastq_file = settings.base_name_to_read_file_name(base_read_file_name);

        string filename = candidate_junction_hash_file_name + ".1.bt2";
        if (!file_exists(filename.c_str())) 
          continue;
    
        /// NEW CODE mapping to junctions with somewhat relaxed parameters
        uint32_t bowtie2_seed_substring_size_junction = trunc(summary.sequence_conversion.reads[settings.read_files[i].base_name()].avg_read_length * 0.3);
        // Check bounds
        bowtie2_seed_substring_size_junction = max<uint32_t>(9, bowtie2_seed_substring_size_junction);
        bowtie2_seed_substring_size_junction = min<uint32_t>(31, bowtie2_seed_substring_size_junction);
        
        string command = "bowtie2 -t -p " + s(settings.num_processors) + " --local " + " -L " + to_string<uint32_t>(bowtie2_seed_substring_size_junction) + " "
        + settings.bowtie2_score_parameters + " " + settings.bowtie2_junction_alignment_reporting_parameters + " " + settings.bowtie2_min_score_junction + " --reorder -x " + candidate_junction_hash_file_name + " -U " + read_fastq_file + " -S " + candidate_junction_sam_file_name;
        
        SYSTEM(command);
        
        settings.track_intermediate_file(settings.alignment_correction_done_file_name, candidate_junction_sam_file_name + "*");
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
		string junction_bam_file_name = settings.junction_bam_file_name;

    //string samtools = settings.ctool("samtools");
    //string command;

    // only run samtools if we are predicting junctions and there were results in the sam file
    // first part of conditional really not necessary @JEB
		if (!file_empty(resolved_junction_sam_file_name.c_str()))
		{
      /*
			command = samtools + " import " + candidate_junction_faidx_file_name + " " + resolved_junction_sam_file_name + " " + junction_bam_unsorted_file_name;
			SYSTEM(command);
			command = samtools + " sort " + junction_bam_unsorted_file_name + " " + junction_bam_prefix;
      SYSTEM(command);
			if (!settings.keep_all_intermediates)
				remove_file(junction_bam_unsorted_file_name.c_str());
			command = samtools + " index " + junction_bam_file_name;
			SYSTEM(command);
       */
      
      samtools_import(candidate_junction_faidx_file_name, resolved_junction_sam_file_name, junction_bam_unsorted_file_name);
      samtools_sort(junction_bam_unsorted_file_name, junction_bam_file_name, settings.num_processors);
      if (!settings.keep_all_intermediates)
        remove_file(junction_bam_unsorted_file_name.c_str());
      samtools_index(junction_bam_file_name);
      
		}

		string resolved_reference_sam_file_name = settings.resolved_reference_sam_file_name;
		string reference_bam_unsorted_file_name = settings.reference_bam_unsorted_file_name;
		string reference_bam_file_name = settings.reference_bam_file_name;

    /*
		command = samtools + " import " + reference_faidx_file_name + " " + resolved_reference_sam_file_name + " " + reference_bam_unsorted_file_name;
    SYSTEM(command);
		command = samtools + " sort " + reference_bam_unsorted_file_name + " " + reference_bam_prefix;
    SYSTEM(command);
		if (!settings.keep_all_intermediates)
			remove_file(reference_bam_unsorted_file_name.c_str());
		command = samtools + " index " + reference_bam_file_name;
    SYSTEM(command);
    */
    
    samtools_import(reference_faidx_file_name, resolved_reference_sam_file_name, reference_bam_unsorted_file_name);
    samtools_sort(reference_bam_unsorted_file_name, reference_bam_file_name, settings.num_processors);
    if (!settings.keep_all_intermediates)
      remove_file(reference_bam_unsorted_file_name.c_str());
    samtools_index(reference_bam_file_name);
    
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
    if (!settings.no_deletion_prediction) {
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
      
      
      // If the fit failed or there was insufficient coverage for some reference sequences
      //    deletion_coverage_propagation_cutoff == -1 for insufficient coverage
      //    nbinom_mean_parameter == 0 for failed fit
      vector<string> failed_fit_seq_ids;
      vector<string> no_coverage_seq_ids;

      for (uint32_t i = 0; i < ref_seq_info.size(); i++) {
        string seq_id = ref_seq_info[i].m_seq_id;
        
        // For junction-only sequences, don't provide these warnings
        if (settings.call_mutations_seq_id_set().count(seq_id)) {
          if (summary.unique_coverage[seq_id].deletion_coverage_propagation_cutoff <= 0) {
            no_coverage_seq_ids.push_back(seq_id);
          } else if (summary.unique_coverage[seq_id].nbinom_mean_parameter == 0) {
            failed_fit_seq_ids.push_back(seq_id);
          }
        }
      }
      
      if (failed_fit_seq_ids.size()) {
        
        string warning_string = "Failed to fit coverage distribution for some reference sequences. This may degrade the quality of predicting mutations from new sequence junctions (JC evidence).";
        if (!settings.deletion_coverage_propagation_cutoff) {
          warning_string += " You may want to set --deletion-coverage-propagation-cutoff to improve the quality of deletion prediction (MC evidence).";
        }
        warning_string += "\n\nFailed fit seq ids: " + join(failed_fit_seq_ids, ", ");
        
        WARN(warning_string);
      }
      
      if (no_coverage_seq_ids.size()) {
        
        string warning_string = "Insufficient coverage to call mutations for some reference sequences. Set either the --targeted-sequencing or --contig-reference option if you want mutations called on these reference sequences.";
        warning_string += "\n\nInsufficient coverage seq ids: " + join(no_coverage_seq_ids, ", ");
        
        WARN(warning_string);
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
        //cout << ref_seq_info[i].m_seq_id << " " << to_string<double>(deletion_propagation_cutoffs.back()) << " " << to_string<double>(deletion_seed_cutoffs.back()) << endl;
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
        settings.polymorphism_precision_decimal,
        settings.polymorphism_precision_places,
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
      
      
      settings.done_step(settings.copy_number_variation_done_file_name);
                       
    } // End of if cnv done file
  
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

    cGenomeDiff ra_mc_gd(settings.ra_mc_genome_diff_file_name);
    test_RA_evidence(ra_mc_gd, ref_seq_info, settings);
    
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

    // Add metadata that will only be in the output.gd file
    
    //#=REFSEQ header lines.
    mpgd.metadata.ref_seqs = settings.all_reference_file_names;
    
    //#=READSEQ header lines.
    mpgd.metadata.read_seqs.resize(settings.read_files.size());
    for (size_t i = 0; i < settings.read_files.size(); i++) {
      mpgd.metadata.read_seqs[i] = settings.read_files[i].file_name();
    }
    
    // These fields will be overwritten by any header GenomeDiff provided
    mpgd.metadata.title = settings.run_name;
    
    // Load metadata from existing header GenomeDiff file
    // This is purposefully done AFTER the previous metadata defs
    // in order to overwrite those (but not others)
    if (settings.header_genome_diff_file_name.size() > 0) {
      cGenomeDiff header_gd(settings.header_genome_diff_file_name);
      mpgd.metadata = header_gd.metadata;
    }
    
    // Fields here and below will write over any values in the header GenomeDiff file
    mpgd.metadata.program = string(PACKAGE_STRING) + " " + string(HG_REVISION);
    mpgd.metadata.command = settings.full_command_line;
    mpgd.metadata.created = Settings::time2string(time(NULL));
    
    // Add additional SUMMARY metadata
    
    mpgd.add_breseq_data("INPUT-READS", to_string(summary.sequence_conversion.original_num_reads));
    mpgd.add_breseq_data("INPUT-BASES", to_string(summary.sequence_conversion.original_num_bases));
    mpgd.add_breseq_data("CONVERTED-READS", to_string(summary.alignment_resolution.total_reads));
    mpgd.add_breseq_data("CONVERTED-BASES", to_string(summary.alignment_resolution.total_bases));
    mpgd.add_breseq_data("MAPPED-READS", to_string(static_cast<int64_t>(summary.alignment_resolution.total_reads_mapped_to_references)));
    mpgd.add_breseq_data("MAPPED-BASES", to_string(summary.alignment_resolution.total_bases_mapped_to_references));
    
    /*
    // Add additional header lines if needed.
    if (settings.add_metadata_to_gd){
      for (storable_map<string, Summary::Coverage>::iterator it = summary.unique_coverage.begin();
            it != summary.unique_coverage.end(); it ++) {
         //Usually needed for gathering breseq data.
       }
    }
    */

    // Write the output.gd file
    cerr << "  Writing final GD file..." << endl;
    mpgd.write(settings.final_genome_diff_file_name);
    
    //Don't reload -- we lose invisible fields that we need
    //cGenomeDiff gd(settings.final_genome_diff_file_name);
    cGenomeDiff gd = mpgd;
    
    // Empty metadata... necessary to keep consistency tests
    // when things like time are being added to output.gd
    gd.metadata = cGenomeDiff::Metadata();
    
    // Write VCF conversion
    cerr << "  Writing final VCF file..." << endl;

    mpgd.write_vcf(settings.output_vcf_file_name, ref_seq_info);
    
    // Write a final JSON file with all summary information
    summary.store(settings.data_summary_file_name);
    
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
		// --- must occur after marking entries no_show
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

  signal(SIGSEGV, seg_fault_handler);
  
  Settings::set_global_paths(argc, argv);
  
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
    Settings::command_line_run_header();
    breseq_default_action(argc, argv); // Gives default usage in this case.
    return -1; 
	}

	// Pass the command to the proper handler
	command = to_upper(command);
  
  
  // gnu standard return version string
  if ( (argc_new == 1) && (command == "--VERSION") ) {
    cout << "breseq " << VERSION << endl;
    return 0;
  }
  
  // Print out our generic header
  Settings::command_line_run_header();
  
  // Sequence COMMANDs:
  if (command == "CONVERT-FASTQ") {
		return do_convert_fastq(argc_new, argv_new);
	} else if (command == "CONVERT-REFERENCE") {
		return do_convert_reference(argc_new, argv_new);
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
  } else if (command == "CNV") {
    return do_copy_number_variation(argc_new, argv_new);
  } else if (command == "COVERAGE-BIAS"){
    return do_coverage_bias(argc_new, argv_new);
  } else if (command == "ERROR_COUNT") {
    return do_error_count(argc_new, argv_new);
  } else if (command == "JUNCTION-POLYMORPHISM") {
    return do_junction_polymorphism(argc_new, argv_new);
  } else if (command == "CL-TABULATE") {
    return do_tabulate_contingency_loci(argc_new, argv_new);
  } else if (command == "CL-SIGNIFICANCE") {
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
