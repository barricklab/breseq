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

// Autoconf header
#include <config.h>

#include "libbreseq/settings.h"

#include "libbreseq/anyoption.h"


using namespace std;

namespace breseq
{
  
  const char* kBreseqAlignmentScoreBAMTag = "X5";
  const char* kBreseqBestAlignmentScoreBAMTag = "X6";

  string Settings::output_divider("====================================================================================");
  
	void cReadFiles::Init(const vector<string>& read_file_names, bool sam_files)
	{
		//clear any existing info
		this->clear();

		uint32_t on_id = 0;
		uint32_t on_error_group = 0;
		uint32_t on_paired_end_group = 0;

		for (vector<string>::const_iterator it = read_file_names.begin(); it < read_file_names.end(); it++)
		{
			cReadFile rf;
			rf.m_original_file_name = *it;

			rf.m_paired_end_group = on_paired_end_group++;
			rf.m_error_group = on_error_group++;
			rf.m_id = on_id++;

			// create base name
			rf.m_base_name = rf.m_original_file_name;
			// - beginning path
			size_t pos = rf.m_base_name.rfind("/");
			if (pos != string::npos) rf.m_base_name.erase(0, pos + 1);
			// - trailing .fastq, must be at the end of the sequence
      if (!sam_files) {
        pos = rf.m_base_name.rfind(".fastq");
        if ((pos != string::npos) && (pos = rf.m_base_name.size() - 6))
          rf.m_base_name.erase(pos);
      } else {
        pos = rf.m_base_name.rfind(".sam");
        if ((pos != string::npos) && (pos = rf.m_base_name.size() - 4))
          rf.m_base_name.erase(pos);
      }
			
      // set up the map for converting base names to fastq file names to be used
      // check to see whether this is a duplicate
      ASSERT(read_file_to_fastq_file_name_map[rf.m_base_name].size() == 0, 
              "Read file provided multiple times:\n1)" + read_file_to_fastq_file_name_map[rf.m_base_name] + "\n2)" + rf.m_original_file_name);
      read_file_to_fastq_file_name_map[rf.m_base_name] = rf.m_original_file_name;
      
			this->push_back(rf);
		}
	}
  
  string cReadFiles::base_name_to_read_file_name(const string& base_name)
  {    
    if (read_file_to_converted_fastq_file_name_map.count(base_name)) 
    {
      return read_file_to_converted_fastq_file_name_map[base_name];
    }
    
    assert(read_file_to_fastq_file_name_map.count(base_name));
    
    return read_file_to_fastq_file_name_map[base_name];
  }
  
  // Set up defaults and build paths without command-line arguments
  Settings::Settings(const string& _base_output_path)
  {
    // things work fine if this is empty ""
    this->base_output_path = _base_output_path;
    
    this->pre_option_initialize();
    // no command line arguments here
		this->post_option_initialize();
	}
  
  Settings::Settings(int argc, char* argv[])
  {
    this->pre_option_initialize(argc, argv);
    
    // We need the path to the executable to locate scripts and our own installed versions of binaries
    // --> get this before parsing the options
    this->bin_path = getExecPath(argv[0]);
    size_t slash_pos = this->bin_path.rfind("/");
    if (slash_pos != string::npos) this->bin_path.erase(slash_pos);

    
    // setup and parse configuration options:
    AnyOption options("Usage: breseq -r reference.gbk [-r reference2.gbk ...] reads1.fastq [reads2.fastq reads3.fastq ...]");
    
    options
		("help,h", "Produce help message showing advanced options", TAKES_NO_ARGUMENT)
    ("verbose,v","Produce verbose output",TAKES_NO_ARGUMENT, ADVANCED_OPTION)
		("output,o", "Path to breseq output", ".")
		("reference,r", "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files. (REQUIRED)")
    ("name,n", "Human-readable name of sample/run for output [empty]", "")
    ("num-processors,j", "Number of processors to use in multithreaded steps", 1)
    ("aligned-sam", "Input files are aligned SAM files, rather than FASTQ files. Junction prediction steps will be skipped.", TAKES_NO_ARGUMENT)
    ;
    
    options.addUsage("Special Reference Sequences", true);
    options
    ("junction-only-reference,s", "File containing reference sequences in GenBank, GFF3, or FASTA format. These references are only used for calling junctions with other reference sequences. An example of appropriate usage is including a transposon sequence not present in a reference genome. Option may be provided multiple times for multiple files.", NULL, ADVANCED_OPTION)
    ("targeted-sequencing,t", "Reference sequences were targeted for ultra-deep sequencing (using pull-downs or amplicons). Do not fit coverage distribution.", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ;
    
    options.addUsage("Junction Options", true);
    options
    ("base-quality-cutoff,b", "Ignore bases with quality scores lower than this value", 3, ADVANCED_OPTION)
    ("quality-score-trim", "Trim the ends of reads past any base with a quality score below --base-quality-score-cutoff.", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("require-match-length", "Only consider alignments that cover this many bases of a read", 0, ADVANCED_OPTION)
    ("require-match-fraction", "Only consider alignments that cover this fraction of a read", 0.9)
    ("deletion-coverage-propagation-cutoff","Value for coverage above which deletions are cutoff. 0 = calculated from coverage distribution", 0, ADVANCED_OPTION)
    ("deletion-coverage-seed-cutoff","Value for coverage below which deletions are seeded", 0, ADVANCED_OPTION)
    ("mutation-score-cutoff", "Log10 E-value cutoff for base-substitution and micro-indel predictions", 10, ADVANCED_OPTION)
    ;
    
    options.addUsage("Junction Options", true);
    options
    ("no-junction-prediction", "Do not predict new sequence junctions", TAKES_NO_ARGUMENT)
    ("junction-alignment-pair-limit", "Only consider this many passed alignment pairs when creating candidate junction sequences", 100000, ADVANCED_OPTION)
    ("junction-score-cutoff", "Maximum negative log10 probability of uneven coverage across a junction breakpoint to accept (0=OFF)", 3.0, ADVANCED_OPTION)
    ("junction-minumum-pos-hash-score", "Minimum number of distinct spanning read start positions required to accept a junction", 3, ADVANCED_OPTION)
    ;
        
    options.addUsage("Polymorphism (Mixed Population) Options", true);
    options
    ("polymorphism-prediction,p", "Predict polymorphic (mixed) mutations", TAKES_NO_ARGUMENT)
    ("polymorphism-no-indels", "Do not predict insertion/deletion polymorphisms", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("polymorphism-reject-homopolymer-length", "Reject polymorphisms predicted in homopolymer repeats with this length or greater (DEFAULT= consensus mode, 3; polymorphism mode, 0) ", "", ADVANCED_OPTION)
    ("polymorphism-score-cutoff", "Log10 E-value cutoff for test of polymorphism vs no polymorphism (DEFAULT= consensus mode, 10; polymorphism mode, 2)", "", ADVANCED_OPTION)
    ("polymorphism-bias-cutoff", "P-value criterion for Fisher's exact test for strand bias AND K-S test for quality score bias (DEFAULT = OFF)", "", ADVANCED_OPTION)
    ("polymorphism-frequency-cutoff", "Only predict polymorphisms where both allele frequencies are > than this value (DEFAULT= consensus mode, 0.1; polymorphism mode, 0.0)", "", ADVANCED_OPTION)
    ("polymorphism-minimum-coverage-each-strand", "Only predict polymorphisms where this many reads on each strand support alternative alleles (DEFAULT= consensus mode, 2; polymorphism mode, 0", "", ADVANCED_OPTION)
    ;
    
    // CNV and Periodicity block
    options.addUsage("Experimental Options (Use at your own risk)", true);
    options
    ("user-junction-gd","User supplied genome diff file of JC evidence to use as candidate junctions and report support for.", "", ADVANCED_OPTION) 
    ("cnv","Do experimental copy number variation prediction",TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("cnv-tile-size", "Tile size for copy number variation prediction", 500, ADVANCED_OPTION)
    ("cnv-ignore-redundant", "Only consider non-redundant coverage when using cnv", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("per-position-file", "Create additional file of per-position aligned bases", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("periodicity", "Finding sum of differences squared of a coverage file", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("periodicity-method", "Which method to use for periodicity", 1, ADVANCED_OPTION)
    ("periodicity-start", "Start of offsets", 1, ADVANCED_OPTION)
    ("periodicity-end", "End of offsets", 2, ADVANCED_OPTION)
    ("periodicity-step", "Increment of offsets", 1, ADVANCED_OPTION)
    ;
    
    options.processCommandArgs(argc, argv);
    
    options.addUsage("");
		options.addUsage("Utility Command Usage: breseq [command] options ...");
    options.addUsage("  Sequence Utility Commands: CONVERT-FASTQ, CONVERT-REFERENCE, GET-SEQUENCE");
    options.addUsage("  Breseq Post-Run Commands: BAM2ALN, BAM2COV");
    options.addUsage("");
    options.addUsage("For help using a utility command, type: breseq [command] ");
    options.addUsage("");
    options.addUsage(output_divider);
    
    // make sure that the other config options are good:
    if (options.count("help"))
    {
      options.printAdvanced();
      exit(-1);
    }
    
    //! Settings: Global Workflow and Output
    
    this->base_output_path = options["output"];
    // Do this so the default printed can be pretty
    if (this->base_output_path == "current directory")
      this->base_output_path = ".";
    
    // Remaining command line items are read files
    // Read sequence file provided?
		if (options.getArgc() == 0)
		{
      options.addUsage("");
      options.addUsage("No read sequence files provided.");
      options.printUsage();
      exit(-1);		
    }
    for (int32_t i = 0; i < options.getArgc(); i++)
    {
      string read_file_name = options.getArgv(i);
      this->read_file_names.push_back(read_file_name);
    }
    this->aligned_sam_mode = options.count("aligned-sam");
    if (aligned_sam_mode) {
      cerr << "Input files are aligned SAM instead of FASTQ (--aligned-sam option)." << endl;
      cerr << "No junction prediction will take place." << endl;
      cerr << output_divider << endl;
    }
    this->read_files.Init(read_file_names, this->aligned_sam_mode);
    
    // Reference sequence provided?
		if (options.count("reference") == 0)
		{
      options.addUsage("");
      options.addUsage("No reference sequences provided (-r).");
      options.printUsage();
      exit(-1);
		}
    this->reference_file_names = from_string<vector<string> >(options["reference"]);
    
    // Important to check for NULL before converting
    if (options.count("junction-only-reference")) {
      this->junction_only_file_names = from_string<vector<string> >(options["junction-only-reference"]);
    }
    
    this->reference_file_names.insert(  this->reference_file_names.end(), 
                                        this->junction_only_file_names.begin(), 
                                        this->junction_only_file_names.end() 
                                      );    
    
    this->user_junction_genome_diff_file_name = options["user-junction-gd"];
    
    this->run_name = options["name"];
    
    this->num_processors = from_string<int32_t>(options["num-processors"]);
    
    this->do_copy_number_variation = options.count("cnv");
    this->copy_number_variation_tile_size = from_string<uint32_t>(options["cnv-tile-size"]);
    this->ignore_redundant_coverage = options.count("cnv-ignore-redundant");
    this->bowtie2 = options.count("bowtie2");
    this->bowtie2 = true; // testing
    this->bowtie2_align = options.count("bowtie2-align");
    this->bowtie2_align = true; // testing
    
    
    this->do_periodicity = options.count("periodicity");
    this->periodicity_method = from_string<uint32_t>(options["periodicity-method"]);
    this->periodicity_start = from_string<uint32_t>(options["periodicity-start"]);
    this->periodicity_end = from_string<uint32_t>(options["periodicity-end"]);
    this->periodicity_step = from_string<uint32_t>(options["periodicity-step"]);
    
    this->verbose = options.count("verbose");
    
    //! Settings: Read Alignment and Candidate Junction Read Alignment

    this->require_match_length = from_string<uint32_t>(options["require-match-length"]);
    this->require_match_fraction = from_string<double>(options["require-match-fraction"]);

    //! Settings: Mutation Identification
    this->base_quality_cutoff = from_string<uint32_t>(options["base-quality-cutoff"]);
    if (options.count("quality-score-trim")) this->quality_score_trim = this->base_quality_cutoff;

    this->deletion_coverage_propagation_cutoff = from_string<double>(options["deletion-coverage-propagation-cutoff"]);
    ASSERT(this->deletion_coverage_propagation_cutoff >= 0, "Argument --deletion-coverage-propagation-cutoff must be > 0")
    this->deletion_coverage_seed_cutoff = from_string<double>(options["deletion-coverage-seed-cutoff"]);
    ASSERT(this->deletion_coverage_propagation_cutoff >= 0, "Argument --deletion-coverage-seed-cutoff must be > 0")
		this->mutation_log10_e_value_cutoff = from_string<double>(options["mutation-score-cutoff"]);
    
    // Junction Prediction
    this->no_junction_prediction = options.count("no-junction-prediction");
    this->maximum_junction_sequence_passed_alignment_pairs_to_consider = from_string<uint64_t>(options["junction-alignment-pair-limit"]);
    this->junction_pos_hash_neg_log10_p_value_cutoff = from_string<double>(options["junction-score-cutoff"]);
    this->minimum_alignment_resolution_pos_hash_score = from_string<uint32_t>(options["junction-minumum-pos-hash-score"]);
    
    //
    // Set the read alignment evidence model
    // Different defaults for the three modes:
    //   Consensus - only consensus base calls (old default mode)
    //   Mixed base - consensus and strong polymorphism calls (new default mode)
    //   Polymorphism prediction - consensus and any mixed calls (--p options)
    //
    
    this->polymorphism_prediction = options.count("polymorphism-prediction");
    if (this->polymorphism_prediction) {
      
      // different default values
      this->polymorphism_frequency_cutoff = 0; // cut off if < X or > 1-X
      this->mixed_base_prediction = false;
      this->polymorphism_reject_homopolymer_length = 0;
      this->polymorphism_log10_e_value_cutoff = 2;
      this->polymorphism_bias_p_value_cutoff = 0.01;
      this->polymorphism_minimum_new_coverage_each_strand = 0;
      this->no_indel_polymorphisms = false;
      
      this->junction_pos_hash_neg_log10_p_value_cutoff = 0; // OFF
    }
    if (this->mixed_base_prediction) {
      // Not calculated for mixed base mode
      ASSERT(!options.count("polymorphism-bias-cutoff"), "Option --polymorphism-bias-cutoff requires --polymorphism-prediction.")

      this->polymorphism_frequency_cutoff = 0.1;
      this->mixed_base_prediction_marginal_frequency_cutoff = 0.5;
      this->no_indel_polymorphisms = false;
      this->polymorphism_reject_homopolymer_length = 3;
      this->polymorphism_log10_e_value_cutoff = 10;
      this->polymorphism_bias_p_value_cutoff = 0;
      this->polymorphism_minimum_new_coverage_each_strand = 2;
      this->no_indel_polymorphisms = false;

    }
    
    // override the default settings 
    if (options.count("polymorphism-frequency-cutoff"))
      this->polymorphism_frequency_cutoff = from_string<double>(options["polymorphism-frequency-cutoff"]);
    if (options.count("polymorphism-no-indels"))
      this->no_indel_polymorphisms = options.count("polymorphism-no-indels");
    if (options.count("polymorphism-reject-homopolymer-length"))
      this->polymorphism_reject_homopolymer_length = from_string<int32_t>(options["polymorphism-reject-homopolymer-length"]);
    if (options.count("polymorphism-score-cutoff"))
      this->polymorphism_log10_e_value_cutoff = from_string<double>(options["polymorphism-score-cutoff"]);
    if (options.count("polymorphism-bias-cutoff"))
      this->polymorphism_bias_p_value_cutoff = from_string<double>(options["polymorphism-bias-cutoff"]);
    if (options.count("polymorphism-minimum-coverage-each-strand"))
      this->polymorphism_minimum_new_coverage_each_strand = from_string<int32_t>(options["polymorphism-minimum-coverage-each-strand"]);
    
    
    this->targeted_sequencing = options.count("targeted-sequencing");
    
    this->print_mutation_identification_per_position_file = options.count("per-position-file");
    
		this->post_option_initialize();
    
    // Log the command line
    time_t stamp_time = time(NULL);
    this->log(ctime(&stamp_time));	
    this->log(this->full_command_line + "\n");
	}
  
  void Settings::command_line_run_header()
  {
    cerr << output_divider << endl;
    fprintf(stderr, "%s  %s   %s\n", PACKAGE_STRING, HG_REVISION, PACKAGE_URL);
    fprintf(stderr, "\n");
    fprintf(stderr, "Authors: Barrick JE, Borges JJ, Colburn GR, Knoester DB, Meyer AG, Reba A, Strand MD\n");
    fprintf(stderr, "Contact: %s\n", PACKAGE_BUGREPORT);
    fprintf(stderr, "\n");
    fprintf(stderr, "%s is free software; you can redistribute it and/or modify it under the\n", PACKAGE_NAME);
    fprintf(stderr, "terms the GNU General Public License as published by the Free Software \n");
    fprintf(stderr, "Foundation; either version 1, or (at your option) any later version.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Copyright (c) 2008-2010 Michigan State University\n");
    fprintf(stderr, "Copyright (c) 2011-2012 The University of Texas at Austin\n");
    cerr << output_divider << endl;
  }

	void Settings::pre_option_initialize(int argc, char* argv[])
	{    
    ////////////////////
    //! Data
    ////////////////////
    
    this->byline = "<b><i>breseq</i></b>&nbsp;&nbsp;version ";
    this->byline += PACKAGE_VERSION;
    this->byline += "&nbsp;&nbsp;";
    this->byline += HG_REVISION;

    this->website = PACKAGE_URL;
    
    // Settings
    // Initialize all variables that do not have default initializers (non-strings mostly) 

    this->arguments = "";
    if (argc > 0) this->full_command_line = argv[0];
    for(int i=1; i<argc; i++)
    {
      if (!this->arguments.empty()) this->arguments += " ";
      this->arguments += argv[i];
    }
    this->full_command_line += " " + this->arguments;
   
    ////////////////////
    //! Settings
    ////////////////////
    
    //! Settings: Global Workflow and Output
    
    //! Options that control which parts of the pipeline to execute
    this->no_read_filtering = false;
    this->no_junction_prediction = false;
		this->no_mutation_prediction = false;
		this->no_deletion_prediction = false;
    this->no_alignment_generation = false;
		this->do_copy_number_variation = false;
		this->do_periodicity = false;
    
    //! DEBUG options
    this->verbose = false;
		this->alignment_read_limit = 0;         
		this->resolve_alignment_read_limit = 0; 
    this->candidate_junction_read_limit = 0;
    this->no_unmatched_reads = false;
    this->keep_all_intermediates = false;

    //! Settings: Read Alignment and Candidate Junction Read Alignment
    this->ssaha2_seed_length = 13;
    this->ssaha2_skip_length = 1;
    this->require_match_length = 0;         
    this->require_match_fraction = 0.9; // CL value overrides
    this->maximum_read_mismatches = -1;

    this->bowtie2 = false;
    this->bowtie2_align = false;
    this->bowtie2_maximum_alignments_to_consider_per_read = 2000;
    
    this->bowtie2_score_parameters = "--ma 1 --mp 3 --np 0 --rdg 2,3 --rfg 2,3";
    this->bowtie2_score_parameters += (bowtie2_maximum_alignments_to_consider_per_read > 0) 
    ? " -k " + to_string(this->bowtie2_maximum_alignments_to_consider_per_read) : " -a";
    
    this->bowtie2_min_score_stringent = "-i S,1,0.25 --score-min L,0,0.9 "; // "-L 22 -i S,1,0.25 --score-min L,0,0.9 ";
    
    // Note: this leaves off -L, since it is set based on read length
    this->bowtie2_min_score_relaxed  = "-i S,1,0.25 --score-min L,6,0.2 "; // "-L 9 -i C,1,0 --score-min L,6,0.2 ";

    this->num_processors = 1;
    
    //! Settings: Candidate Junction Prediction
		this->preprocess_junction_min_indel_split_length = 3;
    this->required_both_unique_length_per_side_fraction = 0.2; 
    this->unmatched_end_length_factor =  1 - this->require_match_fraction;
    this->unmatched_end_minimum_read_length = 50;
    this->maximum_junction_sequence_negative_overlap_length_fraction = 0.4;
    this->maximum_junction_sequence_negative_overlap_length_minimum = 12;
    this->maximum_junction_sequence_positive_overlap_length_fraction = 0.4;
    this->maximum_junction_sequence_positive_overlap_length_minimum = 12;
    this->highly_redundant_junction_ignore_passed_pair_limit = 200;
    this->maximum_junction_sequence_passed_alignment_pairs_to_consider = 1000000;

    // Extra options that are mostly phased out because they are hard nucleotide cutoffs
    this->maximum_junction_sequence_insertion_length = 0;
    this->maximum_junction_sequence_overlap_length = 0;
    this->required_both_unique_length_per_side = 0;
    this->required_one_unique_length_per_side = 0;
    
    this->minimum_candidate_junction_pos_hash_score = 2;
    this->minimum_candidate_junctions = 100;
		this->maximum_candidate_junctions = 5000;
		this->maximum_candidate_junction_length_factor = 0.1;
    this->penalize_negative_junction_overlap = true;

    //! Settings: Alignment Resolution
    this->add_split_junction_sides = true;
    this->junction_pos_hash_neg_log10_p_value_cutoff = 2;
    this->minimum_alignment_resolution_pos_hash_score = 3;

    //! Settings: Mutation Identification
    this->base_quality_cutoff = 3;
    this->quality_score_trim = 0;
    this->deletion_coverage_propagation_cutoff = 0;
    this->deletion_coverage_seed_cutoff = 0;
    this->polymorphism_prediction = false;
    this->mixed_base_prediction = true;
    this->mixed_base_prediction_marginal_frequency_cutoff = 0.5;
    
    this->polymorphism_log10_e_value_cutoff = this->mutation_log10_e_value_cutoff;
		this->polymorphism_bias_p_value_cutoff = 0;
		this->polymorphism_frequency_cutoff = 0.1;
		this->polymorphism_minimum_new_coverage_each_strand = 0;
    this->polymorphism_reject_homopolymer_length = 0;
		this->no_indel_polymorphisms = false;
    
    //! Settings: Output
    this->maximum_reads_to_align = 100;
    this->max_rejected_read_alignment_evidence_to_show = 20;
		this->max_rejected_polymorphisms_to_show = 20;
		this->max_rejected_junction_evidence_to_show = 10;
		this->hide_circular_genome_junctions = true;
    
		this->lenski_format = false;
		this->no_evidence = false;
    this->add_metadata_to_gd = true;
    
    ////////////////////
    //! File Paths
    ////////////////////
    
    this->bin_path = ".";
    
    //ASSERT(this->required_both_unique_length_per_side <= this->ssaha2_seed_length,
    //       "--required-both-unique-length-per-side " + to_string(this->required_both_unique_length_per_side) + " must be less than or equal to --ssaha2-seed-length [" + to_string(this->ssaha2_seed_length) + "].");
	}

	void Settings::post_option_initialize()
	{     
    this->init_installed();
    
    // DATADIR is a preprocessor directive set by Automake in config.h
		this->program_data_path = DATADIR; 
        
    // Unless we are in "make test" mode where this environental variable is defined.
    char * breseq_data_path;
    breseq_data_path = getenv ("BRESEQ_DATA_PATH");
    if (breseq_data_path!=NULL) {
      this->program_data_path = breseq_data_path;
      cerr << "In test mode. Program data path: " << breseq_data_path << endl;
    }
    
    ////////////////////
    //! Settings
    ////////////////////
    
		// problems if there are spaces b/c shell removes quotes before we know about them
		// thus require run names to only use underscores (but when printing output, remove).
    this->print_run_name = substitute(this->run_name, "_", " ");
    
    ////////////////////
    //! File Paths
    ////////////////////
    
		//// '#' replaced with read file name
		//// '@' replaced by seq_id of reference sequence

    //! Paths: Sequence conversion
		this->sequence_conversion_path = "01_sequence_conversion";
		if (this->base_output_path.size() > 0) this->sequence_conversion_path = this->base_output_path + "/" + this->sequence_conversion_path;
		this->sequence_conversion_done_file_name = this->sequence_conversion_path + "/sequence_conversion.done";

		this->converted_fastq_file_name = this->sequence_conversion_path + "/#.converted.fastq";
		this->unwanted_fasta_file_name = this->sequence_conversion_path + "/unwanted.fasta";
		this->reference_trim_file_name = this->sequence_conversion_path + "/@.trims";
		this->sequence_conversion_summary_file_name = this->sequence_conversion_path + "/summary.bin";

    //! Paths: Read alignment
		this->reference_alignment_path = "02_reference_alignment";
		if (this->base_output_path.size() > 0) this->reference_alignment_path = this->base_output_path + "/" + this->reference_alignment_path;
		this->reference_alignment_done_file_name = this->reference_alignment_path + "/alignment.done";

		this->reference_hash_file_name = this->reference_alignment_path + "/reference";
		this->reference_sam_file_name = this->reference_alignment_path + "/#.reference.sam";

    this->stage1_reference_sam_file_name = this->reference_alignment_path + "/#.stage1.sam";
    this->stage1_unmatched_fastq_file_name = this->reference_alignment_path + "/#.stage1.unmatched.fastq";
    this->stage2_reference_sam_file_name = this->reference_alignment_path + "/#.stage2.matched.sam";

    //! Paths: Junction Prediction
		this->candidate_junction_path = "03_candidate_junctions";
		if (this->base_output_path.size() > 0) this->candidate_junction_path = this->base_output_path + "/" + this->candidate_junction_path;

		this->preprocess_junction_done_file_name = this->candidate_junction_path + "/preprocess_junction_alignment.done";
		this->preprocess_junction_best_sam_file_name = this->candidate_junction_path + "/best.sam";
		this->preprocess_junction_split_sam_file_name = this->candidate_junction_path + "/#.split.sam";
    
    this->coverage_junction_done_file_name = this->candidate_junction_path + "/coverage_junction_alignment.done";
		this->coverage_junction_best_bam_unsorted_file_name = this->candidate_junction_path + "/best.unsorted.bam";
		this->coverage_junction_best_bam_file_name = this->candidate_junction_path + "/best.bam";
		this->coverage_junction_best_bam_prefix = this->candidate_junction_path + "/best";
		this->coverage_junction_distribution_file_name = this->candidate_junction_path + "/@.unique_only_coverage_distribution.tab";
		this->coverage_junction_plot_file_name = this->candidate_junction_path + "/@.coverage.pdf";
		this->coverage_junction_summary_file_name = this->candidate_junction_path + "/coverage.summary.bin";
    this->coverage_junction_error_count_summary_file_name = this->candidate_junction_path + "/error_count.summary.bin";

    this->candidate_junction_done_file_name = this->candidate_junction_path + "/candidate_junction.done";
		this->candidate_junction_summary_file_name = this->candidate_junction_path + "/candidate_junction_summary.bin";
		this->candidate_junction_fasta_file_name = this->candidate_junction_path + "/candidate_junction.fasta";
		this->candidate_junction_faidx_file_name = this->candidate_junction_path + "/candidate_junction.fasta.fai";

    //! Paths: Junction Alignment
		this->candidate_junction_alignment_path = "04_candidate_junction_alignment";
		if (this->base_output_path.size() > 0) this->candidate_junction_alignment_path = this->base_output_path + "/" + this->candidate_junction_alignment_path;
		this->candidate_junction_alignment_done_file_name = this->candidate_junction_alignment_path + "/candidate_junction_alignment.done";

		this->candidate_junction_hash_file_name = this->candidate_junction_alignment_path + "/candidate_junction";
		this->candidate_junction_sam_file_name = this->candidate_junction_alignment_path + "/#.candidate_junction.sam";

    //! Paths: Alignment Resolution
		this->alignment_resolution_path = "05_alignment_correction";
		if (this->base_output_path.size() > 0) this->alignment_resolution_path = this->base_output_path + "/" + this->alignment_resolution_path;
		this->alignment_correction_done_file_name = this->alignment_resolution_path + "/alignment_resolution.done";

		this->resolved_reference_sam_file_name = this->alignment_resolution_path + "/reference.sam";
		this->resolved_junction_sam_file_name = this->alignment_resolution_path + "/junction.sam";
		this->alignment_resolution_summary_file_name = this->alignment_resolution_path + "/summary.bin";
		this->jc_genome_diff_file_name = this->alignment_resolution_path + "/jc_evidence.gd";

    //! Paths: BAM conversion
		this->bam_path = "06_bam";
		if (this->base_output_path.size() > 0) this->bam_path = this->base_output_path + "/" + this->bam_path;
		this->bam_done_file_name = this->bam_path + "/bam.done";

		this->reference_bam_unsorted_file_name = this->bam_path + "/reference.unsorted.bam";
		this->junction_bam_unsorted_file_name = this->bam_path + "/junction.unsorted.bam";
		this->junction_bam_prefix = this->bam_path + "/junction";
		this->junction_bam_file_name = this->bam_path + "/junction.bam";

		//! Paths: Error Calibration
		this->error_calibration_path = "07_error_calibration";
		if (this->base_output_path.size() > 0) this->error_calibration_path = this->base_output_path + "/" + this->error_calibration_path;
		this->error_counts_file_name = this->error_calibration_path + "/#.error_counts.tab";
		this->error_rates_done_file_name = this->error_calibration_path + "/error_rates.done";

		this->error_rates_file_name = this->error_calibration_path + "/error_rates.tab";
		this->error_counts_done_file_name = this->error_calibration_path + "/error_counts.done";
		this->coverage_file_name = this->error_calibration_path + "/@.coverage.tab";
		this->unique_only_coverage_distribution_file_name = this->error_calibration_path + "/@.unique_only_coverage_distribution.tab";
		this->error_rates_summary_file_name = this->error_calibration_path + "/summary.bin";
		this->error_rates_base_qual_error_prob_file_name = this->error_calibration_path + "/base_qual_error_prob.#.tab";
		this->plot_error_rates_r_script_file_name = this->program_data_path + "/plot_error_rate.r";
		this->plot_error_rates_fit_r_script_file_name = this->error_calibration_path + "/fit.#.r_script";
		this->plot_error_rates_r_script_log_file_name = this->error_calibration_path + "/#.plot_error_rate.log";

		//! Paths: Mutation Identification
		this->mutation_identification_path = "08_mutation_identification";
		if (this->base_output_path.size() > 0) this->mutation_identification_path = this->base_output_path + "/" + this->mutation_identification_path;

    this->mutation_identification_done_file_name = this->mutation_identification_path + "/mutation_identification.done";
		this->mutation_identification_per_position_file_name = this->mutation_identification_path + "/per_position_file.tab";
		this->complete_coverage_text_file_name = this->mutation_identification_path + "/@.coverage.tab";
		this->genome_error_counts_file_name = this->mutation_identification_path + "/error_counts.tab";
		this->ra_mc_genome_diff_file_name = this->mutation_identification_path + "/ra_mc_evidence.gd";

    this->polymorphism_statistics_done_file_name = this->mutation_identification_path + "/polymorphism_statistics.done";
		this->polymorphism_statistics_input_file_name = this->mutation_identification_path + "/polymorphism_statistics_input.tab";
		this->polymorphism_statistics_output_file_name = this->mutation_identification_path + "/polymorphism_statistics_output.tab";
		this->polymorphism_statistics_r_script_file_name = this->program_data_path + "/polymorphism_statistics.r";
		this->polymorphism_statistics_r_script_log_file_name = this->mutation_identification_path + "/polymorphism_statistics_output.log";
		this->polymorphism_statistics_ra_mc_genome_diff_file_name = this->mutation_identification_path + "/ra_mc_evidence_polymorphism_statistics.gd";

    //! Paths: Copy Number Variation
		this->copy_number_variation_path = "09_copy_number_variation";
    if (this->base_output_path.size() > 0) this->copy_number_variation_path = this->base_output_path + "/" + this->copy_number_variation_path;
    this->copy_number_variation_done_file_name = this->copy_number_variation_path + "/copy_number_variation.done";
    
    this->tiled_complete_coverage_text_file_name = this->copy_number_variation_path + "/@.tiled.tab";
    this->tiled_for_edging_text_file_name = this->copy_number_variation_path + "/@.tiled_for_edging.tab";
    this->ranges_text_file_name = this->copy_number_variation_path + "/@.ranges.tab";
    this->cnv_history_text_file_name = this->copy_number_variation_path + "/@.history.tab";
    this->smoothed_ranges_text_file_name = this->copy_number_variation_path + "/@.smoothed_ranges.tab";
    this->final_cnv_text_file_name = this->copy_number_variation_path + "/@.cnv_final.tab";
    this->copy_number_variation_cn_genome_diff_file_name = this->copy_number_variation_path + "/@.cn_evidence.gd";
    
    this->periodicity_done_file_name = this->copy_number_variation_path + "/@.periodicity.done";
    this->periodicity_table_file_name = this->copy_number_variation_path + "/@.periodicity.tab";
    
    //! Paths: Output
		this->output_path = "output";
		if (this->base_output_path.size() > 0) this->output_path = this->base_output_path + "/" + this->output_path;
		this->output_done_file_name = this->output_path + "/output.done";

		this->log_file_name = this->output_path + "/log.txt";
		this->index_html_file_name = this->output_path + "/index.html";
		this->summary_html_file_name = this->output_path + "/summary.html";
		this->marginal_html_file_name = this->output_path + "/marginal.html";

		this->local_evidence_path = "evidence";
		this->evidence_path = this->output_path + "/" + this->local_evidence_path;
		this->evidence_genome_diff_file_name = this->evidence_path + "/evidence.gd";
		this->final_genome_diff_file_name = this->output_path + "/output.gd";
		this->annotated_genome_diff_file_name = this->evidence_path + "/annotated.gd";
    
		this->local_coverage_plot_path = "evidence";
		this->coverage_plot_path = this->output_path + "/" + this->local_coverage_plot_path;
    this->coverage_plot_r_script_file_name = this->program_data_path + "/plot_coverage.r";    
		this->overview_coverage_plot_file_name = this->coverage_plot_path + "/@.overview.png";

		this->output_calibration_path = this->output_path + "/calibration";
		this->unique_only_coverage_plot_file_name = this->output_calibration_path + "/@.unique_coverage.pdf";
		this->error_rates_plot_file_name = this->output_calibration_path + "/#.error_rates.pdf";

		this->breseq_small_graphic_from_file_name = this->program_data_path + "/breseq_small.png";
		this->breseq_small_graphic_to_file_name = this->output_path + "/" + this->local_evidence_path + "/breseq_small.png";
    
    //! Paths: Data
		this->data_path = "data";
		if (this->base_output_path.size() > 0) this->data_path = this->base_output_path + "/" + this->data_path;

		this->reference_bam_prefix = this->data_path + "/reference";
		this->reference_bam_file_name = this->data_path + "/reference.bam";
		this->reference_fasta_file_name = this->data_path + "/reference.fasta";
		this->reference_faidx_file_name = this->data_path + "/reference.fasta.fai";
		this->reference_gff3_file_name = this->data_path + "/reference.gff3";
		this->unmatched_read_file_name = this->data_path + "/#.unmatched.fastq";

    //! Paths: Experimental
    this->long_pairs_file_name = this->output_path + "/long_pairs.tab";

  }
  
  void Settings::init_reference_sets(cReferenceSequences& ref_seq_info)
  {
    bool debug = false;
    
    set<string> junction_only_file_name_set(junction_only_file_names.begin(), junction_only_file_names.end());
    if (debug) cerr << "Initializing reference sets..." << endl;
    for(cReferenceSequences::iterator it=ref_seq_info.begin(); it!=ref_seq_info.end(); it++) {
      if (debug) cerr << it->get_file_name() << endl;
      if (junction_only_file_name_set.count(it->get_file_name())) {
        this->junction_only_seq_id_set.insert(it->m_seq_id);
        if (debug) cerr << "Junction only: " << it->m_seq_id << endl;
      }
      else {
        this->reference_seq_id_set.insert(it->m_seq_id);
        if (debug) cerr << "Normal: " << it->m_seq_id << endl;
      }
    }
  }


	void Settings::init_installed()
	{

    // Save the path for reference
    char * pPath;
    pPath = getenv("PATH");
    this->installed["path"] = "";
    if (pPath!=NULL)
    {
      this->installed["path"] = pPath;
    }
    
    // For checking path in XCode etc...
    //cerr << this->installed["path"] << endl;
    
    // SAMtools executables - look in the local bin path only
    string samtools_path = this->bin_path + "/samtools";
    if (!file_empty(samtools_path)) {
		  this->installed["samtools"] = samtools_path;
    }
    
    // Unless we are in "make test" mode where this environental variable is defined.
    char* breseq_samtools_path = getenv("BRESEQ_SAMTOOLS_PATH");
    if (breseq_samtools_path !=NULL) {
      if (!file_empty(string(breseq_samtools_path) + "/samtools")) {
        this->installed["samtools"] = string(breseq_samtools_path) + "/samtools";
      }
      cerr << "In test mode. Samtools path: " << breseq_samtools_path << endl;
    }
    
    if (this->installed["samtools"].size() != 0) {
      string version_string = SYSTEM_CAPTURE(this->installed["samtools"], true);
      size_t start_version_pos = version_string.find("Version:");
      if (start_version_pos != string::npos) {
        start_version_pos = version_string.find_first_not_of(" \t\r\n", start_version_pos+8);
        size_t end_version_pos = version_string.find_first_of("\r\n", start_version_pos+1);
        // This would just get the short version number
        //size_t end_version_pos = version_string.find_first_not_of("0123456789.", start_version_pos+1);
        this->installed["samtools_version_string"] = version_string.substr(start_version_pos, end_version_pos - start_version_pos);
      }
      else {
        this->installed["samtools_version_error_message"] = version_string;
      }
    }
    
    // which can return an error message    
    // detect SSAHA2 system-wide install
    this->installed["SSAHA2"] = SYSTEM_CAPTURE("which ssaha2", true);
    this->installed["SSAHA2_version_string"] = "";
    
    if (this->installed["SSAHA2"].size() != 0) {
      string version_string = SYSTEM_CAPTURE(this->installed["SSAHA2"] + " --version", true);
      size_t start_version_pos = version_string.find("version");
      if (start_version_pos != string::npos) {
        start_version_pos = version_string.find_first_not_of(" \t\r\n", start_version_pos+7);
        size_t end_version_pos = version_string.find_first_not_of("0123456789.", start_version_pos+1);
        this->installed["SSAHA2_version_string"] = version_string.substr(start_version_pos, end_version_pos - start_version_pos);
      }
      else {
        this->installed["SSAHA2_version_error_message"] = version_string;
      }
    }
    
    //cerr << "SSAHA2: " << this->installed["SSAHA2"] << " SSAHA2 version: " << this->installed["SSAHA2_version_string"] << endl;
    
    // detect bowtie2 and bowtie2-build system-wide install
    this->installed["bowtie2"] = SYSTEM_CAPTURE("which bowtie2", true);
    this->installed["bowtie2-build"] = SYSTEM_CAPTURE("which bowtie2-build", true);
    this->installed["bowtie2_version"] = "";
    this->installed["bowtie2_version_string"] = "";
    
    if (this->installed["bowtie2"].size() != 0) {
      string version_string = SYSTEM_CAPTURE(this->installed["bowtie2"] + " --version", true);
      size_t start_version_pos = version_string.find("version");
      if (start_version_pos != string::npos) {
        start_version_pos = version_string.find_first_not_of(" \t\r\n", start_version_pos+7);
        size_t end_version_pos = version_string.find_first_not_of("0123456789", start_version_pos+1);
        string new_version_string = version_string.substr(start_version_pos, end_version_pos - start_version_pos);
        
        start_version_pos = version_string.find_first_not_of(".", end_version_pos+1);
        end_version_pos = version_string.find_first_of(" \t\r\n", start_version_pos+1);
        
        new_version_string += "." + version_string.substr(start_version_pos, end_version_pos - start_version_pos);
        this->installed["bowtie2_version_string"] = new_version_string;
        new_version_string = substitute(new_version_string, "-", ".");
        
        // beta counts as another sub version
        if (new_version_string.find("beta") != string::npos) {
          new_version_string = substitute(new_version_string, "beta", "");
        } else {
          new_version_string += ".0";
        }
        
        vector<string> split_version = split(new_version_string, ".");
        uint32_t numerical_version = 0;
        for (vector<string>::iterator it=split_version.begin(); it != split_version.end(); it++) {
          numerical_version *= 1000;
          numerical_version += n(*it);
        }
        this->installed["bowtie2_version"] = s(numerical_version);
        //cout << this->installed["bowtie2_version"] << endl;
      }
      else {
        this->installed["bowtie2_version_error_message"] = version_string;
      }
    }
    
    //cerr << "Bowtie2: " << this->installed["bowtie2"] << " Bowtie2 version: " << this->installed["bowtie2_version_string"] << endl;

    
		this->installed["R"] = SYSTEM_CAPTURE("which R", true);
    this->installed["R_version_string"] = "";
    this->installed["R_version"] = "";
    
		if (this->installed["R"].size() > 0)
		{
			string R_version = SYSTEM_CAPTURE("R --version", true);
      
      // default if output does not match our pattern
      this->installed["R_version"] = "0";
      
      size_t start_version_pos = R_version.find("R version ");
      if (start_version_pos != string::npos)
      {
        start_version_pos+=10;
        size_t end_version_pos = R_version.find(" ", start_version_pos);
        
        if (end_version_pos == string::npos)
          end_version_pos = R_version.size();
        end_version_pos--;
        
        string version_string = R_version.substr(start_version_pos, end_version_pos - start_version_pos + 1);
        this->installed["R_version_string"] = version_string;
        vector<string> split_version_string = split(version_string, ".");
        if (split_version_string.size() == 3)
        {
          uint32_t R_numerical_version = from_string<uint32_t>(split_version_string[0]) * 1000000 + from_string<uint32_t>(split_version_string[1]) * 1000 + from_string<uint32_t>(split_version_string[2]);
          this->installed["R_version"] = to_string(R_numerical_version);
        }
      }
      else {
        this->installed["R_version_error_message"] = R_version;
        
      }
		}

    //cerr << "R: " << this->installed["R"] << " R version: " << this->installed["R_version_string"] << endl;
    
	}

	void Settings::check_installed()
	{
    // Developer's Note
    //
    // If you are running things through a debugger (like in XCode), your $PATH may need to be
    // set to include the paths where you have SSAHA2, Bowtie2, and R installed within your IDE.
    
		bool good_to_go = true;

		if (!this->bowtie2) {
      if (this->installed["SSAHA2"].size() == 0)
      {
        good_to_go = false;
        cerr << "---> ERROR Required executable \"ssaha2\" not found." << endl;
        cerr << "---> See http://www.sanger.ac.uk/resources/software/ssaha2" << endl;
      }
      else if (this->installed.count("SSAHA2_version_error_message")) {
        good_to_go = false;
        cerr << "---> ERROR Could not determine version of installed executable \"SSAHA2\"." << endl;
        cerr << "---> For found executable installed at [" << this->installed["SSAHA2"] << "]" << endl;
        cerr << "---> Commands \"SSAHA2 --version\" returned:" << endl;
        cerr << this->installed["SSAHA2_version_error_message"] << endl;
      }
    }
    
		if (this->bowtie2) 
    {
      if (this->installed["bowtie2"].size() == 0)
      {
        good_to_go = false;
        cerr << "---> ERROR Required executable \"bowtie2\" or \"bowtie2-build\" not found." << endl;
        cerr << "---> See http://bowtie-bio.sourceforge.net/bowtie2" << endl;
      }
      else if (this->installed.count("bowtie2_version_error_message")) {
        good_to_go = false;
        cerr << "---> ERROR Could not determine version of installed executable \"bowtie2\" or \"bowtie2-build\"." << endl;
        cerr << "---> For found executable installed at [" << this->installed["bowtie2"] << "]" << endl;
        cerr << "---> Commands \"bowtie2 --version\" or \"bowtie2-build --version\" returned:" << endl;
        cerr << this->installed["bowtie2_version_error_message"] << endl;
      }
      // version encoded in triplets of numbers
      else if (from_string<uint32_t>(this->installed["bowtie2_version"]) < 2000000007) {
        good_to_go = false;
        cerr << "---> ERROR Required executable \"bowtie2 version 2.0.0-beta7 or later\" not found." << endl;
        cerr << "---> Your version is " << this->installed["bowtie2_version_string"] << "." << endl;
        cerr << "---> See http://bowtie-bio.sourceforge.net/bowtie2" << endl;
      }
      else if ((from_string<uint32_t>(this->installed["bowtie2_version"]) == 2000003000)
           ||  (from_string<uint32_t>(this->installed["bowtie2_version"]) == 2000004000)) {
        good_to_go = false;
        cerr << "---> ERROR \"bowtie2 versions 2.0.3 and 2.0.4 are known to have bugs in" << endl;
        cerr << "---> ERROR \"SAM output that can cause breseq to crash. Please upgrade." << endl;
        cerr << "---> Your version is " << this->installed["bowtie2_version_string"] << "." << endl;
        cerr << "---> See http://bowtie-bio.sourceforge.net/bowtie2" << endl;
      }
      else {
        cout << "---> bowtie2  :: version " << this->installed["bowtie2_version_string"] << " [" << this->installed["bowtie2"] << "]" << endl;
      }
    }
		// R version 2.1 required
		if (this->installed["R"].size() == 0)
		{
			good_to_go = false;
			cerr << "---> ERROR Required executable \"R\" not found." << endl;
			cerr << "---> See http://www.r-project.org" << endl;
		}
    else if (this->installed.count("R_version_error_message")) {
      good_to_go = false;
      cerr << "---> ERROR Could not determine version of installed executable \"R\"." << endl;
      cerr << "---> For found executable installed at [" << this->installed["R"] << "]" << endl;
      cerr << "---> Command \"R --version\" returned:" << endl;
      cerr << this->installed["R_version_error_message"] << endl;
    }
    else if (from_string<uint32_t>(this->installed["R_version"]) < 2001000)
		{
      good_to_go = false;
      cerr << "---> ERROR Required executable \"R version 2.1.0 or later\" not found." << endl;
      cerr << "---> Your version is " << this->installed["R_version_string"] << "." << endl;
      cerr << "---> See http://www.r-project.org" << endl;
    }
    else {
      cout << "---> R        :: version " << this->installed["R_version_string"] << " [" << this->installed["R"] << "]" << endl;
    }

    if (this->installed["samtools"].size() == 0)
    {
      good_to_go = false;
      cerr << "---> ERROR Required executable \"samtools\" not found." << endl;
      cerr << "---> This should have been installed by the breseq installer." << endl;
    }
    else if (this->installed.count("samtools_version_error_message")) {
      good_to_go = false;
      cerr << "---> ERROR Could not determine version of installed executable \"samtools\"." << endl;
      cerr << "---> For found executable installed at [" << this->installed["samtools"] << "]" << endl;
      cerr << "---> Command \"samtools\" returned:" << endl;
      cerr << this->installed["samtools_version_error_message"] << endl;
    }
    else {
      cout << "---> samtools :: version " << this->installed["samtools_version_string"] << " [" << this->installed["samtools"] << "]" << endl;
    }

		if (!good_to_go) exit(0);
	}

	bool Settings::do_step(string done_key, string message)
	{
		string done_file_name = done_key;
		this->done_key_messages[done_key] = message;
    
		if (!file_exists(done_file_name.c_str()))
		{
      cerr << "+++   NOW PROCESSING " << message << endl;
      this->record_start_time(message);
      return true;
		}
    
		cerr << "--- ALREADY COMPLETE " << message << endl;

		ExecutionTime et;
    et.retrieve(done_file_name);    
    // @JEB Should check for errors, such as incomplete reads...

    this->execution_times.push_back(et);

		return false;
	}

	void Settings::done_step(string done_key)
	{
		string done_file_name = done_key;
		string message = this->done_key_messages[done_key];
		this->record_end_time(message);

		// create the done file with timing information
		this->execution_times.back().store(done_file_name);
	}
  
  
  void Settings::log(const string& message)
  {    
    create_path(this->output_path);
    ofstream LOG(this->log_file_name.c_str(), ios::app);
    ASSERT(!LOG.fail(), "Failed to open log file: " + this->log_file_name);
    LOG << message << endl;
    LOG.close();
  }



} // breseq namespace

