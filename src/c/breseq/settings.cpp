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

// Autoconf header
#include <config.h>

#include "libbreseq/settings.h"

#include "libbreseq/anyoption.h"


using namespace std;

namespace breseq
{

	void cReadFiles::Init(const vector<string>& read_file_names)
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
			pos = rf.m_base_name.rfind(".fastq");
			if ((pos != string::npos) && (pos = rf.m_base_name.size() - 6))
			{
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
    AnyOption options("Usage: breseq -r reference.gbk reads1.fastq [reads2.fastq, reads3.fastq...]");
    
    options
		("help,h", "produce advanced help message", TAKES_NO_ARGUMENT)
		// convert to basing everything off the main output path, so we don't have to set so many options
		("output,o", "path to breseq output [current path]", ".")
		("reference,r", "reference sequence in GenBank flatfile format (REQUIRED)")
    ("name,n", "human-readable name of sample/run for output [empty]", "")
    ("no-junction-prediction,j", "do not predict new sequence junctions", TAKES_NO_ARGUMENT)
    ("polymorphism-prediction,p", "predict polymorphic mutations", TAKES_NO_ARGUMENT)
    ("base-quality-cutoff,b", "ignore bases with quality scores lower than this value", 3)
    ("deletion-coverage-propagation-cutoff,u","value for coverage above which deletions are cutoff. 0 = calculated from coverage distribution", 0)
    ("deletion-coverage-seed-cutoff,s","value for coverage below which deletions are seeded", 0)
    ("require-complete-match", "only consider alignments that extend from end to end of a read", TAKES_NO_ARGUMENT)
    ("require-match-length", "only consider alignments that cover this many bases of a read", 0)
    ("require-match-fraction", "only consider alignments that cover this fraction of a read", 0.9)
    ("values-to-gd","",TAKES_NO_ARGUMENT, ADVANCED_OPTION) // @JEB @GRC added in for gathering/analyzing breseq values
    ("cnv","do experimental copy number variation prediction",TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("cnv-tile-size", "tile size for copy number variation prediction", 500, ADVANCED_OPTION)
    ("cnv-ignore-redundant", "only consider non-redundant coverage when using cnv", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
        
    ("verbose,v","",TAKES_NO_ARGUMENT)

    .processCommandArgs(argc, argv);
    
    options.addUsage("");
		options.addUsage("Utility Command Usage: breseq [command] options ...");
    options.addUsage("  Breseq Pipeline Commands: CONVERT_FASTQ, ERROR_COUNT");
    options.addUsage("  Breseq Post-Run Commands: BAM2ALN, BAM2COV");
    options.addUsage("  Genome Diff Commands: APPLY, COMPARE, INTERSECT, NOT_EVIDENCE, SUBTRACT, UNION");
    options.addUsage("  For help using a utility command, type: breseq [command] ");
    
    // make sure that the other csettings settonfig options are good:
    if (options.count("help"))
    {
      options.printAdvanced();
      exit(-1);
    }
    
    //! Settings: Global Workflow and Output
    
    this->base_output_path = options["output"];

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
    this->read_files.Init(read_file_names);
    
    // Reference sequence provided?
		if (options.count("reference") == 0)
		{
      options.addUsage("");
      options.addUsage("No reference sequences provided (-r).");
      options.printUsage();
      exit(-1);
		}
    this->reference_file_names = from_string<vector<string> >(options["reference"]);
    
    this->run_name = options["name"];
    
    this->no_junction_prediction = options.count("no-junction-prediction");
    this->do_copy_number_variation = options.count("cnv");
    this->copy_number_variation_tile_size = from_string<uint32_t>(options["cnv-tile-size"]);
    this->ignore_redundant_coverage = options.count("cnv-ignore-redundant");
    
    this->verbose = options.count("verbose");
    
    //! Settings: Read Alignment and Candidate Junction Read Alignment

    this->require_complete_match = (options.count("require-complete-match") > 0);
    this->require_match_length = from_string<uint32_t>(options["require-match-length"]);
    this->require_match_fraction = from_string<double>(options["require-match-fraction"]);

    //! Settings: Mutation Identification
    
    this->base_quality_cutoff = from_string<uint32_t>(options["base-quality-cutoff"]);

    this->deletion_coverage_propagation_cutoff = from_string<double>(options["deletion-coverage-propagation-cutoff"]);
    ASSERT(this->deletion_coverage_propagation_cutoff >= 0, "Argument --deletion-coverage-propagation-cutoff must be > 0")
    this->deletion_coverage_seed_cutoff = from_string<double>(options["deletion-coverage-seed-cutoff"]);
    ASSERT(this->deletion_coverage_propagation_cutoff >= 0, "Argument --deletion-coverage-seed-cutoff must be > 0")
    
    this->polymorphism_prediction = options.count("polymorphism-prediction");
    if (this->polymorphism_prediction) {
      
      // different default value
      this->polymorphism_frequency_cutoff = 0; // cut off if < X or > 1-X
      this->mixed_base_prediction = false;
    } 
 
    //! Settings: Experimental

    this->add_metadata_to_gd = options.count("values-to-gd");

    
		this->post_option_initialize();
    
    // Log the command line
    time_t stamp_time = time(NULL);
    this->log(ctime(&stamp_time));	
    this->log(this->full_command_line + "\n");
	}
  
  void Settings::command_line_run_header()
  {
    fprintf(stderr, "====================================================================================\n");
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
    fprintf(stderr, "Copyright (c) 2011      The University of Texas at Austin\n");
    fprintf(stderr, "====================================================================================\n");
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
    this->no_junction_prediction = false;
		this->no_mutation_prediction = false;
		this->no_deletion_prediction = false;
    this->no_alignment_generation = false;
		this->do_copy_number_variation = false;
    
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
    this->require_complete_match = false;
    this->require_match_length = 0;         
    this->require_match_fraction = 0.9;
    this->maximum_read_mismatches = -1;

    this->smalt = false;
    this->max_smalt_diff = 0;

    //! Settings: Candidate Junction Prediction
		this->preprocess_junction_min_indel_split_length = 3;
		this->required_both_unique_length_per_side = 5;
		this->required_one_unique_length_per_side = this->ssaha2_seed_length;  
    this->minimum_candidate_junction_pos_hash_score = 2;
    this->maximum_junction_sequence_insertion_length = 20;
    this->maximum_junction_sequence_overlap_length = 20;
    this->minimum_candidate_junctions = 10;
		this->maximum_candidate_junctions = 5000;
		this->maximum_candidate_junction_length_factor = 0.1;

    //! Settings: Alignment Resolution
    this->add_split_junction_sides = true;
    this->junction_pos_hash_neg_log10_p_value_cutoff = 3;

    //! Settings: Mutation Identification
    this->base_quality_cutoff = 3;
    this->deletion_coverage_propagation_cutoff = 0;
    this->deletion_coverage_seed_cutoff = 0;
    this->polymorphism_prediction = false;
    this->mixed_base_prediction = true;
    
		this->mutation_log10_e_value_cutoff = 10;
    this->polymorphism_log10_e_value_cutoff = this->mutation_log10_e_value_cutoff;
		this->polymorphism_bias_p_value_cutoff = 0.05;
		this->polymorphism_frequency_cutoff = 0.1;
		this->polymorphism_coverage_both_strands = 0;
    this->polymorphism_reject_homopolymer_length = 0;
		this->no_indel_polymorphisms = false;
    
    //! Settings: Output
    this->maximum_reads_to_align = 100;
		this->max_rejected_polymorphisms_to_show = 20;
		this->max_rejected_junctions_to_show = 10;
		this->hide_circular_genome_junctions = true;
    
		this->lenski_format = false;
		this->no_evidence = false;
		this->shade_frequencies = false;
    this->add_metadata_to_gd = false;
    
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

		this->error_rates_file_name = this->error_calibration_path + "/#.error_rates.tab";
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
		this->complete_mutations_text_file_name = this->mutation_identification_path + "/@.mutations.tab";
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
    this->ranges_text_file_name = this->copy_number_variation_path + "/@.ranges.tab";
    this->smoothed_ranges_text_file_name = this->copy_number_variation_path + "/@.smoothed_ranges.tab";
    this->copy_number_variation_cn_genome_diff_file_name = this->copy_number_variation_path + "/@.cn_evidence.gd";
    
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
    
    // SAMtools executables - look in the local bin path only
    string test_command = "which " + this->bin_path + "/samtools";
		this->installed["samtools"] = SYSTEM_CAPTURE(test_command, true);
    
		// search first for ssaha2 in the same location as breseq    
    test_command = "which " + this->bin_path + "/ssaha2";
		this->installed["SSAHA2"] = SYSTEM_CAPTURE(test_command, true);
    
    // attempt to fall back on system-wide install
    if (this->installed["SSAHA2"].size() == 0)
      this->installed["SSAHA2"] = SYSTEM_CAPTURE("which ssaha2", true);

    /*
		// check for default names
		this->installed["smalt"] = system("which smalt &>/dev/null") ? "smalt" : "";
		if (this->installed["smalt"].size() == 0)
		{
			this->installed["smalt"] = system("which smalt_i386 &>/dev/null") ? "smalt_i386" : "";
		}
		if (this->installed["smalt"].size() == 0)
		{
			this->installed["smalt"] = system("which smalt_ia64 &>/dev/null") ? "smalt_ia64" : "";
		}
		if (this->installed["smalt"].size() == 0)
		{
			this->installed["smalt"] = system("which smalt_x86_64 &>/dev/null") ? "smalt_x86_64" : "";
		}
		if (this->installed["smalt"].size() == 0)
		{
			this->installed["smalt"] = system("which smalt_MacOSX_i386 &>/dev/null") ? "smalt_MacOSX_i386" : "";
		}
    */
    
		this->installed["R"] = SYSTEM_CAPTURE("which R", true).size() ? "R" : "";
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
        vector<string> split_version_string = split(version_string, ".");
        if (split_version_string.size() == 3)
        {
          uint32_t R_numerical_version = from_string<uint32_t>(split_version_string[0]) * 1000000 + from_string<uint32_t>(split_version_string[1]) * 1000 + from_string<uint32_t>(split_version_string[2]);
          this->installed["R_version"] = to_string(R_numerical_version);
        }
      }
		}

    
	}

	void Settings::check_installed()
	{
    // Developer's Note
    //
    // If you are running things through a debugger (like in XCode), your $PATH may need to be
    // set to include the paths where you have SSAHA2 and R installed within your IDE.
    
		bool good_to_go = true;

		if (!this->smalt && this->installed["SSAHA2"].size() == 0)
		{
			good_to_go = false;
			cerr << "---> ERROR Required executable \"ssaha2\" not found." << endl;
			cerr << "---> See http://www.sanger.ac.uk/resources/software/ssaha2" << endl;
		}

		if (this->smalt && this->installed["smalt"].size() == 0)
		{
			good_to_go = false;
			cerr << "---> ERROR Required executable \"smalt\" not found." << endl;
			cerr << "---> See http://www.sanger.ac.uk/resources/software/smalt/" << endl;
		}

		// R version 2.1 required
		if (this->installed["R"].size() == 0)
		{
			good_to_go = false;
			cerr << "---> ERROR Required executable \"R\" not found." << endl;
			cerr << "---> See http://www.r-project.org" << endl;
		}
    else if (from_string<uint32_t>(this->installed["R_version"]) < 2001000)
		{
      good_to_go = false;
      uint32_t R_numerical_version = from_string<uint32_t>(this->installed["R_version"]);
      string R_version = to_string<uint32_t>(static_cast<uint32_t>(floor(R_numerical_version/1000000))) 
        + "." + to_string<uint32_t>(static_cast<uint32_t>(floor(R_numerical_version%1000000/1000)))
        + "." + to_string<uint32_t>(static_cast<uint32_t>(floor(R_numerical_version%1000)));

      cerr << "---> ERROR Required executable \"R version 2.1.0 or later\" not found." << endl;
      cerr << "---> Your version is " << R_version << endl;
      cerr << "---> See http://www.r-project.org" << endl;
    }

    if (this->installed["samtools"].size() == 0)
    {
      good_to_go = false;
      cerr << "---> ERROR Required executable \"samtools\" not found." << endl;
      cerr << "---> This should have been installed by the breseq installer." << endl;
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

