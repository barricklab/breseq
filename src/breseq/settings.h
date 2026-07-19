/*****************************************************************************

 AUTHORS

   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors

 LICENSE AND COPYRIGHT

   Copyright (c) 2008-2010 Michigan State University
   Copyright (c) 2011-2025 The University of Texas at Austin
   Copyright (c) 2025-     Michigan State University

   breseq is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2, or (at your option) any later version.

   SPDX-License-Identifier: GPL-2.0-or-later

*****************************************************************************/


#ifndef _BRESEQ_SETTINGS_H_
#define _BRESEQ_SETTINGS_H_

#include "common.h"
#include "storable.h"
#include "reference_sequence.h"
#include "ctpl_stl.h"

namespace breseq
{
  // constants used more places than just settings
  extern const int32_t kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation;
  extern const int32_t kBreseq_large_mutation_size_cutoff;
  extern const double kBreseq_ignore_within_this_multiple_of_average_read_length_of_contig_end;
  extern const char* kBreseqAlignmentScoreBAMTag;
  extern const char* kBreseqBestAlignmentScoreBAMTag;

	class ExecutionTime : public Storable {
  public:
		string _message;
		string _name;
		time_t _time;
		string _formatted_time;
		double _time_elapsed; // in seconds
		string _formatted_time_elapsed;
		time_t _time_start;
		string _formatted_time_start;
		time_t _time_end;
		string _formatted_time_end;
    // stores a list of files that need to be cleaned up after later steps
    storable_map<string,storable_vector<string> > _done_key_intermediate_files;
    
    virtual ~ExecutionTime() {};
    
    void serialize(ofstream& f)
    {
      write_to_file(f, _message);
      write_to_file(f, _name);
      write_to_file(f, _time);
      write_to_file(f, _formatted_time);
      write_to_file(f, _time_elapsed);
      write_to_file(f, _formatted_time_elapsed);
      write_to_file(f, _time_start);
      write_to_file(f, _formatted_time_start);
      write_to_file(f, _time_end);
      write_to_file(f, _formatted_time_end);
      _done_key_intermediate_files.serialize(f);
    }
    void deserialize(ifstream& f)
    {
      read_from_file(f, _message);
      read_from_file(f, _name);
      read_from_file(f, _time);
      read_from_file(f, _formatted_time);
      read_from_file(f, _time_elapsed);
      read_from_file(f, _formatted_time_elapsed);
      read_from_file(f, _time_start);
      read_from_file(f, _formatted_time_start);
      read_from_file(f, _time_end);
      read_from_file(f, _formatted_time_end);
      _done_key_intermediate_files.deserialize(f);
    }
	};
	
	// We need to be able to group read files for two reasons
	// 1) They may be paired-end, so we want to map them together
	// 2) They may have the same error rates, so we want to treat them together for error analysis

	struct cReadFile
	{
	public:
		string m_original_file_name;  // the original name provided at the command line
		string m_base_name;           // the original name minus path and .fastq ending (if any)
    string m_converted_file_name; // the name of the converted FASTQ file (if it exists)
		uint32_t m_paired_end_group;  // indicates what file contains paired reads
		uint32_t m_error_group;       // indicates what other read files have the same error rates
		uint32_t m_id;                // index used to refer to this fastq file in BAM

    cReadFile()
    {
      m_paired_end_group = UINT_MAX;
      m_error_group = UINT_MAX;
      m_id = UINT_MAX;
    }

    string file_name()
    {
      if (m_converted_file_name != "") return m_converted_file_name;
      return m_original_file_name;
    }

    string base_name() const { return m_base_name; }
	};

  // A set of one (unpaired) or two (paired) read files that are mapped together.
  // Two files are paired when their base names are identical length and differ in
  // exactly one character position where one has '1' and the other has '2'.
  // The collective name substitutes 'X' at that differing position.
  struct cReadFileSet
  {
  public:
    vector<cReadFile> m_files;  // 1 file (unpaired) or 2 files (paired, R1 at [0], R2 at [1])
    string m_base_name;         // collective name (X-substituted for paired, or individual base name)

    bool is_paired() const { return m_files.size() == 2; }
  };

	class cReadFileSets : public vector<cReadFileSet>
	{
	public:
    map<string,string> read_file_to_fastq_file_name_map;
    map<string,string> read_file_to_converted_fastq_file_name_map;

    cReadFileSets() { };
    cReadFileSets(const vector<string>& read_file_names, bool sam_files, bool allow_paired = false) { Init(read_file_names, sam_files, allow_paired); };
    ~cReadFileSets() { };

    void Init(const vector<string>& read_file_names, bool sam_files, bool allow_paired = false);
    string base_name_to_read_file_name(const string& base_name);

    // Returns base names of all individual files (flat, across all sets).
    vector<string> base_names() const
    {
      vector<string> return_value;
      for (const auto& rfs : *this)
        for (const auto& rf : rfs.m_files)
          return_value.push_back(rf.base_name());
      return return_value;
    }

    // Returns all individual cReadFile objects in a flat list (preserves set order, R1 before R2).
    vector<cReadFile> flat_files() const
    {
      vector<cReadFile> return_value;
      for (const auto& rfs : *this)
        for (const auto& rf : rfs.m_files)
          return_value.push_back(rf);
      return return_value;
    }
	};

  //
  // cReferenceSequenceSettings
  //
  // For grouping reference sequences and
  // flagging some for different treatment
  
  class cReferenceSequenceSettings
  {    
  public:
  
    // Saved directly from command line
    vector<string> m_reference_file_names;         
    vector<string> m_contig_reference_file_names;
    vector<string> m_junction_only_reference_file_names;
        
    // For lookup and iterating certain sets
    set<string> m_junction_only_seq_id_set;     // All seq_ids in junction_only reference files
    set<string> m_contig_seq_id_set;            // All seq_ids in the contig reference files
    set<string> m_call_mutations_seq_id_set;    // Has all NOT in junction_only, mutations called
    
    // For traversing groups that should have their coverage assigned together
    vector< vector<string> > m_seq_ids_by_coverage_group;
    // For looking up the coverage group id of a reference sequence
    map<string, uint32_t> m_seq_id_to_coverage_group_map;
  
    cReferenceSequenceSettings() {};
  
    cReferenceSequenceSettings(
                            cReferenceSequences& ref_seq_info,
                            vector<string>& reference_file_names,
                            vector<string>& contig_reference_file_names,
                            vector<string>& junction_only_reference_file_names  
                            );
  
  };


  //
  // Settings
  //
  // Main global command-line settings object

	class Settings
	{
    
    static string output_divider;
    static string global_bin_path;
    static string global_program_data_path;
    
  public:
    
    //! Multithreading
    static ctpl::thread_pool pool;
    static void sync_threads() { pool.sync(true); }
    
    static std::mutex lock;
    
    ////////////////////
    //! Data
    ////////////////////
    
    string byline;
		string website;
    string full_command_line;
		string arguments;
    
    //! Done file tracking
    string current_step_done_key;
    map<string,string> done_key_messages;
		vector<ExecutionTime> execution_times;
    storable_map<string,storable_vector<string> > done_key_intermediate_files;
    
    //! Read files
    cReadFileSets read_file_sets;
    
    ////////////////////
    //! Settings
    ////////////////////
    
    //! Settings: Global Workflow and Output
    string base_output_path;              // Default = cwd COMMAND-LINE OPTION
    
    //! Options that control which parts of the pipeline to execute
    string custom_run_name;                 // Default = <none> COMMAND-LINE OPTION
    string print_custom_run_name;           // custom_run_name with '_' replaced by ' '
    string genbank_field_for_seq_id;        // Header fields to use for seq_id Default = VERSION, other valid values are LOCUS and ACCESSION
    bool skip_read_filtering;               // Default = false
    bool skip_new_junction_prediction;      // Default = false COMMAND-LINE OPTION
    bool skip_read_alignment_and_missing_coverage_prediction;          // Default = false
    bool skip_missing_coverage_prediction;  // Default = false set to true if targeted_sequencing
    bool predict_copy_number;                    // Default = false COMMAND-LINE OPTION
    bool do_periodicity;                    // Default = false COMMAND-LINE OPTION
    
    //! Settings: Read File Options
    vector<string> read_file_names;             // REQUIRED COMMAND-LINE OPTION
    bool aligned_sam_mode;                      // Default = false COMMAND-LINE OPTION
    double  read_file_coverage_fold_limit;      // Default = 0 (OFF) COMMAND-LINE OPTION
    uint32_t read_file_read_length_min;         // Default = 18 COMMAND-LINE OPTION
    double read_file_max_same_base_fraction;    // Default = 0.9 COMMAND-LINE OPTION
    double read_file_max_N_fraction;            // Default = 0.5 COMMAND-LINE OPTION
    uint32_t read_file_long_read_trigger_length;// Default = 1000 COMMAND-LINE OPTION
    uint32_t read_file_long_read_split_length;  // Default = 200 COMMAND-LINE OPTION
    bool read_file_long_read_distribute_remainder;           // Default = false COMMAND-LINE OPTION
    bool paired_mapping;                                     // Default = true COMMAND-LINE OPTION (disable with --no-paired-mapping)

    // Reference sequences
    vector<string> all_reference_file_names;    // REQUIRED COMMAND-LINE OPTION (filled by below)
    vector<string> normal_reference_file_names; // Default = EMPTY COMMAND-LINE OPTION
    vector<string> contig_reference_file_names; // Default = EMPTY COMMAND-LINE OPTION
    vector<string> junction_only_file_names;    // Default = EMPTY COMMAND-LINE OPTION
    cReferenceSequenceSettings refseq_settings; // This is extra settings data initialized from
                                                // refseqs and how they are provided to the CLI

    //! DEBUG options
    
    // verbose level
    uint32_t verbose;                         // Default = 0 (OFF) COMMAND-LINE OPTION
		uint32_t alignment_read_limit;            // Default = 0 (OFF)
		uint32_t candidate_junction_read_limit;   // Default = 0 (OFF)
    uint32_t resolve_alignment_read_limit;    // Default = 0 (OFF)
    bool     alignment_mask_ref_matches;      // Default false

    //! Don't delete intermediate files
    bool keep_all_intermediates;              // Default = false

    //! Settings: Read Alignment and Candidate Junction Read Alignment
    
    int32_t num_processors;       // Default = 1

    // bowtie2 options. Each of bowtie2_stage1/stage2/junction is a complete, independent
    // command-line string -- match/mismatch/gap scoring, alignment mode, and mapping cutoffs
    // are all bundled together per stage rather than sharing one scoring string across all
    // three, since stage 2/junction can use different scoring/alignment modes than stage 1
    // (see end_to_end below).
    string bowtie2_stage1;          // COMMAND-LINE OPTION
    string bowtie2_stage2;          // COMMAND-LINE OPTION (can be blank to skip step)
    string bowtie2_junction;        // COMMAND-LINE OPTION
    //! Use --end-to-end (instead of --local) for the stage 1 bowtie2 alignment and the
    //! candidate-junction realignment pass. This is now the DEFAULT (Default = true); the
    //! --local-mapping command-line option sets it false to restore --local behavior. When
    //! true, swaps in complete, independent defaults for bowtie2_stage1/bowtie2_junction
    //! (translated to --end-to-end's scoring scale, including --ma 0 -- bowtie2 rejects a
    //! nonzero --ma whenever --score-min can go negative) rather than patching the
    //! --local-oriented pristine defaults at the point bowtie2 is actually invoked.
    //! bowtie2_stage2 is unaffected (always --local, --ma 1). An explicit
    //! --bowtie2-stage1/--bowtie2-junction always wins over the corresponding default.
    bool     end_to_end;            // COMMAND-LINE OPTION (disabled by --local-mapping)
    uint64_t bowtie2_junction_maximum_alignments_to_consider_per_read;       // Default = 2000
    uint64_t bowtie2_genome_maximum_alignments_to_consider_per_read;         // Default = 2000

    //! reads are never included in the BAM alignment file if they fail these guards
    uint32_t minimum_mapping_quality;  // COMMAND-LINE OPTION
		uint32_t require_match_length;     // Default = 0 (OFF) COMMAND-LINE OPTION
		double   require_match_fraction;   // Default = 0.9     COMMAND-LINE OPTION
    //! ignore reads with this many or more mismatches (I+D+MM)
    int32_t  maximum_read_mismatches;     // Default = -1 (OFF)
    
    //! Settings: Candidate Junction Prediction

    //! Threshold (in bases) for splitting an indel found during candidate-junction
    //! preprocessing: an indel already present in a single alignment's own CIGAR that is this
    //! length or longer is split into separate junction-candidate-supporting alignment
    //! records (JC evidence), unless it is entirely a length-change to a reference homopolymer
    //! that was already this long or longer (in which case it is left as RA evidence).
    //! Default = 5.
    int32_t  junction_indel_split_length;
		int32_t required_both_unique_length_per_side;           // Set = junction_minimum_side_match
    double   required_both_unique_length_per_side_fraction; // Default = 0.2 
		int32_t required_one_unique_length_per_side;            // Default = 0 (OFF)
    uint32_t unmatched_end_minimum_read_length;             // Default = 50
    double   unmatched_end_length_factor;                   // Default = 0.0
    
    // The number of mismatches allowed at the end of a read is:
    //   (read_length - unmatched_end_minimum_read_length) * unmatched_end_length_factor 
    
    uint32_t required_junction_read_end_min_coordinate(uint32_t read_length) const {
      int32_t unmatched_end_max_length = static_cast<int32_t>(floor(
        static_cast<double>(static_cast<int32_t>(read_length) - static_cast<int32_t>(unmatched_end_minimum_read_length)) 
        * unmatched_end_length_factor
      ));
      if (unmatched_end_max_length <= 0) 
        return read_length;
      else
        return read_length - unmatched_end_max_length;
    }
    
    bool valid_junction_read_bounds(pair<uint32_t,uint32_t> _read_bounds, uint32_t _read_length) const {      
      return (_read_bounds.first == 1) && (_read_bounds.second >= required_junction_read_end_min_coordinate(_read_length) );
    }
    
    uint32_t maximum_junction_sequence_insertion_length;  // Default = 20
    uint32_t maximum_junction_sequence_overlap_length;    // Default = 20
    
    // The maximum overlap (negative or positive) is 
    //   (read_length - unmatched_end_minimum_read_length) * unmatched_end_length_factor 
    double maximum_junction_sequence_negative_overlap_length_fraction;    // Default = 0.1
    uint32_t maximum_junction_sequence_negative_overlap_length_minimum;   // Default = 12

    double maximum_junction_sequence_positive_overlap_length_fraction;    // Default = 0.4
    uint32_t maximum_junction_sequence_positive_overlap_length_minimum;   // Default = 0
    
    uint32_t highly_redundant_junction_ignore_passed_pair_limit;          // Default = 200
    uint64_t maximum_junction_sequence_passed_alignment_pairs_to_consider;// Default = 100000, 0 = DO NOT LIMIT
    
    // Which candidate junctions do we test?
		uint32_t minimum_candidate_junction_pos_hash_score;   // Default = 2
		uint32_t concordant_pairs_to_make_unique;             // Default = 3 (--paired-mapping side disambiguation)
		uint32_t minimum_candidate_junctions;                 // Default = 10
		uint32_t maximum_candidate_junctions;                 // Default = 5000, 0 = DO NOT LIMIT
    double maximum_candidate_junction_length_factor;      // Default = 0.1, 0 = DO NOT LIMIT
    
    bool penalize_negative_junction_overlap;              // Manually set. True for experimental treatment.
      
    //! Settings: Alignment Resolution
		bool add_split_junction_sides;                        // Default = true (possibly remove this option) 
    uint32_t minimum_alignment_resolution_pos_hash_score; // Default = 2
    double minimum_pr_no_read_start_per_position;         // Default = 0.1
    int32_t junction_minimum_side_match;                  // Default = 1 or 6;
    double junction_pos_hash_neg_log10_p_value_cutoff;    // Default = 2.0, 0 = means don't calculate (to implement)
    bool junction_allow_suboptimal_matches;                       // Default = false
    
    
    //! Settings: Mutation Identification
    
    string user_evidence_genome_diff_file_name; // Default = none COMMAND-LINE OPTION
    
    //! ignore bases below this cutoff for RA evidence (still counted for deletions?)
    uint32_t base_quality_cutoff;                         // Default 3    COMMAND-LINE OPTION
    uint32_t quality_score_trim;                          // Default OFF  COMMAND-LINE OPTION

    //! treated as read numbers if integers >= 1.0 and percentages of average coverage if > 0.0 and < 1.0
    double deletion_coverage_propagation_cutoff;          // Default = calculated COMMAND-LINE OPTION
    double deletion_coverage_seed_cutoff;                 // Default = 0;         COMMAND-LINE OPTION
    
    bool call_mutations_overlapping_missing_coverage;     // Default = false      COMMAND-LINE OPTION

    //! Settings: Discordant Pair (DP) Evidence
    int32_t discordant_pair_seed;                         // Default = 3;         COMMAND-LINE OPTION
    double discordant_pair_skew_cutoff;                   // Default = 3.0;       COMMAND-LINE OPTION

    //! Settings: Soft-Clipping Evidence
    bool     predict_soft_clipping;                       // Default = false     COMMAND-LINE OPTION
    uint32_t soft_clipping_minimum_bases;                 // Default = 12        COMMAND-LINE OPTION
    double   soft_clipping_log10_e_value_cutoff;          // Default = 3.0       COMMAND-LINE OPTION
    
    //! These are mutually exclusive settings (polymorphism prediction overrides mixed_base_prediction)
    bool polymorphism_prediction;                         // Default = false COMMAND-LINE OPTION
    //! Predict not only consensus genotype calls, but test mixed states between them.
    bool mixed_base_prediction;                           // Default = true
    
    //! References are amplicons or ultra-deep sequencing of a small reference
    bool targeted_sequencing;                             // Default = false COMMAND-LINE OPTION
    
    //! Verbose output of bases encountered at each position
    bool print_mutation_identification_per_position_file;

    double mutation_log10_e_value_cutoff;                         // Default = 10
    uint32_t consensus_minimum_variant_coverage;                  // Default = 0
    uint32_t consensus_minimum_total_coverage;                    // Default = 0
    uint32_t consensus_minimum_variant_coverage_each_strand;      // Default = 0
    uint32_t consensus_minimum_total_coverage_each_strand;        // Default = 0
    uint32_t consensus_reject_indel_homopolymer_length;        // Default = 0 (OFF)
    uint32_t consensus_reject_surrounding_homopolymer_length;  // Default = 0 (OFF)

    double consensus_frequency_cutoff;                            // Default = 0.8
    
    double polymorphism_log10_e_value_cutoff;                     // Default = mutation_log10_e_value_cutoff = 10
		double polymorphism_bias_p_value_cutoff;                      // Default = 0.05 for mixed base | 0 (OFF) for polymorphism
		double polymorphism_frequency_cutoff;                         // Default = 0.1 for mixed base | 0.0 for polymorphism
    uint32_t polymorphism_minimum_variant_coverage;               // Default = 0
    uint32_t polymorphism_minimum_total_coverage;                 // Default = 0
    uint32_t polymorphism_minimum_variant_coverage_each_strand;   // Default = 0
		uint32_t polymorphism_minimum_total_coverage_each_strand;     // Default = 0 for mixed base | 2 for polymorphism
		uint32_t polymorphism_reject_indel_homopolymer_length;        // Default = 0 (OFF)
    uint32_t polymorphism_reject_surrounding_homopolymer_length;  // Default = 0 (OFF)
		bool polymorphism_no_indels;                                  // Default = false
    double polymorphism_precision_decimal;                        // Default = not used for mixed base | 0.0000000001 for polymorphism
    uint32_t polymorphism_precision_places;                       // Default = 3 for mixed base | 10 for polymorphism

		
		//! Settings: Copy Number Variation
    uint32_t periodicity_method;
    uint32_t periodicity_start;
    uint32_t periodicity_end;
    uint32_t periodicity_step;
    
    //! Settings: Mutation Prediction
    int32_t size_cutoff_AMP_becomes_INS_DEL_mutation;
      // Default = kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation
    double ignore_within_this_multiple_of_average_read_length_of_contig_end;
      // Default = kBreseq_ignore_within_this_multiple_of_average_read_length_of_contig_end
    
    //! Settings: Output
    //!
    bool output_unmapped_reads;                             // Default = false COMMAND-LINE OPTION
    bool no_evidence_html;                                  // Default = false COMMAND-LINE OPTION
    bool zip_html;                                          // Default = true COMMAND-LINE OPTION (disable with --unzipped-html)
    
    //! Don't include javascript in HTML output (provided for Galaxy integration)
    bool no_javascript;                                     // Default = false COMMAND-LINE OPTION

    uint32_t max_flanking_columns;                          // Default = 100   COMMAND-LINE OPTION
    uint32_t max_displayed_reads;                           // Default = 100   COMMAND-LINE OPTION

    bool no_list_js;                                        // Default = false (use list.js for filtering tables)

    string header_genome_diff_file_name;                    // Default = NONE  COMMAND-LINE OPTION
    string mask_genome_diff_file_name;                      // Default = ""    COMMAND-LINE OPTION
    string mask_mode;                                       // Default = "NONE" COMMAND-LINE OPTION
    uint32_t max_nucleotides_to_show_in_tables;      // Default = 8
    uint32_t max_rejected_read_alignment_evidence_to_show;  // Default = 20
		uint32_t max_rejected_junction_evidence_to_show;        // Default = 10
		bool hide_circular_genome_junctions;                    // Default = true
    //! special output for Blount paper - not implemented in C++!
		bool lenski_format;                                     // Default = false (not a public option!)
    //! don't create any HTML evidence files
		bool no_evidence;                                       // Default = false (rarely used)

    
    //! Settings: Experimental
    
    //! @GRC added in for gathering/analyzing breseq values
    bool add_metadata_to_gd;                              // Default = false COMMAND-LINE OPTION
    bool junction_debug;                                  // Default = false COMMAND-LINE OPTION


    ////////////////////
    //! File Paths
    ////////////////////
    
    //! Paths populated from location of executable
		string bin_path;                  // absolute path to where this binary resides
    string program_data_path;         // path to where R scripts and images reside
		map<string, string> installed;    // hash of paths to programs that we call
    
    //! Paths: Sequence conversion
    string sequence_conversion_path;
    string sequence_conversion_done_file_name;

		string converted_fastq_file_name;
		string unwanted_fasta_file_name;
    string reference_trim_file_name;
		string sequence_conversion_summary_file_name;
    
    //! Paths: Read alignment
		string reference_alignment_path;
    string reference_alignment_done_file_name;
    
		string reference_hash_file_name;

    // Staged alignment to final reference SAM file
    string stage1_reference_sam_file_name;
    string stage1_unmapped_reads_fastq_file_name;
    string stage2_reference_sam_file_name;
    string reference_sam_file_name;
    
    //! Paths: Junction Prediction
    string candidate_junction_path;
    
    string preprocess_junction_done_file_name;
    string preprocess_junction_best_sam_file_name;
		string preprocess_junction_split_sam_file_name;
		string preprocess_junction_split_pair_positions_file_name;

    string paired_mapping_distance_done_file_name;
    string paired_mapping_distance_distribution_file_name;
    string paired_mapping_distance_plot_file_name;
    string paired_mapping_distance_summary_file_name;
    string discordant_pairs_file_name;
    string discordant_pairs_plot_file_name;

    string candidate_junction_done_file_name;
		string coverage_junction_best_bam_unsorted_file_name;
		string coverage_junction_best_bam_file_name;
		string coverage_junction_distribution_file_name;
		string coverage_junction_plot_file_name;
		string coverage_junction_summary_file_name;
    string coverage_junction_error_count_summary_file_name;

    string coverage_junction_done_file_name;
    string candidate_junction_fasta_file_name;
    string candidate_junction_detailed_file_name;
    string candidate_junction_faidx_file_name;
    string candidate_junction_trim_file_name;
		string candidate_junction_summary_file_name;
    
    //! Paths: Junction Alignment
    string candidate_junction_alignment_path;
    string candidate_junction_alignment_done_file_name;
    
    string candidate_junction_hash_file_name;
		string candidate_junction_sam_file_name;
    
    //! Paths: Alignment Resolution
    string alignment_resolution_path;
		string alignment_correction_done_file_name;
    
		string resolved_reference_sam_file_name;
		string resolved_junction_sam_file_name;
		string alignment_resolution_summary_file_name;
		string resolved_paired_mapping_distance_summary_file_name;
		// Per-seq_id interior concordant-pair crossing histogram (crossing<TAB>count), written in resolve
		// as a CSV/tab intermediate; read by the DP score and the crossing plots (# = seq_id).
		string concordant_pair_crossing_distribution_file_name;
    string jc_genome_diff_file_name;
    
    string junction_debug_file_name;
    
    //! Paths: BAM conversion
		string bam_path;
    string bam_done_file_name;

		string reference_bam_unsorted_file_name;
		string junction_bam_unsorted_file_name;
		string junction_bam_file_name;

		//! Paths: Error Calibration
		string error_calibration_path;
		string error_counts_done_file_name;
		string error_rates_done_file_name;
    
		string error_counts_file_name;
		string error_rates_file_name;
		string coverage_file_name;
		string unique_only_coverage_distribution_file_name;
		string error_rates_summary_file_name;
		string error_rates_base_qual_error_prob_file_name;
		string soft_clipping_counts_file_name;
		string soft_clipping_summary_file_name;

		//! Paths: Mutation Identification
		string mutation_identification_path;
    
		string mutation_identification_done_file_name;    
		string mutation_identification_per_position_file_name;
		string complete_coverage_text_file_name;
    string ra_mc_genome_diff_file_name;
    string dp_candidate_regions_file_name;
    string discordant_pair_done_file_name;
    string dp_genome_diff_file_name;


		//! Paths: Copy Number Variation
    string copy_number_variation_path;
    string copy_number_variation_done_file_name;
    string copy_number_evidence_genome_diff_file_name;
    
    string periodicity_table_file_name;
    string periodicity_done_file_name;
    
    //! Paths: Output
    string output_path;
		string output_done_file_name;
    // subsets of output
    string output_mutation_prediction_done_file_name;
    string output_mutation_annotation_done_file_name;
    
		string log_file_name;
		string index_html_file_name;
		string summary_html_file_name;
		string marginal_html_file_name;

    string local_evidence_path;
		string evidence_path;
    string local_html_archive_file_name;
    string html_archive_file_name;
    string local_html_evidence_file_name;
    string html_evidence_file_name;
    string evidence_genome_diff_file_name;
    string final_genome_diff_file_name;
    string preannotated_genome_diff_file_name;
    
		string local_coverage_plot_path;
		string coverage_plot_path;
		string overview_coverage_plot_file_name;
    
		string output_calibration_path;
		string unique_only_coverage_plot_file_name;
		string concordant_pair_crossing_plot_file_name;       // single run-wide reference plot (no token)
		string concordant_pair_crossing_seq_plot_file_name;   // per-seq empirical-vs-projected overlay (# = seq_id)
		string error_rates_plot_file_name;
    
		string breseq_icon_graphic_from_file_name;
		string breseq_icon_graphic_to_file_name;
    
    //! Paths: Data
    string data_path;
    
		string reference_bam_file_name;
    string reference_fasta_file_name;
		string reference_faidx_file_name;
		string reference_gff3_file_name;
    string unmapped_reads_fastq_file_name;

    string data_vcf_file_name;
    string data_genome_diff_file_name;
    string data_annotated_genome_diff_file_name;
    string data_json_summary_file_name;
    
    //! Paths: Experimental
		string long_pairs_file_name;

    
    ////////////////////
    //! Methods
    ////////////////////
    
    // Set up defaults here -- this version does not require command line
		Settings(const string& _base_output_path = ".");
    
    // Constructor for default action
    Settings(int argc, char* argv[]);
    
    static void command_line_run_header();
    

    //! Call this at the very beginning of execution to set up program paths
    void static set_global_paths()
    {
      //debug
      //cout << "DATADIR:" << DATADIR << endl;
      
      global_bin_path = getExecPath();
      
      // 1st choice: use *absolute* DATADIR (set by XCode)
      if (string(DATADIR).substr(0, 1) == "/") {
        global_program_data_path = DATADIR;
        string test_file = global_program_data_path + "/breseq_icon.png";
        ASSERT(file_exists(test_file.c_str()), "Could not find expected data files inside of path set by development environment: " + global_program_data_path);
      }
      
      // 2nd choice: use BRESEQ_DATA_PATH environmental variable
      //             This is automatically set when running tests
      char * breseq_data_path;
      breseq_data_path = getenv ("BRESEQ_DATA_PATH");
      if (breseq_data_path!=NULL) {
        global_program_data_path = breseq_data_path;
        cerr << "Program data path set via BRESEQ_DATA_PATH: " << breseq_data_path << endl;
        
        string test_file = global_program_data_path + "/breseq_icon.png";
        ASSERT(file_exists(test_file.c_str()), "Could not find expected data files inside of data path set by environmental variable BRESEQ_DATA_PATH: " + global_program_data_path + "\nPlease correct this setting.");
      }
      
      // 3rd choice: Build it from the executable path, if we found one.
      //             DATADIR is a relative path (../share/breseq)
      //             This case is the default for binary distributions
      if (global_program_data_path.length() == 0) {
        ASSERT(global_bin_path.length() != 0, string("Failed to automatically detect the location of this executable. To continue, set the BRESEQ_DATA_PATH environmental variable to the 'share/") + PACKAGE_NAME + "' directory of your installation. For example, to '/path/to/executable/../share/'" + PACKAGE_NAME + "'");
        
        global_program_data_path = global_bin_path + "/" + DATADIR;
        
        // Get rid of any trailing slash in breseq_data_path
        if (global_program_data_path[global_program_data_path.length()-1] == '/') {
          global_program_data_path.erase(global_program_data_path.length()-1, 1);
        }
        
        string test_file = global_program_data_path + "/breseq_icon.png";
        ASSERT(file_exists(test_file.c_str()), "Could not find expected data files inside of data path set relative to executable: " + global_program_data_path + "\nPlease, see the installation instructions in the HTML documentation.");
      }
      
      //for debug
      //cout << "Global bin path: " + global_bin_path << endl;
      //cout << "Global program data path: " + global_program_data_path << endl;
    }

    //! Utility functions for getting paths
    static string get_bin_path()
    {
      return global_bin_path;
    }
    
    static string get_program_data_path()
    {
      return global_program_data_path;
    }
    
		//! Utility function to substitute specific details into a generic file name
		static string file_name(const string& _file_name_key, const string& _substitute, const string& _with)
		{
			string s(_file_name_key);
      if (_substitute.size())
        s = substitute(s, _substitute, _with);
			return s;
		}
    
    static string path_file_name(const string& _path, const string& _file_name_key, const string& _substitute, const string& _with)
		{
			string s(path_to_filename(_file_name_key));
      if (_substitute.size())
        s = substitute(s, _substitute, _with);

			return ((_path.size() > 0) ? _path + "/" + s : s);
		}
    
    static string path_file_name(const string& _path, const string& _file_name_key)
		{
			string s(path_to_filename(_file_name_key));
			return ((_path.size() > 0) ? _path + "/" + s : s);
		}
    
		//! Utility function to get relative path from two file locations for making HTML files
    static string relative_path(const string& full_path, const string& base_path) 
    {
      string s(full_path);
      if (s.substr(0, base_path.size()) == base_path)
      {
        s.erase(0,base_path.size());
      }
      
      if (s.at(0) == '/')
        s.erase(0,1);
      
      return s;
    }

		string ctool(string tool_name, bool allow_fail = false)
		{
      if (this->installed[tool_name].size() == 0)
      {
        cerr << "Executable '" << tool_name << "' not found in breseq bin path '" << this->bin_path << "'." << endl;
        if (!allow_fail) exit(-1);
        return "";
      }
			return this->installed[tool_name];
		}
    
    static string time2string(const time_t& _time)
    {
      const struct tm * time_info = localtime(&_time);
      char s[1024];
      strftime(s, 1024, "%H:%M:%S %d %b %Y", time_info);
      return s;
    }
    
    static string elapsedtime2string(double _diff_time)
    {
      stringstream ss;
      
      double t = _diff_time; // in seconds
      uint32_t tm_yday = static_cast<uint32_t>(floor( t / (60*60*24)));
      t -= tm_yday * 60*60*24;
      uint32_t tm_hour = static_cast<uint32_t>(floor( t / (60*60)));
      t -= tm_hour * 60*60;
      uint32_t tm_min = static_cast<uint32_t>(floor( t / (60)));
      t -= tm_min * 60;
      uint32_t tm_sec = static_cast<uint32_t>(floor( t / (1)));      
      if (tm_yday > 0)
      {
        if (ss.str().length() > 0) ss << " ";
        ss << tm_yday;
        ss << ((tm_yday!=1)?" days":" day");
      }

      if (tm_hour > 0)
      {
        if (ss.str().length() > 0) ss << " ";
        ss << tm_hour;
        ss << ((tm_hour!=1)?" hours":" hour");
      }
      
      if (tm_min > 0)
      {
        if (ss.str().length() > 0) ss << " ";
        ss << tm_min;
        ss << ((tm_min!=1)?" minutes":" minute");
      }
      
      //if (tm_sec > 0)
      {
        if (ss.str().length() > 0) ss << " ";
        ss << tm_sec;
        ss << ((tm_sec!=1)?" seconds":" second");
      }

      return ss.str();
    }

    
    void record_start_time(const string& message)
    {
      ExecutionTime ex_time;
      time_t this_time = time(NULL);
      ex_time._time_start = this_time;
      ex_time._formatted_time_start = time2string(this_time);
      ex_time._time_end = 0;
      ex_time._formatted_time_end = "";
      ex_time._time_elapsed = 0;
      ex_time._formatted_time_elapsed = "";
      ex_time._message = message;

      this->execution_times.push_back(ex_time);
    }

		void record_end_time(const string& message)
		{
			uint32_t i = 0;
			while (i < this->execution_times.size())
			{
				if (this->execution_times[i]._message == message)
					break;
				i++;
			}

			if (i >= this->execution_times.size())
			{
				cout << "Did not find matching start time for:" << endl << message << endl;
				ExecutionTime blank;
				this->execution_times.push_back(blank);
			}

			ExecutionTime& ex_time = this->execution_times[i];
			ex_time._message = message;

			time_t this_time = time(NULL);
			ex_time._time_end = this_time;
			ex_time._formatted_time_end = time2string(this_time);

			//if we had a previous time, calculate elapsed
			if (i < this->execution_times.size())
			{
				ex_time._time_elapsed = difftime(ex_time._time_end, ex_time._time_start);
				ex_time._formatted_time_elapsed = elapsedtime2string(ex_time._time_elapsed);
			}
		}
    
    // Convenience accessors for reference file settings.

    void init_reference_sets(cReferenceSequences& ref_seq_info) {
      refseq_settings = cReferenceSequenceSettings(ref_seq_info, normal_reference_file_names, contig_reference_file_names, junction_only_file_names);
    }
    
    bool is_junction_only_reference(string& seq_id) { return refseq_settings.m_junction_only_seq_id_set.count(seq_id); }
    set<string> call_mutations_seq_id_set() const { return refseq_settings.m_call_mutations_seq_id_set; }
    set<string> junction_only_seq_id_set() const { return refseq_settings.m_junction_only_seq_id_set; }
    vector< vector<string> > seq_ids_by_coverage_group() const { return refseq_settings.m_seq_ids_by_coverage_group; }
    uint32_t seq_id_to_coverage_group(const string& _seq_id) const { return refseq_settings.m_seq_id_to_coverage_group_map.find(_seq_id)->second; }

    // Convenience accessors for read file settings.
     
		string base_name_to_read_file_name(string base_name)
		{
      return this->read_file_sets.base_name_to_read_file_name(base_name);
		}

		void check_installed();
    
    bool do_step(const string& done_key, const string& message);
    void set_current_step_done_key(const string& done_key) { current_step_done_key = done_key; }
    string get_current_step_done_key() {return current_step_done_key; };
		void done_step(const string& done_key);
    void track_intermediate_file(const string& done_key, const string& file_path);
    
    void log(const string& message);
        
    ~Settings() { }
    
  private:

		void pre_option_initialize(int argc = 0, char* argv[] = NULL);
		void post_option_initialize();
    void init_installed();
    

  };

class UserOutput
{
  public:
    UserOutput(const string& command, size_t shift = 4)
      : _command(command)
      , _n_shift(shift)
      , _is_shifted(false)
      , _n_default(shift)
    {
      cerr << endl << "*** Begin " << _command << " ***" << endl;
    }
    
    UserOutput(size_t shift = 4)
      : _command("")
      , _n_shift(shift)
      , _is_shifted(false)
      , _n_default(shift)
    { }
    
    ~UserOutput()
    {
      if (_command.size()) {
        cerr << endl << "*** End " << _command << " ***" << endl;
      }
    }

    template<class T> UserOutput& operator<<(const T& value)
    {
      if (!_is_shifted) {
        cerr << string(_n_shift, ' ');
        _is_shifted = true;
      }
      cerr << value;
      return *this;
    }

    UserOutput& operator<<(ostream& (*op_ptr)(ostream&))
    {
      //TODO find way to op_ptr == &std::endl.
      if (_is_shifted) {
        _is_shifted = false; 
      }

      (*op_ptr)(cerr);
      return *this;
    }

    UserOutput& operator()(const string& step, const string& value = "")
    {
      this->end_step();

      cerr << endl << string(_n_shift, ' ') << step;
      if (value.size()) {
        cerr << ": " << value;
        this->end_step();
      } else {
        _n_shift += 4;
      } 
      cerr << endl;
      return *this;
    }

    template<class T> UserOutput& operator()(const string& step, const T& value)
    {
      this->end_step();

      cerr << string(_n_shift, ' ') << step << ": " << value << endl;
      this->end_step();

      return *this;
    }

    void end_step(void)
    {
      _n_shift = _n_default;
    }

  private:
    string        _command;
    size_t        _n_shift;
    bool          _is_shifted;
    size_t        _n_default;
};
  



   
} // breseq namespace

#endif

