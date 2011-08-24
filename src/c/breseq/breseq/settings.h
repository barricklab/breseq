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


#ifndef _BRESEQ_SETTINGS_H_
#define _BRESEQ_SETTINGS_H_

#include "common.h"
#include "storable.h"

#ifndef UNDEFINED
#define UNDEFINED UINT_MAX
#define is_defined(x) (x != UINT_MAX)
#endif

using namespace std;

namespace breseq
{

	class ExecutionTime : public Storable {
  public:
		string _message;
		string _name;
		time_t _time;
		string _formatted_time;
		time_t _time_elapsed;
		string _formatted_time_elapsed;
		time_t _time_start;
		string _formatted_time_start;
		time_t _time_end;
		string _formatted_time_end;
    
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
    }
	};

	class Coverage : public Storable
	{
  public:
		double junction_accept_score_cutoff;
		double deletion_coverage_propagation_cutoff;
		double junction_coverage_cutoff;
		double junction_keep_score_cutoff;
		double nbinom_size_parameter;
		double nbinom_mean_parameter;
		double nbinom_prob_parameter;
		double average;
		double variance;
		double dispersion;
    
    void serialize(ofstream& f)
    {
      write_to_file(f, junction_accept_score_cutoff);
      write_to_file(f, deletion_coverage_propagation_cutoff);
      write_to_file(f, junction_coverage_cutoff);
      write_to_file(f, junction_keep_score_cutoff);
      write_to_file(f, nbinom_size_parameter);
      write_to_file(f, nbinom_mean_parameter);
      write_to_file(f, nbinom_prob_parameter);
      write_to_file(f, average);
      write_to_file(f, variance);
      write_to_file(f, dispersion);
    }
    void deserialize(ifstream& f)
    {
      read_from_file(f, junction_accept_score_cutoff);
      read_from_file(f, deletion_coverage_propagation_cutoff);
      read_from_file(f, junction_coverage_cutoff);
      read_from_file(f, junction_keep_score_cutoff);
      read_from_file(f, nbinom_size_parameter);
      read_from_file(f, nbinom_mean_parameter);
      read_from_file(f, nbinom_prob_parameter);
      read_from_file(f, average);
      read_from_file(f, variance);
      read_from_file(f, dispersion);
    }
	};
	
	class cReferenceSequences;

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
      m_paired_end_group = UNDEFINED;
      m_error_group = UNDEFINED;
      m_id = UNDEFINED;
    }
    
    string file_name()
    {
      if (m_converted_file_name != "") return m_converted_file_name; 
      return m_original_file_name;
    }
    
    string base_name() { return m_base_name; }
	};

	typedef vector<vector<cReadFile> > cReadFileGroup; // unused? delete? @JEB

	class cReadFiles : public vector<cReadFile>
	{
	public:
    map<string,string> read_file_to_fastq_file_name_map;
    map<string,string> read_file_to_converted_fastq_file_name_map;

    
		cReadFiles() { };
		cReadFiles(const vector<string>& read_file_names) { Init(read_file_names); };
		~cReadFiles() { };

		void Init(const vector<string>& read_file_names);
    string base_name_to_read_file_name(const string& base_name);
    vector<string> base_names()
    {
      vector<string> return_value;
      for(vector<cReadFile>::iterator it=this->begin(); it!=this->end(); it++)
      {
        return_value.push_back(it->base_name());
      }
      return return_value;
    }
    
	};

	struct Settings
	{
	public:
		// Set up defaults here
		Settings(const string& _base_output_path = "");
    
    // Constructor for default action
    Settings(int argc, char* argv[]);
    
    //////////////////////////////////
    ////    Settings Variables    ////
    //////////////////////////////////
    
		// Paths populated from location of executable
		string bin_path;                  // absolute path to where this binary resides
    string lib_path;                  // path to where R scripts and images reside
		map<string, string> installed;
    
		// Paths populated relative to output_path
    string output_path;               // command line option

		string sequence_conversion_path;
		string reference_trim_file_name;

		string reference_alignment_path;

		string reference_sam_file_name;
		string reference_fasta_file_name;

		string candidate_junction_path;
		string candidate_junction_fasta_file_name;
		string candidate_junction_faidx_file_name;
		string jc_genome_diff_file_name;
		string preprocess_junction_best_sam_file_name;
		string preprocess_junction_split_sam_file_name;

		string data_path;

		string candidate_junction_alignment_path;
		string candidate_junction_sam_file_name;

		string alignment_correction_path;
		string resolved_reference_sam_file_name;
		string resolved_junction_sam_file_name;

		string reference_features_file_name;

		string unmatched_read_file_name;

		string local_evidence_path;
		string evidence_path;
		string evidence_genome_diff_file_name;
		string final_genome_diff_file_name;

		// Options...

		bool unmatched_reads;
		bool no_unmatched_reads;
		bool add_split_junction_sides;
		bool require_complete_match;

		uint32_t alignment_read_limit;
		uint32_t candidate_junction_read_limit;
		double maximum_candidate_junction_length_factor;
		int32_t maximum_read_length;
		int32_t max_read_mismatches;

		int32_t required_match_length;

		string full_command_line;
		string arguments;
		string predicted_quality_type;
		uint32_t min_quality;
		uint32_t max_quality;
		string run_name;
		uint32_t clean;

		string error_model_method;
		uint32_t base_quality_cutoff;

		////   CandidateJunctions.pm
		bool no_junction_prediction;
		string candidate_junction_score_method;
		uint32_t preprocess_junction_min_indel_split_length;

		//// Scoring to decide which pairs of alignments to the same read to consider
		uint32_t required_extra_pair_total_length;
		uint32_t required_both_unique_length_per_side;
		uint32_t required_one_unique_length_per_side;

		//// Scoring section to choose which ones from list to take
		uint32_t minimum_candidate_junction_pos_hash_score;
		uint32_t minimum_candidate_junction_min_overlap_score;
		uint32_t maximum_inserted_junction_sequence_length;
		uint32_t minimum_candidate_junctions;
		uint32_t maximum_candidate_junctions;

		bool smalt;
		uint32_t max_smalt_diff;

		bool keep_all_intermediates;

		bool no_mutation_prediction;
		bool no_deletion_prediction;
		bool no_alignment_generation;
		uint32_t correction_read_limit;
		bool no_filter_unwanted;
		string unwanted_prefix;
		float mutation_log10_e_value_cutoff;
		float polymorphism_log10_e_value_cutoff;
		float polymorphism_bias_p_value_cutoff;
		float polymorphism_frequency_cutoff;
		uint32_t polymorphism_coverage_both_strands;
		uint32_t polymorphism_reject_homopolymer_length;
		bool no_indel_polymorphisms;
		uint32_t max_rejected_polymorphisms_to_show;
		uint32_t max_rejected_junctions_to_show;
		string byline;
		string website;
		bool strict_polymorphism_prediction;
		uint32_t maximum_read_mismatches;

		string base_output_path;	// main path containing all output

		string converted_fastq_file_name;
		string unwanted_fasta_file_name;
		string sequence_conversion_summary_file_name;
		string sequence_conversion_done_file_name;
		string reference_hash_file_name;
		string reference_alignment_done_file_name;
		string preprocess_junction_done_file_name;
		string coverage_junction_best_bam_unsorted_file_name;
		string coverage_junction_best_bam_file_name;
		string coverage_junction_best_bam_prefix;
		string coverage_junction_distribution_file_name;
		string coverage_junction_plot_file_name;
		string coverage_junction_summary_file_name;
		string coverage_junction_done_file_name;
		string candidate_junction_summary_file_name;
		string candidate_junction_done_file_name;
		string candidate_junction_hash_file_name;
		string candidate_junction_alignment_done_file_name;
		string alignment_correction_summary_file_name;
		string alignment_correction_done_file_name;
		string bam_path;
		string reference_bam_unsorted_file_name;
		string junction_bam_unsorted_file_name;
		string junction_bam_prefix;
		string junction_bam_file_name;
		string bam_done_file_name;
		string reference_bam_prefix;
		string reference_bam_file_name;
		string reference_faidx_file_name;
		string reference_gff3_file_name;

		////// error rates and coverage distribution //////
		string error_calibration_path;
		string error_counts_file_name;
		string error_rates_file_name;
		string error_counts_done_file_name;
		string error_rates_done_file_name;
		string coverage_file_name;
		string unique_only_coverage_distribution_file_name;
		string error_rates_summary_file_name;
		string error_rates_base_qual_error_prob_file_name;
		string plot_error_rates_r_script_file_name;
		string plot_error_rates_fit_r_script_file_name;
		string plot_error_rates_r_script_log_file_name;

		////// mutation identification //////
		string mutation_identification_path;
		string predicted_mutation_file_name;
		string ra_mc_genome_diff_file_name;
		string complete_mutations_text_file_name;
		string complete_coverage_text_file_name;
		string mutation_identification_done_file_name;
		string cnv_coverage_tab_file_name;
		string genome_error_counts_file_name;
		string polymorphism_statistics_input_file_name;

		string polymorphism_statistics_output_file_name;
		string polymorphism_statistics_r_script_file_name;
		string polymorphism_statistics_r_script_log_file_name;
		string polymorphism_statistics_ra_mc_genome_diff_file_name;
		string polymorphism_statistics_done_file_name;

		string output_done_file_name;
		string log_file_name;
		string index_html_file_name;
		string summary_html_file_name;
		string marginal_html_file_name;

		string local_coverage_plot_path;
		string coverage_plot_path;
		string deletions_text_file_name;
		string coverage_plot_file_name;
		string output_calibration_path;
		string unique_only_coverage_plot_file_name;
		string error_rates_plot_file_name;

		// text output files, to be replaced...
		string settings_text_file_name;
		string summary_text_file_name;
		string tiled_coverage_text_file_name;

		string breseq_small_graphic_from_file_name;
		string breseq_small_graphic_to_file_name;

		string long_pairs_file_name;

		uint32_t total_reference_sequence_length;
		uint32_t max_read_length;

    bool junction_prediction; // whether to perform junction prediction step
    
		cReadFiles read_files;
		vector<string> read_file_names;
    vector<string> reference_file_names;

		storable_map<string, Coverage> unique_coverage;

    map<string,string> done_key_messages;
    
		//! TODO @JEB / GRC Setting options needed for HTML outputs
		string print_run_name; 
		bool hide_circular_genome_junctions;
		bool polymorphism_prediction;
		bool lenski_format;
		bool no_evidence;
		bool shade_frequencies;
		bool no_header;
    string html_path(const string& input) const {return "not implemented";}
    bool verbose;

		vector<ExecutionTime> execution_times;
    string time2string(const time_t* timer, const bool& relative);

    //End of settings needed for HTML outputs

		// Utility function to substitute specific details into a generic file name

		static string file_name(const string& file_name_key, const string& substitute, const string& with)
		{
			string s(file_name_key);

			if (substitute.size() > 0)
			{
				size_t pos = s.find(substitute);
				if (pos != string::npos)
				{
					s.replace(pos, 1, with);
				}
			}

			return s;
		}

		// assumes things are in our path for now

		string ctool(string tool_name, bool allow_fail = false)
		{
      if (this->installed[tool_name].size() == 0)
      {
        cerr << "Executable '" << tool_name << "' not found in breseq bin path '" << this->bin_path << "'." << endl;
        if (!allow_fail) exit(-1);
        return "";
      }
			return tool_name;
		}
    
    void record_start_time(const string& message)
    {
      ExecutionTime ex_time;
      time_t this_time = time(NULL);
      ex_time._time_start = this_time;
      ex_time._formatted_time_start = ctime(&this_time);
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

			ExecutionTime ex_time = this->execution_times[i];
			ex_time._message = message;

			time_t this_time = time(NULL);
			ex_time._time_end = this_time;
			ex_time._formatted_time_end = ctime(&this_time);

			//if we had a previous time, calculate elapsed
			if (i < this->execution_times.size())
			{
				time_t time_interval = ex_time._time_end - ex_time._time_start;
				ex_time._time_elapsed = time_interval;
				stringstream ss;
				ss << time_interval;
				ex_time._formatted_time_elapsed = ss.str(); //, 1);
			}
		}

		string base_name_to_read_file_name(string base_name)
		{
      return this->read_files.base_name_to_read_file_name(base_name);
		}

		bool do_step(string done_key, string message);
		void done_step(string done_key);
		void check_installed();

	private:

		void pre_option_initialize();
		void post_option_initialize();
		void init_installed();
	};

	class AnalyzeFastq : public Storable {
  public:
		uint32_t max_read_length;
		uint32_t num_reads;
		uint32_t min_quality_score;
		uint32_t max_quality_score;
		uint32_t num_bases;
		string original_qual_format;
    string quality_format;
		string converted_fastq_name;
    
    AnalyzeFastq() {};
    
    AnalyzeFastq(
                 uint32_t _max_read_length, 
                 uint32_t _num_reads, 
                 uint32_t _min_quality_score, 
                 uint32_t _max_quality_score, 
                 uint32_t _num_bases, 
                 const string& _original_qual_format, 
                 const string& _quality_format,
                 const string& _converted_fastq_name
                )
    : max_read_length(_max_read_length)
    , num_reads(_num_reads)
    , min_quality_score(_min_quality_score)
    , max_quality_score(_max_quality_score)
    , num_bases(_num_bases)
    , original_qual_format(_original_qual_format)
    , quality_format(_quality_format)
    , converted_fastq_name(_converted_fastq_name)
    { }
    
    void serialize(ofstream& f)
    {
      write_to_file(f, max_read_length);
      write_to_file(f, num_reads);
      write_to_file(f, min_quality_score);
      write_to_file(f, max_quality_score);
      write_to_file(f, num_bases);
      write_to_file(f, original_qual_format);
      write_to_file(f, quality_format);
      write_to_file(f, converted_fastq_name);
    }
    void deserialize(ifstream& f)
    {
      read_from_file(f, max_read_length);
      read_from_file(f, num_reads);
      read_from_file(f, min_quality_score);
      read_from_file(f, max_quality_score);
      read_from_file(f, num_bases);
      read_from_file(f, original_qual_format);
      read_from_file(f, quality_format);
      read_from_file(f, converted_fastq_name);
    }
	};

	struct Summary : public Storable
	{
	public:

		struct AlignmentCorrection : public Storable
		{
			map<string, map<string, int32_t> > read_file;

			struct NewJunction
			{
				map<int32_t, int32_t> observed_min_overlap_score_distribution;
				map<int32_t, int32_t> accepted_min_overlap_score_distribution;
				map<int32_t, int32_t> observed_pos_hash_score_distribution;
				map<int32_t, int32_t> accepted_pos_hash_score_distribution;
			} new_junctions;

			void serialize(ofstream& f)
			{
				write_to_file(f, read_file);
				write_to_file(f, new_junctions.observed_min_overlap_score_distribution);
				write_to_file(f, new_junctions.accepted_min_overlap_score_distribution);
				write_to_file(f, new_junctions.observed_pos_hash_score_distribution);
				write_to_file(f, new_junctions.accepted_pos_hash_score_distribution);
			}
			void deserialize(ifstream& f)
			{
				read_from_file(f, read_file);
				read_from_file(f, new_junctions.observed_min_overlap_score_distribution);
				read_from_file(f, new_junctions.accepted_min_overlap_score_distribution);
				read_from_file(f, new_junctions.observed_pos_hash_score_distribution);
				read_from_file(f, new_junctions.accepted_pos_hash_score_distribution);
			}

		} alignment_correction;

		storable_map<string, Coverage> preprocess_coverage;
		storable_map<string, Coverage> unique_coverage;

		struct CandidateJunctionSummaryData : public Storable
		{
			struct Total
			{
				int32_t number;
				int32_t length;
			} total;

			struct Accepted
			{
				int32_t number;
				int32_t length;
				int32_t pos_hash_score_cutoff;
				int32_t min_overlap_score_cutoff;
			} accepted;

			map<int32_t, int32_t> pos_hash_score_distribution;
			map<int32_t, int32_t> min_overlap_score_distribution;

			map<string, map<string, int32_t> > read_file;

      void serialize(ofstream& f)
      {
        write_to_file(f, total);
        write_to_file(f, accepted);
				write_to_file(f, pos_hash_score_distribution);
				write_to_file(f, min_overlap_score_distribution);
				write_to_file(f, read_file);
      }
      
      void deserialize(ifstream& f)
      {
        read_from_file(f, total);
        read_from_file(f, accepted);
				read_from_file(f, pos_hash_score_distribution);
				read_from_file(f, min_overlap_score_distribution);
				read_from_file(f, read_file);
      }

		} candidate_junction;

		struct SequenceConversion : public Storable
		{
			float avg_read_length;
			uint32_t max_qual;
			uint32_t num_reads;
			uint32_t num_bases;
			map<string, string> converted_fastq_name;
			storable_map<string, AnalyzeFastq> reads;
			uint32_t total_reference_sequence_length;
			uint32_t max_read_length;

			void serialize(ofstream& f)
			{
				write_to_file(f, avg_read_length);
        write_to_file(f, max_qual);
				write_to_file(f, num_reads);
				write_to_file(f, num_bases);
				write_to_file(f, converted_fastq_name);
        reads.serialize(f);
        write_to_file(f, total_reference_sequence_length);
				write_to_file(f, max_read_length);
			}
			void deserialize(ifstream& f)
			{
        read_from_file(f, avg_read_length);
        read_from_file(f, max_qual);
				read_from_file(f, num_reads);
				read_from_file(f, num_bases);
				read_from_file(f, converted_fastq_name);
        reads.deserialize(f);
        read_from_file(f, total_reference_sequence_length);
				read_from_file(f, max_read_length);
			}

		} sequence_conversion;

		void serialize(ofstream& f)
		{
      sequence_conversion.serialize(f);
      candidate_junction.serialize(f);
      alignment_correction.serialize(f);
      preprocess_coverage.serialize(f);
      unique_coverage.serialize(f);
    }
    
		void deserialize(ifstream& f)
		{
      sequence_conversion.deserialize(f);
      candidate_junction.deserialize(f);
      alignment_correction.deserialize(f);
      preprocess_coverage.deserialize(f);
      unique_coverage.deserialize(f);
		}
	};

} // breseq namespace

#endif

