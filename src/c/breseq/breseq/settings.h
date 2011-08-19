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

	struct ExecutionTime {
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
	};

	struct Coverage
	{
		double junction_accept_score_cutoff;
		double deletion_coverage_propagation_cutoff;
		double junction_coverage_cutoff;
		double junction_keep_score_cutoff;
		string nbinom_size_parameter;
		string nbinom_mean_parameter;
		string nbinom_prob_parameter;
		uint32_t average;
		string variance;
		string dispersion;
	};
	
	class cReferenceSequences;

	// We need to be able to group read files for two reasons
	// 1) They may be paired-end, so we want to map them together
	// 2) They may have the same error rates, so we want to treat them together for error analysis

	struct cReadFile
	{
	public:
		string m_fastq_file_name;
		string m_base_name;
		uint32_t m_paired_end_group; // indicated what file contains paired reads
		uint32_t m_error_group; // indicates what other read files have the same error rates
		uint32_t m_id; // index used to refer to this fastq file in BAM
	};

	typedef vector<vector<cReadFile> > cReadFileGroup;

	class cReadFiles : public vector<cReadFile>
	{
	public:

		cReadFiles()
		{
		};
		cReadFiles(const vector<string>& read_file_names);

		~cReadFiles()
		{
		};

		void Init(const vector<string>& read_file_names);

	};

	struct Settings
	{
	public:
		// Set up defaults here
		Settings(const string& _base_output_path = "");

		// Fields
		map<string, string> installed;
		string bin_path;

		// Path to files....
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

		string output_path;
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
		string lib_path;
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

		cReadFiles read_structures;
		vector<string> read_files;
		map<string,string> read_file_to_converted_fastq_file;

		map<string, Coverage> unique_coverage;

		//!@GRC Setting options needed for HTML outputs
		string print_run_name; //!< need to set default to "unnamed"
		bool hide_circular_genome_junctions;
		bool polymorphism_prediction;
		bool lenski_format;
		bool no_evidence;
		bool shade_frequencies;
		bool no_header;
                inline string html_path(string input) {return "not implemented";}
                bool verbose;

		vector<ExecutionTime> execution_times;
                string time2string(const time_t* timer, const bool& relative);

                //@GRC End of settings needed for HTML outputs


		// Utility function to substitute specific details into a generic file name

		static string file_name(const string& file_name_key, const string& substitute = "", const string& with = "")
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

		string ctool(string tool_name) const
		{
			//                      my ($self, $tool_name, $allow_fail) = @_;
			//
			//                      if (!$self->{installed}->{$tool_name})
			//                      {
			//                              if ($allow_fail)
			//                              {
			//                                      $self->warn("Executable \"$tool_name\" not found in breseq bin path\"$self->{bin_path}\".");
			//                                      return undef; # couldn't find it, but it's not an error.
			//                              }
			//                              else
			//                              {
			//                                      $self->throw("Executable \"$tool_name\" not found in breseq bin path\"$self->{bin_path}\".");
			//                              }
			//                      }

			return tool_name;
		}

		string create_path(string path_key)
		{
			string path = file_name(path_key);
			//(-e $path) or Breseq::File::Path::make_path($path) or $self->throw("Could not create path \'$path\'.");
			return path;
		}

		string remove_path(string path_key)
		{
			string path = file_name(path_key);
			//(-e $path) and Breseq::File::Path::remove_tree($path) or $self->throw("Could not remove path \'$path\'.");
			return path;
		}

		void record_end_time(string message)
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

		//transparent to whether read trimming is on

		string read_file_to_fastq_file_name(string read_file)
		{
			/*
			if (defined $self->{read_file_to_converted_fastq_file}->{$read_file}) {
				return $self->{read_file_to_converted_fastq_file}->{$read_file};
			}
			$self->throw if (!defined $self->{read_file_to_fastq_file}->{$read_file});
			return $self->{read_file_to_fastq_file}->{$read_file};
			 */
      return "";
		}


		bool do_step(string done_key, string message);
		void done_step(string done_key);
		void check_installed();

	private:

		void pre_option_initialize();
		void post_option_initialize();
		void init_installed();
	};

	struct AnalyzeFastq {
		uint32_t max_read_length;
		uint32_t num_reads;
		uint32_t min_quality_score;
		uint32_t max_quality_score;
		uint32_t num_bases;
		string quality_format;
		string qual_format;
		string original_qual_format;
		string converted_fastq_name;
	};

	struct Summary : Storable
	{
	public:

		struct AlignmentCorrection : Storable
		{
			map<string, map<string, int32_t> > read_file;

			struct NewJunction
			{
				map<int32_t, int32_t> observed_min_overlap_score_distribution;
				map<int32_t, int32_t> accepted_min_overlap_score_distribution;
				map<int32_t, int32_t> observed_pos_hash_score_distribution;
				map<int32_t, int32_t> accepted_pos_hash_score_distribution;
			} new_junctions;

			void store(string filename)
			{
				ofstream outfile(filename.c_str());
				write_to_file(outfile, *this);
				write_to_file(outfile, read_file);
				write_to_file(outfile, new_junctions.observed_min_overlap_score_distribution);
				write_to_file(outfile, new_junctions.accepted_min_overlap_score_distribution);
				write_to_file(outfile, new_junctions.observed_pos_hash_score_distribution);
				write_to_file(outfile, new_junctions.accepted_pos_hash_score_distribution);
				outfile.close();
			}
			void retrieve(string filename)
			{
				ifstream infile(filename.c_str());
				read_from_file(infile, *this);
				read_from_file(infile, read_file);
				read_from_file(infile, new_junctions.observed_min_overlap_score_distribution);
				read_from_file(infile, new_junctions.accepted_min_overlap_score_distribution);
				read_from_file(infile, new_junctions.observed_pos_hash_score_distribution);
				read_from_file(infile, new_junctions.accepted_pos_hash_score_distribution);
				infile.close();
			}

		} alignment_correction;

		map<string, Coverage> preprocess_coverage;
		map<string, Coverage> unique_coverage;

		struct CandidateJunctionSummaryData : Storable
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

			void store(string filename)
			{
				ofstream outfile(filename.c_str());
				write_to_file(outfile, *this);
				write_to_file(outfile, pos_hash_score_distribution);
				write_to_file(outfile, min_overlap_score_distribution);
				write_to_file(outfile, read_file);
				outfile.close();
			}
			void retrieve(string filename)
			{
				ifstream infile(filename.c_str());
				read_from_file(infile, *this);
				read_from_file(infile, pos_hash_score_distribution);
				read_from_file(infile, min_overlap_score_distribution);
				read_from_file(infile, read_file);
				infile.close();
			}

		} candidate_junction;

		struct SequenceConversion : Storable
		{
			float avg_read_length;
			uint32_t max_qual;
			uint32_t num_reads;
			uint32_t num_bases;
			map<string, string> converted_fastq_name;
			map<string, AnalyzeFastq> reads;

			uint32_t total_reference_sequence_length;
			uint32_t max_read_length;
			cReferenceSequences* reference_sequences;

			void store(string filename)
			{
				ofstream outfile(filename.c_str());
				write_to_file(outfile, *this);
				write_to_file(outfile, converted_fastq_name);
				write_to_file(outfile, reads);
				outfile.close();
			}
			void retrieve(string filename)
			{
				ifstream infile(filename.c_str());
				read_from_file(infile, *this);
				read_from_file(infile, converted_fastq_name);
				read_from_file(infile, reads);
				infile.close();
			}

		} sequence_conversion;

		void store(string filename)
		{
			ofstream outfile(filename.c_str());
			write_to_file(outfile, *this);
			write_to_file(outfile, alignment_correction.read_file);
			write_to_file(outfile, alignment_correction.new_junctions.observed_min_overlap_score_distribution);
			write_to_file(outfile, alignment_correction.new_junctions.accepted_min_overlap_score_distribution);
			write_to_file(outfile, alignment_correction.new_junctions.observed_pos_hash_score_distribution);
			write_to_file(outfile, alignment_correction.new_junctions.accepted_pos_hash_score_distribution);
			write_to_file(outfile, preprocess_coverage);
			write_to_file(outfile, unique_coverage);
			write_to_file(outfile, candidate_junction.pos_hash_score_distribution);
			write_to_file(outfile, candidate_junction.min_overlap_score_distribution);
			write_to_file(outfile, candidate_junction.read_file);
			write_to_file(outfile, sequence_conversion.converted_fastq_name);
			write_to_file(outfile, sequence_conversion.reads);
			outfile.close();
		}
		void retrieve(string filename)
		{
			ifstream infile(filename.c_str());
			read_from_file(infile, *this);
			read_from_file(infile, alignment_correction.read_file);
			read_from_file(infile, alignment_correction.new_junctions.observed_min_overlap_score_distribution);
			read_from_file(infile, alignment_correction.new_junctions.accepted_min_overlap_score_distribution);
			read_from_file(infile, alignment_correction.new_junctions.observed_pos_hash_score_distribution);
			read_from_file(infile, alignment_correction.new_junctions.accepted_pos_hash_score_distribution);
			read_from_file(infile, preprocess_coverage);
			read_from_file(infile, unique_coverage);
			read_from_file(infile, candidate_junction.pos_hash_score_distribution);
			read_from_file(infile, candidate_junction.min_overlap_score_distribution);
			read_from_file(infile, candidate_junction.read_file);
			read_from_file(infile, sequence_conversion.converted_fastq_name);
			read_from_file(infile, sequence_conversion.reads);
			infile.close();
		}
	};


	string capture_system(string command);

} // breseq namespace

#endif

