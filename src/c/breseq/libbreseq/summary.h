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


#ifndef _BRESEQ_SUMMARY_H_
#define _BRESEQ_SUMMARY_H_

#include "common.h"
#include "storable.h"

namespace breseq{
  

	class Summary : public Storable
	{
	public:
    
		class AlignmentResolution : public Storable
		{
    public:
      
      // @JEB TODO: these statistics are not implemented
      class ReadFile {
      public:
        uint32_t num_unmatched_reads;
        uint32_t num_unique_reads;
        uint32_t num_repeat_reads;
        
        ReadFile()
        : num_unmatched_reads(0)
        , num_unique_reads(0)
        , num_repeat_reads(0)
        {}
        
        void serialize(ofstream& f)
        {
          write_to_file(f, num_unmatched_reads);
          write_to_file(f, num_unique_reads);
          write_to_file(f, num_repeat_reads);
        }
        void deserialize(ifstream& f)
        {
          read_from_file(f, num_unmatched_reads);
          read_from_file(f, num_unique_reads);
          read_from_file(f, num_repeat_reads);
        }
        
      };
      storable_map<string,ReadFile> read_file;
      uint32_t total_unmatched_reads;
      uint32_t total_unique_reads;
      uint32_t total_repeat_reads;
      int32_t max_sam_base_quality_score; // @JEB only filled in for aligned_sam_mode
      
      // these statistics are implemented
      
      //! map by reference seq_id of number of bases from possible overlap
      //! that will be accepted
      map<string,int32_t> distance_cutoffs; 
      //! map by reference seq_id, then list by non-overlap distance possible
      //! of minimum pos hash score needed for accepting a junction
      storable_map<string, storable_vector<int32_t> > pos_hash_cutoffs;   
      map<int32_t, int32_t> observed_pos_hash_score_distribution;
      map<int32_t, int32_t> accepted_pos_hash_score_distribution;

			void serialize(ofstream& f)
			{
        read_file.serialize(f);
        write_to_file(f, distance_cutoffs);
				pos_hash_cutoffs.serialize(f);
				write_to_file(f, observed_pos_hash_score_distribution);
				write_to_file(f, accepted_pos_hash_score_distribution);
        write_to_file(f, max_sam_base_quality_score);
			}
			void deserialize(ifstream& f)
			{
        read_file.deserialize(f);
        read_from_file(f, distance_cutoffs);
        pos_hash_cutoffs.deserialize(f);
				read_from_file(f, observed_pos_hash_score_distribution);
				read_from_file(f, accepted_pos_hash_score_distribution);
        read_from_file(f, max_sam_base_quality_score);
			}

		} alignment_resolution;
    
    class Coverage : public Storable
    {
    public:
      double deletion_coverage_propagation_cutoff;
      double deletion_coverage_seed_cutoff;
      double nbinom_size_parameter;
      double nbinom_mean_parameter;
      double nbinom_prob_parameter;
      double average;
      double variance;
      double dispersion;
      
      void serialize(ofstream& f)
      {
        write_to_file(f, deletion_coverage_propagation_cutoff);
        write_to_file(f, deletion_coverage_seed_cutoff);
        write_to_file(f, nbinom_size_parameter);
        write_to_file(f, nbinom_mean_parameter);
        write_to_file(f, nbinom_prob_parameter);
        write_to_file(f, average);
        write_to_file(f, variance);
        write_to_file(f, dispersion);
      }
      void deserialize(ifstream& f)
      {
        read_from_file(f, deletion_coverage_propagation_cutoff);
        read_from_file(f, deletion_coverage_seed_cutoff);
        read_from_file(f, nbinom_size_parameter);
        read_from_file(f, nbinom_mean_parameter);
        read_from_file(f, nbinom_prob_parameter);
        read_from_file(f, average);
        read_from_file(f, variance);
        read_from_file(f, dispersion);
      }
    };
    
    class AnalyzeFastq : public Storable {
    public:
      uint32_t max_read_length;
      uint32_t original_reads;
      uint32_t homopolymer_filtered_reads;
      uint32_t N_filtered_reads;
      uint32_t num_reads;         // original_reads = homopolymer_filtered_reads + N_filtered_reads + num_reads
      uint32_t min_quality_score;
      uint32_t max_quality_score;
      uint32_t original_num_bases;
      uint32_t num_bases;
      string original_qual_format;
      string quality_format;
      string converted_fastq_name;
      
      AnalyzeFastq() {};
      
      AnalyzeFastq(
                   uint32_t _max_read_length, 
                   uint32_t _original_reads,
                   uint32_t _homopolymer_filtered_reads,
                   uint32_t _N_filtered_reads,
                   uint32_t _num_reads, 
                   uint32_t _min_quality_score, 
                   uint32_t _max_quality_score, 
                   uint32_t _original_num_bases, 
                   uint32_t _num_bases, 
                   const string& _original_qual_format, 
                   const string& _quality_format,
                   const string& _converted_fastq_name
                   )
      : max_read_length(_max_read_length)
      , original_reads(_original_reads)
      , homopolymer_filtered_reads(_homopolymer_filtered_reads)
      , N_filtered_reads(_N_filtered_reads)
      , num_reads(_num_reads)
      , min_quality_score(_min_quality_score)
      , max_quality_score(_max_quality_score)
      , original_num_bases(_original_num_bases)
      , num_bases(_num_bases)
      , original_qual_format(_original_qual_format)
      , quality_format(_quality_format)
      , converted_fastq_name(_converted_fastq_name)

      { }
      
      void serialize(ofstream& f)
      {
        write_to_file(f, max_read_length);
        write_to_file(f, original_reads);
        write_to_file(f, homopolymer_filtered_reads);
        write_to_file(f, N_filtered_reads);
        write_to_file(f, num_reads);
        write_to_file(f, min_quality_score);
        write_to_file(f, max_quality_score);
        write_to_file(f, original_num_bases);
        write_to_file(f, num_bases);
        write_to_file(f, original_qual_format);
        write_to_file(f, quality_format);
        write_to_file(f, converted_fastq_name);

      }
      void deserialize(ifstream& f)
      {
        read_from_file(f, max_read_length);
        read_from_file(f, original_reads);
        read_from_file(f, homopolymer_filtered_reads);
        read_from_file(f, N_filtered_reads);
        read_from_file(f, num_reads);
        read_from_file(f, min_quality_score);
        read_from_file(f, max_quality_score);
        read_from_file(f, original_num_bases);
        read_from_file(f, num_bases);
        read_from_file(f, original_qual_format);
        read_from_file(f, quality_format);
        read_from_file(f, converted_fastq_name);
      }
    };

		storable_map<string, Coverage> preprocess_coverage;
		storable_map<string, Coverage> unique_coverage;

    class PreprocessAlignments : public Storable
    {
    public:
      int32_t aligned_reads;
			int32_t alignments;
      int32_t alignments_split_on_indels;
      int32_t reads_with_alignments_split_on_indels;
      int32_t split_alignments;
      int32_t reads_with_split_alignments;
      
      PreprocessAlignments()
      : aligned_reads(0)
      , alignments(0)
      , alignments_split_on_indels(0)
      , reads_with_alignments_split_on_indels(0)
      , split_alignments(0)
      , reads_with_split_alignments(0)
      {}
      
      void serialize(ofstream& f)
      {
        write_to_file(f, aligned_reads);
        write_to_file(f, alignments);
				write_to_file(f, alignments_split_on_indels);
        write_to_file(f, reads_with_alignments_split_on_indels);
        write_to_file(f, split_alignments);
				write_to_file(f, reads_with_split_alignments);
      }
      
      void deserialize(ifstream& f)
      {
        read_from_file(f, aligned_reads);
        read_from_file(f, alignments);
				read_from_file(f, alignments_split_on_indels);
        read_from_file(f, reads_with_alignments_split_on_indels);
        read_from_file(f, split_alignments);
				read_from_file(f, reads_with_split_alignments);
      }
		} preprocess_alignments;
    
		class CandidateJunctionSummaryData : public Storable
    {
    public:
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
			} accepted;

			map<int32_t, int32_t> pos_hash_score_distribution;
      
      CandidateJunctionSummaryData()
      {}

      void serialize(ofstream& f)
      {
        write_to_file(f, total);
        write_to_file(f, accepted);
				write_to_file(f, pos_hash_score_distribution);
      }
      
      void deserialize(ifstream& f)
      {
        read_from_file(f, total);
        read_from_file(f, accepted);
				read_from_file(f, pos_hash_score_distribution);
      }

		} candidate_junction;

		class SequenceConversion : public Storable
		{
    public:
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
    
    class ErrorCount : public Storable
    {
    public:
      double no_pos_hash_per_position_pr;
      
      void serialize(ofstream& f)
      {
        write_to_file(f, no_pos_hash_per_position_pr);
      }
      void deserialize(ifstream& f)
      {
        read_from_file(f, no_pos_hash_per_position_pr);
      }
    };
    
    storable_map<string, ErrorCount> preprocess_error_count;


    // Overall functions for all of summary

		void serialize(ofstream& f)
		{
      sequence_conversion.serialize(f);
      candidate_junction.serialize(f);
      alignment_resolution.serialize(f);
      preprocess_coverage.serialize(f);
      unique_coverage.serialize(f);
      preprocess_error_count.serialize(f);
      preprocess_alignments.serialize(f);
    }
    
		void deserialize(ifstream& f)
		{
      sequence_conversion.deserialize(f);
      candidate_junction.deserialize(f);
      alignment_resolution.deserialize(f);
      preprocess_coverage.deserialize(f);
      unique_coverage.deserialize(f);
      preprocess_error_count.deserialize(f);
      preprocess_alignments.deserialize(f);
		}
	};
  
} // breseq namespace

#endif
