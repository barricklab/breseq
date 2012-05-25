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
			}
			void deserialize(ifstream& f)
			{
        read_file.deserialize(f);
        read_from_file(f, distance_cutoffs);
        pos_hash_cutoffs.deserialize(f);
				read_from_file(f, observed_pos_hash_score_distribution);
				read_from_file(f, accepted_pos_hash_score_distribution);
			}

		} alignment_resolution;

		storable_map<string, Coverage> preprocess_coverage;
		storable_map<string, Coverage> unique_coverage;

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
    }
    
		void deserialize(ifstream& f)
		{
      sequence_conversion.deserialize(f);
      candidate_junction.deserialize(f);
      alignment_resolution.deserialize(f);
      preprocess_coverage.deserialize(f);
      unique_coverage.deserialize(f);
      preprocess_error_count.deserialize(f);
		}
	};
  
} // breseq namespace

#endif
