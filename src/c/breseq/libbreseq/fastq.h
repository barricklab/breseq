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

#ifndef _BRESEQ_FASTQ_H_
#define _BRESEQ_FASTQ_H_

#include "libbreseq/common.h"
#include "libbreseq/settings.h"
#include "libbreseq/summary.h"

namespace breseq {

class cAnnotatedSequence;

	/*! Interface for loading and manipulating files in FASTQ format.
   */
  
  /*! normalize_fastq
   
      Main function for this analysis. Prints summary information
      about fastq and converts file to SANGER format if necessary.
   */  
  AnalyzeFastqSummary normalize_fastq(
                                        const string &file_name, 
                                        const string &convert_file_name, 
                                        const uint32_t file_index, 
                                        const int32_t trim_end_on_base_quality, 
                                        const bool filter_reads,
                                        uint64_t current_read_file_bases,
                                        const uint64_t read_file_base_limit,
                                        const uint32_t read_length_min,
                                        const double max_same_base_fraction,
                                        const double max_N_fraction
                                        );
  
  // Utility function for converting FASTQ files between formats
  void convert_fastq(const string &from_file_name, const string &to_file_name, const string &from_format, const string &to_format, bool reverse_complement = false);
  
  /*! Sequence class.
   */
   
  struct cFastqSequence {
    public:
      string   m_name;        //@NAME
      string   m_sequence;    //sequence
      string   m_name_plus;   //+NAME
      string   m_qualities;   //quality score characters
      bool     m_numerical_qualities; // quality scores read were read from numerical format
      uint32_t m_base_counts[base_list_including_N_size]; // number of each base, including N. Used for filtering.
    
    
    cFastqSequence() 
    : m_numerical_qualities(false)
    {}
    
    bool identical(cFastqSequence& seq);
    
    size_t length() { return m_sequence.length(); }
   };

  void fastq_sequence_trim_end_on_base_quality(cFastqSequence& seq, const uint32_t base_quality);


  class cSimFastqSequence {
    public:
      static int32_t SEED_VALUE;

      class GaussianRNG {
        static void box_muller_transform(float* z0, float* z1);  
        static const double PI;

        public:
          GaussianRNG(int mean, int stdev);
          int32_t sample();
          vector<int32_t> samples(int size);

        private:
          float m_mean;
          float m_stdev;

          float m_z0;
          float m_z1;
          float m_store;

      };

      static map<uint32_t, uint32_t> qscore_cumulative_probability_table;
      static map<char, string> random_snp_base_options;

      static char random_insertion_base_options[];

      static char get_random_quality_score(void);
      static char get_random_error_base(const char not_this_base);
      static char get_random_insertion_base(void);

      static bool is_random_error_base(char ascii_qscore);
      static bool is_random_deletion_base(void);
      static bool is_random_insertion_base(void);

      static cFastqSequence simulate(const cAnnotatedSequence& ref_sequence,
                                     uint32_t start_1,
                                     uint32_t read_size,
                                     int8_t strand,
                                     uint32_t id = 0,
                                     bool verbose = false);

      static void simulate_single_ends(const cAnnotatedSequence& sequence,
                                       uint32_t n_reads,
                                       uint32_t read_size, 
                                       string file_name,
                                       bool verbose = false);

      static void simulate_paired_ends(const cAnnotatedSequence& sequence,
                                       uint32_t n_reads,
                                       uint32_t read_size, 
                                       uint32_t mean,
                                       uint32_t stdev,
                                       string pair_1_file_name,
                                       string pair_2_file_name,
                                       bool verbose = false);


  };

  
  /*! Quality score conversion class.
   */
  
  struct cFastqQualityConverter : public vector<uint8_t> {
  public:    
    cFastqQualityConverter(const string &_from_quality_format, const string &_to_quality_format);
    ~cFastqQualityConverter() {};
    
    string from_quality_format;
    string to_quality_format;
    string from_quality_type;
    string to_quality_type;
    int32_t from_chr_offset;
    int32_t to_chr_offset;
    
    void convert_sequence(cFastqSequence &seq);
    
    static string predict_fastq_file_format(const string& file_name, uint64_t& num_original_reads, uint64_t& num_original_bases, uint32_t& read_length_min, uint32_t& read_length_max, uint8_t& min_quality_score, uint8_t& max_quality_score);

  };
   

	/*! File class.
	 */ 
  
  class cFastqFile : public flexgzfstream {
    
  protected:
    uint32_t  m_current_line;
    string    m_file_name;
    
  public:
    // keep track of duplicate read names one after another and append r# to later ones
    bool m_check_for_repeated_read_names;
    string m_last_read_name;
    uint32_t m_repeated_read_name_count;
    
  
    cFastqFile();
    cFastqFile(const string &file_name, ios_base::openmode mode);
    ~cFastqFile() { };
      
    bool read_sequence(cFastqSequence &sequence, cFastqQualityConverter& fqc);
    void write_sequence(const cFastqSequence &sequence);
  };

  /*! General sequence helper function.
	 */ 
  
  extern char reverse_complement_lookup_table[256];
  
  inline string reverse_complement(const string& seq)
	{
		string retval(seq.length(), ' ');
		for (uint32_t i = 0; i < seq.size(); i++)      
      retval[i] = reverse_complement_lookup_table[static_cast<uint8_t>(seq[seq.size() - 1 - i])];
    return retval;
	}
  
  inline char reverse_complement(const char& seq)
	{
		char retval(' ');
    retval = reverse_complement_lookup_table[static_cast<uint8_t>(seq)];
    return retval;
	}
  
  inline cFastqSequence reverse_complement(const cFastqSequence& _seq) {
    cFastqSequence seq = _seq;
    seq.m_sequence = reverse_complement(seq.m_sequence);
    seq.m_qualities = reverse_string(seq.m_qualities);
    seq.m_name += "_RC";
    return seq;
  }
	
} // breseq namespace

#endif
