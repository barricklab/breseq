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

#ifndef _BRESEQ_FASTQ_H_
#define _BRESEQ_FASTQ_H_

#include "common.h"
#include "settings.h"
#include "summary.h"
#include "portable_random.h"

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
                                        const double max_N_fraction,
                                        const uint32_t _long_read_trigger_length,
                                        const uint32_t _long_read_split_length,
                                        const bool _long_read_distribute_remainder,
                                        const uint32_t num_threads = 1
                                        );

  pair<AnalyzeFastqSummary, AnalyzeFastqSummary> normalize_fastq_paired(
                                        const string &r1_file_name,
                                        const string &r1_convert_file_name,
                                        const string &r2_file_name,
                                        const string &r2_convert_file_name,
                                        const uint32_t r1_file_index,
                                        const uint32_t r2_file_index,
                                        const int32_t trim_end_on_base_quality,
                                        const bool filter_reads,
                                        uint64_t current_read_file_bases,
                                        const uint64_t read_file_base_limit,
                                        const uint32_t read_length_min,
                                        const double max_same_base_fraction,
                                        const double max_N_fraction,
                                        const uint32_t long_read_trigger_length,
                                        const uint32_t long_read_split_length,
                                        const bool long_read_distribute_remainder,
                                        const uint32_t num_threads = 1
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
    {
      for (uint8_t b=0; b<base_list_including_N_size; b++) {
        m_base_counts[b] = 0;
      }
    }
    
    bool identical(cFastqSequence& seq);
    
    size_t length() { return m_sequence.length(); }
   };

  void fastq_sequence_trim_end_on_base_quality(cFastqSequence& seq, const uint32_t base_quality);


  //! Per-cycle base quality model for simulated reads.
  //
  //  Deliberately parametric rather than an empirical table: no data files, no new dependencies, and
  //  integer-only so the output is byte-identical on every platform (see portable_random.h).
  //
  //  The mean quality declines QUADRATICALLY with the cycle number, which is the shape real Illumina
  //  data has -- essentially flat over the first half of the read, falling off fastest over the last
  //  few cycles. This replaces a hard-coded cumulative table taken from 2008 GAII runs that had two
  //  problems: it was position-independent (the same distribution at cycle 1 and cycle 100, which is
  //  the least realistic thing a quality model can be), and its marginal implied a ~3.3% per-base
  //  substitution rate, about 3 mismatches in every 100 bp read.
  //
  //  Constant quality needs no special case: set phred_start == phred_end and phred_stdev = 0.
  struct cSimQualityModel {
    uint32_t phred_start;   //!< mean Phred at the first cycle
    uint32_t phred_end;     //!< mean Phred at the last cycle
    uint32_t phred_stdev;   //!< per-base scatter around the cycle mean
    uint32_t phred_min;     //!< floor
    uint32_t phred_max;     //!< ceiling

    cSimQualityModel()
      : phred_start(38), phred_end(28), phred_stdev(4), phred_min(2), phred_max(41) {}

    //! Draw one quality score as a Sanger-encoded (offset 33) character. cycle is 0-based.
    char draw(cPortableRandom& rng, uint32_t cycle, uint32_t read_size) const;
  };

  //! Indel error rates for simulated reads, per million bases. The substitution rate is NOT here:
  //  it derives from each base's own drawn quality (see cSimFastqSequence::is_error_base), which is
  //  what makes the quality string meaningful to breseq's per-quality-bin error calibration.
  struct cSimErrorModel {
    uint32_t insertion_rate_per_million;
    uint32_t deletion_rate_per_million;

    cSimErrorModel() : insertion_rate_per_million(10), deletion_rate_per_million(10) {}  // 1e-5 each
  };

  class cSimFastqSequence {
    public:
      static int32_t SEED_VALUE;

      static map<char, string> random_snp_base_options;

      static char random_insertion_base_options[];

      static char get_random_error_base(cPortableRandom& rng, const char not_this_base);
      static char get_random_insertion_base(cPortableRandom& rng);

      //! P(substitution) = 10^(-phred/10), evaluated from an integer table rather than pow().
      static bool is_error_base(cPortableRandom& rng, uint32_t phred);

      //! Build one read whose footprint in the reference is [left_1, left_1 + read_size - 1].
      //
      //  strand == +1 returns that sequence; strand == -1 returns its reverse complement, so the
      //  read's 5' end sits at left_1 + read_size - 1. Stating the footprint directly (rather than a
      //  strand-dependent buffer origin, as the previous interface did) removes a genuinely
      //  confusing asymmetry: with the old contract a minus-strand read's footprint was actually
      //  [start_1 + read_size, start_1 + 2*read_size - 1], not [start_1, start_1 + read_size - 1].
      //
      //  Simulated deletions consume extra reference bases, so a tail beyond the footprint is
      //  fetched; near the end of a LINEAR sequence that tail is short, and the builder stops
      //  generating deletions rather than reading past it.
      static cFastqSequence simulate_read(const cAnnotatedSequence& ref_sequence,
                                          int32_t left_1,
                                          uint32_t read_size,
                                          int8_t strand,
                                          const cSimQualityModel& quality_model,
                                          const cSimErrorModel& error_model,
                                          cPortableRandom& rng,
                                          const string& read_name,
                                          char& min_quality_char,
                                          bool verbose = false);

      static void simulate_single_ends(const cAnnotatedSequence& sequence,
                                       uint32_t n_reads,
                                       uint32_t read_size,
                                       uint32_t seed,
                                       const cSimQualityModel& quality_model,
                                       const cSimErrorModel& error_model,
                                       string file_name,
                                       bool verbose = false);

      static void simulate_paired_ends(const cAnnotatedSequence& sequence,
                                       uint32_t n_pairs,
                                       uint32_t read_size,
                                       uint32_t fragment_mean,
                                       uint32_t fragment_stdev,
                                       uint32_t seed,
                                       const cSimQualityModel& quality_model,
                                       const cSimErrorModel& error_model,
                                       string pair_1_file_name,
                                       string pair_2_file_name,
                                       bool verbose = false);

      static void simulate_tiled(const cAnnotatedSequence& sequence,
                                 uint32_t read_size,
                                 uint32_t coverage,
                                 const cSimQualityModel& quality_model,
                                 string file_name,
                                 bool verbose);

      //! breseq infers a FASTQ's quality encoding from the LOWEST quality character in the whole
      //  file (see normalize_fastq): a minimum of 59 or more is read as SOLEXA and every score is
      //  then shifted by 31. A simulated library with no low-quality tail therefore gets silently
      //  reinterpreted. Callers pass the minimum character they actually wrote; this warns if it is
      //  high enough to be misread.
      static void warn_if_quality_format_ambiguous(char min_quality_char);
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
    cFastqFile(const string &file_name, ios_base::openmode mode, unsigned int num_threads = 1);
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
