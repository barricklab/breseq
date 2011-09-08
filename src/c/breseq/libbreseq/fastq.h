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

#ifndef _BRESEQ_FASTQ_H_
#define _BRESEQ_FASTQ_H_

#include "common.h"

using namespace std;

namespace breseq {
	
	/*! Interface for loading and manipulating files in FASTQ format.
   */
  
  /*! analyze_fastq
   
      Main function for this analysis. Prints summary information
      about fastq and converts file to SANGER format if necessary.
   */  
  AnalyzeFastq analyze_fastq(const string &file_name, const string &convert_file_name);
  
  // Utility function for converting FASTQ files between formats
  void convert_fastq(const string &from_file_name, const string &to_file_name, const string &from_format, const string &to_format);

  
  /*! Sequence class.
   */
   
  struct cFastqSequence {
    public:
      string m_name;      //@NAME
      string m_sequence;  //sequence
      string m_name_plus; //+NAME
      string m_qualities; //quality score characters
      uint32_t m_num_N_bases; // number of N bases, used for filtering
   }; 
  
  /*! Quality score conversion class.
   */
  
  struct cFastqQualityConverter : public vector<uint8_t> {
  public:
    cFastqQualityConverter(const string &from_quality, const string &to_quality);
    ~cFastqQualityConverter() {};
    
    void convert_sequence(cFastqSequence &seq);
  }; 
   

	/*! File class.
	 */ 
  
  class cFastqFile : public fstream {
    
  protected:
    uint32_t  m_current_line;
    string    m_file_name;
    bool      m_needs_conversion;
    
  public:
    // keep track of duplicate read names one after another and append r# to later ones
    bool m_check_for_repeated_read_names;
    string m_last_read_name;
    uint32_t m_repeated_read_name_count;
    
  
    cFastqFile();
    cFastqFile(const string &file_name, ios_base::openmode mode); 
    ~cFastqFile() {};
      
    bool read_sequence(cFastqSequence &sequence);
    void write_sequence(const cFastqSequence &sequence);
  
    bool needs_conversion() { return m_needs_conversion; }
  };

  /*! General sequence helper function.
	 */ 
  
  extern char reverse_complement_lookup_table[256];
  
  inline string reverse_complement(string seq)
	{
		string retval(seq.length(), ' ');
		for (uint32_t i = 0; i < seq.size(); i++)      
      retval[i] = reverse_complement_lookup_table[static_cast<uint8_t>(seq[seq.size() - 1 - i])];
    return retval;
	}
	
} // breseq namespace

#endif
