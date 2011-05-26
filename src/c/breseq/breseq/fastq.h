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

#include "breseq/common.h"

namespace breseq {
	
	/*! Interface for loading and manipulating files in FASTQ format.
   */
  
  extern const uint8_t k_SANGER_quality_score_offset;
  
  /*! Sequence class.
   */
   
  struct cFastqSequence {
    public:
      std::string m_name;
      std::string m_sequence;
      std::string m_blank;
      std::string m_qualities;
   };   
   

	/*! File class.
	 */ 
  
  class cFastqFile {
    protected:
      
      //! length of longest read
      int32_t    m_max_read_length;

      //highest and lowest scores
      int32_t    m_min_quality_score;    
      int32_t    m_max_quality_score;
    
      //! total number of bases in file
      int64_t    m_total_base;
    
      //! total number of reads in file
      uint64_t   m_total_reads;
    
      //! temporary file name for constructor
      std::string   m_temp_filename;
      std::string   m_quality_format;
    
      //! active fstream that was opened when constructed
      std::fstream     m_file;
      std::fstream     m_file_convert;
      std::fstream     m_temp_file;
    
      std::vector<int>   m_something2sanger;
    
    public:
    
    cFastqFile(const std::string &file_name, std::ios_base::openmode mode, const std::string &temp_path_name); 
    
      void error_in_file_format(int count, int num_reads, int position);
      void check_if_file_opened();
    
      void read_sequence(cFastqSequence &sequence);
      void write_sequence(cFastqSequence &sequence);
    
      void write_summary_file(cFastqSequence &sequence);
    
      void convert_to_sanger(cFastqSequence &sequence);
  };
	
} // breseq namespace

#endif
