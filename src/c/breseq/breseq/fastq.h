/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2010 Michigan State University

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _FASTQ_H_
#define _FASTQ_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <assert.h>

namespace breseq {
	
	/*! Interface for loading and manipulating files in FASTQ format.
   */
  
  
  /*! Sequence class.
   */
   
  struct cFastqSequence {
    public:
      std::string m_name;
      std::string m_sequence;
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
    
      //! active fstream that was opened when constructed
      std::fstream     m_file;
      std::string      m_filename;
      std::ios_base::openmode m_mode;
    
    public:
    
    cFastqFile(const std::string &file_name, std::ios_base::openmode mode); 
    
      void check_if_file_opened();
    
      bool read_sequence(cFastqSequence &sequence);
      void write_sequence(cFastqSequence &sequence);
    
      void write_quality_score_distribution_file(std::string filename);
      void write_summary_file();
  };
	
} // breseq namespace

#endif
