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
      //! The number of bases observed with each quality score
      std::vector<uint32_t> m_quality_score_distribution;
      
      //! length of longest read
      int32_t    m_max_read_length;

      //! length of longest read
      int32_t    m_min_quality_score;    
      int32_t    m_max_quality_score;
    
      //! total number of bases in file
      int64_t    m_total_base;
      
      //! active fstream that was opened when constructed
      std::fstream     m_file;
    
    public:
    cFastqFile(const std::string &file_name, std::ios_base::openmode mode);
    
      void open();
    
      bool read_sequence(cFastqSequence &sequence);
      void write_sequence(cFastqSequence &sequence);
    
      void write_quality_score_distribution_file(std::string filename);
      void write_summary_file(std::string filename);
  };
	
} // breseq namespace

#endif
