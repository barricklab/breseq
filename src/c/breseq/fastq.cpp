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


#include "breseq/fastq.h"

namespace breseq {

  cFastqFile::cFastqFile(const std::string &file_name, std::ios_base::openmode mode)
  {
    m_file.open(file_name.c_str(), mode);
    m_filename          = file_name;
    m_mode              = mode;
    m_total_base        = 0;
    m_max_quality_score = 0;
    m_min_quality_score = INT_MAX;
    m_max_read_length   = 0;
    m_total_reads       = 0;
    
  }
  
  void cFastqFile::check_if_file_opened() {
    assert(m_file.is_open());
  }
  
  bool cFastqFile::read_sequence(cFastqSequence &sequence) {
    
    std::string line, longest_read("");
    uint8_t max_score(0), min_score(UINT8_MAX);
    uint32_t num_reads(0), num_bases(0);
    
    if ( m_file.is_open() ) {
      while( getline( m_file, line )) {
        num_reads++;
        int count(0);
        
        while( count < 4 ) {
          switch (count) {
            case 0:
              sequence.m_name = line;
              getline(m_file, line);
              break;
            case 1:
              sequence.m_sequence = line;
              getline(m_file, line);
              break;
            case 2:
              getline(m_file, line);
              break;
            case 3:
              sequence.m_qualities = line;
              break;
          }
          count++;
        }
        if( sequence.m_sequence.size() > longest_read.size() ) longest_read = line;
        
        num_bases += sequence.m_sequence.size();
        
        for (uint32_t i=0; i<sequence.m_qualities.size(); i++) {
          int this_score(uint8_t(sequence.m_qualities[i]));
          if( this_score > max_score ) max_score = this_score;
          if( this_score < min_score ) min_score = this_score;
        }
        
      }  

    }
    m_total_reads = num_reads;
    m_total_base = num_bases;
    m_max_read_length = longest_read.size();
    m_min_quality_score = min_score-33;
    m_max_quality_score = max_score-33;
    
    return true;
  }

  void cFastqFile::write_sequence(cFastqSequence &sequence) {
    // don't need to implement in first stage
  }
  
  void cFastqFile::write_summary_file() {
    
    // Should have lines like: max_read_length 36
    //                         num_reads 1020200
    //                         min_quality_score 4
    //                         max_quality_score 40
    //                         num_bases 62646176
    std::cout << "max_read_length "   << m_max_read_length   << std::endl;
    std::cout << "num_reads "         << m_total_reads       << std::endl;
    std::cout << "min_quality_score " << m_min_quality_score << std::endl;
    std::cout << "max_quality_score " << m_max_quality_score << std::endl;
    std::cout << "num_bases "         << m_total_base        << std::endl;
  }
  
} // breseq namespace

