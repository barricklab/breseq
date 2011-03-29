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

  cFastqFile::cFastqFile(const std::string &file_name, std::ios_base::openmode mode) {
    
    // open for reading or writing...
  }
  
  bool cFastqFile::read_sequence(cFastqSequence &sequence) {
    
  }

  void cFastqFile::write_sequence(cFastqSequence &sequence) {
    // don't need to implement in first stage
  }
  
  void cFastqFile::write_summary_file(std::string filename) {
    
    // Should have lines like: max_read_length 36
    //                         num_reads 1020200
    //                         min_quality_score 4
    //                         max_quality_score 40
    //                         num_bases 62646176
  }
  
} // breseq namespace

