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

#include "libbreseq/fasta.h"

using namespace std;

namespace breseq {
  
  //constructor
  cFastaFile::cFastaFile(const string &file_name, ios_base::openmode mode) 
    : fstream(file_name.c_str(), mode)
    , m_file_name(file_name)
    , m_current_line()
    , m_current_line_num(0)
    , m_bases_per_line(60)
  {
    _assert(!(*this).fail(), "Failed to open FASTA file: " + m_file_name);
    
    if (mode == ios_base::in) {
      std::getline(*this, m_current_line);
    }
  }

  // read one sequence record from the file
  bool cFastaFile::read_sequence(cFastaSequence &sequence) {
    
    // clear sequence
    sequence.m_name = "";
    sequence.m_sequence = "";
    
    // We're done, no error
    if (this->eof()) return false;
    
    // Current line should begin with >
    assert(m_current_line[0] == '>');
    sequence.m_name = m_current_line.substr(1);
    
    std::getline(*this, m_current_line);
    m_current_line_num++;
    
    while ((m_current_line[0] != '>') && !this->eof()) {
      
      // should really clean the sequence of spaces and extra line returns
      // and possibly also check for valid characters
      sequence.m_sequence += m_current_line;
      
      std::getline(*this, m_current_line);
      m_current_line_num++;
    }
      
    assert(sequence.m_sequence.length() > 0);
    return true;
  }

  void cFastaFile::write_sequence(const cFastaSequence &sequence) {
    
    (*this) << ">" << sequence.m_name << endl;
    
    uint32_t start = 0;
    while (start < sequence.m_sequence.length()) {
      (*this) << sequence.m_sequence.substr(start,m_bases_per_line) << endl;
      start += m_bases_per_line;
    }

  }
  
} // breseq namespace

