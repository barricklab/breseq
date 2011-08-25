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

#ifndef _BRESEQ_CALCULATE_TRIMS_H_
#define _BRESEQ_CALCULATE_TRIMS_H_

#include "common.h"

using namespace std;

namespace breseq {

  void calculate_trims( const string& in_fasta, const string& in_output_path);
  void calculate_trims_1( const string& in_fasta, const string& in_output_path);

  /*! Trim classes
   
	 */ 
  
  struct Trims
	{
    Trims() : L(0), R(0) {}
		uint32_t L;
		uint32_t R;
	};
  
  class SequenceTrims {
  private: 
    uint8_t * trim_data;
    // Trim data stores LEFT trims in 0..m_length-1 and RIGHT trims in m_length..m_length*2-1
    uint32_t m_length;
  public: 
    
    SequenceTrims() : trim_data(NULL), m_length(0) {  };
    
    SequenceTrims(const SequenceTrims& _in) : trim_data(NULL), m_length(0)
    {      
      if ((_in.m_length == 0) && (_in.trim_data == NULL)) return;
      
      m_length = _in.m_length;
      trim_data = new uint8_t[2*m_length];
      memcpy(trim_data, _in.trim_data, 2*m_length);
    }
    
    SequenceTrims(const string& _in_seq);
    
    ~SequenceTrims() { if (trim_data) delete[] trim_data; };
    
    inline void ReadFile(const string& trim_file_name, uint32_t in_seq_length) 
    { 
      // delete any old content
      if (trim_data) delete[] trim_data;
      
      m_length = in_seq_length;
      trim_data = new uint8_t[2*in_seq_length];
      ifstream in_trim_file(trim_file_name.c_str(), ios::out | ios::binary);
      in_trim_file.read(reinterpret_cast<char *>(trim_data), 2*in_seq_length);
      assert(!in_trim_file.fail());
      assert(in_trim_file.gcount() == 2*in_seq_length);
      in_trim_file.close();
    }

    inline void WriteFile(const string& trim_file_name) 
    {       
      ofstream out(trim_file_name.c_str(), ios::out | ios::binary);
      out.write(reinterpret_cast<char *>(trim_data), 2*m_length);
      out.close();
    }
    
    uint8_t left_trim_0(uint32_t pos_0) const
    { 
      assert(pos_0 < m_length); 
      return trim_data[pos_0];
    }
    
    uint8_t right_trim_0(uint32_t pos_0) const
    { 
      assert(pos_0 < m_length); 
      return trim_data[m_length+pos_0];
    }
    
    Trims operator[] (uint32_t pos_0) const
    { 
      assert(pos_0 < m_length); 
      Trims t; 
      t.R = trim_data[pos_0]; 
      t.L = trim_data[m_length+pos_0];  
      return t; 
    }
  };
  
  typedef vector<SequenceTrims> SequenceTrimsList;
  
  
} // breseq

#endif
