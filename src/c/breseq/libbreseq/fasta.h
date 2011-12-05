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

#ifndef _BRESEQ_FASTA_H_
#define _BRESEQ_FASTA_H_

#include "common.h"

using namespace std;

namespace breseq {
	
	/*! Interface for loading and manipulating files in FASTA format.
   */
  
  
  /*! Sequence class.
   */
  class cSequence : public string
  {
    public:
      //! Constructors
      cSequence() : string() {}
      cSequence(const string& value) : string(value) {}
      cSequence(size_t n, const char c) : string(n, c) {}

      //! Public Methods
      inline void operator=(const string &value)
      {
        this->assign(value);
      }

      inline string circularSubStr(const size_t &start_pos, const size_t &n_pos) const
      {
        const size_t  &this_size = this->size();
        const size_t  &max_pos   = start_pos + n_pos;
        const int32_t &size_diff = this_size - max_pos;

        string ret_val;
        if (size_diff >= 0 ) {
          ret_val = this->substr(start_pos, n_pos);
        } else {
          ret_val.append(this->substr(start_pos, this_size + size_diff));
          ret_val.append(this->substr(0, abs(size_diff)));
        }

        return ret_val;
      }

      inline bool isBlank(void)
      {
        size_t size = this->size();
        for (int i = 0; i < size; i++) {
          if ((*this)[i] != ' ') {
            return false;
          }
        }

        return true;
      }
  };

   
  struct cFastaSequence {
    public:
      string m_name;          //>NAME DESCRIPTION 
      string m_description;   //
      cSequence m_sequence;      //sequence ...
   }; 
  

	/*! File class.
	 */ 
  
  class cFastaFile : public fstream {
    
  protected:
    string    m_file_name;
    string    m_current_line;
    uint32_t  m_current_line_num;
    uint32_t  m_bases_per_line;
    
  public:
  
    cFastaFile(const string &file_name, ios_base::openmode mode);
    ~cFastaFile() {};
      
    bool read_sequence(cFastaSequence &sequence);
    void write_sequence(const cFastaSequence &sequence);
    void set_current_line(const string& current_line)
    { m_current_line = current_line;}

  
  };
	
} // breseq namespace

#endif
