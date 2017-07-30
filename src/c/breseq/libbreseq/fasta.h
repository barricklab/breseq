/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2017 The University of Texas at Austin

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
   
  class cFastaSequence {
    
    private:
      string m_name;          //>NAME DESCRIPTION
      string m_description;   //
      string m_sequence;      //sequence ...
    
    public:
      cFastaSequence() {}
      cFastaSequence(const string& _name, const string& _description, const string& _sequence)
          : m_name(_name), m_description(_description), m_sequence(_sequence) {}
      ~cFastaSequence() {}
    
    void set_name(const string& _name) { m_name = _name; }
    void set_description(const string& _description) { m_description = _description; }
    void set_sequence(const string& _sequence) { m_sequence = _sequence; }
    
    const string get_description() const { return m_description; }
    const string get_name() const { return m_name; }
    const string get_sequence() const { return m_sequence; }
    
    // Utility function for checking and correcting bounds
    void correct_coordinate_bounds(int64_t &pos_1) const {
      if (pos_1 < 1) {
        WARN_WITH_BACKTRACE("Coordinate (" + to_string<int64_t>(pos_1) + ") requested for get_sequence_1 is < 1.");
        pos_1 = 1;
      }
      
      if (pos_1 > (int64_t)(m_sequence.length())) {
        WARN_WITH_BACKTRACE("Coordinate (" + to_string<int64_t>(pos_1) + ") requested for get_sequence_1 is > length of sequence (" + to_string(m_sequence.length()) + ") .");
        pos_1 = m_sequence.length();
      }
      
    }
    
    void correct_range_bounds(int64_t &start_1, int64_t &end_1) const {
      
      if (start_1 > end_1) {
        WARN_WITH_BACKTRACE("Start coordinate (" + to_string<int64_t>(start_1) + ") is > end coordinate (" + to_string<int64_t>(end_1) + ") for range.");
        end_1 = start_1;
      }
    }
    
    // We some checking here to make sure we don't throw an out-of-bounds error.
    const string get_sequence_1(int64_t start_1, int64_t end_1) const
    {
      
      // No error, this is sometimes intended
      if ((start_1==0) && (end_1==0)) return "";
      
      correct_coordinate_bounds(start_1);
      correct_coordinate_bounds(end_1);
      correct_range_bounds(start_1, end_1);
      return m_sequence.substr(start_1-1, end_1-start_1+1);
    }
    
    char get_sequence_1(int64_t pos_1) const
    {
      correct_coordinate_bounds(pos_1);
      return m_sequence[pos_1-1];
    }

    
    void replace_sequence_1(int64_t start_1, int64_t end_1, const string &replacement_seq) {
      correct_coordinate_bounds(start_1);
      correct_coordinate_bounds(end_1);
      correct_range_bounds(start_1, end_1);
      m_sequence.replace(start_1-1, end_1-start_1+1, replacement_seq);
    }
    
    void insert_sequence_1(int64_t pos_1, const string &insertion_seq) {
      // Allow a value of zero, whcih means insert at the start of the sequence
      if (pos_1 != 0) correct_coordinate_bounds(pos_1);
      m_sequence.insert(pos_1, insertion_seq);
    }

    
    size_t get_sequence_length() const
    {
      return m_sequence.length();
    }
    

   };
  
  inline ostream& operator<<(ostream& out, const cFastaSequence& fasta_sequence)
  {
    out << ">" <<fasta_sequence.get_name() << endl;
    out << fasta_sequence.get_sequence() << endl;
    return out;
  }


	/*! File class.
	 */ 
  
  class cFastaFile : public fstream {
    
  public:
    string    m_file_name;  
    
  protected:
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
