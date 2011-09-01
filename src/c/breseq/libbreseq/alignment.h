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

#ifndef _BRESEQ_ALIGNMENT_H_
#define _BRESEQ_ALIGNMENT_H_

#include "common.h"

#include "calculate_trims.h"


using namespace std;

namespace breseq {

// pre-defs
class cReferenceSequences;
class pileup;
class bam_alignment;
  
/*! class alignment
    Represents a single alignment within a pileup.
    This is purely a WRAPPER class. It does not allocate memory.
*/
class alignment_wrapper {
  public:
    friend class bam_alignment;
    //! Constructor.
    alignment_wrapper() : _a(NULL) {};
    alignment_wrapper(const bam1_t* a) : _a(a) {};
    ~alignment_wrapper() {};
  
    // Copying is not allowed. Because it would only copy
    // a pointer, which could go stale.

  friend int first_level_fetch_callback(const bam1_t *b, void *data);
  
  protected:
    alignment_wrapper(const alignment_wrapper& _in) { _a = _in._a; }
    alignment_wrapper&  operator =(const alignment_wrapper& _in) { _a = _in._a; return *this; }
  
  public:
    
    //! Does this alignment have any redundancies?
    bool is_redundant() const;

    //! Number of redundancies at this alignment.
    uint32_t redundancy() const;

    //! Is the read aligned to the reverse strand?
    //  Returns 1 if read aligned to bottom strand, 0 if aligned to top strand
    inline uint32_t reference_target_id() const { return _a->core.tid; }
  
    //!Retrieve name of read.
    inline string read_name() const { return bam1_qname(_a);}

    //! Is the read aligned to the reverse strand?
    //  Returns 1 if read aligned to bottom strand, 0 if aligned to top strand
    inline bool reversed() const { return bam1_strand(_a); }
		
    //! Which strand is the read aligned to?
    //  Returns -1 if read aligned to bottom strand, +1 if aligned to top strand
    inline int32_t strand() const { return (bam1_strand(_a) ? -1 : +1); }
    
    //! Retrieve the query sequence (always on top strand).
    inline uint8_t* read_bam_sequence() const { return bam1_seq(_a); }

    //! Retrieve the query sequence (always on top strand).
    string read_char_sequence() const { 
      string s;
      for (uint32_t p=0; p<read_length(); p++) {
        s += basebam2char(read_base_bam_0(p));
      }
      return s;
    }
    
    //! Calculate the total length of the read.
    inline uint32_t read_length() const { return _a->core.l_qseq; }
  
    //! Retrieve the base at a specified position in the read (was 0-indexed)
    //  Methods available for 0-indexed and 1-indexed coordinates.
    inline uint8_t read_base_bam_0(const uint32_t pos) const { assert(pos<read_length()); return bam1_seqi(read_bam_sequence(), pos); }
    inline uint8_t read_base_bam_1(const uint32_t pos) const { assert(pos<=read_length()); return bam1_seqi(read_bam_sequence(), pos-1); }

    inline char read_base_char_0(const uint32_t pos) const { return basebam2char(read_base_bam_0(pos)); }
    inline char read_base_char_1(const uint32_t pos) const { return basebam2char(read_base_bam_1(pos)); }
	
    //! Retrieve the quality score array. Raw quality scores.
    inline uint8_t* read_base_quality_bam_sequence() const { return bam1_qual(_a); }
  
    inline string  read_base_quality_bam_string() const
    {
      string s;
      uint8_t* qscore = read_base_quality_bam_sequence();      
      for (uint32_t j = 0; j < read_length(); j++)
      {
        s += static_cast<char>(*qscore);
        qscore++;
      }
      return s;
    };
  
    //! Retrieve quality scores as a string of SANGER (+33) offset characters.
    inline string read_base_quality_char_string() const
    {
      string s;
      uint8_t* qscore = read_base_quality_bam_sequence();      
      for (uint32_t j = 0; j < read_length(); j++)
      {
        s += static_cast<char>(qscore[j] + 33);
      }
      return s;
    };
  
    //! Retrieve the quality score of a single base. (was 0-indexed)
    //  Methods available for 0-indexed and 1-indexed coordinates.
    inline uint8_t read_base_quality_0(const uint32_t pos) const { assert(pos<read_length()); return bam1_qual(_a)[pos]; }
    inline uint8_t read_base_quality_1(const uint32_t pos) const { assert(pos<=read_length()); return bam1_qual(_a)[pos-1]; }

    //! Retrieve the index of the read file that contained this alignment.
    uint32_t fastq_file_index() const;
    
    //! Has this alignment been trimmed?
    bool is_trimmed() const;
    
    //! Return number of locations on left of sequence to be trimmed
    uint32_t trim_left() const;
    //! Return number of locations on right of sequence to be trimmed
    uint32_t trim_right() const;

    //! Start and end coordinates of the aligned part of the read. (was 1-indexed)
    //! Start is always < End. reversed() tells you which strand the match was on.
    //  Methods available for 0-indexed and 1-indexed coordinates.
    std::pair<uint32_t,uint32_t> query_bounds_0() const;
    void query_bounds_0(uint32_t& start, uint32_t& end) const;
    std::pair<uint32_t,uint32_t> query_bounds_1() const;
    void query_bounds_1(uint32_t& start, uint32_t& end) const;

    //! Reverse start and end coords if on opposite strand
    std::pair<uint32_t,uint32_t> query_stranded_bounds_1() const;
    void query_stranded_bounds_1(uint32_t& start, uint32_t& end) const;

    //! Starting coordinates of the aligned part of the read (was 1-indexed).
    //  Methods available for 0-indexed and 1-indexed coordinates.
    uint32_t query_start_0() const { return query_start_1()-1; };
    uint32_t query_start_1() const;

    //! Ending coordinates of the aligned part of the read (was 1-indexed).
    //  Methods available for 0-indexed and 1-indexed coordinates.    
    uint32_t query_end_0() const { return query_end_1()-1; };
    uint32_t query_end_1() const;
  
    uint32_t query_match_length() const { return query_end_1() - query_start_1() + 1; };

    //! Starting and ending coordinates of the alignment part of the read
    //  on the reference sequence
    uint32_t reference_start_0() const {return _a->core.pos; } ;
    uint32_t reference_start_1() const {return reference_start_0() + 1; };

    uint32_t reference_end_0() const { return reference_end_1() - 1; };
    uint32_t reference_end_1() const { return bam_calend(&_a->core, bam1_cigar(_a)); };

    std::pair<uint32_t,uint32_t> reference_bounds_0() const 
      { return std::make_pair(reference_start_0(),reference_end_0()); };
    void reference_bounds_0(uint32_t& start, uint32_t& end) const
      { start = reference_start_0(); end = reference_end_0(); }
    std::pair<uint32_t,uint32_t> reference_bounds_1() const
      { return std::make_pair(reference_start_1(),reference_end_1()); };
    void reference_bounds_1(uint32_t& start, uint32_t& end) const
      { start = reference_start_1(); end = reference_end_1(); }

  
    uint32_t reference_match_length() const { return reference_end_1() - reference_start_1() + 1; };

  
    int32_t mate_start_0() const {return _a->core.mpos; } ;
    int32_t mate_start_1() const {return mate_start_0() + 1; };

    inline  bool beginning_to_end_match() const
    { return ((!unmapped()) && (query_start_1() == 1) && (query_end_1() == read_length())); }


    //! Number of bases before this position (on read strand)
    //  that are the same base.
    uint32_t base_repeat_0(uint32_t q_pos_0) const;
  
    //! Operations on CIGAR match string
    inline uint32_t* cigar_array() const { return bam1_cigar(_a); }
    inline uint32_t cigar_array_length() const { return _a->core.n_cigar; }
    inline uint32_t cigar_query_length() const { return bam_cigar2qlen(&_a->core, cigar_array()); };

    string cigar_string() const {
      uint32_t* cigar_list = cigar_array();
      stringstream cigar_string_ss;
      
      for (uint32_t j = 0; j < cigar_array_length(); j++) //foreach my $c (@$cigar_list)
      {
        uint32_t op = cigar_list[j] & BAM_CIGAR_MASK;
        uint32_t len = cigar_list[j] >> BAM_CIGAR_SHIFT;
        cigar_string_ss << len << op_to_char[op]; //$cigar_string += $c->[1] + $c->[0];
      }
      return cigar_string_ss.str();
    }
  
    //! Returns a vector of pairs where internal bam format has been changed to chars and ints
    //  first item in pair is operation, second is length
    inline vector<pair<char,uint16_t> > cigar_pair_array() const
    {
      vector<pair<char,uint16_t> > cigar_pair_list;
      uint32_t* cigar_list = cigar_array();
      for (uint32_t i=0; i<cigar_array_length(); i++)
      {
        uint32_t op = cigar_list[i] & BAM_CIGAR_MASK;
        uint32_t len = cigar_list[i] >> BAM_CIGAR_SHIFT;
        cigar_pair_list.push_back(make_pair( op_to_char[op], len));
      }
      return cigar_pair_list;
    }
  
    inline int32_t isize() const { return _a->core.isize; }
    inline uint32_t quality() const { return _a->core.qual; }
    inline uint16_t flag() const { return _a->core.flag; }
    inline uint8_t* aux_get(const char tag[2]) const { return bam_aux_get(_a, tag); }
    
    //! Is this read unmapped?
    inline bool unmapped() const { return flag() & BAM_FUNMAP; }
    
    inline uint32_t aux_get_i(const char tag[2]) const 
    { 
      uint8_t *auxl = aux_get(tag); 
      return (uint32_t)bam_aux2i(auxl); 
    }

	inline bool proper_pair() const
	{
		return ((_a->core.flag & BAM_FPROPER_PAIR) != 0);
	}

	static const char op_to_char[10];

  protected:
    const bam1_t* _a;             //!< Alignment.
};
  
/*! class pileup_alignment
 Represents a single alignment within a pileup.
 Inherits al methods of alignment, and adds additional ones that get information about pileup position.
 This is purely a WRAPPER class. It does not allocate memory.
 */  
    
class pileup_wrapper : public alignment_wrapper
{
public:
  friend class pileup;
  
  pileup_wrapper() : alignment_wrapper(NULL), _p(NULL) {};
  pileup_wrapper(const bam_pileup1_t* p) : alignment_wrapper(p->b), _p(p) {};

  // Copying is not allowed, except for special friend cases. 
  // Because it only copies a pointer, which could go stale.
  
protected:
  pileup_wrapper(const pileup_wrapper& _in) : alignment_wrapper()
    { _a = _in._a; _p = _in._p; }
  pileup_wrapper&  operator =(const pileup_wrapper& _in) 
    { _a = _in._a; _p = _in._p; return *this; }
  
public:
  
  //! Is this alignment a deletion?
  inline bool is_del() const { return _p->is_del; }
  
  //! Indel length; 0 for no indel, positive for ins and negative for del
  inline int indel() const { return _p->indel; }
  
  //! Indel length; 0 for no indel, -1 for this base deleted and +n for number of inserted bases
  inline int on_base_indel() const {
    int on_indel=indel();
    if(on_indel < 0) on_indel = 0;
    if(is_del()) on_indel = -1; 
    return on_indel;
  }
  
  //! Retrieve the base at a specified insert count relative to the current position
  inline uint8_t on_base_bam(int32_t insert_count=0) const { 
    uint32_t pos0 = query_position_0() + insert_count;
    if (on_base_indel() < insert_count) return '.'; // gap character
    return read_base_bam_0(pos0); 
  };
  
  //! Retrieve the position of the alignment in the query (was 0-indexed).
  //  Methods available for 0-indexed and 1-indexed coordinates.
  inline uint32_t query_position_0() const { return _p->qpos; }
  inline uint32_t query_position_1() const { return _p->qpos+1; }

  //! Is this alignment under the current position?
  bool is_alignment_spanning_position()  const {
    if((_p->qpos+1) > 0)  {
      return true;
    } else {
      return false;
    }
  }
  
  //! Has this alignment been trimmed?
  inline bool is_trimmed() const {
    // is our query position in the left-side trimmed region?
    uint8_t *auxl = bam_aux_get(_a,"XL");
    if(auxl) {
      if((query_position_1()) <= (uint32_t)bam_aux2i(auxl)) {
        return true;
      }
    }	
    // is our query position in the right-side trimmed region?
    uint8_t *auxr = bam_aux_get(_a,"XR");
    if(auxr) {
      if((read_length()-(query_position_1())) <= (uint32_t)bam_aux2i(auxr)) {
        return true;
      }
    }
    
    return false;
  }

  
protected:
  const bam_pileup1_t* _p;     //!< Pileup.

};

// This class is provided when we do want to copy the pointer
// Basically, this is only in pileup
// Do not use this class unless you know what you are doing.
class copiable_pileup_wrapper : public pileup_wrapper {

public:
  copiable_pileup_wrapper() : pileup_wrapper() {};
  copiable_pileup_wrapper(const bam_pileup1_t* p) : pileup_wrapper(p) {};
  
  copiable_pileup_wrapper(const copiable_pileup_wrapper& _in) : pileup_wrapper(_in) { }
  copiable_pileup_wrapper&  operator =(const copiable_pileup_wrapper& _in) 
  { pileup_wrapper::operator=(_in); return *this; }
};
  
/*! class bam_alignment
    This class stores the memory of the bam1_t.
    Thus, it is a copiable version of alignment_wrapper
 */  
  class bam_alignment : public bam1_t, public alignment_wrapper
{
public:
  bam_alignment() : bam1_t(), alignment_wrapper(this)
  {
    l_aux = 0;
    data_len = 0;
    m_data = 0;
    data = NULL;
  }

    bam_alignment(const bam1_t& _in) : alignment_wrapper(this)
  {
    l_aux = 0;
    data_len = 0;
    m_data = 0;
    data = NULL;
    bam_copy1(this, &_in);
  }
    
    bam_alignment(const alignment_wrapper& _in)  : alignment_wrapper(this)
    {
        l_aux = 0;
        data_len = 0;
        m_data = 0;
        data = NULL;
        bam_copy1(this, _in._a);
    }
  
  bam_alignment(const bam_alignment& _in) : alignment_wrapper(_in._a)
  {
    l_aux = 0;
    data_len = 0;
    m_data = 0; 
    data = NULL;
    bam_copy1(this, &_in);
  }
  
  ~bam_alignment()
  {
    free(data);
  }
  
};

typedef counted_ptr<bam_alignment> bam_alignment_ptr;
typedef list<bam_alignment_ptr> alignment_list;  
  
// for debugging
inline void print_alignment_list(const alignment_list& alignments)  
{
  cout << ">>> Begin list" << endl;
  for(alignment_list::const_iterator it=alignments.begin(); it != alignments.end(); it++)
  {
//    printf(" %p\n", (*it)->read_name());
  }
  cout << "<<< End list" << endl;
}

  
class tam_file {

public:
  tam_file() : input_tam(NULL) {}
  tam_file(const string& tam_file_name, const string& fasta_file_name, ios_base::openmode mode);
  ~tam_file();
  
  void open_read(const string& tam_file_name, const string& fasta_file_name);
  void open_write(const string& tam_file_name, const string& fasta_file_name);
  
  bool read_alignments(alignment_list& alignments, bool paired = false);
  void write_alignments(
                        int32_t fastq_file_index, 
                        const alignment_list& alignments, 
                        vector<Trims>* trims = NULL,
                        const cReferenceSequences* ref_seq_info = NULL,
                        bool shift_gaps = false
                        );
  
  void write_moved_alignment(const alignment_wrapper& a, const string& rname, uint32_t fastq_file_index, const string& seq_id, int32_t reference_pos, int32_t reference_strand, int32_t reference_overlap, const uint32_t junction_side, int32_t junction_flanking, int32_t junction_overlap, const Trims* trim = NULL);
  void write_split_alignment(uint32_t min_indel_split_len, const alignment_wrapper& a);

  inline const char* target_name(const alignment_wrapper& a)
  {
    int32_t tid = a.reference_target_id();
    assert (tid < bam_header->n_targets);
    return bam_header->target_name[tid];
  }
  
  bam_header_t* bam_header;
  
protected:
  tamFile input_tam;                // used for input
  ofstream output_tam;              // used for output
  bam_alignment_ptr last_alignment; // contains alignment* last_alignment
};
	
}

#endif
