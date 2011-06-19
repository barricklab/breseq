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

using namespace std;

namespace breseq {

  /*! Represents a single alignment within a pileup.
 */
class alignment {
  public:
    //! Constructor.
    alignment(const bam_pileup1_t* p);
    alignment(const bam1_t* a);
    alignment();
    ~alignment();
  
    
    //! Does this alignment have any redundancies?
    bool is_redundant() const;

    //! Number of redundancies at this alignment.
    uint32_t redundancy() const;

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

    //! Retrieve the base at a specified insert count relative to the current position
    inline uint8_t on_base_bam(int32_t insert_count=0) const { 
      uint32_t pos0 = query_position_0() + insert_count;
      if (on_base_indel() < insert_count) return '.'; // gap character
      return read_base_bam_0(pos0); 
    };
    
    //! Calculate the total length of the read.
    inline uint32_t read_length() const { return _a->core.l_qseq; }
  
    //! Retrieve the base at a specified position in the read (was 0-indexed)
    //  Methods available for 0-indexed and 1-indexed coordinates.
    inline uint8_t read_base_bam_0(const uint32_t pos) const { assert((pos>=0) && (pos<read_length())); return bam1_seqi(read_bam_sequence(), pos); }
    inline uint8_t read_base_bam_1(const uint32_t pos) const { assert((pos>0) && (pos<=read_length())); return bam1_seqi(read_bam_sequence(), pos-1); }

    //! Retrieve the position of the alignment in the query (was 0-indexed).
    //  Methods available for 0-indexed and 1-indexed coordinates.
    inline uint32_t query_position_0() const { return _p->qpos; }
    inline uint32_t query_position_1() const { return _p->qpos+1; }
	
    //! Retrieve the quality score array.
    inline uint8_t* read_base_quality_sequence() const { return bam1_qual(_a); }

    //! Retrieve the quality score of a single base. (was 0-indexed)
    //  Methods available for 0-indexed and 1-indexed coordinates.
    inline uint8_t read_base_quality_0(const uint32_t pos) const { assert((pos>=0) && (pos<read_length())); return bam1_qual(_a)[pos]; }
    inline uint8_t read_base_quality_1(const uint32_t pos) const { assert((pos>0) && (pos<=read_length())); return bam1_qual(_a)[pos-1]; }

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

    uint32_t reference_end_0() const;
    uint32_t reference_end_1() const {return reference_end_0() + 1; };

    uint32_t reference_match_length() const { return reference_end_1() - reference_start_1() + 1; };

  
    int32_t mate_start_0() const {return _a->core.mpos; } ;
    int32_t mate_start_1() const {return mate_start_0() + 1; };

    inline  bool beginning_to_end_match() const
    { return ((query_start_1() == 1) && (query_end_1() == read_length())); }


    //! Number of bases before this position (on read strand)
    //  that are the same base.
    uint32_t base_repeat_0(uint32_t q_pos_0) const;
    
    //! Is this alignment under the current position?
    bool is_alignment_spanning_position()  const {
      if((_p->qpos+1) > 0)  {
        return true;
      } else {
        return false;
      }
    }
  
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
  
  
	inline string qseq() const {
	    string seq(_a->core.l_qseq, ' ');
	    for (int32_t i = 0; i < _a->core.l_qseq; i++)
			seq[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(_a),i)];
	    return seq;
	}
	inline int32_t isize() const { return _a->core.isize; }
	inline uint8_t quality() const { return _a->core.qual; }
	inline uint16_t flag() const { return _a->core.flag; }
	inline uint8_t* aux_get(const char tag[2]) const { return bam_aux_get(_a, tag); }
  
  //! Is this read unmapped?
  inline bool unmapped() { return flag() & BAM_FUNMAP; }
  
  inline uint32_t aux_get_i(const char tag[2]) const 
  { 
    uint8_t *auxl = aux_get(tag); 
    return (uint32_t)bam_aux2i(auxl); 
  } 

  protected:
    const bam_pileup1_t* _p; //!< Pileup.
    const bam1_t* _a; //!< Alignment.
  
    static const char op_to_char[10];
};
  
typedef vector<alignment> alignment_list;  
  
class tam_file {

public:
  tam_file(const string& tam_file_name, const string& fasta_file_name, ios_base::openmode mode);
  ~tam_file();
  
  void open_read(const string& tam_file_name, const string& fasta_file_name);
  void open_write(const string& tam_file_name, const string& fasta_file_name);
  
  bool read_alignments(alignment_list& alignments, bool paired = false);
  void write_alignments(int32_t fastq_file_index, alignment_list& alignments, vector<Trim>* trims = NULL);
  void write_split_alignment(uint32_t min_indel_split_len, const alignment& a);
  
protected:
  bam_header_t* bam_header;
  tamFile input_tam;                // used for input
  ofstream output_tam;              // used for output
  vector<bam1_t*> loaded_alignments; // contains alignment* last_alignment
  
  void free_loaded_alignments();
};
	
}

#endif
