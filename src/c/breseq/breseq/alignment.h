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

namespace breseq {

	/*! Represents a single alignment within a pileup.
	 */
	class alignment {
	public:
		//! Constructor.
		alignment(const bam_pileup1_t* p);

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
		inline uint32_t reversed() const { return bam1_strand(_a); }
		
    //! Which strand is the read aligned to?
    //  Returns -1 if read aligned to bottom strand, +1 if aligned to top strand
    inline uint32_t strand() const { return (bam1_strand(_a) ? -1 : +1); }
    
		//! Retrieve the query sequence (always on top strand).
		inline uint8_t* query_bam_sequence() const { return bam1_seq(_a); }

    //! Retrieve the query sequence (always on top strand).
    std::string query_char_sequence() const { 
      std::string s;
      for (uint32_t p=0; p<query_length(); p++) {
        s += basebam2char(query_base_bam_0(p));
      }
      return s;
    }

    //! Retrieve the base at a specified insert count relative to the current position
    inline uint8_t on_base_bam(int32_t insert_count=0) const { 
      uint32_t pos0 = query_position_0() + insert_count;
      if (on_base_indel() < insert_count) return '.'; // gap character
      return query_base_bam_0(pos0); 
    };
    
    //! Retrieve the base at a specified position in the read (was 0-indexed)
    //  Methods available for 0-indexed and 1-indexed coordinates.
		inline uint8_t query_base_bam_0(const uint32_t pos) const { assert((pos>=0) && (pos<query_length())); return bam1_seqi(query_bam_sequence(), pos); }
		inline uint8_t query_base_bam_1(const uint32_t pos) const { assert((pos>0) && (pos<=query_length())); return bam1_seqi(query_bam_sequence(), pos-1); }

		//! Retrieve the position of the alignment in the query (was 0-indexed).
    //  Methods available for 0-indexed and 1-indexed coordinates.
		inline uint32_t query_position_0() const { return _p->qpos; }
		inline uint32_t query_position_1() const { return _p->qpos+1; }

		//! Calculate the total length of the query.
		uint32_t query_length() const;
		
		//! Retrieve the quality score array.
		inline uint8_t* quality_scores() const { return bam1_qual(_a); }

		//! Retrieve the quality score of a single base. (was 0-indexed)
    //  Methods available for 0-indexed and 1-indexed coordinates.
		inline uint8_t quality_base_0(const uint32_t pos) const { assert((pos>=0) && (pos<query_length())); return bam1_qual(_a)[pos]; }
		inline uint8_t quality_base_1(const uint32_t pos) const { assert((pos>0) && (pos<=query_length())); return bam1_qual(_a)[pos-1]; }

		//! Retrieve the index of the read file that contained this alignment.
		uint32_t fastq_file_index() const;
		
		//! Has this alignment been trimmed?
		bool is_trimmed() const;
		
		//! Start and end coordinates of the aligned part of the read. (was 1-indexed)
    //! Start is always < End. reversed() tells you which strand the match was on.
    //  Methods available for 0-indexed and 1-indexed coordinates.
		std::pair<uint32_t,uint32_t> query_bounds_0() const;
		std::pair<uint32_t,uint32_t> query_bounds_1() const;

		//! Starting coordinates of the aligned part of the read (was 1-indexed).    
    //  Methods available for 0-indexed and 1-indexed coordinates.
		uint32_t query_start_0() const { return query_start_1()-1; };
    uint32_t query_start_1() const;

		//! Ending coordinates of the aligned part of the read (was 1-indexed).    
    //  Methods available for 0-indexed and 1-indexed coordinates.    
    uint32_t query_end_0() const { return query_end_1()-1; };
		uint32_t query_end_1() const;
    
    //! Number of bases before this position (on read strand)
    //  that are the same base.
    uint32_t base_repeat_0(uint32_t q_pos_0) const;
    
	protected:
		const bam_pileup1_t* _p; //!< Pileup.
		const bam1_t* _a; //!< Alignment.
	};
	
}

#endif
