#ifndef _BRESEQ_ALIGNMENT_H_
#define _BRESEQ_ALIGNMENT_H_

#include <utility>
#include <bam.h>

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
		int32_t redundancy() const;
		
		//! Is this alignment a deletion?
		inline bool is_del() const { return _p->is_del; }
		
		//! Indel length; 0 for no indel, positive for ins and negative for del
		inline int indel() const { return _p->indel; }

		//! Which strand are we on?
		inline uint32_t strand() const { return bam1_strand(_a); }
		
		//! Retrieve the query sequence.
		inline uint8_t* query_sequence() const { return bam1_seq(_a); }
		
		//! Retrieve the position of the alignment in the query.
		inline int32_t query_position() const { return _p->qpos; }
		
		//! Calculate the length of this query out to the last non-clip, non-skip.
		int32_t query_length() const;
		
		//! Retrieve the quality score array.
		inline uint8_t* quality_scores() const { return bam1_qual(_a); }
		
		//! Retrieve the index of the read file that contained this alignment.
		int32_t fastq_file_index() const;
		
		//! Has this alignment been trimmed?
		bool is_trimmed() const;
		
		//! Start and end coordinates of the aligned part of the read (1-indexed). @dk: verify indexing
		std::pair<int32_t,int32_t> query_bounds() const;

	protected:
		const bam_pileup1_t* _p; //!< Pileup.
		const bam1_t* _a; //!< Alignment.
	};
	
}

#endif
