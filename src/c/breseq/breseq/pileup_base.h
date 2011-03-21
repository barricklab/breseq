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

#ifndef _BRESEQ_PILEUP_BASE_H_
#define _BRESEQ_PILEUP_BASE_H_

#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <string>
#include <map>
#include <vector>
#include <sam.h>
#include <faidx.h>


namespace breseq {
	
	// pre-decs
	class pileup;
	int first_level_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data);
	
	//! Helper struct to manage a single reference sequence.
	struct reference_sequence {
	public:		
		reference_sequence(const std::string& fasta_filename, const std::string& target);
		~reference_sequence();
		
		faidx_t* _ref; //!< FAI file handle.
		char* _seq; //!< Reference sequence (ascii).
		int _len; //<! Length of reference sequence.
		
	private:
		// not allowed:
		reference_sequence(const reference_sequence& that);
		reference_sequence& operator=(const reference_sequence& that);
	};
	
	
	/*! Class to assist in developing pileup-related functionality.
	 */
	class pileup_base {
	public:
		//! Type for a list of reference sequences.
		typedef std::vector<boost::shared_ptr<reference_sequence> > refseq_list_t;
		
		//! Constructor.
		pileup_base(const std::string& bam, const std::string& fasta);

		//! Destructor.
		virtual ~pileup_base();
		
		//! Retrieve the reference sequence for the given target and fai index.
		char* get_refseq(uint32_t target) const;

		//! Retrieve the name of the given target.
		const char* target_name(uint32_t target) const { return _bam->header->target_name[target]; }

		//! Retrieve the name of the given target.
		const uint32_t target_length(uint32_t target) const { return _refs[target]->_len; }

    char reference_base_char_1(uint32_t target, uint32_t pos1) const  { return get_refseq(target)[pos1-1]; } ;
		
		//! Do the pileup; will trigger callback for each alignment.
		void do_pileup();
		
		//! Pileup callback.
		virtual void callback(const pileup& p) = 0;

		//! Called after the pileup has completed.
		virtual void at_end(uint32_t tid, uint32_t seqlen) { }

	protected:
		friend int first_level_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data);

		samfile_t* _bam; //!< BAM file handle.
		refseq_list_t _refs; //!< Reference sequences.
		boost::optional<uint32_t> _last_tid; //!< The "last target" for which the first-level-callback was called.
	};
	
} // breseq

#endif
