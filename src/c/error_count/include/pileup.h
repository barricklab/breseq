#ifndef _BRESEQ_PILEUP_H_
#define _BRESEQ_PILEUP_H_

#include <boost/shared_ptr.hpp>
#include <string>
#include <map>
#include <vector>
#include <sam.h>
#include <faidx.h>

namespace breseq {

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
		pileup_base(const std::string& bam, const std::vector<std::string>& fastas);

		//! Destructor.
		virtual ~pileup_base();
		
		//! Retrieve the reference sequence for the given target and fai index.
		char* get_refseq(int target, int idx);
		
		//! Run the pileup; will trigger callback for each alignment.
		void pileup();
		
		//! Alignment callback.
		virtual int callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile) = 0;
		
	protected:
		friend int first_level_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data);

		samfile_t* _bam; //!< BAM file handle.
		std::map<int,refseq_list_t> _refs; //!< Reference sequences.
	};	
	
} // breseq

#endif
