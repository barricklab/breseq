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

#ifndef _BRESEQ_CONTINGENCY_LOCI_H_
#define _BRESEQ_CONTINGENCY_LOCI_H_

#include "breseq/common.h"
#include "breseq/pileup_base.h"
#include "breseq/annotated_sequence.h"

using namespace std;

namespace breseq {
	      
	/*! Analyze contingency loci.
	 
	 */
	void analyze_contingency_loci(const string& bam,
                        const string& fasta,
                        const string& output
                        );
	
	
  // Structure to hold information about repeats
	struct homopolymer_repeat {

		homopolymer_repeat() {
			bzero(this,sizeof(homopolymer_repeat));
		}
		
    string    seq_id;
    uint32_t  start;
    uint32_t  length;
		char      base;
	};
	
  typedef vector<homopolymer_repeat> homopolymer_repeat_list;

  // Function to identify repeats
  void identify_homopolymer_repeats(homopolymer_repeat_list& hr, const cReferenceSequences& ref_seqs);
  
	/*! Error-counting class.
	 
	 This class is used by the above identify_mutations() function in order to count errors.
	 */
	class contingency_loci_pileup : public pileup_base {
	public:

		
		//! Constructor.
		contingency_loci_pileup(
                              const string& bam,
															const string& fasta
                            );
				
		//! Destructor.
		virtual ~contingency_loci_pileup() {};		
		
    //! Main function that analyzes a single repeat locus...
    //  ... should either return or output results!
    void analyze_contingency_locus(const string& region_of_interest);

		//! Called for each alignment.
		virtual void fetch_callback(const alignment& a);
    
		
	protected:

    // Add variables that keep track of distribution while fetch_callback is called....

	};  
  
} // breseq namespace

#endif
