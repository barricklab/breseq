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

#include "libbreseq/common.h"
#include "libbreseq/pileup_base.h"
#include "libbreseq/annotated_sequence.h"

using namespace std;

namespace breseq {
	      
    //Unnecessary stuff that I need to graph stuff
    
    vector<int> readIndices();
    

	void analyze_contingency_loci(const string& bam,
                                const string& fasta,
                                const vector<string>& ref_seq_file_names,
                                const string& output,
                                const string& loci,
                                int strict
                                );
	

	
  // Structure to hold information about repeats
	struct homopolymer_repeat {
        /*
		homopolymer_repeat() {
			bzero(this,sizeof(homopolymer_repeat));
		}*/
		
        string    seq_id;
        uint32_t  start;
        uint32_t  length;
        char      base;
        vector<double> freqs;
	};
    
    struct repeat_stats {
        string region;
        vector<double> freqs;
        
        repeat_stats( string r ){
            region = r;
        }
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
															const string& fasta,
                              const string& loci,
                              int string
                            );
				
		//! Destructor.
		virtual ~contingency_loci_pileup() {};		
		
    //! Main function that analyzes a single repeat locus...
    //  ... should either return or output results!
    void analyze_contingency_locus(const string& region_of_interest);

		//! Called for each alignment.
		virtual void fetch_callback(const alignment_wrapper& a);
    
    void printStats(const string& output, cReferenceSequences& ref_seq_info);
    
    void readIndices( vector<int>& indices, vector<string>&  names, const string& loci);
        
        
	protected:
    // These are used to store information for each run of analyze_contingency_locus
    vector<repeat_stats> repeats;
    homopolymer_repeat current_region;
    string fastaf;
    tam_file tf;
    int strict;
    
    vector<int> indices;
    vector<string> names;
	};  
    
  void writeTAM( const string& tam_file_name, const string& fasta, contingency_loci_pileup clp, alignment_list al );

  
} // breseq namespace


#endif
