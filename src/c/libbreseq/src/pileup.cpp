#include <iostream>

#include "breseq/pileup.h"


/*! Constructor.
 */
breseq::pileup::pileup(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, breseq::pileup_base& pb)
: _tid(tid)
, _pos(pos)
, _pb(pb) {
	
	// build our alignment objects:
	reserve(static_cast<std::size_t>(n));
	for(int i=0; i<n; ++i) {
		push_back(alignment(&pile[i]));
	}		
}


/*! Retrieve the reference sequence for this pileup.
 */
char* breseq::pileup::reference_sequence() const { 
	return _pb.get_refseq(_tid); 
}
