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

#include <iostream>

#include "breseq/pileup.h"


/*! Constructor.
 */
breseq::pileup::pileup(uint32_t tid, uint32_t pos_1, int n, const bam_pileup1_t *pile, breseq::pileup_base& pb)
: _tid(tid)
, _pos_1(pos_1)
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
