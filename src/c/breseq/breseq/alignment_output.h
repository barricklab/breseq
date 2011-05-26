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

#ifndef _BRESEQ_ALIGNMENT_OUTPUT_H_
#define _BRESEQ_ALIGNMENT_OUTPUT_H_

#include "common.h"

using namespace std;

namespace breseq {

	/*! This class is a factory for generating HTML alignments
	 */
	class alignment_output {
  private:
    uint32_t  m_maximum_to_align;
    
	public:
		//! Constructor.
		alignment_output(uint32_t in_maximum_to_align)
    : m_maximum_to_align(in_maximum_to_align)
    {} ;
    
    /// ... add in all of the functions from Perl alignment_output.pm ...

    //! Output an HTML alignment.
    string html_alignment(string bam, string fasta, string region);
    
  };
	
}

#endif
