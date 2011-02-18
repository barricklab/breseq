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

// Note, this is essentially a direct port of Perl tabulate coverage function.

#ifndef _BRESEQ_TABULATE_COVERAGE_H_
#define _BRESEQ_TABULATE_COVERAGE_H_

#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <boost/optional.hpp>

#include "breseq/pileup_base.h"
#include "breseq/pileup.h"

namespace breseq {
	  
	void tabulate_coverage( const std::string& bam,
													const std::string& fasta,
                          const std::string& output,
                          const std::string& region,
                          const uint32_t downsample );
  
  
	class tabulate_coverage_pileup : public breseq::pileup_base {
	public:		
		//! Constructor.
		tabulate_coverage_pileup(const std::string& bam, const std::string& fasta, const std::string& output);
		
		//! Destructor.
		virtual ~tabulate_coverage_pileup();		
		
		//! Called for each alignment.
		virtual void callback(const pileup& p);
    
		//! Called at end of fragment.
    void at_end(uint32_t tid, uint32_t seqlen);

	
  protected:
    std::ofstream _output_table;
    uint32_t _last_position_processed;
	};
  
  
  
} // breseq namespace

#endif
