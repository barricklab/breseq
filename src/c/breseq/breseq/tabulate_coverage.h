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

// Note, this is essentially a direct port of Perl tabulate coverage function.

#ifndef _BRESEQ_TABULATE_COVERAGE_H_
#define _BRESEQ_TABULATE_COVERAGE_H_

#include <boost/optional.hpp>

#include "common.h"
#include "pileup_base.h"
#include "pileup.h"

namespace breseq {
	  
	void tabulate_coverage( const std::string& bam,
													const std::string& fasta,
                          const std::string& output,
                          const std::string& region,
                          const uint32_t downsample,
                          const std::string& read_start_output,
                          const std::string& gc_output  );
  
  
	class tabulate_coverage_pileup : public breseq::pileup_base {
	public:		
		//! Constructor.
		tabulate_coverage_pileup(const std::string& bam, const std::string& fasta, const std::string& output,
      const std::string& read_begin_output, const std::string& gc_output);
		
		//! Destructor.
		virtual ~tabulate_coverage_pileup();		
		
		//! Called for each alignment.
		virtual void pileup_callback(const pileup& p);
    
		//! Called at end of fragment.
    void at_target_end(const uint32_t tid);

	
  protected:
    std::ofstream m_output_table, m_read_begin_output, m_gc_output;
    
    std::map<std::string,uint32_t> m_read_begin_top_bins;
    std::map<std::string,uint32_t> m_read_begin_bot_bins;
    std::map<std::string,uint32_t> m_ref_begin_top_bins;
    std::map<std::string,uint32_t> m_ref_begin_bot_bins;
    
    std::vector<uint32_t> m_gc_content_bins;
    
	};
  
  
  
} // breseq namespace

#endif
