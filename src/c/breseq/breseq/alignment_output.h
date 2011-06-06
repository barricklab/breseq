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
#include "pileup_base.h"
#include "pileup.h"
#include "alignment.h"


using namespace std;

namespace breseq {

	/*! This class is a FACTORY for generating HTML alignments
	 */

class alignment_output_pileup : public pileup_base {
  public:
    //! Constructor.
    alignment_output_pileup(const std::string& bam, const std::string& fasta);

    //! Destructor.
    virtual ~alignment_output_pileup();

    //! Called for each alignment.
    virtual void callback(const pileup& aligned_reference);

    //! Called at end of fragment.
    //void at_end(uint32_t tid, uint32_t seqlen);
};

class alignment_output {
  private:
    uint32_t  m_maximum_to_align;
    alignment_output_pileup m_alignment_output_object;
    string m_aligned_reference;
    map<string, string> m_aligned_reads;


  public:
    //! Constructor.
    alignment_output(string bam, string fasta, uint32_t in_maximum_to_align);
    //! Output an HTML alignment.
    string html_alignment(const string region);


};




}//end namespace breseq
#endif
