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

#include "breseq/alignment_output.h"
#include "breseq/common.h"
#include "breseq/pileup.h"
#include "assert.h"
#include "boost/optional.hpp"

using namespace std;

namespace breseq {

alignment_output_pileup::alignment_output_pileup(const string& bam, const string& fasta)
  : pileup_base(bam, fasta) { }

alignment_output_pileup::~alignment_output_pileup()  {}


alignment_output::alignment_output(string bam, string fasta, uint32_t in_maximum_to_align)
  : m_maximum_to_align(in_maximum_to_align)
  , m_alignment_output_object(bam, fasta) { }

string alignment_output::html_alignment(const string region) {
  m_alignment_output_object.do_pileup(region);

}



/*! Called for each alignment. pileup& p comes from samtools
 */

void alignment_output_pileup::callback(const pileup& p) {



}
	// tabula rasa

} // namespace breseq

