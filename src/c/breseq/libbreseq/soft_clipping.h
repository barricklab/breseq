/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2022 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_SOFT_CLIPPING_H_
#define _BRESEQ_SOFT_CLIPPING_H_

#include "common.h"

#include "reference_sequence.h"
#include "genome_diff.h"
#include "pileup_base.h"


using namespace std;

namespace breseq {

bool is_read_soft_clipped(
                          Settings& settings,
                          Summary& summary,
                          cReferenceSequences& ref_seq_info,
                          const string&  output_file_name
                          );


void analyze_soft_clipping(
                           const vector<string>& bam_file_names,
                           const string& fasta_file_name,
                           const string& output_file_name,
                           const uint32_t minimum_clipped_bases
                           );
}

#endif
