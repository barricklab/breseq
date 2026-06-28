/*****************************************************************************

 AUTHORS

   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors

 LICENSE AND COPYRIGHT

   Copyright (c) 2008-2010 Michigan State University
   Copyright (c) 2011-2025 The University of Texas at Austin
   Copyright (c) 2025-     Michigan State University

   breseq is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2, or (at your option) any later version.

   SPDX-License-Identifier: GPL-2.0-or-later

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

void tabulate_soft_clipping_counts(
                                   const Settings& settings,
                                   Summary& summary,
                                   const vector<string>& bam_file_names,
                                   const string& fasta_file_name,
                                   const cReferenceSequences& ref_seq_info
                                   );
}

#endif
