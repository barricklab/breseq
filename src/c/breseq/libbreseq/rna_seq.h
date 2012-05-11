/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the
  terms the GNU General Public License as published by the Free Software
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_RNA_SEQ_H_
#define _BRESEQ_RNA_SEQ_H_

#include "common.h"

#include "reference_sequence.h"
#include "settings.h"

using namespace std;

namespace breseq {

	class RNASeq
	{
	public:

    static void tam_to_gene_counts(cReferenceSequences& ref_seq_info, const string& fasta_file_name, const vector<string>& tam_files, const string& out_tab_file, bool verbose = false);
    
	}; // class RNASeq

} // namespace breseq

#endif
