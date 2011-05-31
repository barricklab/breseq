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

#ifndef _BRESEQ_RESOLVE_ALIGNMENTS_H_
#define _BRESEQ_RESOLVE_ALIGNMENTS_H_

#include "breseq/common.h"
#include "breseq/settings.h"

using namespace std;

namespace breseq {
  
  extern const string k_junction_name_separator;
  
  struct junction_info_side {
    string    m_seq_id;
    uint32_t  m_position;
    int8_t   m_strand;
    bool      m_redundant;
  };
  
  struct junction_info {
    junction_info_side  m_side_1;
    junction_info_side  m_side_2;
    uint32_t m_alignment_overlap;
    string m_unique_read_sequence;
    uint32_t m_flanking_left;
    uint32_t m_flanking_right;    
  };
  
  void junction_name_split(junction_info& ji, const string& junction_name);
  
  void resolve_alignments( 
                          const bool junction_prediction,
                          const string &reference_fasta,
                          const string &junction_fasta,
                          const string &reference_sam_path,
                          const string &junction_sam_path,
                          const string &resolved_path,
                          const string &data_path,
                          const string &features_file,
                          const cReadFiles &read_files,
                          const uint32_t max_read_length
                          );
	
}

#endif
