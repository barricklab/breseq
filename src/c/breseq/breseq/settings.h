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


#ifndef _BRESEQ_SETTINGS_H_
#define _BRESEQ_SETTINGS_H_

#include "breseq/common.h"

namespace breseq {
	  
  // We need to be able to group read files for two reasons
  // 1) They may be paired-end, so we want to map them together
  // 2) They may have the same error rates, so we want to treat them together for error analysis

  struct cReadFile {
  public:
    string m_fastq_file_name;
    string m_base_name;
    uint32_t m_paired_end_group;    // indicated what file contains paired reads
    uint32_t m_error_group;         // indicates what other read files have the same error rates
    uint32_t m_id;                  // index used to refer to this fastq file in BAM
  };
  
  typedef vector<vector<cReadFile> > cReadFileGroup;
  
  
  class cReadFiles : public vector<cReadFile> {
    
  protected:
    
  public:
    cReadFiles(const vector<string>& read_file_names);
    ~cReadFiles();
    
  };
  
  
} // breseq namespace

#endif
