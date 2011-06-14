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

#include "breseq/settings.h"

using namespace std;

namespace breseq {

  cReadFiles::cReadFiles(const vector<string>& read_file_names) {

    Init(read_file_names);
  }
  
  void cReadFiles::Init(const vector<string>& read_file_names) {
    
    //clear any existing info
    (*this).clear();
    
    uint32_t on_id = 0;
    uint32_t on_error_group =0;
    uint32_t on_paired_end_group = 0;
    
    for (vector<string>::const_iterator it = read_file_names.begin(); it < read_file_names.end(); it++) {
      
      cReadFile rf;
      rf.m_fastq_file_name = *it;
      
      rf.m_paired_end_group = on_paired_end_group++; 
      rf.m_error_group = on_error_group++;         
      rf.m_id = on_id++;
      
      
      // create base name
      rf.m_base_name = rf.m_fastq_file_name;
      // - beginning path
      size_t pos = rf.m_base_name.rfind("/");
      if (pos != string::npos) rf.m_base_name.erase(0,pos+1);
      // - trailing .fastq
      pos = rf.m_base_name.rfind(".fastq");
      if (pos == rf.m_base_name.size() - 6) {
        rf.m_base_name.erase(pos);
      }
      
      (*this).push_back(rf);
    }
  }
  
}
