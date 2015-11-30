/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2008-2010 Michigan State University
 Copyright (c) 2011-2012 The University of Texas at Austin
 
 breseq is free software; you can redistribute it and/or modify it under the
 terms the GNU General Public License as published by the Free Software
 Foundation; either version 1, or (at your option) any later version.
 
 *****************************************************************************/

#ifndef _BRESEQ_SAMTOOLS_COMMANDS_H_
#define _BRESEQ_SAMTOOLS_COMMANDS_H_

#include "common.h"

using namespace std;

namespace breseq {

  
// This class takes care of allocating a C style argument list
// when given a C++ string vector.
class cli_arguments : public vector<string> {
  
  char* _argv[255];
  int _argc;

public:
  
  cli_arguments() : _argc(0) {}
  
  ~cli_arguments() {
    delete_argv();
  }
  
  void delete_argv() {
    for(int i=0; i<_argc; i++) {
      delete[] _argv[i];
    }
    _argc = 0;
  }
  
  // Note that this object maintains the memory allocation
  // so do not use the returned argv  after this object
  // is destroyed!
  void get_args(int* arg_c, char **argv[] ) {
    
    // These are from C library getopt
    // they must be reset between calls to samtools
    optind = 1;
    opterr = 0;
    optopt = 0;
    
    delete_argv();
    
    _argc = this->size();
    
    for(int i=0; i<_argc; i++) {
      _argv[i] = new char[(*this)[i].size()+1];
      strcpy(_argv[i], (*this)[i].c_str());
    }
    
    *arg_c = _argc;
    *argv = _argv;
  }
  
  const string as_string() const{
    return join(*this, " ");
  }
  
};
  
void samtools_index(const string& bam_file_name);

void samtools_sort(
                   const string& unsorted_bam_file_name,
                   const string& sorted_bam_prefix
                   );

void samtools_import(
                     const string& faidx_file_name,
                     const string& sam_file_name,
                     const string& output_bam_file_name
                     );

void samtools_faidx(const string& fasta_file_name);

}


#endif
