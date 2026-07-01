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
  
void samtools_index(const string& bam_file_name
                    );

void samtools_sort(
                   const string& unsorted_bam_file_name,
                   const string& sorted_bam_file_name,
                   const uint32_t num_threads
                   );

void samtools_faidx(const string& fasta_file_name);

// Split a BAM produced by paired-end bowtie2 into two per-file BAMs.
// Records whose read name begins with "<r1_file_index>:" go to r1_output_bam;
// those beginning with "<r2_file_index>:" go to r2_output_bam.
// Records appear in the output in the same order as in the input (no sort needed
// when bowtie2 was run with --reorder).
void split_bam_by_file_index(
    const string& input_bam,
    uint32_t r1_file_index, const string& r1_output_bam,
    uint32_t r2_file_index, const string& r2_output_bam
);

// Filter a paired-end BAM by per-read alignment score and split into per-file BAMs and FASTQs.
// For each record, the file_index is parsed from the QNAME prefix before ':'.
// A read passes if it is mapped (FLAG & 0x4 == 0) AND its AS:i tag value >=
//   score_intercept + score_slope * seq_length  (default 0.8*len, matching --score-min L,0,0.8).
// Passing reads are written to the appropriate per-file BAM (paired flags preserved).
// Failing reads (unmapped or low AS) are written as FASTQ to the appropriate .fastq.gz file.
// Pass an empty string for r1_output_fastq / r2_output_fastq to discard failing reads.
void filter_and_split_paired_bam(
    const string& input_bam,
    uint32_t r1_file_index,
    const string& r1_output_bam,
    const string& r1_output_fastq,
    uint32_t r2_file_index,
    const string& r2_output_bam,
    const string& r2_output_fastq,
    double score_slope = 0.8,
    double score_intercept = 0.0
);

}


#endif
