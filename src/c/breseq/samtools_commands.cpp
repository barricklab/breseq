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

#include "libbreseq/samtools_commands.h"

using namespace std;

namespace breseq {


void samtools_index(const string& bam_file_name) {
  string cmd = "samtools index " + double_quote(bam_file_name)
             + " " + double_quote(bam_file_name + ".bai");
  cerr << "[samtools] " << cmd << endl;
  SYSTEM(cmd, false, false, false);
}

void samtools_sort(const string& unsorted_bam_file_name, const string& sorted_bam_file_name, uint32_t const num_threads) {
  string cmd = "samtools sort --threads " + to_string(num_threads)
             + " -o " + double_quote(sorted_bam_file_name)
             + " -T " + double_quote(sorted_bam_file_name)
             + " " + double_quote(unsorted_bam_file_name);
  cerr << "[samtools] " << cmd << endl;
  SYSTEM(cmd, false, false, false);
}

void samtools_faidx(const string& fasta_file_name) {
  string cmd = "samtools faidx " + double_quote(fasta_file_name);
  cerr << "[samtools] " << cmd << endl;
  SYSTEM(cmd, false, false, false);
}

}
