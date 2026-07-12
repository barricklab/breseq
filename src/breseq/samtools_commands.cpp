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

#include "samtools_commands.h"

using namespace std;

namespace breseq {


void samtools_index(const string& bam_file_name) {
  string cmd = "samtools index " + double_quote(bam_file_name)
             + " " + double_quote(bam_file_name + ".bai");
  SYSTEM(cmd, false, false, false);
}

void samtools_sort(const string& unsorted_bam_file_name, const string& sorted_bam_file_name, uint32_t const num_threads) {
  // samtools sort's --threads/-@ is the number of *additional* worker threads on
  // top of the CPU-active main merge thread, so total busy threads = num_threads+1.
  // Pass num_threads-1 (min 0 = single-threaded) so total matches the -j budget.
  uint32_t const sort_threads = (num_threads > 0) ? num_threads - 1 : 0;
  string cmd = "samtools sort --threads " + to_string(sort_threads)
             + " -o " + double_quote(sorted_bam_file_name)
             + " -T " + double_quote(sorted_bam_file_name)
             + " " + double_quote(unsorted_bam_file_name);
  cerr << "[system] " << cmd << endl;

  string stderr_tmp = sorted_bam_file_name + ".stderr.tmp";
  int ret = system((cmd + " 2>" + double_quote(stderr_tmp)).c_str());

  ifstream err_file(stderr_tmp.c_str());
  string line;
  while (getline(err_file, line)) {
    if (ret != 0 || line.find("[bam_sort_core]") == string::npos) {
      cerr << line << "\n";
    }
  }
  err_file.close();
  remove(stderr_tmp.c_str());

  ASSERT(ret == 0, "Error running command:\n[system] " + cmd + "\nResult code: " + to_string(ret));
}

void samtools_faidx(const string& fasta_file_name) {
  string cmd = "samtools faidx " + double_quote(fasta_file_name);
  SYSTEM(cmd, false, false, false);
}

}
