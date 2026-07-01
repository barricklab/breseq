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

#include "libbreseq/samtools_commands.h"
#include <zlib.h>

using namespace std;

namespace breseq {


void samtools_index(const string& bam_file_name) {
  string cmd = "samtools index " + double_quote(bam_file_name)
             + " " + double_quote(bam_file_name + ".bai");
  SYSTEM(cmd, false, false, false);
}

void samtools_sort(const string& unsorted_bam_file_name, const string& sorted_bam_file_name, uint32_t const num_threads) {
  string cmd = "samtools sort --threads " + to_string(num_threads)
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

void split_bam_by_file_index(
    const string& input_bam,
    uint32_t r1_file_index, const string& r1_output_bam,
    uint32_t r2_file_index, const string& r2_output_bam
) {
  samFile* in = hts_open(input_bam.c_str(), "rb");
  ASSERT(in, "Could not open BAM for reading: " + input_bam);
  bam_hdr_t* hdr = sam_hdr_read(in);
  ASSERT(hdr, "Could not read BAM header from: " + input_bam);

  samFile* out1 = hts_open(r1_output_bam.c_str(), "wb");
  ASSERT(out1, "Could not open BAM for writing: " + r1_output_bam);
  ASSERT(sam_hdr_write(out1, hdr) == 0, "Failed to write BAM header to: " + r1_output_bam);

  samFile* out2 = hts_open(r2_output_bam.c_str(), "wb");
  ASSERT(out2, "Could not open BAM for writing: " + r2_output_bam);
  ASSERT(sam_hdr_write(out2, hdr) == 0, "Failed to write BAM header to: " + r2_output_bam);

  bam1_t* b = bam_init1();
  while (sam_read1(in, hdr, b) >= 0) {
    const char* qname = bam_get_qname(b);
    const char* colon = strchr(qname, ':');
    ASSERT(colon, string("Malformed read name (missing ':'): ") + qname);
    uint32_t file_index = 0;
    for (const char* p = qname; p < colon; p++)
      file_index = file_index * 10 + static_cast<uint32_t>(*p - '0');
    if (file_index == r1_file_index) {
      ASSERT(sam_write1(out1, hdr, b) >= 0, "Failed to write alignment to: " + r1_output_bam);
    } else if (file_index == r2_file_index) {
      ASSERT(sam_write1(out2, hdr, b) >= 0, "Failed to write alignment to: " + r2_output_bam);
    }
  }
  bam_destroy1(b);
  bam_hdr_destroy(hdr);
  hts_close(in);
  hts_close(out1);
  hts_close(out2);
}

void filter_and_split_paired_bam(
    const string& input_bam,
    uint32_t r1_file_index,
    const string& r1_output_bam,
    const string& r1_output_fastq,
    uint32_t r2_file_index,
    const string& r2_output_bam,
    const string& r2_output_fastq,
    double score_slope,
    double score_intercept
) {
  samFile* in = hts_open(input_bam.c_str(), "rb");
  ASSERT(in, "Could not open BAM for reading: " + input_bam);
  bam_hdr_t* hdr = sam_hdr_read(in);
  ASSERT(hdr, "Could not read BAM header from: " + input_bam);

  samFile* out1 = hts_open(r1_output_bam.c_str(), "wb");
  ASSERT(out1, "Could not open BAM for writing: " + r1_output_bam);
  ASSERT(sam_hdr_write(out1, hdr) == 0, "Failed to write BAM header to: " + r1_output_bam);

  samFile* out2 = hts_open(r2_output_bam.c_str(), "wb");
  ASSERT(out2, "Could not open BAM for writing: " + r2_output_bam);
  ASSERT(sam_hdr_write(out2, hdr) == 0, "Failed to write BAM header to: " + r2_output_bam);

  gzFile fq1 = NULL;
  gzFile fq2 = NULL;
  if (!r1_output_fastq.empty()) {
    fq1 = gzopen(r1_output_fastq.c_str(), "wb");
    ASSERT(fq1, "Could not open FASTQ for writing: " + r1_output_fastq);
  }
  if (!r2_output_fastq.empty()) {
    fq2 = gzopen(r2_output_fastq.c_str(), "wb");
    ASSERT(fq2, "Could not open FASTQ for writing: " + r2_output_fastq);
  }

  bam1_t* b = bam_init1();
  while (sam_read1(in, hdr, b) >= 0) {
    const char* qname = bam_get_qname(b);
    const char* colon = strchr(qname, ':');
    ASSERT(colon, string("Malformed read name (missing ':'): ") + qname);
    uint32_t file_index = 0;
    for (const char* p = qname; p < colon; p++)
      file_index = file_index * 10 + static_cast<uint32_t>(*p - '0');

    bool passes = false;
    if (!(b->core.flag & BAM_FUNMAP)) {
      uint8_t* as_tag = bam_aux_get(b, "AS");
      if (as_tag) {
        int32_t as_val = bam_aux2i(as_tag);
        double threshold = score_intercept + score_slope * static_cast<double>(b->core.l_qseq);
        passes = (as_val >= static_cast<int32_t>(threshold));
      }
    }

    if (passes) {
      samFile* out = (file_index == r1_file_index) ? out1 : out2;
      ASSERT(sam_write1(out, hdr, b) >= 0, string("Failed to write BAM alignment: ") + qname);
    } else {
      gzFile fq = (file_index == r1_file_index) ? fq1 : fq2;
      if (fq) {
        int len = b->core.l_qseq;
        const uint8_t* bseq  = bam_get_seq(b);
        const uint8_t* bqual = bam_get_qual(b);

        string seq(len, 'N');
        string qual(len, '!');
        for (int i = 0; i < len; i++) {
          seq[i]  = seq_nt16_str[bam_seqi(bseq, i)];
          qual[i] = static_cast<char>(bqual[i] + 33);
        }

        if (b->core.flag & BAM_FREVERSE) {
          string rc_seq(len, 'N');
          string rc_qual(len, '!');
          for (int i = 0; i < len; i++) {
            char c = seq[len - 1 - i];
            switch (c) {
              case 'A': rc_seq[i] = 'T'; break;
              case 'T': rc_seq[i] = 'A'; break;
              case 'C': rc_seq[i] = 'G'; break;
              case 'G': rc_seq[i] = 'C'; break;
              default:  rc_seq[i] = 'N'; break;
            }
            rc_qual[i] = qual[len - 1 - i];
          }
          seq  = rc_seq;
          qual = rc_qual;
        }

        gzprintf(fq, "@%s\n%s\n+\n%s\n", qname, seq.c_str(), qual.c_str());
      }
    }
  }

  bam_destroy1(b);
  bam_hdr_destroy(hdr);
  hts_close(in);
  hts_close(out1);
  hts_close(out2);
  if (fq1) gzclose(fq1);
  if (fq2) gzclose(fq2);
}

}
