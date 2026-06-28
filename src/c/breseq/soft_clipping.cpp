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

#include "libbreseq/soft_clipping.h"

#include "libbreseq/alignment.h"
#include "libbreseq/identify_mutations.h"
#include "libbreseq/reference_sequence.h"
#include "libbreseq/output.h"
#include "libbreseq/settings.h"
#include "libbreseq/summary.h"

using namespace std;

namespace breseq {

void analyze_soft_clipping(
                           const vector<string>& bam_file_names,
                           const string& fasta_file_name,
                           const string& output_file_name,
                           const uint32_t minimum_clipped_bases
                           )
{
  // Create the output path if necessary
  string path = path_to_dirname(output_file_name);
  create_path(path);
  
  ofstream out_file(output_file_name);
  out_file << join(make_vector<string>("read_name")("seq_id")("position")("direction")("num_bases")("strand"), ",") << endl;
  
  for(vector<string>::const_iterator it=bam_file_names.begin(); it != bam_file_names.end(); it++ ) {
    string bam_file_name = *it;
    cerr << "    Processing " << bam_file_name << endl;
    bam_file final_bam_file(bam_file_name, fasta_file_name, ios::in);
    
    alignment_list al;
    
    while (final_bam_file.read_alignments(al, false)) {
      
      uint32_t num_alignments = al.size();
      
      // We don't want any with multiple alignments...
      // this should have been resolved or they are uninformative.
      //if (num_alignments > 1) continue;
      
      for (alignment_list::iterator it = al.begin(); it != al.end(); it++) {
        bam_alignment& a = *it->get();
        
        // is it mapped?
        if (a.unmapped()) continue;
        
        // check if it has our mark for being multiply mapped and skip if so
        uint32_t num_equivalent_alignments;
        if (a.aux_get_i("X1", num_equivalent_alignments)) {
          if (num_equivalent_alignments > 1) continue;
        }
        
        // Check to see if this is one side of a junction match
        // in which case we need to ignore the soft trimming on one side
        // But, we don't know which side... so
        int32_t ignore_side = -1;
        uint32_t junction_side;
        // either 1 or 2, this tells us which side to ignore
        if (a.aux_get_i("XJ", junction_side)) {
          ignore_side = junction_side;
        }
        
        // NOTE: a.read_length() can be zero when there are multiple alignments
        //       for minimap2. To avoid these problems, we need to use the CIGAR string to calculate the length
        uint32_t read_length = a.cigar_query_length();
        
        uint32_t query_begin_soft_clipping = a.query_start_1() - 1;
        uint32_t query_end_soft_clipping = read_length - a.query_end_1();
        
        // debugging
        //if (a.read_name() == "e176333e-29ea-4179-be3c-890e0e65dc60") {
        //  cout << "DEBUG!" << endl;
        //}
        
        //cout << a.read_name() << endl;
        //cout << "  Read Length/Strand: " << read_length << "  " << a.strand() << endl;
        //cout << "  Reference: " << a.reference_start_1() << "-" << a.reference_end_1() << endl;
        //cout << "  Query: " << a.query_start_1() << "-" << a.query_end_1() << endl;
        //cout << "  Clipping: " << query_begin_soft_clipping << "   " << query_end_soft_clipping << endl;
        
        
        if (query_begin_soft_clipping > minimum_clipped_bases) {
          
          uint32_t clipping_coord = a.reference_start_1();
          int32_t clipping_direction = -1;
          uint32_t num_bases = query_begin_soft_clipping;
          
          //cout << "PRINTED BEGIN" << endl;
          out_file << join(
                           make_vector<string>
                           (a.read_name())
                           (final_bam_file.target_name(a))
                           (to_string(clipping_coord))
                           (to_string(clipping_direction))
                           (to_string(num_bases))
                           (to_string(a.strand())),
                           ",") << endl;
        }
        
        if (query_end_soft_clipping > minimum_clipped_bases) {
          
          uint32_t clipping_coord = a.reference_end_1();
          int32_t clipping_direction = +1;
          uint32_t num_bases = query_end_soft_clipping;
          
          //cout << "PRINTED END" << endl;
          out_file << join(
                           make_vector<string>
                           (a.read_name())
                           (final_bam_file.target_name(a))
                           (to_string(clipping_coord))
                           (to_string(clipping_direction))
                           (to_string(num_bases))
                           (to_string(a.strand())),
                           ",") << endl;
          
        }
        
      } // end for loop
      
    } // end BAM file loop
    } // end while loop
    
  } // end function
 
// Computes per-position soft-clipping counts and genome-wide totals.
// Boundary positions (1 and seq_length) are excluded from both counts and the null model,
// since reads always soft-clip at reference edges.
void tabulate_soft_clipping_counts(
                                   const Settings& settings,
                                   Summary& summary,
                                   const vector<string>& bam_file_names,
                                   const string& fasta_file_name,
                                   const cReferenceSequences& ref_seq_info
                                   )
{
  uint32_t minimum_clipped_bases = settings.soft_clipping_minimum_bases;

  // Build seq_length lookup from ref_seq_info
  map<string, uint32_t> seq_lengths;
  for (uint32_t i = 0; i < ref_seq_info.size(); i++) {
    seq_lengths[ref_seq_info[i].m_seq_id] = ref_seq_info[i].m_length;
  }

  // Per-position soft-clip event counts (numerator only).
  //
  // Direction/strand semantics (correct for reads on EITHER strand, because BAM stores
  // CIGAR/SEQ in forward-reference orientation):
  //   direction -1 : a *leading* CIGAR soft-clip (query_begin_soft_clipping) at reference_start_1.
  //                  This captures both a top-strand read clipped at its 5' beginning AND a
  //                  bottom-strand read clipped at its 3' end -- both produce a leading clip and
  //                  match toward higher reference coordinates.
  //   direction +1 : a *trailing* CIGAR soft-clip (query_end_soft_clipping) at reference_end_1.
  //                  Captures a top-strand 3' clip and a bottom-strand 5' clip.
  map<string, map<uint32_t, uint32_t>> clipped_neg; // direction = -1, indexed by reference_start_1
  map<string, map<uint32_t, uint32_t>> clipped_pos; // direction = +1, indexed by reference_end_1

  // Per-sequence "read-through" coverage, accumulated as difference arrays. A read only counts
  // toward a position's denominator if it reads through that position by at least
  // minimum_clipped_bases aligned bases on the relevant side -- i.e. it extends as far past the
  // position as a clipped read must be clipped to count. This makes the numerator and denominator
  // symmetric and keeps reads that barely reach a position out of the denominator.
  //   direction -1 (leading-clip wall at p): a read reads through on the low side if it matches at
  //       least min_bases below p, i.e. reference_start_1 <= p - min_bases. It contributes over
  //       [reference_start_1 + min_bases, reference_end_1].
  //   direction +1 (trailing-clip wall at p): a read reads through on the high side if it matches at
  //       least min_bases above p, i.e. reference_end_1 >= p + min_bases. It contributes over
  //       [reference_start_1, reference_end_1 - min_bases].
  map<string, vector<int32_t> > cov_minus_diff; // direction -1 read-through coverage
  map<string, vector<int32_t> > cov_plus_diff;  // direction +1 read-through coverage

  uint64_t total_clipped_read_ends = 0;  // total clip events, both directions (a read clipped at
                                         // both ends counts twice)

  for (vector<string>::const_iterator bam_it = bam_file_names.begin(); bam_it != bam_file_names.end(); bam_it++) {
    string bam_file_name = *bam_it;
    bam_file final_bam_file(bam_file_name, fasta_file_name, ios::in);

    alignment_list al;

    while (final_bam_file.read_alignments(al, false)) {

      for (alignment_list::iterator it = al.begin(); it != al.end(); it++) {
        bam_alignment& a = *it->get();

        if (a.unmapped()) continue;

        // Skip multiply-mapped reads
        uint32_t num_equivalent_alignments;
        if (a.aux_get_i("X1", num_equivalent_alignments)) {
          if (num_equivalent_alignments > 1) continue;
        }

        // Skip junction-side reads
        uint32_t junction_side;
        if (a.aux_get_i("XJ", junction_side)) continue;

        uint32_t read_length = a.cigar_query_length();
        uint32_t query_begin_soft_clipping = a.query_start_1() - 1;
        uint32_t query_end_soft_clipping = read_length - a.query_end_1();

        string seq_id = final_bam_file.target_name(a);
        uint32_t seq_length = seq_lengths.count(seq_id) ? seq_lengths[seq_id] : 0;
        if (seq_length == 0) continue;

        uint32_t ref_start = a.reference_start_1();
        uint32_t ref_end   = a.reference_end_1();

        // Accumulate direction-specific read-through coverage (difference arrays). The read only
        // contributes where it extends at least minimum_clipped_bases past the position.
        vector<int32_t>& cdm = cov_minus_diff[seq_id];
        if (cdm.empty()) cdm.resize(seq_length + 2, 0);
        if (ref_start + minimum_clipped_bases <= ref_end) {
          cdm[ref_start + minimum_clipped_bases] += 1;
          cdm[ref_end + 1] -= 1;
        }

        vector<int32_t>& cdp = cov_plus_diff[seq_id];
        if (cdp.empty()) cdp.resize(seq_length + 2, 0);
        if (ref_end >= ref_start + minimum_clipped_bases) {
          cdp[ref_start] += 1;
          cdp[ref_end - minimum_clipped_bases + 1] -= 1;
        }

        // Direction -1: leading soft-clip at reference_start_1
        if (query_begin_soft_clipping >= minimum_clipped_bases) {
          bool at_boundary = (ref_start <= 1 || ref_start >= seq_length);
          if (!at_boundary) {
            clipped_neg[seq_id][ref_start]++;
            total_clipped_read_ends++;
          }
        }

        // Direction +1: trailing soft-clip at reference_end_1
        if (query_end_soft_clipping >= minimum_clipped_bases) {
          bool at_boundary = (ref_end <= 1 || ref_end >= seq_length);
          if (!at_boundary) {
            clipped_pos[seq_id][ref_end]++;
            total_clipped_read_ends++;
          }
        }

      } // end for each alignment in group
    } // end while read_alignments
  } // end for each bam file

  // Prefix-sum the difference arrays into per-position read-through coverage, and sum the
  // read-through opportunities over non-boundary positions (both directions).
  map<string, vector<uint32_t> > cov_minus;
  map<string, vector<uint32_t> > cov_plus;
  uint64_t total_spanning_read_bases = 0;
  for (map<string, uint32_t>::iterator seq_it = seq_lengths.begin(); seq_it != seq_lengths.end(); seq_it++) {
    const string& seq_id = seq_it->first;
    uint32_t seq_length = seq_it->second;
    if (cov_minus_diff.count(seq_id) == 0) continue; // no reads on this sequence

    vector<int32_t>& cdm = cov_minus_diff[seq_id];
    vector<int32_t>& cdp = cov_plus_diff[seq_id];
    vector<uint32_t>& cm = cov_minus[seq_id];
    vector<uint32_t>& cp = cov_plus[seq_id];
    cm.resize(cdm.size(), 0);
    cp.resize(cdp.size(), 0);

    int32_t running_m = 0, running_p = 0;
    for (uint32_t p = 1; p < cdm.size(); p++) {
      running_m += cdm[p];
      running_p += cdp[p];
      cm[p] = (running_m > 0) ? static_cast<uint32_t>(running_m) : 0;
      cp[p] = (running_p > 0) ? static_cast<uint32_t>(running_p) : 0;
      bool at_boundary = (p <= 1 || (seq_length > 0 && p >= seq_length));
      if (!at_boundary) total_spanning_read_bases += cm[p] + cp[p];
    }
  }

  // Store genome-wide totals in summary.
  // p0 = baseline probability that a read reaching a position is clipped there
  //    = total clip events / (total clip events + total read-through opportunities).
  uint64_t total_opportunities = total_clipped_read_ends + total_spanning_read_bases;
  summary.soft_clipping.total_spanning_read_bases = total_spanning_read_bases;
  summary.soft_clipping.total_clipped_read_ends = total_clipped_read_ends;
  summary.soft_clipping.soft_clipping_rate = (total_opportunities > 0)
    ? static_cast<double>(total_clipped_read_ends) / static_cast<double>(total_opportunities)
    : 0.0;

  // Write per-position counts file. Only positions with clip events are written; total_count is
  // the number of reads tested at that position = clipped reads + read-through reads.
  ofstream out(settings.soft_clipping_counts_file_name.c_str());
  out << "seq_id\tposition\tdirection\tclipped_count\ttotal_count\n";

  for (map<string, map<uint32_t, uint32_t>>::iterator seq_it = clipped_neg.begin(); seq_it != clipped_neg.end(); seq_it++) {
    const string& seq_id = seq_it->first;
    const vector<uint32_t>& cm = cov_minus[seq_id];
    for (map<uint32_t, uint32_t>::iterator pos_it = seq_it->second.begin(); pos_it != seq_it->second.end(); pos_it++) {
      uint32_t p = pos_it->first;
      uint32_t clipped = pos_it->second;
      uint32_t total = clipped + ((p < cm.size()) ? cm[p] : 0);
      out << seq_id << "\t" << p << "\t-1\t" << clipped << "\t" << total << "\n";
    }
  }

  for (map<string, map<uint32_t, uint32_t>>::iterator seq_it = clipped_pos.begin(); seq_it != clipped_pos.end(); seq_it++) {
    const string& seq_id = seq_it->first;
    const vector<uint32_t>& cp = cov_plus[seq_id];
    for (map<uint32_t, uint32_t>::iterator pos_it = seq_it->second.begin(); pos_it != seq_it->second.end(); pos_it++) {
      uint32_t p = pos_it->first;
      uint32_t clipped = pos_it->second;
      uint32_t total = clipped + ((p < cp.size()) ? cp[p] : 0);
      out << seq_id << "\t" << p << "\t1\t" << clipped << "\t" << total << "\n";
    }
  }

  out.close();

  cerr << "  Soft-clipping summary: " << total_clipped_read_ends << " clip events over "
       << total_opportunities << " read opportunities (rate=" << summary.soft_clipping.soft_clipping_rate << ")" << endl;
}

} // namespace breseq

