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


#include "candidate_junctions.h"

#include "calculate_trims.h"
#include "fastq.h"

using namespace std;

namespace breseq {

  /*! Utility for sorting alignments by the number of mismatches
   */
  
  class alignment_score_map_class : public map<bam_alignment*,uint32_t> {
  public:
    inline bool operator() (alignment_list::iterator a1, alignment_list::iterator a2) { return (*this)[a1->get()] > (*this)[a2->get()]; } 
  } alignment_score_map;
  
  
  bool sort_by_alignment_score (counted_ptr<bam_alignment>& a1, counted_ptr<bam_alignment>& a2) { return alignment_score_map[a1.get()] > alignment_score_map[a2.get()]; }  
  
  
  /*! Filter a list of alignments to only those that are eligible for mapping
   *  and return the BEST SCORE, which is currently the length of the aligned region of the read
   *  minus the number of base mismatches and minus the number of indels.
   *
   *  In junction_mode, negative overlap of the junction sequence is also subtracted from the score.
   *    and returns all alignments above the minimum.
   *
   *  @JEB: 2019-12-04 change to only return best-scoring junction alignments!
   *
   *  Otherwise, returns only alignments tied for best.
   */
  
 uint32_t eligible_read_alignments(
                                   const Settings& settings,
                                   const cReferenceSequences& ref_seq_info,
                                   alignment_list& alignments,
                                   bool keep_suboptimal_matches,
                                   int32_t min_match_score
                                   )
  {
    bool verbose = false;
    
    if (alignments.size() == 0) return 0;
    
    // Require read to be mapped and get out of here if it is not!
    for (alignment_list::iterator it = alignments.begin(); it != alignments.end();)
    {
      if (it->get()->unmapped())
        it = alignments.erase(it);
      else
        it++;
    }
    if (alignments.size() == 0) return 0;
    
    // JEB 2018-05-05
    // We must find the best alignment score **BEFORE** applying the
    // guard on minimum length of the read being aligned. Otherwise,
    // we can throw out the match with the best alignment score if it doesn't
    // extend all the way to the end of the read and be left with a worse
    // match that happens to extend further and pass but has many mismatches.
    
    uint32_t read_length = alignments.front()->read_length();
    alignment_score_map.clear();
    for (alignment_list::iterator it = alignments.begin(); it != alignments.end();)
    {
      bam_alignment* ap = it->get(); // we are saving the pointer value as the map key
        
      int32_t as;
      
      
      // IMPORTANT NOTE @JEB 2016-12-09
      //
      // We would like to use the 'AS' output by bowtie2, but it isn't consistent.
      // For equal alignments it still gives different scores even with --ignore-quals
      // at times when it matches on different strands.
      //
      // If this gets ironed out someday...
      //
      
      as = alignment_score(*ap, ref_seq_info);
       
      // We only record positive scores
      as = (as < 0) ? 0 : as;
      ap->aux_set("AS", 'I', sizeof(int32_t), (uint8_t*)&as);

      
      // Only keep ones with a minimum score
      if (as < min_match_score) {
        it = alignments.erase(it);
      }
      else {
        ap->aux_set(kBreseqAlignmentScoreBAMTag, 'I', sizeof(uint32_t), (uint8_t*)&as);
        alignment_score_map[ap] = as;
        it++;
      }
      
    }
    
    alignments.sort(sort_by_alignment_score);
    
    if (verbose)
    {
      for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++)
      {
        bam_alignment& a = *(it->get());
        cerr << a.query_start_1() << "-" << a.query_end_1() << " ";
        cerr << alignment_score_map[it->get()] << "\n";
      }
    }
      
    if (alignments.size() == 0) return 0;
   
    // how many reads share the best score?
    uint32_t last_best(0);
    uint32_t best_score = alignment_score_map[alignments.front().get()];
    
    // Return only the best matches
    // This code must be executed to set kBreseqBestAlignmentScoreBAMTag
    
    for (alignment_list::iterator it = alignments.begin()++; it != alignments.end(); it++)
    {
      bam_alignment* ap = it->get();
      uint8_t is_best_score = (alignment_score_map[ap] == best_score) ? 1 : 0;
      
      ap->aux_set(kBreseqBestAlignmentScoreBAMTag, 'C', sizeof(uint8_t), (uint8_t*)&is_best_score);
      
      if (is_best_score)
        last_best++;
    }
    
    // Default is to truncate only to ties for best
    if (!keep_suboptimal_matches)
      alignments.resize(last_best);
    
    // Require a minimum length of the read to be mapped
    // Must be applied AFTER sorting and accepting based on score
    for (alignment_list::iterator it = alignments.begin(); it != alignments.end();)
    {
      if ( !test_read_alignment_requirements(settings, ref_seq_info, *(it->get())) )
        it = alignments.erase(it);
      else
        it++;
    }
    
    if (verbose)
    {
      for (alignment_list::iterator it = alignments.begin(); it != alignments.end(); it++)
      {
        bam_alignment& a = *(it->get());
        cerr << a.query_start_1() << "-" << a.query_end_1() << endl;
      }
    }
    
    // Note that the score we return is higher for better matches...
    return best_score;
  }
  
  /*! Test the requirements for an alignment to be counted.
   *
   *  Returns whether the alignment passes.
   */
  bool test_read_alignment_requirements(const Settings& settings, const cReferenceSequences& ref_seq_info, const alignment_wrapper& a)
  {
    if (a.unmapped()) return false;

    // Credit soft-clipped bases that extend past a reference boundary:
    // they couldn't align because the reference ran out, not because
    // the read is a poor match.
    uint32_t effective_match_length = a.query_match_length();
    uint32_t seq_length = ref_seq_info[a.reference_target_id()].m_length;
    if (a.reference_start_1() <= 1)
      effective_match_length += a.query_start_1() - 1;
    if (a.reference_end_1() >= seq_length)
      effective_match_length += a.read_length() - a.query_end_1();

    if (effective_match_length < settings.require_match_length)
      return false;

    // Later: Move to being a guard during mutation identification and not here.
    if (a.mapping_quality() < settings.minimum_mapping_quality)
      return false;

    if (effective_match_length < settings.require_match_fraction * static_cast<double>(a.read_length()))
      return false;

    if (settings.maximum_read_mismatches >= 0)
    {
      int32_t mismatches = alignment_mismatches(a, ref_seq_info);
      if (mismatches > settings.maximum_read_mismatches)
        return false;
    }

    return true;
  }

  // Splits unaligned reads out of a SAM file
  void PreprocessAlignments::split_mapped_and_unmapped_alignments(
                                                                    uint32_t fastq_file_index,
                                                                    string fasta_file_name, 
                                                                    string input_sam_file_name, 
                                                                    string mapped_sam_file_name,
                                                                    string unmapped_reads_fastq_file_name
                                                                    )
  {

    bam_file mapped(mapped_sam_file_name, fasta_file_name, ios::out);
    ofstream unmapped(unmapped_reads_fastq_file_name.c_str());

    bam_file in(input_sam_file_name, fasta_file_name, ios::in);

    alignment_list al;
    while (in.read_alignments(al, false)) {
      if (al.front()->unmapped()) {
        unmapped << "@" << al.front()->read_name() << endl << al.front()->read_char_sequence() << endl << "+" << endl << al.front()->read_base_quality_char_string() << endl;
      } else {
        mapped.write_alignments(fastq_file_index, al);
      }
    }
  }
  void PreprocessAlignments::split_matched_alignments(uint32_t fastq_file_index,
                                                      string fasta_file_name, 
                                                      string input_sam_file_name,
                                                      string matched_sam_file_name) 
  {
    bam_file in(input_sam_file_name, fasta_file_name, ios::in);
    bam_file matched(matched_sam_file_name, fasta_file_name, ios::out);

    alignment_list al;
    while (in.read_alignments(al, false)) {
      if (!al.front()->unmapped()) {
        matched.write_alignments(fastq_file_index, al);
      }
    }

    return;
  }
  
  static void bam_to_read_index(const char* qname, int64_t& index_a, int64_t& index_b) {
    const char* colon = strchr(qname, ':');
    ASSERT(colon, "Malformed qname (missing ':'): " + string(qname));
    index_a = 0;
    const char* p = colon + 1;
    while (*p && *p != 'S') index_a = index_a * 10 + (*p++ - '0');
    index_b = 0;
    if (*p == 'S') { ++p; while (*p) index_b = index_b * 10 + (*p++ - '0'); }
  }

  //! Processes one read's worth of alignments (as grouped by bam_file::read_alignments)
  //! for candidate junction preprocessing: writes split-supporting alignments to PSAM, and
  //! writes the best-scoring alignments to BSAM, updating summary statistics.
  //! Returns false if candidate_junction_read_limit was hit (caller should stop).
  bool PreprocessAlignments::preprocess_alignment_list(
                                                        const Settings& settings,
                                                        Summary& summary,
                                                        const cReferenceSequences& ref_seq_info,
                                                        bam_file& BSAM,
                                                        bam_file& PSAM,
                                                        alignment_list& alignments,
                                                        uint32_t& i,
                                                        bool& wrote_to_split
                                                        )
  {
    wrote_to_split = false;

    if (++i % 100000 == 0) {
      ostringstream progress_message;
      progress_message << "    ALIGNED READ:" << setw(12) << right << i;
      print_progress_line(progress_message.str());
    }

    summary.preprocess_alignments.aligned_reads++;
    summary.preprocess_alignments.alignments += alignments.size();

    // for testing...
    if (settings.candidate_junction_read_limit != 0 && i > settings.candidate_junction_read_limit) return false;

    if (alignments.front()->unmapped()) return true;

    // write remaining split-supporting alignments for candidate junction identification
    wrote_to_split = write_junction_candidate_alignments(settings, summary, PSAM, alignments);

    // write best alignments
    int32_t best_score = eligible_read_alignments(settings, ref_seq_info, alignments);
    (void) best_score;
    BSAM.write_alignments(0, alignments, NULL);

    return true;
  }

  // Merges two BAM files (stage1 and stage2 alignments for one read file) so that
  // read order matches the original FASTQ file, writing the result to
  // output_sam_file_name. Requires reads to be renamed as we expect from
  // 01_sequence_conversion: file_num:read_num. If do_preprocess is set, also feeds
  // each mapped read's alignments to preprocess_alignment_list (writing BSAM/PSAM)
  // as part of the same pass.
  void PreprocessAlignments::merge_two_sam_files(
                                                 Settings& settings,
                                                 Summary& summary,
                                                 const cReferenceSequences& ref_seq_info,
                                                 const string& input_sam_file_name_1,
                                                 const string& input_sam_file_name_2,
                                                 const string& output_sam_file_name,
                                                 bool do_preprocess,
                                                 bam_file& BSAM,
                                                 bam_file& PSAM,
                                                 uint32_t& i
                                                 )
  {
    string reference_fasta_file_name = settings.reference_fasta_file_name;

    bam_file in1(input_sam_file_name_1, reference_fasta_file_name, ios_base::in);
    bam_file in2(input_sam_file_name_2, reference_fasta_file_name, ios_base::in);

    samFile* out = hts_open(output_sam_file_name.c_str(), "wb");
    ASSERT(out, "Could not open output BAM file: " + output_sam_file_name);
    ASSERT(sam_hdr_write(out, in1.bam_header) == 0, "Failed to write BAM header to: " + output_sam_file_name);

    alignment_list group1, group2;
    bool not_done_1 = in1.read_alignments(group1, false);
    bool not_done_2 = in2.read_alignments(group2, false);

    int64_t index_1a = -1, index_1b = -1, index_2a = -1, index_2b = -1;

    if (not_done_1) bam_to_read_index(group1.front()->read_name().c_str(), index_1a, index_1b);
    if (not_done_2) bam_to_read_index(group2.front()->read_name().c_str(), index_2a, index_2b);

    bool keep_going = true;
    while (keep_going && (not_done_1 || not_done_2)) {
      bool take_1 = not_done_1 && (
        (index_1a < index_2a) || (index_1a == index_2a && index_1b < index_2b) || !not_done_2);

      alignment_list& group = take_1 ? group1 : group2;
      bam_hdr_t* hdr = take_1 ? in1.bam_header : in2.bam_header;

      if (!group.front()->unmapped()) {
        if (do_preprocess) {
          // Split alignments with their own large indel into JC-candidate-supporting pieces.
          // Splitting reads an alignment's own CIGAR and never mutates group.
          int32_t split_cutoff = settings.junction_indel_split_length;
          split_alignments_on_indels(settings, summary, ref_seq_info, PSAM, split_cutoff, group);
        }
        for (alignment_list::iterator it = group.begin(); it != group.end(); it++) {
          ASSERT(sam_write1(out, hdr, it->get()) >= 0, "Failed to write alignment to: " + output_sam_file_name);
        }
        if (do_preprocess) {
          bool wrote_to_split_unused = false;
          keep_going = preprocess_alignment_list(settings, summary, ref_seq_info, BSAM, PSAM, group, i, wrote_to_split_unused);
        }
      }

      if (take_1) {
        not_done_1 = in1.read_alignments(group1, false);
        if (not_done_1) bam_to_read_index(group1.front()->read_name().c_str(), index_1a, index_1b);
        else { index_1a = index_1b = numeric_limits<int64_t>::max(); }
      } else {
        not_done_2 = in2.read_alignments(group2, false);
        if (not_done_2) bam_to_read_index(group2.front()->read_name().c_str(), index_2a, index_2b);
        else { index_2a = index_2b = numeric_limits<int64_t>::max(); }
      }
    }

    hts_close(out);
  }

  // Reads an already-complete reference_sam_file_name (written directly by a single
  // bowtie2 pass) and feeds each read's alignments to preprocess_alignment_list,
  // writing BSAM/PSAM. Does not rewrite reference_sam_file_name.
  void PreprocessAlignments::preprocess_one_sam_file(
                                                     Settings& settings,
                                                     Summary& summary,
                                                     const cReferenceSequences& ref_seq_info,
                                                     const string& reference_sam_file_name,
                                                     bam_file& BSAM,
                                                     bam_file& PSAM,
                                                     uint32_t& i
                                                     )
  {
    bam_file tam(reference_sam_file_name, settings.reference_fasta_file_name, ios_base::in);

    alignment_list alignments;
    while (tam.read_alignments(alignments, false))
    {
      if (!alignments.front()->unmapped()) {
        int32_t split_cutoff = settings.junction_indel_split_length;
        split_alignments_on_indels(settings, summary, ref_seq_info, PSAM, split_cutoff, alignments);
      }
      bool wrote_to_split_unused = false;
      if (!preprocess_alignment_list(settings, summary, ref_seq_info, BSAM, PSAM, alignments, i, wrote_to_split_unused))
        break;
    }
  }

  // For a matched read pair (R1's and R2's alignment summaries for the same read number),
  // finds the (R1 alignment, R2 alignment) combination sharing a reference sequence (tid)
  // with the smallest outermost-coordinate distance. Orientation letters are assigned by
  // 5'-end order (the mate whose 5' end is lower first), not by R1/R2 identity or leftmost
  // mapped coordinate; "RR" is folded to "FF" since they represent the same relative orientation. Reference
  // sequence circularity is ignored (always a simple linear distance). Returns false if no
  // alignment combination shares a tid (caller should not add a row for this read pair).
  static bool best_pair_orientation_and_distance(
                                                 const vector<AlignmentSummary>& r1_alignments,
                                                 const vector<AlignmentSummary>& r2_alignments,
                                                 string& orientation_out,
                                                 int64_t& distance_out
                                                 )
  {
    bool found = false;
    int64_t best_distance = 0;
    string best_orientation;

    for (vector<AlignmentSummary>::const_iterator a = r1_alignments.begin(); a != r1_alignments.end(); a++) {
      for (vector<AlignmentSummary>::const_iterator b = r2_alignments.begin(); b != r2_alignments.end(); b++) {
        if (a->tid != b->tid) continue;

        // Order by each mate's 5' end (forward -> start_1, reverse -> end_1), NOT by leftmost
        // mapped coordinate. Leftmost coordinate ties for short/overlapping fragments and would
        // mislabel normal FR fragments as RF; the 5' end gives the correct read geometry regardless
        // of overlap.
        const uint32_t a_5p = a->reversed ? a->end_1 : a->start_1;
        const uint32_t b_5p = b->reversed ? b->end_1 : b->start_1;
        const AlignmentSummary& lower = (a_5p <= b_5p) ? *a : *b;
        const AlignmentSummary& higher = (a_5p <= b_5p) ? *b : *a;

        string orientation;
        orientation += lower.reversed ? 'R' : 'F';
        orientation += higher.reversed ? 'R' : 'F';
        if (orientation == "RR") orientation = "FF";

        int64_t distance = static_cast<int64_t>(max(a->end_1, b->end_1)) - static_cast<int64_t>(min(a->start_1, b->start_1));

        if (!found || distance < best_distance) {
          found = true;
          best_distance = distance;
          best_orientation = orientation;
        }
      }
    }

    if (found) {
      orientation_out = best_orientation;
      distance_out = best_distance;
    }
    return found;
  }

  // A single logical read-alignment source, in increasing read-number order, backed by
  // either one already-complete BAM file or a 2-way merge of two BAM files (stage1 and
  // stage2 of the same read file) by read number -- i.e. the same take_1 bookkeeping as
  // merge_two_sam_files, but exposed as a step-able cursor so two of these (one per mate)
  // can be advanced in lockstep by a caller without buffering either file in memory.
  class MergedReadStream {
  public:
    // Two-stage: merge file1 (stage1) and file2 (stage2) by read number.
    MergedReadStream(const string& file1, const string& file2, const string& reference_fasta_file_name)
      : in1(file1, reference_fasta_file_name, ios_base::in),
        in2(file2, reference_fasta_file_name, ios_base::in)
    {
      not_done_1 = in1.read_alignments(group1, false);
      not_done_2 = in2.read_alignments(group2, false);
      if (not_done_1) bam_to_read_index(group1.front()->read_name().c_str(), index_1a, index_1b);
      if (not_done_2) bam_to_read_index(group2.front()->read_name().c_str(), index_2a, index_2b);
    }

    // Single-stage: one already-complete file; in2 stays default-constructed/unopened and
    // not_done_2 stays permanently false, so take_1() always resolves to stream 1.
    MergedReadStream(const string& file1, const string& reference_fasta_file_name)
      : in1(file1, reference_fasta_file_name, ios_base::in)
    {
      not_done_1 = in1.read_alignments(group1, false);
      not_done_2 = false;
      if (not_done_1) bam_to_read_index(group1.front()->read_name().c_str(), index_1a, index_1b);
    }

    bool has_current() const { return not_done_1 || not_done_2; }
    int64_t current_read_number() const { return take_1() ? index_1a : index_2a; }
    alignment_list& current_group() { return take_1() ? group1 : group2; }
    bam_hdr_t* header() { return in1.bam_header; }

    // Two-stage: true when the current group came from the stage1 file (in1); false = stage2
    // (in2). Single-stage: always true (there is no separate stage2 file; not_done_2 is
    // permanently false, so take_1() always resolves to stream 1).
    bool current_is_stage1() const { return take_1(); }

    void advance() {
      if (take_1()) {
        not_done_1 = in1.read_alignments(group1, false);
        if (not_done_1) bam_to_read_index(group1.front()->read_name().c_str(), index_1a, index_1b);
        else { index_1a = index_1b = numeric_limits<int64_t>::max(); }
      } else {
        not_done_2 = in2.read_alignments(group2, false);
        if (not_done_2) bam_to_read_index(group2.front()->read_name().c_str(), index_2a, index_2b);
        else { index_2a = index_2b = numeric_limits<int64_t>::max(); }
      }
    }

  private:
    bool take_1() const {
      return not_done_1 && ((index_1a < index_2a) || (index_1a == index_2a && index_1b < index_2b) || !not_done_2);
    }

    bam_file in1, in2; // in2 default-constructed (unopened) in the 1-file ctor -- never read
    alignment_list group1, group2;
    bool not_done_1, not_done_2;
    int64_t index_1a = -1, index_1b = -1, index_2a = -1, index_2b = -1;
  };

  // Extracts one read group's alignment summaries -- used to snapshot the group BEFORE it
  // is handed to preprocess_alignment_list, which can erase entries from (or fully empty)
  // the list via eligible_read_alignments.
  static vector<AlignmentSummary> summarize_alignment_group(const alignment_list& group)
  {
    vector<AlignmentSummary> summaries;
    summaries.reserve(group.size());
    for (alignment_list::const_iterator it = group.begin(); it != group.end(); it++) {
      summaries.push_back(AlignmentSummary{ (*it)->reference_target_id(), (*it)->reference_start_1(), (*it)->reference_end_1(), (*it)->reversed() });
    }
    return summaries;
  }

  // Synchronously traverses two MergedReadStreams (one per mate) by read number, always
  // advancing whichever side is behind. Writes each mate's current group to its own output
  // BAM (if out1/out2 non-NULL) and preprocesses it into its own PSAM (sharing BSAM and i),
  // then -- for every read number seen as mapped in both mates at the same time -- writes
  // one row to csv with the best (smallest-distance) orientation/distance for that pair.
  void PreprocessAlignments::merge_join_paired_read_streams(
                                             Settings& settings,
                                             Summary& summary,
                                             const cReferenceSequences& ref_seq_info,
                                             MergedReadStream& r1_stream,
                                             MergedReadStream& r2_stream,
                                             samFile* out1,
                                             samFile* out2,
                                             const string& output_sam_file_name_r1,
                                             const string& output_sam_file_name_r2,
                                             map<pair<string, int64_t>, uint32_t>& distance_counts,
                                             bool do_preprocess,
                                             bam_file& BSAM,
                                             bam_file& PSAM1,
                                             bam_file& PSAM2,
                                             uint32_t& i,
                                             ofstream* pair_pos_csv_r1,
                                             ofstream* pair_pos_csv_r2
                                             )
  {
    // Writes one sidecar row: the split read's name, then a 4-tuple (seq_id,F/R,start_1,end_1)
    // per alignment of its (not-split) mate. Consumed by identify_candidate_junctions to pin a
    // redundant junction side to the copy the mate maps concordantly next to.
    auto write_pair_positions_row = [](ofstream* csv, const string& read_name,
                                       const vector<AlignmentSummary>& mate, bam_hdr_t* mate_header) {
      if (!csv) return;
      *csv << read_name;
      for (vector<AlignmentSummary>::const_iterator it = mate.begin(); it != mate.end(); it++) {
        *csv << "," << mate_header->target_name[it->tid]
             << "," << (it->reversed ? 'R' : 'F')
             << "," << it->start_1
             << "," << it->end_1;
      }
      *csv << endl;
    };

    bool keep_going = true;
    while (keep_going && (r1_stream.has_current() || r2_stream.has_current())) {
      bool r1_avail = r1_stream.has_current();
      bool r2_avail = r2_stream.has_current();
      bool matched = r1_avail && r2_avail && (r1_stream.current_read_number() == r2_stream.current_read_number());
      bool process_r1 = r1_avail && (matched || !r2_avail || r1_stream.current_read_number() < r2_stream.current_read_number());
      bool process_r2 = r2_avail && (matched || !r1_avail || r2_stream.current_read_number() < r1_stream.current_read_number());

      bool r1_mapped = false, r2_mapped = false;
      bool r1_is_stage1 = false, r2_is_stage1 = false;
      vector<AlignmentSummary> r1_summary, r2_summary;
      bool r1_wrote_split = false, r2_wrote_split = false;
      string r1_read_name, r2_read_name;

      if (process_r1) {
        alignment_list& group = r1_stream.current_group();
        r1_is_stage1 = r1_stream.current_is_stage1();
        r1_mapped = !group.front()->unmapped();
        if (r1_mapped) {
          r1_read_name = group.front()->read_name();
          if (do_preprocess) {
            int32_t split_cutoff = settings.junction_indel_split_length;
            if (split_alignments_on_indels(settings, summary, ref_seq_info, PSAM1, split_cutoff, group) > 0)
              r1_wrote_split = true;
          }
          if (out1) {
            for (alignment_list::iterator it = group.begin(); it != group.end(); it++) {
              ASSERT(sam_write1(out1, r1_stream.header(), it->get()) >= 0, "Failed to write alignment to: " + output_sam_file_name_r1);
            }
          }
          if (matched) r1_summary = summarize_alignment_group(group);
          if (do_preprocess) {
            bool wrote_to_split = false;
            keep_going = preprocess_alignment_list(settings, summary, ref_seq_info, BSAM, PSAM1, group, i, wrote_to_split) && keep_going;
            if (wrote_to_split) r1_wrote_split = true;
          }
        }
      }

      if (process_r2) {
        alignment_list& group = r2_stream.current_group();
        r2_is_stage1 = r2_stream.current_is_stage1();
        r2_mapped = !group.front()->unmapped();
        if (r2_mapped) {
          r2_read_name = group.front()->read_name();
          if (do_preprocess) {
            int32_t split_cutoff = settings.junction_indel_split_length;
            if (split_alignments_on_indels(settings, summary, ref_seq_info, PSAM2, split_cutoff, group) > 0)
              r2_wrote_split = true;
          }
          if (out2) {
            for (alignment_list::iterator it = group.begin(); it != group.end(); it++) {
              ASSERT(sam_write1(out2, r2_stream.header(), it->get()) >= 0, "Failed to write alignment to: " + output_sam_file_name_r2);
            }
          }
          if (matched) r2_summary = summarize_alignment_group(group);
          if (do_preprocess) {
            bool wrote_to_split = false;
            keep_going = preprocess_alignment_list(settings, summary, ref_seq_info, BSAM, PSAM2, group, i, wrote_to_split) && keep_going;
            if (wrote_to_split) r2_wrote_split = true;
          }
        }
      }

      // Paired-mapping sidecar: for a matched pair where exactly one mate is a split read, record
      // the OTHER (not-split, but possibly multiply-mapped) mate's positions, keyed by the split
      // read's own name so identify_candidate_junctions can look them up when building candidates.
      if (matched) {
        if (r1_wrote_split && !r2_wrote_split && !r2_summary.empty())
          write_pair_positions_row(pair_pos_csv_r1, r1_read_name, r2_summary, r2_stream.header());
        if (r2_wrote_split && !r1_wrote_split && !r1_summary.empty())
          write_pair_positions_row(pair_pos_csv_r2, r2_read_name, r1_summary, r1_stream.header());
      }

      // Only tabulate distance/orientation stats for pairs where both mates' alignments came
      // from stage1 -- stage1's stringent --score-min threshold means a heavily soft-clipped
      // (e.g. junction-spanning) read fails it outright and is only recovered by stage2's
      // relaxed pass, so restricting to stage1 excludes that contamination from the distribution.
      if (matched && r1_mapped && r2_mapped && r1_is_stage1 && r2_is_stage1) {
        string orientation;
        int64_t distance;
        if (best_pair_orientation_and_distance(r1_summary, r2_summary, orientation, distance)) {
          distance_counts[make_pair(orientation, distance)]++;
        }
      }

      if (process_r1) r1_stream.advance();
      if (process_r2) r2_stream.advance();
    }
  }

  // Paired analogue of merge_two_sam_files: synchronously traverses R1's and R2's
  // stage1+stage2 merges by read number, writing each mate's own output BAM exactly as
  // merge_two_sam_files would (sharing BSAM and i), and additionally writes a condensed
  // (orientation,distance,count) histogram CSV -- one row per distinct combination, not
  // per read pair -- named after the read file set's collective base name in
  // settings.data_path (so it survives the run).
  void PreprocessAlignments::merge_two_sets_of_paired_sam_files(
                                                 Settings& settings,
                                                 Summary& summary,
                                                 const cReferenceSequences& ref_seq_info,
                                                 const cReadFileSet& read_file_set,
                                                 const string& input_sam_file_name_1_r1,
                                                 const string& input_sam_file_name_2_r1,
                                                 const string& output_sam_file_name_r1,
                                                 const string& input_sam_file_name_1_r2,
                                                 const string& input_sam_file_name_2_r2,
                                                 const string& output_sam_file_name_r2,
                                                 bool do_preprocess,
                                                 bam_file& BSAM,
                                                 bam_file& PSAM1,
                                                 bam_file& PSAM2,
                                                 uint32_t& i
                                                 )
  {
    MergedReadStream r1_stream(input_sam_file_name_1_r1, input_sam_file_name_2_r1, settings.reference_fasta_file_name);
    MergedReadStream r2_stream(input_sam_file_name_1_r2, input_sam_file_name_2_r2, settings.reference_fasta_file_name);

    samFile* out1 = hts_open(output_sam_file_name_r1.c_str(), "wb");
    ASSERT(out1, "Could not open output BAM file: " + output_sam_file_name_r1);
    ASSERT(sam_hdr_write(out1, r1_stream.header()) == 0, "Failed to write BAM header to: " + output_sam_file_name_r1);

    samFile* out2 = hts_open(output_sam_file_name_r2.c_str(), "wb");
    ASSERT(out2, "Could not open output BAM file: " + output_sam_file_name_r2);
    ASSERT(sam_hdr_write(out2, r2_stream.header()) == 0, "Failed to write BAM header to: " + output_sam_file_name_r2);

    ofstream pair_pos_csv_r1, pair_pos_csv_r2;
    if (do_preprocess && settings.paired_mapping) {
      string pp1 = Settings::file_name(settings.preprocess_junction_split_pair_positions_file_name, "#", read_file_set.m_files[0].m_base_name);
      string pp2 = Settings::file_name(settings.preprocess_junction_split_pair_positions_file_name, "#", read_file_set.m_files[1].m_base_name);
      pair_pos_csv_r1.open(pp1.c_str());
      pair_pos_csv_r2.open(pp2.c_str());
      settings.track_intermediate_file(settings.candidate_junction_done_file_name, pp1);
      settings.track_intermediate_file(settings.candidate_junction_done_file_name, pp2);
    }

    map<pair<string, int64_t>, uint32_t> distance_counts;
    merge_join_paired_read_streams(settings, summary, ref_seq_info, r1_stream, r2_stream,
                                    out1, out2, output_sam_file_name_r1, output_sam_file_name_r2,
                                    distance_counts, do_preprocess, BSAM, PSAM1, PSAM2, i,
                                    pair_pos_csv_r1.is_open() ? &pair_pos_csv_r1 : NULL,
                                    pair_pos_csv_r2.is_open() ? &pair_pos_csv_r2 : NULL);

    hts_close(out1);
    hts_close(out2);

    string csv_file_name = Settings::file_name(settings.paired_mapping_distance_distribution_file_name, "#", read_file_set.m_base_name);
    ofstream csv(csv_file_name.c_str());
    ASSERT(csv.good(), "Could not open file for writing read-pair mapping statistics: " + csv_file_name);
    csv << "orientation,distance,count" << endl;
    for (map<pair<string, int64_t>, uint32_t>::const_iterator it = distance_counts.begin(); it != distance_counts.end(); it++) {
      csv << it->first.first << "," << it->first.second << "," << it->second << endl;
    }
    csv.close();
  }

  // Paired analogue of preprocess_one_sam_file: synchronously traverses R1's and R2's
  // already-complete reference_sam_file_name by read number (sharing BSAM and i, writing
  // nothing new -- single-stage mode never rewrites reference_sam_file_name), and
  // additionally writes the same read-pair mapping statistics CSV as
  // merge_two_sets_of_paired_sam_files. Only called when do_preprocess is true (matching
  // preprocess_one_sam_file's existing guard).
  void PreprocessAlignments::preprocess_one_set_of_paired_sam_files(
                                                 Settings& settings,
                                                 Summary& summary,
                                                 const cReferenceSequences& ref_seq_info,
                                                 const cReadFileSet& read_file_set,
                                                 const string& reference_sam_file_name_r1,
                                                 const string& reference_sam_file_name_r2,
                                                 bam_file& BSAM,
                                                 bam_file& PSAM1,
                                                 bam_file& PSAM2,
                                                 uint32_t& i
                                                 )
  {
    MergedReadStream r1_stream(reference_sam_file_name_r1, settings.reference_fasta_file_name);
    MergedReadStream r2_stream(reference_sam_file_name_r2, settings.reference_fasta_file_name);

    ofstream pair_pos_csv_r1, pair_pos_csv_r2;
    if (settings.paired_mapping) {
      string pp1 = Settings::file_name(settings.preprocess_junction_split_pair_positions_file_name, "#", read_file_set.m_files[0].m_base_name);
      string pp2 = Settings::file_name(settings.preprocess_junction_split_pair_positions_file_name, "#", read_file_set.m_files[1].m_base_name);
      pair_pos_csv_r1.open(pp1.c_str());
      pair_pos_csv_r2.open(pp2.c_str());
      settings.track_intermediate_file(settings.candidate_junction_done_file_name, pp1);
      settings.track_intermediate_file(settings.candidate_junction_done_file_name, pp2);
    }

    map<pair<string, int64_t>, uint32_t> distance_counts;
    merge_join_paired_read_streams(settings, summary, ref_seq_info, r1_stream, r2_stream,
                                    NULL, NULL, reference_sam_file_name_r1, reference_sam_file_name_r2,
                                    distance_counts, true /* only called when do_preprocess is true */,
                                    BSAM, PSAM1, PSAM2, i,
                                    pair_pos_csv_r1.is_open() ? &pair_pos_csv_r1 : NULL,
                                    pair_pos_csv_r2.is_open() ? &pair_pos_csv_r2 : NULL);

    string csv_file_name = Settings::file_name(settings.paired_mapping_distance_distribution_file_name, "#", read_file_set.m_base_name);
    ofstream csv(csv_file_name.c_str());
    ASSERT(csv.good(), "Could not open file for writing read-pair mapping statistics: " + csv_file_name);
    csv << "orientation,distance,count" << endl;
    for (map<pair<string, int64_t>, uint32_t>::const_iterator it = distance_counts.begin(); it != distance_counts.end(); it++) {
      csv << it->first.first << "," << it->first.second << "," << it->second << endl;
    }
    csv.close();
  }

  /*! PreprocessAlignments::merge_sort_and_preprocess_alignments
   *
   *  For each read file, merges the stage1/stage2 alignment BAMs (if two-stage
   *  alignment was used) so that read order matches the original FASTQ file,
   *  writing reference_sam_file_name. If new junction prediction is enabled
   *  (!settings.skip_new_junction_prediction), this is done in the same pass as
   *  preprocessing for candidate junction identification: writing one SAM file
   *  of partial read alignments that could support junctions
   *  (preprocess_junction_split_sam_file_name) and another SAM file of the best
   *  read matches to the reference genome for a preliminary analysis of coverage
   *  (preprocess_junction_best_sam_file_name).
   */
  void PreprocessAlignments::merge_sort_and_preprocess_alignments(Settings& settings, Summary& summary, const cReferenceSequences& ref_seq_info)
  {
    bool do_preprocess = !settings.skip_new_junction_prediction;
    string reference_fasta_file_name = settings.reference_fasta_file_name;

    if (do_preprocess) {
      cout << "Preprocessing alignments." << endl;
      create_path(settings.candidate_junction_path);
    }

    // BSAM includes best matches as they are merged from all alignment files
    bam_file BSAM;
    if (do_preprocess) {
      BSAM.open_write(settings.preprocess_junction_best_sam_file_name, reference_fasta_file_name);
    }

    uint32_t i = 0;
    for (const auto& rfs : settings.read_file_sets)
    {
      if (!rfs.is_paired())
      {
        // ---- UNPAIRED: unchanged from before, just nested under the set loop ----
        for (const auto& read_file : rfs.m_files)   // always exactly 1 file here
        {
          string reference_sam_file_name = Settings::file_name(settings.reference_sam_file_name, "#", read_file.m_base_name);

          bam_file PSAM;
          if (do_preprocess) {
            end_progress_line();
            cerr << "  READ FILE::" << read_file.m_base_name << endl;
            string preprocess_junction_split_sam_file_name = Settings::file_name(settings.preprocess_junction_split_sam_file_name, "#", read_file.m_base_name);
            PSAM.open_write(preprocess_junction_split_sam_file_name, reference_fasta_file_name);
            settings.track_intermediate_file(settings.candidate_junction_done_file_name, preprocess_junction_split_sam_file_name);
          }

          if (settings.bowtie2_stage2.size() != 0) {
            // Two-stage alignment: merge stage1 + stage2 BAMs into reference_sam_file_name
            string stage1_reference_sam_file_name = Settings::file_name(settings.stage1_reference_sam_file_name, "#", read_file.m_base_name);
            string stage2_reference_sam_file_name = Settings::file_name(settings.stage2_reference_sam_file_name, "#", read_file.m_base_name);

            merge_two_sam_files(settings, summary, ref_seq_info, stage1_reference_sam_file_name, stage2_reference_sam_file_name,
                                 reference_sam_file_name, do_preprocess, BSAM, PSAM, i);
          } else if (do_preprocess) {
            // Single-stage alignment: reference_sam_file_name was already written directly by bowtie2
            preprocess_one_sam_file(settings, summary, ref_seq_info, reference_sam_file_name, BSAM, PSAM, i);
          }

          if (do_preprocess) {
            ostringstream progress_message;
            progress_message << "    ALIGNED READ:" << setw(12) << right << i;
            print_progress_line(progress_message.str());
          }

          settings.track_intermediate_file(settings.alignment_correction_done_file_name, reference_sam_file_name);
        }
      }
      else
      {
        // ---- PAIRED: process R1 then R2 (same order/effects as two flat-loop iterations
        // would have had), plus compute/write the read-pair mapping statistics CSV ----
        const cReadFile& read_file_r1 = rfs.m_files[0];
        const cReadFile& read_file_r2 = rfs.m_files[1];

        string reference_sam_file_name_r1 = Settings::file_name(settings.reference_sam_file_name, "#", read_file_r1.m_base_name);
        string reference_sam_file_name_r2 = Settings::file_name(settings.reference_sam_file_name, "#", read_file_r2.m_base_name);

        bam_file PSAM1, PSAM2;
        if (do_preprocess) {
          end_progress_line();
          cerr << "  READ FILE::" << read_file_r1.m_base_name << endl;
          string split1 = Settings::file_name(settings.preprocess_junction_split_sam_file_name, "#", read_file_r1.m_base_name);
          PSAM1.open_write(split1, reference_fasta_file_name);
          settings.track_intermediate_file(settings.candidate_junction_done_file_name, split1);

          cerr << "  READ FILE::" << read_file_r2.m_base_name << endl;
          string split2 = Settings::file_name(settings.preprocess_junction_split_sam_file_name, "#", read_file_r2.m_base_name);
          PSAM2.open_write(split2, reference_fasta_file_name);
          settings.track_intermediate_file(settings.candidate_junction_done_file_name, split2);
        }

        if (settings.bowtie2_stage2.size() != 0) {
          string stage1_r1 = Settings::file_name(settings.stage1_reference_sam_file_name, "#", read_file_r1.m_base_name);
          string stage2_r1 = Settings::file_name(settings.stage2_reference_sam_file_name, "#", read_file_r1.m_base_name);
          string stage1_r2 = Settings::file_name(settings.stage1_reference_sam_file_name, "#", read_file_r2.m_base_name);
          string stage2_r2 = Settings::file_name(settings.stage2_reference_sam_file_name, "#", read_file_r2.m_base_name);

          merge_two_sets_of_paired_sam_files(settings, summary, ref_seq_info, rfs,
                                              stage1_r1, stage2_r1, reference_sam_file_name_r1,
                                              stage1_r2, stage2_r2, reference_sam_file_name_r2,
                                              do_preprocess, BSAM, PSAM1, PSAM2, i);
        } else if (do_preprocess) {
          preprocess_one_set_of_paired_sam_files(settings, summary, ref_seq_info, rfs,
                                                  reference_sam_file_name_r1, reference_sam_file_name_r2,
                                                  BSAM, PSAM1, PSAM2, i);
        }

        if (do_preprocess) {
          ostringstream progress_message;
          progress_message << "    ALIGNED READ:" << setw(12) << right << i;
          print_progress_line(progress_message.str());
        }

        settings.track_intermediate_file(settings.alignment_correction_done_file_name, reference_sam_file_name_r1);
        settings.track_intermediate_file(settings.alignment_correction_done_file_name, reference_sam_file_name_r2);
      }
    }

    if (do_preprocess) {
      end_progress_line();
      cerr << "  Summary... " << endl
           << "  Aligned reads:                          " << setw(12) << right << summary.preprocess_alignments.aligned_reads << endl
           << "  Read alignments:                        " << setw(12) << right << summary.preprocess_alignments.alignments << endl
           << "  Alignments split on indels:             " << setw(12) << right << summary.preprocess_alignments.alignments_split_on_indels << endl
           << "  Reads with alignments split on indels:  " << setw(12) << right << summary.preprocess_alignments.reads_with_alignments_split_on_indels << endl
           << "  Split alignments:                       " << setw(12) << right << summary.preprocess_alignments.split_alignments << endl
           << "  Reads with split alignments:            " << setw(12) << right << summary.preprocess_alignments.reads_with_split_alignments << endl
      ;
    }
  }

  //! Returns true if a D/I op of the given length at ref_pos, with the given actual content
  //! (deleted reference bases for a D, inserted read bases for an I), is a pure length-change
  //! to a reference homopolymer that was already min_homopolymer_length bases or longer: the
  //! indel's own content must be entirely one repeated base, and the existing reference
  //! homopolymer run of that base -- for a deletion, including the deleted span itself, since
  //! those bases are genuinely present in the reference; for an insertion, counting only the
  //! flanking reference bases, since the inserted content isn't reference at all -- must meet
  //! the threshold.
  static bool indel_is_pure_homopolymer_change(
                                                const cAnnotatedSequence& ref_seq,
                                                bool is_deletion,
                                                int32_t ref_pos,
                                                int32_t len,
                                                const string& indel_bases,
                                                int32_t min_homopolymer_length
                                                )
  {
    char c = indel_bases[0];
    for (size_t k = 1; k < indel_bases.size(); k++) {
      if (indel_bases[k] != c) return false;
    }

    int32_t left_pos = ref_pos - 1;
    int32_t left_run = 0;
    while ((left_pos >= 1) && (ref_seq.get_sequence_1(left_pos) == c)) { left_run++; left_pos--; }

    int32_t right_pos = is_deletion ? (ref_pos + len) : ref_pos;
    int32_t right_run = 0;
    while ((right_pos <= static_cast<int32_t>(ref_seq.get_sequence_length())) && (ref_seq.get_sequence_1(right_pos) == c)) { right_run++; right_pos++; }

    int32_t existing_homopolymer_length = left_run + right_run + (is_deletion ? len : 0);
    return existing_homopolymer_length >= min_homopolymer_length;
  }

  /*! Split alignments interrupted by indels into their segments and write to a SAM file.
	 */

  uint32_t PreprocessAlignments::split_alignments_on_indels(const Settings& settings, Summary& summary, const cReferenceSequences& ref_seq_info, bam_file& PSAM, int32_t min_indel_split_len, const alignment_list& alignments)
  {
    //##
    //## @JEB: Note that this may affect the order of alignments in the SAM file. This has
    //## consequences for the AlignmentCorrection step, which may assume highest scoring
    //## alignments are first in the TAM file.
    //##
    //## So only USE the split alignment file for creating a list of candidate junctions.
    //##

    assert(min_indel_split_len >= 0);

    // Keeps track of which to split
    alignment_list untouched_alignments;
    alignment_list split_alignments;

    untouched_alignments.read_base_quality_char_string = alignments.read_base_quality_char_string;
    untouched_alignments.read_base_quality_char_string_reversed = alignments.read_base_quality_char_string_reversed;
    uint32_t alignments_written = 0;

    for(alignment_list::const_iterator it = alignments.begin(); it != alignments.end(); it++) {

      const bam_alignment& a = **it;
      vector<pair<char,uint16_t> > cigar_list = a.cigar_pair_char_op_array();
      const cAnnotatedSequence& ref_seq = ref_seq_info[a.reference_target_id()];
      string read_seq = a.read_char_sequence();

      bool do_split = false;
      int32_t ref_pos = static_cast<int32_t>(a.reference_start_1());
      uint32_t query_pos = a.query_start_1();

      for (size_t i = 0; i < cigar_list.size(); i++) {
        char op = cigar_list[i].first;
        uint16_t len = cigar_list[i].second;
        if (op == 'S') continue;

        if (((op == 'D') || (op == 'I')) && (len >= static_cast<uint16_t>(min_indel_split_len))) {
          bool is_deletion = (op == 'D');
          string indel_bases = is_deletion
              ? ref_seq.get_sequence_1(ref_pos, ref_pos + len - 1)
              : read_seq.substr(query_pos - 1, len);
          if (!indel_is_pure_homopolymer_change(ref_seq, is_deletion, ref_pos, len, indel_bases, min_indel_split_len)) {
            do_split = true;
          }
        }

        if ((op == 'M') || (op == 'D')) ref_pos += len;
        if ((op == 'M') || (op == 'I')) query_pos += len;
      }

      if (do_split) {
        split_alignments.push_back( *it );

      }
      else {
        // Check guard showing that the main alignment is good and thus we don't look at other alignments
        if ( (*it)->query_match_length() >= alignments.front()->read_length() - settings.required_both_unique_length_per_side)
          return 0;

        untouched_alignments.push_back( *it );
      }
    }
    (void) untouched_alignments;

    // Untouched alignments are deliberately NOT written here (unlike the original version of
    // this function): write_junction_candidate_alignments, called later on this same
    // alignments list, is solely responsible for those, so writing them here too would just
    // duplicate them as JC candidates.
    for(alignment_list::iterator it = split_alignments.begin(); it != split_alignments.end(); ++it) {
      PSAM.write_split_alignment(min_indel_split_len, **it, alignments, ref_seq_info);
      alignments_written += 2;
    }

    // record statistics
    if (split_alignments.size()>0) summary.preprocess_alignments.reads_with_alignments_split_on_indels++;
    summary.preprocess_alignments.alignments_split_on_indels +=split_alignments.size() ;

    if (alignments_written > 0) summary.preprocess_alignments.reads_with_split_alignments++;
    summary.preprocess_alignments.split_alignments += alignments_written;

    return alignments_written;
  }

  //! Writes a read's alignments that could support a junction candidate to
  //! PSAM: skipped entirely if one alignment already covers nearly the whole read (nothing
  //! left to usefully pair), otherwise all remaining alignments are written as long as there
  //! is more than one (it takes at least two to make a candidate junction).
  bool PreprocessAlignments::write_junction_candidate_alignments(const Settings& settings, Summary& summary, bam_file& PSAM, const alignment_list& alignments)
  {
    (void) summary;

    for(alignment_list::const_iterator it = alignments.begin(); it != alignments.end(); it++) {
      // Guard showing that a main alignment is already good - it covers so much of the match that it should
      // be considered definitive and not used for junction candidates, so there's nothing left to do.
      // NOTE: The read will still be re-aligned to candidate junctions and could match them better later
      // but it shouldn't help those alternative junctions get seeded => this can lead to odd candidates.
      if (
          ( (*it)->query_match_length() >= alignments.front()->read_length() - settings.required_both_unique_length_per_side) ||
          ( (*it)->query_match_length() >= alignments.front()->read_length() * ( 1 - settings.required_both_unique_length_per_side_fraction) )
         )
        return false;
    }

    // Don't write if there is only one alignment to be written,
    // it takes at least two to make a candidate junction.
    if (alignments.size() > 1) {
      PSAM.write_alignments(0, alignments);
      return true;
    }
    return false;
  }

	/*! Predicts candidate junctions
	 */
  // Parses one #.split_pair_positions.csv line into the split read's name + its mate's alignments.
  // Each line: read_name,seq_id,F|R,start_1,end_1[,seq_id,F|R,start_1,end_1]...
  // Returns false for an empty/malformed line (no complete 4-tuple).
  static bool parse_split_pair_positions_line(const string& line, string& read_name, vector<MateAln>& mates)
  {
    mates.clear();
    if (line.empty()) return false;
    vector<string> f = split(line, ",");
    if (f.size() < 5) return false; // need name + at least one 4-tuple
    read_name = f[0];
    for (size_t k = 1; k + 3 < f.size(); k += 4) {
      MateAln m;
      m.seq_id = f[k];
      m.reversed = (f[k+1] == "R");
      m.start_1 = from_string<int32_t>(f[k+2]);
      m.end_1 = from_string<int32_t>(f[k+3]);
      mates.push_back(m);
    }
    return true;
  }

  // Tests whether the read whose junction-flanking piece is `piece` would form a concordant
  // pair with ANY of its mate's alignments (--paired-mapping). The piece's reference position
  // is projected back to the read's imagined 5' start (extending along the read direction by the
  // read bases 5' of the piece), so the mate distance is a full-fragment length directly
  // comparable to the fitted distance_cutoff -- not a junction-truncated partial. The orientation
  // and distance convention exactly mirrors best_pair_orientation_and_distance().
  static bool mate_concordant_with_piece(bam_alignment& piece,
                                         const cReferenceSequences& ref_seq_info,
                                         const vector<MateAln>& mate_positions,
                                         const string& majority_orientation,
                                         double distance_cutoff)
  {
    uint32_t qs, qe;
    piece.query_stranded_bounds_1(qs, qe);
    int32_t offset_5p = static_cast<int32_t>(qs) - 1; // read bases 5' of this piece on its strand
    int32_t read_len = static_cast<int32_t>(piece.read_length());
    bool rev = piece.reversed();

    // Virtual full-read alignment anchored at the read's imagined 5' start in this locus.
    int32_t v_start, v_end;
    if (!rev) {
      v_start = static_cast<int32_t>(piece.reference_start_1()) - offset_5p;
      v_end = v_start + read_len - 1;
    } else {
      v_end = static_cast<int32_t>(piece.reference_end_1()) + offset_5p;
      v_start = v_end - read_len + 1;
    }

    const string& piece_seq_id = ref_seq_info[piece.reference_target_id()].m_seq_id;

    for (vector<MateAln>::const_iterator m = mate_positions.begin(); m != mate_positions.end(); m++) {
      if (m->seq_id != piece_seq_id) continue;

      int32_t a_5p = rev ? v_end : v_start;
      int32_t b_5p = m->reversed ? m->end_1 : m->start_1;
      bool a_lower = (a_5p <= b_5p);
      bool lower_rev = a_lower ? rev : m->reversed;
      bool higher_rev = a_lower ? m->reversed : rev;

      string orientation;
      orientation += lower_rev ? 'R' : 'F';
      orientation += higher_rev ? 'R' : 'F';
      if (orientation == "RR") orientation = "FF";

      int64_t distance = static_cast<int64_t>(max(v_end, m->end_1)) - static_cast<int64_t>(min(v_start, m->start_1));

      if ((orientation == majority_orientation) && (distance <= static_cast<int64_t>(distance_cutoff)))
        return true;
    }
    return false;
  }

	void CandidateJunctions::identify_candidate_junctions(const Settings& settings, Summary& summary, const cReferenceSequences& ref_seq_info)
	{
		(void)summary; // TO DO: save statistics
    bool verbose = false;
    int32_t read_length_max = summary.sequence_conversion.read_length_max;
    
    /// Load all of the user-defined junctions
    
    map<string,cDiffEntry> user_defined_junctions = load_user_junctions(settings, summary, ref_seq_info);
    
    // hash by junction sequence
    SequenceToKeyToJunctionCandidateMap candidate_junctions;
    
    // shortcut to summary data for this step
    CandidateJunctionSummary& hcs(summary.candidate_junction);
    
    uint32_t i = 0;
    uint64_t passed_alignment_pairs_considered = 0;
    uint64_t total_reads_examined = 0;
    vector<cReadFile> flat_read_files_cj2 = settings.read_file_sets.flat_files();
    for (uint32_t j = 0; j < flat_read_files_cj2.size(); j++)
    {
      cReadFile read_file = flat_read_files_cj2[j];

      string read_file_name = read_file.m_base_name;
      end_progress_line();
      cerr << "  READ FILE::" << read_file_name << endl;

      // Decide which input SAM file we are using...

      string reference_sam_file_name = Settings::file_name(settings.preprocess_junction_split_sam_file_name, "#", flat_read_files_cj2[j].m_base_name);

      bam_file tam(reference_sam_file_name, settings.reference_fasta_file_name, ios_base::in);
      alignment_list alignments;

      // --paired-mapping: stream this read file's mate-position sidecar in lockstep with the split
      // BAM (both are emitted in the same read-number order; the sidecar is a subset holding only
      // reads whose mate was not itself split), and resolve the fitted majority orientation /
      // concordant-distance cutoff for its read-file set, so we can pin a redundant junction side to
      // the copy its mate maps concordantly next to.
      ifstream sidecar_in;
      string pending_name;
      vector<MateAln> pending_mates;
      bool has_pending = false;
      int64_t pending_a = -1, pending_b = -1;
      string majority_orientation;
      double distance_cutoff = 0.0;
      if (settings.paired_mapping) {
        string pp_file = Settings::file_name(settings.preprocess_junction_split_pair_positions_file_name, "#", read_file.m_base_name);
        sidecar_in.open(pp_file.c_str());
        for (const auto& rfs : settings.read_file_sets) {
          if (!rfs.is_paired()) continue;
          bool contains = false;
          for (const auto& rf : rfs.m_files) if (rf.m_base_name == read_file.m_base_name) contains = true;
          if (!contains) continue;
          PairedMappingDistanceDistributionSummaries::const_iterator pit = summary.preliminary_paired_mapping_distance_distribution.find(rfs.m_base_name);
          if (pit != summary.preliminary_paired_mapping_distance_distribution.end()) {
            majority_orientation = pit->second.majority_orientation;
            distance_cutoff = pit->second.distance_cutoff;
          }
          break;
        }
      }

      // Advance the sidecar cursor to the next parseable line (sets has_pending / pending_*).
      auto advance_pending = [&]() {
        string line;
        has_pending = false;
        while (getline(sidecar_in, line)) {
          if (parse_split_pair_positions_line(line, pending_name, pending_mates)) {
            bam_to_read_index(pending_name.c_str(), pending_a, pending_b);
            has_pending = true;
            return;
          }
        }
      };
      if (sidecar_in.is_open()) advance_pending(); // prime the first line

      while (tam.read_alignments(alignments, false))
      {
        if (alignments.size() == 0)
          break;

        if (++i % 10000 == 0) {
          ostringstream progress_message;
          progress_message << "    ALIGNED READ:" << setw(12) << right << i << " CANDIDATE JUNCTIONS:" << setw(12) << right << candidate_junctions.size();
          print_progress_line(progress_message.str());
        }

        // for testing...
        if (settings.candidate_junction_read_limit != 0 && i > settings.candidate_junction_read_limit)
          break;

        // Advance the sidecar cursor in lockstep and look up this split read's mate positions
        // (present only when its mate was not itself split, so the sidecar is a subset of these
        // reads). Advance past any sidecar reads that precede the current one, then match by name.
        const vector<MateAln>* mate_positions = NULL;
        if (settings.paired_mapping && !majority_orientation.empty() && (distance_cutoff > 0.0) && has_pending) {
          int64_t cur_a, cur_b;
          bam_to_read_index(alignments.front()->read_name().c_str(), cur_a, cur_b);
          while (has_pending && ((pending_a < cur_a) || ((pending_a == cur_a) && (pending_b < cur_b))))
            advance_pending();
          if (has_pending && (pending_name == alignments.front()->read_name()))
            mate_positions = &pending_mates;
        }

        // pass back how many were considered
        passed_alignment_pairs_considered += alignments_to_candidate_junctions(settings, summary, ref_seq_info, candidate_junctions, alignments, mate_positions, majority_orientation, distance_cutoff);
        
        if ((settings.maximum_junction_sequence_passed_alignment_pairs_to_consider != 0) && (passed_alignment_pairs_considered >= settings.maximum_junction_sequence_passed_alignment_pairs_to_consider))
          break;
      }

      {
        ostringstream progress_message;
        progress_message << "    ALIGNED READ:" << setw(12) << right << i << " CANDIDATE JUNCTIONS:" << setw(12) << right << candidate_junctions.size();
        print_progress_line(progress_message.str());
      }

      if ((settings.maximum_junction_sequence_passed_alignment_pairs_to_consider != 0) && (passed_alignment_pairs_considered >= settings.maximum_junction_sequence_passed_alignment_pairs_to_consider)) {
        break;
      }
      
    }
    
    summary.candidate_junction.passed_alignment_pairs_considered = passed_alignment_pairs_considered;
    end_progress_line();
    cerr << "  Passed alignment pairs examined: " << passed_alignment_pairs_considered << endl;
    if ( (settings.maximum_junction_sequence_passed_alignment_pairs_to_consider != 0) && (passed_alignment_pairs_considered >= settings.maximum_junction_sequence_passed_alignment_pairs_to_consider) ) {
      cerr << "  WARNING: Reached limit of " << settings.maximum_junction_sequence_passed_alignment_pairs_to_consider << " passed alignment pairs." << endl;
      cerr << "  Specify a greater value for --junction-alignment-pair-limit for more thorough junction prediction." << endl;
    }
    ///
    // --paired-mapping: a junction side is reassigned to a specific repeat copy (marked unique)
    // only if at least settings.concordant_pairs_to_make_unique concordant read pairs agree on it.
    // Convert the accumulated per-side tally into the pair_unique flag the merge below acts on.
    // (With the feature off every count is 0, so pair_unique stays false and merging is unchanged.)
    ///
    if (settings.paired_mapping) {
      for (SequenceToKeyToJunctionCandidateMap::iterator seq_it = candidate_junctions.begin(); seq_it != candidate_junctions.end(); seq_it++)
        for (KeyToJunctionCandidateMap::iterator key_it = seq_it->second.begin(); key_it != seq_it->second.end(); key_it++)
          for (int32_t s = 0; s < 2; s++)
            key_it->second->sides[s].pair_unique =
              (static_cast<uint32_t>(key_it->second->sides[s].pair_unique_count) >= settings.concordant_pairs_to_make_unique);
    }

    ///
		// Merge all junctions with the same exact sequence
    //   * They are hashed together for speed in this comparison
		////

    // New list consisting of merged junction candidates
    list<JunctionCandidatePtr> junction_candidate_list;
    
		for (SequenceToKeyToJunctionCandidateMap::iterator same_seq_it = candidate_junctions.begin(); same_seq_it != candidate_junctions.end(); same_seq_it++)
    {
      KeyToJunctionCandidateMap& kcjm = same_seq_it->second;
      
      KeyToJunctionCandidateMap::iterator it = kcjm.begin();
      
      // compare each pair 
      JunctionCandidatePtr one_junction = it->second;

      for (++it; it != kcjm.end(); ++it)
      {
        JunctionCandidatePtr *jcp1 = &one_junction;
        JunctionCandidatePtr *jcp2 = &it->second;

        merge_candidate_junctions(jcp1, jcp2);
        one_junction = *jcp1;
      }
      
      if (one_junction->pos_hash_score() >= settings.minimum_candidate_junction_pos_hash_score) junction_candidate_list.push_back(one_junction);
    }
    
    ////
		// Second round of merging, attempt to merge each item in the list with _every_ other one
    //   * This takes care of cases where one is a subsequence of another - we keep the shortest one
		////
    
    for (list<JunctionCandidatePtr>::iterator it1 = junction_candidate_list.begin(); it1 != junction_candidate_list.end();)
    {
      list<JunctionCandidatePtr>::iterator it2 = it1;
      bool deleted_it1 = false;
      
      for (++it2; it2 != junction_candidate_list.end(); )
      {
        bool deleted_it2 = false;
        JunctionCandidatePtr *jcp1 = &(*it1);
        JunctionCandidatePtr *jcp2 = &(*it2);
        
        if (verbose)
        {
          cout << "Comparing:" << endl;
          cout << (*jcp1)->junction_key() << endl;
          cout << (*jcp1)->sequence << endl;
          cout << (*jcp2)->junction_key() << endl;
          cout << (*jcp2)->sequence << endl;
        }
        
        
        bool merged = merge_candidate_junctions(jcp1, jcp2);
        if (merged)
        {
          if (verbose) cout << "MERGED" << endl;
          
          if (*jcp1 == *it1)
          {
            it2 = junction_candidate_list.erase(it2);
            deleted_it2 = true;
          }
          else
          {
            it1 = junction_candidate_list.erase(it1);
            deleted_it1 = true;
            break;
          }
        }
        
        if (!deleted_it2) it2++;
      }
      if (!deleted_it1) it1++;
    }

		////
		//  Combine hash into a list, retaining only one item (the best pos_hash score) for each unique sequence 
    //    (and also its reverse complement)
		////
    
		PosHashScoreDistribution observed_pos_hash_score_distribution;
    
		vector<JunctionCandidate> combined_candidate_junctions;
    
		map<string, int32_t> handled_seq;
		map<string, JunctionCandidate> ids_to_print;
    
		// Map is pre-sorted by CandidateJunction::Sorter
		for (list<JunctionCandidatePtr>::iterator it = junction_candidate_list.begin(); it != junction_candidate_list.end(); ++it)
		{
      JunctionCandidatePtr& jcp = *it;
      JunctionCandidate& best_candidate_junction = *jcp;
      
			string junction_seq = best_candidate_junction.sequence;
      string rc_junction_seq = best_candidate_junction.reverse_complement_sequence;

      if (verbose) cout << "Handling " << junction_seq << endl;
      
			// These shouldn't be necessary checks, but keeping to detect unintended errors
      
      // We already handled the sequence or its reverse complement, skip.
      if (handled_seq.count(junction_seq)) continue;
      if (handled_seq.count(rc_junction_seq)) continue;

      string junction_id = best_candidate_junction.junction_key();
            
			handled_seq[junction_seq]++;
			handled_seq[rc_junction_seq]++;
            
			// Save the score in the distribution
			observed_pos_hash_score_distribution.add_score(best_candidate_junction.pos_hash_score());
      
			// Check minimum requirements
			if (best_candidate_junction.pos_hash_score() < settings.minimum_candidate_junction_pos_hash_score)
				continue;
      
      if (verbose) cout << "  best candidate junction:" << endl << best_candidate_junction.junction_key() << endl;
      
			// Make sure it isn't a duplicate junction id -- this should NEVER happen and causes downstream problem.
			// ---> Begin sanity check
			if (ids_to_print.count(junction_id))
			{
        JunctionCandidate& ccj = ids_to_print[junction_id];
        
				cout << "Attempt to create junction candidate with duplicate id: " << junction_id << endl;
        
				cout << "==Existing junction==" << endl;      
        cout << "  id: " << ccj.junction_key() << endl;
        cout << "  pos_hash_score: " << ccj.pos_hash_score() << endl;
        cout << "  pos_hash_score: " << ccj.pos_hash_score() << endl;
        cout << "  seq: " << ccj.sequence << endl;
        cout << "  rc_seq: " << ccj.reverse_complement_sequence << endl;   
        
				cout << "==New junction==" << endl;
        cout << "  id: " << best_candidate_junction.junction_key() << endl;
        cout << "  pos_hash_score: " << best_candidate_junction.pos_hash_score() << endl;
        cout << "  pos_hash_score: " << best_candidate_junction.pos_hash_score() << endl;
        cout << "  seq: " << best_candidate_junction.sequence << endl;
        cout << "  rc_seq: " << best_candidate_junction.reverse_complement_sequence << endl;  
        
				assert (best_candidate_junction.sequence == ids_to_print[junction_id].sequence);
				exit(-1);
			}
      // <--- End sanity check

			ids_to_print[best_candidate_junction.junction_key()] = best_candidate_junction;
      
			combined_candidate_junctions.push_back(best_candidate_junction);
		}
    
		sort(combined_candidate_junctions.begin(), combined_candidate_junctions.end(), JunctionCandidate::sort_by_scores_and_seq_length);
    
    if (verbose)
    {
      for (vector<JunctionCandidate>::iterator it=combined_candidate_junctions.begin(); it < combined_candidate_junctions.end(); it++ )
      {
        cout << "ID: " << it->junction_key() << endl;
        cout << "  pos_hash_score: " << it->pos_hash_score() << endl;
        cout << "  seq: " << it->sequence << endl;
        cout << "  rc_seq: " << it->reverse_complement_sequence << endl;
      }
    }
    
		///
		// Limit the number of candidate junctions that we print by:
		//   (1) A maximum number of candidate junctions
		//   (2) A maximum length of the sequences in candidate junctions
    //   (3) But take at least some minimum despite these
		///
    
		cerr << "  Taking top candidate junctions..." << endl;
    
		// adding up the lengths might be too time-consuming to be worth it...
		int32_t total_cumulative_cj_length = 0;
		int32_t total_candidate_junction_number = combined_candidate_junctions.size();
		for (uint32_t j = 0; j < combined_candidate_junctions.size(); j++)
			total_cumulative_cj_length += combined_candidate_junctions[j].sequence.size();
    
		uint32_t cumulative_cj_length = 0;
		int32_t lowest_accepted_pos_hash_score = 0;
    
		// Right now we limit the candidate junctions to have a length no longer than the reference sequence times some factor.
		uint32_t cj_length_limit = static_cast<uint32_t>(summary.sequence_conversion.total_reference_sequence_length * settings.maximum_candidate_junction_length_factor);
		uint32_t maximum_candidate_junctions = settings.maximum_candidate_junctions;
		uint32_t minimum_candidate_junctions = settings.minimum_candidate_junctions;
    
		fprintf(stderr, "  Minimum number to keep: %7d \n", minimum_candidate_junctions);
    if (maximum_candidate_junctions==0)
      fprintf(stderr, "  Maximum number to keep: no limit \n");
    else
      fprintf(stderr, "  Maximum number to keep: %7d \n", maximum_candidate_junctions);
    
    if (cj_length_limit==0)
      fprintf(stderr, "  Maximum length to keep: no limit \n");
    else
      fprintf(stderr, "  Maximum length to keep: %7d bases\n", cj_length_limit);
    
		cerr << "    Initial: Number = " << total_candidate_junction_number << ", Cumulative Length = " << total_cumulative_cj_length << " bases" << endl;
    
		if (combined_candidate_junctions.size() > 0)
		{
			vector<JunctionCandidate> remaining_ids;
			vector<JunctionCandidate> list_in_waiting;
			int32_t add_cj_length = 0;

			i = 0;
			uint32_t current_pos_hash_score = combined_candidate_junctions[i].pos_hash_score();
      
			// Check to make sure that adding everything from the last iteration doesn't put us over any limits...
			uint32_t new_number = remaining_ids.size() + list_in_waiting.size();
			uint32_t new_length = cumulative_cj_length + add_cj_length;
			while (	( new_number <= minimum_candidate_junctions ) || ( ((cj_length_limit == 0) || (new_length <= cj_length_limit)) && ((maximum_candidate_junctions == 0) || (new_number <= maximum_candidate_junctions))) )
			{
				// OK, add everything from the last iteration
				cumulative_cj_length += add_cj_length;
				remaining_ids.reserve(remaining_ids.size() + list_in_waiting.size());
				remaining_ids.insert(remaining_ids.end(), list_in_waiting.begin(), list_in_waiting.end());
        
				lowest_accepted_pos_hash_score = current_pos_hash_score;
        
				// Zero out what we will add
				add_cj_length = 0;
				list_in_waiting.clear();

				// Check to make sure we haven't exhausted the list
				if (i >= combined_candidate_junctions.size()) break;
        
        // Grab the next chunk with the same score
				current_pos_hash_score = combined_candidate_junctions[i].pos_hash_score();
				while (
               (i < combined_candidate_junctions.size())
               && (combined_candidate_junctions[i].pos_hash_score() == current_pos_hash_score)
               )
				{
					JunctionCandidate c = combined_candidate_junctions[i];

          list_in_waiting.push_back(c);
          add_cj_length += c.sequence.size();
					i++;
				}
        
				new_number = remaining_ids.size() + list_in_waiting.size();
				new_length = cumulative_cj_length + add_cj_length;
        
        fprintf(stderr, "      Testing Pos Hash Score = %4d, Num = %6d, Length = %6d\n", current_pos_hash_score, int32_t(list_in_waiting.size()), add_cj_length);
			}
			combined_candidate_junctions = remaining_ids;
		}
    
		int32_t accepted_candidate_junction_number = combined_candidate_junctions.size();
		cerr << "    Accepted: Number = " << accepted_candidate_junction_number << ", Pos Hash Score >= " << lowest_accepted_pos_hash_score << ", Cumulative Length = " << cumulative_cj_length << " bases" << endl;
    
		// Save summary statistics
		hcs.total_number = total_candidate_junction_number;
		hcs.total_length = total_cumulative_cj_length;
    
		hcs.accepted_number = accepted_candidate_junction_number;
		hcs.accepted_length = cumulative_cj_length;
		hcs.accepted_pos_hash_score_cutoff = lowest_accepted_pos_hash_score;
		hcs.pos_hash_score_distribution = observed_pos_hash_score_distribution;
    
    ///
    // Mark junctions as user defined and append ones that are only user-defined
    ///
    for (uint32_t j = 0; j < combined_candidate_junctions.size(); j++)
		{
			JunctionCandidate& junction = combined_candidate_junctions[j];
      string junction_key = junction.junction_key(false); // false here is important for not including redundant tags
      if (user_defined_junctions.count(junction_key)) {
        cout << "Found as user junction: " << junction.junction_key() << endl;
        junction.user_defined = true;
        user_defined_junctions.erase(junction_key);
      }
		}
    for (map<string,cDiffEntry>::iterator it = user_defined_junctions.begin(); it != user_defined_junctions.end(); it++) {
      cDiffEntry& user_junction = it->second;
      JunctionInfo user_junction_info(user_junction);
      user_junction_info.user_defined = true;
      string junction_sequence = construct_junction_sequence(ref_seq_info, user_junction, read_length_max, false);
      JunctionCandidate new_jc(user_junction_info, junction_sequence);
      combined_candidate_junctions.push_back(new_jc);
    }
    
		///
		// Print out the candidate junctions, sorted by the lower coordinate, higher coord, then number
		///
    
		sort(combined_candidate_junctions.begin(), combined_candidate_junctions.end(), JunctionCandidate::sort_by_ref_seq_coord);

    // Load reference genome trims, so we can reject candidate junctions whose non-redundant
    // side(s) or own junction sequence are guaranteed to have zero confidently-matched register
    // positions once trimmed -- these would otherwise cause a division by zero later, when
    // assign_one_junction_read_counts() normalizes read counts by possible_overlap_registers.
    // See junction_read_counter::count_confident_overlap_registers() in resolve_alignments.cpp
    // for the same math applied to real reads later in the pipeline.
    map<string, SequenceTrims> reference_trims_by_seq_id;
    for (cReferenceSequences::const_iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); it++) {
      string this_file_name = Settings::file_name(settings.reference_trim_file_name, "@", it->m_seq_id);
      reference_trims_by_seq_id[it->m_seq_id].ReadFile(this_file_name, it->m_length);
    }
    uint32_t read_length_avg = static_cast<uint32_t>(round(summary.sequence_conversion.read_length_avg));
    uint32_t num_junctions_rejected_zero_registers = 0;

    cFastaFile out(settings.candidate_junction_fasta_file_name, ios_base::out);
    ofstream detailed;
    if (settings.junction_debug) {
      detailed.open(settings.candidate_junction_detailed_file_name.c_str());
    }

		for (uint32_t j = 0; j < combined_candidate_junctions.size(); j++) {

			JunctionCandidate& junction = combined_candidate_junctions[j];

      // Would a non-redundant side, or the junction itself, have zero confidently-matched
      // register positions after trimming? If so, this junction is guaranteed to produce a
      // zero possible_overlap_registers denominator downstream -- reject it now, before any
      // reads get remapped to it. The single-point window used per side (rather than the full
      // corrected window resolve_alignments.cpp computes later, which isn't determinable this
      // early) is a safe upper bound: every later correction (overlap-split, continuation,
      // minimum-side-match) only widens that window, so if even this narrower check is already
      // zero, the real one downstream is mathematically guaranteed to be zero too.
      bool reject_zero_registers = false;
      for (int32_t side_index = 0; (side_index < 2) && !reject_zero_registers; side_index++) {
        JunctionSide& side = junction.sides[side_index];
        if (side.redundant) continue; // redundant sides don't factor into the denominator
        const cAnnotatedSequence& ref_seq = ref_seq_info[side.seq_id];
        uint32_t seq_length = static_cast<uint32_t>(ref_seq.m_length);

        // For a circular replicon, a read can wrap around the origin -- our (linear) trim
        // model has no way to represent that, so near either end it would underestimate the
        // register count (potentially to zero) even though real reads can and do land there.
        // Skip the check rather than risk a false rejection in that region.
        bool near_circular_boundary = ref_seq.is_circular()
          && ((static_cast<uint32_t>(side.position) <= read_length_avg)
              || (static_cast<uint32_t>(side.position) > seq_length - read_length_avg));
        if (near_circular_boundary) continue;

        if (count_confident_overlap_registers(reference_trims_by_seq_id[side.seq_id], seq_length, side.position, side.position, read_length_avg) == 0) {
          reject_zero_registers = true;
        }
      }
      if (!reject_zero_registers) {
        SequenceTrims junction_trims(junction.sequence);
        int32_t junction_window_start = junction.flanking_left;
        int32_t junction_window_end = junction.flanking_left + abs(junction.alignment_overlap) + 1;
        if (count_confident_overlap_registers(junction_trims, static_cast<uint32_t>(junction.sequence.size()), junction_window_start, junction_window_end, read_length_avg) == 0) {
          reject_zero_registers = true;
        }
      }

      if (reject_zero_registers) {
        num_junctions_rejected_zero_registers++;
        if (settings.junction_debug) {
          detailed << junction.junction_key() << "\tREJECTED: zero possible overlap registers after trimming" << endl;
        }
        continue;
      }

      cFastaSequence seq(junction.junction_key(), "", junction.sequence); //= { junction.junction_key(), "", junction.sequence };
			out.write_sequence(seq);

      // write to detailed file
      if (settings.junction_debug) {
      detailed << seq;
        for (vector<JunctionCandidatePtr>::iterator it = junction.merged_from.begin(); it != junction.merged_from.end(); it++) {
          JunctionCandidate& merged_junction = **it;
          detailed << merged_junction.junction_key() << "\t" << merged_junction.sequence << "\n";
        }
        detailed << "======" << endl;
      }
		}
		out.close();

    if (num_junctions_rejected_zero_registers > 0) {
      cerr << "  Rejected " << num_junctions_rejected_zero_registers << " candidate junction(s) with zero possible overlap registers after trimming (would cause division by zero in frequency calculation)." << endl;
    }

		summary.candidate_junction = hcs;
	}
  
  map<string,cDiffEntry> 
  CandidateJunctions::load_user_junctions(
                                              const Settings& settings,
                                              const Summary& summary,
                                              const cReferenceSequences& ref_seq_info
                                              )
  {
    map<string,cDiffEntry> user_defined_junctions;
   
    // File must exist for us to process
    if (settings.user_evidence_genome_diff_file_name == "")
      return user_defined_junctions;
    
    cGenomeDiff gd(settings.user_evidence_genome_diff_file_name);

    int32_t read_length_max = summary.sequence_conversion.read_length_max;
    
    diff_entry_list_t _entry_list = gd.get_list(make_vector<gd_entry_type>(JC));
    for (diff_entry_list_t::iterator it = _entry_list.begin(); it != _entry_list.end(); it++)
    {
      cDiffEntry& user_junction = **it;

      // debug
      /*
      if (user_junction[SIDE_1_POSITION] == "588495") {
        cout << "problem" << endl;
      }
       */
      //cout << user_junction.as_string() << endl;
      
      // set initial flanking lengths, these may be reduced by construct_junction_sequence
      user_junction["flanking_left"] = to_string<int32_t>(read_length_max);
      user_junction["flanking_right"] = to_string<int32_t>(read_length_max);
      

      // Fix the overlap...
      normalize_junction_overlap(ref_seq_info, user_junction);
      
      // @JEB 2018-02-09 Construct the junction sequence here purely to correct the
      // 'flanking_left' and 'flanking_right' fields so that we will have the correct junction key
      // for comparing to the junctions that are predicted in this sample during merging.
      // Failure to do so can result in duplicate junctions. Must do this after fixing overlap!
      string throwaway_sequence = construct_junction_sequence(ref_seq_info, user_junction, read_length_max);
      
      JunctionInfo junction_info(user_junction);
      user_defined_junctions[junction_info.junction_key(false)] = user_junction;
      
      //cout << user_junction.as_string() << endl;
    }
    return user_defined_junctions; 
  }

  
// CandidateJunctions::merge_candidate_junctions
//
//   Determined whether two junctions are equivalent (i.e., one is a subsequence of the other)
//   If they are, it merges theminto one, according to these criteria: 
//      1) Shorter one if they differ in length (this implies more overlap)
//      2) Favoring one with the two sides on the same reference sequence
//      3) Favoring the one with the closest reference coordinates
  
  bool CandidateJunctions::merge_candidate_junctions(JunctionCandidatePtr*& jcp1, JunctionCandidatePtr*& jcp2)
  {
    bool verbose = false;
    
    
    uint32_t merged_strand = 0; // 0 for not merged, 1 for sequence, 2 for reverse complement of sequence
    
    JunctionCandidate& jc1 = **jcp1;
    JunctionCandidate& jc2 = **jcp2;

    // Determine whether one is a subsequence of the other (including on the opposite strand)
    if (jc1.sequence.size() > jc2.sequence.size())
    {
      if (jc1.sequence.find(jc2.sequence) != string::npos)
        merged_strand = 1;
      else if (jc1.reverse_complement_sequence.find(jc2.sequence) != string::npos )
        merged_strand = 2;
    }
    else
    {
      if (jc2.sequence.find(jc1.sequence) != string::npos)
        merged_strand = 1;
      else if (jc2.sequence.find(jc1.reverse_complement_sequence) != string::npos)
        merged_strand = 2;
    }
    
    if (!merged_strand) return false;
    
    if (verbose) cout << ((merged_strand == 1) ? "Merged same strand" : "Merged reverse complement") << endl;
    
    JunctionCandidatePtr* merge_into_p = NULL;
    JunctionCandidatePtr* merge_from_p = NULL;
    
    // this is a rather complicated compare function to favor shorter sequences and those with close coords

    // This < comparison of junctions favors...
    //   1) the shorter sequence
    //   2) both sides on the same reference sequence
    //   3) the smallest coordinate on side 1

    // --paired-mapping: if exactly one candidate has a pair-confirmed (pinned) side, make it the
    // representative so its pinned coordinate -- not the arbitrary smallest-coord tie-break --
    // becomes the junction's name/coordinate. pair_unique is only ever set under --paired-mapping,
    // so with the feature off both flags are false and this reduces to the original jc2 < jc1.
    bool jc1_pinned = jc1.sides[0].pair_unique || jc1.sides[1].pair_unique;
    bool jc2_pinned = jc2.sides[0].pair_unique || jc2.sides[1].pair_unique;

    if (jc1_pinned != jc2_pinned)
    {
      if (jc1_pinned) { merge_into_p = jcp1; merge_from_p = jcp2; }
      else            { merge_into_p = jcp2; merge_from_p = jcp1; }
    }
    else if ( jc2 < jc1 )
    {
      merge_into_p = jcp1;
      merge_from_p = jcp2;
    }
    else
    {
      merge_into_p = jcp2;
      merge_from_p = jcp1;
    }
    
    JunctionCandidate& merge_into = **merge_into_p;
    JunctionCandidate& merge_from = **merge_from_p;
    
    // Add the merged junctions to our list
    if (merge_into.merged_from.size() == 0) merge_into.merged_from.push_back(*merge_into_p);
    merge_into.merged_from.push_back(*merge_from_p);
    
    if (verbose) cout << "Merging into:" << merge_into.junction_key() << endl;
    if (verbose) cout << merge_into.sequence << endl;
    if (verbose) cout << "Merging from:" << merge_from.junction_key() << endl;
    if (verbose) cout << merge_from.sequence << endl;
    
    // Carry over redundancy that was already assigned (based on different coord matching, but sequence being the same)
    if (merged_strand == 2)
    {
      merge_into.sides[0].redundant = merge_from.sides[1].redundant && merge_into.sides[0].redundant;
      merge_into.sides[1].redundant = merge_from.sides[0].redundant && merge_into.sides[1].redundant;
    }
    else // merged == 1
    {
      merge_into.sides[0].redundant = merge_from.sides[0].redundant && merge_into.sides[0].redundant;
      merge_into.sides[1].redundant = merge_from.sides[1].redundant && merge_into.sides[1].redundant;
    }
    
    // If one of the sides is not equivalent (in terms of reference coordinate and strand), 
    // then mark that side as redundant
    
    if (verbose) cout << "Merging into carryover:" << merge_into.junction_key() << endl;
    for (uint32_t into_side = 0; into_side < 2; into_side++) 
    {
      bool found = false; // whether we found an identical coordinate
      
      for (uint32_t from_side = 0; from_side < 2; from_side++) 
      {
        if (merge_into.sides[into_side] == merge_from.sides[from_side])
          found = true;
        
        // There may sometimes be equivalent sides due to overlap with different descriptions
        // NC_005966__3079290__-1__NC_005966__3079300__1__0____101__101
        // NC_005966__3079295__-1__NC_005966__3079300__1__5____101__101
        
        // it would really be better to normalize both junctions here if that code were reliable
        
        if ((merge_into.alignment_overlap >= 0) && (merge_from.alignment_overlap >= 0)) {
          
          JunctionSide shifted_merge_into(merge_into.sides[into_side]);
          JunctionSide shifted_merge_from(merge_from.sides[from_side]);
          
          shifted_merge_into.position += merge_into.alignment_overlap * shifted_merge_into.strand;
          shifted_merge_from.position += merge_from.alignment_overlap * shifted_merge_from.strand;

          
          if (shifted_merge_into == shifted_merge_from)
            found = true;
        }
        
      }
      
      // we did not find an equivalent side, meaning this side must have multiple descriptions == redundant
      if (!found)
      {
        // --paired-mapping: a pair-confirmed side is pinned to a specific copy, so keep it unique
        // instead of marking it redundant -- UNLESS the merge_from candidate's corresponding side
        // is ALSO pinned but to a different coordinate (two reads disagree on the copy), in which
        // case the evidence is contradictory and we fall back to redundant. pair_unique is only set
        // under --paired-mapping, so this whole guard is inert with the feature off.
        bool keep_unique = merge_into.sides[into_side].pair_unique;
        if (keep_unique)
        {
          uint32_t corresponding_from_side = (merged_strand == 2) ? (1 - into_side) : into_side;
          if (merge_from.sides[corresponding_from_side].pair_unique &&
              (merge_from.sides[corresponding_from_side].pair_confirmed_position != merge_into.sides[into_side].pair_confirmed_position))
            keep_unique = false;
        }

        if (!keep_unique)
        {
          merge_into.sides[into_side].redundant = true;
          merge_into.sides[into_side].pair_unique = false;
          merge_into.sides[into_side].pair_confirmed_position = 0;
          if (verbose) cout << "Marking side " << into_side << " as redundant." << endl;
        }
      }
    }
    
    jcp1 = merge_into_p;
    jcp2 = merge_from_p;

    return true;
  }


  /*! Computes reference-side breakpoint geometry for a pair of partial alignments of the
   *  same read, correcting for mismatches near the junction, for use by candidate-junction
   *  construction (CandidateJunctions::alignment_pair_to_candidate_junction).
   *
   *  a1/a2 are reordered internally (via the q1/q2 pointer members, not by copying or
   *  swapping the alignment objects) so that q1 is read-earlier; the caller's a1/a2
   *  themselves are never modified.
   */
  AlignmentPairGeometry::AlignmentPairGeometry(const cReferenceSequences& ref_seq_info, bam_alignment& a1, bam_alignment& a2)
  {
    // First, sort matches by their order in the query
    q1 = &a1;
    q2 = &a2;
    q1->query_stranded_bounds_1(q1_start, q1_end);
    q2->query_stranded_bounds_1(q2_start, q2_end);

    // Reverse the coordinates to be consistently such that 1 refers to lowest...
    if (q2_start < q1_start)
    {
      swap(q1_start, q2_start);
      swap(q1_end, q2_end);
      swap(q1, q2);
    }

    // create hash key and store information about the location of this hit
    hash_strand_1 = q1->reversed();
    hash_seq_id_1 = ref_seq_info[q1->reference_target_id()].m_seq_id;
    const cAnnotatedSequence * ref_seq_1(&ref_seq_info[q1->reference_target_id()]);

    hash_strand_2 = !q2->reversed();
    hash_seq_id_2 = ref_seq_info[q2->reference_target_id()].m_seq_id;
    const cAnnotatedSequence * ref_seq_2(&ref_seq_info[q2->reference_target_id()]);

    // how much overlap is there between the two matches?
    // positive if the middle sequence can match either side of the read
    // negative if there is sequence in the read NOT matched on either side
    overlap = -1 * (q2_start - q1_end - 1);

    //
    // OVERLAP MISMATCH CORRECTION
    //
    // If there are mismatches in one or the other read in the overlap region
    // then we need to adjust the coordinates. Why? All sequences that we
    // retrieve are from the reference sequence and there are two choices
    // for where to extract this non-necessarily identical sequence!

    // save these as variables, because we may have to adjust them
    r1_start = q1->reference_start_1();
    r1_end = q1->reference_end_1();
    r2_start = q2->reference_start_1();
    r2_end = q2->reference_end_1();

    // Adjust the overlap in cases where there is a mismatch within the overlap region
    int32_t overlap_in_reference = 0;

    if (q1->strand() == q2->strand()) {

      if ((r2_start >=  r1_start) && (r2_start <=  r1_end)) {
        if (q1->strand() == +1)
          overlap_in_reference = r1_end - r2_start + 1;
        else
          overlap_in_reference = r2_end - r1_start + 1;
      }
    }

    if ((overlap > 0) || (overlap_in_reference > 0))
    {
      overlap = max(overlap, overlap_in_reference);

      int32_t q1_move, q2_move, r1_move, r2_move;
      q1->num_matches_from_end(ref_seq_info, false, overlap, q1_move, r1_move);
      q2->num_matches_from_end(ref_seq_info, true, overlap, q2_move, r2_move);

      if (q1_move >= 0)
      {
        // change where it ENDS
        q1_end -= q1_move;
        if (q1->reversed())
          r1_start += r1_move;
        else
          r1_end -= r1_move;
      }

      if (q2_move >= 0)
      {
        // change where it STARTS
        q2_start += q2_move;
        if (!q2->reversed())
          r2_start += r2_move;
        else
          r2_end -= r2_move;
      }

      // JEB 2013-10-12
      // We might re-check that they still have the required amount of unique length if they moved.
      // For now this is checked only on the original alignments that have not had their overlap corrected.

      //re-calculate the overlap
      overlap = -1 * (q2_start - q1_end - 1);
    }

    //Recalculate the overlap in the reference
    if (q1->strand() == q2->strand()) {

      if ((r2_start >=  r1_start) && (r2_start <=  r1_end)) {
        if (q1->strand() == +1)
          overlap_in_reference = r1_end - r2_start + 1;
        else
          overlap_in_reference = r2_end - r1_start + 1;
      }
    }

    // create hash coords AFTER overlap adjustment
    hash_coord_1 = (hash_strand_1) ? r1_start : r1_end;
    hash_coord_2 = (hash_strand_2) ? r2_start : r2_end;

    // Further correction for zero overlap.
    //
    // If there are multiple ways the two sides could have been aligned...shift them
    // over so as much is included in the lower reference coordinate side as possible
    //
    // This case can arise for reads matching the same reference bases in SSAHA2
    // or after correcting for mismatches in the overlap (?)

    if (overlap == 0)
    {
      int32_t lower_coord_side = (hash_coord_1 < hash_coord_2) ? -1 : +1;
      int32_t move_r1_pos = (hash_strand_1) ? -1 : +1;
      int32_t move_r2_pos = (!hash_strand_2) ? -1 : +1;
      move_r1_pos *= lower_coord_side;
      move_r2_pos *= lower_coord_side;

      uint32_t test_r1_pos = (hash_strand_1) ? r1_start : r1_end;
      uint32_t test_r2_pos = (hash_strand_2) ? r2_start : r2_end;

      if (lower_coord_side == -1)
        test_r2_pos += move_r2_pos;
      else
        test_r1_pos += move_r1_pos;

      if ( (test_r1_pos >= 1) && (test_r1_pos <= ref_seq_1->get_sequence_length() )
        && (test_r2_pos >= 1) && (test_r2_pos <= ref_seq_2->get_sequence_length() ) )
      {
        string test_r1_char;
        string test_r2_char;

        test_r1_char = ref_seq_1->get_sequence_1(test_r1_pos);
        if (hash_strand_1) test_r1_char = reverse_complement(test_r1_char);
        test_r2_char = ref_seq_2->get_sequence_1(test_r2_pos);
        if (!hash_strand_2) test_r2_char = reverse_complement(test_r2_char);

        while (test_r1_char == test_r2_char)
        {
          test_r1_pos += move_r1_pos;
          test_r2_pos += move_r2_pos;

          if (! (
                 (test_r1_pos >= 1) && (test_r1_pos <= ref_seq_1->get_sequence_length())
                 && (test_r2_pos >= 1) && (test_r2_pos <= ref_seq_2->get_sequence_length())
                 ) )
          {
            test_r1_pos -= move_r1_pos;
            test_r2_pos -= move_r2_pos;
            break;
          }

          test_r1_char = ref_seq_1->get_sequence_1(test_r1_pos);
          if (hash_strand_1) test_r1_char = reverse_complement(test_r1_char);
          test_r2_char = ref_seq_2->get_sequence_1(test_r2_pos);
          if (!hash_strand_2) test_r2_char = reverse_complement(test_r2_char);
        }

        // backtrack by one
        if (lower_coord_side == -1)
          test_r2_pos -= move_r2_pos;
        else
          test_r1_pos -= move_r1_pos;

        hash_coord_1 = test_r1_pos;
        hash_coord_2 = test_r2_pos;
      }
    }

    // these are the positions of the beginning and end of the read, across the junction
    // query 1 is the start of the read, which is why we hash by this coordinate
    // (it is less likely to be shifted by a nucleotide or two by base errors)
    read_begin_coord = (hash_strand_1) ? r1_end : r1_start;
  }

	bool CandidateJunctions::alignment_pair_to_candidate_junction(
                                                            const Settings& settings,
                                                            Summary& summary,
                                                            const cReferenceSequences& ref_seq_info,
                                                            AlignmentPair& ap,
                                                            JunctionCandidatePtr& returned_junction_candidate,
                                                            const vector<MateAln>* mate_positions,
                                                            const string& majority_orientation,
                                                            double distance_cutoff
                                                            )
	{
    (void) settings;
    bool verbose = false;

    //if (ap.a1.read_name() == "1:2096442") {
    //   verbose = true;
    //}

    // clear the return value
    returned_junction_candidate = JunctionCandidatePtr(NULL);
    

		// set up local settings
		int32_t flanking_length = summary.sequence_conversion.read_length_max;

    bam_alignment& a1 = ap.a1;
    bam_alignment& a2 = ap.a2;

		// Method
		//
		// Hash junctions by a key showing the inner coordinate of the read.
		// and the direction propagating from that position in the reference sequence.
		// Prefer the lower coordinate side of the read for main hash.
		//
		// REL606__1__1__REL606__4629812__-1__0
		// means the junction sequence is 36-1 + 4629812-4629777 from the reference sequence
		//
		// On the LEFT side: -1 means this is highest coord of alignment, junction seq begins at lower coord
		//                    1 means this is lowest coord of alignment, junction seq begins at higher coord
		// On the RIGHT side:-1 means this is highest coord of alignment, junction seq continues to lower coord
		//                    1 means this is lowest coord of alignment, junction seq continues to higher coord

    // Reference-side geometry (breakpoint coordinates, corrected for mismatches near the
    // junction) is computed by AlignmentPairGeometry. Copy into local mutable variables,
    // since the rest of this function adjusts them further (overlap normalization, strand
    // canonicalization).
    AlignmentPairGeometry geom(ref_seq_info, a1, a2);

    bool hash_strand_1 = geom.hash_strand_1;
    bool hash_strand_2 = geom.hash_strand_2;
    string hash_seq_id_1 = geom.hash_seq_id_1;
    string hash_seq_id_2 = geom.hash_seq_id_2;
    int32_t hash_coord_1 = geom.hash_coord_1;
    int32_t hash_coord_2 = geom.hash_coord_2;
    int32_t overlap = geom.overlap;
    int32_t r1_start = geom.r1_start, r1_end = geom.r1_end;
    int32_t r2_start = geom.r2_start, r2_end = geom.r2_end;
    int32_t read_begin_coord = geom.read_begin_coord;
    bam_alignment& q1 = *geom.q1;
    uint32_t q1_end = geom.q1_end;

    // --paired-mapping: is each side's flanking piece concordant with the (not-split) mate?
    // pin_c1/pin_c2 track side-1 (geom.q1) / side-2 (geom.q2) BEFORE the strand-canonicalization
    // swap below; they are swapped in lockstep with the coords so they stay attached to the
    // right final side. The read-start piece is geometrically across the junction from the mate
    // and won't be concordant; the read-end piece pins its side to the copy the mate sits next to.
    bool pin_c1 = false, pin_c2 = false;
    if (mate_positions && !majority_orientation.empty() && (distance_cutoff > 0.0)) {
      pin_c1 = mate_concordant_with_piece(*geom.q1, ref_seq_info, *mate_positions, majority_orientation, distance_cutoff);
      pin_c2 = mate_concordant_with_piece(*geom.q2, ref_seq_info, *mate_positions, majority_orientation, distance_cutoff);
    }

		// Calculate an offset that only applies if the overlap is positive (sequence is shared between the two ends)
		int32_t overlap_offset = (overlap > 0) ? overlap : 0;
		if (verbose)
			cout << "Overlap offset: " << overlap_offset << endl;
     
    ////
		// Create the sequence of the candidate junction
    ////
		    
    // Get the unique read sequence... (coords need to be on correct strand)
    
    string unique_read_seq_string = "";
    if (overlap < 0) {
			unique_read_seq_string = q1.read_char_stranded_sequence_1(q1_end + 1, q1_end - overlap);
      if (verbose) {
        cout << "Unique read sequence: " << unique_read_seq_string << endl;
      }
    }
    
    // Construct a sequence from a temporary DiffEntry
    cDiffEntry jc(JC);
    jc[SIDE_1_SEQ_ID] = to_string(hash_seq_id_1);
    jc[SIDE_2_SEQ_ID] = to_string(hash_seq_id_2);
    jc[SIDE_1_POSITION] = to_string(hash_coord_1);
    jc[SIDE_2_POSITION] = to_string(hash_coord_2);
    jc[SIDE_1_STRAND] = hash_strand_1 ? "1" : "-1";
    jc[SIDE_2_STRAND] = hash_strand_2 ? "1" : "-1";
    jc["overlap"] = to_string(overlap);
    jc["unique_read_sequence"] = unique_read_seq_string;
    
    // Why do we still need to normalize overlap? It turns out having an N in the read in the overlap region can cause
    // too much overlap to be removed and a user junction comes out having a diff sequence because of this. So, this
    // call will add back that overlap by looking only at the reference
    if (verbose) cout << "Normalizing overlap:\n" << jc.as_string() << endl;
    normalize_junction_overlap(ref_seq_info,jc);
    
    // If we changed, we need to back everything out to continue the code normally
    if (overlap != from_string<int32_t>(jc["overlap"])) {
      // Seq_ids won't change...
      hash_coord_1 = from_string<int32_t>(jc[SIDE_1_POSITION]);
      hash_coord_2 = from_string<int32_t>(jc[SIDE_2_POSITION]);
      hash_strand_1 = (jc[SIDE_1_STRAND] == "1");
      hash_strand_2 = (jc[SIDE_2_STRAND] == "1");
      overlap = from_string<int32_t>(jc["overlap"]);
      unique_read_seq_string = jc["unique_read_sequence"];
    }
    if (verbose) cout << "Converted to:\n" << jc.as_string() << endl;

    string junction_seq_string = construct_junction_sequence(ref_seq_info, jc, flanking_length);
    
    // Save these return values - we do not allow coords to be adjusted within construct_junction_sequence
    int32_t flanking_left = from_string<int32_t>(jc["flanking_left"]);
    int32_t flanking_right = from_string<int32_t>(jc["flanking_right"]);
    
		// create hash coords after adjustment for overlap
		if (!hash_strand_1)
			r1_end -= overlap_offset;
		else //reversed
			r1_start += overlap_offset;

		if (!hash_strand_2)
			r2_end -= overlap_offset;
		else //reversed
			r2_start += overlap_offset;

    ///
    //  Important: Here is where we choose the strand of the junction sequence
    ///
    
    // At this point, we don't know whether a certain side is redundant.
    //
		// Rules...
    // If they are both on the same
		if ( hash_seq_id_1.compare(hash_seq_id_2) > 0 || ((hash_seq_id_1.compare(hash_seq_id_2) == 0) && (hash_coord_2 < hash_coord_1)) )
		{
			swap(hash_coord_1, hash_coord_2);
			swap(hash_strand_1, hash_strand_2);
			swap(hash_seq_id_1, hash_seq_id_2);
			swap(flanking_left, flanking_right);
			swap(pin_c1, pin_c2); // keep pair-concordance attached to its (now swapped) side

			junction_seq_string = reverse_complement(junction_seq_string);
			unique_read_seq_string = reverse_complement(unique_read_seq_string);
		}
    
    // create junction candidate pointer
		JunctionCandidate* candidate_junction_ptr = new JunctionCandidate(
          JunctionInfo(
                       JunctionSide(hash_seq_id_1,	hash_coord_1,	hash_strand_1 ? +1 : -1), // note conversion of strand 0/1 to -1/+1
                       JunctionSide(hash_seq_id_2,	hash_coord_2,	hash_strand_2 ? +1 : -1), // note conversion of strand 0/1 to -1/+1
                       overlap, 				
                       unique_read_seq_string,
                       flanking_left,
                       flanking_right
                       ), 
          junction_seq_string
    );
    candidate_junction_ptr->read_begin_hash[read_begin_coord]++;

    // Provisionally pin the side whose flanking piece is concordant with the mate. Only when
    // EXACTLY ONE side is concordant (the other piece is across the junction on a different
    // locus); if both or neither pass, this read gives no unambiguous placement. The "exactly one
    // candidate per read" filter is applied by the caller (alignments_to_candidate_junctions).
    if (pin_c1 != pin_c2) {
      int32_t s = pin_c1 ? 0 : 1;
      candidate_junction_ptr->sides[s].pair_unique_count = 1;   // this read = one concordant pair
      candidate_junction_ptr->sides[s].pair_confirmed_position = candidate_junction_ptr->sides[s].position;
    }

    // set the return value (which takes control of the allocated pointer)
    returned_junction_candidate = JunctionCandidatePtr(candidate_junction_ptr);

		if (verbose)
		{
      string junction_id = candidate_junction_ptr->junction_key();
			cout << "READ ID: " << a1.read_name() << endl;
			cout << "JUNCTION ID: " << junction_id << endl;
		}

		ASSERT(junction_seq_string.size() > 0, "Junction sequence not found."); 
		ASSERT(junction_seq_string.size() == flanking_left + flanking_right + static_cast<uint32_t>(abs(overlap)),
           "Incorrect junction sequence length for " + candidate_junction_ptr->junction_key() + "\n" + junction_seq_string);
    
		return true;
	}

	uint64_t CandidateJunctions::alignments_to_candidate_junctions(
                                                             const Settings& settings,
                                                             Summary& summary, const
                                                             cReferenceSequences& ref_seq_info,
                                                             SequenceToKeyToJunctionCandidateMap& candidate_junctions,
                                                             alignment_list& alignments,
                                                             const vector<MateAln>* mate_positions,
                                                             const string& majority_orientation,
                                                             double distance_cutoff
                                                             )
	{
    (void)summary;
		bool verbose = false;

		if (verbose)
		{
			cout << endl << "###########################" << endl;
			cout << alignments.front()->read_name();
			cout << endl << "###########################" << endl;
		}

		// Must  have multiple matches to support a new junction.
		if (alignments.size() <= 1)
			return 0;

    ////
    // Split the read alignments into two lists
    //  (1) list1 contains all matches starting at the beginning of the read
    //  (2) list2 contains all other matches
    // Accepted pairs must have a member from each list
    // Requiring (1) saves a number of comparisons and gets rid of a lot of bad matches.
    ////

    alignment_list list1, list2;
    
		if (verbose)
		{
			cout << alignments.front()->read_name() << endl;
			cout << "Total matches: " << alignments.size() << endl;
		}

    const uint32_t unmatched_end_min_coord = settings.required_junction_read_end_min_coordinate(alignments.front()->read_length());
    
    for (alignment_list::iterator it=alignments.begin(); it != alignments.end(); it++)
		{
			bam_alignment_ptr a = *it;
      
			uint32_t a_start, a_end;
      a->query_stranded_bounds_1(a_start, a_end);

			if (verbose) cout << "(" << a_start << ", " << a_end << ")" << endl;
          
      // In order to pass later guards have to match at least a minimum amount of the read
      uint32_t min_match_length = static_cast<uint32_t>(ceil(a->read_length() * settings.required_both_unique_length_per_side_fraction));
      if (a_end - a_start + 1 < min_match_length)
        continue;
      
			if (a_start == 1) {
				list1.push_back(a);
        if (verbose) cout << "  List 1" << endl;
      }
			else if (a_end >= unmatched_end_min_coord) {
				list2.push_back(a);
        if (verbose) cout << "  List 2" << endl;
      }
    }
    
		// The first match in this category is the longest
		if (verbose)
		{
			cout << "  List1: " << list1.size() << endl;
			cout << "  List2: " << list2.size() << endl;
		}


    ////
		// Try adding together each pair of matches to make a junction, by looking at read coordinates
    //   Only keep pairs that cover the maximum number of query bases encountered
    ////
    
    int32_t max_end_to_end_length = 0;
		vector<AlignmentPair> passed_pair_list;

    for (alignment_list::iterator it1 = list1.begin(); it1 != list1.end(); it1++)
		{
			bam_alignment_ptr& a1 = *it1;

      for (alignment_list::iterator it2 = list2.begin(); it2 != list2.end(); it2++)
			{
				bam_alignment_ptr& a2 = *it2;

				// constructing the alignment pair calculates statistics about their overlap
        // and tests the guards that are in Settings.
        
        AlignmentPair ap(*a1, *a2, settings);

				if (ap.pass)
        {
          // Use the end_to_end_length as a kind of score
          // matches with greater end_to_end_length trump those with smaller
          if (ap.end_to_end_length > max_end_to_end_length)
          {
            max_end_to_end_length = ap.end_to_end_length;
            passed_pair_list.clear();
          }

          passed_pair_list.push_back(ap);
        }
      }
		}
    
    // Ignore matches that predict many highly redundant junctions!!
    if (passed_pair_list.size() > settings.highly_redundant_junction_ignore_passed_pair_limit)
      return 0;
    
    // Build all of this read's candidates first, so the --paired-mapping "exactly one
    // pair-concordant candidate disambiguates" rule can be applied across the whole set before
    // any are stored (each passing pair may pin a side via mate concordance).
    vector<JunctionCandidatePtr> this_read_candidates;
		for (uint32_t i = 0; i < passed_pair_list.size(); i++)
		{
			AlignmentPair& ap = passed_pair_list[i];

      JunctionCandidatePtr new_junction_ptr;
			bool passed = alignment_pair_to_candidate_junction(settings, summary, ref_seq_info, ap, new_junction_ptr,
                                                         mate_positions, majority_orientation, distance_cutoff);
			if (!passed) continue;
      this_read_candidates.push_back(new_junction_ptr);
    }

    // Exactly-one rule: the mate disambiguates only if it pins a SINGLE reference coordinate.
    // When both junction sides are repeats, one spanning read produces several candidates that all
    // pin the same side to the same coordinate (they differ only in the other, still-ambiguous
    // side) -- that is ONE distinct pin and must promote. But if the mate is concordant with two
    // DIFFERENT copies (different coordinates), the placement is genuinely ambiguous and no pin
    // survives. So key the decision on the set of distinct pinned positions, not the candidate count.
    if (mate_positions) {
      set<int32_t> pinned_positions;
      for (uint32_t k = 0; k < this_read_candidates.size(); k++)
        for (int32_t s = 0; s < 2; s++)
          if (this_read_candidates[k]->sides[s].pair_unique_count > 0)
            pinned_positions.insert(this_read_candidates[k]->sides[s].pair_confirmed_position);
      if (pinned_positions.size() != 1) {
        for (uint32_t k = 0; k < this_read_candidates.size(); k++) {
          for (int32_t s = 0; s < 2; s++) {
            this_read_candidates[k]->sides[s].pair_unique_count = 0;
            this_read_candidates[k]->sides[s].pair_confirmed_position = 0;
          }
        }
      }
    }

    // Store / update the candidates.
    for (uint32_t k = 0; k < this_read_candidates.size(); k++)
		{
      JunctionCandidatePtr new_junction_ptr = this_read_candidates[k];
      JunctionCandidate& new_junction = *new_junction_ptr;
      if (verbose) cout << "Testing junction: " << new_junction_ptr->junction_key() << endl << new_junction_ptr->sequence << endl;

			string junction_id = new_junction.junction_key();
      if (verbose) cout << junction_id << endl;

      if ((candidate_junctions.count(new_junction.sequence) == 0) || (candidate_junctions[new_junction.sequence].count(junction_id) == 0)) {
        // this is a new combination of sequence and id
        if (verbose) cout << "New saved junction " << new_junction.sequence << " " << junction_id << endl;
        candidate_junctions[new_junction.sequence][junction_id] = new_junction_ptr;
      }
      else
      {
        // update score of existing junction
        JunctionCandidate& cj = *candidate_junctions[new_junction.sequence][junction_id];
        cj.read_begin_hash[new_junction.read_begin_hash.begin()->first]++;
        // Same key => identical coordinates, so accumulate this read's concordant-pair count onto
        // the stored candidate (each supporting read contributes +1). A later threshold pass in
        // identify_candidate_junctions converts the tally into the pair_unique flag.
        for (int32_t s = 0; s < 2; s++) {
          if (new_junction.sides[s].pair_unique_count > 0) {
            cj.sides[s].pair_unique_count += new_junction.sides[s].pair_unique_count;
            cj.sides[s].pair_confirmed_position = new_junction.sides[s].pair_confirmed_position;
          }
        }
        if (verbose) cout << "Updating score of existing " << new_junction.sequence << " " << junction_id << endl
          << "Pos: " << new_junction.read_begin_hash.begin()->first << " Score: " << cj.pos_hash_score()  << endl;
      }

    } // end passed pair list

    return passed_pair_list.size();
	}

  
  // Extends and moves overlap so that it is uniformly on one side.
  // Use before getting sequence on user-input junctions and junctions
  // predicted by other programs such as TopHat.
  
  void CandidateJunctions::normalize_junction_overlap (
         const cReferenceSequences& ref_seq_info,
         cDiffEntry& jc
         )
  {
    //ASSERT(false, "CandidateJunctions::normalize_junction_overlap not fully tested");
    bool verbose = false;

    // We need to back out any alignment_overlap that was resolved to change
    //   to overlap=0 at the end of a previous breseq run
    int32_t alignment_overlap = jc.entry_exists("alignment_overlap") ? from_string<int32_t>(jc["alignment_overlap"]) : 0;
    int32_t overlap = from_string<int32_t>(jc["overlap"]);

    // alignment_overlap can be negative, in which case both sides are zero and some checks fail
    if (alignment_overlap != overlap) {
    
      int32_t side_1_overlap = jc.entry_exists("side_1_overlap") ? from_string<int32_t>(jc["side_1_overlap"]) : 0;
      int32_t side_2_overlap = jc.entry_exists("side_2_overlap") ? from_string<int32_t>(jc["side_2_overlap"]) : 0;
      ASSERT(alignment_overlap == side_1_overlap + side_2_overlap, "overlap != side_1_overlap + side_2_overlap" + jc.as_string() );
      
      // This number means this side was shifted this much in assigning overlap to each side
      //   if it was all assigned to this side, then we don't need to move it
      int32_t shift_side_2 = alignment_overlap - side_2_overlap;
      if (shift_side_2 > 0) {
        jc["side_2_position"] = to_string<int32_t>( from_string<int32_t>(jc["side_2_position"]) - from_string<int32_t>(jc["side_2_strand"]) * shift_side_2  );
      }
      
      int32_t shift_side_1 = alignment_overlap - side_1_overlap;
      if (shift_side_1 > 0) {
        jc["side_1_position"] = to_string<int32_t>( from_string<int32_t>(jc["side_1_position"]) - from_string<int32_t>(jc["side_1_strand"]) * shift_side_1  );
      }
    }
    
    if (alignment_overlap > 0) {
      overlap = alignment_overlap;
    }
    
    // This section catches cases where zero overlap was assigned because it
    // was converted from a prediction by a different program (e.g. TopHat2)
    // it will also expand overlap if it was too small to begin with...
    
    // Do nothing if there is negative overlap (unique junction sequence).
    if (overlap < 0)
      return;
    
    int32_t hash_strand_1 = from_string<int32_t>(jc["side_1_strand"]);
    int32_t hash_strand_2 = from_string<int32_t>(jc["side_2_strand"]);
    
    int32_t hash_coord_1 = from_string<int32_t>(jc["side_1_position"]);
    int32_t hash_coord_2 = from_string<int32_t>(jc["side_2_position"]);
    
    int32_t overlap_offset = max(0, overlap);
    
    const cAnnotatedSequence& ref_seq_1 = ref_seq_info[jc["side_1_seq_id"]];
    const cAnnotatedSequence& ref_seq_2 = ref_seq_info[jc["side_2_seq_id"]];
    

    // overlap may have been resolved and we need to shift coords and correct for that...
    
    // First shift things as far to the left as possible,
    // then shift to the right to count up the overlap
    
    // reverse direction
    int32_t reverse_overlap = 0;
    {
      
      int32_t test_pos_1 = hash_coord_1 + overlap_offset * hash_strand_1;
      int32_t test_pos_2 = hash_coord_2 - hash_strand_2;
      
      //@JEB: notice the minus sign added to the first strand in get_stranded_sequence_1
      //      due to how junction directions are labeled
      while ( (test_pos_1 >= 1) && (test_pos_2 >= 1)
              && (static_cast<uint32_t>(test_pos_1) <= ref_seq_1.get_sequence_length())
              && (static_cast<uint32_t>(test_pos_2) <= ref_seq_2.get_sequence_length())
              && (ref_seq_1.get_stranded_sequence_1(-hash_strand_1, test_pos_1) == ref_seq_2.get_stranded_sequence_1(hash_strand_2, test_pos_2)) ) {
        test_pos_1 += hash_strand_1;
        test_pos_2 -= hash_strand_2;
        reverse_overlap++;
        if (test_pos_1 < 1) break;
        if (test_pos_2 < 1) break;
        if (static_cast<uint32_t>(test_pos_1) > ref_seq_1.get_sequence_length()) break;
        if (static_cast<uint32_t>(test_pos_2) > ref_seq_2.get_sequence_length()) break;
      }
    }
    
    // forward direction
    int32_t forward_overlap = 0;
    {
      int32_t test_pos_1 = hash_coord_1 - hash_strand_1;
      int32_t test_pos_2 = hash_coord_2 + overlap_offset * hash_strand_2;
      // most of these check for remaining in bounds
      while ( (test_pos_1 >= 1) && (test_pos_2 >= 1)
              && (static_cast<uint32_t>(test_pos_1) <= ref_seq_1.get_sequence_length())
              && (static_cast<uint32_t>(test_pos_2) <= ref_seq_2.get_sequence_length())
              && (ref_seq_1.get_stranded_sequence_1(-hash_strand_1, test_pos_1) == ref_seq_2.get_stranded_sequence_1(hash_strand_2, test_pos_2)) ) {
        test_pos_1 -= hash_strand_1;
        test_pos_2 += hash_strand_2;
        forward_overlap++;
      }
    }
    
    hash_coord_2 -= hash_strand_2 * reverse_overlap;
    hash_coord_1 -= hash_strand_1 * forward_overlap;
    
    overlap_offset = overlap_offset + forward_overlap + reverse_overlap;
    if (verbose) {
      cout << "Adjusted for overlap:" << endl;
      cout << " Hash coord 1:" << hash_coord_1 << endl;
      cout << " Hash coord 2:" << hash_coord_2 << endl;
      cout << " Reverse overlap:" << reverse_overlap << endl;
      cout << " Forward overlap:" << forward_overlap << endl;
      cout << " Overlap offset:" << overlap_offset << endl;
    }
    
    // save new values
    jc["side_1_position"] = to_string<int32_t>(hash_coord_1);
    jc["side_2_position"] = to_string<int32_t>(hash_coord_2); 
    jc["overlap"] = to_string<int32_t>(overlap_offset);     
  }
  
  //
  // string CandidateJunctions::construct_junction_sequence
  //
	// This is a Swiss-Army Knife function for reconstructing a junction sequence
  // from the parameters that fully describe it.
  //
  // It is used:
  //   1) In the main pipeline to output candidate junctions for re-alignment (exclusive mode).
  //   2) To define specific junctions to look for in a sample (exclusive mode).
  //   3) When comparing junction predictions from different tools (inclusive mode). 
  //
  // In the input cDiffEntry of 'JC' type, these fields must be provided:
  //  'side_1_seq_id', 'side_1_position', 'side_1_strand',
  //  'side_2_seq_id', 'side_2_position', 'side_2_strand',
  //  'overlap' must be defined, but can be corrected from zero to the true value
  //  'unique_junction_sequence' may be optionally provided...
  //  If it is provided, then overlap will be set to be negative its length.
  //
  // The flanking sequence length can be applied two different ways:
  //   Such that it includes overlap (needed for comparing junctions between programs)
  //   Such that it excludes overlap (needed for junction prediction within breseq)
  //
  // Exclusive means that if there are 3 bases of overlap and the flanking length is 50, 
  // (with no unique junction sequence) then the output sequence has a length of 103.
  //
  // Inclusive means that if there are 3 bases of overlap and the flanking length is 50, 
  // (with no unique junction sequence) then the output sequence has a length of 97.
  //
  // The diff entry is returned with a corrected flanking length if we ran into the edge of a sequence
  
  string CandidateJunctions::construct_junction_sequence( 
    const cReferenceSequences& ref_seq_info,
    cDiffEntry& jc,
    int32_t flanking_length,
    bool inclusive_overlap
    )
  {
    bool verbose = false;
    // set up local settings
    
    if (verbose) cout << "Original Junction Diff Entry:\n" << jc << endl;
        
    int32_t hash_strand_1 = from_string<int32_t>(jc["side_1_strand"]);
    int32_t hash_strand_2 = from_string<int32_t>(jc["side_2_strand"]);
    
    int32_t hash_coord_1 = from_string<int32_t>(jc["side_1_position"]);
    int32_t hash_coord_2 = from_string<int32_t>(jc["side_2_position"]);
        
    int32_t overlap = from_string<int32_t>(jc["overlap"]);
    int32_t overlap_offset = max(0, overlap);
    
    ASSERT(ref_seq_info.seq_id_exists(jc["side_1_seq_id"]), "Reference seq ID not found:" + jc["side_1_seq_id"]);
    ASSERT(ref_seq_info.seq_id_exists(jc["side_2_seq_id"]), "Reference seq ID not found:" + jc["side_2_seq_id"]);

    const cAnnotatedSequence& ref_seq_1 = ref_seq_info[jc["side_1_seq_id"]];
    const cAnnotatedSequence& ref_seq_2 = ref_seq_info[jc["side_2_seq_id"]];
    
    ////
    // Create the sequence of the candidate junction
    ////

    string junction_seq_string = "";
    
    if (inclusive_overlap) flanking_length -= overlap_offset;

    // first end - includes or excludes >0 overlap
    int32_t flanking_left = flanking_length;
    
    if (hash_strand_1 == -1) { // alignment is not reversed
      
      // start_pos is in 1-based coordinates
      int32_t start_pos = hash_coord_1 - (flanking_left - 1) - overlap_offset;
      if (start_pos < 1) {
        if (verbose)
          cout << "START POS 1: " << start_pos << " < 0" << endl;
        flanking_left += start_pos - 1;
        flanking_left = max(0, flanking_left);
        start_pos = 1;
      }
      //cout << "number 1:" << endl;
      //cout << start_pos + 1 << " " << start_pos + flanking_left + overlap_offset << endl;
      
      if (flanking_left + overlap_offset > 0) {
        string add_seq = ref_seq_1.get_sequence_1(start_pos, start_pos + flanking_left + overlap_offset - 1);        
        if (verbose) cout << "1F: " << add_seq << endl;
        junction_seq_string += add_seq;
      }
   
    } else { // alignment is reversed
      
      // end_pos is in 1-based coordinates
      int32_t end_pos = hash_coord_1 + (flanking_left - 1) + overlap_offset;
      if (end_pos > ref_seq_1.m_length) {
        if (verbose) cout << "END POS 1: " << end_pos << " > length" << endl;
        flanking_left -= end_pos - ref_seq_1.m_length;
        flanking_left = max(0, flanking_left);
        end_pos = ref_seq_1.m_length;
      }
      
      //cout << "number 1:" << endl;
      //cout << end_pos - (flanking_left + overlap_offset) + 1 << " " << end_pos << endl;
      if (flanking_left + overlap_offset > 0) {
        string add_seq = ref_seq_1.get_sequence_1(end_pos - (flanking_left + overlap_offset) + 1, end_pos);
        add_seq = reverse_complement(add_seq);
        if (verbose) cout << "1R: " << add_seq << endl;
        junction_seq_string += add_seq;
      }
    }
    
    // Add any unique junction sequence that was only in the read
    // and NOT present in the reference genome (overlap < 0)
    if (jc.count("unique_read_sequence") > 0)
      junction_seq_string += jc["unique_read_sequence"];
    
    if (verbose) cout << "Unique junction sequence: " << jc["unique_read_sequence"] << endl;
    
    // second end - added without overlapping sequence
    int32_t flanking_right = flanking_length;
    
    if (hash_strand_2 == 1) //alignment is not reversed
    {
      // end_pos is in 1-based coordinates
      int32_t end_pos = hash_coord_2 + (flanking_right - 1) + overlap_offset;
      if (end_pos > ref_seq_2.m_length) {
        if (verbose)
          cout << "END POS 2: " << end_pos << " > length" << endl;
        flanking_right -= (end_pos - ref_seq_2.m_length);
        flanking_right = max(0, flanking_right);
        end_pos = ref_seq_2.m_length;
      }
      //string add_seq = ref_seq_2.substr(end_pos - flanking_right, flanking_right);
      //cout << "number 2:" << endl;
      //cout << end_pos - flanking_right << " " << end_pos - 1 << endl;
      if (flanking_right > 0) {
        string add_seq = ref_seq_2.get_sequence_1(end_pos - flanking_right + 1, end_pos);
        if (verbose) cout << "2F: " << add_seq << endl;
        junction_seq_string += add_seq;
      }
    }
    else // alignment is reversed
    {
      // start_pos is in 1-based coordinates
      int32_t start_pos = hash_coord_2 - (flanking_right - 1) - overlap_offset;
      if (start_pos < 1)
      {
        if (verbose) cout << "START POS 2: " << start_pos << " < 0" << endl;
        flanking_right += start_pos - 1;
        flanking_right = max(0, flanking_right);
        start_pos = 1;
      }
      
      if (flanking_right > 0) {
        string add_seq = ref_seq_2.get_sequence_1(start_pos, start_pos + flanking_right - 1);
        add_seq = reverse_complement(add_seq);
        if (verbose) cout << "2R: " << add_seq << endl;
        junction_seq_string += add_seq;
      }
    }
  
    // Check the length - debug code
    /* @JEB This check assumes that we did not subtract off flanking earlier!!
    if (inclusive_overlap && ( static_cast<int32_t>(junction_seq_string.size()) != flanking_left + flanking_right - abs(overlap))) {
      stringstream s;
      s << jc << endl;
      ERROR( s.str() + "Incorrect junction sequence length: "  +  to_string(junction_seq_string.size()));
    }
    else 
    */
    if (static_cast<int32_t>(junction_seq_string.size()) != flanking_left + flanking_right + abs(overlap)) {
      stringstream s;
      s << jc << endl;
      ERROR( s.str() + "Incorrect junction sequence length: "  +  to_string(junction_seq_string.size()) + "\nExpected junction sequence length: " + to_string(flanking_left + flanking_right + abs(overlap)));
    }
    
    //cout << hash_coord_1 << " " << hash_strand_1 << " " << hash_coord_2 << " " << hash_strand_2 << endl;
    //cout << junction_seq_string << endl;
        
    jc["flanking_left"] = to_string<int32_t>(flanking_left);
    jc["flanking_right"] = to_string<int32_t>(flanking_right);
    
    if (verbose) cout << "Returned Junction Diff Entry:\n" << jc << endl;
    
    return junction_seq_string;
  }
  
  
  /*!
   * Calculates statistics and tests whether it passes required conditions in settings
   */
  
  AlignmentPair::AlignmentPair(bam_alignment& _a1, bam_alignment& _a2, const Settings &settings)
  : a1(_a1), a2(_a2), hash_coord(0)
    , a1_unique_start(0), a1_unique_end(0), a1_unique_length(0)
    , a2_unique_start(0), a2_unique_end(0), a2_unique_length(0)
    , end_to_end_length(0), union_length(0), intersection_length(0)
    , pass(false)
  {
    calculate_union_and_unique();
    pass = test(settings);
  }
  
  /*!
   * Calculates statistics about how the alignments overlap
   */
  
  void AlignmentPair::calculate_union_and_unique()
  {
    bool verbose = false;
    
    uint32_t a1_start, a1_end;
    a1.query_stranded_bounds_1(a1_start, a1_end);
    
    uint32_t a2_start, a2_end;
    a2.query_stranded_bounds_1(a2_start, a2_end);
    
    // move adjust so order is a1-a2 on read sequence
    if (a1_start > a2_start)
    {
      swap(a1, a2);
      swap(a1_start, a2_start);
      swap(a1_end, a2_end);
    }

    // this is the coord of the first base in the read
    // (since we required this to match)
    uint32_t a1_reference_start, a1_reference_end;
    a1.reference_stranded_bounds_1(a1_reference_start, a1_reference_end);
    hash_coord = a1_reference_start;
    ASSERT(a1_start == 1, "Read does not match from beginning.");
    
		int32_t intersection_start = a2_start;
		int32_t intersection_end = a1_end;
		this->intersection_length = intersection_end - intersection_start + 1;
    
    int32_t union_start = a1_start;
		int32_t union_end = a2_end;
    this->end_to_end_length = union_end - union_start + 1;
    // Note: last term subtracts missing bases when the two alignments don't overlap in the middle
		this->union_length = end_to_end_length - max(0, -intersection_length);
    
    a1_unique_start = a1_start;
    a1_unique_end = min(a2_start - 1, a1_end);
    this->a1_unique_length = max(a1_unique_end - a1_unique_start + 1, 0);
    
    a2_unique_start = max(a1_end + 1, a2_start);
    a2_unique_end = a2_end;
    this->a2_unique_length = max(a2_unique_end - a2_unique_start + 1, 0);
    
    if (verbose)
    {
      cout << " Read: " << a1.read_name() << endl;
      cout << "=== Match1: " << a1_start << "-" << a1_end << "   Match2: " << a2_start << "-" << a2_end << endl;
			cout << "    Union: " << this->union_length << "   Intersection: " << this->intersection_length << endl;
      
      cout << "    Unique length 1: " << this->a1_unique_length << " Unique length 2:"<< this->a2_unique_length << endl;
    }
  }
  
  /*!
   * Calculates whether the read pair passes required conditions in settings
   */
  
  bool AlignmentPair::test(const Settings& settings)
  {
    int32_t intersection_length_negative = -min(0, intersection_length);
    int32_t intersection_length_positive = max(0, intersection_length);
    
    int32_t  scaled_maximum_junction_sequence_negative_overlap_length_fraction = settings.maximum_junction_sequence_negative_overlap_length_minimum +
      static_cast<int32_t>(floor(static_cast<double>(a1.read_length() - settings.maximum_junction_sequence_negative_overlap_length_minimum) * settings.maximum_junction_sequence_negative_overlap_length_fraction));

    int32_t  scaled_maximum_junction_sequence_positive_overlap_length_fraction = settings.maximum_junction_sequence_positive_overlap_length_minimum +
    static_cast<int32_t>(floor(static_cast<double>(a1.read_length() - settings.maximum_junction_sequence_positive_overlap_length_minimum) * settings.maximum_junction_sequence_positive_overlap_length_fraction));

    
		//// Require negative overlap (inserted unique sequence length) to be less than some value
		if (intersection_length_negative > scaled_maximum_junction_sequence_negative_overlap_length_fraction)
			return false;
    
		if (settings.maximum_junction_sequence_insertion_length &&
        (intersection_length_negative > static_cast<int32_t>(settings.maximum_junction_sequence_insertion_length)))
			return false;
    
    //// Require positive overlap (shared by both ends) to be less than some value
    if (intersection_length_positive > scaled_maximum_junction_sequence_positive_overlap_length_fraction)
			return false;
    
    if (settings.maximum_junction_sequence_overlap_length && 
        (intersection_length_positive > static_cast<int32_t>(settings.maximum_junction_sequence_overlap_length)))
			return false;
    
		//// Require both ends to extend a certain minimum length outside of the overlap

    // This can be an absolute number or a fraction of the read length
    int32_t scaled_required_both_unique_length_per_side 
      = static_cast<int32_t>(ceil(static_cast<double>(a1.read_length()) * settings.required_both_unique_length_per_side_fraction));

		if (a1_unique_length < static_cast<int32_t>(settings.required_both_unique_length_per_side))
			return false;
        
    if (a1_unique_length < scaled_required_both_unique_length_per_side)
      return false;
    
		if (a2_unique_length < static_cast<int32_t>(settings.required_both_unique_length_per_side))
			return false;
    
    if (a2_unique_length < scaled_required_both_unique_length_per_side)
      return false;

		//// Require one end to extend a higher minimum length outside of the overlap
		if ((a1_unique_length <  static_cast<int32_t>(settings.required_one_unique_length_per_side))
        && (a2_unique_length <  static_cast<int32_t>(settings.required_one_unique_length_per_side)))
			return false;
    
		//// Test all of the normal criteria for counting a match to the reference
		if (end_to_end_length < static_cast<int32_t>(settings.require_match_length))
			return false;
    
    if (end_to_end_length < settings.require_match_fraction * static_cast<double>(a1.read_length()) )
      return false;

    return true;
  }


  
} // namespace breseq
