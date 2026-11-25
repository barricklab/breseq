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

#include "fastq.h"

#include <zlib.h>
#include "reference_sequence.h"

using namespace std;

namespace breseq {

  // How much of a FASTQ file, as stored on disk, estimate_fastq_file_bases() reads
  static const uint64_t FASTQ_ESTIMATE_COMPRESSED_BYTES = 16 * 1024 * 1024;

  uint64_t estimate_fastq_file_bases(const string &file_name)
  {
    uint64_t file_size = 0;
    {
      ifstream f(file_name.c_str(), ios::binary | ios::ate);
      ASSERT(!f.fail(), "Could not open file: " + file_name);
      file_size = static_cast<uint64_t>(f.tellg());
    }
    if (file_size == 0) return 0;

    // gzopen reads a file that is not gzipped as it is, and gzoffset gives the position in the
    // file on disk either way, which is what the scaling below needs.
    gzFile gz = gzopen(file_name.c_str(), "rb");
    ASSERT(gz != NULL, "Could not open file: " + file_name);

    vector<char> buffer(1 << 20);
    uint64_t bases = 0;
    uint64_t consumed = 0;
    uint8_t line_in_record = 0;   // 0 name, 1 sequence, 2 plus, 3 qualities
    bool reached_end = false;
    while (true) {
      int n = gzread(gz, buffer.data(), static_cast<unsigned>(buffer.size()));
      if (n <= 0) { reached_end = true; break; }
      for (int i = 0; i < n; i++) {
        const char c = buffer[i];
        if (c == '\n') line_in_record = (line_in_record + 1) % 4;
        else if ((line_in_record == 1) && (c != '\r')) bases++;
      }
      consumed = static_cast<uint64_t>(gzoffset(gz));
      if (consumed >= FASTQ_ESTIMATE_COMPRESSED_BYTES) break;
    }
    gzclose(gz);

    if (reached_end || (consumed == 0) || (consumed >= file_size)) return bases;
    return static_cast<uint64_t>(static_cast<double>(bases) * static_cast<double>(file_size) / static_cast<double>(consumed));
  }
  
  /*
   normalize_fastq
   
   correct common errors in input fastq and normalize to standard SANGER format

   */

  // Defined with the rest of the mate-pair junction adapter code, below
  namespace {
    size_t trim_mate_pair_junction_adapter(cFastqSequence& seq);
    void warn_mate_pair_junction_adapter(const string& file_name, const double fraction, const string& adapter_name, const bool trimming);
  }

  AnalyzeFastqSummary normalize_fastq(
                                        const string &file_name,
                                        const string &convert_file_name,
                                        const uint32_t file_index,
                                        const int32_t trim_end_on_base_quality,
                                        const bool filter_reads,
                                        const double read_keep_fraction,
                                        const uint32_t _read_length_min,
                                        const double _max_same_base_fraction,
                                        const double _max_N_fraction,
                                        const uint32_t _long_read_trigger_length,
                                        const uint32_t _long_read_split_length,
                                        const bool _long_read_distribute_remainder,
                                        const uint32_t constant_base_quality_replacement,
                                        const bool trim_mate_pair_junction_adapters,
                                        const uint32_t num_threads,
                                        const uint32_t long_read_pair_distance,
                                        const string &lp1_convert_file_name,
                                        const string &lp2_convert_file_name,
                                        AnalyzeFastqSummary* lp1_summary,
                                        AnalyzeFastqSummary* lp2_summary
                                        )
{
    cerr << "    Converting/filtering FASTQ file..." << endl;

    // EXPERIMENTAL: synthetic read pairs from long reads (see fastq.h)
    const bool make_long_read_pairs = (long_read_pair_distance != 0) && (_long_read_split_length != 0)
                                      && !lp1_convert_file_name.empty() && !lp2_convert_file_name.empty();
    // Stride between the two mates of a pair, in pieces. The outer span of a pair is
    // (stride + 1) pieces, so this is the stride whose span is closest to the requested distance.
    uint32_t long_read_pair_stride = 1;
    if (make_long_read_pairs) {
      uint32_t spanned_pieces = (long_read_pair_distance + _long_read_split_length / 2) / _long_read_split_length;
      long_read_pair_stride = (spanned_pieces > 2) ? (spanned_pieces - 1) : 1;
    }

    // Set up maps between formats
    map<string,uint8_t> format_to_chr_offset;
    format_to_chr_offset["SANGER"] = 33;
    format_to_chr_offset["SOLEXA"] = 64;
    format_to_chr_offset["ILLUMINA_1.3+"] = 64;

    // Honor that zero means no splitting
    uint32_t long_read_trigger_length = _long_read_trigger_length;
    if (long_read_trigger_length == 0) {
      long_read_trigger_length = numeric_limits<uint32_t>::max();
    }
    uint32_t long_read_split_length = _long_read_split_length;
    if (long_read_split_length == 0) {
      long_read_split_length = numeric_limits<uint32_t>::max();
    }
    bool long_read_distribute_remainder = _long_read_distribute_remainder;

    // Same thresholds as cFastqQualityConverter::predict_fastq_file_format(), but
    // operating on an already-known min raw quality score.
    auto predict_format_from_min_quality_score = [&](uint8_t min_q) -> string {
      string format = "SANGER";
      if (min_q >= format_to_chr_offset["SOLEXA"] - 5) {
        format = "SOLEXA";
      }
      if (min_q >= format_to_chr_offset["ILLUMINA_1.3+"]) {
        format = "ILLUMINA_1.3+";
      }
      return format;
    };

    // Per-read processing: quality convert, trim, split, filter, write, and update stats.
    struct AttemptResult {
      uint64_t num_original_reads = 0;
      uint64_t num_original_bases = 0;
      uint8_t  overall_min_quality_score = 255;
      uint8_t  overall_max_quality_score = 0;
      uint64_t num_reads = 0;
      uint64_t num_bases = 0;
      uint64_t num_filtered_too_short_reads = 0;
      uint64_t num_filtered_too_short_bases = 0;
      uint64_t num_junction_adapter_trimmed_reads = 0;
      uint64_t num_junction_adapter_trimmed_bases = 0;
      uint64_t num_filtered_same_base_reads = 0;
      uint64_t num_filtered_same_base_bases = 0;
      uint64_t num_filtered_too_many_N_reads = 0;
      uint64_t num_filtered_too_many_N_bases = 0;
      uint32_t read_length_min = numeric_limits<uint32_t>::max();
      uint32_t read_length_max = 0;
      bool     file_has_split_reads = false;
      uint64_t num_filtered_coverage_limit_reads = 0;
      uint64_t num_filtered_coverage_limit_bases = 0;
      // Synthetic long-read pairs: [0] is mate 1, [1] is mate 2. NOT included in num_reads/num_bases.
      uint64_t lp_num_reads[2] = {0, 0};
      uint64_t lp_num_bases[2] = {0, 0};
      uint32_t lp_read_length_min[2] = {numeric_limits<uint32_t>::max(), numeric_limits<uint32_t>::max()};
      uint32_t lp_read_length_max[2] = {0, 0};
    };

    const char constant_base_quality_replacement_char = static_cast<char>(format_to_chr_offset["SANGER"] + constant_base_quality_replacement);

    // Decided from the first reads, below, before any pass runs
    bool is_mate_pair_library = false;
    bool trim_junction_adapter = false;
    double junction_adapter_fraction = 0.0;
    string junction_adapter_name;

    auto run_pass = [&](const string& pass_quality_format,
                         const bool replace_base_qualities,
                         std::function<bool(cFastqSequence&)> get_next_sequence) -> AttemptResult
    {
      AttemptResult r;

      cFastqQualityConverter fqc(pass_quality_format, "SANGER");
      cFastqFile output_fastq_file(convert_file_name.c_str(), fstream::out, num_threads);
      unique_ptr<cFastqFile> lp1_fastq_file, lp2_fastq_file;
      if (make_long_read_pairs) {
        lp1_fastq_file.reset(new cFastqFile(lp1_convert_file_name.c_str(), fstream::out, num_threads));
        lp2_fastq_file.reset(new cFastqFile(lp2_convert_file_name.c_str(), fstream::out, num_threads));
      }
      // Pieces of the current long read that passed the filters, by piece index (pairing mode only)
      vector<cFastqSequence> kept_pieces;
      vector<bool> piece_was_kept;

      // Started afresh in each pass, so a second pass with a corrected quality format keeps the same reads
      cReadSubsampler subsampler(read_keep_fraction);
      uint32_t on_read = 1;
      cFastqSequence original_sequence;

      while (get_next_sequence(original_sequence)) {

        r.num_original_reads++;
        r.num_original_bases += original_sequence.length();

        // track raw (pre-conversion) min/max quality score
        for (uint32_t i=0; i<original_sequence.m_qualities.size(); i++) {
          uint8_t this_score = static_cast<uint8_t>(original_sequence.m_qualities[i]);
          r.overall_min_quality_score = min(r.overall_min_quality_score, this_score);
          r.overall_max_quality_score = max(r.overall_max_quality_score, this_score);
        }

        // --limit-fold-coverage: a read that is not chosen still counts toward the original totals
        // and the quality format above, and still takes its number, so that read names give the
        // read's position in the input whatever the limit was.
        if (!subsampler.keep()) {
          r.num_filtered_coverage_limit_reads++;
          r.num_filtered_coverage_limit_bases += original_sequence.length();
          on_read++;
          continue;
        }

        // truncate second name
        original_sequence.m_name_plus = "";

        // fastq quality convert
        fqc.convert_sequence(original_sequence);

        // Every base in the file has the same score: it carries no information (see fastq.h).
        // Replaced before trimming, which would otherwise act on the uninformative scores.
        if (replace_base_qualities) {
          original_sequence.m_qualities.assign(original_sequence.m_qualities.size(), constant_base_quality_replacement_char);
        }

        // Mate-pair library: cut the read at the junction adapter (see fastq.h)
        if (trim_junction_adapter) {
          size_t removed = trim_mate_pair_junction_adapter(original_sequence);
          if (removed) {
            r.num_junction_adapter_trimmed_reads++;
            r.num_junction_adapter_trimmed_bases += removed;
          }
        }

        // trim bad quality scores from the end
        if (trim_end_on_base_quality) {
          fastq_sequence_trim_end_on_base_quality(original_sequence, trim_end_on_base_quality);
        }

        // New loop that enables us to split long reads

        // Uniformly name, to prevent problems drawing alignments
        // and allows us to know the input order when merge sorting later
        original_sequence.m_name = to_string(file_index) + ":" + to_string(on_read++);

        // Decide on a per-read basis whether this read is subject to splitting
        uint32_t effective_split_length = (original_sequence.m_sequence.length() < long_read_trigger_length)
                                            ? numeric_limits<uint32_t>::max()
                                            : long_read_split_length;
        if (effective_split_length != numeric_limits<uint32_t>::max()) {
          r.file_has_split_reads = true;
        }

        // Always have at least one piece and only add new pieces if it is over the length (hence minus 1)
        uint32_t num_split_read_pieces = 1 + (original_sequence.m_sequence.length()-1) / effective_split_length;

        size_t chunk_start_0, chunk_end_0;
        double chunk_size = effective_split_length;
        bool read_was_split = num_split_read_pieces > 1;

        if (read_was_split) {
          if (!long_read_distribute_remainder) {
            // One fewer pieces because we ignore the last if not distributing the remainder
            num_split_read_pieces--;
          } else {
            // Update the chunk size (can be fractional) if we are distributing the remainder
            chunk_size = static_cast<double>(original_sequence.m_sequence.length()) / static_cast<double>(num_split_read_pieces);
          }
        }

        // Counts and writes one read to a given output. Both mates of a synthetic pair keep the long
        // read's strand: the pairs are FF, which the pair-evidence code reads through pair_geometry.
        auto emit_single = [&](const cFastqSequence& seq) {
          r.num_reads++;
          r.num_bases += seq.m_sequence.length();
          r.read_length_min = min<size_t>(seq.m_sequence.length(), r.read_length_min);
          r.read_length_max = max<size_t>(seq.m_sequence.length(), r.read_length_max);
          output_fastq_file.write_sequence(seq);
        };
        auto emit_mate = [&](const cFastqSequence& piece, uint32_t mate_index, const string& mate_name) {
          cFastqSequence seq = piece;
          seq.m_name = mate_name;
          r.lp_num_reads[mate_index]++;
          r.lp_num_bases[mate_index] += seq.m_sequence.length();
          r.lp_read_length_min[mate_index] = min<size_t>(seq.length(), r.lp_read_length_min[mate_index]);
          r.lp_read_length_max[mate_index] = max<size_t>(seq.length(), r.lp_read_length_max[mate_index]);
          ((mate_index == 0) ? lp1_fastq_file : lp2_fastq_file)->write_sequence(seq);
        };

        const bool pair_this_read = make_long_read_pairs && read_was_split;
        if (pair_this_read) {
          kept_pieces.assign(num_split_read_pieces, cFastqSequence());
          piece_was_kept.assign(num_split_read_pieces, false);
        }

        for (uint32_t i=0; i<num_split_read_pieces; i++) {

          cFastqSequence on_sequence;
          chunk_start_0 = ceil(i * chunk_size);
          chunk_end_0 = ceil((i+1) * chunk_size) - 1;

          // Suffix name with chunk number
          on_sequence.m_name = original_sequence.m_name;
          if (read_was_split) on_sequence.m_name += "S" + to_string(i+1);

          // Find correct chunk of bases and quals
          on_sequence.m_sequence = original_sequence.m_sequence.substr(chunk_start_0, chunk_end_0 - chunk_start_0 + 1);
          on_sequence.m_qualities = original_sequence.m_qualities.substr(chunk_start_0, chunk_end_0 - chunk_start_0 + 1);

          if ( filter_reads ) {

            // Discard sequences that are too short
            if ( _read_length_min && (on_sequence.length() < _read_length_min) ) {
              r.num_filtered_too_short_reads++;
              r.num_filtered_too_short_bases += on_sequence.length();
              continue;
            }

            // If doing tests that require, copy over info that is normally calculated on reading
            if (_max_N_fraction || _max_same_base_fraction) {

              if (read_was_split) {
                // Have to recalculate these since they are only done on reading the FASTQ right now
                for(size_t j=0; j<on_sequence.length(); j++) {
                  on_sequence.m_base_counts[basechar2index(on_sequence.m_sequence[j])]++;
                }
              } else {
                // Copy over if we didn't split the read
                for (uint8_t b=0; b<base_list_including_N_size; b++) {
                  on_sequence.m_base_counts[b] = original_sequence.m_base_counts[b];
                }
              }
            }

            // Discard sequences that are 50% or more N.
            if ( _max_N_fraction ) {
              if (_max_N_fraction * static_cast<double>(on_sequence.length()) <= on_sequence.m_base_counts[base_list_N_index]) {
                r.num_filtered_too_many_N_reads++;
                r.num_filtered_too_many_N_bases += on_sequence.length();
                continue;
              }
            }

            // Ignore heavily homopolymer reads, as these are a common type of machine error
            // Discard sequences that are 90% or more of a single base or N.
            bool same_base_filtered = false;
            if (_max_same_base_fraction) {
              for (uint8_t b=0; b<base_list_including_N_size; b++) {
                if ((_max_same_base_fraction * static_cast<double>(on_sequence.length())) <=
                    static_cast<double>(on_sequence.m_base_counts[b] + on_sequence.m_base_counts[base_list_N_index]) ) {
                  same_base_filtered = true;
                  break;
                }
              }
            }

            if (same_base_filtered)  {
              r.num_filtered_same_base_reads++;
              r.num_filtered_same_base_bases += on_sequence.length();
              continue;
            }

          } // end filter read block

          if (pair_this_read) {
            // Held back until we know which pieces of this read have a partner
            kept_pieces[i] = on_sequence;
            piece_was_kept[i] = true;
          } else {
            emit_single(on_sequence);
          }
        }

        // Pair up the pieces of a long read. Pieces are taken in blocks of 2*stride: within a
        // block, piece i pairs with piece i+stride, so each piece belongs to at most one pair and
        // every pair spans the same distance (for 4 pieces and a stride of 2: 1-3 and 2-4). A
        // final partial block pairs whatever still has a partner that far downstream. Anything
        // left over -- no partner, or a partner that was filtered -- is written as a single.
        if (pair_this_read) {
          vector<int32_t> partner(num_split_read_pieces, -1);
          for (uint32_t block_start = 0; block_start < num_split_read_pieces; block_start += 2 * long_read_pair_stride) {
            for (uint32_t i = block_start; (i < block_start + long_read_pair_stride) && (i + long_read_pair_stride < num_split_read_pieces); i++) {
              uint32_t j = i + long_read_pair_stride;
              if (piece_was_kept[i] && piece_was_kept[j]) {
                partner[i] = static_cast<int32_t>(j);
                partner[j] = static_cast<int32_t>(i);
              }
            }
          }

          // Both mates are named for the pair's FIRST piece, so that they share one (read, piece)
          // index -- that tuple is what re-joins them after alignment -- and both files stay in
          // increasing order of it. The prefixes are the two files' indices.
          for (uint32_t i=0; i<num_split_read_pieces; i++) {
            if (!piece_was_kept[i]) continue;
            if (partner[i] == -1) {
              emit_single(kept_pieces[i]);
            } else if (partner[i] > static_cast<int32_t>(i)) {
              string pair_suffix = original_sequence.m_name.substr(original_sequence.m_name.find(':')) + "S" + to_string(i+1);
              emit_mate(kept_pieces[i], 0, to_string(file_index + 1) + pair_suffix);
              emit_mate(kept_pieces[partner[i]], 1, to_string(file_index + 2) + pair_suffix);
            }
          }
        }
      }

      return r;
    };

    // ---- Sniff up to FASTQ_FORMAT_SNIFF_READS reads to predict the quality format ----
    const size_t FASTQ_FORMAT_SNIFF_READS = 10000;
    const int32_t prelim_offset = 64;

    cFastqFile input_fastq_file(file_name.c_str(), fstream::in, num_threads);
    cFastqQualityConverter prelim_fqc("ILLUMINA_1.3+", "SANGER");

    vector<cFastqSequence> buffered_reads;
    buffered_reads.reserve(FASTQ_FORMAT_SNIFF_READS);

    uint8_t sniffed_min_quality_score = 255;
    uint8_t sniffed_max_quality_score = 0;

    // The read type is judged from the same reads, as they are in the file: before the names
    // are replaced, the qualities converted or replaced, and long reads split.
    cReadTypeDetector read_type_detector;
    bool sniffed_numerical_qualities = false;
    cMatePairJunctionAdapterCounter junction_adapter_counter;

    {
      cFastqSequence seq;
      while ( (buffered_reads.size() < FASTQ_FORMAT_SNIFF_READS) && input_fastq_file.read_sequence(seq, prelim_fqc) ) {
        junction_adapter_counter.add(seq.m_sequence);
        for (uint32_t i=0; i<seq.m_qualities.size(); i++) {
          uint8_t this_score = static_cast<uint8_t>(seq.m_qualities[i]);
          sniffed_min_quality_score = min(sniffed_min_quality_score, this_score);
          sniffed_max_quality_score = max(sniffed_max_quality_score, this_score);
        }
        read_type_detector.add(seq.m_name, seq.m_name_comment, seq.m_qualities, seq.length());
        sniffed_numerical_qualities = sniffed_numerical_qualities || seq.m_numerical_qualities;
        buffered_reads.push_back(seq);
        seq = cFastqSequence();
      }
    }

    string quality_format = predict_format_from_min_quality_score(sniffed_min_quality_score);

    // A mate-pair library, judged from the same reads (see fastq.h). Set before the first pass,
    // which captures it, and not revisited: a second pass must trim the same reads.
    junction_adapter_fraction = junction_adapter_counter.fraction(junction_adapter_name);
    is_mate_pair_library = (junction_adapter_fraction >= kMatePairJunctionAdapterFraction);
    trim_junction_adapter = is_mate_pair_library && trim_mate_pair_junction_adapters;

    // Adjustment for buffered reads that were read in numerical-quality format using
    // the preliminary ILLUMINA_1.3+ offset (64), but should have used the chosen format's offset.
    int32_t numerical_quality_offset_adjustment = prelim_offset - static_cast<int32_t>(format_to_chr_offset[quality_format]);

    size_t buffered_idx = 0;
    auto get_buffered_then_stream = [&](cFastqSequence& seq) -> bool {
      if (buffered_idx < buffered_reads.size()) {
        seq = std::move(buffered_reads[buffered_idx++]);
        if (seq.m_numerical_qualities && numerical_quality_offset_adjustment != 0) {
          for (uint32_t i=0; i<seq.m_qualities.size(); i++) {
            seq.m_qualities[i] = static_cast<char>(static_cast<int32_t>(static_cast<uint8_t>(seq.m_qualities[i])) - numerical_quality_offset_adjustment);
          }
        }
        return true;
      }
      cFastqQualityConverter fqc_for_offset(quality_format, "SANGER");
      return input_fastq_file.read_sequence(seq, fqc_for_offset);
    };

    // Whether every base has one quality score is decided from the same reads and checked
    // against the whole file below, so that a file whose scores are informative is read only once
    bool replace_base_qualities = (constant_base_quality_replacement != 0) && (sniffed_min_quality_score == sniffed_max_quality_score);

    AttemptResult result = run_pass(quality_format, replace_base_qualities, get_buffered_then_stream);

    // ---- Check whether the format predicted from the sniff matches the format
    // derived from the actual min/max quality scores seen over the whole file ----
    string corrected_quality_format = predict_format_from_min_quality_score(result.overall_min_quality_score);
    const bool file_has_constant_base_quality = (result.num_original_bases > 0) && (result.overall_min_quality_score == result.overall_max_quality_score);
    const bool corrected_replace_base_qualities = (constant_base_quality_replacement != 0) && file_has_constant_base_quality;

    unique_ptr<cFastqFile> fresh_input_fastq_file;

    if ((corrected_quality_format != quality_format) || (corrected_replace_base_qualities != replace_base_qualities)) {

      if (corrected_quality_format != quality_format) {
        cerr << "    Warning: Quality score format predicted from the first " << FASTQ_FORMAT_SNIFF_READS << " reads (" << quality_format
             << ") did not match the format determined from the entire file (" << corrected_quality_format << ")." << endl;
      }
      if (corrected_replace_base_qualities != replace_base_qualities) {
        cerr << "    Warning: Every base in the first " << FASTQ_FORMAT_SNIFF_READS << " reads had the same quality score, but the entire file has more than one." << endl;
      }
      cerr << "    Re-converting FASTQ file using corrected quality score format..." << endl;

      quality_format = corrected_quality_format;
      replace_base_qualities = corrected_replace_base_qualities;

      fresh_input_fastq_file.reset(new cFastqFile(file_name.c_str(), fstream::in, num_threads));
      cFastqQualityConverter fresh_read_fqc(quality_format, "SANGER");

      auto get_fresh = [&](cFastqSequence& seq) -> bool {
        return fresh_input_fastq_file->read_sequence(seq, fresh_read_fqc);
      };

      result = run_pass(quality_format, replace_base_qualities, get_fresh);
    }

    // Every pass reads the whole file, so the original totals already include the reads that
    // --limit-fold-coverage did not choose.
    uint64_t num_filtered_coverage_limit_reads = result.num_filtered_coverage_limit_reads;
    uint64_t num_filtered_coverage_limit_bases = result.num_filtered_coverage_limit_bases;

    uint64_t num_original_reads = result.num_original_reads;
    uint64_t num_original_bases = result.num_original_bases;

    // Convert the overall min/max raw quality scores to SANGER for reporting
    cFastqQualityConverter final_fqc(quality_format, "SANGER");
    cFastqSequence min_max_sequence;
    min_max_sequence.m_qualities.append(1, (char)result.overall_min_quality_score);
    min_max_sequence.m_qualities.append(1, (char)result.overall_max_quality_score);
    final_fqc.convert_sequence(min_max_sequence);
    uint8_t min_quality_score = (uint8_t)min_max_sequence.m_qualities[0] - format_to_chr_offset["SANGER"];
    uint8_t max_quality_score = (uint8_t)min_max_sequence.m_qualities[1] - format_to_chr_offset["SANGER"];

    // The converted reads carry the replacement, and stage 07 sizes its error table from these
    int32_t replaced_constant_base_quality = -1;
    if (replace_base_qualities) {
      replaced_constant_base_quality = min_quality_score;
      min_quality_score = static_cast<uint8_t>(constant_base_quality_replacement);
      max_quality_score = static_cast<uint8_t>(constant_base_quality_replacement);
    }

    uint32_t width_for_reads = to_string(num_original_reads).size();
    uint32_t width_for_bases = to_string(num_original_bases).size();

    cerr << "    Original base quality format: " << quality_format << " New format: SANGER"<< endl;
    if (replace_base_qualities) {
      warn_constant_base_quality(file_name, replaced_constant_base_quality, constant_base_quality_replacement);
    }
    cerr << "    Original reads: " << num_original_reads << " bases: "<< num_original_bases << endl;
    if (is_mate_pair_library) warn_mate_pair_junction_adapter(file_name, junction_adapter_fraction, junction_adapter_name, trim_junction_adapter);
    if (trim_junction_adapter) {
      cerr << "    Mate-pair junction adapter trimmed from reads: " << result.num_junction_adapter_trimmed_reads
           << " bases removed: " << result.num_junction_adapter_trimmed_bases << endl;
    }

    if (num_filtered_coverage_limit_reads) {
      cerr << "    Filtered reads: " << setw(width_for_reads) << num_filtered_coverage_limit_reads;
      cerr << " bases: "<< setw(width_for_bases) << num_filtered_coverage_limit_bases;
      cerr << " (coverage limit option)" << endl;
    }

    // Then report splitting because other filters happened on split reads
    if (result.file_has_split_reads) {
      uint64_t total_split_reads = result.num_reads + result.num_filtered_same_base_reads + result.num_filtered_too_many_N_reads;
      uint64_t total_split_bases = result.num_bases + result.num_filtered_same_base_bases + result.num_filtered_too_many_N_bases;

      cerr << "    >> Long reads split to " << (long_read_distribute_remainder ? "≤": "exactly ");
      cout << long_read_split_length << " bases" << (long_read_distribute_remainder ? "" : " (extra bases discarded)") << endl;
      cerr << "    >> Split reads: " << setw(width_for_reads) << total_split_reads;
      cerr << " bases: "<< setw(width_for_bases) << total_split_bases << endl;
    }

    if (filter_reads) {
      if (result.num_filtered_too_many_N_reads + result.num_filtered_same_base_reads + result.num_filtered_too_short_reads + num_filtered_coverage_limit_reads == 0) {
        cerr << "    Filtered reads: none" << endl;
      } else {

        if (result.num_filtered_too_short_reads) {
          cerr << "    Filtered reads: " << setw(width_for_reads) << result.num_filtered_too_short_reads;
          cerr << " bases: "<< setw(width_for_bases) << result.num_filtered_too_short_bases;
          cerr << " (<" << _read_length_min << " bases long)" << endl;
        }
        if (result.num_filtered_too_many_N_reads) {
          string percentage = formatted_double(100 * _max_N_fraction, 1).to_string();
          cerr << "    Filtered reads: " << setw(width_for_reads) << result.num_filtered_too_many_N_reads;
          cerr << " bases: "<< setw(width_for_bases) << result.num_filtered_too_many_N_bases;
          cerr << " (≥" << percentage << "% bases N)" << endl;
        }
        if (result.num_filtered_same_base_reads) {
          string percentage = formatted_double(100 * _max_same_base_fraction, 1).to_string();
          cerr << "    Filtered reads: " << setw(width_for_reads) << result.num_filtered_same_base_reads;
          cerr << " bases: "<< setw(width_for_bases) << result.num_filtered_same_base_bases;
          cerr << " (≥" << percentage << "% same base)" << endl;
        }
      }
    }
    cerr << "    Analyzed reads: " << setw(width_for_reads) << result.num_reads << " bases: " << setw(width_for_bases) << result.num_bases << endl;

    double read_length_avg = (result.num_reads > 0) ? static_cast<double>(result.num_bases) / static_cast<double>(result.num_reads) : 0;

    AnalyzeFastqSummary retval(
                                 result.read_length_min,
                                 result.read_length_max,
                                 read_length_avg,
                                 num_original_reads,
                                 result.num_filtered_too_short_reads,
                                 result.num_filtered_same_base_reads,
                                 result.num_filtered_too_many_N_reads,
                                 num_filtered_coverage_limit_reads,
                                 result.num_reads,
                                 min_quality_score,
                                 max_quality_score,
                                 num_original_bases,
                                 result.num_bases,
                                 result.file_has_split_reads,
                                 quality_format,
                                 "SANGER",
                                 convert_file_name
                                 );
    retval.replaced_constant_base_quality = replaced_constant_base_quality;
    retval.junction_adapter_fraction = junction_adapter_fraction;
    retval.junction_adapter_name = junction_adapter_name;
    retval.junction_adapter_trimmed = trim_junction_adapter;
    retval.num_junction_adapter_trimmed_reads = result.num_junction_adapter_trimmed_reads;
    retval.num_junction_adapter_trimmed_bases = result.num_junction_adapter_trimmed_bases;
    // Numerical qualities were read with the preliminary offset whatever the format turned out to be
    read_type_detector.detect(sniffed_numerical_qualities ? static_cast<uint8_t>(prelim_offset) : format_to_chr_offset[quality_format],
                              retval.read_type, retval.read_type_source);

    // Synthetic long-read pairs get a summary per mate file. The original-read and filtered
    // counts stay with the file the reads came from (retval), so nothing is counted twice.
    if (make_long_read_pairs) {
      cerr << "    >> Long-read pairs: " << result.lp_num_reads[0] << " (mates " << (long_read_pair_stride + 1) * _long_read_split_length
           << " bases apart, outer span) bases: " << (result.lp_num_bases[0] + result.lp_num_bases[1]) << endl;
      AnalyzeFastqSummary* lp_summaries[2] = {lp1_summary, lp2_summary};
      const string lp_convert_file_names[2] = {lp1_convert_file_name, lp2_convert_file_name};
      for (uint32_t m = 0; m < 2; m++) {
        if (lp_summaries[m] == NULL) continue;
        *lp_summaries[m] = AnalyzeFastqSummary(
                                 (result.lp_num_reads[m] > 0) ? result.lp_read_length_min[m] : 0,
                                 result.lp_read_length_max[m],
                                 (result.lp_num_reads[m] > 0) ? static_cast<double>(result.lp_num_bases[m]) / static_cast<double>(result.lp_num_reads[m]) : 0,
                                 0, 0, 0, 0, 0,
                                 result.lp_num_reads[m],
                                 min_quality_score,
                                 max_quality_score,
                                 0,
                                 result.lp_num_bases[m],
                                 true,
                                 quality_format,
                                 "SANGER",
                                 lp_convert_file_names[m]
                                 );
        lp_summaries[m]->replaced_constant_base_quality = replaced_constant_base_quality;
        lp_summaries[m]->read_type = retval.read_type;
        lp_summaries[m]->read_type_source = retval.read_type_source;
      }
    }
    return retval;
  }
  
  pair<AnalyzeFastqSummary, AnalyzeFastqSummary> normalize_fastq_paired(
                                        const string &r1_file_name,
                                        const string &r1_convert_file_name,
                                        const string &r2_file_name,
                                        const string &r2_convert_file_name,
                                        const uint32_t r1_file_index,
                                        const uint32_t r2_file_index,
                                        const int32_t trim_end_on_base_quality,
                                        const bool filter_reads,
                                        const double read_keep_fraction,
                                        const uint32_t _read_length_min,
                                        const double _max_same_base_fraction,
                                        const double _max_N_fraction,
                                        const uint32_t _long_read_trigger_length,
                                        const uint32_t _long_read_split_length,
                                        const bool _long_read_distribute_remainder,
                                        const uint32_t constant_base_quality_replacement,
                                        const bool trim_mate_pair_junction_adapters,
                                        const uint32_t num_threads
                                        )
  {
    (void)_long_read_split_length;
    (void)_long_read_distribute_remainder;

    cerr << "    Converting/filtering paired FASTQ files..." << endl;

    map<string,uint8_t> format_to_chr_offset;
    format_to_chr_offset["SANGER"] = 33;
    format_to_chr_offset["SOLEXA"] = 64;
    format_to_chr_offset["ILLUMINA_1.3+"] = 64;

    auto predict_format_from_min_quality_score = [&](uint8_t min_q) -> string {
      string fmt = "SANGER";
      if (min_q >= format_to_chr_offset["SOLEXA"] - 5) fmt = "SOLEXA";
      if (min_q >= format_to_chr_offset["ILLUMINA_1.3+"]) fmt = "ILLUMINA_1.3+";
      return fmt;
    };

    // Per-file stats accumulated inside run_pass_paired
    struct FileResult {
      uint64_t num_original_reads = 0;
      uint64_t num_original_bases = 0;
      uint8_t  overall_min_quality_score = 255;
      uint8_t  overall_max_quality_score = 0;
      uint64_t num_reads = 0;
      uint64_t num_bases = 0;
      uint64_t num_filtered_too_short_reads = 0;
      uint64_t num_filtered_too_short_bases = 0;
      uint64_t num_junction_adapter_trimmed_reads = 0;
      uint64_t num_junction_adapter_trimmed_bases = 0;
      uint64_t num_filtered_same_base_reads = 0;
      uint64_t num_filtered_same_base_bases = 0;
      uint64_t num_filtered_too_many_N_reads = 0;
      uint64_t num_filtered_too_many_N_bases = 0;
      uint32_t read_length_min = numeric_limits<uint32_t>::max();
      uint32_t read_length_max = 0;
    };

    struct PairResult {
      FileResult r1, r2;
      // --limit-fold-coverage: pairs not chosen, counted per file
      uint64_t num_filtered_coverage_limit_reads[2] = {0, 0};
      uint64_t num_filtered_coverage_limit_bases[2] = {0, 0};
    };

    const char constant_base_quality_replacement_char = static_cast<char>(format_to_chr_offset["SANGER"] + constant_base_quality_replacement);

    // replace_base_qualities_r1/_r2: every base of that file has one quality score (see fastq.h)
    // Decided from the first reads of R1, below, before any pass runs (see fastq.h)
    bool is_mate_pair_library = false;
    bool trim_junction_adapter = false;
    double junction_adapter_fraction = 0.0;
    string junction_adapter_name;

    auto run_pass_paired = [&](
        const string& pass_quality_format,
        const bool replace_base_qualities_r1,
        const bool replace_base_qualities_r2,
        std::function<bool(cFastqSequence&)> get_next_r1,
        std::function<bool(cFastqSequence&)> get_next_r2) -> PairResult
    {
      PairResult pr;

      cFastqQualityConverter fqc(pass_quality_format, "SANGER");
      // Two output streams are open at once here, each with its own compressor
      // thread pool. Split the -j budget so the two pools sum to ~num_threads
      // rather than 2*num_threads. (At num_threads==1 this rounds up to 1 each,
      // i.e. 2 total; avoiding that last overshoot would require the two streams
      // to share a single pool.)
      unsigned int const per_stream_threads = std::max<unsigned int>(1, num_threads / 2);
      cFastqFile out_r1(r1_convert_file_name.c_str(), fstream::out, per_stream_threads);
      cFastqFile out_r2(r2_convert_file_name.c_str(), fstream::out, per_stream_threads);

      // One decision per pair, so mates are always kept or dropped together. Started afresh in each
      // pass, so a second pass with a corrected quality format keeps the same pairs.
      cReadSubsampler subsampler(read_keep_fraction);
      uint32_t on_read = 1;
      bool warned_long_reads = false;
      uint32_t long_read_trigger = (_long_read_trigger_length == 0) ? numeric_limits<uint32_t>::max() : _long_read_trigger_length;

      cFastqSequence orig_r1, orig_r2;

      while (get_next_r1(orig_r1)) {
        if (!get_next_r2(orig_r2)) {
          cerr << "  Warning: R2 file ended before R1 in paired FASTQ processing. Files may have different read counts." << endl;
          break;
        }

        // Count originals
        pr.r1.num_original_reads++;
        pr.r1.num_original_bases += orig_r1.length();
        pr.r2.num_original_reads++;
        pr.r2.num_original_bases += orig_r2.length();

        // Track raw quality scores (pre-conversion) for format detection
        for (uint32_t i = 0; i < orig_r1.m_qualities.size(); i++) {
          uint8_t q = static_cast<uint8_t>(orig_r1.m_qualities[i]);
          pr.r1.overall_min_quality_score = min(pr.r1.overall_min_quality_score, q);
          pr.r1.overall_max_quality_score = max(pr.r1.overall_max_quality_score, q);
        }
        for (uint32_t i = 0; i < orig_r2.m_qualities.size(); i++) {
          uint8_t q = static_cast<uint8_t>(orig_r2.m_qualities[i]);
          pr.r2.overall_min_quality_score = min(pr.r2.overall_min_quality_score, q);
          pr.r2.overall_max_quality_score = max(pr.r2.overall_max_quality_score, q);
        }

        // --limit-fold-coverage: a pair that is not chosen still counts toward the original totals
        // and the quality format above, and still takes its number (see normalize_fastq).
        if (!subsampler.keep()) {
          pr.num_filtered_coverage_limit_reads[0]++;
          pr.num_filtered_coverage_limit_bases[0] += orig_r1.length();
          pr.num_filtered_coverage_limit_reads[1]++;
          pr.num_filtered_coverage_limit_bases[1] += orig_r2.length();
          on_read++;
          continue;
        }

        orig_r1.m_name_plus = "";
        orig_r2.m_name_plus = "";

        fqc.convert_sequence(orig_r1);
        fqc.convert_sequence(orig_r2);

        // Replaced before trimming, which would otherwise act on the uninformative scores
        if (replace_base_qualities_r1) orig_r1.m_qualities.assign(orig_r1.m_qualities.size(), constant_base_quality_replacement_char);
        if (replace_base_qualities_r2) orig_r2.m_qualities.assign(orig_r2.m_qualities.size(), constant_base_quality_replacement_char);

        // Mate-pair library: cut each read at the junction adapter (see fastq.h). A read cut
        // shorter than --read-min-length takes its mate with it, below, like any short read.
        if (trim_junction_adapter) {
          size_t removed = trim_mate_pair_junction_adapter(orig_r1);
          if (removed) { pr.r1.num_junction_adapter_trimmed_reads++; pr.r1.num_junction_adapter_trimmed_bases += removed; }
          removed = trim_mate_pair_junction_adapter(orig_r2);
          if (removed) { pr.r2.num_junction_adapter_trimmed_reads++; pr.r2.num_junction_adapter_trimmed_bases += removed; }
        }

        if (trim_end_on_base_quality) {
          fastq_sequence_trim_end_on_base_quality(orig_r1, trim_end_on_base_quality);
          fastq_sequence_trim_end_on_base_quality(orig_r2, trim_end_on_base_quality);
        }

        // Name reads with their respective file indices so names remain unique
        orig_r1.m_name = to_string(r1_file_index) + ":" + to_string(on_read);
        orig_r2.m_name = to_string(r2_file_index) + ":" + to_string(on_read);
        on_read++;

        if (!warned_long_reads && (orig_r1.length() >= long_read_trigger || orig_r2.length() >= long_read_trigger)) {
          cerr << "    Warning: Long read(s) detected in paired FASTQ files. Long-read splitting is not supported in paired mode; reads will be processed without splitting." << endl;
          warned_long_reads = true;
        }

        // Filter: check R1 first, then R2. If either fails, skip both and
        // increment the triggered filter count in both R1 and R2 stats.
        if (filter_reads) {

          if (_read_length_min && (orig_r1.length() < _read_length_min)) {
            pr.r1.num_filtered_too_short_reads++;
            pr.r1.num_filtered_too_short_bases += orig_r1.length();
            pr.r2.num_filtered_too_short_reads++;
            pr.r2.num_filtered_too_short_bases += orig_r2.length();
            continue;
          }
          if (_read_length_min && (orig_r2.length() < _read_length_min)) {
            pr.r1.num_filtered_too_short_reads++;
            pr.r1.num_filtered_too_short_bases += orig_r1.length();
            pr.r2.num_filtered_too_short_reads++;
            pr.r2.num_filtered_too_short_bases += orig_r2.length();
            continue;
          }

          if (_max_N_fraction) {
            if (_max_N_fraction * static_cast<double>(orig_r1.length()) <= orig_r1.m_base_counts[base_list_N_index]) {
              pr.r1.num_filtered_too_many_N_reads++;
              pr.r1.num_filtered_too_many_N_bases += orig_r1.length();
              pr.r2.num_filtered_too_many_N_reads++;
              pr.r2.num_filtered_too_many_N_bases += orig_r2.length();
              continue;
            }
            if (_max_N_fraction * static_cast<double>(orig_r2.length()) <= orig_r2.m_base_counts[base_list_N_index]) {
              pr.r1.num_filtered_too_many_N_reads++;
              pr.r1.num_filtered_too_many_N_bases += orig_r1.length();
              pr.r2.num_filtered_too_many_N_reads++;
              pr.r2.num_filtered_too_many_N_bases += orig_r2.length();
              continue;
            }
          }

          if (_max_same_base_fraction) {
            bool r1_same_base = false;
            for (uint8_t b = 0; b < base_list_including_N_size; b++) {
              if (_max_same_base_fraction * static_cast<double>(orig_r1.length()) <=
                  static_cast<double>(orig_r1.m_base_counts[b] + orig_r1.m_base_counts[base_list_N_index])) {
                r1_same_base = true;
                break;
              }
            }
            if (r1_same_base) {
              pr.r1.num_filtered_same_base_reads++;
              pr.r1.num_filtered_same_base_bases += orig_r1.length();
              pr.r2.num_filtered_same_base_reads++;
              pr.r2.num_filtered_same_base_bases += orig_r2.length();
              continue;
            }

            bool r2_same_base = false;
            for (uint8_t b = 0; b < base_list_including_N_size; b++) {
              if (_max_same_base_fraction * static_cast<double>(orig_r2.length()) <=
                  static_cast<double>(orig_r2.m_base_counts[b] + orig_r2.m_base_counts[base_list_N_index])) {
                r2_same_base = true;
                break;
              }
            }
            if (r2_same_base) {
              pr.r1.num_filtered_same_base_reads++;
              pr.r1.num_filtered_same_base_bases += orig_r1.length();
              pr.r2.num_filtered_same_base_reads++;
              pr.r2.num_filtered_same_base_bases += orig_r2.length();
              continue;
            }
          }

        } // end filter block

        // Both reads pass all filters; write them
        pr.r1.num_reads++;
        pr.r1.num_bases += orig_r1.length();
        pr.r2.num_reads++;
        pr.r2.num_bases += orig_r2.length();

        pr.r1.read_length_min = min<size_t>(orig_r1.length(), pr.r1.read_length_min);
        pr.r1.read_length_max = max<size_t>(orig_r1.length(), pr.r1.read_length_max);
        pr.r2.read_length_min = min<size_t>(orig_r2.length(), pr.r2.read_length_min);
        pr.r2.read_length_max = max<size_t>(orig_r2.length(), pr.r2.read_length_max);

        out_r1.write_sequence(orig_r1);
        out_r2.write_sequence(orig_r2);
      }

      return pr;
    };

    // ---- Sniff up to FASTQ_FORMAT_SNIFF_READS reads from R1 to predict the quality format ----
    const size_t FASTQ_FORMAT_SNIFF_READS = 10000;
    const int32_t prelim_offset = 64;

    cFastqFile r1_input(r1_file_name.c_str(), fstream::in, num_threads);
    cFastqQualityConverter prelim_fqc("ILLUMINA_1.3+", "SANGER");

    vector<cFastqSequence> buffered_reads;
    buffered_reads.reserve(FASTQ_FORMAT_SNIFF_READS);

    uint8_t sniffed_min_quality_score = 255;
    uint8_t sniffed_max_quality_score = 0;

    // Judged from R1 alone, like the quality format: the two files of a pair are one library
    cReadTypeDetector read_type_detector;
    bool sniffed_numerical_qualities = false;
    cMatePairJunctionAdapterCounter junction_adapter_counter;

    {
      cFastqSequence seq;
      while ((buffered_reads.size() < FASTQ_FORMAT_SNIFF_READS) && r1_input.read_sequence(seq, prelim_fqc)) {
        junction_adapter_counter.add(seq.m_sequence);
        for (uint32_t i = 0; i < seq.m_qualities.size(); i++) {
          uint8_t q = static_cast<uint8_t>(seq.m_qualities[i]);
          sniffed_min_quality_score = min(sniffed_min_quality_score, q);
          sniffed_max_quality_score = max(sniffed_max_quality_score, q);
        }
        read_type_detector.add(seq.m_name, seq.m_name_comment, seq.m_qualities, seq.length());
        sniffed_numerical_qualities = sniffed_numerical_qualities || seq.m_numerical_qualities;
        buffered_reads.push_back(seq);
        seq = cFastqSequence();
      }
    }

    string quality_format = predict_format_from_min_quality_score(sniffed_min_quality_score);

    int32_t numerical_quality_offset_adjustment = prelim_offset - static_cast<int32_t>(format_to_chr_offset[quality_format]);

    // A mate-pair library, judged from R1 like the quality format: both files are one library
    junction_adapter_fraction = junction_adapter_counter.fraction(junction_adapter_name);
    is_mate_pair_library = (junction_adapter_fraction >= kMatePairJunctionAdapterFraction);
    trim_junction_adapter = is_mate_pair_library && trim_mate_pair_junction_adapters;

    size_t buffered_idx = 0;
    auto get_buffered_then_stream_r1 = [&](cFastqSequence& seq) -> bool {
      if (buffered_idx < buffered_reads.size()) {
        seq = std::move(buffered_reads[buffered_idx++]);
        if (seq.m_numerical_qualities && numerical_quality_offset_adjustment != 0) {
          for (uint32_t i = 0; i < seq.m_qualities.size(); i++) {
            seq.m_qualities[i] = static_cast<char>(static_cast<int32_t>(static_cast<uint8_t>(seq.m_qualities[i])) - numerical_quality_offset_adjustment);
          }
        }
        return true;
      }
      cFastqQualityConverter fqc_for_offset(quality_format, "SANGER");
      return r1_input.read_sequence(seq, fqc_for_offset);
    };

    cFastqFile r2_input(r2_file_name.c_str(), fstream::in, num_threads);
    cFastqQualityConverter r2_initial_fqc(quality_format, "SANGER");
    auto get_r2_stream = [&](cFastqSequence& seq) -> bool {
      return r2_input.read_sequence(seq, r2_initial_fqc);
    };

    // Only R1 is sniffed, so R2 is first assumed to be like it; both are checked over the whole files
    // below, and converted again if either guess was wrong.
    bool replace_base_qualities_r1 = (constant_base_quality_replacement != 0) && (sniffed_min_quality_score == sniffed_max_quality_score);
    bool replace_base_qualities_r2 = replace_base_qualities_r1;

    PairResult result = run_pass_paired(quality_format, replace_base_qualities_r1, replace_base_qualities_r2, get_buffered_then_stream_r1, get_r2_stream);

    // ---- Check format from R1's actual min quality score (same logic as normalize_fastq) ----
    string corrected_quality_format = predict_format_from_min_quality_score(result.r1.overall_min_quality_score);
    auto file_has_constant_base_quality = [](const FileResult& fr) -> bool {
      return (fr.num_original_bases > 0) && (fr.overall_min_quality_score == fr.overall_max_quality_score);
    };
    const bool corrected_replace_base_qualities_r1 = (constant_base_quality_replacement != 0) && file_has_constant_base_quality(result.r1);
    const bool corrected_replace_base_qualities_r2 = (constant_base_quality_replacement != 0) && file_has_constant_base_quality(result.r2);
    const bool replacement_guess_wrong = (corrected_replace_base_qualities_r1 != replace_base_qualities_r1)
                                      || (corrected_replace_base_qualities_r2 != replace_base_qualities_r2);

    unique_ptr<cFastqFile> fresh_r1, fresh_r2;

    if ((corrected_quality_format != quality_format) || replacement_guess_wrong) {
      if (corrected_quality_format != quality_format) {
        cerr << "    Warning: Quality score format predicted from the first " << FASTQ_FORMAT_SNIFF_READS << " reads (" << quality_format
             << ") did not match the format determined from the entire file (" << corrected_quality_format << ")." << endl;
      }
      if (replacement_guess_wrong) {
        cerr << "    Warning: Whether every base has the same quality score, judged from the first " << FASTQ_FORMAT_SNIFF_READS
             << " reads of R1, did not hold for the entire R1 and R2 files." << endl;
      }
      cerr << "    Re-converting FASTQ files using corrected quality score format..." << endl;

      quality_format = corrected_quality_format;
      replace_base_qualities_r1 = corrected_replace_base_qualities_r1;
      replace_base_qualities_r2 = corrected_replace_base_qualities_r2;

      fresh_r1.reset(new cFastqFile(r1_file_name.c_str(), fstream::in, num_threads));
      fresh_r2.reset(new cFastqFile(r2_file_name.c_str(), fstream::in, num_threads));

      cFastqQualityConverter fresh_r1_fqc(quality_format, "SANGER");
      cFastqQualityConverter fresh_r2_fqc(quality_format, "SANGER");

      auto get_fresh_r1 = [&](cFastqSequence& seq) -> bool {
        return fresh_r1->read_sequence(seq, fresh_r1_fqc);
      };
      auto get_fresh_r2 = [&](cFastqSequence& seq) -> bool {
        return fresh_r2->read_sequence(seq, fresh_r2_fqc);
      };

      result = run_pass_paired(quality_format, replace_base_qualities_r1, replace_base_qualities_r2, get_fresh_r1, get_fresh_r2);
    }

    // Every pass reads both files to the end, so the original totals already include the pairs
    // that --limit-fold-coverage did not choose.
    uint64_t num_filtered_coverage_limit_reads_r1 = result.num_filtered_coverage_limit_reads[0];
    uint64_t num_filtered_coverage_limit_bases_r1 = result.num_filtered_coverage_limit_bases[0];
    uint64_t num_filtered_coverage_limit_reads_r2 = result.num_filtered_coverage_limit_reads[1];
    uint64_t num_filtered_coverage_limit_bases_r2 = result.num_filtered_coverage_limit_bases[1];

    // Convert raw min/max quality scores to SANGER for reporting
    cFastqQualityConverter final_fqc(quality_format, "SANGER");

    auto build_summary = [&](const FileResult& fr,
                              uint64_t ncov_reads, uint64_t ncov_bases,
                              const string& convert_name,
                              const bool replaced_base_qualities) -> AnalyzeFastqSummary {
      cFastqSequence mm_seq;
      mm_seq.m_qualities.append(1, (char)fr.overall_min_quality_score);
      mm_seq.m_qualities.append(1, (char)fr.overall_max_quality_score);
      final_fqc.convert_sequence(mm_seq);
      uint8_t min_q = (uint8_t)mm_seq.m_qualities[0] - format_to_chr_offset["SANGER"];
      uint8_t max_q = (uint8_t)mm_seq.m_qualities[1] - format_to_chr_offset["SANGER"];

      // The converted reads carry the replacement, and stage 07 sizes its error table from these
      int32_t replaced_constant_base_quality = -1;
      if (replaced_base_qualities) {
        replaced_constant_base_quality = min_q;
        min_q = static_cast<uint8_t>(constant_base_quality_replacement);
        max_q = static_cast<uint8_t>(constant_base_quality_replacement);
      }

      // The original totals already include the reads --limit-fold-coverage did not choose
      (void)ncov_bases;
      uint64_t num_orig_reads = fr.num_original_reads;
      uint64_t num_orig_bases = fr.num_original_bases;
      double avg = (fr.num_reads > 0) ? static_cast<double>(fr.num_bases) / static_cast<double>(fr.num_reads) : 0.0;

      AnalyzeFastqSummary afs(
          fr.read_length_min,
          fr.read_length_max,
          avg,
          num_orig_reads,
          fr.num_filtered_too_short_reads,
          fr.num_filtered_same_base_reads,
          fr.num_filtered_too_many_N_reads,
          ncov_reads,
          fr.num_reads,
          min_q,
          max_q,
          num_orig_bases,
          fr.num_bases,
          false,   // reads_were_split: not supported in paired mode
          quality_format,
          "SANGER",
          convert_name
      );
      afs.replaced_constant_base_quality = replaced_constant_base_quality;
      afs.junction_adapter_fraction = junction_adapter_fraction;
      afs.junction_adapter_name = junction_adapter_name;
      afs.junction_adapter_trimmed = trim_junction_adapter;
      afs.num_junction_adapter_trimmed_reads = fr.num_junction_adapter_trimmed_reads;
      afs.num_junction_adapter_trimmed_bases = fr.num_junction_adapter_trimmed_bases;
      return afs;
    };

    AnalyzeFastqSummary s_r1 = build_summary(result.r1, num_filtered_coverage_limit_reads_r1, num_filtered_coverage_limit_bases_r1, r1_convert_file_name, replace_base_qualities_r1);
    AnalyzeFastqSummary s_r2 = build_summary(result.r2, num_filtered_coverage_limit_reads_r2, num_filtered_coverage_limit_bases_r2, r2_convert_file_name, replace_base_qualities_r2);
    read_type_detector.detect(sniffed_numerical_qualities ? static_cast<uint8_t>(prelim_offset) : format_to_chr_offset[quality_format],
                              s_r1.read_type, s_r1.read_type_source);
    s_r2.read_type = s_r1.read_type;
    s_r2.read_type_source = s_r1.read_type_source;

    // Print combined summary for paired files
    uint64_t total_orig_reads = s_r1.num_original_reads + s_r2.num_original_reads;
    uint64_t total_orig_bases = s_r1.num_original_bases + s_r2.num_original_bases;
    uint64_t total_reads = s_r1.num_reads + s_r2.num_reads;
    uint64_t total_bases = s_r1.num_bases + s_r2.num_bases;

    uint32_t width_for_reads = to_string(total_orig_reads).size();
    uint32_t width_for_bases = to_string(total_orig_bases).size();

    cerr << "    Original base quality format: " << quality_format << " New format: SANGER" << endl;
    if (replace_base_qualities_r1) warn_constant_base_quality(r1_file_name, s_r1.replaced_constant_base_quality, constant_base_quality_replacement);
    if (replace_base_qualities_r2) warn_constant_base_quality(r2_file_name, s_r2.replaced_constant_base_quality, constant_base_quality_replacement);
    cerr << "    Original reads: " << total_orig_reads << " bases: " << total_orig_bases << endl;
    if (is_mate_pair_library) warn_mate_pair_junction_adapter(r1_file_name, junction_adapter_fraction, junction_adapter_name, trim_junction_adapter);
    if (trim_junction_adapter) {
      cerr << "    Mate-pair junction adapter trimmed from reads: " << (s_r1.num_junction_adapter_trimmed_reads + s_r2.num_junction_adapter_trimmed_reads)
           << " bases removed: " << (s_r1.num_junction_adapter_trimmed_bases + s_r2.num_junction_adapter_trimmed_bases) << endl;
    }

    uint64_t total_cov_filtered = num_filtered_coverage_limit_reads_r1 + num_filtered_coverage_limit_reads_r2;
    if (total_cov_filtered) {
      cerr << "    Filtered reads: " << setw(width_for_reads) << total_cov_filtered;
      cerr << " bases: " << setw(width_for_bases) << (num_filtered_coverage_limit_bases_r1 + num_filtered_coverage_limit_bases_r2);
      cerr << " (coverage limit option)" << endl;
    }

    if (filter_reads) {
      uint64_t total_short = result.r1.num_filtered_too_short_reads + result.r2.num_filtered_too_short_reads;
      uint64_t total_N = result.r1.num_filtered_too_many_N_reads + result.r2.num_filtered_too_many_N_reads;
      uint64_t total_same = result.r1.num_filtered_same_base_reads + result.r2.num_filtered_same_base_reads;

      if (total_short + total_N + total_same + total_cov_filtered == 0) {
        cerr << "    Filtered reads: none" << endl;
      } else {
        if (total_short) {
          cerr << "    Filtered reads: " << setw(width_for_reads) << total_short;
          cerr << " bases: " << setw(width_for_bases) << (result.r1.num_filtered_too_short_bases + result.r2.num_filtered_too_short_bases);
          cerr << " (<" << _read_length_min << " bases long)" << endl;
        }
        if (total_N) {
          string percentage = formatted_double(100 * _max_N_fraction, 1).to_string();
          cerr << "    Filtered reads: " << setw(width_for_reads) << total_N;
          cerr << " bases: " << setw(width_for_bases) << (result.r1.num_filtered_too_many_N_bases + result.r2.num_filtered_too_many_N_bases);
          cerr << " (≥" << percentage << "% bases N)" << endl;
        }
        if (total_same) {
          string percentage = formatted_double(100 * _max_same_base_fraction, 1).to_string();
          cerr << "    Filtered reads: " << setw(width_for_reads) << total_same;
          cerr << " bases: " << setw(width_for_bases) << (result.r1.num_filtered_same_base_bases + result.r2.num_filtered_same_base_bases);
          cerr << " (≥" << percentage << "% same base)" << endl;
        }
      }
    }
    cerr << "    Analyzed reads: " << setw(width_for_reads) << total_reads << " bases: " << setw(width_for_bases) << total_bases << endl;

    return {s_r1, s_r2};
  }

  // ---- Mate-pair junction adapters (see fastq.h) ----

  const vector<cMatePairJunctionAdapter>& mate_pair_junction_adapters()
  {
    static vector<cMatePairJunctionAdapter> adapters;
    if (adapters.empty()) {
      cMatePairJunctionAdapter nextera;
      nextera.name = "Nextera mate-pair junction adapter";
      nextera.sequence = "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG";
      adapters.push_back(nextera);
      cMatePairJunctionAdapter cre_lox;
      cre_lox.name = "Cre-lox mate-pair junction adapter";
      cre_lox.sequence = "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG";
      adapters.push_back(cre_lox);
    }
    return adapters;
  }

  namespace {

    // One strand of one adapter, as it may be read
    struct JunctionAdapterPattern {
      size_t adapter_index;
      string sequence;
      string first_half;
      string second_half;
    };

    // Both strands of every adapter (one, for an adapter that is its own reverse complement)
    const vector<JunctionAdapterPattern>& junction_adapter_patterns()
    {
      static vector<JunctionAdapterPattern> patterns;
      if (patterns.empty()) {
        const vector<cMatePairJunctionAdapter>& adapters = mate_pair_junction_adapters();
        for (size_t i = 0; i < adapters.size(); i++) {
          vector<string> strands;
          strands.push_back(adapters[i].sequence);
          string rc = reverse_complement(adapters[i].sequence);
          if (rc != adapters[i].sequence) strands.push_back(rc);
          for (size_t s = 0; s < strands.size(); s++) {
            JunctionAdapterPattern p;
            p.adapter_index = i;
            p.sequence = strands[s];
            p.first_half = p.sequence.substr(0, p.sequence.size() / 2);
            p.second_half = p.sequence.substr(p.sequence.size() / 2);
            patterns.push_back(p);
          }
        }
      }
      return patterns;
    }

    // First position at which pattern matches sequence with at most max_mismatches; npos if none
    size_t find_with_mismatches(const string& sequence, const string& pattern, const uint32_t max_mismatches)
    {
      const size_t pattern_length = pattern.size();
      if (sequence.size() < pattern_length) return string::npos;
      const size_t last = sequence.size() - pattern_length;
      for (size_t i = 0; i <= last; i++) {
        uint32_t mismatches = 0;
        size_t j = 0;
        for (; j < pattern_length; j++) {
          if ((sequence[i + j] != pattern[j]) && (++mismatches > max_mismatches)) break;
        }
        if (j == pattern_length) return i;
      }
      return string::npos;
    }
  }

  int32_t exact_mate_pair_junction_adapter(const string& sequence)
  {
    const vector<JunctionAdapterPattern>& patterns = junction_adapter_patterns();
    for (size_t p = 0; p < patterns.size(); p++)
      if (sequence.find(patterns[p].sequence) != string::npos) return static_cast<int32_t>(patterns[p].adapter_index);
    return -1;
  }

  size_t mate_pair_junction_adapter_position(const string& sequence)
  {
    const vector<JunctionAdapterPattern>& patterns = junction_adapter_patterns();
    size_t cut = string::npos;

    for (size_t p = 0; p < patterns.size(); p++) {
      const JunctionAdapterPattern& pattern = patterns[p];

      // The first half: the adapter starts here
      size_t pos = find_with_mismatches(sequence, pattern.first_half, kJunctionAdapterMaxMismatches);
      if ((pos != string::npos) && ((cut == string::npos) || (pos < cut))) cut = pos;

      // The second half: the adapter started one first half before it, even if that half is too
      // damaged to recognize or the read begins inside it
      pos = find_with_mismatches(sequence, pattern.second_half, kJunctionAdapterMaxMismatches);
      if (pos != string::npos) {
        size_t start = (pos >= pattern.first_half.size()) ? pos - pattern.first_half.size() : 0;
        if ((cut == string::npos) || (start < cut)) cut = start;
      }
    }
    if (cut != string::npos) return cut;

    // The read ends partway into the first half of an adapter
    for (size_t p = 0; p < patterns.size(); p++) {
      const string& first_half = patterns[p].first_half;
      for (size_t length = first_half.size() - 1; length >= kJunctionAdapterMinEndMatch; length--) {
        if ((sequence.size() >= length) && (sequence.compare(sequence.size() - length, length, first_half, 0, length) == 0)) {
          size_t start = sequence.size() - length;
          if ((cut == string::npos) || (start < cut)) cut = start;
          break;
        }
      }
    }
    return cut;
  }

  bool contains_mate_pair_junction_adapter_kmer(const string& sequence)
  {
    if (sequence.size() < kJunctionAdapterKmer) return false;
    const vector<JunctionAdapterPattern>& patterns = junction_adapter_patterns();
    for (size_t p = 0; p < patterns.size(); p++) {
      const string& adapter = patterns[p].sequence;
      for (size_t i = 0; i + kJunctionAdapterKmer <= adapter.size(); i++) {
        if (sequence.find(adapter.c_str() + i, 0, kJunctionAdapterKmer) != string::npos) return true;
      }
    }
    return false;
  }

  namespace {
    // Cuts a read at a junction adapter. Returns the number of bases removed (0 = no adapter).
    size_t trim_mate_pair_junction_adapter(cFastqSequence& seq)
    {
      size_t cut = mate_pair_junction_adapter_position(seq.m_sequence);
      if (cut == string::npos) return 0;
      size_t removed = seq.m_sequence.size() - cut;
      seq.m_sequence.erase(cut);
      seq.m_qualities.erase(cut);
      // The N and same-base filters work from these, which are counted as the read is read
      memset(seq.m_base_counts, 0, sizeof(seq.m_base_counts));
      for (size_t j = 0; j < seq.m_sequence.size(); j++) seq.m_base_counts[basechar2index(seq.m_sequence[j])]++;
      return removed;
    }

    void warn_mate_pair_junction_adapter(const string& file_name, const double fraction, const string& adapter_name, const bool trimming)
    {
      WARN("Read file " + file_name + "\nlooks like a MATE-PAIR library: "
           + formatted_double(100.0 * fraction, 1).to_string() + "% of its first reads contain the\n"
           + adapter_name + ".\n"
           "A read through that adapter joins the two ends of a fragment that is several kb long, so it\n"
           "is not one piece of the genome. "
           + (trimming ? string("breseq will cut these reads at the adapter and keep the part\nbefore it.")
                       : string("These reads are being used as they are, because\n--no-mate-pair-adapter-trimming was given."))
           + " For better results, process the reads with a tool made for mate-pair libraries\n"
           "before running breseq.");
    }
  }

  // ---- Read types (see fastq.h) ----

  const char* const kReadTypeIllumina     = "illumina";
  const char* const kReadTypeNanopore     = "nanopore";
  const char* const kReadTypeNanoporeFast = "nanopore-fast";
  const char* const kReadTypeNanoporeHAC  = "nanopore-hac";
  const char* const kReadTypePacBio       = "pacbio";
  const char* const kReadTypePacBioCLR    = "pacbio-clr";
  const char* const kReadTypePacBioHiFi   = "pacbio-hifi";
  const char* const kReadTypeUnknown      = "unknown";

  const char* const kReadTypeSourceSpecified = "specified";
  const char* const kReadTypeSourceHeader    = "header";
  const char* const kReadTypeSourceQuality   = "quality";

  vector<string> read_type_names()
  {
    return make_vector<string>(kReadTypeIllumina)(kReadTypeNanopore)(kReadTypeNanoporeFast)(kReadTypeNanoporeHAC)
      (kReadTypePacBio)(kReadTypePacBioCLR)(kReadTypePacBioHiFi)(kReadTypeUnknown);
  }

  string canonical_read_type(const string& name)
  {
    // Synonyms are brefito's (LONG_READ_TYPE_SYNONYMS and its illumina-* layout names). breseq
    // works out pairing from the file names, so every Illumina layout is just "illumina".
    static map<string,string> synonyms;
    if (synonyms.empty()) {
      vector<string> names = read_type_names();
      for (vector<string>::const_iterator it = names.begin(); it != names.end(); it++) synonyms[*it] = *it;
      synonyms["ont"]          = kReadTypeNanopore;
      synonyms["ont-fast"]     = kReadTypeNanoporeFast;
      synonyms["ont-raw"]      = kReadTypeNanoporeFast;
      synonyms["nanopore-raw"] = kReadTypeNanoporeFast;
      synonyms["ont-hac"]      = kReadTypeNanoporeHAC;
      synonyms["ont-hq"]       = kReadTypeNanoporeHAC;
      synonyms["nanopore-hq"]  = kReadTypeNanoporeHAC;
      // Super accuracy is not told apart from high accuracy
      synonyms["nanopore-sup"] = kReadTypeNanoporeHAC;
      synonyms["ont-sup"]      = kReadTypeNanoporeHAC;
      synonyms["illumina-se"]     = kReadTypeIllumina;
      synonyms["illumina-pe"]     = kReadTypeIllumina;
      synonyms["illumina-r1"]     = kReadTypeIllumina;
      synonyms["illumina-r2"]     = kReadTypeIllumina;
      synonyms["illumina-pair"]   = kReadTypeIllumina;
      synonyms["illumina-paired"] = kReadTypeIllumina;
    }
    map<string,string>::const_iterator it = synonyms.find(to_lower(name));
    return (it == synonyms.end()) ? "" : it->second;
  }

  bool is_nanopore_read_type(const string& read_type)
  {
    return (read_type == kReadTypeNanopore) || (read_type == kReadTypeNanoporeFast)
        || (read_type == kReadTypeNanoporeHAC);
  }

  bool is_long_read_type(const string& read_type)
  {
    return is_nanopore_read_type(read_type) || (read_type == kReadTypePacBio)
        || (read_type == kReadTypePacBioCLR) || (read_type == kReadTypePacBioHiFi);
  }

  string read_type_platform(const string& read_type)
  {
    // The values the SAM specification lists for @RG PL
    if (read_type == kReadTypeIllumina) return "ILLUMINA";
    if (is_nanopore_read_type(read_type)) return "ONT";
    if (is_long_read_type(read_type)) return "PACBIO";
    return "";
  }

  string read_type_source_description(const string& read_type_source)
  {
    if (read_type_source == kReadTypeSourceSpecified) return "specified on the command line";
    if (read_type_source == kReadTypeSourceHeader) return "detected from read names";
    if (read_type_source == kReadTypeSourceQuality) return "guessed from read lengths and base qualities";
    return "not determined";
  }

  void sample_fastq_file_for_read_type(const string& file_name, cReadTypeDetector& detector)
  {
    cFastqFile input_fastq_file(file_name.c_str(), fstream::in);
    cFastqQualityConverter fqc("SANGER", "SANGER");
    cFastqSequence seq;
    while ((detector.num_reads() < kReadTypeSampleReads) && input_fastq_file.read_sequence(seq, fqc))
      detector.add(seq.m_name, seq.m_name_comment, seq.m_qualities, seq.length());
  }

  void resolve_read_type(const string& file_name, const string& specified_read_type,
                         string& read_type, string& read_type_source)
  {
    if (specified_read_type.empty()) {
      cerr << "    Read type: " << read_type;
      if (!read_type_source.empty()) cerr << " (" << read_type_source_description(read_type_source) << ")";
      cerr << endl;
      return;
    }

    // Different platforms are a contradiction whatever the evidence. Within a platform only read
    // names are worth contradicting the user with, and "nanopore"/"pacbio" contradict no subtype.
    const string specified_platform = read_type_platform(specified_read_type);
    const string detected_platform = read_type_platform(read_type);
    const bool specified_is_generic = (specified_read_type == kReadTypeNanopore) || (specified_read_type == kReadTypePacBio);
    const bool disagrees = !detected_platform.empty() && !specified_platform.empty()
      && ( (detected_platform != specified_platform)
        || ((read_type_source == kReadTypeSourceHeader) && !specified_is_generic && (read_type != specified_read_type)) );

    cerr << "    Read type: " << specified_read_type << " (" << read_type_source_description(kReadTypeSourceSpecified) << ")" << endl;
    if (disagrees) {
      WARN("Read file " + file_name + "\nwas given the read type " + specified_read_type + ", but its reads look like "
           + read_type + " (" + read_type_source_description(read_type_source) + ").\n"
           "The read type given on the command line will be used.");
    }

    read_type = specified_read_type;
    read_type_source = kReadTypeSourceSpecified;
  }

  string split_read_type_prefix(const string& read_file_argument, string& read_type)
  {
    read_type = "";
    size_t colon_pos = read_file_argument.find(":");
    if ((colon_pos == string::npos) || (colon_pos == 0) || (colon_pos + 1 == read_file_argument.size()))
      return read_file_argument;

    string prefix_type = canonical_read_type(read_file_argument.substr(0, colon_pos));
    if (prefix_type.empty()) return read_file_argument;

    // A file that really is named "illumina:x.fastq" is that file
    if (path_exists(read_file_argument)) return read_file_argument;

    read_type = prefix_type;
    return read_file_argument.substr(colon_pos + 1);
  }

  string strip_read_type_prefix(const string& read_file_argument)
  {
    string read_type;
    return split_read_type_prefix(read_file_argument, read_type);
  }

  cReadTypeDetector::cReadTypeDetector()
  : m_num_reads(0)
  {
    memset(m_header_votes, 0, sizeof(m_header_votes));
    memset(m_quality_chr_counts, 0, sizeof(m_quality_chr_counts));
  }

  namespace {

    bool all_digits(const string& s, size_t start = 0, size_t end = string::npos)
    {
      if (end == string::npos) end = s.size();
      if (start >= end) return false;
      for (size_t i = start; i < end; i++)
        if (!isdigit(static_cast<unsigned char>(s[i]))) return false;
      return true;
    }

    // m64011_190830_220126/101/ccs            HiFi (also .../ccs/fwd, by-strand)
    // m54119U_210220_034433/66/0_11274        subread: CLR
    // m140913_050931_42139_c1007..._s1_p0/9/0_8372   RS II subread
    // Returns 0 if this is not a PacBio name, 1 for a subread, 2 for a CCS read.
    int pacbio_name_kind(const string& token)
    {
      if ((token.size() < 8) || (token[0] != 'm') || !isdigit(static_cast<unsigned char>(token[1]))) return 0;
      size_t slash_1 = token.find('/');
      if (slash_1 == string::npos) return 0;
      if (token.find('_') > slash_1) return 0;
      size_t slash_2 = token.find('/', slash_1 + 1);
      if (slash_2 == string::npos) return 0;
      if (!all_digits(token, slash_1 + 1, slash_2)) return 0;
      size_t slash_3 = token.find('/', slash_2 + 1);
      string last = token.substr(slash_2 + 1, (slash_3 == string::npos) ? string::npos : slash_3 - slash_2 - 1);
      if (last == "ccs") return 2;
      size_t underscore = last.find('_');
      if ((underscore != string::npos) && all_digits(last, 0, underscore) && all_digits(last, underscore + 1)) return 1;
      return 0;
    }

    // HWI-ST1234:100:C0ABCACXX:6:1101:1234:2100      CASAVA 1.8+ (7 fields; 8 with a UMI)
    // HWUSI-EAS100R:6:73:941:1973#0/1                older (5 fields)
    // The last four numeric fields are lane:tile:x:y in both.
    bool is_illumina_name(const string& token)
    {
      vector<string> fields = split(token, ":");
      size_t first_numeric;
      if ((fields.size() == 7) || (fields.size() == 8)) first_numeric = 3;
      else if (fields.size() == 5) first_numeric = 1;
      else return false;
      if (fields[0].empty()) return false;
      for (size_t i = first_numeric; i < first_numeric + 4; i++) {
        string field = fields[i];
        if (i == first_numeric + 3) {
          // y carries the old-style "#index/mate" or "/mate" suffix
          size_t end = field.find_first_of("#/");
          if (end != string::npos) field.erase(end);
        }
        if (!all_digits(field)) return false;
      }
      return true;
    }

    // 8-4-4-4-12 hexadecimal digits: what every ONT basecaller names a read
    //
    // The SRA and ENA pass it on in several forms, all seen in real downloads: with a mate suffix
    // (<uuid>/1), with the hyphens turned into underscores, and as <uuid>;<uuid> for a duplex read.
    bool is_uuid(const string& token)
    {
      if (token.size() < 36) return false;
      for (size_t i = 0; i < 36; i++) {
        if ((i == 8) || (i == 13) || (i == 18) || (i == 23)) {
          if ((token[i] != '-') && (token[i] != '_')) return false;
        } else if (!isxdigit(static_cast<unsigned char>(token[i]))) return false;
      }
      if (token.size() == 36) return true;
      if ((token[36] == '/') && all_digits(token, 37)) return true;
      if (token[36] == ';') return is_uuid(token.substr(37));
      return false;
    }

  } // anonymous namespace

  cReadTypeDetector::HeaderEvidence cReadTypeDetector::classify_header(const string& _name, const string& _comment)
  {
    // The FASTQ reader splits the header at the first space, but tags written by samtools fastq -T
    // follow the name after a TAB: split there too.
    string name = _name;
    string comment = _comment;
    size_t tab_pos = name.find('\t');
    if (tab_pos != string::npos) {
      comment = name.substr(tab_pos + 1) + (comment.empty() ? "" : " " + comment);
      name.erase(tab_pos);
    }

    // A dorado (or late guppy) model name ends in the accuracy and version: ..._400bps_hac@v4.3.0.
    // It is in basecall_model_version_id=, or in an RG:Z: tag when the FASTQ came from a BAM.
    // A super accuracy model is reported as high accuracy: the two are one read type.
    if (comment.find("_sup@v") != string::npos) return kNanoporeHAC;
    if (comment.find("_hac@v") != string::npos) return kNanoporeHAC;
    if (comment.find("_fast@v") != string::npos) return kNanoporeFast;

    // Any other ONT header: the accuracy is left to the quality scores
    if ((comment.find("runid=") != string::npos) || (comment.find("basecall_model_version_id=") != string::npos)
        || (comment.find("ch:i:") != string::npos))
      return kNanopore;

    // The name, then each field after it (SRA: "@SRR1.1 <original name> length=100")
    string token = name;
    size_t start = 0;
    while (true) {
      if (!token.empty()) {
        int pacbio_kind = pacbio_name_kind(token);
        if (pacbio_kind == 2) return kPacBioHiFi;
        if (pacbio_kind == 1) return kPacBioCLR;
        if (is_illumina_name(token)) return kIllumina;
        if (is_uuid(token)) return kNanopore;
      }
      if (start >= comment.size()) break;
      size_t end = comment.find_first_of(" \t", start);
      if (end == string::npos) end = comment.size();
      token = comment.substr(start, end - start);
      start = end + 1;
    }

    return kNone;
  }

  void cReadTypeDetector::add(const string& name, const string& comment, const string& qualities, const size_t length)
  {
    m_num_reads++;
    m_header_votes[classify_header(name, comment)]++;
    m_read_lengths.push_back(static_cast<uint32_t>(length));
    for (size_t i = 0; i < qualities.size(); i++)
      m_quality_chr_counts[static_cast<uint8_t>(qualities[i])]++;
  }

  uint32_t cReadTypeDetector::read_length_percentile(const double fraction) const
  {
    if (m_read_lengths.empty()) return 0;
    sort(m_read_lengths.begin(), m_read_lengths.end());
    size_t index = static_cast<size_t>(fraction * static_cast<double>(m_read_lengths.size() - 1) + 0.5);
    return m_read_lengths[min(index, m_read_lengths.size() - 1)];
  }

  double cReadTypeDetector::mean_quality(const uint8_t quality_chr_offset) const
  {
    double total_error = 0;
    uint64_t num_bases = 0;
    for (uint32_t c = 0; c < 256; c++) {
      if (m_quality_chr_counts[c] == 0) continue;
      double quality = static_cast<double>(c) - static_cast<double>(quality_chr_offset);
      total_error += static_cast<double>(m_quality_chr_counts[c]) * pow(10.0, -quality / 10.0);
      num_bases += m_quality_chr_counts[c];
    }
    if ((num_bases == 0) || (total_error <= 0)) return 0;
    return -10.0 * log10(total_error / static_cast<double>(num_bases));
  }

  double cReadTypeDetector::fraction_bases_at_least_quality(const uint8_t quality_chr_offset, const uint32_t quality) const
  {
    uint64_t num_bases = 0, num_at_least = 0;
    for (uint32_t c = 0; c < 256; c++) {
      num_bases += m_quality_chr_counts[c];
      if (c >= static_cast<uint32_t>(quality_chr_offset) + quality) num_at_least += m_quality_chr_counts[c];
    }
    return (num_bases == 0) ? 0 : static_cast<double>(num_at_least) / static_cast<double>(num_bases);
  }

  uint32_t cReadTypeDetector::num_distinct_qualities() const
  {
    uint32_t n = 0;
    for (uint32_t c = 0; c < 256; c++) if (m_quality_chr_counts[c]) n++;
    return n;
  }

  void cReadTypeDetector::detect(const uint8_t quality_chr_offset, string& read_type, string& read_type_source) const
  {
    read_type = kReadTypeUnknown;
    read_type_source = "";
    if (m_num_reads == 0) return;

    // ---- Headers: trusted when most of the sampled reads agree ----
    HeaderEvidence best = kNone;
    for (uint32_t e = kNone + 1; e < kNumHeaderEvidence; e++)
      if (m_header_votes[e] > m_header_votes[best] || (best == kNone)) best = static_cast<HeaderEvidence>(e);
    const bool header_decides = (best != kNone) && (2 * m_header_votes[best] > m_num_reads);

    // An ONT header without a model name still settles the platform; so does everything below
    const bool nanopore_by_header = header_decides && (best == kNanopore);

    if (header_decides && !nanopore_by_header) {
      read_type_source = kReadTypeSourceHeader;
      switch (best) {
        case kIllumina:     read_type = kReadTypeIllumina; break;
        case kNanoporeFast: read_type = kReadTypeNanoporeFast; break;
        case kNanoporeHAC:  read_type = kReadTypeNanoporeHAC; break;
        case kPacBioCLR:    read_type = kReadTypePacBioCLR; break;
        case kPacBioHiFi:   read_type = kReadTypePacBioHiFi; break;
        default: break;
      }
      return;
    }

    // ---- Lengths and qualities ----
    read_type_source = kReadTypeSourceQuality;

    if (!nanopore_by_header) {

      // Short reads: nearly all under 700 bases, which clears merged 2x300 pairs, and none of
      // 1000 or more. The second condition is for nanopore libraries of short fragments, which
      // are mostly short reads but always have some long ones (SRR14995061: median 153 bases,
      // longest of the sample 6400), where a short-read platform has a hard maximum. A nanopore
      // file with no read that long and no usable names is still called illumina.
      if ( (read_length_percentile(0.99) < kReadTypeLongReadLength)
        && (read_length_percentile(1.0) < kReadTypeLongestShortRead) ) {
        read_type = kReadTypeIllumina;
        return;
      }

      // Every base one score: there are no real qualities. That is how PacBio subreads come
      // from Sequel instruments on (all 7 sampled Sequel / Sequel II CLR runs).
      const uint32_t num_distinct = num_distinct_qualities();
      if (num_distinct <= 1) {
        read_type = kReadTypePacBioCLR;
        return;
      }

      // PacBio RS II subreads do have qualities, but only Q0-Q15 (all 5 sampled RS II runs), where
      // even the worst nanopore sample seen (mean Q5) has 5% of its bases at Q20 or better.
      if (fraction_bases_at_least_quality(quality_chr_offset, kReadTypeCLRMaximumQuality) < kReadTypeCLRMaximumQualityFraction) {
        read_type = kReadTypePacBioCLR;
        return;
      }

      // HiFi consensus qualities run to the top of the scale (Q93), which no single-pass
      // basecaller reaches; newer instruments bin them into a handful of values instead, and a
      // nanopore file has dozens of distinct values.
      if ( (fraction_bases_at_least_quality(quality_chr_offset, kReadTypeHiFiQuality) >= kReadTypeHiFiQualityFraction)
        || (num_distinct <= kReadTypeBinnedQualityValues) ) {
        read_type = kReadTypePacBioHiFi;
        return;
      }
    }

    // Nanopore: the basecalling accuracy shows in the mean quality
    const double q = mean_quality(quality_chr_offset);
    read_type = (q < kReadTypeNanoporeHACMeanQuality) ? kReadTypeNanoporeFast : kReadTypeNanoporeHAC;
  }

  void warn_constant_base_quality(const string& file_name, const uint32_t original_quality, const uint32_t replacement_quality)
  {
    string original = (original_quality == kMissingBamBaseQuality)
                      ? "no base quality scores (QUAL is '*')"
                      : "the same base quality score (Q" + to_string(original_quality) + ")";
    WARN("Every base in read file " + file_name + "\nhas " + original + ". These scores carry no information:\n"
         "PacBio CLR reads downloaded from the SRA, for example, are all Q0. breseq will treat them as\n"
         "uninformative and give every base in this file quality score Q" + to_string(replacement_quality) + " instead.\n"
         "Its error rate is calibrated from the reads, as for any quality score.\n"
         "Set the score with --constant-base-quality-replacement, or keep the original scores\n"
         "with --no-constant-base-quality-replacement.");
  }

  // converts a sequence file
  void convert_fastq(const string &from_file_name, const string &to_file_name, const string &from_format, const string &to_format, bool _reverse_complement)
  {
    cFastqFile input_fastq_file(from_file_name.c_str(), ios::in);
    cFastqFile output_fastq_file(to_file_name.c_str(), ios::out);

    cFastqQualityConverter fqc(from_format, to_format);

    cFastqSequence on_sequence;
    while (input_fastq_file.read_sequence(on_sequence, fqc)) 
    {
      fqc.convert_sequence(on_sequence);
      if (_reverse_complement) 
        on_sequence = reverse_complement(on_sequence);
      output_fastq_file.write_sequence(on_sequence);
    }
    
  }
  
  void fastq_sequence_trim_end_on_base_quality(cFastqSequence& seq, const uint32_t base_quality)
  {
    for (uint32_t i=0; i<seq.m_qualities.size(); i++) {
      if (static_cast<uint8_t>(seq.m_qualities[i]-33) < base_quality) {
        seq.m_sequence.resize(i);
        seq.m_qualities.resize(i);
        break; 
      }
    }
  }
  
  bool cFastqSequence::identical(cFastqSequence& seq)
  {
    return ( (this->m_sequence == seq.m_sequence) && (this->m_qualities == seq.m_qualities) );
  }


  // constructor
  cFastqQualityConverter::cFastqQualityConverter(const string &_from_quality_format, const string &_to_quality_format)
  {
    // Set up maps between formats
    map<string,uint8_t> format_to_chr_offset;
    format_to_chr_offset["SANGER"] = 33;
    format_to_chr_offset["SOLEXA"] = 64;
    format_to_chr_offset["ILLUMINA_1.3+"] = 64;
    
    map<string,string> format_to_quality_type;
    format_to_quality_type["SANGER"] = "PHRED";
    format_to_quality_type["SOLEXA"] = "SOLEXA";
    format_to_quality_type["ILLUMINA_1.3+"] = "PHRED";
    
    from_quality_format = _from_quality_format;
    to_quality_format = _to_quality_format;
    
    // check what we asked for is valid...
    ASSERT(format_to_chr_offset.count(from_quality_format), 
           "Unknown FASTQ quality score format: " + from_quality_format + "\nValid choices are 'SANGER', 'SOLEXA', 'ILLUMINA_1.3+', 'NUMERICAL'");
    ASSERT(format_to_chr_offset.count(to_quality_format), 
           "Unknown FASTQ quality score format: " + to_quality_format + "\nValid choices are 'SANGER', 'SOLEXA', 'ILLUMINA_1.3+', 'NUMERICAL'");

    
    from_quality_type = format_to_quality_type[from_quality_format];
    to_quality_type = format_to_quality_type[to_quality_format];

    from_chr_offset = format_to_chr_offset[from_quality_format];
    to_chr_offset = format_to_chr_offset[to_quality_format];

    
    this->resize(256);
    for (uint16_t i = 0; i<=255; i++) {
      (*this)[i] = 0;
    }
    
    for (uint16_t from_chr = 0; from_chr<=255; from_chr++) {

      int32_t from_quality = from_chr - from_chr_offset;
      
      // Calculate the probability of error
      double probability_of_error;
      
      if (from_quality_type == "SOLEXA") {
        probability_of_error = 1 / (1+pow(10,(double)from_quality/10));
      } else if (from_quality_type == "PHRED") {
        probability_of_error = pow(10,-(double)from_quality/10);
      } else {
        cerr << "Unknown base quality score type: " << from_quality_type << endl;
        exit(-1);
      }
      
      //Convert back to quality score
      int32_t to_quality;
            
      if (to_quality_type == "SOLEXA") {
        to_quality = static_cast<uint32_t>(round(10 * log((1-probability_of_error)/probability_of_error) / log(10)));
      } else if (to_quality_type == "PHRED") {
        to_quality = static_cast<uint32_t>(round(-10 * log(probability_of_error) / log(10)));
      } else {
        cerr << "Unknown base quality score type: " << to_quality_type << endl;
        exit(-1);
      }
            
      int16_t to_chr = to_quality + to_chr_offset;
      
      // May be out of range
      if ((to_chr < 0) || (to_chr > 255)) continue;

      (*this)[(uint8_t)from_chr] = (uint8_t)to_chr;
      
      // Debug
      //cerr << from_chr << " => " << to_chr << endl;
    }     
    
  }

  void cFastqQualityConverter::convert_sequence(cFastqSequence &seq) {
    
    for(uint32_t i=0; i < seq.m_qualities.size(); i++)
    {
      seq.m_qualities[i] = (*this)[seq.m_qualities[i]];
    }
  }
  
  string cFastqQualityConverter::predict_fastq_file_format(const string& file_name, uint64_t& num_original_reads, uint64_t& num_original_bases, uint32_t& read_length_min, uint32_t& read_length_max, uint8_t& min_quality_score, uint8_t& max_quality_score)
  {
  // Initialize the input variables!
    num_original_reads = 0;
    num_original_bases = 0;
    read_length_min = numeric_limits<uint32_t>::max();
    read_length_max = 0;
    min_quality_score = 255;
    max_quality_score = 0;
    
  // Set up maps between formats
  map<string,uint8_t> format_to_chr_offset;
  format_to_chr_offset["SANGER"] = 33;
  format_to_chr_offset["SOLEXA"] = 64;
  format_to_chr_offset["ILLUMINA_1.3+"] = 64;
    
  cFastqFile input_fastq_file(file_name.c_str(), fstream::in);
  input_fastq_file.m_check_for_repeated_read_names = true;
  
  cFastqSequence on_sequence;
  cFastqQualityConverter prelim_fqc("ILLUMINA_1.3+", "SANGER");
  
  while (input_fastq_file.read_sequence(on_sequence, prelim_fqc)) {
    
    //increment read number
    num_original_reads++;
    
    //check sequence length
    read_length_min = min<uint32_t>(read_length_min, on_sequence.m_sequence.size());
    read_length_max = max<uint32_t>(read_length_max, on_sequence.m_sequence.size());

    
    //add current sequence length to number of bases
    num_original_bases += on_sequence.m_sequence.size();
      
      //iterate through sequence grabbing the associated scores
    for (uint32_t i=0; i<on_sequence.m_qualities.size(); i++) {
      int this_score(uint8_t(on_sequence.m_qualities[i]));
      if( this_score > max_quality_score ) max_quality_score = this_score;
        if( this_score < min_quality_score ) min_quality_score = this_score;
    }
  }
  
  // Default is SANGER
  string quality_format = "SANGER";
  
  // Typical range: (-5, 40) + 64
  if (min_quality_score >= format_to_chr_offset["SOLEXA"] - 5) {
    quality_format = "SOLEXA";
  } 
  // Typical range:  (0, 40) + 64
  if (min_quality_score >= format_to_chr_offset["ILLUMINA_1.3+"]) {
    quality_format = "ILLUMINA_1.3+";
  }
    
    return quality_format;
  }
  
  //constructor
  cFastqFile::cFastqFile()
    : flexgzfstream()
    , m_current_line(0)
    , m_file_name("")
    , m_check_for_repeated_read_names(false)
    , m_last_read_name(""), m_repeated_read_name_count(0)
  {
  }
 
  
  cFastqFile::cFastqFile(const string &file_name, std::ios_base::openmode mode, unsigned int num_threads)
    : flexgzfstream(file_name.c_str(), mode, num_threads)
    , m_current_line(0)
    , m_file_name(file_name)
    , m_check_for_repeated_read_names(false)
    , m_last_read_name("")
    , m_repeated_read_name_count(0)
  {
  }

  // read one sequence record from the file
  bool cFastqFile::read_sequence(cFastqSequence &sequence, cFastqQualityConverter& fqc) {
    
    // We're done, no error
    if (m_stream->eof())
     return false; 

    
    uint32_t count = 0;
    string line;
        
    memset(sequence.m_base_counts, 0, sizeof(sequence.m_base_counts));
    
    // get the next four lines
    while (count < 4) {
      breseq::getline(*m_stream, line);
      
      m_current_line++;
      
      // Didn't get a first line, then we ended correctly
      if (m_stream->eof()) {
        if (count == 0) {
          return false;
        } else {
          uint32_t last_valid_line = static_cast<uint32_t>(floor((m_current_line-1)/4.0) * 4);
          fprintf(stderr, "Incomplete FASTQ sequence record found at end of file.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line-1);
          fprintf(stderr, "You may be able to repair this damage and salvage the reads before this point with the command:\n");
          fprintf(stderr, "  head -n %u %s > new.fastq\n", last_valid_line, m_file_name.c_str());
          fprintf(stderr, "Then use \"new.fastq\" as input.\n");
          exit(-1);
        }
      }
      
      // Skip empty lines
      if (line.size()==0) continue;
      
      switch (count) {
        case 0:
          if( line[0] != '@' ) {
            fprintf(stderr, "FASTQ sequence record does not begin with @NAME line.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line);
            exit(-1);
          }
          sequence.m_name = line.substr(1,string::npos);
          
          // Delete any sequence name information after the first space...
          // Necessary for scrubbing SRA FASTQs, for example.
          { // block to keep inside this switch case
            // ...but keep it aside: it is where the basecaller model or the original read
            // name of an SRA download is, which cReadTypeDetector reads the read type from.
            size_t space_pos = sequence.m_name.find(" ");
            if (space_pos != string::npos)
            {
              sequence.m_name_comment.assign(sequence.m_name, space_pos + 1, string::npos);
              sequence.m_name.erase(space_pos);
            }
            else
            {
              sequence.m_name_comment.clear();
            }
          }
          
          // some SRA files have identical read names, we don't like this...
          if (m_check_for_repeated_read_names)
          {
            string original_read_name = sequence.m_name;
            if (m_last_read_name == sequence.m_name)
            {
              m_repeated_read_name_count++;
              sequence.m_name += "r" + to_string(m_repeated_read_name_count);
            }
            else
            {
              m_repeated_read_name_count = 0;
            }
            m_last_read_name = original_read_name;
          }
          
          break;
          
        case 1:
          sequence.m_sequence = line;
          
          for (uint32_t i=0; i<sequence.m_sequence.size(); i++) {
            
            // convert to uppercase and require
            // reformatting if this was necessary
            switch (sequence.m_sequence[i]) {
                
              case 'A':
              case 'T':
              case 'C':
              case 'G':
                break;
                
              case 'N':
                break;
                
              case 'a':
                sequence.m_sequence.replace(i,1,1,'A');
                break;
                
              case 't':
                sequence.m_sequence.replace(i,1,1,'T');
                break;
                
              case 'c':
                sequence.m_sequence.replace(i,1,1,'C');
                break;
                
              case 'g':
                sequence.m_sequence.replace(i,1,1,'G');
                break;

              case 'n':
                sequence.m_sequence.replace(i,1,1,'N');
                break;
              
              // all other characters converted to 'N'
              default :
                sequence.m_sequence.replace(i,1,1,'N');

            }

            // keep a count of the number of each base for detecting homopolymeric reads
            sequence.m_base_counts[basechar2index(sequence.m_sequence[i])]++;
            
            if(sequence.m_sequence[i] != 'A' && 
               sequence.m_sequence[i] != 'T' && 
               sequence.m_sequence[i] != 'G' && 
               sequence.m_sequence[i] != 'C' && 
               sequence.m_sequence[i] != 'N') {
              
              fprintf(stderr, "FASTQ sequence character not allowed %c.\nSequence: %s\nFile %s\nLine: %d\n", 
                      sequence.m_sequence[i], sequence.m_sequence.c_str(), m_file_name.c_str(), m_current_line);
              exit(-1);
            }
          }
          
          break;
        case 2:
          
          //Only need to see if the first character is a +
          if( line[0] != '+' ) {
            fprintf(stderr, "FASTQ sequence record does not contain +NAME line.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line);
            exit(-1);
          }
          // Could optionally check to see if the name after the + was either absent or identical to the earlier name
          sequence.m_name_plus = line.substr(1,string::npos);

          break;
        case 3:
          
          if (sequence.m_sequence.size() == line.size()) {
            sequence.m_qualities = line;
          } else if ((line.find_first_of(" ") != string::npos) && (line.find_first_not_of(" -0123456789\t") == string::npos)) {
            
            sequence.m_numerical_qualities = true;
            vector<string> numerical_qualities(split(line, " "));
            if( sequence.m_sequence.size() != numerical_qualities.size() ) {
              fprintf(stderr, "FASTQ sequence record has different SEQUENCE and numerical QUALITY lengths.\nFile %s\nLine: %d\n", m_file_name.c_str(), m_current_line);
              exit(-1);
            }
            
            // convert the qualities to characters with the Illumina offset (which keeps things from being negative)
            sequence.m_qualities = "";
            for(vector<string>::iterator it = numerical_qualities.begin(); it != numerical_qualities.end(); it++)
            {
              // use of uint16_t is on purpose to force proper conversion @JEB
              sequence.m_qualities += static_cast<char>(from_string<int16_t>(*it)) + fqc.from_chr_offset;
            }
          } else {
            ERROR("FASTQ QUALITY line length does not match SEQUENCE length.\nFile: " + m_file_name + " Line: " + to_string(m_current_line) + "\nSequence:     " + sequence.m_sequence + "\nQuality Line: " + line);
          }

          break;
      }
      
      count++;
  }
    
    return true;
  }

  void cFastqFile::write_sequence(const cFastqSequence &sequence) {
    *m_stream << "@" << sequence.m_name << endl;
    *m_stream << sequence.m_sequence << endl;
    *m_stream << "+" << sequence.m_name_plus << endl;
    *m_stream << sequence.m_qualities << endl;
  }


  // Default seed is the wall clock, so an unseeded run is NOT reproducible. That is deliberate --
  // it preserves the historical behavior -- but do_simulate_reads prints the resolved seed so any
  // run can be reproduced after the fact, and every test passes --seed explicitly.
  int32_t cSimFastqSequence::SEED_VALUE = time(NULL);

  map<char, string> cSimFastqSequence::random_snp_base_options =
  make_map<char, string>
    ('A', "TCG")
    ('T', "ACG")
    ('C', "ATG")
    ('G', "ATC")
    ('N', "ACTG");

  char cSimFastqSequence::random_insertion_base_options[] =
  {'A', 'C', 'T', 'G'};

  // 1e9 * 10^(-k/10) for k = 0..9. Every Phred score's error probability is one of these ten
  // mantissas divided by a power of ten, so the whole scale is exact integer arithmetic. This
  // replaces a pow(10, -q/10.0) call: libm is not required to be correctly rounded, and a last-ulp
  // difference between glibc and Apple libm would change which bases become errors and so change
  // the simulated reads between platforms. See portable_random.h.
  static const uint32_t kPhredMantissaPerBillion[10] = {
    1000000000u, 794328235u, 630957344u, 501187234u, 398107171u,
     316227766u, 251188643u, 199526231u, 158489319u, 125892541u
  };

  static uint32_t phred_error_rate_per_billion(uint32_t phred)
  {
    uint32_t decades = phred / 10;
    uint32_t rate = kPhredMantissaPerBillion[phred % 10];
    for (uint32_t i = 0; i < decades; i++) rate /= 10;
    return rate;
  }

  bool cSimFastqSequence::is_error_base(cPortableRandom& rng, uint32_t phred)
  {
    return rng.next_event_per_billion(phred_error_rate_per_billion(phred));
  }

  char cSimFastqSequence::get_random_error_base(cPortableRandom& rng, const char not_this_base)
  {
    ASSERT(random_snp_base_options.count(not_this_base), "Unexpected base: " + to_string(not_this_base));

    const string& options = random_snp_base_options[not_this_base];
    return options[rng.next_below(static_cast<uint32_t>(options.size()))];
  }

  char cSimFastqSequence::get_random_insertion_base(cPortableRandom& rng)
  {
    return base_char_list[rng.next_below(4)];
  }

  char cSimQualityModel::draw(cPortableRandom& rng, uint32_t cycle, uint32_t read_size) const
  {
    // Real Illumina quality is roughly flat over most of the read and falls off fastest over the
    // last few cycles, so the decline is quadratic in the cycle rather than linear. Integer
    // throughout: the division truncates identically everywhere.
    uint32_t n = (read_size > 1) ? (read_size - 1) : 1;
    uint32_t drop = (phred_start > phred_end) ? (phred_start - phred_end) : 0;
    uint32_t shift = static_cast<uint32_t>(
      (static_cast<uint64_t>(drop) * cycle * cycle) / (static_cast<uint64_t>(n) * n));

    int32_t center = static_cast<int32_t>(phred_start) - static_cast<int32_t>(shift);
    int32_t q = rng.next_gaussian(center, static_cast<int32_t>(phred_stdev));

    if (q < static_cast<int32_t>(phred_min)) q = static_cast<int32_t>(phred_min);
    if (q > static_cast<int32_t>(phred_max)) q = static_cast<int32_t>(phred_max);

    return static_cast<char>(q + 33);   // Sanger offset
  }

  void cSimFastqSequence::warn_if_quality_format_ambiguous(char min_quality_char)
  {
    // normalize_fastq classifies a file from the minimum quality character in it: >= 59 is taken to
    // be SOLEXA and every score is shifted by 31. A library with no low-quality tail is therefore
    // silently reinterpreted, which is very hard to diagnose downstream.
    if (min_quality_char >= 59) {
      WARN("The lowest quality score written was '" + to_string(min_quality_char) +
           "' (Phred " + to_string(static_cast<int32_t>(min_quality_char) - 33) + "). breseq infers a "
           "FASTQ's quality encoding from the lowest score in the file and will misread this one as "
           "SOLEXA rather than SANGER. Lower --quality-end or --quality-min so that the file "
           "contains at least one score below Phred 26.");
    }
  }

  cFastqSequence cSimFastqSequence::simulate_read(const cAnnotatedSequence& ref_sequence,
                                                  int32_t  left_1,
                                                  uint32_t read_size,
                                                  int8_t   strand,
                                                  const cSimQualityModel& quality_model,
                                                  const cSimErrorModel&   error_model,
                                                  cPortableRandom& rng,
                                                  const string& read_name,
                                                  char& min_quality_char,
                                                  bool verbose)
  {
    cFastqSequence ret_val;
    ret_val.m_name = read_name;
    ret_val.m_name_plus = "[strand]:" + to_string(static_cast<int32_t>(strand)) +
                          "\t[left_1]:" + to_string(left_1);

    /*! Algorithm:
        1) Fetch the read's reference footprint [left_1, left_1+read_size-1], plus a tail of up to
           read_size more bases in the direction the read is built, so that simulated deletions have
           somewhere to draw replacement bases from. Reverse-complement the whole buffer for the
           minus strand, which puts the read's 5' end first.
        2) Walk the buffer, deciding per base whether it is a deletion, an insertion, a substitution
           error, or a normal base.
        3) Give every emitted base a quality score from the per-cycle model, and derive the
           substitution probability from that score.
     */

    // Buffer origin: the plus-strand read starts at its footprint; the minus-strand read's buffer
    // must START read_size earlier, because after reverse-complementing, the footprint ends up at
    // the FRONT. Fetching an equal tail in each case keeps the deletion slack symmetric.
    int32_t seq_len = static_cast<int32_t>(ref_sequence.get_sequence_length());
    int32_t fetch_start = (strand == -1) ? (left_1 - static_cast<int32_t>(read_size)) : left_1;

    // On a circular sequence get_sequence_1_start_size wraps; on a linear one it silently CLAMPS and
    // returns a shorter string. The old code indexed the returned buffer blindly and so could read
    // past the end of it for any read near a linear sequence's end -- both undefined behavior and a
    // reproducibility hazard, since the bytes past the end are allocator-dependent. Clamp here
    // instead, and let the loop below cope with a short buffer.
    string ref_segment;
    if (ref_sequence.is_circular()) {
      ref_segment = ref_sequence.get_sequence_1_start_size(fetch_start, 2 * read_size);
    } else {
      int32_t lo = max<int32_t>(1, fetch_start);
      int32_t hi = min<int32_t>(seq_len, fetch_start + 2 * static_cast<int32_t>(read_size) - 1);
      if (lo <= hi) ref_segment = ref_sequence.get_sequence_1(lo, hi);
    }

    if (strand == -1) {
      ref_segment = reverse_complement(ref_segment);
    }

    ASSERT(ref_segment.size() >= read_size,
           "Simulated read at " + to_string(left_1) + " would extend past the end of sequence " +
           ref_sequence.m_seq_id + ". The caller must place reads so their whole footprint fits.");

    string verbose_errors(read_size, ' ');
    string verbose_deletions(ref_segment.size(), ' ');
    string verbose_insertions(ref_segment.size(), ' ');

    ret_val.m_sequence.resize(read_size);
    ret_val.m_qualities.resize(read_size);

    size_t index_to_assign = 0;       // position in the simulated read
    size_t index_in_ref_segment = 0;  // position in the reference buffer

    while (index_to_assign < read_size) {

      // A deletion consumes a reference base without emitting one, so it is only possible while
      // there is still buffer left to consume. Near the end of a linear sequence there may not be.
      bool can_delete = (index_in_ref_segment + 1 < ref_segment.size());

      if (can_delete && rng.next_event_per_billion(error_model.deletion_rate_per_million * 1000)) {
        verbose_deletions[index_in_ref_segment] = ref_segment[index_in_ref_segment];
        ++index_in_ref_segment;
        continue;
      }

      char quality_char = quality_model.draw(rng, static_cast<uint32_t>(index_to_assign), read_size);
      if (quality_char < min_quality_char) min_quality_char = quality_char;
      uint32_t phred = static_cast<uint32_t>(quality_char) - 33;

      if (rng.next_event_per_billion(error_model.insertion_rate_per_million * 1000)) {
        // An insertion emits a base without consuming one from the reference.
        char base_to_insert = cSimFastqSequence::get_random_insertion_base(rng);
        ret_val.m_sequence[index_to_assign] = base_to_insert;
        ret_val.m_qualities[index_to_assign] = quality_char;
        verbose_insertions[index_in_ref_segment] = base_to_insert;
        ++index_to_assign;
        continue;
      }

      char ref_base = ref_segment[index_in_ref_segment];

      if (cSimFastqSequence::is_error_base(rng, phred)) {
        char new_base = cSimFastqSequence::get_random_error_base(rng, ref_base);
        ret_val.m_sequence[index_to_assign] = new_base;
        verbose_errors[index_to_assign] = new_base;
      } else {
        ret_val.m_sequence[index_to_assign] = ref_base;
      }
      ret_val.m_qualities[index_to_assign] = quality_char;

      ++index_to_assign;
      ++index_in_ref_segment;
    }

    if (verbose) {
      if (verbose_deletions.find_first_not_of(' ')  != string::npos ||
          verbose_insertions.find_first_not_of(' ') != string::npos ||
          verbose_errors.find_first_not_of(' ')     != string::npos) {
        printf("\tVerbose output for simulated read    :  %s\n", ret_val.m_name.c_str());
        printf("\tReference segment                    :  %s\n", ref_segment.c_str());
        if (verbose_errors.find_first_not_of(' ')     != string::npos)
          printf("\tSimulated errors                     :  %s\n", verbose_errors.c_str());
        if (verbose_deletions.find_first_not_of(' ')  != string::npos)
          printf("\tSimulated DEL                        :  %s\n", verbose_deletions.c_str());
        if (verbose_insertions.find_first_not_of(' ') != string::npos)
          printf("\tSimulated INS                        :  %s\n", verbose_insertions.c_str());
        printf("\tFinal simulated read sequence        :  %s\n", ret_val.m_sequence.c_str());
        printf("\tFinal simulated read quality scores  :  %s\n", ret_val.m_qualities.c_str());
        printf("\n");
      }
    }

    return ret_val;
  }

  void cSimFastqSequence::simulate_single_ends(const cAnnotatedSequence& sequence,
                                               uint32_t n_reads,
                                               uint32_t read_size,
                                               uint32_t seed,
                                               const cSimQualityModel& quality_model,
                                               const cSimErrorModel& error_model,
                                               string file_name,
                                               bool verbose)
  {
    cFastqFile out(file_name.c_str(), ios_base::out);

    int32_t seq_len = static_cast<int32_t>(sequence.get_sequence_length());
    char min_quality_char = 126;

    for (uint32_t i = 0; i < n_reads; ++i) {
      // One independent stream per read, so a read depends only on its own index. Adding a draw
      // elsewhere then cannot shift every subsequent read.
      cPortableRandom rng(seed, i);

      int8_t strand = rng.next_coin() ? +1 : -1;

      // Place the whole read inside the sequence. The old code drew a start anywhere in the
      // sequence and let the fetch clamp, which produced short reads (and out-of-bounds indexing).
      int32_t left_1 = sequence.is_circular()
                         ? static_cast<int32_t>(rng.next_in_range_1(1, seq_len))
                         : static_cast<int32_t>(rng.next_in_range_1(1, seq_len - read_size + 1));

      cFastqSequence read = simulate_read(sequence, left_1, read_size, strand,
                                          quality_model, error_model, rng,
                                          "READ-" + to_string(i + 1), min_quality_char, verbose);
      out.write_sequence(read);

      if (verbose && ((i + 1) % 10000 == 0 || (i + 1) == n_reads)) {
        ostringstream progress_message;
        progress_message << "\tREAD: " << setw(12) << right << (i + 1) << "/" << n_reads;
        print_progress_line(progress_message.str());
      }
    }

    warn_if_quality_format_ambiguous(min_quality_char);
  }

  void cSimFastqSequence::simulate_paired_ends(const cAnnotatedSequence& sequence,
                                               uint32_t n_pairs,
                                               uint32_t read_size,
                                               uint32_t fragment_mean,
                                               uint32_t fragment_stdev,
                                               uint32_t seed,
                                               const cSimQualityModel& quality_model,
                                               const cSimErrorModel& error_model,
                                               string pair_1_file_name,
                                               string pair_2_file_name,
                                               bool verbose)
  {
    cFastqFile pair_1_out(pair_1_file_name.c_str(), ios_base::out);
    cFastqFile pair_2_out(pair_2_file_name.c_str(), ios_base::out);

    int32_t seq_len = static_cast<int32_t>(sequence.get_sequence_length());
    char min_quality_char = 126;

    for (uint32_t i = 0; i < n_pairs; ++i) {
      cPortableRandom rng(seed, i);

      // Fragment length. Clamping at read_size (rather than at 2*read_size) is what allows the mates
      // to OVERLAP, which is a perfectly normal library -- tests/data/tmv_plasmid is 190 bp
      // fragments read 2 x 151. It also removes the unsigned underflow the old code had when the
      // drawn length came out below 2*read_size.
      int32_t L = rng.next_gaussian(static_cast<int32_t>(fragment_mean),
                                    static_cast<int32_t>(fragment_stdev));
      L = max<int32_t>(L, static_cast<int32_t>(read_size));
      L = min<int32_t>(L, seq_len);

      // Leftmost base of the fragment. On a linear sequence the whole fragment must fit inside, just
      // as a library prepared from linear DNA would; a circular one may wrap.
      int32_t left = sequence.is_circular()
                       ? static_cast<int32_t>(rng.next_in_range_1(1, seq_len))
                       : static_cast<int32_t>(rng.next_in_range_1(1, seq_len - L + 1));

      // FR ("innie"): each mate reads inward from one end of the fragment. Which physical end is
      // read 1 is a coin flip, so both strands are equally represented in R1 -- previously every R1
      // was a plus-strand read. This does NOT change the pair's orientation label, which breseq
      // derives from 5'-end order rather than R1/R2 identity (see
      // best_pair_orientation_and_distance in candidate_junctions.cpp).
      bool r1_is_left = rng.next_coin();
      int32_t left_read_left_1  = left;
      int32_t right_read_left_1 = left + L - static_cast<int32_t>(read_size);

      string name = "READ-" + to_string(i + 1);

      cFastqSequence r1 = simulate_read(sequence,
                                        r1_is_left ? left_read_left_1 : right_read_left_1,
                                        read_size, r1_is_left ? +1 : -1,
                                        quality_model, error_model, rng,
                                        name + "/1", min_quality_char, verbose);
      cFastqSequence r2 = simulate_read(sequence,
                                        r1_is_left ? right_read_left_1 : left_read_left_1,
                                        read_size, r1_is_left ? -1 : +1,
                                        quality_model, error_model, rng,
                                        name + "/2", min_quality_char, verbose);

      // ALWAYS in this order: breseq pairs the two files positionally, reading them in lockstep.
      pair_1_out.write_sequence(r1);
      pair_2_out.write_sequence(r2);

      if (verbose && ((i + 1) % 10000 == 0 || (i + 1) == n_pairs)) {
        ostringstream progress_message;
        progress_message << "\tPAIR: " << setw(12) << right << (i + 1) << "/" << n_pairs;
        print_progress_line(progress_message.str());
      }
    }

    warn_if_quality_format_ambiguous(min_quality_char);
  }

// Simulates perfectly tiled reads - assuming a circular genome
//
// If coverage is 2 x read length (or greater), then every possible read is simulated
// on top and bottom strands.
//
// if it is less, then every ceiling(2 x length / coverage) bases, reads on both strands
// are simulated on each strand for every start position.
//
// This assumes a circular genome
void cSimFastqSequence::simulate_tiled(const cAnnotatedSequence& sequence,
                                             uint32_t read_size,
                                             uint32_t coverage,
                                             const cSimQualityModel& quality_model,
                                             string file_name,
                                             bool verbose)
{
  (void) verbose;

  // Tiled reads are deterministic by construction -- fixed spacing, no RNG, no simulated errors --
  // so this mode reproduces regardless of the seed. Qualities come from the shared model rather than
  // a hard-coded constant, but with the flat settings this mode uses (start == end, stdev 0) the
  // model consumes no randomness and emits the same character for every base, as before.
  cPortableRandom unused_rng(0, 0);

  cFastqFile out(file_name.c_str(), ios_base::out);
  vector<int8_t>strands = make_vector<int8_t>(-1)(+1);
  
  uint32_t spacing = ceil( (2.0 * read_size) / coverage );
  spacing = max<uint32_t>(spacing, 0);
  cout << "Read length : " << read_size << endl;
  cout << "Coverage    : " << coverage << endl;
  cout << "Spacing     : " << spacing << endl;
  
  uint32_t read_index = 1;
  for (uint32_t start_1 = 1; start_1 <= sequence.get_sequence_length(); start_1+=spacing) {
    for (vector<int8_t>::iterator strand_it = strands.begin(); strand_it != strands.end(); strand_it++) {
      int8_t strand = *strand_it;
      
      cFastqSequence read;
      sprintf(read.m_name, "READ-%i", read_index++);
      sprintf(read.m_name_plus, "[strand]:%i\t[start_1]:%u", strand, start_1);
      
      read.m_sequence = sequence.get_sequence_1_start_size(start_1, read_size);
      if (strand == -1) {
        read.m_sequence = reverse_complement(read.m_sequence);
      }
      
      read.m_qualities.resize(read_size);
      for (uint32_t c = 0; c < read_size; c++)
        read.m_qualities[c] = quality_model.draw(unused_rng, c, read_size);

      out.write_sequence(read);
    }
  }

return;
}

// Reverse complement and also uppercase
// Convert most characters to 'N'. Might want to give errors on non-printable characters
  char reverse_complement_lookup_table[256] = {
/*  0*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 16*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 32*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 48*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 64*/    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 80*/    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/* 96*/    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*112*/    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*128*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*144*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*160*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*176*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*192*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*208*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*224*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
/*240*/    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
  };
  
} // breseq namespace

