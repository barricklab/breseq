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

#include "soft_clipping.h"

#include "alignment.h"
#include "identify_mutations.h"
#include "reference_sequence.h"
#include "output.h"
#include "settings.h"
#include "summary.h"

using namespace std;

namespace breseq {

// Clipped tails are packed 3 bits per base into a uint64_t, so at most 21 bases
// of a tail participate in the consensus. That is well beyond the default clip
// threshold (12) and keeps memory at 8 bytes per clip event rather than ~32 for
// a std::string -- on a 4 Mb genome at 150x there are ~200,000 clip events.
const uint32_t kSoftClippingMaxConsensusBases = 21;

// 21 bases x 3 bits fills bits 0-62, leaving the top bit free. The read's strand rides
// there so that per-strand agree counts cost nothing beyond the tail itself.
const uint64_t kSoftClippingTailStrandBit = (1ULL << 63);

// Base <-> 3-bit code. 4 means "not A/C/G/T"; it never wins a consensus column
// and always counts as a mismatch.
static inline uint32_t base_to_code(char b)
{
  switch (b) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
  }
}

static inline char code_to_base(uint32_t c)
{
  static const char* bases = "ACGT";
  return (c < 4) ? bases[c] : 'N';
}

// column <= 20, so the shift is at most 60 and the 3-bit mask reads bits 60-62 --
// kSoftClippingTailStrandBit (bit 63) is never part of any base code.
static inline uint32_t tail_code_at(uint64_t packed, uint32_t column)
{
  return static_cast<uint32_t>((packed >> (3 * column)) & 0x7ULL);
}

static inline bool tail_is_reversed(uint64_t packed)
{
  return (packed & kSoftClippingTailStrandBit) != 0;
}

/*
 * Per-column consensus over the clipped tails at one (seq_id, position, direction), and the
 * number of reads agreeing with it.
 *
 * Columns are indexed by distance from the clip wall, so reads align column-wise regardless
 * of direction and regardless of how much each read was clipped. A column whose bases are all
 * non-ACGT becomes 'N' and is dropped from both sides of the match count. Ties are broken by
 * lowest base code (A < C < G < T) so results are deterministic and goldens reproduce exactly.
 *
 * A read agrees if it matches the consensus in at least base_fraction of the informative
 * columns (rounded down). Reads clipped for uninteresting reasons -- adapter read-through,
 * low-quality tails -- disagree with each other, while reads spanning a single real breakpoint
 * all carry the same donor sequence.
 *
 * Note what this deliberately does NOT report: a reference copy of a mobile element receives
 * reads from every copy in the sample, so its boundary carries several distinct donor
 * sequences at once (ADP1 position 2151230 has 290 clipped reads in seven near-equal groups of
 * 47/38/37/36/36/36/21). The per-column plurality over such a mixture is a chimera matching no
 * actual read, so agree_count collapses to ~0 and the position is discarded. That is intended:
 * repeat-element boundaries are not the target of this evidence type, and junction (JC)
 * evidence reports them properly. An earlier version clustered the tails and used the largest
 * group specifically to keep these loci; it was removed because they are not wanted. On every
 * genuine single-donor breakpoint the two approaches agree exactly (one donor => one group =>
 * the group consensus is the per-column consensus), and on noisy positions the simpler rule
 * here is the more conservative one.
 *
 * The agreeing reads are also counted separately by the strand of the read they came from
 * (agree_count_out == agree_forward_out + agree_reverse_out, always). This is what the
 * strand test in add_sc_evidence() runs on, and it is the discriminator that matters most on
 * real Illumina data: the dominant false positive is a dark-cycle poly-G tail, which is always
 * the 3' end of the read, so it produces direction +1 clips only from forward-strand reads and
 * direction -1 clips only from reverse-strand reads. Measured over 29 LTEE clones, 993 of 1040
 * accepted SC calls had every clipped read on one strand, while every plausible real breakpoint
 * was strand-balanced. Note that a poly-G tail is stored reference-forward, so it reads as
 * poly-C for direction -1 -- the tails agree with each other perfectly and the consensus test
 * above cannot see anything wrong with them.
 *
 * base_fraction == 0 turns the whole test off: every clipped read counts and no sequence
 * is reported. consensus_out is returned in column order; the caller reverses it for
 * direction -1 so the stored sequence is always reference-forward.
 */
static void compute_clipped_tail_consensus(
                                           const vector<uint64_t>& tails,
                                           uint32_t K,
                                           double base_fraction,
                                           string& consensus_out,
                                           uint32_t& agree_count_out,
                                           uint32_t& agree_forward_out,
                                           uint32_t& agree_reverse_out
                                           )
{
  consensus_out.clear();
  agree_count_out = static_cast<uint32_t>(tails.size());
  agree_forward_out = 0;
  agree_reverse_out = 0;
  for (vector<uint64_t>::const_iterator it = tails.begin(); it != tails.end(); it++) {
    if (tail_is_reversed(*it)) agree_reverse_out++; else agree_forward_out++;
  }
  if (tails.empty()) return;

  // With the consensus test off every clipped read counts, and the strand split above is
  // already the split of all of them.
  if (base_fraction <= 0.0) return;

  // Per-column plurality consensus over every clipped tail at this position.
  consensus_out.resize(K, 'N');
  vector<uint32_t> consensus_codes(K, 4);
  uint32_t informative_columns = 0;

  for (uint32_t col = 0; col < K; col++) {
    uint32_t tally[4] = {0, 0, 0, 0};
    for (vector<uint64_t>::const_iterator it = tails.begin(); it != tails.end(); it++) {
      uint32_t c = tail_code_at(*it, col);
      if (c < 4) tally[c]++;
    }
    uint32_t best = 4, best_count = 0;
    for (uint32_t c = 0; c < 4; c++) {
      if (tally[c] > best_count) { best_count = tally[c]; best = c; }  // strict > keeps A<C<G<T tie-break
    }
    consensus_codes[col] = best;
    consensus_out[col] = code_to_base(best);
    if (best < 4) informative_columns++;
  }

  if (informative_columns == 0) {
    consensus_out.clear();
    agree_count_out = 0;
    agree_forward_out = 0;
    agree_reverse_out = 0;
    return;
  }

  // Rounded DOWN, so the default 0.95 over 12 compared bases requires 11 matches,
  // i.e. tolerates exactly one mismatched base.
  uint32_t required_matches =
    static_cast<uint32_t>(floor(base_fraction * static_cast<double>(informative_columns)));
  if (required_matches < 1) required_matches = 1;

  agree_count_out = 0;
  agree_forward_out = 0;
  agree_reverse_out = 0;
  for (vector<uint64_t>::const_iterator it = tails.begin(); it != tails.end(); it++) {
    uint32_t matches = 0;
    for (uint32_t col = 0; col < K; col++) {
      if (consensus_codes[col] >= 4) continue;      // uninformative column
      if (tail_code_at(*it, col) == consensus_codes[col]) matches++;
    }
    if (matches >= required_matches) {
      agree_count_out++;
      if (tail_is_reversed(*it)) agree_reverse_out++; else agree_forward_out++;
    }
  }
}

// A position is only testable where a spanning read is geometrically possible. A read must have
// minimum_clipped_bases of aligned reference on BOTH sides of p to count toward the denominator
// (see tabulate_soft_clipping_counts below), so the first and last minimum_clipped_bases positions
// of a contig can never have one -- a clip event there would be scored against nothing and report
// frequency 1.000. Position 1 and seq_length are always excluded (reads always soft-clip at
// reference edges), which is what this reduces to when minimum_clipped_bases is 0 or 1.
//
// This predicate MUST be applied identically to the clip events (numerator) and to the spanning
// coverage (denominator), or the identity
//     total_opportunities == sum over (position, direction) of n_i
// that the dispersion estimate depends on stops holding.
static inline bool sc_position_is_testable(uint32_t p, uint32_t seq_length, uint32_t minimum_clipped_bases)
{
  if (seq_length == 0) return false;
  uint32_t margin = (minimum_clipped_bases > 1) ? minimum_clipped_bases : 1;
  if (p <= margin) return false;
  if (static_cast<uint64_t>(p) + margin > seq_length) return false;
  return true;
}

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

  // Packed clipped tails, keyed identically to clipped_neg/clipped_pos. Columns run
  // outward from the clip wall in both cases (see the pushes below), so reads from
  // either direction and either strand are directly column-comparable. BAM stores SEQ
  // reference-forward, so no reverse-complementing is needed.
  const uint32_t consensus_bases = min(minimum_clipped_bases, kSoftClippingMaxConsensusBases);
  map<string, map<uint32_t, vector<uint64_t> > > tails_neg;
  map<string, map<uint32_t, vector<uint64_t> > > tails_pos;

  // Per-sequence "spanning" coverage, accumulated as a difference array. A read counts toward the
  // denominator at position p only if it reads through p with at least minimum_clipped_bases of
  // ALIGNED reference on BOTH sides:
  //
  //     reference_start_1 + min_bases <= p <= reference_end_1 - min_bases
  //
  // Because both sides are required the condition does not depend on clip direction, so one array
  // serves both tests.
  //
  // Requiring only the anchor side (min_bases below p for direction -1, as this once did) lets a
  // read whose alignment stops a few bases past p vote that the reference is contiguous through p
  // -- but a read terminating right next to p is far likelier to be evidence OF the breakpoint
  // being tested. A mobile-element insertion makes this quantitative: it creates two junctions a
  // target-site duplication apart, so reads supporting one junction end 1-9 bases from the other
  // and were being counted against it. Measured on ADP1 (Himar1 insertions), position 2968150
  // direction -1 had 71 clip events against a 64-read denominator in which NOT ONE read extended
  // 12 aligned bases past the position; 44 of the 64 carried their own qualifying clip, 42 of
  // those one base away at 2968151. Every fixed insertion measured came out near 0.53 (2968150:
  // 0.526, 2311060: 0.553, 600249: 0.542) -- the paired junction roughly doubles the denominator
  // and halves the frequency. Requiring both sides puts all three at 1.000, and costs background
  // positions only ~10% of their denominator (186 -> 164, 218 -> 198); p0 rises by the same ~10%,
  // so the test stays calibrated.
  //
  // The test is deliberately geometric rather than "does this read carry its own clip near p":
  // 20 of those 64 contaminating reads had no qualifying clip at all, because the aligner extended
  // a few mismatched bases into the donor sequence or clipped fewer than min_bases. Alignment
  // geometry is robust to that choice; the clip annotation is not.
  //
  // NOTE: a read with a large internal deletion (CIGAR D/N) is credited with spanning positions it
  // does not actually align to. Pre-existing behavior of the difference-array approach, but this
  // is now the only denominator, so it matters more.
  //
  // Kept split by the strand of the read, because that is the reference the strand test compares
  // the clipped reads against: coverage itself is not always 50/50, and a locally strand-skewed
  // pileup must not be read as evidence of a strand-skewed clip. The two arrays sum to what a
  // single array held before, so every genome-wide total below is unchanged.
  map<string, vector<int32_t> > spanning_diff_fw;
  map<string, vector<int32_t> > spanning_diff_rv;

  uint64_t total_clipped_read_ends = 0;  // total clip events, both directions (a read clipped at
                                         // both ends counts twice)
  uint64_t total_agreeing_clipped_read_ends = 0;  // subset whose tail matches its position consensus
  uint64_t total_strand_pure_agreeing_clipped_read_ends = 0;  // ...of those, the ones at positions
                                         // where every agreeing clipped read was on one strand
  uint64_t total_tested_positions = 0;   // N: (position, direction) pairs with n_i > 0

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

        // Accumulate spanning coverage (difference array). The read contributes only over the
        // stretch where it has minimum_clipped_bases of aligned reference on BOTH sides, which
        // needs an aligned span of at least 2*min_bases + 1 reference bases. The guard is
        // evaluated in 64 bits so it cannot wrap; inside it ref_end - min_bases >= ref_start +
        // min_bases >= 1, so neither index underflows, and the high index is at most
        // seq_length + 1, which is in bounds for a seq_length + 2 array.
        const bool read_reversed = a.reversed();
        const uint64_t strand_bit = read_reversed ? kSoftClippingTailStrandBit : 0ULL;

        vector<int32_t>& sd_fw = spanning_diff_fw[seq_id];
        vector<int32_t>& sd_rv = spanning_diff_rv[seq_id];
        if (sd_fw.empty()) sd_fw.resize(seq_length + 2, 0);
        if (sd_rv.empty()) sd_rv.resize(seq_length + 2, 0);
        vector<int32_t>& sd = read_reversed ? sd_rv : sd_fw;
        if (static_cast<uint64_t>(ref_end) >=
            static_cast<uint64_t>(ref_start) + 2ULL * minimum_clipped_bases) {
          sd[ref_start + minimum_clipped_bases] += 1;
          sd[ref_end - minimum_clipped_bases + 1] -= 1;
        }

        // Direction -1: leading soft-clip at reference_start_1
        if (query_begin_soft_clipping >= minimum_clipped_bases) {
          if (sc_position_is_testable(ref_start, seq_length, minimum_clipped_bases)) {
            clipped_neg[seq_id][ref_start]++;
            total_clipped_read_ends++;

            // Column 0 is the clipped base adjacent to the wall (just before
            // query_start_1); columns then run backwards, away from the reference.
            // Bit 63 carries the read's strand; see kSoftClippingTailStrandBit.
            uint64_t packed = strand_bit;
            for (uint32_t col = 0; col < consensus_bases; col++) {
              uint32_t q = a.query_start_1() - 1 - col;   // >= 1, guaranteed by the clip length test
              packed |= static_cast<uint64_t>(base_to_code(a.read_base_char_1(q))) << (3 * col);
            }
            tails_neg[seq_id][ref_start].push_back(packed);
          }
        }

        // Direction +1: trailing soft-clip at reference_end_1
        if (query_end_soft_clipping >= minimum_clipped_bases) {
          if (sc_position_is_testable(ref_end, seq_length, minimum_clipped_bases)) {
            clipped_pos[seq_id][ref_end]++;
            total_clipped_read_ends++;

            // Column 0 is the clipped base adjacent to the wall (just after query_end_1);
            // columns run forward, already reference-forward.
            // Bit 63 carries the read's strand; see kSoftClippingTailStrandBit.
            uint64_t packed = strand_bit;
            for (uint32_t col = 0; col < consensus_bases; col++) {
              uint32_t q = a.query_end_1() + 1 + col;     // <= read_length, guaranteed by the clip length test
              packed |= static_cast<uint64_t>(base_to_code(a.read_base_char_1(q))) << (3 * col);
            }
            tails_pos[seq_id][ref_end].push_back(packed);
          }
        }

      } // end for each alignment in group
    } // end while read_alignments
  } // end for each bam file

  // Prefix-sum the difference array into per-position spanning coverage, and total the
  // read-through opportunities over testable positions. The same spanning count is the
  // denominator of BOTH direction tests at a position, so each position contributes it twice.
  // Keeping the factor of two here is what preserves
  //
  //     total_opportunities == sum over (position, direction) of n_i
  //
  // which the Pearson chi-square expansion below depends on, and it keeps
  // summary.total_spanning_read_bases in the same units as previous versions.
  map<string, vector<uint32_t> > spanning_fw;
  map<string, vector<uint32_t> > spanning_rv;
  uint64_t total_spanning_read_bases = 0;
  uint64_t total_spanning_read_bases_fw = 0;
  uint64_t total_spanning_read_bases_rv = 0;
  for (map<string, uint32_t>::iterator seq_it = seq_lengths.begin(); seq_it != seq_lengths.end(); seq_it++) {
    const string& seq_id = seq_it->first;
    uint32_t seq_length = seq_it->second;
    if (spanning_diff_fw.count(seq_id) == 0) continue; // no reads on this sequence

    vector<int32_t>& sd_fw = spanning_diff_fw[seq_id];
    vector<int32_t>& sd_rv = spanning_diff_rv[seq_id];
    vector<uint32_t>& sp_fw = spanning_fw[seq_id];
    vector<uint32_t>& sp_rv = spanning_rv[seq_id];
    sp_fw.resize(sd_fw.size(), 0);
    sp_rv.resize(sd_rv.size(), 0);

    int32_t running_fw = 0;
    int32_t running_rv = 0;
    for (uint32_t p = 1; p < sd_fw.size(); p++) {
      running_fw += sd_fw[p];
      running_rv += sd_rv[p];
      sp_fw[p] = (running_fw > 0) ? static_cast<uint32_t>(running_fw) : 0;
      sp_rv[p] = (running_rv > 0) ? static_cast<uint32_t>(running_rv) : 0;
      uint32_t sp = sp_fw[p] + sp_rv[p];
      if (sc_position_is_testable(p, seq_length, minimum_clipped_bases)) {
        total_spanning_read_bases += 2ULL * static_cast<uint64_t>(sp);
        total_spanning_read_bases_fw += 2ULL * static_cast<uint64_t>(sp_fw[p]);
        total_spanning_read_bases_rv += 2ULL * static_cast<uint64_t>(sp_rv[p]);
        // N counts (position, direction) pairs with n_i > 0. A non-zero spanning count makes
        // both directions testable; positions with clips but no spanning reads are added in
        // the loop below, which only fires when read_through == 0, so nothing double counts.
        if (sp > 0) total_tested_positions += 2;
      }
    }
  }

  if (total_spanning_read_bases == 0) {
    WARN("No read spans any reference position with " + to_string(minimum_clipped_bases) +
         " aligned bases on both sides. --soft-clipping-minimum-bases is too large relative to "
         "the read lengths in this sample; soft-clipping (SC) evidence will not be meaningful.");
  }

  // Resolve each clipped position's tail consensus. This is done once, up front, because
  // both the null-rate estimate and the per-position rows below need agree_count.
  struct clipped_position {
    string   seq_id;
    uint32_t position;
    int32_t  direction;
    uint32_t clipped_count;
    uint32_t total_count;
    uint32_t agree_count;
    uint32_t agree_count_fw;   // agree_count split by the strand of the clipped read;
    uint32_t agree_count_rv;   //   agree_count_fw + agree_count_rv == agree_count
    uint32_t spanning_fw;      // read-through reads at this position, by strand;
    uint32_t spanning_rv;      //   spanning_fw + spanning_rv == total_count - clipped_count
    string   consensus_tail;   // always stored reference-forward
  };
  vector<clipped_position> positions;
  positions.reserve(
    [&]() -> size_t {
      size_t n = 0;
      for (map<string, map<uint32_t, uint32_t> >::iterator it = clipped_neg.begin(); it != clipped_neg.end(); it++) n += it->second.size();
      for (map<string, map<uint32_t, uint32_t> >::iterator it = clipped_pos.begin(); it != clipped_pos.end(); it++) n += it->second.size();
      return n;
    }());

  for (int32_t direction = -1; direction <= 1; direction += 2) {
    map<string, map<uint32_t, uint32_t> >& clipped = (direction < 0) ? clipped_neg : clipped_pos;
    map<string, map<uint32_t, vector<uint64_t> > >& tails = (direction < 0) ? tails_neg : tails_pos;

    for (map<string, map<uint32_t, uint32_t> >::iterator seq_it = clipped.begin(); seq_it != clipped.end(); seq_it++) {
      const string& seq_id = seq_it->first;
      // Same denominator for both directions.
      const vector<uint32_t>& cov_fw = spanning_fw[seq_id];
      const vector<uint32_t>& cov_rv = spanning_rv[seq_id];
      map<uint32_t, vector<uint64_t> >& seq_tails = tails[seq_id];

      for (map<uint32_t, uint32_t>::iterator pos_it = seq_it->second.begin(); pos_it != seq_it->second.end(); pos_it++) {
        clipped_position cp;
        cp.seq_id = seq_id;
        cp.position = pos_it->first;
        cp.direction = direction;
        cp.clipped_count = pos_it->second;
        cp.spanning_fw = (cp.position < cov_fw.size()) ? cov_fw[cp.position] : 0;
        cp.spanning_rv = (cp.position < cov_rv.size()) ? cov_rv[cp.position] : 0;
        uint32_t read_through = cp.spanning_fw + cp.spanning_rv;
        cp.total_count = cp.clipped_count + read_through;

        compute_clipped_tail_consensus(seq_tails[cp.position], consensus_bases,
                                       settings.soft_clipping_consensus_base_fraction,
                                       cp.consensus_tail, cp.agree_count,
                                       cp.agree_count_fw, cp.agree_count_rv);

        // Columns run outward from the wall. For a leading clip that is backwards
        // relative to the reference, so reverse it: the stored sequence is then always
        // "what should have continued the reference" reading in the forward direction.
        if (direction < 0) reverse(cp.consensus_tail.begin(), cp.consensus_tail.end());

        total_agreeing_clipped_read_ends += cp.agree_count;
        // Diagnostic only (reported in the SC gates table): how much of the agreeing clip
        // population sits at positions that saw only one strand. A clean library runs a few
        // percent; a run dominated by dark-cycle poly-G runs most of the way to 100%, which is
        // the single number that says "the SC calls in this run are an artifact".
        if ((cp.agree_count > 0) && ((cp.agree_count_fw == 0) || (cp.agree_count_rv == 0)))
          total_strand_pure_agreeing_clipped_read_ends += cp.agree_count;
        // A position with clips but no read-through was not counted in the prefix-sum loop.
        if (read_through == 0) total_tested_positions++;

        positions.push_back(cp);
      }
    }
  }

  // Store genome-wide totals in summary.
  // The raw rate keeps its historical meaning; the null rate the test actually uses is
  // the *agreeing* rate, since agree_count is the numerator being tested.
  uint64_t total_opportunities = total_clipped_read_ends + total_spanning_read_bases;
  summary.soft_clipping.total_spanning_read_bases = total_spanning_read_bases;
  // The genome-wide strand split of the read-through population. add_sc_evidence() falls back to
  // this as the expected strand ratio at a position with no read-through of its own, which is
  // exactly the frequency == 1.000 case -- otherwise those positions have an empty contingency
  // row and the strand test silently has no power where the clip count is highest.
  summary.soft_clipping.total_spanning_read_bases_forward = total_spanning_read_bases_fw;
  summary.soft_clipping.total_spanning_read_bases_reverse = total_spanning_read_bases_rv;
  summary.soft_clipping.total_clipped_read_ends = total_clipped_read_ends;
  summary.soft_clipping.total_agreeing_clipped_read_ends = total_agreeing_clipped_read_ends;
  summary.soft_clipping.total_strand_pure_agreeing_clipped_read_ends = total_strand_pure_agreeing_clipped_read_ends;
  summary.soft_clipping.soft_clipping_rate = (total_opportunities > 0)
    ? static_cast<double>(total_clipped_read_ends) / static_cast<double>(total_opportunities)
    : 0.0;

  double p0_raw = (total_opportunities > 0)
    ? static_cast<double>(total_agreeing_clipped_read_ends) / static_cast<double>(total_opportunities)
    : 0.0;
  double p0 = p0_raw;
  if ((settings.soft_clipping_minimum_rate > 0.0) && (p0 < settings.soft_clipping_minimum_rate)) {
    p0 = settings.soft_clipping_minimum_rate;
  }

  /*
   * Single-pass over-dispersion estimate.
   *
   * Pearson chi-square for the quasi-binomial, expanded so that zero-clip positions
   * need no separate term (their k=0 contribution is just n_i*p0/(1-p0), which the
   * same formula produces):
   *
   *   chi2 = sum_i (k_i - n_i p0)^2 / (n_i p0 (1-p0))
   *        = [ S_k2_over_n - 2*p0*S_k + p0^2 * S_n ] / (p0 * (1-p0))
   *
   * so only S_k and S_k2_over_n (over clipped positions) plus S_n and N (over all
   * positions) are needed -- no extra pass over the BAM or the difference arrays.
   *
   * Positions at or above the trim frequency are excluded from the fit: a genuine
   * breakpoint must not define the background it is judged against. Measured on real
   * data this matters a lot -- with the real breakpoints included the estimate is
   * inflated several-fold, and excluding them recovers the same null whether or not
   * they are present.
   */
  double S_k = 0.0, S_k2_over_n = 0.0;
  double S_n = static_cast<double>(total_opportunities);
  double S_n_trimmed = 0.0;
  uint64_t trimmed_positions = 0;

  for (vector<clipped_position>::const_iterator it = positions.begin(); it != positions.end(); it++) {
    if (it->total_count == 0) continue;
    double k = static_cast<double>(it->agree_count);
    double n = static_cast<double>(it->total_count);
    if ((settings.soft_clipping_dispersion_trim_frequency > 0.0) &&
        (k / n >= settings.soft_clipping_dispersion_trim_frequency)) {
      trimmed_positions++;
      S_n_trimmed += n;
      continue;
    }
    S_k += k;
    S_k2_over_n += k * k / n;
  }

  uint64_t N_prime = (total_tested_positions > trimmed_positions)
    ? (total_tested_positions - trimmed_positions) : 0;
  double S_n_prime = S_n - S_n_trimmed;
  double mean_n_prime = (N_prime > 0) ? S_n_prime / static_cast<double>(N_prime) : 0.0;

  double phi = 0.0;
  double rho = 0.0;
  if ((N_prime > 1) && (mean_n_prime > 1.0) && (p0 > 0.0) && (p0 < 1.0)) {
    double chi2 = (S_k2_over_n - 2.0 * p0 * S_k + p0 * p0 * S_n_prime) / (p0 * (1.0 - p0));
    phi = chi2 / static_cast<double>(N_prime - 1);
    rho = (phi - 1.0) / (mean_n_prime - 1.0);
  }
  double rho_raw = rho;
  if (!std::isfinite(rho) || (rho < 0.0)) rho = 0.0;          // degenerate -> plain binomial
  if ((settings.soft_clipping_maximum_dispersion > 0.0) &&
      (rho > settings.soft_clipping_maximum_dispersion)) {
    rho = settings.soft_clipping_maximum_dispersion;
  }
  if (rho < 1e-6) rho = 0.0;                                  // numerically identical to binomial

  summary.soft_clipping.soft_clipping_null_rate = p0;
  summary.soft_clipping.soft_clipping_dispersion = rho;
  summary.soft_clipping.soft_clipping_pearson_phi = phi;
  summary.soft_clipping.soft_clipping_tested_positions = N_prime;
  summary.soft_clipping.soft_clipping_trimmed_positions = trimmed_positions;
  summary.soft_clipping.soft_clipping_mean_tested_reads = mean_n_prime;

  // Write per-position counts file. Only positions with clip events are written; total_count is
  // the number of reads tested at that position = clipped reads + read-through reads.
  //
  // The leading '#' line records the parameters this file was tabulated under. The whole
  // tabulation sits behind the error-count done-file, so a stale file is the most likely way
  // to silently get wrong results; add_sc_evidence() checks these and errors out on a mismatch.
  ofstream out(settings.soft_clipping_counts_file_name.c_str());
  // The leading format token exists so that a counts file written by an older binary fails
  // loudly rather than being silently reinterpreted. Bump it whenever the meaning of a column
  // or of the denominator changes, not just when a column is added.
  out << "#sc_format=3"
      << "\tsoft_clipping_minimum_bases=" << minimum_clipped_bases
      << "\tconsensus_base_fraction=" << settings.soft_clipping_consensus_base_fraction
      << "\tnull_rate=" << p0
      << "\tdispersion=" << rho << "\n";
  out << "seq_id\tposition\tdirection\tclipped_count\ttotal_count\tagree_count\tclipped_sequence"
         "\tagree_count_forward\tagree_count_reverse\tspanning_forward\tspanning_reverse\n";

  for (vector<clipped_position>::const_iterator it = positions.begin(); it != positions.end(); it++) {
    out << it->seq_id << "\t" << it->position << "\t" << it->direction << "\t"
        << it->clipped_count << "\t" << it->total_count << "\t" << it->agree_count << "\t"
        << (it->consensus_tail.empty() ? "." : it->consensus_tail) << "\t"
        << it->agree_count_fw << "\t" << it->agree_count_rv << "\t"
        << it->spanning_fw << "\t" << it->spanning_rv << "\n";
  }

  out.close();

  cerr << "  Soft-clipping summary: " << total_clipped_read_ends << " clip events ("
       << total_agreeing_clipped_read_ends << " agreeing with their position consensus) over "
       << total_opportunities << " read opportunities" << endl;
  if (total_agreeing_clipped_read_ends > 0) {
    cerr << "    " << total_strand_pure_agreeing_clipped_read_ends << " ("
         << (100.0 * static_cast<double>(total_strand_pure_agreeing_clipped_read_ends)
                   / static_cast<double>(total_agreeing_clipped_read_ends))
         << "%) of the agreeing clip events are at positions that saw only one read strand" << endl;
  }
  cerr << "    null rate p0 = " << p0;
  if (p0 != p0_raw) cerr << " (raised from " << p0_raw << " by --soft-clipping-minimum-rate)";
  cerr << endl;
  cerr << "    dispersion: phi = " << phi << ", rho = " << rho_raw;
  if (rho != rho_raw) cerr << " (used " << rho << " after clamping)";
  cerr << endl;
  cerr << "    fit over " << N_prime << " positions (mean " << mean_n_prime
       << " reads), " << trimmed_positions << " trimmed at frequency >= "
       << settings.soft_clipping_dispersion_trim_frequency << endl;
}

} // namespace breseq

