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

#ifndef _BRESEQ_HOMOLOGOUS_DELETION_H_
#define _BRESEQ_HOMOLOGOUS_DELETION_H_

#include "common.h"
#include "reference_sequence.h"
#include "genome_diff.h"
#include "alignment.h"
#include "pileup_base.h"
#include "pileup.h"

namespace breseq {

  //! Machinery for locating a deletion that occurred between two near-identical copies of an
  //  UNANNOTATED repeat -- paralogous genes, an rRNA operon, any homology that carries no
  //  repeat_region feature.
  //
  //  Such a deletion leaves no split-read junction at all: a read crossing the breakpoint aligns just
  //  as well to either copy, so candidate junction identification has nothing to find and
  //  predictMCplusJCtoDEL cannot convert the missing coverage into a mutation. What survives in the
  //  genome is a single HYBRID copy -- the left copy up to the crossover, the right copy after it.
  //
  //  The signal is in the read alignments themselves. Let d be the register (the deletion size), so
  //  that reference position x in the left copy is homologous to x+d in the right copy, and let p be
  //  the crossover. At a DISCRIMINATING column -- one where the two copies actually differ -- the
  //  hybrid carries the left copy's base if x < p and the right copy's base if x >= p. Because a read
  //  that cannot be told apart between the copies is placed at BOTH of them, pooling the pile at x
  //  with the pile at x+d recovers that single series at full depth from the surviving flank. Reads
  //  that span the breakpoint are not required, which is what makes this work for short single-end
  //  data.
  //
  //  It also means the test is self-nulling. If both copies are still present, every discriminating
  //  column reads ~50/50 and nothing is called. If a deletion removed one whole copy but was not
  //  mediated by the homology, every column calls the same allele and there is no flip.
  //
  //  IMPORTANT: this evidence lives in reads that identify_mutations deliberately throws away. A read
  //  that ties between the two copies is marked redundant, and base calling skips redundant reads, so
  //  no RA evidence is ever emitted for it. With long reads enough crossing reads break the tie and
  //  land uniquely (which is why the RA-based path works on 150 bp data), but with 36 bp reads
  //  essentially all of them tie. The counter below therefore goes to the BAM and counts redundant
  //  placements too.

  //! A/C/G/T counts at one reference column, INCLUDING redundantly placed reads.
  struct base_tally {
    uint32_t n[4];      //!< indexed by basechar2index(): A C G T
    uint32_t total;     //!< n[0..3] summed, i.e. bases that passed the quality gate

    base_tally() : total(0) { n[0] = n[1] = n[2] = n[3] = 0; }
  };
  typedef map<int32_t, base_tally> base_tally_map_t;

  //! One homologous column pair, pooled across both copies.
  struct homology_column {
    int32_t  position;      //!< 1-indexed, in LEFT-copy coordinates
    char     left_base;     //!< seq[position]
    char     right_base;    //!< seq[position + d]; differs from left_base by construction
    uint32_t left_allele;   //!< pooled count of left_base  at position and at position + d
    uint32_t right_allele;  //!< pooled count of right_base at position and at position + d
    uint32_t total;         //!< pooled total at both columns
    int32_t  call;          //!< -1 hybrid shows the left allele, +1 the right allele, 0 uncalled

    homology_column()
    : position(0), left_base('N'), right_base('N')
    , left_allele(0), right_allele(0), total(0), call(0) {}
  };
  typedef vector<homology_column> homology_columns_t;

  //! The alleles one read placement carries, in ascending LEFT-copy coordinate.
  typedef vector<pair<int32_t, char> > read_allele_list_t;
  typedef map<string, read_allele_list_t> read_alleles_t;

  //! Counts reference-column bases from data/reference.bam, INCLUDING redundantly placed reads.
  //
  //  Deliberately different from identify_mutations_pileup in three ways, each of which matters here:
  //
  //   1. Redundancy is neither filtered nor down-weighted. A redundantly placed read is genuine
  //      evidence about the hybrid at BOTH homologous columns, and those columns are pooled anyway.
  //      This is the whole point -- identify_mutations discards exactly these reads.
  //   2. Trimmed bases are not excluded. Read-end trimming exists to stop breseq calling a variant
  //      from a read end whose placement among repeat copies is ambiguous, but that ambiguity is
  //      precisely the quantity being measured here, and pooling both copies makes it harmless.
  //      Excluding trimmed bases would delete most of the signal inside a long, highly identical block.
  //   3. Mapping quality is not filtered, matching identify_mutations (whose MAPQ guard is commented
  //      out). Redundant placements have low MAPQ by construction, so a MAPQ filter would just be the
  //      redundancy filter again.
  //
  //  Base quality IS gated, at the same settings.base_quality_cutoff bar identify_mutations uses.
  //
  //  Note the pileup returns at most 8000 reads per column. That is far above these coverages, and it
  //  degrades gracefully in any case: a truncated pile is still a fair sample of the allele ratio.
  class homology_base_counter : public pileup_base {
  public:
    homology_base_counter(const string& bam, const string& fasta, uint8_t base_quality_cutoff)
    : pileup_base(bam, fasta), m_position_offset(0)
    , m_region_start(0), m_region_end(0), m_base_quality_cutoff(base_quality_cutoff) {}

    //! Accumulate one region. Safe to call repeatedly; tallies accumulate until clear().
    void count_region(const string& seq_id, int32_t start_1, int32_t end_1);

    void clear() { m_tallies.clear(); m_read_alleles.clear(); m_record_positions.clear(); m_position_offset = 0; }

    //! Also record per-read alleles at these reference positions, shifting them by -offset so they
    //  come back in left-copy coordinates. Set before each count_region() call.
    void set_read_recording(const set<int32_t>& positions, int32_t offset)
      { m_record_positions = positions; m_position_offset = offset; }

    const base_tally_map_t& tallies() const { return m_tallies; }
    const read_alleles_t& read_alleles() const { return m_read_alleles; }

    virtual void pileup_callback(const pileup& p);

  protected:
    base_tally_map_t m_tallies;
    read_alleles_t   m_read_alleles;
    set<int32_t>     m_record_positions;
    int32_t          m_position_offset;
    int32_t          m_region_start;
    int32_t          m_region_end;
    uint8_t          m_base_quality_cutoff;
  };

  //! Fraction of positions x in [lo,hi] with seq[x] == seq[x+d]. Returns 0 if the range is empty
  //  after clamping both x and x+d into the sequence.
  double homology_identity(const cAnnotatedSequence& seq, int32_t d, int32_t lo, int32_t hi);

  //! Best register d > 0 aligning the neighborhood of either anchor to its homologous partner.
  //
  //  Scores each candidate d as max(identity around anchor_left, identity around anchor_right - d):
  //  either MC boundary may have landed outside the homologous block, so requiring both to agree
  //  loses real cases. Returns 0 unless the winner clears min_identity AND beats the best register at
  //  least margin_exclusion away by min_margin -- between paralogs that are identical in long
  //  stretches everywhere, a high score on its own means nothing.
  int32_t find_homologous_register(const cAnnotatedSequence& seq,
                                   int32_t anchor_left, int32_t anchor_right, int32_t d_seed,
                                   int32_t search_half_width, int32_t window,
                                   double min_identity, double min_margin, int32_t margin_exclusion,
                                   double& best_identity, double& runner_up_identity);

  //! Grow [block_start, block_end] (LEFT-copy coordinates) outward from seed while each successive
  //  window of `step` bases keeps identity >= min_identity at register d.
  void homologous_block_extent(const cAnnotatedSequence& seq, int32_t d, int32_t seed,
                               int32_t step, double min_identity, int32_t max_extent,
                               int32_t& block_start, int32_t& block_end);

  //! Positions x in [block_start, block_end] where seq[x] != seq[x+d], ascending.
  void discriminating_positions(const cAnnotatedSequence& seq, int32_t d,
                                int32_t block_start, int32_t block_end,
                                vector<int32_t>& positions);

  //! Pool the BAM tallies at x and x+d into a per-column allele call.
  void build_homology_columns(const cAnnotatedSequence& seq, int32_t d,
                              const vector<int32_t>& positions,
                              const base_tally_map_t& tallies,
                              uint32_t min_column_coverage, double min_allele_fraction,
                              homology_columns_t& columns);

  //! Same, but filled from RA evidence rather than the BAM, for callers that have no BAM (gdtools).
  //  Only positions where an RA exists are informative; everything else stays uncalled, which
  //  reproduces the bracketing the RA-only implementation used to do.
  void build_homology_columns_from_ra(const cAnnotatedSequence& seq, int32_t d,
                                      const vector<int32_t>& positions,
                                      const map<int32_t, string>& ra_base,
                                      homology_columns_t& columns);

  //! Locate the single left->right flip in the called columns, choosing the split that misclassifies
  //  the fewest of them. The crossover lies in (last_left, first_right]. Returns false unless both
  //  sides carry at least min_columns_each_side calls and the flip is clean to within max_violations.
  bool find_crossover_interval(const homology_columns_t& columns,
                               int32_t min_columns_each_side, int32_t max_violations,
                               int32_t& last_left, int32_t& first_right,
                               int32_t& n_left, int32_t& n_right, int32_t& n_violations);

  //! Fraction of read placements whose own alleles form a single left->right step.
  //
  //  This is the direct statement that no sequencing error interrupts the repeat structure: a read
  //  that shows the right copy's allele and then the left copy's allele further along cannot have
  //  come from the hybrid. Only placements covering at least two discriminating columns are counted;
  //  informative_reads returns how many those were, so the caller can tell "all consistent" from
  //  "nothing to check".
  double read_allele_consistency(const cAnnotatedSequence& seq, int32_t d,
                                 const read_alleles_t& read_alleles,
                                 int32_t& informative_reads);

  //! DEL with seq_id/position/size/frequency, plus the between=<f1>,<f2> tag naming the two
  //  homologous features the deletion sits between.
  cDiffEntry make_homologous_deletion(cAnnotatedSequence& seq, int32_t position, int32_t size,
                                      const vector<string>& evidence_ids,
                                      bool polymorphism_prediction);

} // namespace breseq

#endif
