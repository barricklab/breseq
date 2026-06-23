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

#ifndef _BRESEQ_IDENTIFY_REPEATS_H_
#define _BRESEQ_IDENTIFY_REPEATS_H_

#include "common.h"

using namespace std;

namespace breseq {

// One matching pair of approximately repeated regions found by a self-vs-self
// MUMmer alignment. seq_id_1/start_1/end_1 are always read in the forward
// orientation; seq_id_2/start_2/end_2 are normalized to start_2 < end_2 with
// strand_2 recording whether that side matched the reverse complement.
struct cRepeatMatch {
  string   seq_id_1;
  uint32_t start_1;
  uint32_t end_1;
  uint32_t length_1;
  string   seq_id_2;
  uint32_t start_2;
  uint32_t end_2;
  uint32_t length_2;
  int8_t   strand_2;
  double   percent_identity;
};

// One region that is a member of a repeat group: a set of regions in which
// every pair was found to match each other (a clique in the match graph).
// group_size == 1 means this region matched at least one other region, but
// its connected component of matches was not a clique, so (for now -- see
// group_repeats()) it was left ungrouped rather than merged with the others.
// min_percent_identity is the weakest pairwise identity within the group, and
// is only meaningful when group_size > 1.
struct cRepeatRegion {
  uint32_t group_id;
  uint32_t group_size;
  string   seq_id;
  uint32_t start;
  uint32_t end;
  uint32_t length;
  double   min_percent_identity;
};

// Runs a self-vs-self nucmer alignment of the given reference sequences,
// filters matches by length and percent identity, condenses them into
// groups of mutually-matching regions, and writes the result as a CSV.
void identify_repeats(
                      const vector<string>& reference_file_names,
                      const string& output_file_name,
                      const uint32_t minimum_length,
                      const double minimum_identity
                      );

// Parses the tab-delimited output of:
//   show-coords -r -T -l -H <delta_file>
// into one cRepeatMatch per line, before any filtering/deduplication.
vector<cRepeatMatch> parse_show_coords_file(const string& coords_file_name);

// Groups regions from a filtered, deduplicated match list into cliques: a
// group is only formed when every pair of its members appears in matches.
// Connected components of the match graph that are not themselves cliques
// are left ungrouped (each region becomes its own group of size 1) rather
// than guessed at -- this is the simple, non-fuzzy first pass.
vector<cRepeatRegion> group_repeats(const vector<cRepeatMatch>& matches);

void write_repeats_csv(const vector<cRepeatRegion>& regions, const string& file_name);
vector<cRepeatRegion> read_repeats_csv(const string& file_name);

}

#endif
