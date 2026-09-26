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

#ifndef _BRESEQ_ORIGINAL_COORDINATES_H_
#define _BRESEQ_ORIGINAL_COORDINATES_H_

#include "common.h"

namespace breseq {

/* Where a base of a mutated ("applied") sequence came from in the original sequence.
 *
 * A base that the original sequence contains has offset 0 and position = its coordinate there.
 * A base of newly inserted sequence has no original coordinate of its own: position is then the
 * coordinate of the last original base before the insertion (0 when there is none) and offset >= 1
 * says that this is the offset-th inserted base after it.
 */
struct cOriginalCoordinate {
  int32_t position;
  int32_t offset;
  int8_t  strand;     // +1 or -1 for an original base (-1 inside an applied inversion); 0 if inserted

  cOriginalCoordinate() : position(0), offset(0), strand(0) {}
  cOriginalCoordinate(int32_t _position, int32_t _offset, int8_t _strand)
    : position(_position), offset(_offset), strand(_strand) {}

  bool is_inserted() const { return offset > 0; }
};

/* Correspondence between the coordinates of a sequence that has had mutations applied to it and
 * the coordinates of the sequence they were applied to.
 *
 * Kept as a sorted list of blocks that together cover the whole applied sequence. The map is
 * inactive (empty) until reset_identity() is called, and every edit is then a no-op, so a sequence
 * that is not being tracked pays nothing. The three sequence-editing primitives of
 * cAnnotatedSequence call replace(), insert() and invert() with the same arguments they apply to
 * the bases, so every mutation type is covered without knowing what a mutation is.
 */
class cOriginalCoordinateMap {
public:

  struct cBlock {
    int32_t applied_start;     // 1-based coordinate in the applied sequence
    int32_t length;
    int32_t original_start;    // original coordinate of the FIRST applied base of the block; 0 = inserted sequence
    int8_t  strand;            // +1: original increases along the block; -1: decreases (inverted)

    cBlock() : applied_start(0), length(0), original_start(0), strand(1) {}
    cBlock(int32_t _applied_start, int32_t _length, int32_t _original_start, int8_t _strand)
      : applied_start(_applied_start), length(_length), original_start(_original_start), strand(_strand) {}

    bool    is_inserted() const { return original_start == 0; }
    int32_t applied_end() const { return applied_start + length - 1; }
    // Original coordinate of the i-th base of the block (0-based index)
    int32_t original_at(int32_t i) const { return is_inserted() ? 0 : original_start + strand * i; }
    int32_t original_end() const { return original_at(length - 1); }
  };

  cOriginalCoordinateMap() : m_active(false) {}

  bool is_active() const { return m_active; }

  // Start tracking: every base maps to itself
  void reset_identity(int32_t length);
  void clear() { m_blocks.clear(); m_active = false; }

  // The applied bases start_1..end_1 are replaced by new_length bases. The first min(old, new)
  // bases keep their original coordinates (a substituted base is still that base of the original),
  // the rest of the old range is deleted, and any surplus new bases are inserted sequence.
  void replace(int32_t start_1, int32_t end_1, int32_t new_length);

  // length bases of inserted sequence are added AFTER applied position pos_1 (0 = at the start)
  void insert(int32_t pos_1, int32_t length);

  // The applied bases start_1..end_1 are reverse complemented in place
  void invert(int32_t start_1, int32_t end_1);

  cOriginalCoordinate original(int32_t applied_pos_1) const;

  int32_t applied_length() const;
  const vector<cBlock>& blocks() const { return m_blocks; }

  // Appends the block list to a file; one line per block, tab-delimited:
  //   seq_id  applied_start  applied_end  original_start  original_end  strand
  // with '.' in the original and strand columns of inserted sequence.
  static void write_header(ostream& out);
  void write(ostream& out, const string& seq_id) const;
  // Adds one block read back from such a line (activates the map)
  void add_block(const cBlock& block);

private:
  vector<cBlock> m_blocks;
  bool m_active;

  // Index of the block containing applied_pos_1
  size_t find_block(int32_t applied_pos_1) const;
  // Makes a block start at applied_pos_1 (splitting one if needed) and returns its index, which is
  // m_blocks.size() when applied_pos_1 is one past the end of the sequence.
  size_t split_at(int32_t applied_pos_1);
  // Recomputes applied_start for every block from index onward
  void renumber_from(size_t index);
  // Merges adjacent blocks that describe one continuous run
  void coalesce();
};

} // namespace breseq

#endif
