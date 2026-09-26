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

#include "original_coordinates.h"

namespace breseq {

void cOriginalCoordinateMap::reset_identity(int32_t length)
{
  m_blocks.clear();
  if (length > 0) m_blocks.push_back(cBlock(1, length, 1, +1));
  m_active = true;
}

int32_t cOriginalCoordinateMap::applied_length() const
{
  return m_blocks.empty() ? 0 : m_blocks.back().applied_end();
}

size_t cOriginalCoordinateMap::find_block(int32_t applied_pos_1) const
{
  ASSERT((applied_pos_1 >= 1) && (applied_pos_1 <= applied_length()),
         "Applied coordinate " + to_string(applied_pos_1) + " is outside of the mapped sequence (1-" + to_string(applied_length()) + ").");
  // First block whose start is past the position, then step back one
  size_t lo = 0, hi = m_blocks.size();
  while (lo < hi) {
    size_t mid = (lo + hi) / 2;
    if (m_blocks[mid].applied_start <= applied_pos_1) lo = mid + 1;
    else hi = mid;
  }
  return lo - 1;
}

size_t cOriginalCoordinateMap::split_at(int32_t applied_pos_1)
{
  if (applied_pos_1 == applied_length() + 1) return m_blocks.size();
  size_t i = find_block(applied_pos_1);
  cBlock& b = m_blocks[i];
  if (b.applied_start == applied_pos_1) return i;

  int32_t left_length = applied_pos_1 - b.applied_start;
  cBlock right(applied_pos_1, b.length - left_length, b.is_inserted() ? 0 : b.original_at(left_length), b.strand);
  b.length = left_length;
  m_blocks.insert(m_blocks.begin() + i + 1, right);
  return i + 1;
}

void cOriginalCoordinateMap::renumber_from(size_t index)
{
  int32_t next_start = (index == 0) ? 1 : m_blocks[index - 1].applied_end() + 1;
  for (size_t i = index; i < m_blocks.size(); i++) {
    m_blocks[i].applied_start = next_start;
    next_start += m_blocks[i].length;
  }
}

void cOriginalCoordinateMap::coalesce()
{
  if (m_blocks.size() < 2) return;
  vector<cBlock> merged;
  merged.reserve(m_blocks.size());
  merged.push_back(m_blocks[0]);
  for (size_t i = 1; i < m_blocks.size(); i++) {
    cBlock& a = merged.back();
    const cBlock& b = m_blocks[i];
    bool continuous;
    if (a.is_inserted() || b.is_inserted()) {
      continuous = a.is_inserted() && b.is_inserted();
    } else {
      continuous = (a.strand == b.strand) && (a.original_at(a.length) == b.original_start);
    }
    if (continuous) {
      a.length += b.length;
    } else {
      merged.push_back(b);
    }
  }
  m_blocks.swap(merged);
}

void cOriginalCoordinateMap::replace(int32_t start_1, int32_t end_1, int32_t new_length)
{
  if (!m_active) return;
  ASSERT(start_1 <= end_1, "Replaced range is empty: " + to_string(start_1) + "-" + to_string(end_1));
  int32_t old_length = end_1 - start_1 + 1;

  if (new_length < old_length) {
    // The bases after the ones that keep their coordinates are gone
    size_t i = split_at(start_1 + new_length);
    size_t j = split_at(end_1 + 1);
    m_blocks.erase(m_blocks.begin() + i, m_blocks.begin() + j);
    renumber_from(i);
    coalesce();
  } else if (new_length > old_length) {
    insert(end_1, new_length - old_length);
  }
}

void cOriginalCoordinateMap::insert(int32_t pos_1, int32_t length)
{
  if (!m_active || (length <= 0)) return;
  ASSERT((pos_1 >= 0) && (pos_1 <= applied_length()), "Insertion after applied coordinate " + to_string(pos_1) + " is outside of the mapped sequence.");
  size_t i = split_at(pos_1 + 1);
  m_blocks.insert(m_blocks.begin() + i, cBlock(pos_1 + 1, length, 0, +1));
  renumber_from(i + 1);
  coalesce();
}

void cOriginalCoordinateMap::invert(int32_t start_1, int32_t end_1)
{
  if (!m_active) return;
  ASSERT(start_1 <= end_1, "Inverted range is empty: " + to_string(start_1) + "-" + to_string(end_1));
  size_t i = split_at(start_1);
  size_t j = split_at(end_1 + 1);
  std::reverse(m_blocks.begin() + i, m_blocks.begin() + j);
  for (size_t k = i; k < j; k++) {
    cBlock& b = m_blocks[k];
    if (b.is_inserted()) continue;
    // The block is now read from its other end
    b.original_start = b.original_end();
    b.strand = -b.strand;
  }
  renumber_from(i);
  coalesce();
}

cOriginalCoordinate cOriginalCoordinateMap::original(int32_t applied_pos_1) const
{
  ASSERT(m_active, "Original coordinates were requested from a sequence that is not tracking them.");
  size_t i = find_block(applied_pos_1);
  const cBlock& b = m_blocks[i];
  if (!b.is_inserted()) {
    return cOriginalCoordinate(b.original_at(applied_pos_1 - b.applied_start), 0, b.strand);
  }
  // Inserted: anchor on the last original base before this run of inserted sequence
  size_t k = i;
  while ((k > 0) && m_blocks[k - 1].is_inserted()) k--;
  if (k == 0) {
    return cOriginalCoordinate(0, applied_pos_1, 0);
  }
  const cBlock& anchor = m_blocks[k - 1];
  return cOriginalCoordinate(anchor.original_end(), applied_pos_1 - anchor.applied_end(), 0);
}

void cOriginalCoordinateMap::write_header(ostream& out)
{
  out << "#seq_id\tapplied_start\tapplied_end\toriginal_start\toriginal_end\tstrand" << endl;
}

void cOriginalCoordinateMap::write(ostream& out, const string& seq_id) const
{
  for (vector<cBlock>::const_iterator it = m_blocks.begin(); it != m_blocks.end(); it++) {
    out << seq_id << "\t" << it->applied_start << "\t" << it->applied_end() << "\t";
    if (it->is_inserted()) {
      out << ".\t.\t.";
    } else {
      out << it->original_start << "\t" << it->original_end() << "\t" << ((it->strand == -1) ? "-" : "+");
    }
    out << endl;
  }
}

void cOriginalCoordinateMap::add_block(const cBlock& block)
{
  ASSERT(block.applied_start == applied_length() + 1,
         "Original coordinate blocks are not contiguous at applied coordinate " + to_string(block.applied_start) + ".");
  m_blocks.push_back(block);
  m_active = true;
}

} // namespace breseq
