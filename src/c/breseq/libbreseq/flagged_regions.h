/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2008-2010 Michigan State University
 Copyright (c) 2011-2012 The University of Texas at Austin
 
 breseq is free software; you can redistribute it and/or modify it under the
 terms the GNU General Public License as published by the Free Software
 Foundation; either version 1, or (at your option) any later version.
 
 *****************************************************************************/

#ifndef _BRESEQ_FLAGGED_REGIONS_H_
#define _BRESEQ_FLAGGED_REGIONS_H_

#include "common.h"

namespace breseq {

// Predefs
class cReferenceSequences;
class cReferenceCoordinate;
  
/* Helper class for a set of start-end regions
 *
 * Handles loading from MUMmer results,
 * padding, merging, and testing for overlap
 */
  
class cFlaggedRegions  {
public:
  typedef pair<uint32_t, uint32_t>  region_t;
  typedef set<region_t>             regions_t;
  
  cFlaggedRegions()
  :m_regions() {
    return;
  }
  
  //! I/O.
  cFlaggedRegions& read_mummer(const string& file_path, cReferenceSequences& ref_seq);
  void write(const string& file_path) const;
  void print() const;
  
  //! Add region to be marked.
  cFlaggedRegions& flag_region(const string& seq_id, const uint32_t _start_1, const uint32_t _end_1);
  
  //! Extend the ends of every region and remove small regions if necessary
  void add_padding_to_ends(const int32_t padding, const int32_t minimum_size = 0);
  
  //! Extend the ends of every region and merge if necessary
  void merge_within_distance(const uint32_t merge_distance);
  
  //! Remove overlapping regions, adds segments if partial overlapping occurs.
  //cFlaggedRegions& unflag_region(uint32_t start_1, uint32_t end_1 = 0);
  
  //! Tests if start_1 to end_1 spans over a marked region.
  bool is_flagged(const string& seq_id, const uint32_t _start_1, const uint32_t _end_1 = 0) const;
  bool is_flagged(const string& seq_id, const cReferenceCoordinate& _start_1, const cReferenceCoordinate& _end_1) const;
  
  //! Test if two regions overlap each other.
  bool overlaps(const string& seq_id_1, const region_t& region_1, const string& seq_id_2, const region_t& region_2) const;
  
  //! Test if a position is within a region.
  bool overlaps(const string& seq_id_1, const uint32_t pos_1, const string& seq_id_2, const region_t& region) const;
  
  //! Return overlapping regions, defaults to all regions.
  regions_t regions_that_overlap(const string& seq_id, const cReferenceCoordinate& _start_1, const cReferenceCoordinate& _end_1) const;
  
  //! Return containing regions, must completely lie within the region
  bool is_contained(const string& seq_id, const cReferenceCoordinate& start_1, const cReferenceCoordinate& end_1) const;
  regions_t regions_that_contain(const string& seq_id, const cReferenceCoordinate& _start_1, const cReferenceCoordinate& _end_1) const;

  //! return all regions for the specified seq_id
  regions_t all_regions(const string& seq_id) {
      return (m_regions.count(seq_id)) ? m_regions.at(seq_id) : regions_t();
  }
  
  //! Remove regions.
  cFlaggedRegions& remove(const string& seq_id, const regions_t& regions);
  
  list<string> get_seq_ids() const
  {
    list<string> return_list;
    for (map<string,regions_t>::const_iterator it = m_regions.begin(); it != m_regions.end(); it++) {
      return_list.push_back(it->first);
    }
    return return_list;
  }
  
protected:
  map<string,regions_t> m_regions;
};

}
#endif