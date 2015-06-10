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

#include "libbreseq/flagged_regions.h"

#include "libbreseq/reference_sequence.h"

namespace breseq {
  
cFlaggedRegions& cFlaggedRegions::read_mummer(const string& file_path, cReferenceSequences& ref_seq_info) {
  /*
   * Example command:
   * $ mummer -maxmatch -b -c -l 36 REL606.fna REL606.fna > exclude
   *
   * Input file format:
   > SeqID
   1544819    3595894       266
   3398753    4156036      200
   2712737    2713011        38
   ......
   > SeqID Reverse
   Start1     Start2    Length
   1544819    3595894       266
   3398753    4156036       200
   2712737    2713011        38
   *
   * Columns are: Start1     Start2    Length
   * In the first block, the match goes forward from both coords
   * In the second block, the match goes up forward the first coord and downward from the second coord
   */
  
  ifstream in(file_path.c_str());
  assert(in);
  
  bool on_reverse_block = false;
  
  string current_ref_seq_id;
  
  string line = "";
  while (getline(in, line)) {
    const vector<string>& tokens = split_on_whitespace(line);

    ASSERT(tokens.size() >= 2, "Unexpected number of items on line:\n" + line);
    uint32_t first  = from_string<uint32_t>(tokens[0]);
    uint32_t second = from_string<uint32_t>(tokens[1]);

    if (tokens.size() == 2) {
      ASSERT(tokens[0] == ">", "Expected line matching '> seq_id':\n" + line);
      current_ref_seq_id = tokens[1];
      continue;
    } else {
      ASSERT(tokens.size() >= 3, "Unexpected number of items on line:\n" + line);
      
      if ((tokens[0] == ">") && (tokens[2] == "Reverse")) {
        current_ref_seq_id = tokens[1];
        on_reverse_block = true;
        //cerr << "Found reverse block" << endl;
        continue;
      }
    }
    
    uint32_t size   = from_string<uint32_t>(tokens[2]);
    
    //ignore match to self!
    if (size==ref_seq_info.get_sequence_length(current_ref_seq_id))
      continue;
    
    if (!on_reverse_block) {
      this->flag_region(current_ref_seq_id, first,  first  + size - 1);
      this->flag_region(current_ref_seq_id, second, second + size - 1);
      string seq1 = ref_seq_info.get_sequence_1(current_ref_seq_id, first,  first  + size - 1);
      string seq2 = ref_seq_info.get_sequence_1(current_ref_seq_id, second, second + size - 1);
      ASSERT(seq1 == seq2, "Problem with line in MUMmer output file. Not repeat:\n" + line);
    } else {
      this->flag_region(current_ref_seq_id, first,  first  + size - 1);
      this->flag_region(current_ref_seq_id, (second + 1) - size, second);
      string seq1 = ref_seq_info.get_sequence_1(current_ref_seq_id, first,  first + size - 1);
      string seq2 = ref_seq_info.get_sequence_1(current_ref_seq_id, (second + 1) - size, second);
      seq2 = reverse_complement(seq2);
      ASSERT(seq1 == seq2, "Problem with line in MUMmer output file. Not repeat:\n" + line);
    }

  }
  
  return *this;
}


void cFlaggedRegions::write(const string& file_path) const {
  ofstream out(file_path.c_str());
  assert(out);
  
  out << "SeqID" << '\t' << "Start1" << '\t' << "End1" << endl;
  
  for (map<string,regions_t>::const_iterator itm = m_regions.begin(); itm != m_regions.end(); ++itm) {

    for (regions_t::iterator itr = itm->second.begin(); itr != itm->second.end(); ++itr) {
      out << itm->first << '\t' << itr->first << '\t' << itr->second << endl;
    }
  }
  
  out.close();
  
  return;
}

void cFlaggedRegions::print() const {
  
  cout << "SeqID" << '\t' << "Start1" << '\t' << "End1" << endl;

  for (map<string,regions_t>::const_iterator itm = m_regions.begin(); itm != m_regions.end(); ++itm) {
    for (regions_t::iterator itr = itm->second.begin(); itr != itm->second.end(); ++itr) {
      cout << itm->first << '\t' << itr->first << '\t' << itr->second << endl;
    }
  }
  return;
}

bool cFlaggedRegions::overlaps(const string& seq_id_1, const uint32_t pos_1, const string& seq_id_2, const region_t& region) const {
  return ((seq_id_1==seq_id_2) && (region.first <= pos_1) && (pos_1 <= region.second));
}

bool cFlaggedRegions::overlaps(const string& seq_id_1, const region_t& region_1, const string& seq_id_2, const region_t& region_2) const {
  
  if (seq_id_1 != seq_id_2)
    return false;
  
  region_t min = region_1 < region_2 ? region_1 : region_2;
  region_t max = region_1 > region_2 ? region_1 : region_2;
  
  return ((min.second >= max.first) && (max.second >= min.first));
}

  
void cFlaggedRegions::merge_within_distance(const uint32_t merge_distance)
{
  // Before merging!
  //cout << "Before merging:" << endl;
  //this->print();
  
  // save regions then clear current values
  for (map<string,regions_t>::const_iterator itm = m_regions.begin(); itm != m_regions.end(); ++itm) {
    const string& seq_id = itm->first;
    regions_t regions = this->all_regions(seq_id);
    this->m_regions[seq_id].clear();
    
    uint32_t last_begin = 0;
    uint32_t last_end = 0;

    for(regions_t::iterator it=regions.begin(); it != regions.end(); it++)
    {
      if (last_begin != 0) {
        if (last_end + merge_distance >= it->first) {
          last_end = it->second;
          continue;
        }
      
        flag_region(seq_id, last_begin, last_end);
      }
      
      last_begin = it->first;
      last_end = it->second;
    }
    
    // add the last one
    if (last_begin != 0) {
      flag_region(seq_id, last_begin, last_end);
    }
  }
  
  // After merging!
  //cout << endl << "After merging:" << endl;
  //this->print();
}
  

// Adds padding to both ends of every region, merging if necessary
// Padding can be negative!
void cFlaggedRegions::add_padding_to_ends(const int32_t padding, const int32_t minimum_size)
{
  // Before padding!
  //cout << "Before padding:" << endl;
  //this->print();
  
  // save regions then clear current values
  for (map<string,regions_t>::const_iterator itm = m_regions.begin(); itm != m_regions.end(); ++itm) {
    const string& seq_id = itm->first;
    regions_t regions = this->all_regions(seq_id);
    this->m_regions[seq_id].clear();
    for(regions_t::iterator it=regions.begin(); it != regions.end(); it++)
    {
      uint32_t new_start = (static_cast<int32_t>(it->first) > padding) ? static_cast<uint32_t>(static_cast<int32_t>(it->first) - padding) : 1;
      uint32_t new_end = static_cast<uint32_t>(static_cast<int32_t>(it->second) + padding); // we don't know anything about the ends of the sequence
      
      if (new_end >= new_start + minimum_size - 1) {
        flag_region(seq_id, new_start, new_end);
      }
    }
  }
  
  // After padding!
  //cout << endl << "After padding:" << endl;
  //this->print();
}

//<! Flags the specified region
cFlaggedRegions& cFlaggedRegions::flag_region(const string& seq_id, const uint32_t _start_1, const uint32_t _end_1) {
  
  uint32_t start_1 = _start_1;
  uint32_t end_1 = _end_1;
  ASSERT( start_1 <= end_1, "[start_1]: " + to_string<uint32_t>(start_1) + " is greater than [end_1]: " + to_string<uint32_t>(end_1));
  
  //cout << "adding:" << start_1 << " " << end_1 << endl;
  
  const regions_t& regions = this->regions_that_overlap(seq_id, start_1, end_1);
  
  //Remove and merge overlapping regions if needed.
  if (regions.size()) {
    start_1 = min(start_1, regions.begin()->first);
    end_1 = max(end_1, regions.rbegin()->second);
    this->remove(seq_id, regions);
  }
  
  //cout << "merged:" << start_1 << " " << end_1 << endl;
  
  m_regions[seq_id].insert(make_pair(start_1, end_1));
  
  return *this;
}

/* @JEB: Not tested
 cFlaggedRegions& cFlaggedRegions::unflag_region(uint32_t start_1, uint32_t end_1) {
 end_1 = end_1 == 0 ? start_1 : end_1;
 ASSERT(start_1 <= end_1, "[start_1]: " + s(start_1) + " is greater than [end_1]: " +s(end_1));
 
 regions_t regions = this->regions(start_1, end_1);
 
 this->remove(regions);
 
 //Add front segment if a region was partially overlapping.
 region_t front = *regions.begin();
 if (front.first < start_1) {
 this->flag_region(front.first, start_1 - 1);
 }
 
 //Add back segment if a region was partially overlapping.
 region_t back = *regions.rbegin();
 if (end_1 < back.second) {
 this->flag_region(end_1 + 1, back.second );
 }
 
 return *this;
 }
 */

  
//<! Returns whether the specified region is flagged (i.e., overlaps any region at all)
bool cFlaggedRegions::is_flagged(const string& seq_id, const cReferenceCoordinate& start_1, const cReferenceCoordinate& end_1) const {
  
  ASSERT(start_1 <= end_1, "[start_1]: " + s(start_1.get_position()) + ":" + s(start_1.get_insert_position()) + " is greater than [end_1]: " + s(end_1.get_position()) + ":" + s(end_1.get_insert_position()) );
  
  return this->regions_that_overlap(seq_id, start_1, end_1).size();
}

//<! Returns whether the specified region is flagged (i.e., overlaps any region at all)
bool cFlaggedRegions::is_flagged(const string& seq_id, const uint32_t start_1, const uint32_t end_1) const {
  
  return this->is_flagged(seq_id, cReferenceCoordinate(start_1), cReferenceCoordinate(end_1));
}

//<! Removes the specified regions
cFlaggedRegions& cFlaggedRegions::remove(const string& seq_id, const regions_t& regions) {
  for (regions_t::iterator it = regions.begin(); it != regions.end(); ++it) {
    m_regions[seq_id].erase(*it);
  }
  return *this;
}

//<! Returns all of the regions that overlap the specified coordinates
cFlaggedRegions::regions_t cFlaggedRegions::regions_that_overlap(const string& seq_id, const cReferenceCoordinate& start_1, const cReferenceCoordinate& end_1) const {
  
  // create empty item and return if seq_id does not exist
  if (!m_regions.count(seq_id))
    return regions_t();
  
  ASSERT(start_1 <= end_1, "[start_1]: " + s(start_1.get_position()) + ":" + s(start_1.get_insert_position()) + " is greater than [end_1]: " + s(end_1.get_position()) + ":" + s(end_1.get_insert_position()) );
  
  // Return all regions by default
  if ((start_1 == cReferenceCoordinate(0,0)) && (end_1 == cReferenceCoordinate(0,0))) {
    return m_regions.at(seq_id);
  }

  // find all overlapping ones
  
  // Should be something clever we can do here to use binary search
  //regions_t::iterator it_lower = m_regions.at(seq_id).lower_bound( make_pair<uint32_t, uint32_t>(start_1, 0) );
  //if (it_lower != m_regions.at(seq_id).begin())
  //  it_lower--;
  
  regions_t overlapping;
  
  for (regions_t::iterator it = m_regions.at(seq_id).begin(); it != m_regions.at(seq_id).end(); it++) {
    
    if ((start_1 >= cReferenceCoordinate(it->first)) && (start_1 <= cReferenceCoordinate(it->second))) {
      overlapping.insert(*it);
    }
    else if ((cReferenceCoordinate(it->first) >= start_1) && (cReferenceCoordinate(it->first) <= end_1)) {
      overlapping.insert(*it);

    }
  }
  
  return overlapping;
}
  
  
//<! Returns whether the specified region is contained within a region (i.e., completely overlaps any region at all)
bool cFlaggedRegions::is_contained(const string& seq_id, const cReferenceCoordinate& start_1, const cReferenceCoordinate& end_1) const {
  
  ASSERT(start_1 <= end_1, "[start_1]: " + s(start_1.get_position()) + ":" + s(start_1.get_insert_position()) + " is greater than [end_1]: " + s(end_1.get_position()) + ":" + s(end_1.get_insert_position()) );
  
  return this->regions_that_contain(seq_id, start_1, end_1).size();
}
  
//<! Returns all of the regions that overlap the specified coordinates
cFlaggedRegions::regions_t cFlaggedRegions::regions_that_contain(const string& seq_id, const cReferenceCoordinate& start_1, const cReferenceCoordinate& end_1) const {
  
  // create empty item and return if seq_id does not exist
  if (!m_regions.count(seq_id))
    return regions_t();
  
  ASSERT(start_1 <= end_1, "[start_1]: " + s(start_1.get_position()) + ":" + s(start_1.get_insert_position()) + " is greater than [end_1]: " + s(end_1.get_position()) + ":" + s(end_1.get_insert_position()) );
  
  // Return all regions by default
  if ((start_1 == cReferenceCoordinate(0,0)) && (end_1 == cReferenceCoordinate(0,0))) {
    return m_regions.at(seq_id);
  }
  
  // find all containing ones
  regions_t contained;
  
  for (regions_t::iterator it = m_regions.at(seq_id).begin(); it != m_regions.at(seq_id).end(); it++) {
    
    if ((start_1 >= cReferenceCoordinate(it->first)) && (end_1 <= cReferenceCoordinate(it->second))) {
      contained.insert(*it);
    }
  }
  
  return contained;
}


}//namespace bresesq