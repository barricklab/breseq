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

#include "mutation_predictor.h"

#include "output.h"
#include "dp_evidence.h"
#include "homologous_deletion.h"
#include "identify_mutations.h"
#include "resolve_alignments.h"

using namespace std;

namespace breseq {

  cReferenceSequences MutationPredictor::ref_seq_info;
  
	MutationPredictor::MutationPredictor(cReferenceSequences& _ref_seq_info)
	{
		 ref_seq_info = _ref_seq_info;
	}

	// Private methods
	cFeatureLocation* MutationPredictor::within_repeat(string seq_id, int32_t position)
	{
		cFeatureLocationList& repeat_list = ref_seq_info[seq_id].m_repeat_locations;
    cFeatureLocation* repeat= NULL;
    
    // by returning the last one we encounter that we are inside, 
    // we get the inner repeat in nested cases
    for(cFeatureLocationList::iterator it = repeat_list.begin(); it != repeat_list.end(); it++) {
      cFeatureLocation& test_repeat = *it;
			if ((test_repeat.get_start_1() <= position) && (position <= test_repeat.get_end_1()))
				repeat = &test_repeat;
    }
		return repeat;
	}

	bool MutationPredictor::sort_by_hybrid(const counted_ptr<cDiffEntry>& a, const counted_ptr<cDiffEntry>& b)
	{
    int32_t a_pos = n(a->entry_exists("_side_1_is") ? (*a)["side_2_position"] : (*a)["side_1_position"]);
		int32_t b_pos = n(b->entry_exists("_side_1_is") ? (*b)["side_2_position"] : (*b)["side_1_position"]);

		int32_t a_seq_order = (a->entry_exists("_side_1_is") ? ref_seq_info.seq_order[(*a)["side_2_seq_id"]] : ref_seq_info.seq_order[(*a)["side_1_seq_id"]]);
		int32_t b_seq_order = (b->entry_exists("_side_1_is") ? ref_seq_info.seq_order[(*b)["side_2_seq_id"]] : ref_seq_info.seq_order[(*b)["side_1_seq_id"]]);

		uint32_t a_reject_order = a->number_reject_reasons();
		uint32_t b_reject_order = b->number_reject_reasons();

		// sort by seq_id, position, fewer reject reasons, then score (highest to lowest)
    
    if (a_seq_order != b_seq_order) 
      return (a_seq_order < b_seq_order);
    if (a_pos != b_pos) 
      return (a_pos < b_pos);
    if (a_reject_order != b_reject_order) 
      return (a_reject_order < b_reject_order);

    if (a->entry_exists("pos_hash_score") && b->entry_exists("pos_hash_score")) {
      return (n((*a)["pos_hash_score"]) > n((*b)["pos_hash_score"]));
    }
    
    return false;
  }

	bool MutationPredictor::sort_by_reject_score(const counted_ptr<cDiffEntry>& a, const counted_ptr<cDiffEntry>& b)
	{
		uint32_t a_reject_order = a->number_reject_reasons();
		uint32_t b_reject_order = b->number_reject_reasons();

		// sort by seq_id, position, fewer reject reasons, then score (highest to lowest)
    if (a_reject_order != b_reject_order) 
      return (a_reject_order < b_reject_order);
    
    if (a->entry_exists("pos_hash_score") && b->entry_exists("pos_hash_score")) {
      return (n((*a)["pos_hash_score"]) > n((*b)["pos_hash_score"]));
    }
    
    return false;
	}

	// Look at SNPs and small indels predicted by read alignments.
	// Make sure they are sorted by position.
	bool MutationPredictor::sort_by_pos(const counted_ptr<cDiffEntry>& a, const counted_ptr<cDiffEntry>& b)
	{
		return (
			((*a)["seq_id"] != (*b)["seq_id"]) ? ((*a)["seq_id"] < (*b)["seq_id"]) :
			((*a)["position"] != (*b)["position"]) ? (n((*a)["position"]) < n((*b)["position"])) :
			(n((*a)["insert_position"]) < n((*b)["insert_position"]))
		);
	}

  /*
	 Title   : prepare_junctions
	 Function: Adds fields to junction items in preparation for mutation prediction
             And mark junctions that are near the ends of contigs
   */
  void MutationPredictor::prepare_junctions(Settings& settings, Summary& summary, cGenomeDiff& gd)
  {
    (void) settings;
    (void) summary;
    
    // For all that follows, we need information about repeat_regions overlapping the sides of junctions
    vector<gd_entry_type> jc_types = make_vector<gd_entry_type>(JC);
		diff_entry_list_t jc = gd.get_list(jc_types);
    
    const int32_t max_distance_to_repeat = 50;
    
		for (diff_entry_list_t::iterator jc_it=jc.begin(); jc_it!=jc.end(); jc_it++)
		{
			cDiffEntry& j = **jc_it;
      
			j["_side_1_read_side"] = "-1";
			j["_side_2_read_side"] = "1";
      
			for (uint32_t side = 1; side <= 2; side++)
			{
				string side_key = "side_" + s(side);
        int32_t this_max_distance_to_repeat = max_distance_to_repeat;
				cFeatureLocation* is = ref_seq_info.find_closest_repeat_region_boundary(
                                                                                  n(j[side_key + "_position"]),
                                                                                  ref_seq_info[j[side_key + "_seq_id"]].m_repeats,
                                                                                  this_max_distance_to_repeat,
                                                                                  n(j[side_key + "_strand"])
                                                                                  );
				if (is != NULL)
				{
					j["_" + side_key + "_is"] = "1";
					j["_" + side_key + "_is_start"] = s(is->get_start_1());
					j["_" + side_key + "_is_end"] = s(is->get_end_1());
          j["_" + side_key + "_is_name"] = (*(is->get_feature()))["name"];
          j["_" + side_key + "_is_strand"] = s(is->get_strand());
          j["_" + side_key + "_is_distance"] = s(this_max_distance_to_repeat);
				}
				
				j[side_key + "_annotate_key"] = (j.entry_exists("_" + side_key + "_is_start") || (j.entry_exists(side_key + "_redundant") && n(j[side_key + "_redundant"]))) ? "repeat" : "gene";
			}
      
			// by default, we are sorted by this coord
			j["_unique_interval"] = "side_1";
      
			// Determine which side of the junction is the IS and which is unique
			// these point to the correct initial interval...
			if (j.entry_exists("_side_1_is"))
			{
        //cout << n(j["_side_1_is_start"]) << " " << n(j["_side_1_is_end"]) << " " << n(j["side_1_position"]) << endl;
        
				if (abs(n(j["_side_1_is_start"]) - n(j["side_1_position"])) <= max_distance_to_repeat)
				{
					j["_is_interval"] = "side_1";
					j["_is_interval_closest_side_key"] = "start";
					j["_unique_interval"] = "side_2";
				}
				else if (abs(n(j["_side_1_is_end"]) - n(j["side_1_position"])) <= max_distance_to_repeat)
				{
					j["_is_interval"] = "side_1";
					j["_is_interval_closest_side_key"] = "end";
					j["_unique_interval"] = "side_2";
				}
			}
      
      if (j.entry_exists("_side_2_is"))
			{
        //cout << n(j["_side_2_is_start"]) << " " << n(j["_side_2_is_end"]) << " " << n(j["side_2_position"]) << endl;
        
				if (abs(n(j["_side_2_is_start"]) - n(j["side_2_position"])) <= max_distance_to_repeat)
				{
					j["_is_interval"] = "side_2";
					j["_is_interval_closest_side_key"] = "start";
					j["_unique_interval"] = "side_1";
				}
				else if (abs(n(j["_side_2_is_end"]) - n(j["side_2_position"])) <= max_distance_to_repeat)
				{
					j["_is_interval"] = "side_2";
					j["_is_interval_closest_side_key"] = "end";
					j["_unique_interval"] = "side_1";
				}
			}
      
      // add in the rest of the unique interval information
      j["_unique_interval_seq_id"] = j[j["_unique_interval"] + "_seq_id"];
      j["_unique_interval_position"] = j[j["_unique_interval"] + "_position"];
      j["_unique_interval_strand"] = j[j["_unique_interval"] + "_strand"];
    }
    
  }
  
  void MutationPredictor::predictMCplusJCtoDEL(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc)
  {
    (void) summary;
    bool verbose = false; // for debugging
    
    int32_t read_length_max = summary.sequence_conversion.read_length_max;
    
		// DEL prediction:
		// (1) there is a junction that exactly crosses the deletion boundary
		// (2) there is no junction, but both ends of the deletion are in repeat sequences
		// (3) there is a junction between unique sequence and a repeat element
    
    if (verbose)
      cout << "DEL PREDICTION" << endl;
    
    for(diff_entry_list_t::iterator mc_it = mc.begin(); mc_it != mc.end(); mc_it++)
    {
      cDiffEntry& mc_item = **mc_it;
      
      if (verbose)
        cout << mc_item << endl;
      
			if (mc_item.entry_exists("reject"))
			  continue;
      
			// set up generic deletion item
			cDiffEntry mut;
      mut._type = DEL;
      mut._evidence = make_vector<string>(mc_item._id);
      int32_t size = n(mc_item["end"]) - n(mc_item["start"]) + 1;
			mut
      ("seq_id", mc_item["seq_id"])
      ("position", mc_item["start"])
      ("size", s(size));
			;
      
      if (settings.polymorphism_prediction) {
        mut[FREQUENCY] = "1";
      }
      
			///
			// (0) this is a deletion of an entire fragment
			///     
      
      uint32_t tid = ref_seq_info.seq_id_to_index(mut[SEQ_ID]); 
      if ( (n(mut[POSITION]) == 1) && (n(mut[POSITION]) + n(mut["size"]) - 1 == static_cast<int32_t>(ref_seq_info[tid].m_length)) )
      {
        gd.add(mut);
        continue;
      }
      
      
			///
			// (1) there is a junction that exactly crosses the deletion boundary 
			///
      ///
      if (verbose)
        cout << "(1) Checking for JC with same boundaries as MC." << endl;
      
      bool done = false;
			for(diff_entry_list_t::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) //JC
			{
				cDiffEntry& jc_item = **jc_it;
        
        if (verbose)
          cout << "  " << jc_item << endl;
        
				if (jc_item["side_1_seq_id"] != mut["seq_id"] || jc_item["side_2_seq_id"] != mut["seq_id"])  {
          continue;  
        }
        
        // Due to overlap resolution, this can change from the time when we sorted this way... *Ugh* .. 
        // so we fix it here so that the lower coordinate part of the junction is first.
        
        int32_t side_1_position = n(jc_item["side_1_position"]);
        int32_t side_2_position = n(jc_item["side_2_position"]);
        int32_t side_1_strand = n(jc_item["side_1_strand"]);
        int32_t side_2_strand = n(jc_item["side_2_strand"]);
        
        if (side_2_position < side_1_position) {
          swap(side_1_position, side_2_position);
          swap(side_1_strand, side_2_strand);
        }
        
				if (
            (side_1_position == n(mut["position"])-1)
            && (side_1_strand == -1)
            && (side_2_position == n(mut["position"])+n(mut["size"]))
            && (side_2_strand == +1)
            )
				{
          
          //it's possible that one or both sides are in repeat elements
          cFeatureLocation* r1_pointer = within_repeat(jc_item["side_1_seq_id"], side_1_position);
          cFeatureLocation* r2_pointer = within_repeat(jc_item["side_2_seq_id"], side_2_position);
          
          // one repeat cases where the end matches up exactly
          if (r1_pointer) 
          {
            // must match up to an end of the repeat
            if (side_1_position == static_cast<int32_t>(r1_pointer->get_start_1())
                || (side_1_position == static_cast<int32_t>(r1_pointer->get_end_1())))
            {
              mut["mediated"] = (*r1_pointer->get_feature())["name"];
            }
          }
          // if it didn't match, then check possibility of a second repeat
          if (!mut.count("mediated") && r2_pointer) 
          {
            // must match up to an end of the repeat
            if ((side_2_position == static_cast<int32_t>(r2_pointer->get_start_1())
                 || (side_2_position) == static_cast<int32_t>(r2_pointer->get_end_1())))
            {
              mut["mediated"] = (*r2_pointer->get_feature())["name"];
            }
          }    
          
					mut._evidence.push_back(jc_item._id);
          
          // If there is unique sequence in the junction, then it is actually a SUB
          JunctionInfo ji(jc_item["key"]);
          if (ji.unique_read_sequence.size() > 0) {
            mut._type = SUB;
            mut[NEW_SEQ] = ji.unique_read_sequence;
          }
          
          jc_it = jc.erase(jc_it);
          
					gd.add(mut);
          mc_it = mc.erase(mc_it); // iterator is now past the erased element
          mc_it--;                //We just removed the current jc, do not iterate.
          
          if (verbose)
            cout << "**** Junction precisely matching deletion boundary found ****\n";
          done = true;
          
          break; // out of jc item loop;
				}
			}
      if (done) continue; // to next mc item
      
      
      cFeatureLocation* r1_pointer = within_repeat(mut["seq_id"], n(mut["position"]));
      cFeatureLocation* r2_pointer = within_repeat(mut["seq_id"], n(mut["position"]) + n(mut["size"]));
      
			///
			// (2) there is no junction, but both ends of the deletion are in different copies of the same repeat sequence
			///
      ///
      if (verbose)
        cout << "(2) Checking for MC between two repeats." << endl;
      
			// Then we will adjust the coordinates to remove...
			if  (
           (r1_pointer != r2_pointer) 
           && (r1_pointer != NULL) 
           && (r2_pointer != NULL) 
           && ((*r1_pointer->get_feature())["name"] == (*r2_pointer->get_feature())["name"])
           )
			{
				cFeatureLocation& r1 = *r1_pointer, r2 = *r2_pointer;
        
				// there may be more evidence that one or the other is deleted...
				int32_t r1_overlap_end = n(mc_item["start"]) + n(mc_item["start_range"]);
				if (r1_overlap_end > r1.get_end_1())
					r1_overlap_end = r1.get_end_1();
				int32_t r1_overlap = r1_overlap_end - n(mc_item["start"]) + 1;
        
				int32_t r2_overlap_start = n(mc_item["end"]) - n(mc_item["end_range"]);
				if (r2_overlap_start < r1.get_start_1())
					r2_overlap_start = r2.get_start_1();
				int32_t r2_overlap = n(mc_item["end"]) - r2_overlap_start + 1;
        
				// it may be really close...defined by read length of genome in which case
				uint32_t slop_distance = read_length_max;
        
				// prefer to delete the second copy
				if ( (static_cast<uint32_t>(abs(r1_overlap - r2_overlap)) <= slop_distance) || (r2_overlap > r1_overlap ))
				{
					mut["position"] = s(r1.get_end_1() + 1);
					mut["size"] = s(r2.get_end_1() - r1.get_end_1());
				}
				else // delete the first copy
				{
					mut["position"] = s(r1.get_start_1());
          mut["size"] = s(r2.get_start_1() - r1.get_start_1());
				}
        
        // @JEB 2014-01-07
        // It's possible for this to be in the SAME copy of the element,
        // in which case the deletion size here is zero bases
        if (n(mut["size"]) > 0) {
         
          // @JEB 2016-03-26
          // Its possible that NOW after this adjustment the MC ends exactly match a JC item
          // This is the loop to find that JC item!
          for(diff_entry_list_t::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) //JC
          {
            cDiffEntry& jc_item = **jc_it;
            
            if (jc_item["side_1_seq_id"] != mut["seq_id"] || jc_item["side_2_seq_id"] != mut["seq_id"])  {
              continue;
            }
            
            // Due to overlap resolution, this can change from the time when we sorted this way... *Ugh* ..
            // so we fix it here so that the lower coordinate part of the junction is first.
            
            int32_t side_1_position = n(jc_item["side_1_position"]);
            int32_t side_2_position = n(jc_item["side_2_position"]);
            int32_t side_1_strand = n(jc_item["side_1_strand"]);
            int32_t side_2_strand = n(jc_item["side_2_strand"]);
            
            if (side_2_position < side_1_position) {
              swap(side_1_position, side_2_position);
              swap(side_1_strand, side_2_strand);
            }
            
            if (
                (side_1_position == n(mut["position"])-1)
                && (side_1_strand == -1)
                && (side_2_position == n(mut["position"])+n(mut["size"]))
                && (side_2_strand == +1)
                )
            {
              mut._evidence.push_back(jc_item._id);
              
              // If there is unique sequence in the junction, then it is actually a SUB
              JunctionInfo ji(jc_item["key"]);
              if (ji.unique_read_sequence.size() > 0) {
                mut._type = SUB;
                mut[NEW_SEQ] = ji.unique_read_sequence;
              }
              
              // do not re-use the junction
              jc_it = jc.erase(jc_it);
              
              if (verbose)
                cout << "**** Junction matching MC deletion between repeats found ****\n";
              done = true;
              
              break; // out of jc item loop;
            }
          } // End of JC loop
          
          // remember the name of the element
          mut["between"] = (*r1.get_feature())["name"];
          gd.add(mut);
          
          mc_it = mc.erase(mc_it); // iterator is now past the erased element
          mc_it--;                //We just removed the current jc, do not iterate.
          
          if (verbose)
            cout << "**** Ends of junction in copies of same repeat element ****\n";
          continue; // to next mc_item
        }
			}
      
      // Both MC ends inside IS elements but case (2) did not match
      if (r1_pointer != NULL && r2_pointer != NULL)
        continue; // to next mc_item

      ///
      // (2b) Both boundaries are in unique sequence, but each is connected by a JC to the
      //      same IS element — predict a MOB (IS insertion, duplication_size=0) at the left
      //      boundary, plus an IS-mediated DEL for the MC region.
      ///
      if (r1_pointer == NULL && r2_pointer == NULL) {
        if (verbose)
          cout << "(2b) Checking for MC flanked by two IS-element junctions." << endl;

        int32_t needed_left  = n(mut["position"]) - 1;              // MC.start - 1
        int32_t needed_right = n(mut["position"]) + n(mut["size"]); // MC.end   + 1

        cDiffEntry* jc_left  = NULL;
        cDiffEntry* jc_right = NULL;
        diff_entry_list_t::iterator jc_left_it, jc_right_it;

        for (diff_entry_list_t::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) {
          cDiffEntry& j = **jc_it;
          if (!j.entry_exists("_is_interval")) continue;
          if (j[j["_unique_interval"] + "_seq_id"] != mut["seq_id"]) continue;

          int32_t uniq_pos    = n(j[j["_unique_interval"] + "_position"]);
          int32_t uniq_strand = n(j[j["_unique_interval"] + "_strand"]);

          if (jc_left == NULL && uniq_pos == needed_left && uniq_strand == -1) {
            jc_left    = &j;
            jc_left_it = jc_it;
          } else if (jc_right == NULL && uniq_pos == needed_right && uniq_strand == +1) {
            jc_right    = &j;
            jc_right_it = jc_it;
          }
          if (jc_left != NULL && jc_right != NULL) break;
        }

        if (jc_left != NULL && jc_right != NULL) {
          string is_name  = (*jc_left) ["_" + (*jc_left) ["_is_interval"] + "_is_name"];
          string is_name2 = (*jc_right)["_" + (*jc_right)["_is_interval"] + "_is_name"];

          if (is_name == is_name2) {
            // MOB strand from JC_left (unique_interval_strand == -1)
            int32_t is_strand = -(
              n((*jc_left)[(*jc_left)["_is_interval"] + "_strand"])
              * n((*jc_left)["_" + (*jc_left)["_is_interval"] + "_is_strand"])
              * n((*jc_left)[(*jc_left)["_unique_interval"] + "_strand"])
            );

            // MOB: IS inserted at left boundary with no TSD; evidence = JC_left only
            cDiffEntry mut_mob;
            mut_mob._type = MOB;
            mut_mob._evidence.push_back(jc_left->_id);
            mut_mob
              ("seq_id",                         mut["seq_id"])
              ("repeat_name",                    is_name)
              ("strand",                         s(is_strand))
              ("duplication_size",               "0")
              ("indeterminate_duplication_size", "1")
              ("position",                       s(n(mut["position"]) - 1));

            // In polymorphism mode every mutation must carry a frequency (checked at the end of
            // predict()). It is 1 here, like the deletion this case predicts alongside: the case only
            // fires because there is missing coverage, and MC is not predicted unless the region is
            // fully absent -- so the insertion that replaced it is fixed, whatever the junction's own
            // fitted frequency says.
            if (settings.polymorphism_prediction) mut_mob[FREQUENCY] = "1";

            // DEL: evidence = MC (already in mut._evidence) + JC_right
            mut["mediated"] = is_name;
            mut._evidence.push_back(jc_right->_id);

            // Erase later iterator first to keep earlier one valid (std::list guarantee)
            jc.erase(jc_right_it);
            jc.erase(jc_left_it);

            gd.add(mut_mob);
            gd.add(mut);

            mc_it = mc.erase(mc_it);
            mc_it--;

            if (verbose)
              cout << "**** Double IS-mediated deletion: MOB + DEL predicted ****\n";
            continue;
          }
        }
        continue; // case (3) requires exactly one non-NULL r-pointer
      }

			///
			// (3) there is a junction between unique sequence and a repeat element
			///
      ///
      if (verbose)
        cout << "(3) Checking for MC between unique sequence and repeat." << endl;
     
			cFeatureLocation& r = (r1_pointer != NULL) ? *r1_pointer : *r2_pointer;
			int32_t redundant_deletion_side = (r1_pointer != NULL) ? -1 : +1;
			int32_t unique_deletion_strand = -redundant_deletion_side;
			int32_t needed_coord = (r1_pointer != NULL)
      ? n(mut["position"]) + n(mut["size"])
      : n(mut["position"]) - 1;
      
			for(diff_entry_list_t::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) //JUNCTION
			{
				cDiffEntry& j = **jc_it;
        
        if (verbose)
          cout << j << endl;
        
				if (!j.entry_exists("_is_interval")) continue;
        
				if (verbose)
					cout << "Check 1: " << j[j["_unique_interval"] + "_seq_id"] << " ne " << mut["seq_id"] << endl;
				if (j[j["_unique_interval"] + "_seq_id"] != mut["seq_id"])
					continue;
				if (verbose)
					cout << "Pass 1" << endl;
        
				// check type of IS
				if (verbose)
					cout << "Check 2: " << (*r.get_feature())["name"] << " ne " << j["_" + j["_is_interval"] + "_is_name"] << endl;
				if ((*r.get_feature())["name"] != j["_" + j["_is_interval"] + "_is_name"])
					continue;
        
				if (verbose) cout << "Pass 2" << endl;
        
				// check that the unique side matches coordinate
				if (verbose)
					cout << "Check 3: " << j[j["_unique_interval"] + "_position"] << " != " << needed_coord << endl;
				if (n(j[j["_unique_interval"] + "_position"]) != needed_coord)
					continue;
				if (verbose)
					cout << "Pass 3" << endl;
        
        
				// check that IS is on the right strand
				if (verbose)
					cout << "Check 4: " << redundant_deletion_side << " * " << to_string(r.get_strand()) << " != " << j[j["_is_interval"] + "_strand"] << " * " << j["_" + j["_is_interval"] + "_is_strand"] << endl;
				if ( (redundant_deletion_side * r.get_strand()) != (n(j[j["_is_interval"] + "_strand"]) * n(j["_" + j["_is_interval"] + "_is_strand"])) )
					continue;
				if (verbose)
					cout << "Pass 4" << endl;
        
				// check that the unique side is on the right strand
				if (verbose)
					cout << "Check 5: " << unique_deletion_strand << " != " << j[j["_unique_interval"] + "_strand"] << endl;
				if ( unique_deletion_strand != n(j[j["_unique_interval"] + "_strand"]) )
					continue;
				if (verbose)
					cout << "Pass 5" << endl;
        
        
				// need to adjust the non-unique coords // mut has the MC coordinates at this point
				if (redundant_deletion_side == -1)
				{          
					uint32_t move_dist = r.get_end_1() + 1 - n(mut["position"]) + n(j["_" + j["_is_interval"] + "_is_distance"]);
					mut["position"] = s(n(mut["position"]) + move_dist);
					mut["size"] = s(n(mut["size"]) - move_dist);
				}
				else
				{
					int32_t move_dist = (n(mut["position"]) + n(mut["size"]) - 1) - (r.get_start_1() - 1) + n(j["_" + j["_is_interval"] + "_is_distance"]);
					mut["size"] = s(n(mut["size"]) - move_dist);
				}
        
        // Don't predict zero length deletions!
        if (n(mut["size"]) <= 0)
          continue;
        
				// OK, we're good!
				mut["mediated"] = (*r.get_feature())["name"];
				mut._evidence.push_back(j._id);
				jc.erase(jc_it);
				gd.add(mut);
        
        mc_it = mc.erase(mc_it); // iterator is now past the erased element
        mc_it--;                //We just removed the current jc, do not iterate.
        
        if (verbose)
          cout << "**** Junction with repeat element corresponding to deletion boundaries found ****\n";
        
        break; // done looking at jc_items
			}
		}

    
    // DEL prediction (separate case):
    // (4) the reference is circular, there is missing coverage at one or both ends,
    //     AND there is junction connecting those ends
    
    if (verbose)
      cout << "(4) Checking for MC crossing circular genome ends." << endl;
    
    for(diff_entry_list_t::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) //JC
    {
      cDiffEntry& jc_item = **jc_it;
      
      if (verbose)
        cout << jc_item << endl;
      
      // They have to be for the same seq_id
      if (jc_item["side_1_seq_id"] != jc_item["side_2_seq_id"] )  {
        continue;
      }
      string& seq_id = jc_item["side_1_seq_id"];
      int32_t seq_length = static_cast<int32_t>(ref_seq_info[seq_id].get_sequence_length());
      
      // This seq_id must be circular
      if (!ref_seq_info[seq_id].is_circular())
        continue;
      
      // They have to connect across the origin of the sequence
      // ---> assumes position_1 is lower
      
      int32_t side_1_position = n(jc_item[SIDE_1_POSITION]);
      int32_t side_2_position = n(jc_item[SIDE_2_POSITION]);
      int32_t side_1_strand = n(jc_item[SIDE_1_STRAND]);
      int32_t side_2_strand = n(jc_item[SIDE_2_STRAND]);
      
      // There should not be any overlap
      int32_t overlap = n(jc_item[OVERLAP]);
      if (overlap > 0) continue;
      
      //This is not necessarily true
      //ASSERT(side_1_position <= side_2_position, "Junction has side_1_position > side_2_position\n" + jc_item.as_string());
      
      if (! ((side_1_strand == +1) &&  (side_2_strand == -1)) ) {
        continue;
      }
      
      // Now we have to find matching MC items
      cDiffEntry* start_seq_mc(NULL);
      cDiffEntry* end_seq_mc(NULL);
      
      for(diff_entry_list_t::iterator mc_it = mc.begin(); mc_it != mc.end(); mc_it++)
      {
        cDiffEntry& mc_item = **mc_it;
        
        if (verbose)
          cout << mc_item << endl;
        
        if (mc_item.entry_exists("reject"))
          continue;
        
        if (mc_item[SEQ_ID] != seq_id)
          continue;
        
        if (n(mc_item[START]) - n(mc_item[START_RANGE]) <= 1)
          start_seq_mc = &mc_item;
        
        if (n(mc_item[END]) + n(mc_item[END_RANGE]) >= seq_length)
          end_seq_mc = &mc_item;
      }
      
      
      if (verbose) {
        if (start_seq_mc)
          cout << "START_SEQ_MC" << endl << *start_seq_mc << endl;
        if (end_seq_mc)
          cout << "END_SEQ_MC" << endl << *end_seq_mc << endl;
      }
      
      // We have to have found MC on both ends unless the junctions are flush to the ends
      if (!start_seq_mc && (side_1_position != 1))
          continue;
      if (!end_seq_mc && (side_2_position != seq_length))
          continue;
      
      // Did we find both missing coverage pieces and they match up with the junction?
      // Then we get to create a mutation! It will extend past the end of the fragment.

      // Note: some additional slop could be added here
      int32_t start_seq_mc_end_start_range_1(0), start_seq_mc_end_end_range_1(0);
      if (start_seq_mc) {
        start_seq_mc_end_start_range_1 = n((*start_seq_mc)[END]);
        start_seq_mc_end_end_range_1 = start_seq_mc_end_start_range_1 + n((*start_seq_mc)[END_RANGE]);
      }
      
      int32_t end_seq_mc_start_start_range_1(0), end_seq_mc_start_end_range_1(0);
      if (end_seq_mc) {
        end_seq_mc_start_end_range_1 = n((*end_seq_mc)[START]);
        end_seq_mc_start_start_range_1 = end_seq_mc_start_end_range_1 - n((*end_seq_mc)[START_RANGE]);
      }
      
      if (verbose)
        cout << "Passed Checks" << endl;
      
      if (
          (!start_seq_mc || ((side_1_position-1 >= start_seq_mc_end_start_range_1) && (side_1_position-1 <= start_seq_mc_end_end_range_1)) )
       && (!end_seq_mc || ((side_2_position+1 >= end_seq_mc_start_start_range_1) && (side_2_position+1 <= end_seq_mc_start_end_range_1)) )
          )
      {
        cDiffEntry mut;
        mut._evidence = make_vector<string>(jc_item._id);
        
        int32_t size = (side_1_position - 1) + (seq_length - side_2_position);
        int32_t position = side_2_position+1;
        // Necessary for case where the deletion is only at the start
        if (position > seq_length) position -= seq_length;
        
        // If there is unique read sequence, then it is a substitution
        if ( jc_item.entry_exists(UNIQUE_READ_SEQUENCE)) {
          
          if (verbose)
            cout << "Predicting SUB" << endl;
          
          mut._type = SUB;

          // If the size is zero then it wil be marked as a normal circular read later,
          // in predictJCtoINSorSUBorDEL
          if (size == 0) continue;
          
          mut
          (SEQ_ID, seq_id)
          (POSITION, s(position))
          (SIZE, s(size))
          (NEW_SEQ, jc_item[UNIQUE_READ_SEQUENCE])
          ;
          
        // Otherwise it is a deletion
        } else {
          
          if (verbose)
            cout << "Predicting DEL" << endl;
          
          mut._type = DEL;

          // If the size is zero then it wil be marked as a normal circular read later,
          // in predictJCtoINSorSUBorDEL
          if (size == 0) continue;
          
          mut
          (SEQ_ID, seq_id)
          (POSITION, s(position))
          (SIZE, s(size));
          ;
        }
        
        if (settings.polymorphism_prediction) {
          mut[FREQUENCY] = "1";
        }
        
        if (start_seq_mc) mut._evidence.push_back(start_seq_mc->_id);
        if (end_seq_mc) mut._evidence.push_back(end_seq_mc->_id);
        
        gd.add(mut);
        
        // Not really necessary to delete
        
        if (verbose)
          cout << "Deleting JC" << endl;
        
        jc_it = jc.erase(jc_it);
        jc_it--;
      }
    }
    
  }
  
  void MutationPredictor::predictJCplusJCtoMOB(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc)
  {
    (void)summary;
    bool verbose = false;
    
    // Built through the shared factory so this function's counters cannot drift out of step with
    // the ones in assign_junction_read_counts(). That mattered: this function re-counts entries
    // that were already counted there, and it runs later, so a counter configured differently
    // here silently overwrites the earlier results -- with no error and no failing assertion.
    counted_ptr<junction_read_counter> reference_jrc, junction_jrc;
    make_junction_read_counters(settings, reference_jrc, junction_jrc);

    
    // Prepare the lists
    for(diff_entry_list_t::iterator it = jc.begin(); it != jc.end(); it++) //JC
		{
			cDiffEntry& j = **it;
      
			// Junction isn't near an IS. Move on.
			if (!j.entry_exists("_is_interval")) continue;
      
			// There is no overlap to correct. Move on.
			if (n(j["overlap"]) <= 0) continue;
      
			// The following code implies n(j["overlap"]) > 0
      
			/// first, adjust the repetitive sequence boundary to get as close to the IS as possible
			int32_t is_interval_position = n(j[j["_is_interval"] + "_position"]);
			uint32_t move_dist = abs(is_interval_position - n(j[j["_is_interval"] + "_is_" + j["_is_interval_closest_side_key"]])); //  + "_position" ?
			uint32_t overlap = n(j["overlap"]);
			if (move_dist > overlap) move_dist = overlap;
      
			is_interval_position += n(j[j["_is_interval"] + "_strand"]) * move_dist;
			j[j["_is_interval"] + "_position"] = s(is_interval_position);
			overlap -= move_dist;
      
			/// second, adjust the unique sequence side with any remaining overlap
			int32_t unique_interval_position = n(j[j["_unique_interval"] + "_position"]);
			unique_interval_position += n(j[j["_unique_interval"] + "_strand"]) * overlap;
			j[j["_unique_interval"] + "_position"] = s(unique_interval_position);
      
			j["overlap"] = "0";
		}
    
    // This sorts by seq_id matched then by position of unique coordinate side
		jc.sort(MutationPredictor::sort_by_hybrid);
    
    // @JEB: Notice that we don't advance the iterator here. This happens at the end of the JC1 loop 
    // (which is never skipped) because we sometimes delete entries in ways that are hard for list iterators to deal with
    for(diff_entry_list_t::iterator jc1_it = jc.begin(); jc1_it != jc.end(); ) //JC1
		{
      bool jc1_erased = false;
			cDiffEntry& j1 = **jc1_it;
     
      if (verbose) {
        cout << j1.as_string() << endl;
      }
      
			// Compile a list of the next possibilities within a certain length of bases
      vector<diff_entry_ptr_t> j2_list;
      
      // start looking at the next JC entry
      diff_entry_list_t::iterator jc2_it = jc1_it;
      for(jc2_it++ ; jc2_it != jc.end(); jc2_it++) //JC2
			{
				cDiffEntry& j2 = **jc2_it;
        
				// must be close together in real coords
				if ( (j1[j1["_unique_interval"] + "_seq_id"] != j2[j2["_unique_interval"] + "_seq_id"])
            || (abs(n(j1[j1["_unique_interval"] + "_position"]) - n(j2[j2["_unique_interval"] + "_position"])) > 20 ) )
					break;
        
				if ( (!j1.entry_exists("_is_interval") || !j2.entry_exists("_is_interval"))
            || (j1["_" + j1["_is_interval"] + "_is_name"] != j2["_" + j2["_is_interval"] + "_is_name"]) )
					continue;
        
				j2_list.push_back(*jc2_it);
			}
      if (verbose)
        cout << "Size of J2 list: " << j2_list.size() << endl;
      
			//sort the $j2_list by reject reason and score
      
			sort(j2_list.begin(), j2_list.end(), sort_by_reject_score);
      
			// We need to go through all with the same coordinate (or within a certain coordinate stretch)
			// because sometimes a failed junction will be in between the successful junctions
      for(size_t i=0; i<j2_list.size(); i++)
			{        
				cDiffEntry& j2 = *(j2_list[i]);
        
        if (verbose) 
        {
          cout << "Sorted: == J1 ==" << endl;
          for(map<string,string>::iterator it=j1.begin(); it!=j1.end(); it++)
          {
            cout << it->first << " = " << it->second << endl; 
          }
          cout << "Sorted: == J2 ==" << endl;
          for(map<string,string>::iterator it=j2.begin(); it!=j2.end(); it++)
          {
            cout << it->first << " = " << it->second << endl; 
          }
        }
        
        // If both are rejects, then don't keep (one _can_ be a reject!)
        if (j1.is_rejected_and_not_user_defined() && j2.is_rejected_and_not_user_defined())
        {
          if (verbose)
          {
            cout << "Both are rejected junctions. No prediction." << endl;
          }
          continue;
        }
        
        
				// positive overlap should be resolved by now
				assert(n(j1["overlap"]) <= 0);
				assert(n(j2["overlap"]) <= 0);
        
				// the first unique coords are going into the IS element
				int32_t uc1_strand = n(j1[j1["_unique_interval"] + "_strand"]);
				int32_t uc2_strand = n(j2[j2["_unique_interval"] + "_strand"]);
				if (uc1_strand != -uc2_strand) continue;
        
				// What strand is the IS on relative to the top strand of the genome
				int32_t is1_strand = - (n(j1[j1["_is_interval"] + "_strand"]) * n(j1["_" + j1["_is_interval"] + "_is_strand"]) * n(j1[j1["_unique_interval"] + "_strand"]));
				int32_t is2_strand = - (n(j2[j2["_is_interval"] + "_strand"]) * n(j2["_" + j2["_is_interval"] + "_is_strand"]) * n(j2[j2["_unique_interval"] + "_strand"]));

				// Create the mutation, with evidence
        
        if (verbose) {
          cout << "Predicting mutation from junctions:" << endl << j1 << endl << j2 << endl;
        }
        
				cDiffEntry mut;
        mut._type = MOB;
				mut._evidence.push_back(j1._id);
				mut._evidence.push_back(j2._id);
				mut
        ("seq_id", j1[j1["_unique_interval"] + "_seq_id"])
				;
				mut["_start"] = (uc1_strand == -1) ? j2[j2["_unique_interval"] + "_position"] : j1[j1["_unique_interval"] + "_position"];
				mut["_end"] = (uc1_strand == -1) ? j1[j1["_unique_interval"] + "_position"] : j2[j2["_unique_interval"] + "_position"];
				mut["repeat_name"] = j1["_" + j1["_is_interval"] + "_is_name"];
        
				mut["position"] = mut["_start"]; // - 1; //position is the first duplicated base...
				mut["duplication_size"] = s(n(mut["_end"]) - n(mut["_start"]) + 1);
        
				// ok, we're actually missing a base of the reference...
				if (n(mut["duplication_size"]) < 0)
				{
					mut["position"] = s(n(mut["position"]) + n(mut["duplication_size"]));
				}
        // @JEB 2014-01-05
        // Special case of no duplicated bases, we shift the position by one backward
        // so that the MOB is inserted AFTER this position (rather than BEFORE).
        else if (n(mut["duplication_size"]) == 0)
        {
          mut["position"] = s(n(mut["_start"])-1);
        }
        
				// get any unique junction sequence
        string j1_unique_read_sequence;
        if (j1.entry_exists("key")) {
          JunctionInfo j1i(j1["key"]);
           j1_unique_read_sequence = j1i.unique_read_sequence;
        }
        
        string j2_unique_read_sequence;
        if (j2.entry_exists("key")) {
          JunctionInfo j2i(j2["key"]);
          j2_unique_read_sequence = j2i.unique_read_sequence;
        }
        
				// _gap_left and _gap_right also refer to the top strand of the genome
        
				mut["_ins_start"] = "";
				mut["_ins_end"] = "";
        
				mut["_del_start"] = "0";
				mut["_del_end"] = "0";
        
        uint32_t start_1 = 0, end_1 = 0, pos_1 = 0;
        
        
        ///////////////////////////////////////////////////////
        // Figuring out bases *inserted* next to the IS element
        ///////////////////////////////////////////////////////
        
        // j1_not_flush_seq is for an insertion
        
				string j1_not_flush_seq = "";
				if (n(j1[j1["_is_interval"] + "_strand"]) == -1)
				{
					mut["_gap_left"] = s(n(j1[j1["_is_interval"] + "_position"]) - n(j1["_" + j1["_is_interval"] + "_is_end"]));
					if (n(mut["_gap_left"]) > 0)
          {
            start_1 = n(j1["_" + j1["_is_interval"] + "_is_end"]) + 1;
            end_1   = n(j1[j1["_is_interval"] + "_position"]);
						j1_not_flush_seq = ref_seq_info.get_sequence_1 (
                                                            j1[j1["_is_interval"] + "_seq_id"],
                                                            start_1,
                                                            end_1
                                                            );
					}
				}
				else
				{
					mut["_gap_left"] = s(n(j1["_" + j1["_is_interval"] + "_is_start"]) - n(j1[j1["_is_interval"] + "_position"]));          
					if (n(mut["_gap_left"]) > 0)
					{
            start_1 = n(j1[j1["_is_interval"] + "_position"]);
            end_1   = n(j1["_" + j1["_is_interval"] + "_is_start"]) - 1;
            j1_not_flush_seq = ref_seq_info.get_sequence_1 (
                                                            j1[j1["_is_interval"] + "_seq_id"],
                                                            start_1,
                                                            end_1
                                                            );
					}
				}
        
        if (n(mut["_gap_left"]) < 0)
				{
					mut["_del_start"] = s(abs(n(mut["_gap_left"])));
				}
        
        
        if (verbose)
          cout << "J1 NF:" << j1_not_flush_seq << " U:" << j1_unique_read_sequence << endl;
        
        if (n(j1["_" + j1["_is_interval"] + "_read_side"]) != n(j1[j1["_is_interval"] + "_strand"]))
        {
          j1_not_flush_seq = reverse_complement(j1_not_flush_seq);
        }
        
        if (n(j1["_" + j1["_is_interval"] + "_read_side"]) == -1)
        {
          mut["_ins_start"] = j1_not_flush_seq + j1_unique_read_sequence;
        }
        else
        {
          mut["_ins_start"] = j1_unique_read_sequence + j1_not_flush_seq;
        }
        
        
				string j2_not_flush_seq = "";
				if (n(j2[j2["_is_interval"] + "_strand"]) == -1)
				{
					mut["_gap_right"] = s(n(j2[j2["_is_interval"] + "_position"]) - n(j2["_" + j2["_is_interval"] + "_is_end"]));
					if (n(mut["_gap_right"]) > 0)
					{
            start_1 = n(j2["_" + j2["_is_interval"] + "_is_end"]) + 1;
            end_1   = n(j2[j2["_is_interval"] + "_position"]);
            j2_not_flush_seq = ref_seq_info.get_sequence_1 (
                                                            j1[j2["_is_interval"] + "_seq_id"],
                                                            start_1,
                                                            end_1
                                                            );
					}
				}
				else
				{
					mut["_gap_right"] = s(n(j2["_" + j2["_is_interval"] + "_is_start"]) - n(j2[j2["_is_interval"] + "_position"]));          
					if (n(mut["_gap_right"]) > 0)
					{
            start_1 = n(j2[j2["_is_interval"] + "_position"]);
            end_1   = n(j2["_" + j2["_is_interval"] + "_is_start"]) - 1;
            j2_not_flush_seq = ref_seq_info.get_sequence_1 (
                                                            j1[j2["_is_interval"] + "_seq_id"],
                                                            start_1,
                                                            end_1
                                                            );
					}
				}
        
        if (n(mut["_gap_right"]) < 0)
				{
					mut["_del_end"] = s(abs(n(mut["_gap_right"])));
				}
        
        if (verbose)
          cout << "J2 NF:" << j2_not_flush_seq << " U:" << j2_unique_read_sequence << endl;
        
        if ( n(j2["_" + j2["_is_interval"] + "_read_side"]) * n(j2[j2["_is_interval"] + "_strand"]) == -1)
        {
          j2_not_flush_seq = reverse_complement(j2_not_flush_seq);
        }
        
        if (n(j2["_" + j2["_is_interval"] + "_read_side"]) == -1)
        {
          mut["_ins_end"] = j2_not_flush_seq + j2_unique_read_sequence;
        }
        else
        {
          mut["_ins_end"] = j2_unique_read_sequence + j2_not_flush_seq;
        }
        
        
				// At this point any added junction sequences are on the strand as you would see them in the alignment.
				// we may need to reverse complement. Note that we never swap sides because these changes are with respect to
        // the reference and do not change if the mobile element is inserted in the reverse orientation.
        
				if (verbose)
					cout << mut["_gap_left"] << " :: " << mut["_gap_right"] << endl;
        
				if ( n(j1[j1["_unique_interval"] + "_strand"]) != n(j1["_" + j1["_unique_interval"] + "_read_side"]))
				{
					if (verbose) cout << "RC left" << endl;
					mut["_ins_start"] = reverse_complement(mut["_ins_start"]);
				}
        
				if ( n(j2[j2["_unique_interval"] + "_strand"]) != n(j2["_" + j2["_unique_interval"] + "_read_side"]))
				{
					if (verbose) cout << "RC right" << endl;
					mut["_ins_end"] = reverse_complement(mut["_ins_end"]);
				}
        
				//// Check for ambiguous insertion direction!
				// Sometimes a strand will be assigned just because there is a 50-50 chance of getting the correct sides of the IS.
				// We need to actually check the sequence on each side of the repeat element on the end as far in as the maximum overlap on that side.
				// Use:			{max_left"] {max_right"]
				// Retrieve sequence on unique side and compare to sequence on the other side of a repeat element
        
        uint32_t j1_not_flush_length = j1_not_flush_seq.size();
				uint32_t j2_not_flush_length = j2_not_flush_seq.size();
				uint32_t max_not_flush_length = max(j1_not_flush_length, j2_not_flush_length);
        
        if (verbose) {
					cout << "J1 not flush length: " << j1_not_flush_length << endl;
					cout << "J2 not flush length: " << j2_not_flush_length << endl;
					cout << "Max not flush length: " << max_not_flush_length << endl;
				}
        
        int32_t j1_is_overlap_length(0);
        if (j1.entry_exists("_is_interval") && j1.entry_exists("_" + j1["_is_interval"] + "_read_side") && j1.entry_exists("max_left") && j1.entry_exists("max_right")) {
          j1_is_overlap_length = (n(j1["_" + j1["_is_interval"] + "_read_side"]) == -1) ? n(j1["max_left"]) : n(j1["max_right"]);
        }

        int32_t j2_is_overlap_length(0);
        if (j2.entry_exists("_is_interval") && j2.entry_exists("_" + j2["_is_interval"] + "_read_side") && j2.entry_exists("max_left") && j2.entry_exists("max_right")) {
          j2_is_overlap_length = (n(j2["_" + j2["_is_interval"] + "_read_side"]) == -1) ? n(j2["max_left"]) : n(j2["max_right"]);
        }
      
				if (verbose) {
					cout << "J1 IS overlap length: " << j1_is_overlap_length << endl;
					cout << "J2 IS overlap length: " << j2_is_overlap_length << endl;
				}
        
        
        bool j1_is_ambiguous = false;
          
        // what are the actual sequences of this length at the end of the IS elements?
        start_1 = n(j1["_" + j1["_is_interval"] + "_is_start"]);
        end_1   = start_1 + j1_is_overlap_length - 1;
        string j1_left_is_sequence;
        if (end_1 >= start_1) {
          j1_left_is_sequence = ref_seq_info.get_sequence_1 (
                                                                  j1[j1["_is_interval"] + "_seq_id"],
                                                                  start_1,
                                                                  end_1
                                                                  );
        }
        
        end_1   = n(j1["_" + j1["_is_interval"] + "_is_end"]);
        start_1 = end_1 - j1_is_overlap_length - 1;
        string j1_right_is_sequence;
        if (end_1 >= start_1) {
          j1_right_is_sequence = ref_seq_info.get_sequence_1 (
                                                                   j1[j1["_is_interval"] + "_seq_id"],
                                                                   start_1,
                                                                   end_1
                                                                   );
        }
        j1_right_is_sequence = reverse_complement(j1_right_is_sequence);
        
        if (verbose) {
          cout << "J1 LEFT : " << j1_left_is_sequence << endl;
          cout << "J1 RIGHT: " << j1_right_is_sequence << endl;
        }
        
        // believe the direction if the sequences are different
        j1_is_ambiguous = (j1_left_is_sequence == j1_right_is_sequence);

        
        bool j2_is_ambiguous = false;
          
        start_1 = n(j2["_" + j2["_is_interval"] + "_is_start"]);
        end_1   = start_1 +j2_is_overlap_length - 1;
        string j2_left_is_sequence;
        if (end_1 >= start_1) {
          j2_left_is_sequence = ref_seq_info.get_sequence_1 (
                                                                  j2[j2["_is_interval"] + "_seq_id"],
                                                                  start_1,
                                                                  end_1
                                                                  );
        }
        
        end_1   = n(j2["_" + j2["_is_interval"] + "_is_end"]);
        start_1 = end_1 - j2_is_overlap_length - 1;
        string j2_right_is_sequence;
        if (end_1 >= start_1) {
          j2_right_is_sequence = ref_seq_info.get_sequence_1 (
                                                                   j2[j2["_is_interval"] + "_seq_id"],
                                                                   start_1,
                                                                   end_1
                                                                   );
        }
        j2_right_is_sequence = reverse_complement(j2_right_is_sequence);
        
        // believe the direction if the sequences are different
        j2_is_ambiguous = (j2_left_is_sequence == j2_right_is_sequence);
        
        if (verbose) {
          cout << "J2 LEFT : " << j2_left_is_sequence << endl;
          cout << "J2 RIGHT: " << j2_right_is_sequence << endl;
          
          //cout << "J1 IS matched length " << j1_is_overlap_length << ": " << j1_is_seq_matched << endl;
          //cout << "J2 IS matched length " << j2_is_overlap_length << ": " << j2_is_seq_matched << endl;
        }
				
				// if the matched IS element sequences are the same then the direction is AMBIGUOUS
				int32_t is_strand = 0;
				if (j1_is_ambiguous && j2_is_ambiguous)
				{
					if (verbose) cout << "AMBIGUOUS strand for mobile element insertion" << endl;
				}
				else if (j1_is_ambiguous)
				{
					is_strand = is2_strand;
				}
				else if (j2_is_ambiguous)
				{
					is_strand = is1_strand;
				}
				else // neither is ambiguous and the strands don't agree (this is not a simple IS insertion)
				{
					if (is1_strand != is2_strand) continue;
          is_strand = is1_strand;
				}
				mut["strand"] = s(is_strand);
        
				////
				//// We are still not checking for a case where one junction side extends far enough to uniquely define the
				//// side of the IS, but the other side does not (giving it the wrong strand).
				////
        
				// Finally, do this AFTER checking for the IS-matched sequences...
				// side_1 of the junction may be the left side, rather than the right side of the insertion, if so we need
        // to swap coords (but not reverse complement anything
				if (uc1_strand == +1)
				{
					if (verbose) cout << "reverse right and left" << endl;
          
          swap(mut["_ins_start"], mut["_ins_end"]);
          swap(mut["_del_start"], mut["_del_end"]);          
				}
        
        
				// only transfer the hidden _keys to normal keys that will be printed if they are different from 0
        if (mut.entry_exists("_del_start") && (mut["_del_start"] != "0")) mut["del_start"] = mut["_del_start"];
        if (mut.entry_exists("_del_end")   && (mut["_del_end"] != "0"))   mut["del_end"] = mut["_del_end"];
        
        if (mut.entry_exists("_ins_start") && (mut["_ins_start"].length() != 0)) mut["ins_start"] = mut["_ins_start"];
        if (mut.entry_exists("_ins_end")   && (mut["_ins_end"].length() != 0))   mut["ins_end"] = mut["_ins_end"];
        
				if (verbose)
					cout << mut["_gap_left"] << " :: " << mut["_gap_right"] << endl;
        
        // print out everything
        if (verbose)
        {
          cout << "== J1 ==" << endl << j1 << endl;
          for(map<string,string>::iterator it=j1.begin(); it!=j1.end(); it++)
          {
            cout << it->first << " = " << it->second << endl; 
          }
          
          cout << "== J2 ==" << endl;
          for(map<string,string>::iterator it=j2.begin(); it!=j2.end(); it++)
          {
            cout << it->first << " = " << it->second << endl; 
          }
          
          cout << "== Mut ==" << endl;
          for(map<string,string>::iterator it=mut.begin(); it!=mut.end(); it++)
          {
            cout << it->first << " = " << it->second << endl; 
          }
        }
        
        // @JEB 08-06-13 
        // Reassign read counts with required overlap to avoid counting reads that go into the
        // target site duplication from counting against the IS element (giving a non-100% value).
        // This value must be positive...
        
        int32_t duplication_require_overlap = max(n(mut["duplication_size"]), 0);
        
        // Snapshot before recounting, so the consensus-mode rejection path below can put the
        // entries back by assignment instead of running the whole pileup a second time to undo
        // itself. Recounting is not free: each call re-fetches every read overlapping three
        // windows across two BAMs.
        cDiffEntry j1_before_offset(j1);
        cDiffEntry j2_before_offset(j2);

        if (verbose) cerr << "Before 1:" << endl << j1 << endl;
        j1["read_count_offset"] = to_string(duplication_require_overlap);
        assign_one_junction_read_counts(settings, summary, j1, reference_jrc, junction_jrc);
        if (verbose) cerr << "After 1:" << endl << j1 << endl;
        
        if (verbose) cerr << "Before 2:" << endl << j2 << endl;
        j2["read_count_offset"] = to_string(duplication_require_overlap);
        assign_one_junction_read_counts(settings, summary, j2, reference_jrc, junction_jrc);
        if (verbose) cerr << "After 2:" << endl << j2 << endl;
        
        if (!settings.polymorphism_prediction) { // consensus mode
          // at this point the prediction type (CONSENSUS/POLYMORPHISM) should be honored
          // in determining whether to include this mutation.
          // Both junctions must be consensus. This previously tested j1 twice and never looked
          // at j2, so a MOB whose second junction was not consensus was accepted regardless.
          if ( (j1[PREDICTION]!="consensus") || (j2[PREDICTION]!="consensus") ) {
            
            // Restore the pre-offset entries. This used to re-run the counting to undo itself.
            j1 = j1_before_offset;
            j2 = j2_before_offset;
            continue;
          }
          
        } else { // polymorphism mode
          // @JEB 01-04-13 updated 09-22-17
          // Calculate a frequency for the mobile element insertion from the reads supporting the new and old junctions on each side
          
          // ref_read_count_1 is the total reads supporting either old junction (side 1/2) in the reference normalized by the number of sides
          // new_read_count_1 a1 is total counts supporting old junctions (added together)
          // function returns false if it cannot determine the counts for/against, which happens when both ref junctions are in repeats
          double new_read_count_1, total_read_count_1;
          bool j1_has_ref_reads = gd.read_counts_for_entry(j1, new_read_count_1, total_read_count_1);
          
          double new_read_count_2, total_read_count_2;
          bool j2_has_ref_reads = gd.read_counts_for_entry(j2, new_read_count_2, total_read_count_2);

          double frequency;
          if (j1_has_ref_reads || j2_has_ref_reads) {
            frequency = (new_read_count_1 + new_read_count_2) / (total_read_count_1 + total_read_count_2);
            mut[FREQUENCY] = formatted_double(frequency, settings.polymorphism_precision_places, true).to_string();
          } else {
            // Can't calculate a frequency if no sides of the junction fall in unique sequence
            mut[FREQUENCY] = "NA";
          }
          
          if (mut[FREQUENCY] == "NA") {
            //Don't predict mutations that have no frequency
            continue;
          }
          
          // Determine consensus vs. polymorphism from exact (Clopper-Pearson) confidence bounds on
          // the combined read counts, matching how RA and JC evidence are classified: the lower
          // bound gates "is this a real variant at all", the upper bound gates "is it fixed rather
          // than polymorphic". Testing the high side with the lower bound would require
          // ln(alpha)/ln(cutoff) reads even at 100% frequency and would drop real mutations for
          // want of depth.
          double combined_new = new_read_count_1 + new_read_count_2;
          double combined_total = total_read_count_1 + total_read_count_2;
          double freq_lower = binomial_frequency_lower_bound(combined_new, combined_total);
          double freq_upper = binomial_frequency_upper_bound(combined_new, combined_total);

          if ((settings.polymorphism_frequency_cutoff > 0.0)
              && (freq_lower < settings.polymorphism_frequency_cutoff)) {
            mut.add_reject_reason("FREQUENCY_CUTOFF");
            // @JEB 08-08-13 we might want to keep the mutation as rejected. This discards completely.
            continue;
          }
          if ((settings.consensus_frequency_cutoff > 0.0)
              && (freq_upper >= settings.consensus_frequency_cutoff)) {
            mut[FREQUENCY] = "1";
          } else if (!settings.polymorphism_prediction) {
            // Consensus mode reports only fixed differences; an intermediate frequency it cannot
            // call fixed is still not a consensus mutation.
            mut.add_reject_reason("FREQUENCY_CUTOFF");
            continue;
          }
        }
        
        // @JEB 12-22-12
        // Add link to missing coverage evidence that corresponds to the deleted region for negative overlap          
        if ( n(mut["duplication_size"]) < 0 ) {
          
          for(diff_entry_list_t::iterator mc_it = mc.begin(); mc_it != mc.end(); mc_it++) {
            cDiffEntry& mc_item = **mc_it;
            
            if  (  ( mc_item[SEQ_ID] == mut[SEQ_ID])
                 && ( n(mc_item[START]) == n(mut[POSITION]) ) 
                 && ( -n(mut["duplication_size"]) == n(mc_item[END]) - n(mc_item[START]) + 1)
                 )
            {
              mut._evidence.push_back(mc_item._id);
              mc.erase(mc_it);
              break; // If you don't break here, you may be past the end of the array
            }
          }
          
        }
        
        // Remove the two JC evidence items that we used from the list.
        // --> be sure to remove JC2 first so we don't invalidate JC1
        diff_entry_list_t::iterator jc2_it(jc1_it);
        do  {
          ASSERT(jc2_it != jc.end(), "Could not find 2nd junction used to predict MOB.");
          jc2_it++;
        } while (**jc2_it != j2);
          
        jc.erase(jc2_it);
        jc1_it = jc.erase(jc1_it); // iterator is now past element erased
        jc1_erased = true; // this tells us to not increment the iterator
        
				gd.add(mut);
				break; // next JC1
			}
      
      // because this is a list, we only advance if we didn't delete the current entry
      if (!jc1_erased) jc1_it++;
		}

    // @JEB 2019-11-29
    // We need to assign "unique=X" tags when the same MOB mutations are predicted
    // because two IS element copies can have the same name but different sequences
    // and these two copies can insert in the exact same way!
    
    diff_entry_list_t mobs_to_add_ids;
    diff_entry_list_t mobs = gd.get_list(make_vector<gd_entry_type>(MOB));
    
    diff_entry_list_t::iterator it2;
    for (diff_entry_list_t::iterator it1 = mobs.begin(); it1 != mobs.end(); it1++) {
      
      it2 = it1;
      it2++;
      
      if ((it2 != mobs.end()) && (**it1 == **it2)) {
        mobs_to_add_ids.push_back(*it2);
      }
    }

    // Now give them unique IDs that are numbered
    uint32_t mob_id = 1;
    for (diff_entry_list_t::iterator it = mobs_to_add_ids.begin(); it != mobs_to_add_ids.end(); it++) {
      cDiffEntry & de = **it;
      de["unique"] = "mob" + to_string(mob_id);
      mob_id++;
    }
    
  }
  
  // Deletions between two homologous copies of a sequence -- paralogous genes, an rRNA operon, any
  // repeat that is not an annotated repeat_region -- leave NO split-read junction: a read crossing the
  // breakpoint aligns just as well to either copy, so there is nothing for candidate-junction
  // identification to find. predictMCplusJCtoDEL therefore drops them (its case (2) requires annotated
  // repeats and everything else requires a JC), and the mutation is lost even though the coverage
  // evidence is unambiguous.
  //
  // What survives such a deletion is a single HYBRID copy: the left copy up to the crossover, the right
  // copy after it. Every step below follows from that. The register d (which is also the deletion size)
  // is the offset that aligns the two copies, and it can be read straight off the reference sequence,
  // seeded by the MC boundaries. The crossover is then the one place where the reads stop showing the
  // left copy's bases and start showing the right copy's, measured at the columns where the two copies
  // actually differ. See homologous_deletion.h for why pooling the two copies' piles recovers that
  // series at full depth, and why the reads carrying it are invisible as RA evidence.
  //
  // A discordant pair, when the library is paired and one is present, gives a better seed than the MC
  // boundaries and is worth recording as evidence -- but it is not required, and nothing else here
  // depends on it.
  void MutationPredictor::predictMCtoDELbyHomology(Settings& settings, Summary& summary, cGenomeDiff& gd,
                                                   diff_entry_list_t& mc, diff_entry_list_t& dp, diff_entry_list_t& ra)
  {
    if (settings.skip_homologous_deletion_prediction) return;
    if (mc.empty()) return;

    int32_t read_length = static_cast<int32_t>(summary.sequence_conversion.read_length_avg + 0.5);
    if (read_length <= 0) return;

    // A DP side is placed at the innermost aligned end of its supporting pairs and an MC boundary is
    // where coverage falls off, so against a homologous boundary they disagree by up to about one read.
    const int32_t kMCtoDPSlop = read_length;

    // Window for measuring whether a candidate register really is a homologous alignment. Deliberately
    // narrow and NOT scaled by read length: at +/-75 the true register scores 0.95-0.98 on the REL606
    // colanic-acid deletion while the best wrong one scores 0.40, and widening it runs the measurement
    // off the end of the homologous block into unrelated sequence (at +/-150 one of these MCs, whose
    // boundary sits only 50 bp inside the block, drops to 0.80 and would be rejected).
    const int32_t kRegisterWindow = 75;

    // The MC boundaries can be wrong by as much as the homologous block is long, because that is the
    // stretch over which redundancy suppresses unique coverage. The worst real case seen is 1181 bp;
    // this leaves room for blocks up to rRNA-operon scale. Scoring one register is ~300 byte
    // comparisons, so a scan this wide is still only single-digit milliseconds per MC.
    const int32_t kRegisterSearchHalfWidth = max(2 * read_length, 5000);

    // A high identity on its own means nothing between paralogs that are identical in long stretches
    // everywhere; what marks the true register is that nothing else comes close. Note this deliberately
    // does NOT use settings.consensus_frequency_cutoff, which is 0.50 in consensus mode -- low enough
    // to admit wrong registers that measure 0.49.
    const double kMinRegisterIdentity = 0.85;
    const double kMinRegisterMargin = 0.25;
    const int32_t kRegisterMarginExclusion = 2 * read_length;

    const int32_t kBlockStep = 50;
    const double kBlockIdentity = 0.80;
    const int32_t kBlockMaxExtent = 20000;
    // Anything a single read can span is a junction, and belongs to predictMCplusJCtoDEL.
    const int32_t kMinBlockLength = max(read_length, 50);

    // The false-positive gate. Under the deletion hypothesis only one copy exists in the sample and a
    // column reads ~1.0; with both copies still present it reads ~0.5.
    const double kMinAlleleFraction = 0.80;
    // Two calls could come from a single mismapped pile; three spanning more than the typical distance
    // between discriminating columns could not.
    const int32_t kMinColumnsEachSide = 3;
    const double kMinReadConsistency = 0.90;

    // The reads that carry this signal are exactly the ones base calling discards as redundant, so the
    // BAM is the primary source. RA evidence can still bracket the crossover when the reads happened to
    // place uniquely, which is the only thing available to callers that have no BAM at all (gdtools).
    bool have_bam = file_exists(settings.reference_bam_file_name.c_str())
                 && file_exists(settings.reference_fasta_file_name.c_str());
    homology_base_counter* counter = NULL;

    // Observed base per position, for every RA that called something different from the reference. At a
    // deletion between near-identical copies these are not mutations at all: they are the surviving
    // hybrid showing the OTHER copy's base where reads mapped to this one, which is exactly the signal
    // that locates the crossover.
    map<string, map<int32_t, string> > ra_base;
    for (diff_entry_list_t::iterator it = ra.begin(); it != ra.end(); it++) {
      cDiffEntry& r = **it;
      if (!r.entry_exists(NEW_BASE) || !r.entry_exists(SEQ_ID) || !r.entry_exists(POSITION)) continue;
      if (r[NEW_BASE].size() != 1) continue;          // only substitutions are informative here
      if (r.entry_exists(INSERT_POSITION) && (n(r[INSERT_POSITION]) != 0)) continue;
      // Only FIXED differences say anything about which copy the hybrid inherited. In polymorphism mode
      // the same list also carries variants at any frequency down to polymorphism_frequency_cutoff, and
      // a minority allele that happens to match the other copy would place the crossover wherever it
      // sits. Require consensus support, the same bar used to call a base at all.
      if (r.entry_exists(FREQUENCY)
          && (from_string<double>(r[FREQUENCY]) < settings.consensus_frequency_cutoff)) continue;
      ra_base[r[SEQ_ID]][n(r[POSITION])] = r[NEW_BASE];
    }

    for (diff_entry_list_t::iterator mc_it = mc.begin(); mc_it != mc.end(); mc_it++) {
      cDiffEntry& mc_item = **mc_it;
      const string& seq_id = mc_item[SEQ_ID];
      int32_t mc_start = n(mc_item["start"]), mc_end = n(mc_item["end"]);

      // --- find a DP with deletion geometry at these boundaries -------------------------------------
      // side_1 strand -1 keeps the flank at coords <= p (the left flank survives); side_2 strand +1
      // keeps coords >= p (the right flank survives). That is a deletion of what lies between them.
      cDiffEntry* dp_item = NULL;
      for (diff_entry_list_t::iterator dp_it = dp.begin(); dp_it != dp.end(); dp_it++) {
        cDiffEntry& d = **dp_it;
        if (d.entry_exists(REJECT)) continue;
        if (d[SIDE_1_SEQ_ID] != seq_id || d[SIDE_2_SEQ_ID] != seq_id) continue;
        if (n(d[SIDE_1_STRAND]) != -1 || n(d[SIDE_2_STRAND]) != +1) continue;
        if (abs(n(d[SIDE_1_POSITION]) - mc_start) > kMCtoDPSlop) continue;
        if (abs(n(d[SIDE_2_POSITION]) - mc_end) > kMCtoDPSlop) continue;
        dp_item = &d;
        break;
      }

      cAnnotatedSequence& seq = ref_seq_info[seq_id];
      int32_t seq_len = static_cast<int32_t>(seq.m_length);

      // --- find the homologous register ---------------------------------------------------------------
      // The deletion joins x to x + d, so the two copies are the same sequence offset by d. That is a
      // property of the reference alone, and the MC gives a seed good to a few hundred bases -- a DP,
      // when there is one, gives a better one. Neither the MC boundaries nor the longest run of
      // identical bases IS the answer: coverage falls off gradually through homology, and between two
      // 95%-identical paralogs long identity runs are everywhere, most of them nowhere near the
      // crossover. What identifies the register is that nothing else comes close to it.
      int32_t d_seed = mc_end - mc_start + 1;
      if (dp_item != NULL)
        d_seed = n((*dp_item)[SIDE_2_POSITION]) - n((*dp_item)[SIDE_1_POSITION]);
      if (d_seed <= 0) continue;

      double best_identity = 0.0, runner_up_identity = 0.0;
      int32_t best_d = find_homologous_register(seq, mc_start, mc_end, d_seed,
                                                kRegisterSearchHalfWidth, kRegisterWindow,
                                                kMinRegisterIdentity, kMinRegisterMargin,
                                                kRegisterMarginExclusion,
                                                best_identity, runner_up_identity);
      if (best_d == 0) continue;

      // --- measure the homologous block ---------------------------------------------------------------
      // Either MC boundary may have fallen outside the block, so grow from both and keep the better.
      int32_t block_start = 0, block_end = 0;
      homologous_block_extent(seq, best_d, mc_start, kBlockStep, kBlockIdentity, kBlockMaxExtent,
                              block_start, block_end);
      {
        int32_t alt_start = 0, alt_end = 0;
        homologous_block_extent(seq, best_d, mc_end - best_d, kBlockStep, kBlockIdentity, kBlockMaxExtent,
                                alt_start, alt_end);
        if ((alt_end - alt_start) > (block_end - block_start)) {
          block_start = alt_start; block_end = alt_end;
        }
      }
      int32_t block_length = block_end - block_start + 1;
      if (block_length < kMinBlockLength) continue;
      if ((block_start < 1) || (block_end + best_d > seq_len)) continue;
      // The two copies have to be separate stretches of the genome. If they overlap, this is a tandem
      // array rather than a pair of copies, the block and its image are the same reference positions,
      // and pooling them would count the same reads twice. Tandem repeats are handled elsewhere.
      if (best_d < block_length) continue;

      // The MC boundaries can only be wrong by as much as the block is long -- that is the stretch over
      // which redundancy suppresses unique coverage -- plus about a read at each end.
      if (abs(d_seed - best_d) > block_length + 2 * read_length) continue;

      vector<int32_t> positions;
      discriminating_positions(seq, best_d, block_start, block_end, positions);
      if (positions.empty()) continue;

      // --- read out which copy the hybrid carries at each discriminating column ------------------------
      homology_columns_t columns;
      int32_t last_left = 0, first_right = 0, n_left = 0, n_right = 0, n_violations = 0;
      bool found = false;

      if (have_bam) {
        if (counter == NULL)
          counter = new homology_base_counter(settings.reference_bam_file_name,
                                              settings.reference_fasta_file_name,
                                              settings.base_quality_cutoff);

        set<int32_t> left_positions, right_positions;
        for (size_t i = 0; i < positions.size(); i++) {
          left_positions.insert(positions[i]);
          right_positions.insert(positions[i] + best_d);
        }

        counter->clear();
        counter->set_read_recording(left_positions, 0);
        counter->count_region(seq_id, block_start, block_end);
        counter->set_read_recording(right_positions, best_d);
        counter->count_region(seq_id, block_start + best_d, block_end + best_d);

        // Scale the per-column floor to the run's own depth, so this works at 10x and at 100x. The
        // surviving flank contributes essentially full coverage to whichever side of the flip a column
        // falls on; well below a fifth of average there is noise or mismapping, not the hybrid.
        uint32_t min_column_coverage = 3;
        CoverageSummaries::const_iterator cov = summary.unique_coverage.find(seq_id);
        if (cov != summary.unique_coverage.end())
          min_column_coverage = max(min_column_coverage,
                                    static_cast<uint32_t>(0.20 * cov->second.average));

        build_homology_columns(seq, best_d, positions, counter->tallies(),
                               min_column_coverage, kMinAlleleFraction, columns);

        // A real crossover is a perfect step. Allow one column of slop for an error or a homopolymer
        // artifact, and 5% once the block is long enough for that to be the looser bound.
        int32_t n_called = 0;
        for (size_t i = 0; i < columns.size(); i++)
          if (columns[i].call != 0) n_called++;

        found = find_crossover_interval(columns, kMinColumnsEachSide,
                                        max(1, static_cast<int32_t>(0.05 * n_called)),
                                        last_left, first_right, n_left, n_right, n_violations);

        // A read that shows the right copy's allele and then the left copy's further along cannot have
        // come from the hybrid. If that keeps happening, whatever is here is not one clean crossover.
        if (found) {
          int32_t informative_reads = 0;
          double consistency = read_allele_consistency(seq, best_d, counter->read_alleles(),
                                                       informative_reads);
          if (consistency < kMinReadConsistency) found = false;
        }
      }

      // Reads that placed uniquely still leave RA evidence, which is all a caller with no BAM has, and
      // it can bracket the crossover from a handful of calls. That is a much thinner basis than the
      // columns above, so only lean on it when there is no BAM to consult, or when a discordant pair
      // has independently confirmed the deletion geometry.
      if (!found && (!have_bam || (dp_item != NULL))) {
        map<string, map<int32_t, string> >::const_iterator rb = ra_base.find(seq_id);
        if (rb == ra_base.end()) continue;
        build_homology_columns_from_ra(seq, best_d, positions, rb->second, columns);
        found = find_crossover_interval(columns, 1, 0,
                                        last_left, first_right, n_left, n_right, n_violations);
      }
      if (!found) continue;

      // The crossover lies in (last_left, first_right]; every position in there yields identical
      // sequence, so the START is ambiguous by that much while the SIZE is exact -- both endpoints
      // slide together. Emit the rightmost consistent start, matching the right-shift convention
      // normalize_to_sequence() applies to DEL, which is also the one choice in that interval that
      // does not need to be qualified.
      int32_t position = first_right;
      if ((position < 1) || (position + best_d - 1 > seq_len)) continue;

      vector<string> evidence_ids = make_vector<string>(mc_item._id);
      if (dp_item != NULL) evidence_ids.push_back(dp_item->_id);

      gd.add(make_homologous_deletion(seq, position, best_d, evidence_ids,
                                      settings.polymorphism_prediction));
      mc_it = mc.erase(mc_it);
      mc_it--;
    }

    if (counter != NULL) delete counter;
  }

  void MutationPredictor::predictJCtoINSorSUBorDEL(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc, bool use_redundant_sides)
  {
    (void)summary;
    (void)mc;
    bool verbose = false;
    
    // variables pulled from the summary
    // JEB 2022-11-19 changed to use kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation below
    // int32_t read_length_avg = static_cast<int32_t>(summary.sequence_conversion.read_length_avg);
    
    // @JEB 03-09-2014 Changed this section to produce INS and DEL instead of AMP,
    //                 when the mutation size is less than or equal to threshold in settings. 
		///
    for(diff_entry_list_t::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) //JC
		{
			cDiffEntry& j = **jc_it;
      
      if (verbose) {
        cout << j.as_string() << endl;
      }
      //cout << j["side_1_seq_id"] << " " << j["side_2_seq_id"] << endl;
      //cout << j["side_1_strand"] << " " << j["side_2_strand"] << endl;
      
      // must be on same sequence
      if ((j["side_1_seq_id"] != j["side_2_seq_id"]))
        continue;
      
      // must be in same orientation (implies strands are opposite)
      int32_t side_1_strand = n(j["side_1_strand"]);
      int32_t side_2_strand = n(j["side_2_strand"]);
      
      // to be safe about making predictions, neither can be ambiguous, unless we're in the pass
      // when we check for ignoring indels
      if (!use_redundant_sides) {
        if ( (n(j["side_1_redundant"]) == 1) || (n(j["side_2_redundant"]) == 1) )
          continue;
      }
      
      if (side_1_strand == side_2_strand)
				continue;
      
			string seq_id = j["side_1_seq_id"];
      int32_t side_1_position = n(j["side_1_position"]);
      int32_t side_2_position = n(j["side_2_position"]);
      
			// We can assume that the lower coordinate will be first since this is NOT a deletion
			// (which would be handled above)
      // By this point any positive overlap should have been resolved.
			ASSERT(n(j["overlap"]) <= 0, "Non-zero overlap in junction when predicting INS/SUB.");
      
			// mutation will always be after this position
			// Special case of circular chromosome
			if ( (side_1_position == 1) && ( side_2_position== ref_seq_info[ref_seq_info.seq_id_to_index(j["side_2_seq_id"])].m_length ) )
			{
				j[IGNORE] = "CIRCULAR_CHROMOSOME";
				continue;
			}
      
      // If we are predicting a very big insertion (longer than read length), 
      // it is likely spurious. Require other evidence to convert to a mutation.
      
      // Decide whether we are reverse complementing unique sequence
      // based on the original left side coords!! -- before swap below.
      bool reverse_complement_unique_sequence = (side_1_strand == 1);
      
      // * We change everything here so that the side with the lower coordinate is first
      // Due to overlap resolution, this can change from the time when we sorted this way... *Ugh* .. so we fix it here.
      if (side_2_position < side_1_position) {
        swap(side_1_position, side_2_position);
        swap(side_1_strand, side_2_strand);
      }
      int32_t size = side_2_position - side_1_position + 1;
      
      // This implies a deletion, and not counting the endpoint nucleotides (hence -1 instead of +1)
      // Note that this is a negative number!
      if (side_1_strand == -1 )
        size = side_1_position - side_2_position + 1;
      
      // Insertion or deletion must be smaller than this setting to be predicted by this evidence alone.
      // @JEB 2022-11-19 changed from being dependent on read length to this constant
      if (abs(size) > kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation)
        continue;
      
      // At this point, we have committed to adding a mutation...
      cDiffEntry mut;
      mut._evidence = make_vector<string>(j._id);
    
      if (size < 0) // this is a deletion!
      {

        if (!j.entry_exists("unique_read_sequence"))
        {    
          //'DEL' if there is no read-only sequence
          //
          // Example (reverse_complement = false)
          //   REL606 2103888 -1 REL606 2103911 1
          //    side_1_pos, side_1_strand    side_2_pos, side_2_strand    size = 22
          //   Output
          //     DEL . . REL606 2103889 22
          //
          // Example (reverse_complement = true)
          //   REL606 2103911 1 REL606 2103888 -1 
          //    side_1_pos, side_1_strand    side_2_pos, side_2_strand    size = 22
          //   Output
          //     DEL . . REL606 2103889 22
          
          mut._type = DEL;
          
          mut
          ("seq_id", seq_id)
          ("position", s(side_1_position+1))
          ("size", s(-size))
          ;
        }
        else {
          //'SUB' = deletion with some unique sequence inserted in its place
          //
          // Example (reverse_complement = false)
          //   REL606 2103888 -1 REL606 2103911 1 unique_read_sequence=TTTTT
          //    side_1_pos, side_1_strand    side_2_pos, side_2_strand    size = 22
          //   Output
          //     SUB . . REL606 2103889 22 TTTTT
          //
          // Example (reverse_complement = true)
          //   REL606 2103911 1 REL606 2103888 -1 REL606 unique_read_sequence=AAAAA
          //    side_1_pos, side_1_strand    side_2_pos, side_2_strand    size = 22
          //   Output
          //     SUB . . REL606 2103889 22 TTTTT
          
          mut._type = SUB;
          string new_seq = j["unique_read_sequence"];  
          if (reverse_complement_unique_sequence) new_seq = reverse_complement(new_seq);
          
          mut
          ("seq_id", seq_id)
          ("position", s(side_1_position+1))
          ("size", s(-size))
          ("new_seq", new_seq)
          ;
        }
      } 

      else // (size >= 0)
      {		 
        // 'INS' otherwise
        //
        // Example (reverse_complement = false)
        //   REL606 2103888 1 REL606 2103911 -1
        //    side_1_pos, side_1_strand    side_2_pos, side_2_strand    size = 24
        //   Output
        //     INS . . REL606 2103911 CAGCCAGCCAGCCAGCCAGCCAGC
        //
        // Example (reverse_complement_unique_sequence = true)
        //   REL606 2103888 1 REL606 2103911 -1
        //    side_1_pos, side_1_strand    side_2_pos, side_2_strand    size = 24
        //   Output
        //     INS . . REL606 2103911 CAGCCAGCCAGCCAGCCAGCCAGC
        //
        // Example with unique_read_sequence (reverse_complement = false)
        //   REL606 2103911 -1 REL606 2103888 1 unique_read_sequence=CAA
        //    side_1_pos, side_1_strand    side_2_pos, side_2_strand    size = 24
        //   Output
        //     INS . . REL606 2103911 CAA|CAGCCAGCCAGCCAGCCAGCCAGC 
        //
        // Example with unique_read_sequence (reverse_complement = true)
        //   REL606 2103888 1 REL606 2103911 -1 unique_read_sequence=TTG
        //    side_1_pos, side_1_strand    side_2_pos, side_2_strand    size = 24
        //   Output
        //     INS . . REL606 2103911 CAA|CAGCCAGCCAGCCAGCCAGCCAGC
        //
        // Example with unique_read_sequence 
        //   2108470	-1	REL606	2108471	1	-8	unique_read_sequence=CAGCCAGC
        //    side_1_pos, side_1_strand    side_2_pos, side_2_strand    size = 0

        
        mut._type = INS;
        
        string ins_seq;
      
        // Special case of size = 0 means that we don't add any reference sequence positions
        if (size != 0)
          ins_seq = ref_seq_info.get_sequence_1(seq_id, side_1_position, side_2_position);
      
        if (j.entry_exists("unique_read_sequence")) {
          string unique_seq = j["unique_read_sequence"];  
          if (reverse_complement_unique_sequence) unique_seq = reverse_complement(unique_seq);
          ins_seq = unique_seq + ins_seq;
        }
        
        mut		
        ("seq_id", seq_id)		
        ("position", (size != 0) ? s(side_2_position) : s(side_1_position))	
        // previously was side_2_position
        //
        // If there is a size, then we are really adding some repeated bases
        //  ===========0123
        //              123=========
        //  and we insert after position 3
        //
        //  If there was not a size, then we have completely novel sequence in between
        //
        //   ==========01
        //                 23==========
        //  and we insert after position 1
        //
        /// If there was both a size and unique sequence, then we still insert after position 3
        //   ==========0123
        //                 AB
        //                   123==========
        //
        ("new_seq", ins_seq)		
        ;		
      }
      
      // If we are in polymorphism mode, propagate the frequency from JC evidence to mutation.
      // mutation_frequency() re-applies the consensus snap: evidence keeps its fitted estimate in
      // [frequency] and records the call in [prediction], while a mutation reports a consensus
      // junction as exactly 1.
      if (settings.polymorphism_prediction)
        mut[FREQUENCY] = j.mutation_frequency();
      
      // Finally, add it
      gd.add(mut);
		}

  }
  
  void MutationPredictor::predictRAtoSNPorDELorINSorSUB(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& ra, diff_entry_list_t& mc)
  {
    (void)summary;
    (void)mc;
    bool verbose = false;
    
		///
		// Ignore RA that overlap DEL, MC unless they are user
		// They are due to low spurious coverage in deleted regions or replicate bases changed!
		///
    
    vector<gd_entry_type> del_types = make_vector<gd_entry_type>(DEL)(INS)(SUB);
    diff_entry_list_t del = gd.get_list(del_types);
    
    // Don't add deleted flags if we are in targeted sequencing mode
    if (!settings.targeted_sequencing) {
      
      for(diff_entry_list_t::iterator ra_it = ra.begin(); ra_it != ra.end(); ra_it++) //RA
      {
        cDiffEntry& ra_item = **ra_it;
        
        bool next_ra = false;
        
        // don't remove these, no matter what!
        if (ra_item.entry_exists(USER_DEFINED))
          continue;
        
        for(diff_entry_list_t::iterator del_it = del.begin(); del_it != del.end(); del_it++) //DEL
        {
          cDiffEntry& del_item = **del_it;
          
          if (ra_item.located_within(del_item))
          {
            ra_item["deleted"] = "1";
            next_ra = true;
            break;
          }
        }
        if (next_ra) continue;
        
        // @JEB Unless an option is supplied, mark anything that overlaps MC as 'deleted' so it will be ignored
        if (!settings.call_mutations_overlapping_missing_coverage) {
          for(diff_entry_list_t::iterator mc_it = mc.begin(); mc_it != mc.end(); mc_it++) //MC
          {
            cDiffEntry& mc_item = **mc_it;
            
            if (ra_item.located_within(mc_item))
            {
              ra_item["deleted"] = "1";
              break;
            }
          }
        }
        
      }
    }
    
		ra.sort(MutationPredictor::sort_by_pos);
    
		///
		// Gather together read alignment mutations that occur next to each other
		// ...unless they are polymorphisms
		///
    
		bool first_time = true;
		cDiffEntry mut;
		vector<cDiffEntry> muts;
    
    
    for(diff_entry_list_t::iterator ra_it = ra.begin(); ra_it != ra.end(); ra_it++) //RA
    {
      cDiffEntry& item = **ra_it;
      
      ///DEBUG
      //cout << item.as_string() << endl;
      
      string ra_seq_id = item[SEQ_ID];
      int32_t ra_position = from_string<int32_t>(item[POSITION]);
      string ra_ref_base;
      if (item.entry_exists(INSERT_POSITION) && (from_string<int32_t>(item[INSERT_POSITION]) != 0)) {
        ra_ref_base = ".";
      } else {
        ra_ref_base = ref_seq_info.get_sequence_1(ra_seq_id, ra_position, ra_position);
      }
      string ra_new_base = (item[MAJOR_BASE] == ra_ref_base) ? item[MINOR_BASE] : item[MAJOR_BASE];
      
      // Whether this position was called fixed or mixed. Read from [prediction] rather than by
      // testing [frequency] against 1: evidence stores its fitted estimate there, so a fixed call
      // reads back as 0.982 rather than 1 and a frequency test would see every consensus SNP as a
      // polymorphism. The mutation built below re-snaps via mutation_frequency().
      bool ra_is_consensus = (item[PREDICTION] == "consensus");
      
      // Sometimes a SNP might be called in a deleted area because the end was wrong,
			// but it was corrected using a junction. (This catches this case.)
			if ( (!item.entry_exists(USER_DEFINED)) && (item.entry_exists("reject") || item.entry_exists("deleted")) )
			  continue;
      
      // If we are predicting mixed bases and not polymorphisms, then don't create
      // mutations for mixed frequency predictions (leave them as unassigned RA evidence)
      if (!settings.polymorphism_prediction && !ra_is_consensus) {
        continue;
      }
      
			bool same = false;
			if (!first_time)
			{
        // @JEB: 2015-11-25 comments enable joining discontinguous INS evidence
        // which can lead to multiple insertions at same position+insert_position
				if ( ((mut["end"] == item["position"]) && (n(mut["insert_end"]) + 1 == n(item["insert_position"])))
						|| ((n(mut["end"]) + 1 == n(item["position"])) && (item["insert_position"] == "0")) )
					same = true;
        
        // This code is only safe if every mutation has a frequency
        if (settings.polymorphism_prediction) {
          if ( !ra_is_consensus || (mut[FREQUENCY] != "1") //don't join polymorphisms
              || (mut[SEQ_ID] != item[SEQ_ID]) )
            same = false;
        }
			}
      
			if (!same)
			{
				if (!first_time) muts.push_back(mut);
				first_time = false;
				cDiffEntry new_mut;
        new_mut._evidence = make_vector<string>(item._id);
				new_mut
        ("seq_id", item["seq_id"])
        ("position", item["position"])
        ("start", item["position"])
        ("end", item["position"])
        ("insert_start", item["insert_position"])
        ("insert_end", item["insert_position"])
        ("ref_seq", (ra_ref_base != ".") ? ra_ref_base : "")
        ("new_seq", (ra_new_base != ".") ? ra_new_base : "")
				;
        
        if (settings.polymorphism_prediction) {
          new_mut[FREQUENCY] = item.entry_exists(FREQUENCY) ? item.mutation_frequency() : "1";
        }
				mut = new_mut;
			}
			else
			{
				mut
        ("insert_end", item["insert_position"])
        ("end", item["position"])
				;
        mut["ref_seq"] += (ra_ref_base != ".") ? ra_ref_base : "";
				mut["new_seq"] += (ra_new_base != ".") ? ra_new_base : "";
				mut._evidence.push_back(item._id);
			}
		}
		//don't forget the last one
		if (!first_time) muts.push_back(mut);
    
		///
		// Finally, convert these items into the fields needed for the various types of mutations
		///
    
    cDiffEntry last_mut;
		for (uint32_t i = 0; i < muts.size(); i++) {
      
      if (i>0) last_mut = mut;
			mut = muts[i];
      
      // insertion and amplification
			if (mut["ref_seq"].size() == 0)
			{
        mut._type = INS;
				// unused fields
				mut.erase("ref_seq");
        
        if (settings.polymorphism_prediction) {
          // This is a special case to keep ordering of multiple inserted bases after 
          // the same base (without it the order is unknown in poly mode
          
          ASSERT( (mut[FREQUENCY]=="1") || (mut["insert_start"] == mut["insert_end"]), "Polymorphism has incorrectly merged INS mutations.");
          string debug_ins_pos = mut["insert_start"];           
          mut["insert_position"] = mut["insert_start"];
        } else { // CONSENSUS mode
          

          // We have a problem sometimes with certain non-adjacent columns passing the
          // score threshold for example, due to differences in error rates
          // the C columns might pass with insert_positions 3, 6, 9.
          //
          //   *  *  *
          // TTCTTCTTC
          
#define NO_DISCONTIGUOUS_INSERTS
          
#ifdef NO_DISCONTIGUOUS_INSERTS
          //
          // Do not accept predictions that don't start at an insert position of 1
          if (n(mut["insert_start"]) != 1) {
            continue;
          }
#else
          // Number insert_positions for INS mutations continuously
          // (starting with implicit nothing = 1) for failing the next condition...
          if ( (last_mut._type == INS) && (last_mut[POSITION] == mut[POSITION]) ) {
            //cout << last_mut.as_string() << endl << mut.as_string() << endl;
            
            uint32_t last_insert_position = (last_mut.entry_exists(INSERT_POSITION) ? n(last_mut[INSERT_POSITION]) : 1);
            last_insert_position++;
            mut["insert_position"] = s(last_insert_position);
          }
#endif
        }
			}
			// deletion
			else if (mut["new_seq"].size() == 0)
			{
        mut._type = DEL;
				mut["size"] = s(n(mut["end"]) - n(mut["start"]) + 1);
        
				// unused fields
				mut.erase("new_seq");
				mut.erase("ref_seq");
			}
			// block substitution
			else if ((mut["ref_seq"].size() > 1) || (mut["new_seq"].size() > 1))
			{
        // This loop will go through the RA evidence for this SUB, and find the lowest
        // RA that has a ref_base of A,C,T, or G.        
        int32_t iRefPos = -1;
        for(diff_entry_list_t::iterator ra_it = ra.begin(); ra_it != ra.end(); ra_it++)
        {
          cDiffEntry& ra_evid = **ra_it;
          
          if((find(mut._evidence.begin(), mut._evidence.end(), ra_evid._id) != mut._evidence.end()) && (iRefPos < 0 || iRefPos > from_string<int32_t>(ra_evid["position"])) && ra_evid["ref_base"] != ".")
            iRefPos = from_string<int32_t>(ra_evid["position"]);
        }
        
        if(iRefPos > -1)mut["position"] = to_string(iRefPos);
        
        mut._type = SUB;
				mut["size"] = s(mut["ref_seq"].size());
				mut.erase("ref_seq");
			}
			//snp
			else
			{
				mut.erase("ref_seq");
        mut._type = SNP;
			}
      
			mut.erase("start");
			mut.erase("end");
			mut.erase("insert_start");
			mut.erase("insert_end");
      
			gd.add(mut);
		}

  }

  void MutationPredictor::filter_JC_indel_evidence(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc)
  {
    (void) summary;
    
    // Filter out JC that predict hompolymer indels and/or polymorphisms
    // Also remove junctions that predicted for nanopore data
    
    cGenomeDiff junction_indel_gd;
    predictJCtoINSorSUBorDEL(settings, summary, junction_indel_gd, jc, mc, true);
    normalize_and_annotate_tandem_repeat_mutations(settings, summary, junction_indel_gd);

      
    vector<gd_entry_type> ins_del_types = make_vector<gd_entry_type>(INS)(DEL);
    diff_entry_list_t jc_prediction = junction_indel_gd.get_list(ins_del_types);
    for (diff_entry_list_t::iterator it=jc_prediction.begin(); it != jc_prediction.end(); it++) {
      diff_entry_ptr_t& mut = *it;
      
      // There should be only one evidence item and it better be a JC
      ASSERT(mut->_evidence.size() == 1, "Expected only one JC evidence item supporting INS/DEL mutation.");
      diff_entry_ptr_t ev = gd.find_by_id(mut->_evidence[0]);
      ASSERT(ev->_type == JC, "Expected only one JC evidence item supporting INS/DEL mutation.");
      
      // Figure out if it is a polymorphism or not...
      //cout << mut->as_string() << endl;
      //cout << ev->as_string() << endl;
      
      // Only check repeats of length one
      if (settings.consensus_reject_indel_homopolymer_length || settings.polymorphism_reject_indel_homopolymer_length) {
        
        if (mut->entry_exists("repeat_length") && ((*mut)["repeat_length"] == "1")) {
          
          ASSERT(mut->entry_exists("repeat_ref_copies"), "Expected 'repeat_ref_copies' field for JC evidence for INS/DEL.");
          uint32_t repeat_ref_copies = n((*mut)["repeat_ref_copies"]);
          
          ASSERT(ev->entry_exists("prediction"), "Expected 'prediction' field for JC evidence for INS/DEL.");
          
          bool reject = false;
          if ((*ev)["prediction"] == "consensus") {
            reject = (repeat_ref_copies >= settings.consensus_reject_indel_homopolymer_length);
          } else if ((*ev)["prediction"] == "polymorphism") {
            reject = (repeat_ref_copies >= settings.polymorphism_reject_indel_homopolymer_length);
          } else if ((*ev)["prediction"] != "unknown") {
            ERROR("Expected 'prediction' field for JC evidence for INS/DEL to be 'consensus', 'polymorphism', or 'unknown'\n"+ev->as_string());
          }
          
          // Add reject reasons to the JC entries, so that they will not be used below to re-predict these mutations
          if (reject) {
            ev->add_reject_reason("INDEL_HOMOPOLYMER");
          }
        }
      }
      
      // If the 'polymorphism_no_indels' option is set:
      // also hide any polymorphic indels (prediction=polymorphism)
      // and any with undetermined frequencies (prediction=unknown).
      // It's possible we should add SUB mutation here?
      
      if (settings.polymorphism_no_indels) {
        ASSERT(ev->entry_exists("prediction"), "Expected 'prediction' field for JC evidence for INS/DEL.");
        if ( ((*ev)["prediction"] == "polymorphism") || ((*ev)["prediction"] == "unknown") )
          ev->add_reject_reason("POLYMORPHIC_INDEL");
      }
    }
  }


  
  void MutationPredictor::remove_mutations_on_deleted_reference_sequences(Settings& settings, Summary& summary, cGenomeDiff& gd)
  {
    (void)settings;
    
    for (map<string,CoverageSummary>::iterator it=summary.unique_coverage.begin(); it!=summary.unique_coverage.end(); it++) {
      
      // This is how we know it was deleted
      if (it->second.deletion_coverage_propagation_cutoff >= 0.0) continue;
      gd.remove_mutations_on_deleted_reference_sequence(it->first, ref_seq_info[it->first].get_sequence_length());

    }
  }


  void MutationPredictor::ignore_evidence_near_contig_ends(Settings& settings, Summary& summary, cGenomeDiff& gd)
  {
    (void) summary;
    
    int32_t contig_end_distance = settings.ignore_within_this_multiple_of_average_read_length_of_contig_end
      * static_cast<int32_t>(round(summary.sequence_conversion.read_length_avg));
    
    // RA evidence
    vector<gd_entry_type> ra_types = make_vector<gd_entry_type>(RA);
    diff_entry_list_t ra = gd.get_list(ra_types);
    
    for (diff_entry_list_t::iterator ra_it=ra.begin(); ra_it!=ra.end(); ra_it++)
    {
      cDiffEntry& ra = **ra_it;
      
      if (!ref_seq_info[ra[SEQ_ID]].m_is_contig)
        continue;
      
      int32_t seq_length = ref_seq_info.get_sequence_length(ra[SEQ_ID]);

      bool near_contig_end = false;
      if (n(ra[POSITION]) <= contig_end_distance)
        near_contig_end = true;
          
      if (n(ra[POSITION]) >= seq_length-contig_end_distance)
        near_contig_end = true;
      
      if (near_contig_end) {
        ra[IGNORE] = "CONTIG_END";
      }
    }
    
    // MC evidence
    vector<gd_entry_type> mc_types = make_vector<gd_entry_type>(MC);
    diff_entry_list_t mc = gd.get_list(mc_types);
    
    for (diff_entry_list_t::iterator mc_it=mc.begin(); mc_it!=mc.end(); mc_it++)
    {
      cDiffEntry& mc = **mc_it;
      
      if (!ref_seq_info[mc[SEQ_ID]].m_is_contig)
        continue;
      
      int32_t seq_length = ref_seq_info.get_sequence_length(mc[SEQ_ID]);

      bool near_contig_end = false;
      if (n(mc[START]) <= contig_end_distance)
        near_contig_end = true;
          
      if (n(mc[END]) + n(mc[END_RANGE]) >= seq_length-contig_end_distance)
        near_contig_end = true;
      
      if (near_contig_end) {
        mc[IGNORE] = "CONTIG_END";
      }
    }
    
    // JC evidence
    vector<gd_entry_type> jc_types = make_vector<gd_entry_type>(JC);
    diff_entry_list_t jc = gd.get_list(jc_types);
    
    for (diff_entry_list_t::iterator jc_it=jc.begin(); jc_it!=jc.end(); jc_it++)
    {
      cDiffEntry& j = **jc_it;
      //cout << j.as_string() << endl;
      
      if (!ref_seq_info[j[SIDE_1_SEQ_ID]].m_is_contig)
        continue;
      if (!ref_seq_info[j[SIDE_2_SEQ_ID]].m_is_contig)
        continue;
      
      int32_t side_1_seq_length = ref_seq_info.get_sequence_length(j[SIDE_1_SEQ_ID]);
      int32_t side_2_seq_length = ref_seq_info.get_sequence_length(j[SIDE_2_SEQ_ID]);

      bool side_1_near_contig_end = false;
      if (n(j[SIDE_1_POSITION]) <= contig_end_distance)
        side_1_near_contig_end = true;
          
      if (n(j[SIDE_1_POSITION]) >= side_1_seq_length-contig_end_distance)
        side_1_near_contig_end = true;

      bool side_2_near_contig_end = false;
      if (n(j[SIDE_2_POSITION]) <= contig_end_distance)
        side_2_near_contig_end = true;
          
      if (n(j[SIDE_2_POSITION]) >= side_2_seq_length-contig_end_distance)
        side_2_near_contig_end = true;
      
      if (side_1_near_contig_end && side_2_near_contig_end) {
        j[IGNORE] = "CONTIG_END";
      }
    }

    // SC evidence
    diff_entry_list_t sc = gd.get_list(make_vector<gd_entry_type>(SC));

    for (diff_entry_list_t::iterator sc_it = sc.begin(); sc_it != sc.end(); sc_it++) {
      cDiffEntry& s = **sc_it;

      if (!ref_seq_info[s[SEQ_ID]].m_is_contig)
        continue;

      int32_t seq_length = ref_seq_info.get_sequence_length(s[SEQ_ID]);

      bool near_contig_end = false;
      if (n(s[POSITION]) <= contig_end_distance)
        near_contig_end = true;

      if (n(s[POSITION]) >= seq_length - contig_end_distance)
        near_contig_end = true;

      if (near_contig_end) {
        s[IGNORE] = "CONTIG_END";
      }
    }

  }


  void MutationPredictor::remove_mutations_near_contig_ends(Settings& settings, Summary& summary, cGenomeDiff& gd)
  {
    
    (void) summary;
    
    int32_t contig_end_distance = settings.ignore_within_this_multiple_of_average_read_length_of_contig_end
      * static_cast<int32_t>(round(summary.sequence_conversion.read_length_avg));
    
    diff_entry_list_t* mutable_entry_list = gd.get_mutable_list_ptr();

    diff_entry_list_t::iterator it=mutable_entry_list->begin();
    bool advance_it = true;
    
    while (it != mutable_entry_list->end())
    {
      advance_it = true;

      cDiffEntry& de = **it;
      if (de.is_mutation()) {
        
        if (ref_seq_info[de[SEQ_ID]].m_is_contig) {
          
          int32_t seq_length = ref_seq_info.get_sequence_length(de[SEQ_ID]);
          
          bool near_contig_end = false;
          if (de.get_reference_coordinate_start() <= contig_end_distance)
            near_contig_end = true;
          
          if (de.get_reference_coordinate_end() >= seq_length-contig_end_distance)
            near_contig_end = true;
          
          
          if (near_contig_end) {
            it = gd.remove(it);
            advance_it = false;
          }
        }
      }

      if (advance_it) ++it;
    }

  }


  // Remove SC (soft-clipping) evidence that is already explained by nearby
  // evidence. For an SC entry at POSITION clipped in direction STRAND, look at the
  // reference "clip window" of soft_clipping_minimum_bases bases immediately in the
  // clip direction (POSITION+1..POSITION+N for STRAND>0, POSITION-N..POSITION-1 for
  // STRAND<0). If any predicted mutation OR rejected RA evidence on the same seq_id
  // overlaps that window, the soft-clipping is redundant and the SC entry is
  // removed entirely. Rejected RA is included because a variant was detected at
  // that position but not promoted to a mutation, which still explains the clip.
  void MutationPredictor::remove_soft_clipping_near_mutations(Settings& settings, cGenomeDiff& gd)
  {
    const int32_t N = static_cast<int32_t>(settings.soft_clipping_minimum_bases);

    // Collect the reference spans that can explain a soft-clipping event:
    // predicted mutations plus rejected RA evidence. Copies of the pointer lists
    // stay valid while we erase SC entries from the underlying entry list below.
    struct explaining_span { string seq_id; int32_t start; int32_t end; };
    vector<explaining_span> spans;

    diff_entry_list_t muts = gd.mutation_list();
    for (diff_entry_list_t::iterator mit = muts.begin(); mit != muts.end(); mit++) {
      cDiffEntry& de = **mit;
      explaining_span span;
      span.seq_id = de[SEQ_ID];
      span.start  = de.get_reference_coordinate_start().get_position();
      span.end    = de.get_reference_coordinate_end().get_position();
      spans.push_back(span);
    }

    diff_entry_list_t ra = gd.get_list(make_vector<gd_entry_type>(RA));
    for (diff_entry_list_t::iterator rit = ra.begin(); rit != ra.end(); rit++) {
      cDiffEntry& de = **rit;
      if (!de.is_rejected()) continue;   // only standard "reject" field counts
      explaining_span span;
      span.seq_id = de[SEQ_ID];
      span.start  = de.get_reference_coordinate_start().get_position();
      span.end    = de.get_reference_coordinate_end().get_position();
      spans.push_back(span);
    }

    diff_entry_list_t* mutable_entry_list = gd.get_mutable_list_ptr();

    diff_entry_list_t::iterator it = mutable_entry_list->begin();
    while (it != mutable_entry_list->end())
    {
      cDiffEntry& sc = **it;
      if (sc._type != SC) { ++it; continue; }

      const string& seq_id = sc[SEQ_ID];
      int32_t position = n(sc[POSITION]);
      int32_t strand = n(sc[STRAND]);

      int32_t window_start, window_end;
      if (strand > 0) {
        window_start = position + 1;
        window_end   = position + N;
      } else {
        window_start = position - N;
        window_end   = position - 1;
      }

      bool explained = false;
      for (vector<explaining_span>::iterator sp = spans.begin(); sp != spans.end(); sp++) {
        if (sp->seq_id != seq_id) continue;
        if (sp->start <= window_end && sp->end >= window_start) {
          explained = true;
          break;
        }
      }

      if (explained) {
        it = gd.remove(it);
      } else {
        ++it;
      }
    }
  }




  // Normalizes INS/DEL mutations by shifting them to the highest coordinates possible
  // by repeat units.
  
  // Adds additional fields "repeat_length", "repeat_ref_copies", and "repeat_new_copies"
  // to any repeat that is above a certain threshold length in the original sequence (by default 6 bp).
  
  // There are some difficulties dealing with deletions that have been split into two
  // e.g. deletion of one G and then another G from a larger G repeat, for now we skip those!
  // Insertions don't have the same problem. We could correct the final repeat number and
  // so on, but it's not obvious how to do this: set it to the final state for all mutations)
  // or set it to the state after each mutation and make that the original state for the next one.
  
  // * 2016-11-26
  // Fixed incorrect overlapping mutations after shifts possible in breseq output
  //
  // For example CCAAAAGT => CCAAATGT will be predicted as two mutations (DEL and SNP) initially
  // like this: CC(-A)AA(A=>T)GT. These two events will be shifted to the same position if we are
  // not careful about looking at the next mutation in the list. The best annotation is a SUB
  // like this: CCAA(AA=>T)GT
  
  // * 2018-08-01
  // There are a lot of times in polymorphism mode where this can go awry if we shift mutations
  // (The nsAA phage samples break the code when there are lots of insertions and also substitutions
  // in the same homopolymer repeat, for example. The same base can be both deleted and substituted!)
  // The only safe way to deal with this will be local re-alignment. For now removing shifts in
  // polymorphism mode. We *do* still annotate them as being in repeats.
  
  void MutationPredictor::normalize_and_annotate_tandem_repeat_mutations(Settings& settings, Summary& summary, cGenomeDiff& gd)
  { 
    (void)settings;
    (void) summary;

    // Sort to absolutely guarantee reproducibility
    // and having mutations in the proper order for merges
    gd.sort();
    
    uint32_t minimum_tandem_repeat_length = 5;
    uint32_t minimum_repeat_ref_copies = 2;
    uint32_t maximum_repeat_sequence_length_to_show = 40;
    
    diff_entry_list_t test_muts = gd.mutation_list();
    
    // Add additional fields for INS or DEL mutations that are in tandem repeats
    
    // We iterate in reverse order, because our normalize functions are moving
    // mutations to the HIGHEST coordinates possible when they are in tandem repeats.
    const cDiffEntry* last_mut(NULL);
    
    for(diff_entry_list_t::reverse_iterator it = test_muts.rbegin(); it != test_muts.rend(); it++) {
      
      cDiffEntry& mut = **it;
      
      // If we are annotated as 'within', do not shift the coordinates,
      // because the reference sequence is actually different by the time we are applied
      
      // We are still potentially in danger of doing the wrong thing here,
      // because a mutation could be applied only after one with the 'before' tag, making the shift
      // incorrect. So, setting 'no_normalize' is an out that can be used.
      
      if (mut.entry_exists("within") || mut.entry_exists("no_normalize") ) goto next_mutation;
      
      if (mut._type == INS) {
        
        int32_t size = mut["new_seq"].size();
        int32_t position = from_string<int32_t>(mut["position"]);
        int32_t insert_position = mut.entry_exists("insert_position") ? from_string<int32_t>(mut["insert_position"]) : 1;

        string mutation_sequence = mut["new_seq"];
        uint32_t repeat_unit_size(0);
        string repeat_unit_sequence;
        find_repeat_unit(mutation_sequence, repeat_unit_size, repeat_unit_sequence);
        
        uint32_t num_repeat_units = size / repeat_unit_size;
        
        // Note shift to +1 to get to where the first unit of a repeat would be for INS
        uint32_t original_num_repeat_units = find_original_num_repeat_units(ref_seq_info[mut["seq_id"]], position+1, repeat_unit_sequence);
        
        // Begin consensus mode shifting of coordinates ---->
        if (!settings.polymorphism_prediction) {
          normalizeINSposition(ref_seq_info[mut["seq_id"]], mut, repeat_unit_sequence);
          int32_t new_position = from_string<int32_t>(mut["position"]);
          string new_mutation_sequence = mut["new_seq"];
          int32_t new_insert_position = mut.entry_exists("insert_position") ? from_string<int32_t>(mut["insert_position"]) : 1;
          
          // Did we get shifted into the position of the next mutation? Then back off
          // Note: We don't do this with converted AMPs as this creates problems (an insertion within them can shift their position)
          if (last_mut && (mut[SEQ_ID] == last_mut->get(SEQ_ID)) && (mut.get_reference_coordinate_end() >= last_mut->get_reference_coordinate_start()) && !mut.entry_exists("_was_AMP")) {
            // The position of this insert mutation should be one before the mutation, unless it is another INS,
            // In the INS case, we need to properly update all of the insert positions
            
            if (last_mut->_type == INS) {
              
              // Don't directly change the position in here... it gets assigned afterward
              int32_t assign_insert_position = 1;
              new_position = n(last_mut->get(POSITION));
              new_insert_position = assign_insert_position;
              
              diff_entry_list_t::reverse_iterator it_ins = it;
              it_ins--;
              cDiffEntry* ins_mut = (*it_ins).get();
              while (   (ins_mut->_type == INS)
                     && (n((*ins_mut)["position"]) == new_position)
                    )
              {
                assign_insert_position++;
                (*ins_mut)["insert_position"] = s(assign_insert_position);
                it_ins--;
                if (it_ins == test_muts.rend()) break;
                ins_mut = (*it_ins).get();
              }
              
            } else {
              new_position = n(last_mut->get(POSITION)) - 1;
            }
          }
          
          // repeat info may have changed, so reload
          position = new_position;
          insert_position = new_insert_position;
          mutation_sequence = new_mutation_sequence;
          find_repeat_unit(mutation_sequence, repeat_unit_size, repeat_unit_sequence);
          num_repeat_units = size / repeat_unit_size;
          
          // Note shift to +1 to get to where the first unit of a repeat would be for INS
          original_num_repeat_units = find_original_num_repeat_units(ref_seq_info[mut["seq_id"]], position+1, repeat_unit_sequence);

          // save normalized position even if we aren't a repeat
          mut["position"] = to_string<int32_t>(position);
          mut["insert_position"] = to_string<int32_t>(insert_position);
          mut["new_seq"] = mutation_sequence;
          
        } // <-------  End of consensus mode shifting of positions
        
        if (original_num_repeat_units * repeat_unit_size < minimum_tandem_repeat_length)
          goto next_mutation;
        
        if (original_num_repeat_units + num_repeat_units < minimum_repeat_ref_copies)
          goto next_mutation;
        
        mut.erase("repeat_seq");
        if (repeat_unit_size <= maximum_repeat_sequence_length_to_show)
          mut["repeat_seq"] = repeat_unit_sequence;
        mut["repeat_length"] = s(repeat_unit_size);
        mut["repeat_ref_copies"] = s(original_num_repeat_units);
        mut["repeat_new_copies"] = s(original_num_repeat_units + num_repeat_units);
      }
      
      else if (mut._type == DEL) {
        
        int32_t size = from_string<int32_t>(mut["size"]);
        
        int32_t position = from_string<int32_t>(mut["position"]);
        
        string mutation_sequence = ref_seq_info.get_sequence_1_start_size(mut["seq_id"], position, size);
        
        uint32_t repeat_unit_size;
        string repeat_unit_sequence;
        find_repeat_unit(mutation_sequence, repeat_unit_size, repeat_unit_sequence);
        
        find_repeat_unit(mutation_sequence, repeat_unit_size, repeat_unit_sequence);
        uint32_t num_repeat_units = size / repeat_unit_size;
        
        // Note shift to +size to get to the spot immediately past all reference repeats
        uint32_t original_num_repeat_units = find_original_num_repeat_units(ref_seq_info[mut["seq_id"]], position+size, repeat_unit_sequence);
        
        // Begin consensus mode shifting of coordinates ---->
        if (!settings.polymorphism_prediction) {
        
          normalizeDELposition(ref_seq_info[mut["seq_id"]], mut, repeat_unit_sequence);
          
          // Did we get shifted into the position of the next mutation? Then back off.
          if (last_mut && !gd.applied_before_id(last_mut->_id, mut._id) && (mut[SEQ_ID] == last_mut->get(SEQ_ID)) && (mut.get_reference_coordinate_end() >= last_mut->get_reference_coordinate_start())) {
            // The position of this insert mutation should be as many bases as the
            // deletion is long before the next mutation
            mut["position"] = s(n(last_mut->get("position")) - size);
          }
          
          // Normalize may actually change the sequence used for the repeat... so call again here.
          position = from_string<int32_t>(mut["position"]);
          mutation_sequence = ref_seq_info.get_sequence_1_start_size(mut["seq_id"], position, size);

          find_repeat_unit(mutation_sequence, repeat_unit_size, repeat_unit_sequence);
          num_repeat_units = size / repeat_unit_size;
          
          // Note shift to +size to get to the spot immediately past all reference repeats
          original_num_repeat_units = find_original_num_repeat_units(ref_seq_info[mut["seq_id"]], position+size, repeat_unit_sequence);
        } // <-------  End of consensus mode shifting of positions

        if (original_num_repeat_units * repeat_unit_size < minimum_tandem_repeat_length)
          goto next_mutation;
        
        if (original_num_repeat_units < minimum_repeat_ref_copies)
          goto next_mutation;
        
        // We only count if there is at least one unit remaining...
        if (original_num_repeat_units - num_repeat_units == 0)
          goto next_mutation;
        
        mut.erase("repeat_seq");
        if (repeat_unit_size <= maximum_repeat_sequence_length_to_show)
          mut["repeat_seq"] = repeat_unit_sequence;
        mut["repeat_length"] = s(repeat_unit_size);
        mut["repeat_ref_copies"] = s(original_num_repeat_units);
        mut["repeat_new_copies"] = s(original_num_repeat_units - num_repeat_units);
      }
      
      // always save this mutation as the last one
      next_mutation: {
        last_mut = &**it;
      }
    }
  }
  
  void MutationPredictor::normalize_INS_to_AMP(Settings& settings, Summary& summary, cGenomeDiff& gd) {
    (void) summary;
    
    int32_t AMP_size_cutoff = settings.size_cutoff_AMP_becomes_INS_DEL_mutation;
    
    // Convert all AMP to INS
    //   so that INS/DEL normalization can take care of them
    diff_entry_list_t mut_list = gd.mutation_list();
    
    // Convert some INS back to AMP
    //   because they are too big to treat as INS
    for(diff_entry_list_t::iterator it=mut_list.begin(); it!=mut_list.end(); it++) {
      cDiffEntry& mut = **it;
      if (mut._type != INS)
        continue;
      
      if (!mut.entry_exists("repeat_length"))
        continue;
      
      int32_t unit_size = from_string<int32_t>(mut["repeat_length"]);
      
      // bail if the repeat length is not long enough
      // even if it could be if we broke into sub-repeats
      if (unit_size <= AMP_size_cutoff)
        continue;
      
      int32_t new_copy_number = from_string<uint32_t>(mut["repeat_new_copies"]);
      
      mut._type = AMP;
      int32_t pos = from_string<uint32_t>(mut[POSITION]) - unit_size * from_string<uint32_t>(mut["repeat_ref_copies"]) + 1;
      // note shift back up of position by one is correct
      // because INS are after this position, but AMP start at this position
      mut["position"] = to_string<int32_t>(pos);
      mut["new_copy_number"] = mut["repeat_new_copies"];
      mut["size"] = mut["repeat_length"];
      
      // delete all the repeat information...
      mut.erase(NEW_SEQ);
      mut.erase("insert_position");
      mut.erase("repeat_length");
      mut.erase("repeat_sequence");
      mut.erase("repeat_new_copies");
      mut.erase("repeat_ref_copies");
    }
  }

  
	/*
	 Title   : predict
	 Usage   : $mp->predict();
	 Function: Predicts mutations from evidence in a GenomeDiff and adds them to it
	 Returns :

	*/
	void MutationPredictor::predict(Settings& settings, Summary& summary, cGenomeDiff& gd)
	{
    bool verbose = false; // for debugging
    
    // Setup
    vector<gd_entry_type> jc_types = make_vector<gd_entry_type>(JC);
    vector<gd_entry_type> mc_types = make_vector<gd_entry_type>(MC);

    
    // This checks for using a bad reference and warns
    {
      vector<gd_entry_type> ev_types = make_vector<gd_entry_type>(MC)(RA)(JC)(CN);
      diff_entry_list_t ev = gd.get_list(ev_types);
      
      ev.remove_if(cDiffEntry::field_exists(REJECT));
      ev.remove_if(cDiffEntry::field_exists(IGNORE));

      // Also remove consensus reject RA, which may be not rejected in CONSENSUS mode
      // There are many of these for MinION data.
      ev.remove_if(cDiffEntry::field_exists(CONSENSUS_REJECT));
      
      uint64_t num_evidence_items = ev.size();
      uint64_t total_ref_seq_length = summary.sequence_conversion.total_reference_sequence_length;
      double maximum_sequence_divergence = static_cast<double>(num_evidence_items)*100/static_cast<double>(total_ref_seq_length);
      if ( (num_evidence_items > 0) && (maximum_sequence_divergence > 0.2) ) {
        WARN("Large number of differences detected between the sample and the reference sequence (" + to_string<uint64_t>(num_evidence_items) + " evidence items, suggesting approximately " +  to_string(formatted_double(maximum_sequence_divergence,2)) + "% sequence divergence). If this is unexpected, check that you are using the closest available reference sequence for this sample (for example, the correct strain of a bacterial species). Mutation prediction can become less accurate with too much divergence from the reference sequence. It may also take a long time to predict mutations and generate output files.");
      }
    }
    

    // Do some filtering if certain options are in play
    if (settings.polymorphism_no_indels || settings.consensus_reject_indel_homopolymer_length || settings.polymorphism_reject_indel_homopolymer_length) {
      cGenomeDiff junction_indel_gd;
      diff_entry_list_t jc = gd.get_list(jc_types);
      jc.remove_if(cDiffEntry::ignored_and_not_user_defined());

      diff_entry_list_t mc = gd.get_list(mc_types);
      mc.remove_if(cDiffEntry::ignored_and_not_user_defined());
      
      filter_JC_indel_evidence(settings, summary, gd, jc, mc);
    }
    
    ignore_evidence_near_contig_ends(settings, summary, gd);
    
		///
		//  Preprocessing of JC evidence
    //  NB: This call is likely redundant in the normal pipeline, but preserved here so that 
    //  predict can be a stand-alone call that is not dependent on other processing.
		///
    
    cerr << "  Preparing junctions..." << endl;
    prepare_junctions(settings, summary, gd);

    ///
    // Create master lists of evidence that will be culled as they are used
    // by functions below.
    ///
    
    // Do not use rejected junctions unless they are user-defined
		diff_entry_list_t jc = gd.get_list(jc_types);
    jc.remove_if(cDiffEntry::ignored_and_not_user_defined());

    // Do not use rejected missing coverage evidence
		diff_entry_list_t mc = gd.get_list(mc_types);
    mc.remove_if(cDiffEntry::ignored_and_not_user_defined());
    mc.remove_if(cDiffEntry::rejected_and_not_user_defined());

    // Copy number evidence, present only under --predict-copy-number. Everything that consumes it
    // below is skipped when this is empty, so a run without the option is completely unaffected.
    diff_entry_list_t cn = gd.get_list(make_vector<gd_entry_type>(CN));

    ///
    // evidence JC + JC = MOB mutation
    ///
    
    cerr << "  Predicting mobile element insertions..." << endl;
    predictJCplusJCtoMOB(settings, summary, gd, jc, mc);
    
    // Don't use rejected junctions for predicting DEL, INS, SUB, *just* MOB
    // because frequencies of rejected ones could be increased by recounting
    // after adjusting for overlap in this case!
    jc.remove_if(cDiffEntry::rejected_and_not_user_defined());
    
		///
		// evidence MC + JC => DEL mutation
		///
    cerr << "  Predicting large deletions..." << endl;
    predictMCplusJCtoDEL(settings, summary, gd, jc, mc);

    ///
    // evidence MC => DEL mutation (deletion between homologous copies, where no JC can exist)
    ///
    // Runs after the JC cases so split-read evidence, which locates a breakpoint to the base, always
    // wins; this only sees the MC items those cases left unexplained.
    cerr << "  Predicting large deletions between homologous sequences..." << endl;
    {
      diff_entry_list_t dp = gd.get_list(make_vector<gd_entry_type>(DP));
      diff_entry_list_t ra_for_crossover = gd.get_list(make_vector<gd_entry_type>(RA));
      predictMCtoDELbyHomology(settings, summary, gd, mc, dp, ra_for_crossover);
    }

    ///
    // evidence JC + CN, and CN + repeat => AMP mutations
    ///
    //
    // Deliberately ahead of the polymorphism filter below: an amplification is usually still growing
    // in the population, so its junction is exactly the sub-consensus kind that filter removes -- as
    // its own comment says. Both are no-ops without CN evidence.
    if (!cn.empty()) {
      cerr << "  Predicting amplifications from junctions and copy number..." << endl;
      predictJCplusCNtoAMP(settings, summary, gd, jc, cn);
      predictCNplusRepeatToAMP(settings, summary, gd, jc, cn);
    }

		///
		// evidence JC => INS, SUB, DEL mutations
    ///
    
    // Do not use polymorphic junctions here in consensus mode
    // -- needs to be relaxed in the future for predicting amplifications!
    // IMPORTANT: they are used above because other evidence implies consensus
    // or predicting the mutation can re-adjust frequencies
    //
    // @JEB 11-12-2015 We don't add a reject reason to these, as they
    // still look like valid junctions even if their frequency is different
    
    if (!settings.polymorphism_prediction) {
      jc.remove_if(cDiffEntry::field_equals(PREDICTION, "polymorphism"));
    }
    cerr << "  Predicting small indels and substitutions from junctions..." << endl;
    predictJCtoINSorSUBorDEL(settings, summary, gd, jc, mc);
    
		///
		// evidence RA => SNP, DEL, INS, SUB
		///
    vector<gd_entry_type> ra_types = make_vector<gd_entry_type>(RA);
    diff_entry_list_t ra = gd.get_list(ra_types);
    ra.remove_if(cDiffEntry::rejected_and_not_user_defined());
    ra.remove_if(cDiffEntry::ignored_and_not_user_defined());
    
    cerr << "  Predicting small indels and substitutions from alignments..." << endl;
    predictRAtoSNPorDELorINSorSUB(settings, summary, gd, ra, mc);
    
    ///
    //  Check for completely deleted reference sequences
    //  we need to be careful and not predict any mutations on these
    //  as they can mess up shifting
    //  @JEB: It would be better to never add these to the GenomeDiff in the first
    //        place, but for now we strip them out at the end, leaving only
    //        the deletion of the entire fragment and no other muts on that reference
    ///
    remove_mutations_on_deleted_reference_sequences(settings, summary, gd);
    
    ///
    // Check for mutations predicted by a JC and by a combination of RA evidence
    ///
     cerr << "  Reconciling mutation predictions..." << endl;
    gd.reconcile_mutations_predicted_two_ways();
    
    cerr << "  Making final adjustments to mutations..." << endl;
    normalize_and_annotate_tandem_repeat_mutations(settings, summary, gd);
    
    // Combine INS/DEL mutations that have been shifted to be adjacent to other mutations.
    // But NOT in polymorphism mode where we are uncertain if they are different mutations!
    if (!settings.polymorphism_prediction) {
      combine_newly_adjacent_mutations(gd);
    }
    
    ///
		// mutation INS => mutation AMP
		///
    normalize_INS_to_AMP(settings, summary, gd);
    
    ///
    // Check for certain kinds of overlap that need 'before' or 'within' fields to resolve
    ///
    if (!settings.polymorphism_prediction) {
      assign_before_within_to_mutations(gd);
    }
    
    ///////////////////////////////////////////////////////
    // Check to be sure the "frequency" field is present //
    // as appropriate in consensus/polymorphism mode.    //
    ///////////////////////////////////////////////////////
    
    {
      diff_entry_list_t check_mut_list = gd.mutation_list();
      for (diff_entry_list_t::iterator it=check_mut_list.begin(); it != check_mut_list.end(); it++) {
        cDiffEntry& de = **it;
        if (settings.polymorphism_prediction) {
          ASSERT(de.entry_exists(FREQUENCY) && !de[FREQUENCY].empty(), "Expected polymorphism field [frequency] not found for mutation.\n" + de.as_string() );
        } else {
          ASSERT(!de.entry_exists(FREQUENCY) || !de[FREQUENCY].empty(), "Field [frequency] not expected in consensus mode for mutation.\n" + de.as_string() );
        }
      }
    }
    
    // Need to remove any mutations that overlap contig ends
    // This catches SUB from JC evidence.
    remove_mutations_near_contig_ends(settings, summary, gd);

    // Remove SC evidence that is redundant because a predicted mutation already
    // explains the reference bases in the clip direction.
    if (settings.predict_soft_clipping) {
      remove_soft_clipping_near_mutations(settings, gd);
    }
    // 32  189  51 bp→32 bp  57.0%  intergenic (–/+224)  – / ← INEPCCPE_03363  –/tRNA‑Ala(tgc)

    
    
    // In consensus mode, polish up the insert_position fields of INS predictions
    // so that if there is only one, then insert_position=1 is removed.
    if (!settings.polymorphism_prediction) {
      diff_entry_list_t ins_list = gd.get_list(make_vector<gd_entry_type>(INS));
      cDiffEntry* last_ins(NULL);
      for (diff_entry_list_t::iterator it=ins_list.begin(); it != ins_list.end(); it++) {
        cDiffEntry* ins = it->get();
        
        if (
            last_ins
            && ((*ins)[POSITION] != (*last_ins)[POSITION])
            && last_ins->entry_exists(INSERT_POSITION)
            && ( n(last_ins->get(INSERT_POSITION)) == 1)
            )
        {
          (*last_ins)["_dont_print_insert_position"] = "1";
          //last_ins->erase(INSERT_POSITION);
        }
        last_ins = ins;
      }
      
      // Catch the last INS mutation
      if (
          last_ins
          && last_ins->entry_exists(INSERT_POSITION)
          && ( n(last_ins->get(INSERT_POSITION)) == 1)
          )
      {
        (*last_ins)["_dont_print_insert_position"] = "1";
        //last_ins->erase(INSERT_POSITION);
      }

    }

    // Drop any discordant-pair (DP) item that a pair-distance (PD) item already describes, BEFORE
    // the DP attachment below, so a superseded DP is never attached to a mutation on its way out.
    remove_DP_superseded_by_PD(settings, summary, gd);

    // Attach a matching discordant-pair (DP) evidence item to each mutation whose JC it matches.
    add_matching_DP_evidence(gd);

    // Same for pair-distance (PD) evidence.
    add_matching_PD_evidence(gd);

    // Same for copy number (CN) evidence, matched by span against a deletion rather than by
    // breakpoint against a junction.
    add_matching_CN_evidence(gd);

    // Also combine DP with MOBs by their unique side (IS sides may be on different copies); a DP that
    // localizes a specific non-consensus IS copy sets the MOB's mob_region.
    combine_DP_with_MOB_by_unique_side(gd);

	}

  // Canonical, side-order-independent key for a junction's two sides (seq_id|position|strand each),
  // so a DP and a JC describing the same breakpoint compare equal regardless of which side is side_1.
  static string junction_side_key(cDiffEntry& e)
  {
    string a = e[SIDE_1_SEQ_ID] + "|" + e[SIDE_1_POSITION] + "|" + e[SIDE_1_STRAND];
    string b = e[SIDE_2_SEQ_ID] + "|" + e[SIDE_2_POSITION] + "|" + e[SIDE_2_STRAND];
    return (a < b) ? (a + "#" + b) : (b + "#" + a);
  }

  // How far apart a DP's unique side and a JC's unique side may sit and still be the same breakpoint.
  // A DP side is placed at the innermost aligned end of its supporting read pairs, so it lands at or
  // just short of the split-read coordinate rather than exactly on it. Same value the MOB unique-side
  // match below has always used.
  static const int32_t kDPJCMatchSlop = 20;

  // Do side `ap` of a and side `bp` of b describe the same breakpoint -- same sequence and strand,
  // position within slop?
  static bool same_junction_side(cDiffEntry& a, const string& ap, cDiffEntry& b, const string& bp, int32_t slop)
  {
    return (a[ap + "_seq_id"] == b[bp + "_seq_id"])
        && (n(a[ap + "_strand"]) == n(b[bp + "_strand"]))
        && (abs(n(a[ap + "_position"]) - n(b[bp + "_position"])) <= slop);
  }

  // Split a junction-like entry (JC or DP) into its unique side and its repeat (IS) side. False unless
  // exactly one side is annotated "repeat" -- with neither or both there is no unique side to match on.
  static bool split_repeat_unique_sides(cDiffEntry& e, string& unique_prefix, string& repeat_prefix)
  {
    string ak1 = e.entry_exists("side_1_annotate_key") ? e["side_1_annotate_key"] : "";
    string ak2 = e.entry_exists("side_2_annotate_key") ? e["side_2_annotate_key"] : "";
    if      (ak1 == "repeat" && ak2 != "repeat") { repeat_prefix = "side_1"; unique_prefix = "side_2"; }
    else if (ak2 == "repeat" && ak1 != "repeat") { repeat_prefix = "side_2"; unique_prefix = "side_1"; }
    else return false;
    return true;
  }

  // An IS-mediated junction, matched by its UNIQUE side alone: each entry has exactly one repeat side
  // and one unique side, and the unique sides agree. The repeat sides are then allowed to sit on
  // different copies of the element -- a DP's IS side is chosen by the per-locus copy vote in
  // resolve_alignments, which has no reason to pick the same copy the split reads did, so requiring
  // both sides would miss every IS-mediated event. This is the same rule combine_DP_with_MOB_by_unique_side
  // already applies to MOB; here it reaches the DEL/AMP/... that an IS-mediated junction also produces.
  // The window stays tight: a DP whose unique side is further off than this is fixed at the source, by
  // the unique-side JC snap in dp_evidence.cpp, not by widening the match here.
  static bool same_breakpoint_by_unique_side(cDiffEntry& j, cDiffEntry& d)
  {
    string ju, jr, du, dr;
    if (!split_repeat_unique_sides(j, ju, jr)) return false;
    if (!split_repeat_unique_sides(d, du, dr)) return false;
    return same_junction_side(j, ju, d, du, kDPJCMatchSlop);
  }

  // Attach a discordant-pair (DP) evidence item as ALSO supporting a mutation when it describes the
  // same breakpoint as a JC that already supports that mutation -- e.g. a DP that snapped onto that
  // junction's split-read coordinates. This makes the mutation show a "DP" evidence link and drops the
  // DP from the "unassigned discordant pair evidence" table, both via the native evidence= machinery
  // (in_evidence_list / filter_used_as_evidence).
  //
  // Two ways to be the same breakpoint:
  //  (1) both sides match exactly (seq_id/position/strand, either side order);
  //  (2) it is IS-mediated -- each entry has exactly one repeat side -- and the UNIQUE sides match.
  //      An IS-mediated DEL is the common case: the JC and the DP agree on the unique boundary, but
  //      their IS sides sit on different copies of the element, so (1) alone never fires.
  void MutationPredictor::add_matching_DP_evidence(cGenomeDiff& gd)
  {
    // Non-rejected DP evidence only (rejected DP live in the marginal table, not the unassigned one,
    // so they are not promoted). Indexed by exact side key for (1); the list is walked for (2).
    map<string, vector<string> > dp_by_key;
    vector<diff_entry_ptr_t> dps;
    diff_entry_list_t dp_list = gd.get_list(make_vector<gd_entry_type>(DP));
    for (diff_entry_list_t::iterator it = dp_list.begin(); it != dp_list.end(); it++) {
      cDiffEntry& dp = **it;
      if (dp.entry_exists(REJECT)) continue;
      dp_by_key[junction_side_key(dp)].push_back(dp._id);
      dps.push_back(*it);
    }
    if (dps.empty()) return;

    diff_entry_list_t muts = gd.mutation_list();
    for (diff_entry_list_t::iterator it = muts.begin(); it != muts.end(); it++) {
      cDiffEntry& mut = **it;
      // Snapshot the current evidence ids (we append to mut._evidence inside the loop).
      vector<string> evidence_now = mut._evidence;
      for (vector<string>::iterator ev = evidence_now.begin(); ev != evidence_now.end(); ev++) {
        diff_entry_ptr_t e = gd.find_by_id(*ev);
        if (e.get() == NULL || e->_type != JC) continue;

        // (1) exact both-sides match
        map<string, vector<string> >::const_iterator m = dp_by_key.find(junction_side_key(*e));
        if (m != dp_by_key.end()) {
          for (vector<string>::const_iterator dpid = m->second.begin(); dpid != m->second.end(); dpid++) {
            if (find(mut._evidence.begin(), mut._evidence.end(), *dpid) == mut._evidence.end())
              mut._evidence.push_back(*dpid);
          }
        }

        // (2) IS-mediated: unique sides match, IS sides may be different copies
        for (vector<diff_entry_ptr_t>::iterator d = dps.begin(); d != dps.end(); d++) {
          if (!same_breakpoint_by_unique_side(*e, **d)) continue;
          if (find(mut._evidence.begin(), mut._evidence.end(), (*d)->_id) == mut._evidence.end())
            mut._evidence.push_back((*d)->_id);
        }
      }
    }
  }

  // Attach a pair-distance (PD) evidence item as ALSO supporting a mutation when it describes the
  // same breakpoint as a JC that already supports it -- the direct analogue of add_matching_DP_evidence
  // above, and it matters most for the case PD snapped onto that very junction, where the two agree
  // to the base.
  void MutationPredictor::add_matching_PD_evidence(cGenomeDiff& gd)
  {
    map<string, vector<string> > pd_by_key;
    vector<diff_entry_ptr_t> pds;
    diff_entry_list_t pd_list = gd.get_list(make_vector<gd_entry_type>(PD));
    for (diff_entry_list_t::iterator it = pd_list.begin(); it != pd_list.end(); it++) {
      cDiffEntry& pd = **it;
      if (pd.entry_exists(REJECT)) continue;
      pd_by_key[junction_side_key(pd)].push_back(pd._id);
      pds.push_back(*it);
    }
    if (pds.empty()) return;

    diff_entry_list_t muts = gd.mutation_list();
    for (diff_entry_list_t::iterator it = muts.begin(); it != muts.end(); it++) {
      cDiffEntry& mut = **it;
      vector<string> evidence_now = mut._evidence;
      for (vector<string>::iterator ev = evidence_now.begin(); ev != evidence_now.end(); ev++) {
        diff_entry_ptr_t e = gd.find_by_id(*ev);
        if (e.get() == NULL || e->_type != JC) continue;

        map<string, vector<string> >::const_iterator m = pd_by_key.find(junction_side_key(*e));
        if (m == pd_by_key.end()) continue;
        for (vector<string>::const_iterator pdid = m->second.begin(); pdid != m->second.end(); pdid++) {
          if (find(mut._evidence.begin(), mut._evidence.end(), *pdid) == mut._evidence.end())
            mut._evidence.push_back(*pdid);
        }
      }
    }
  }

  // Attach a copy number (CN) item as ALSO supporting a deletion that already explains it.
  //
  // Prediction never uses CN to CALL a deletion -- that always comes from MC/JC/RA -- but a region
  // the HMM called at copy number zero across a deletion breseq already predicted is describing that
  // same event. Leaving it unattached reports it a second time under "unassigned copy number
  // evidence", as though it were an observation nothing accounted for.
  //
  // The match is made in TILES, not in percent. A CN boundary can only land on a multiple of CNery's
  // window (tile_size), while a deletion's boundaries come from base-resolution MC/JC/RA evidence, so
  // the two never agree exactly and the disagreement is bounded by the tiling rather than by the size
  // of the event. Percent overlap is the wrong statistic here: it makes a short deletion look like a
  // poor match purely because the tile is a large fraction of it. In the LTEE data a 407 bp deletion
  // pairs with a CN region that overhangs it by 81 bp -- under half a tile, but only 80% overlap,
  // while every large deletion trivially clears 99%.
  //
  // So the CN has to lie INSIDE the deletion, allowing one tile of slop at each end. Deliberately
  // one-directional: the deletion does not have to be mostly covered, which is what lets one long
  // deletion collect several CN regions when the HMM breaks it up, and lets a CN region sitting well
  // inside a longer deletion attach at all.
  //
  // Only copy number zero qualifies. One copy is the baseline and is never emitted; two or more is
  // the AMP paths' input (predictJCplusCNtoAMP and predictCNplusRepeatToAMP both require copies >= 2)
  // -- which is also why no CN item can end up cited by both an AMP and a DEL.
  //
  // Deletions at or below one tile are skipped: CNery cannot resolve an event smaller than its own
  // window, so a CN region overlapping a 1 bp deletion is a coincidence of position and not a
  // measurement of it. Without this a large amplified region would attach itself to any incidental
  // single-base deletion inside it.
  void MutationPredictor::add_matching_CN_evidence(cGenomeDiff& gd)
  {
    // The same two filters output.cpp applies when building the unassigned copy number table, so an
    // item hidden from that table here is exactly one it would otherwise have shown. Neither is set
    // on CN by the current pipeline; they are here so the two stay in step if that changes.
    vector<diff_entry_ptr_t> cns;
    diff_entry_list_t cn_list = gd.get_list(make_vector<gd_entry_type>(CN));
    for (diff_entry_list_t::iterator it = cn_list.begin(); it != cn_list.end(); it++) {
      cDiffEntry& cn = **it;
      if (cn.entry_exists(REJECT)) continue;
      if (cn.entry_exists(IGNORE)) continue;
      // copy_number is a positional field on CN, so it is always present.
      if (from_string<int32_t>(cn[COPY_NUMBER]) != 0) continue;
      cns.push_back(*it);
    }
    // Leaves a run without --predict-copy-number bit-identical.
    if (cns.empty()) return;

    diff_entry_list_t dels = gd.get_list(make_vector<gd_entry_type>(DEL));
    for (diff_entry_list_t::iterator it = dels.begin(); it != dels.end(); it++) {
      cDiffEntry& mut = **it;

      int32_t del_size  = from_string<int32_t>(mut[SIZE]);
      int32_t del_start = from_string<int32_t>(mut[POSITION]);
      int32_t del_end   = del_start + del_size - 1;

      for (vector<diff_entry_ptr_t>::iterator cn_it = cns.begin(); cn_it != cns.end(); cn_it++) {
        cDiffEntry& cn = **cn_it;
        if (cn[SEQ_ID] != mut[SEQ_ID]) continue;

        // Read start/end directly -- get_reference_coordinate_start/end() are not implemented for CN.
        int32_t cn_start = from_string<int32_t>(cn[START]);
        int32_t cn_end   = from_string<int32_t>(cn[END]);

        // tile_size is written by CNEvidence on every CN item it creates, but it is not one of the
        // type's required fields, so a CN read back from a hand-made genome diff may not have it.
        // Falling back to zero slop degrades this to strict containment rather than asserting.
        int32_t tile = cn.entry_exists("tile_size") ? from_string<int32_t>(cn["tile_size"]) : 0;

        if (del_size <= tile) continue;               // below what CNery can resolve
        if (cn_start < del_start - tile) continue;    // hangs off the start by more than one tile
        if (cn_end   > del_end   + tile) continue;    // ... or off the end

        // A CN region abutting a boundary can satisfy the two tests above while overlapping the
        // deletion barely or not at all, so require most of it to actually be inside.
        int32_t overlap = min(del_end, cn_end) - max(del_start, cn_start) + 1;
        if (2 * overlap < cn_end - cn_start + 1) continue;

        // Append, never insert: output.cpp reads in_evidence_list().front() for INS mutations.
        if (find(mut._evidence.begin(), mut._evidence.end(), cn._id) == mut._evidence.end())
          mut._evidence.push_back(cn._id);
      }
    }
  }

  // Remove every discordant-pair (DP) item that a pair-distance (PD) item already describes.
  //
  // The two see the same event from opposite ends of the same distribution. When a deletion is large
  // enough, the tail of its spanning pairs crosses the discordant cutoff and DP seeds on that tail;
  // PD meanwhile uses the whole population. So the DP is built from a handful of the most extreme
  // pairs while the PD is built from all of them, and the PD is better in every respect that matters
  // -- it localizes from more evidence, it can snap to a junction, and it reports the event's size,
  // which DP cannot. Keeping both would list one event twice, with the weaker entry's coordinates
  // disagreeing with the stronger one's.
  //
  // The match window is deliberately much wider than the slop used for JC comparisons. A DP side is
  // placed at the innermost aligned edge of its few supporting pairs, which sits short of the true
  // breakpoint by however much of the fragment those pairs left unsequenced -- routinely tens to
  // hundreds of bases, and by construction never past it. One paired-mapping distance is the natural
  // scale of that error, and two distinct events closer together than that are not separable by pair
  // evidence in any case.
  // The CN entry whose span covers most of [start, end], or NULL. CN boundaries are tiled, so an
  // interval derived from base-resolution evidence never lines up with them exactly; requiring the
  // majority of the interval to be inside is the useful test rather than containment.
  static cDiffEntry* cn_covering(diff_entry_list_t& cn, const string& seq_id, int32_t start, int32_t end)
  {
    cDiffEntry* best = NULL;
    int32_t best_overlap = 0;
    int32_t length = end - start + 1;
    if (length <= 0) return NULL;

    for (diff_entry_list_t::iterator it = cn.begin(); it != cn.end(); it++) {
      cDiffEntry& c = **it;
      if (c[SEQ_ID] != seq_id) continue;
      int32_t overlap = min(end, from_string<int32_t>(c[END]))
                      - max(start, from_string<int32_t>(c[START])) + 1;
      if (overlap <= 0) continue;
      if (2 * overlap < length) continue;
      if (overlap > best_overlap) { best_overlap = overlap; best = &c; }
    }
    return best;
  }

  // Orientation, in reference coordinates, of the NEW element copy a junction reports. This is the
  // same quantity a MOB's "strand" field carries and the same one cFeatureLocation::get_strand()
  // gives for an annotated copy, so all three are directly comparable.
  //
  // Both prefixes carry their trailing underscore ("side_1_"). The triple product is not obvious:
  // the element's own strand says which way the copy the junction READ THROUGH sits, and the two
  // junction side strands then carry that over to where the new copy landed. predictMCplusJCtoDEL's
  // two-junction case fills a MOB's strand with exactly this, and so do both callers below, so an
  // orientation test and the MOB it is tested against cannot disagree about which way round the
  // element sits.
  static int32_t new_repeat_copy_strand(cDiffEntry& j, const string& unique_side, const string& is_side)
  {
    return -( n(j[is_side + "strand"])
            * n(j["_" + is_side + "is_strand"])
            * n(j[unique_side + "strand"]) );
  }

  // Build the MOB implied by a single IS junction, with the target-site duplication left undetermined.
  //
  // One junction fixes where the element landed and which way round it is, but the duplication size
  // needs the junction on the OTHER side of the element to measure against -- and at an amplification
  // boundary that junction lies inside a copy of the element, so it is redundant and gets discarded
  // before prediction. Same situation predictMCplusJCtoDEL's two-junction IS case is in, and the same
  // answer: emit duplication_size 0 flagged indeterminate rather than invent a number.
  //
  // Returns false if the junction was never annotated with an IS interval.
  static bool mob_from_single_IS_junction(cDiffEntry& j, bool polymorphism_prediction, cDiffEntry& mob)
  {
    if (!j.entry_exists("_is_interval")) return false;

    const string& is_side     = j["_is_interval"];
    const string& unique_side = j["_unique_interval"];

    int32_t is_strand = new_repeat_copy_strand(j, unique_side + "_", is_side + "_");

    mob._type = MOB;
    mob._evidence.push_back(j._id);
    mob
      (SEQ_ID,                           j[unique_side + "_seq_id"])
      (REPEAT_NAME,                      j["_" + is_side + "_is_name"])
      (STRAND,                           to_string<int32_t>(is_strand))
      ("duplication_size",               "0")
      ("indeterminate_duplication_size", "1")
      (POSITION,                         j[unique_side + "_position"]);

    if (polymorphism_prediction) mob[FREQUENCY] = "1";
    return true;
  }

  // Nearest annotated copy of one repeat family to `position`, or NULL if none is within `slop`.
  // Filtering by family is the point: the closest element to a boundary is often not the family the
  // junction names (at one real boundary an IS150 sits 423 bp away and the IS186 actually involved
  // 2104 bp away).
  //
  // `direction` -- -1 for a copy expected at or before `position`, +1 at or after -- is a PREFERENCE
  // that breaks ties, not a filter. It used to reject every copy on the far side outright, which
  // contradicts the reason `slop` exists at all: `position` here is a CN boundary, and the HMM places
  // those only to within a few of its own tiles, so the element that really bounds the amplified unit
  // can land either side of where the boundary was drawn. That is not hypothetical -- when CNery
  // halved its default window (200 bp overlapping -> 100 bp non-overlapping) one LTEE clone's CN
  // boundary moved 1.7 kb short of the IS copy bounding it, well inside `slop`, and the hard test
  // silently dropped a 26.6 kb amplification that had been called for years.
  //
  // Distance is measured to the element itself rather than to one chosen edge, so a copy on the far
  // side is not penalized by its own length: measuring a downstream element from its far edge added
  // up to a full element's width to the distance and could push a genuine match past `slop`.
  //
  // `required_strand` (+1/-1, or 0 to accept either) is the orientation the caller needs. Filtering
  // HERE rather than rejecting the answer afterwards is deliberate: it returns the nearest copy that
  // could actually have produced the event, instead of returning the nearest copy and then refusing
  // an event a slightly more distant copy explains.
  static cFeatureLocation* closest_repeat_of_family(cAnnotatedSequence& seq, const string& family,
                                                    int32_t position, int32_t direction, int32_t slop,
                                                    int32_t required_strand = 0)
  {
    cFeatureLocation* best = NULL;
    int32_t best_d = slop + 1;
    bool best_preferred = false;

    for (cFeatureLocationList::iterator it = seq.m_repeat_locations.begin(); it != seq.m_repeat_locations.end(); it++) {
      cFeatureLocation& r = *it;
      if ((*r.get_feature())["name"] != family) continue;
      if (required_strand && (r.get_strand() != required_strand)) continue;

      // Distance from the point to the interval; zero when the point is inside the element.
      int32_t d;
      if (position < r.get_start_1())    d = r.get_start_1() - position;
      else if (position > r.get_end_1()) d = position - r.get_end_1();
      else                               d = 0;
      if (d > slop) continue;

      bool preferred = (d == 0)
                    || ((direction < 0) ? (r.get_end_1() <= position) : (r.get_start_1() >= position));

      // Nearest wins; a copy on the expected side breaks a tie.
      if ((d < best_d) || ((d == best_d) && preferred && !best_preferred)) {
        best_d = d; best = &r; best_preferred = preferred;
      }
    }
    return best;
  }

  static bool mc_fragment_start_lt(const cDiffEntry* a, const cDiffEntry* b)
  {
    return from_string<int32_t>((*const_cast<cDiffEntry*>(a))[START])
         < from_string<int32_t>((*const_cast<cDiffEntry*>(b))[START]);
  }

  void MutationPredictor::predictJCplusCNtoAMP(Settings& settings, Summary& summary, cGenomeDiff& gd,
                                               diff_entry_list_t& jc, diff_entry_list_t& cn)
  {
    (void)settings; (void)summary;
    if (cn.empty()) return;

    for (diff_entry_list_t::iterator jc_it = jc.begin(); jc_it != jc.end(); ) {
      cDiffEntry& j = **jc_it;

      // Same geometry predictJCtoINSorSUBorDEL derives, kept deliberately identical so the two agree
      // on what "points back upstream" means.
      if (j["side_1_seq_id"] != j["side_2_seq_id"]) { jc_it++; continue; }
      if ((n(j["side_1_redundant"]) == 1) || (n(j["side_2_redundant"]) == 1)) { jc_it++; continue; }

      int32_t side_1_strand = n(j["side_1_strand"]);
      int32_t side_2_strand = n(j["side_2_strand"]);
      if (side_1_strand == side_2_strand) { jc_it++; continue; }

      string seq_id = j["side_1_seq_id"];
      int32_t side_1_position = n(j["side_1_position"]);
      int32_t side_2_position = n(j["side_2_position"]);

      if (side_2_position < side_1_position) {
        swap(side_1_position, side_2_position);
        swap(side_1_strand, side_2_strand);
      }
      // Lower side on the minus strand is a deletion, not a duplication.
      if (side_1_strand == -1) { jc_it++; continue; }

      int32_t size = side_2_position - side_1_position + 1;

      // Below the cutoff this is predictJCtoINSorSUBorDEL's INS to make, not ours. Taking only the
      // junctions it drops is what keeps that path's output unchanged.
      if (size <= kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation) { jc_it++; continue; }

      cDiffEntry* cn_item = cn_covering(cn, seq_id, side_1_position, side_2_position);
      if (cn_item == NULL) { jc_it++; continue; }
      int32_t copies = from_string<int32_t>((*cn_item)[COPY_NUMBER]);
      if (copies < 2) { jc_it++; continue; }

      // Right-shift to the highest coordinates that describe the same duplication, which is what
      // gdtools NORMALIZE enforces so that two runs reporting one event compare equal. A tandem block
      // [p, p+size-1] is indistinguishable from [p+1, p+size] exactly when ref[p] == ref[p+size], so
      // walk while that holds. Junction coordinates are not normalized to begin with: on the LTEE
      // rnk-citT amplification this moves 625889 to 625890, which is both where NORMALIZE lands it
      // and what the curated genome diff records.
      //
      // Capped at one full period -- a perfectly repeating block would otherwise walk a whole unit
      // and describe the same duplication all over again.
      int32_t position = side_1_position;
      int32_t seq_length = static_cast<int32_t>(ref_seq_info[seq_id].m_length);
      for (int32_t moved = 0; moved < size; moved++) {
        if (position + size > seq_length) break;
        if (ref_seq_info.get_sequence_1(seq_id, position, position)
            != ref_seq_info.get_sequence_1(seq_id, position + size, position + size)) break;
        position++;
      }

      cDiffEntry mut(AMP);
      mut._evidence = make_vector<string>(j._id)(cn_item->_id);
      mut[SEQ_ID] = seq_id;
      mut[POSITION] = to_string<int32_t>(position);
      mut[SIZE] = to_string<int32_t>(size);
      mut[NEW_COPY_NUMBER] = to_string<int32_t>(copies);
      gd.add(mut);

      jc_it = jc.erase(jc_it);
    }
  }

  void MutationPredictor::predictCNplusRepeatToAMP(Settings& settings, Summary& summary, cGenomeDiff& gd,
                                                   diff_entry_list_t& jc, diff_entry_list_t& cn)
  {
    (void)settings; (void)summary;
    if (cn.empty()) return;

    // How far from a CN boundary the pre-existing copy may sit. The HMM only resolves to its tile
    // size and the amplified unit does not have to end exactly on the element, so this is loose:
    // measured distances on real events are 223, 423 and 2104 bp.
    const int32_t kRepeatBoundarySlop = 3000;
    // How far a MOB or junction may sit from a CN boundary and still be that boundary. The anchor
    // itself is base-resolution; what is loose is the CN edge, which the HMM only places to within a
    // few of its tiles. Measured on real events: 26, 73, 755 and 1195 bp. The family match and the
    // requirement that the partner side lie in a repeat are what actually keep this specific.
    const int32_t kAnchorSlop = 2000;

    diff_entry_list_t mob_list = gd.get_list(make_vector<gd_entry_type>(MOB));

    for (diff_entry_list_t::iterator ci = cn.begin(); ci != cn.end(); ci++) {
      cDiffEntry& c = **ci;
      int32_t copies = from_string<int32_t>(c[COPY_NUMBER]);
      if (copies < 2) continue;

      const string& seq_id = c[SEQ_ID];
      int32_t cn_start = from_string<int32_t>(c[START]);
      int32_t cn_end   = from_string<int32_t>(c[END]);

      if (ref_seq_info.seq_id_to_index(seq_id) >= ref_seq_info.size()) continue;
      cAnnotatedSequence& seq = ref_seq_info[seq_id];

      // Try each end as the NEW copy in turn; the other end then has to carry the annotated one.
      for (int end_choice = 0; end_choice < 2; end_choice++) {
        int32_t annotated_at = (end_choice == 0) ? cn_end : cn_start;
        int32_t anchor_at    = (end_choice == 0) ? cn_start : cn_end;

        // The new copy is not in the reference, so only a MOB or a junction reports it -- and it is
        // that anchor, not the nearest annotation, which says WHICH family recombined.
        bool mob_anchor = false;
        int32_t anchor_position = 0;
        // Which way round the NEW copy sits, in reference coordinates. An amplification arises by
        // recombination between two copies of the element, and only DIRECT repeats -- copies in the
        // same orientation -- give a tandem duplication; inverted copies give an inversion instead.
        // So this is what the copy at the other end has to match.
        int32_t anchor_strand = 0;
        string anchor_id, family;
        cDiffEntry* anchor_jc = NULL;
        diff_entry_ptr_t anchor_mob_ptr;

        for (diff_entry_list_t::iterator mi = mob_list.begin(); mi != mob_list.end(); mi++) {
          cDiffEntry& m = **mi;
          if (m[SEQ_ID] != seq_id) continue;
          int32_t p = from_string<int32_t>(m[POSITION]);
          if (abs(p - anchor_at) > kAnchorSlop) continue;
          mob_anchor = true; anchor_position = p; anchor_id = m._id; family = m[REPEAT_NAME];
          anchor_strand = from_string<int32_t>(m[STRAND]);
          anchor_mob_ptr = *mi;
          break;
        }

        if (!mob_anchor) {
          for (diff_entry_list_t::iterator ji = jc.begin(); ji != jc.end(); ji++) {
            cDiffEntry& j = **ji;
            for (int side = 1; side <= 2; side++) {
              string me    = "side_" + to_string(side) + "_";
              string other = "side_" + to_string(3 - side) + "_";
              if (j[me + "seq_id"] != seq_id) continue;
              int32_t p = n(j[me + "position"]);
              if (abs(p - anchor_at) > kAnchorSlop) continue;
              // The partner side has to sit in a repeat family -- that is the evidence a new copy of
              // one landed here. Read it from the _side_N_is_* fields prepare_junctions computes, NOT
              // from side_N_gene_name: gene names are filled in by annotate_mutations, which runs
              // after prediction, so at this point they are empty.
              if (!j.entry_exists("_" + other + "is")) continue;
              if (j["_" + other + "is_name"].empty()) continue;
              anchor_position = p; anchor_id = j._id; family = j["_" + other + "is_name"];
              // `me` is the unique side (it sits at the CN boundary), `other` the repeat side.
              anchor_strand = new_repeat_copy_strand(j, me, other);
              anchor_jc = &j;
              break;
            }
            if (!anchor_id.empty()) break;
          }
        }
        if (anchor_id.empty() || family.empty()) continue;

        // The other end. Usually a pre-existing copy of the same family, annotated in the reference,
        // and required to be in the SAME orientation as the new copy -- see anchor_strand above.
        cFeatureLocation* rep = closest_repeat_of_family(seq, family, annotated_at,
                                                         (end_choice == 0) ? -1 : 1,
                                                         kRepeatBoundarySlop, anchor_strand);

        // ...but both ends can be NEW copies, when a transposition dropped one element and the
        // amplification then recombined it against a second new copy. Nothing in the reference marks
        // either end, so the far end has to be anchored the same way this one was.
        //
        // NOTE: the MOB minted for that far end below is used only to establish that the junction
        // really does place an element there, and to read its orientation. It is never added to the
        // genome diff -- the insertion itself is reported (or not) by the ordinary MOB paths.
        bool far_mob_is_new = false;
        int32_t far_position = 0;
        string far_id;
        diff_entry_ptr_t far_mob_ptr;
        cDiffEntry* far_jc = NULL;

        if (rep == NULL) {
          for (diff_entry_list_t::iterator mi = mob_list.begin(); mi != mob_list.end(); mi++) {
            cDiffEntry& m = **mi;
            if (m[SEQ_ID] != seq_id) continue;
            if (m[REPEAT_NAME] != family) continue;
            // Direct repeats only, as at the annotated end.
            if (anchor_strand && (from_string<int32_t>(m[STRAND]) != anchor_strand)) continue;
            int32_t p = from_string<int32_t>(m[POSITION]);
            if (abs(p - annotated_at) > kAnchorSlop) continue;
            far_position = p; far_id = m._id; far_mob_ptr = *mi;
            break;
          }

          // Only a junction there: the element is real but its second junction lies inside a copy of
          // itself and was discarded, so mint the MOB it implies with the duplication size left
          // undetermined rather than leaving the insertion unreported.
          if (far_id.empty()) {
            for (diff_entry_list_t::iterator ji = jc.begin(); ji != jc.end(); ji++) {
              cDiffEntry& j2 = **ji;
              if (&j2 == anchor_jc) continue;
              for (int side = 1; side <= 2; side++) {
                string me    = "side_" + to_string(side) + "_";
                string other = "side_" + to_string(3 - side) + "_";
                if (j2[me + "seq_id"] != seq_id) continue;
                int32_t p = n(j2[me + "position"]);
                if (abs(p - annotated_at) > kAnchorSlop) continue;
                if (!j2.entry_exists("_" + other + "is")) continue;
                if (j2["_" + other + "is_name"] != family) continue;
                // A fresh entry each time: the minting appends to _evidence, so a candidate rejected
                // below would leave its junction id behind for the next one to inherit.
                cDiffEntry candidate;
                if (!mob_from_single_IS_junction(j2, settings.polymorphism_prediction, candidate)) continue;
                // Direct repeats only. The minted strand comes from the same triple product
                // anchor_strand does, so the two are on one scale and comparable.
                if (anchor_strand && (from_string<int32_t>(candidate[STRAND]) != anchor_strand)) continue;
                far_position = p; far_mob_is_new = true; far_jc = &j2;
                break;
              }
              if (far_mob_is_new) break;
            }
          }
          if (far_id.empty() && !far_mob_is_new) continue;
        }

        // Two shapes, and which one applies turns on whether the reference already contains a copy of
        // the element inside the amplified interval. Both were checked against the curated LTEE
        // genome diffs (barricklab/LTEE-Ecoli), which agree with the coordinates derived here to the
        // base in both cases.
        //
        //  (a) A pre-existing ANNOTATED copy at one end. Run the block from the new copy through the
        //      end of that element, so the element is inside the block and repeating it reproduces
        //      one between every pair of copies. Nothing else is needed -- no mediated=, no extra MOB.
        //      Curated form: plain "AMP REL606 2749275 26603 2".
        //
        //  (b) BOTH copies new. The reference has no element in the interval, so the block is pure
        //      unique sequence and the element has to come from mediated=, which lays one down with
        //      each new copy. The MOB that reports the first copy anchors the low end, and the AMP is
        //      placed within= its target-site duplication so the two are ordered: the insertion
        //      happened first, then recombination amplified the region between it and the second copy.
        //      Curated form: "MOB ... IS186 1 8" plus
        //      "AMP REL606 3517306 108140 2 mediated=IS186 mediated_strand=-1 within=<mob>:2".
        //      within= only means something when the MOB's duplication size is known; a MOB minted
        //      from a lone junction has an indeterminate one and therefore no second copy to sit in.
        //
        // Only evidence items go in _evidence -- a MOB is a mutation, and a non-evidence id is
        // silently dropped from the list, so the tie to it is carried by within=.
        int32_t position, last;
        string mediated_family;
        int32_t mediated_strand = 0;
        diff_entry_ptr_t order_after_mob;
        vector<string> evidence = make_vector<string>(c._id);
        if (anchor_jc != NULL) evidence.push_back(anchor_jc->_id);

        if (rep != NULL) {
          if (end_choice == 0) {        // annotated copy at the high-coordinate end
            position = anchor_position;
            last     = rep->get_end_1();
          } else {                      // annotated copy at the low-coordinate end
            position = rep->get_start_1();
            last     = anchor_position;
          }
        } else {
          // Needs a MOB to anchor and to order against; a pair of bare junctions says an element is
          // involved but not that it was ever placed.
          if (!mob_anchor) continue;
          // Written with the anchoring MOB at the low end, matching the curated layout.
          if (far_position <= anchor_position) continue;

          position = anchor_position;
          last     = far_position;
          mediated_family = family;
          mediated_strand = from_string<int32_t>((*anchor_mob_ptr)[STRAND]);
          if (far_jc != NULL) evidence.push_back(far_jc->_id);

          // Order the amplification after the insertion it recombined against, when there are two
          // copies of the target site to speak of.
          if (!anchor_mob_ptr->entry_exists("indeterminate_duplication_size")
              && (from_string<int32_t>((*anchor_mob_ptr)["duplication_size"]) > 0))
            order_after_mob = anchor_mob_ptr;
        }

        int32_t size = last - position + 1;
        if (size <= kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation) continue;

        cDiffEntry mut(AMP);
        mut._evidence = evidence;
        mut[SEQ_ID] = seq_id;
        mut[POSITION] = to_string<int32_t>(position);
        mut[SIZE] = to_string<int32_t>(size);
        mut[NEW_COPY_NUMBER] = to_string<int32_t>(copies);

        if (!mediated_family.empty()) {
          mut[MEDIATED] = mediated_family;
          mut[MEDIATED_STRAND] = to_string<int32_t>(mediated_strand);
          if (order_after_mob.get() != NULL)
            mut["within"] = order_after_mob->_id + ":2";
        }

        gd.add(mut);
        break; // this CN region is accounted for
      }
    }
  }

  // Join MC fragments that one CN copy-number-0 region spans. See the header for why.
  //
  // Runs before mutation prediction (from breseq_cmdline, on the merged evidence genome diff) rather
  // than inside predict(), so the joined MC is what gets written to evidence.gd and shown in the
  // HTML. Doing it later would leave the on-disk evidence and the report disagreeing with the
  // deletion that was predicted from it.
  void merge_MC_fragments_spanned_by_CN(cGenomeDiff& gd, uint32_t max_island)
  {
    diff_entry_list_t cn_list = gd.get_list(make_vector<gd_entry_type>(CN));
    if (cn_list.empty()) return;
    diff_entry_list_t mc_list = gd.get_list(make_vector<gd_entry_type>(MC));
    if (mc_list.size() < 2) return;

    set<string> absorbed;

    for (diff_entry_list_t::iterator ci = cn_list.begin(); ci != cn_list.end(); ci++) {
      cDiffEntry& cn = **ci;

      // copy_number is a POSITIONAL field on CN, so it is always present; only a call of zero says
      // the region is deleted outright, which is the only thing that licenses bridging an island.
      if (from_string<int32_t>(cn[COPY_NUMBER]) != 0) continue;

      const string& seq_id = cn[SEQ_ID];
      int32_t cn_start = from_string<int32_t>(cn[START]);
      int32_t cn_end   = from_string<int32_t>(cn[END]);

      // Fragments that are mostly inside this CN region. Overlap rather than containment: a fragment
      // routinely runs past the CN boundary, because the HMM only resolves to its tile size while the
      // MC is called per base (in the real manB/cpsG case the first fragment starts ~1 kb before the
      // CN does). Requiring the majority of the fragment to lie inside keeps a large MC that merely
      // clips the edge of the region from being dragged in.
      vector<cDiffEntry*> frags;
      for (diff_entry_list_t::iterator mi = mc_list.begin(); mi != mc_list.end(); mi++) {
        cDiffEntry& mc = **mi;
        if (mc[SEQ_ID] != seq_id) continue;
        if (absorbed.count(mc._id)) continue;

        int32_t mc_start = from_string<int32_t>(mc[START]);
        int32_t mc_end   = from_string<int32_t>(mc[END]);
        int32_t overlap  = min(mc_end, cn_end) - max(mc_start, cn_start) + 1;
        if (overlap <= 0) continue;
        if (2 * overlap < (mc_end - mc_start + 1)) continue;

        frags.push_back(&mc);
      }
      if (frags.size() < 2) continue;

      sort(frags.begin(), frags.end(), mc_fragment_start_lt);

      // Merge maximal runs of fragments separated by small islands, so one oversized gap costs only
      // the fragments across it rather than the whole region.
      size_t run_start = 0;
      for (size_t i = 1; i <= frags.size(); i++) {
        bool breaks_run = true;
        if (i < frags.size()) {
          int32_t island = from_string<int32_t>((*frags[i])[START])
                         - from_string<int32_t>((*frags[i - 1])[END]) - 1;
          breaks_run = (island > static_cast<int32_t>(max_island));
        }
        if (!breaks_run) continue;

        if (i - run_start >= 2) {
          cDiffEntry& first = *frags[run_start];
          cDiffEntry& last  = *frags[i - 1];

          // The joined entry keeps the outer extents and each end's own coverage numbers; the inner
          // boundaries the islands created are exactly what we are discarding.
          first[END] = last[END];
          if (last.entry_exists(END_RANGE))         first[END_RANGE]         = last[END_RANGE];
          if (last.entry_exists(RIGHT_INSIDE_COV))  first[RIGHT_INSIDE_COV]  = last[RIGHT_INSIDE_COV];
          if (last.entry_exists(RIGHT_OUTSIDE_COV)) first[RIGHT_OUTSIDE_COV] = last[RIGHT_OUTSIDE_COV];

          for (size_t k = run_start + 1; k < i; k++) absorbed.insert(frags[k]->_id);
        }
        run_start = i;
      }
    }

    if (absorbed.empty()) return;

    diff_entry_list_t* all = gd.get_mutable_list_ptr();
    for (diff_entry_list_t::iterator it = all->begin(); it != all->end(); ) {
      if (((*it)->_type == MC) && (absorbed.count((*it)->_id) > 0))
        it = gd.remove(it);
      else
        it++;
    }
  }

  void MutationPredictor::remove_DP_superseded_by_PD(Settings& settings, Summary& summary, cGenomeDiff& gd)
  {
    (void)settings;

    diff_entry_list_t pd_list = gd.get_list(make_vector<gd_entry_type>(PD));
    if (pd_list.empty()) return;
    diff_entry_list_t dp_list = gd.get_list(make_vector<gd_entry_type>(DP));
    if (dp_list.empty()) return;

    bool inner3p = true;
    double D = 0.0, pair_median = 0.0;
    if (!paired_library_params(summary, inner3p, D, pair_median)) return;
    int32_t window = max(1, static_cast<int32_t>(pair_median));

    // Each doomed DP remembers which PD displaced it, so a mutation that was resting on the DP can
    // be handed the PD instead. predictMCplusDPtoDEL builds a deletion out of an MC plus a DP, and
    // simply dropping that DP would leave the deletion standing on half the evidence it was
    // predicted from, with nothing in the output explaining where the other half went.
    map<string, string> superseded_by;
    for (diff_entry_list_t::iterator di = dp_list.begin(); di != dp_list.end(); di++) {
      cDiffEntry& dp = **di;
      for (diff_entry_list_t::iterator pi = pd_list.begin(); pi != pd_list.end(); pi++) {
        cDiffEntry& pd = **pi;
        if (pd.entry_exists(REJECT)) continue;
        if (!same_junction_side(dp, "side_1", pd, "side_1", window)) continue;
        if (!same_junction_side(dp, "side_2", pd, "side_2", window)) continue;
        superseded_by[dp._id] = pd._id;
        break;
      }
    }
    if (superseded_by.empty()) return;

    // Substitute in every mutation's evidence list first, so nothing is left pointing at an entry
    // that no longer exists.
    diff_entry_list_t muts = gd.mutation_list();
    for (diff_entry_list_t::iterator it = muts.begin(); it != muts.end(); it++) {
      cDiffEntry& mut = **it;
      vector<string> updated;
      for (size_t k = 0; k < mut._evidence.size(); k++) {
        map<string, string>::const_iterator s = superseded_by.find(mut._evidence[k]);
        const string& id = (s == superseded_by.end()) ? mut._evidence[k] : s->second;
        if (find(updated.begin(), updated.end(), id) == updated.end()) updated.push_back(id);
      }
      mut._evidence = updated;
    }

    diff_entry_list_t* all = gd.get_mutable_list_ptr();
    for (diff_entry_list_t::iterator it = all->begin(); it != all->end(); ) {
      if (((*it)->_type == DP) && (superseded_by.count((*it)->_id) > 0))
        it = gd.remove(it);
      else
        it++;
    }
  }

  // Repeat-copy region "seq:start-end" that (position, strand) sits in, if it belongs to repeat family
  // repeat_name; else "". Used to map a JC/DP IS side to its specific copy.
  static string mob_copy_region_for(cReferenceSequences& ref_seq_info, const string& seq_id, int32_t pos,
                                     int32_t strand, const string& repeat_name)
  {
    int32_t md = 50;
    cFeatureLocation* is = cReferenceSequences::find_closest_repeat_region_boundary(pos, ref_seq_info[seq_id].m_repeats, md, strand, true);
    if (is == NULL || is->get_feature()->SafeGet("name") != repeat_name) return "";
    return seq_id + ":" + to_string(is->get_start_1()) + "-" + to_string(is->get_end_1());
  }

  void MutationPredictor::combine_DP_with_MOB_by_unique_side(cGenomeDiff& gd)
  {
    // Collect non-rejected DP items, recording their unique side and IS side. The IS side is the one
    // whose annotate_key is "repeat" (NOT the redundant flag -- a side that uniquely matched a variant
    // copy is a repeat side but is flagged redundant=0).
    struct DPitem { string id; string useq; int32_t upos; int32_t ustr; string iseq; int32_t ipos; int32_t istr; int32_t count; };
    vector<DPitem> dps;
    diff_entry_list_t dp_list = gd.get_list(make_vector<gd_entry_type>(DP));
    for (diff_entry_list_t::iterator it = dp_list.begin(); it != dp_list.end(); it++) {
      cDiffEntry& dp = **it;
      if (dp.entry_exists(REJECT)) continue;
      string ak1 = dp.entry_exists("side_1_annotate_key") ? dp["side_1_annotate_key"] : "";
      string ak2 = dp.entry_exists("side_2_annotate_key") ? dp["side_2_annotate_key"] : "";
      string is, un;
      if      (ak1 == "repeat" && ak2 != "repeat") { is = "side_1"; un = "side_2"; }
      else if (ak2 == "repeat" && ak1 != "repeat") { is = "side_2"; un = "side_1"; }
      else continue;  // need exactly one IS (repeat) side and one unique side
      DPitem d; d.id = dp._id;
      d.iseq = dp[is + "_seq_id"]; d.ipos = n(dp[is + "_position"]); d.istr = n(dp[is + "_strand"]);
      d.useq = dp[un + "_seq_id"]; d.upos = n(dp[un + "_position"]); d.ustr = n(dp[un + "_strand"]);
      d.count = dp.entry_exists("discordant_count") ? n(dp["discordant_count"]) : 0;
      dps.push_back(d);
    }
    if (dps.empty()) return;

    // A DP whose IS side lands on a different copy than the MOB's IS side won't match by both sides, so
    // match by the UNIQUE side only, within a small window and the same strand.
    const int32_t SLOP = 20;
    diff_entry_list_t mobs = gd.get_list(make_vector<gd_entry_type>(MOB));
    for (diff_entry_list_t::iterator it = mobs.begin(); it != mobs.end(); it++) {
      cDiffEntry& mut = **it;
      string repeat_name = mut["repeat_name"];
      string consensus = ref_seq_info.repeat_family_sequence(repeat_name, 1, NULL, NULL, NULL, false);

      // Best NON-consensus (variant) copy that any supporting evidence localizes, driving mob_region:
      //  priority 3 = a supporting JC's IS side (uniquely matched, redundant=0 -> not canonicalized)
      //  priority 2 = a matched DP's IS side (tie-broken by discordant_count)
      // A JC IS side with redundant=1 is canonicalized to the consensus copy, so it never localizes a
      // variant. Higher priority wins; within a priority, higher tie wins.
      int best_prio = 0; int32_t best_tie = -1; string best_region;
      // Returns true and records the region when it is a non-consensus variant of this family.
      // (defined as a lambda so both the JC and DP loops can call it.)
      auto consider_variant = [&](const string& region, int prio, int32_t tie) {
        if (region.empty()) return;
        string reg = region;  // repeat_family_sequence takes a non-const string*
        string specific = ref_seq_info.repeat_family_sequence(repeat_name, 1, &reg, NULL, NULL, false);
        if (specific.empty() || specific == consensus) return;   // consensus copy -> no override
        if (prio > best_prio || (prio == best_prio && tie > best_tie)) {
          best_prio = prio; best_tie = tie; best_region = region;
        }
      };

      // The MOB's unique insertion loci come from its supporting JC(s)' unique side (annotate_key !=
      // "repeat"). Gather those AND any variant copy a JC IS side uniquely localizes.
      vector<pair<string, pair<int32_t, int32_t> > > jc_unique;  // (seq_id, (position, strand))
      for (vector<string>::iterator ev = mut._evidence.begin(); ev != mut._evidence.end(); ev++) {
        diff_entry_ptr_t e = gd.find_by_id(*ev);
        if (e.get() == NULL || e->_type != JC) continue;
        cDiffEntry& j = *e;
        string ak1 = j.entry_exists("side_1_annotate_key") ? j["side_1_annotate_key"] : "";
        string ak2 = j.entry_exists("side_2_annotate_key") ? j["side_2_annotate_key"] : "";
        string is, un;
        if      (ak1 == "repeat" && ak2 != "repeat") { is = "side_1"; un = "side_2"; }
        else if (ak2 == "repeat" && ak1 != "repeat") { is = "side_2"; un = "side_1"; }
        else continue;
        jc_unique.push_back(make_pair(j[un + "_seq_id"], make_pair(n(j[un + "_position"]), n(j[un + "_strand"]))));
        // A uniquely-matched (redundant=0) JC IS side is the strongest variant localization.
        bool is_redundant = j.entry_exists(is + "_redundant") && n(j[is + "_redundant"]);
        if (!is_redundant)
          consider_variant(mob_copy_region_for(ref_seq_info, j[is + "_seq_id"], n(j[is + "_position"]), n(j[is + "_strand"]), repeat_name), 3, 0);
      }
      if (jc_unique.empty()) continue;

      // Attach every DP whose unique side matches a JC unique side; a matched DP that localizes a
      // variant contributes a (weaker) mob_region candidate.
      for (vector<DPitem>::iterator d = dps.begin(); d != dps.end(); d++) {
        bool match = false;
        for (size_t k = 0; k < jc_unique.size(); k++) {
          if (d->useq == jc_unique[k].first
              && abs(d->upos - jc_unique[k].second.first) <= SLOP
              && d->ustr == jc_unique[k].second.second) { match = true; break; }
        }
        if (!match) continue;

        if (find(mut._evidence.begin(), mut._evidence.end(), d->id) == mut._evidence.end())
          mut._evidence.push_back(d->id);

        consider_variant(mob_copy_region_for(ref_seq_info, d->iseq, d->ipos, d->istr, repeat_name), 2, d->count);
      }

      if (!best_region.empty()) mut["mob_region"] = best_region;
    }
  }

  // Returns the size of the minimum repeated subseqence
  void MutationPredictor::find_repeat_unit(string& mutation_sequence, uint32_t& repeat_size, string& repeat_sequence)
  {
    // Find the shortest sub-repeat.
    repeat_size = mutation_sequence.size();
    
    // there's no need to test any size where two copies of the sub-repeat are longer than the sequence
    for (uint32_t i= static_cast<uint32_t>(trunc(mutation_sequence.size() / 2)); i>0; i--) {
      
      // Must evenly divide the sequence
      if (repeat_size % i != 0) continue;
      
      string unit = mutation_sequence.substr(0, i);
      bool is_repeat = true;
      uint32_t test_pos = i; // change to zero offset right past first repeat unit
      
      while(test_pos + i <= mutation_sequence.size()) {
        // compare new unit to original unit
        string test = mutation_sequence.substr(test_pos, i);
        if (unit != test) {
          is_repeat = false;
          break;
        }
        test_pos += i;
      }
      
      if (is_repeat) repeat_size = i;
    }
    
    repeat_sequence = mutation_sequence.substr(0, repeat_size);
  }
  
  void MutationPredictor::normalizeINSposition(cAnnotatedSequence& ref_seq, cDiffEntry& de, string& repeat_sequence)
  {

    int32_t test_position = 1 + from_string<int32_t>(de[POSITION]);   
    // The offset of 1 ensures we are on the first base of a possible repeat unit (relative to the reference).

    // attempt to move forward by repeat units at a time from the current position
    while ( test_position + repeat_sequence.size() - 1 <= ref_seq.get_sequence_length()) {
      
      string test_sequence = ref_seq.get_sequence_1(test_position, test_position + repeat_sequence.size() - 1);
      
      if (test_sequence != repeat_sequence) break;
      test_position += repeat_sequence.size();
    }
        
    string new_position = to_string<int32_t>(test_position - 1);
    
    if (new_position != de[POSITION]) {
      de["_original_aligned_position"] = de[POSITION]; // Save the original position for marking in alignments
      de[POSITION] = new_position;
    }
    
    // We still may need to move a fraction of the repeat to deal with 
    // equivalent cases being initially described in different ways:
    //
    // ATCG ATCG ATCG ATCG +(ATCG) AT  - INITIAL
    // ATCG ATCG ATCG ATCG  AT +(CGAT) - PREFERRED 
    
    // Recall that position is where it is inserted after
    int32_t size = repeat_sequence.size();
    if (size > 1) {
      int32_t test_size = size - 1;
      int32_t test_pos_ref = from_string<int32_t>(de[POSITION]) + 1;
      int32_t test_pos_ins = 0; // 0-indexed

      while (test_size > 0) {
        
        string test_seq_1 = ref_seq.get_sequence_1(test_pos_ref, test_pos_ref+test_size-1);
        string test_seq_2 = repeat_sequence.substr(0, test_size);
        
        if (test_seq_1 == test_seq_2) {
          // If not already shifted
          if (!de.entry_exists("_original_aligned_position"))
            de["_original_aligned_position"] = de[POSITION];           
          de[POSITION] = to_string<int32_t>(test_pos_ref+test_size-1);
          string new_seq = de[NEW_SEQ];
          new_seq=  new_seq.substr(test_size) + new_seq.substr(0, test_size);
          de[NEW_SEQ] = new_seq;
          
          break;
        }
        
        test_size--;
      }
    }

  }
  
  // repeat_unit_sequence is the smallest repeated unit
  //
  
  void MutationPredictor::normalizeDELposition(cAnnotatedSequence& ref_seq, cDiffEntry& de, string& repeat_unit_sequence)
  {
    bool verbose = false;
    
    if (verbose)
      cout << de << endl;
    
    // Don't shift mediated or between mutations
    if (de.entry_exists("mediated") || de.entry_exists("between")) return;
    
    uint32_t mutation_size = from_string<uint32_t>(de[SIZE]);
    int32_t original_position = from_string<int32_t>(de[POSITION]);
    int32_t test_position(original_position);
    int32_t new_position(original_position);
    
    // Offset the test position to where the next repeat would fall.
    // Ex: initial del    GG
    //     genome       CTGGGTAAGCTAG
    //     mutation pos   ^
    //     test pos         *
    //     final del       GG
    
    test_position += mutation_size;
    
    // The offset of 1 ensures we are on the first base of a possible repeat unit (relative to the reference).
    
    // attempt to move forward by repeat units at a time from the current position
    // check to be sure the entire new repeat that we might shift mutation to is in bounds!!
    while ( test_position + repeat_unit_sequence.size() - 1 <= ref_seq.get_sequence_length()) {
      
      string test_sequence = ref_seq.get_sequence_1(test_position, test_position + repeat_unit_sequence.size() - 1);
      if (test_sequence == repeat_unit_sequence) {
        // New positon: remember, needs to be back at the start of the whole mutation
        test_position += repeat_unit_sequence.size();
        new_position = test_position - mutation_size;
      } else {
        break;
      }
    }
    
    if (new_position != original_position) {
      de["_original_aligned_position"] = de[POSITION]; // Save the original position for marking in alignments
      de[POSITION] = to_string<int32_t>(new_position);
    }
    
    // This section checks for rare cases where a DEL can be written in two ways 
    // because there is a repeat separated by unique bases that are also deleted
    //                     * test_position
    // REF : ATAA TCGCCAGC G TCGCCAGC ACTG
    // DEL1:      DDDDDDDD D
    // DEL2:               D DDDDDDDD
    //
    // REF : ATAA TCGCCAGC AGC TCGCCAGC ACTG
    // DEL1:      DDDDDDDD DDD 
    // DEL2:               DDD DDDDDDDD   
    //
    // REF : ATAA C CGCCAGCAGCTCGCCAGC C ACTG
    // DEL1:      D DDDDDDDDDDDDDDDDDD
    // DEL2:        DDDDDDDDDDDDDDDDDD D
    
    int32_t size = from_string<int32_t>(de[SIZE]);
    if (size > 1) {
      int32_t test_size = size;
      int32_t test_pos_1 = from_string<int32_t>(de[POSITION]);
      int32_t test_pos_2 = test_pos_1 + test_size;      
      test_size--;
      
      while (test_size > 0) {
        
        if (test_pos_2 + test_size - 1 > static_cast<int32_t>(ref_seq.get_sequence_length()))
          break;
        
        string test_seq_1 = ref_seq.get_sequence_1(test_pos_1, test_pos_1+test_size-1);
        string test_seq_2 = ref_seq.get_sequence_1(test_pos_2, test_pos_2+test_size-1);

        if (test_seq_1 == test_seq_2) {
          // If not already shifted
          if (!de.entry_exists("_original_aligned_position"))
              de["_original_aligned_position"] = de[POSITION];           
          de[POSITION] = to_string<int32_t>(test_pos_1+test_size);
          break;
        }
        test_size--;
      }
    }
  }
  
  
  ///
  // Expects position to be the first position of the repeat we are looking for
  // due to earlier normalization (shifting) of these mutations.
  //
  // For DEL it is the first deleted base (mutation position).
  //
  //          ***
  // ABCABCABCABC    => 4 repeats in reference
  //          ^
  //
  // For INS this is the position one past where the insertion starts (mutation position+1).
  //
  //          ABC
  // ABCABCABC   XYZ => 3 repeats in reference
  //             ^
  //
  ///
  
  uint32_t MutationPredictor::find_original_num_repeat_units(cAnnotatedSequence& ref_seq, int32_t position, string& repeat_sequence)
  {
    uint32_t num_repeat_units = 0;
    
    int32_t test_position = position;
    
    // count repeats backwards
    while ( test_position - repeat_sequence.size() > 1) {
      
      test_position -= repeat_sequence.size();
      
      // Safety valve for overflowing sequence
      if (test_position <= 1) break;
      
      string test_sequence = ref_seq.get_sequence_1_start_size(test_position, repeat_sequence.size());
      
      if (test_sequence != repeat_sequence) break;
      
      num_repeat_units++;
    }
    
    return num_repeat_units;
  }

  
  // We should only have to combine INS-SNP, DEL-SNP, and INS-DEL, DEL-INS combinations
  void MutationPredictor::combine_newly_adjacent_mutations(cGenomeDiff& gd)
  {
    diff_entry_list_t muts = gd.mutation_list();
    diff_entry_list_t new_muts;
    
    bool any_changes(false);
    cDiffEntry* last_mut(NULL);
    for (diff_entry_list_t::iterator it=muts.begin(); it != muts.end(); it++) {
      
      cDiffEntry& mut = **it;
      
      // Require both mutations to be one of these four types
      if ((mut._type != SNP) && (mut._type != DEL) && (mut._type != INS) && (mut._type != SUB)) {
        last_mut = NULL;
        continue;
      }
      
      // If they are directly adjacent on the same SEQ_ID (add one to the position)
      // Notice that we require muts to be of different types, this is to
      // prevent cases of two INS mutations that should have been combined previous
      // and would spuriously be combined based on our checking of positions here.
      if (last_mut && ((*last_mut)[SEQ_ID] == mut[SEQ_ID]) && (last_mut->_type != mut._type)) {
        if (mut.get_reference_coordinate_start() - last_mut->get_reference_coordinate_end() <= cReferenceCoordinate(1, 0) ) {
          
          //cout << "Combining:\n" + last_mut->as_string() + "\n" + mut.as_string() + "\n";
          
          any_changes = true;
          
          // We can only make a SUB mutation
          cDiffEntry* new_mut = new cDiffEntry(SUB);
          (*new_mut)[SEQ_ID] = mut[SEQ_ID];
          
          new_mut->_evidence.insert(new_mut->_evidence.end(), mut._evidence.begin(), mut._evidence.end());
          new_mut->_evidence.insert(new_mut->_evidence.end(), last_mut->_evidence.begin(), last_mut->_evidence.end());
          (*new_mut)[POSITION] = (last_mut->_type != INS) ? (*last_mut)[POSITION] : mut[POSITION];
          
          int32_t size(0);
          string new_seq;

          if (last_mut->_type == SNP) {
            size++;
            new_seq += (*last_mut)[NEW_SEQ];
          } else if (last_mut->_type == DEL) {
            size += from_string<int32_t>((*last_mut)[SIZE]);
          } else if (last_mut->_type == INS) {
            new_seq += (*last_mut)[NEW_SEQ];
          }
          
          if (mut._type == SNP) {
            size++;
            new_seq += mut[NEW_SEQ];
          } else if (mut._type == DEL) {
            size += from_string<int32_t>(mut[SIZE]);
          } else if (mut._type == INS) {
            new_seq += mut[NEW_SEQ];
          }
          
          (*new_mut)[SIZE] = to_string(size);
          (*new_mut)[NEW_SEQ] = new_seq;
          
          //cout << "New:\n" + new_mut->as_string() + "\n";
          
          // Erase both old entries
          it--;           // start on last_mut
          it = muts.erase(it); // advances to mut
          it = muts.erase(it); // now points past mut
          
          // Insert new mutation before the specified position
          diff_entry_ptr_t new_mut_p(new_mut);
          muts.insert(it, new_mut_p);
          it--;           // now go back to point at the inserted mutation
        }
      }
      
      last_mut = &**it;
    }
    
    
    
    // Sorts and fixes ID gaps
    if (any_changes) gd.reassign_unique_ids();
  }
  
  
  
  // Currently, only adds before tag to any small mutation predicted that is within the
  // duplication region of a MOB in consensus mode as 'before'
  
  void MutationPredictor::assign_before_within_to_mutations(cGenomeDiff& gd)
  {
    diff_entry_list_t muts = gd.mutation_list();
    
    // Since they are sorted, we only have to look around the MOB in the list, for efficiency.
    // IMPORTANT: MOBs with the same position get sorted after SNP, INS, DEL, etc., so we must check backwards.
    
    for (diff_entry_list_t::iterator it_mob=muts.begin(); it_mob != muts.end(); it_mob++) {
      
      cDiffEntry& mob = **it_mob;
      if (mob._type != MOB) continue;
      
      // Find the beginning
      diff_entry_list_t::iterator it_test_start = muts.begin();
      if (it_mob != muts.begin()) {
        it_test_start=it_mob;
        it_test_start--;
        while (it_test_start!=muts.begin()) {
          cDiffEntry& test = **it_test_start;
          if (test[SEQ_ID] != mob[SEQ_ID]) break;
          if (test.get_reference_coordinate_start() < mob.get_reference_coordinate_start()) break;
          it_test_start--;
        }
      }
      
      // Find the end
      diff_entry_list_t::iterator it_test_end = muts.end();
      if (it_mob != muts.end()) {
        it_test_end=it_mob;
        it_test_end++;
        while (it_test_end!=muts.end()) {
          cDiffEntry& test = **it_test_end;
          if (test[SEQ_ID] != mob[SEQ_ID]) break;
          if (test.get_reference_coordinate_start() > mob.get_reference_coordinate_end()) break;
          it_test_end++;
        }
      }
      
      if (it_test_start != it_test_end) {
        for (diff_entry_list_t::iterator it_test=it_test_start; it_test != it_test_end; it_test++) {
          if (it_test == it_mob) continue;
          cDiffEntry& test = **it_test;
          
          // Must be entirely within
          if ( test.located_within(mob) ) {
            test[BEFORE] = to_string(mob._id);
          }
        }
      }
      
    }
  }
  
  
  string BaseSubstitutionEffects::separator = ".";
  
  vector<base_char> BaseSubstitutionEffects::base_char_list = make_vector<base_char>
  ('A')('T')('C')('G')
  ;
  
  map<base_char, uint8_t> BaseSubstitutionEffects::base_char_to_base_index = make_map<base_char, uint8_t>
  ('A',0)('T',1)('C',2)('G',3)
  ;
  
  vector<string>  BaseSubstitutionEffects::base_change_list = make_vector<string>
  ("A.G")("T.C")("A.C")("T.G")("A.T")("T.A")("G.A")("C.T")("G.T")("C.A")("G.C")("C.G")
  ;
  
  map<string,uint8_t>  BaseSubstitutionEffects::base_change_to_base_pair_change_index = make_map<string,uint8_t>
  ("A.G",0)("A.C",1)("A.T",2)("G.A",3)("G.T",4)("G.C",5)
  ("T.C",0)("T.G",1)("T.A",2)("C.T",3)("C.A",4)("C.G",5)
  ;
  
  vector<string>  BaseSubstitutionEffects::base_pair_change_list = make_vector<string>
  ("AT.GC")("AT.CG")("AT.TA")("CG.TA")("CG.AT")("CG.GC")
  ;
  
  map<string,string>  BaseSubstitutionEffects::base_change_to_base_pair_change = make_map<string,string>
  ("A.G","AT.GC")("A.C","AT.CG")("A.T","AT.TA")
  ("C.T","CG.TA")("C.A","CG.AT")("C.G","CG.GC")
  ("T.C","AT.GC")("T.G","AT.CG")("T.A","AT.TA")
  ("G.A","CG.TA")("G.T","CG.AT")("G.C","CG.GC")
  ;
  
  vector<string> BaseSubstitutionEffects::base_change_type_list = make_vector<string> 
  ("INTERGENIC")("NONCODING")("PSEUDOGENE")("SYNONYMOUS")("NONSYNONYMOUS")("NONSENSE")("NO_CHANGE")("TOTAL")
  ;
 
   vector<string>  BaseSubstitutionEffects::base_type_list = make_vector<string>
  ("INTERGENIC")("NONCODING")("PSEUDOGENE")("PROTEIN")("TOTAL")
  ;  
  
  vector<string> BaseSubstitutionEffectCounts::base_pair_change_count_list = make_vector<string> 
  ("AT.GC")("AT.CG")("AT.TA")("CG.TA")("CG.AT")("CG.GC")("TOTAL")
  ;
  
  vector<string> BaseSubstitutionEffectCounts::base_change_type_count_list = make_vector<string> 
  ("INTERGENIC")("NONCODING")("PSEUDOGENE")("SYNONYMOUS")("NONSYNONYMOUS")("NONSENSE")("TOTAL")
  ;
  
  map<BaseSubstitutionEffect,BaseType>  BaseSubstitutionEffects::snp_type_to_base_type = make_map<BaseSubstitutionEffect,BaseType>
  (intergenic_base_substitution,intergenic_base)
  (pseudogene_base_substitution,pseudogene_base)
  (noncoding_base_substitution,noncoding_base)
  (synonymous_coding_base_substitution,protein_base)
  (nonsynonymous_coding_base_substitution,protein_base)
  (nonsense_coding_base_substitution,protein_base)
  (unknown_coding_base_substitution,protein_base)
  ;  
  
  
  void BaseSubstitutionEffects::initialize_from_sequence(cReferenceSequences& ref_seq_info) 
  {    
    bool count_synonymous_stop_codons = true;
    bool verbose = false;
    
    map<string,string> codon_synonymous_changes;
    map<string,string> codon_nonsynonymous_changes;
    map<string,string> codon_nonsense_changes;
    
    map<string,string> nonsynonymous_mutations;
    map<string,string> synonymous_mutations;
    
    uint32_t total_num_synonymous_changes = 0;
    uint32_t total_num_nonsynonymous_changes = 0;
    uint32_t total_num_nonsense_changes = 0;
    uint32_t total_codon_nt_positions = 0;
    uint32_t total_nt_position = 0;
    
    map<char,uint32_t> total_bases = make_map<char,uint32_t>('A',0)('T',0)('C',0)('G',0);
    
    uint32_t total_codons = 0;
    uint32_t total_orfs = 0;
    
    // Load sequence
    for(vector<cAnnotatedSequence>::iterator it=ref_seq_info.begin(); it!=ref_seq_info.end(); ++it) {
      cAnnotatedSequence& seq = *it;
      for(size_t i=1; i<=seq.get_sequence_length(); ++i) {
        char base = seq.get_sequence_1(i);
        if ((base != 'A') && (base != 'T') && (base != 'C') && (base != 'G'))
          cerr << "WARNING: Nonstandard base in sequence:" << base << "\n"; 
        else
          total_bases[base]++;
      }
      
      total_nt_position += seq.get_sequence_length();
      
      // Allocate the entire thing as default intergenic    
      SequenceBaseSubstitutionEffects& seq_bse = m_bse[seq.m_seq_id];
      seq_bse.resize(seq.get_sequence_length()*4, intergenic_base_substitution);
      
      // But set the bases that are no change to no_change
      for (uint32_t this_location_0 = 0; this_location_0 < seq.get_sequence_length(); ++this_location_0)
        for (size_t b=0; b<BaseSubstitutionEffects::base_char_list.size(); b++)
          if (BaseSubstitutionEffects::base_char_list[b] == seq.get_sequence_1(this_location_0+1) )
            seq_bse[this_location_0*4+b] = max(seq_bse[this_location_0*4+b], no_change_base_substitution);
      
      SequenceBaseCDSStrands& seq_bcs = m_bcs[seq.m_seq_id];
      seq_bcs.resize(seq.get_sequence_length(), no_CDS);
      
      for(cSequenceFeatureList::iterator it2=seq.m_features.begin(); it2!=seq.m_features.end(); ++it2) {
        cSequenceFeature& f = **it2;
        
        // Catches tRNA/rRNA/pseudogenes...
        // Things within introns are still called INTERGENIC
        if (f["type"] == "gene") {
          for (cFeatureLocationList::iterator it3=f.m_locations.begin(); it3!=f.m_locations.end(); ++it3) {
            cFeatureLocation& loc = *it3;
            
            if (verbose) cout << f.SafeGet("name") << " " << loc.get_start_1() << " " << loc.get_end_1() << " " << loc.get_strand() << endl;
            
            for (int32_t this_location_1=loc.get_start_1(); this_location_1<=loc.get_end_1(); this_location_1++) {
              int32_t this_location_0 = this_location_1-1;
              for (size_t b=0; b<BaseSubstitutionEffects::base_char_list.size(); b++)
                seq_bse[this_location_0*4+b] = max(seq_bse[this_location_0*4+b], f.m_pseudo ? pseudogene_base_substitution : noncoding_base_substitution);
            }
          }
        }
        
        // Remainder is only for coding sequences (and not pseudogenes)
        //
        if ((f["type"] != "CDS") || f.m_pseudo)
          continue;
        
        // We cannot have both an indeterminate start and end for a CDS
        // This should never happen because they are marked pseudo
        ASSERT(!(f.start_is_indeterminate() && f.end_is_indeterminate()), "CDS with indetermine start and end cannot be translated:" + f["locus_tag"]);
        
        
        // initialize gene structure
        cGeneFeature gene(f);
        total_orfs++;
        
        string gene_nt_sequence = gene.get_nucleotide_sequence(seq);
        
        // The position within a codon... indexed to start at 0.
        size_t indeterminate_codon_pos_offset_0 = 0;
        
        // Add padding to put us in-frame if we have an indeterminate start
        if (gene.start_is_indeterminate()) {
          indeterminate_codon_pos_offset_0 =  (3 - gene_nt_sequence.length() % 3) % 3;
          gene_nt_sequence = repeat_char('N', indeterminate_codon_pos_offset_0) + gene_nt_sequence;
        }
        if (gene.end_is_indeterminate()) {
          gene_nt_sequence += repeat_char('N', gene_nt_sequence.length() % 3);
        }
        
        // The position within a codon (bounds: 0-2)
        size_t on_codon_pos_0 = indeterminate_codon_pos_offset_0;
        // Positition withing gene_nt_sequence
        size_t on_nt_pos_0 = on_codon_pos_0;
        
        cFeatureLocationList& sub_locations = gene.m_locations;
        vector<int32_t> this_codon_locations_0(3, numeric_limits<int32_t>::max()); // 0-indexed
        vector<int8_t> this_codon_strands(3, 0);

        // This length includes any incomplete codons
        // at the end or beginning due to indeterminate codons
        // because we use it to bound our traversal of the sequence
        uint32_t total_amino_acid_length = gene_nt_sequence.size() / 3;
        
        // The number of the amino acid / codon indexed to start at 1
        uint32_t on_codon_number_1 = 1;
        string this_codon = "NNN";
        
        cFeatureLocationList::iterator it3=f.m_locations.begin();
        cFeatureLocation& loc = *it3;
        int8_t strand = loc.get_strand();
        int32_t pos_1 = loc.get_strand_aware_initial_position_1();
        
        // Decremented until zero to determine when to move to next location
        size_t location_position_count_down = loc.get_end_1() - loc.get_start_1() + 1;

        while (on_codon_number_1 <= total_amino_acid_length) {
          
          int32_t pos_0 = pos_1 - 1;
          
          // Prevent going out of range at the end of indeterminate genes
          if ((pos_1 >=1) && (pos_1 <= static_cast<int32_t>(seq.get_sequence_length()))) {
            //// Remember the strand of the gene overlapping this position
            if (seq_bcs[pos_0] == conflict) {
              // do nothing
            }
            // Don't count if we have genes on both strands overlapping same nucleotide
            else if (seq_bcs[pos_0] != no_CDS) {
              if ((loc.get_strand() == +1) && (seq_bcs[pos_0] == reverse) )
                seq_bcs[pos_0] = conflict;
              if ((loc.get_strand() == -1) && (seq_bcs[pos_0] == forward) )
                seq_bcs[pos_0] = conflict;
            }
            else {
              seq_bcs[pos_0] = (loc.get_strand() == +1 ? forward : reverse);
            }
            
            //// Handle codon synonymous/nonsynonymous changes
            this_codon_locations_0[on_codon_pos_0] = pos_0;
            this_codon[on_codon_pos_0] = gene_nt_sequence[on_nt_pos_0];
            this_codon_strands[on_codon_pos_0] = strand;
          }
          
          on_codon_pos_0++;
          
          // The codon is filled, now make all mutations and assign to proper nucleotides
          if (on_codon_pos_0 == 3) {
            
            // The adjustment to codon number is so that we don't count
            // the first codon of an indeterminate start as a start codon!
            char original_amino_acid = cReferenceSequences::translate_codon(this_codon, gene.translation_table, ( gene.start_is_indeterminate() && (on_codon_number_1 == 1) ) ? 2 : on_codon_number_1, gene.get_locus_tag());
            
            for (int32_t test_codon_index=0; test_codon_index<3; test_codon_index++) {
              
              // Only do locations that are within range
              if (this_codon_locations_0[test_codon_index] != numeric_limits<int32_t>::max()) {
              
                // Check range to avoid problems writing to random memory
                ASSERT( (this_codon_locations_0[test_codon_index] >= 0) && (this_codon_locations_0[test_codon_index] < static_cast<int32_t>(seq.get_sequence_length())), "Position within gene is out of range: " + gene.get_locus_tag() );
                
                for (size_t b=0; b<BaseSubstitutionEffects::base_char_list.size(); b++) {
                  
                  char mut_base = BaseSubstitutionEffects::base_char_list[b];
                  
                  // We have to complement the base we are changing in the codon if this
                  // part of the reading frame was on the reverse genomic strand!
                  if (this_codon_strands[test_codon_index] == -1)
                    mut_base = complement_base_char(mut_base);
                  
                  string test_codon = this_codon;
                  test_codon[test_codon_index] = mut_base;
                  
                  char mut_amino_acid = cReferenceSequences::translate_codon(test_codon, gene.translation_table, ( gene.start_is_indeterminate() && (on_codon_number_1 == 1) ) ? 2 : on_codon_number_1, gene.get_locus_tag());
                  
                  // We are testing whether we defined this to avoid going out of position due to
                  // indeterminate coordinates
                  
                  
                  if ((mut_amino_acid == '?') || (original_amino_acid == '?'))
                    seq_bse[this_codon_locations_0[test_codon_index]*4+b] = max(seq_bse[this_codon_locations_0[test_codon_index]*4+b], unknown_coding_base_substitution);
                  else if (mut_amino_acid == original_amino_acid)
                    seq_bse[this_codon_locations_0[test_codon_index]*4+b] = max(seq_bse[this_codon_locations_0[test_codon_index]*4+b], synonymous_coding_base_substitution);
                  else if (mut_amino_acid == '*')
                    seq_bse[this_codon_locations_0[test_codon_index]*4+b] = max(seq_bse[this_codon_locations_0[test_codon_index]*4+b], nonsense_coding_base_substitution);
                  else
                    seq_bse[this_codon_locations_0[test_codon_index]*4+b] = max(seq_bse[this_codon_locations_0[test_codon_index]*4+b], nonsynonymous_coding_base_substitution);
                }
              }
            }
            
            // reset
            on_codon_pos_0 = 0;
            this_codon = "NNN";
            on_codon_number_1++;
            vector<int8_t> this_codon_strands(3, 0);
            this_codon_locations_0 = vector<int32_t>(3, numeric_limits<int32_t>::max()); // 0-indexed
          }
          
          // Move to the next location
          location_position_count_down--;
          on_nt_pos_0++;
        
          // we have to allow the countdown to go past the end of the
          // last location for indeterminate end features
          if (location_position_count_down == 0) {
            it3++;
            if (it3 != f.m_locations.end()) {
              loc = *it3;
              strand = loc.get_strand();
              location_position_count_down = loc.get_end_1() - loc.get_start_1() + 1;
              pos_1 = loc.get_strand_aware_initial_position_1();
              continue; // skips moving position below
            }
          }
          
          // otherwise increment
          pos_1 += strand;
          
        } // end location within gene loop

      
      } // end feature loop
      
      if (verbose) {
        for(size_t i=1; i<=seq.get_sequence_length(); ++i) {
          const size_t i_0 = i-1;
          char base = seq.get_sequence_1(i);
          cout << i << "\t" << base << "\t" << seq_bcs[i_0] << "\t"
          << seq_bse[i_0*4+0] << "\t" << seq_bse[i_0*4+1] << "\t" << seq_bse[i_0*4+2] << "\t" << seq_bse[i_0*4+3] << endl;
        }
      }
      
    } // end sequence loop
  }    
  
  void BaseSubstitutionEffectCounts::initialize_possible_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse)
  {
    vector<string> seq_ids = ref_seq_info.seq_ids();
    for (vector<string>::iterator seq_id_it = seq_ids.begin(); seq_id_it != seq_ids.end(); ++seq_id_it) {
      for (uint32_t i=1; i <= ref_seq_info[*seq_id_it].get_sequence_length(); i++) {
        //cout << *seq_id_it << " " << i << endl;
        change_position_1_possible_totals(ref_seq_info, bse, *seq_id_it, i, +1);
      }
    }
  }
   
  void BaseSubstitutionEffectCounts::change_position_1_possible_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse, string seq_id, uint32_t pos_1, int32_t inc)
  {
    BaseSubstitutionEffectPositionInfo pos_info = bse.position_info_1(ref_seq_info, seq_id, pos_1);
   
    // counts of the total number of bases in the genome in category
    m_base_counts[pos_info.m_base_type]["nt"] += inc;
    if ( (pos_info.m_base_char == 'G') || (pos_info.m_base_char == 'C') )
      m_base_counts[pos_info.m_base_type]["gc"] += inc;
    else
      m_base_counts[pos_info.m_base_type]["at"] += inc;
        
    // counts of the total number of base changes in the genome in category
    for (size_t this_base_index = 0; this_base_index < BaseSubstitutionEffects::base_char_list.size(); ++this_base_index) {
      
      base_char this_base_char = BaseSubstitutionEffects::base_char_list[this_base_index];
      
      // Don't try to count if the base is the same -- lookup of base pair change will fail
      if (pos_info.m_base_char == this_base_char) 
        continue;
      
      string base_key = pos_info.m_base_char + BaseSubstitutionEffects::separator + this_base_char;
      string base_pair_key = BaseSubstitutionEffects::base_change_to_base_pair_change[base_key];
      
      string base_change_type_key = BaseSubstitutionEffects::base_change_type_list[pos_info.m_base_substitution_effect[this_base_index]];
      
      m_possible_base_pair_change_counts[base_change_type_key][base_pair_key] += inc;
      m_possible_base_pair_change_counts[base_change_type_key]["TOTAL"] += inc;
      
      m_possible_base_pair_change_counts["TOTAL"][base_pair_key] += inc;
      m_possible_base_pair_change_counts["TOTAL"]["TOTAL"] += inc;
    }
   
  }
   
  void BaseSubstitutionEffectCounts::change_position_1_observed_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse, string seq_id, uint32_t pos_1, string new_base, int32_t inc)
  {
    ASSERT(new_base.size()==1,"Unexpected base string size.");
    base_char new_base_char = new_base[0];
    BaseSubstitutionEffectPositionInfo pos_info = bse.position_info_1(ref_seq_info, seq_id, pos_1);
    string base_change = pos_info.m_base_char + BaseSubstitutionEffects::separator + new_base;
    
    size_t new_base_index = BaseSubstitutionEffects::base_char_to_base_index[new_base_char];
    string base_change_type_key = BaseSubstitutionEffects::base_change_type_list[pos_info.m_base_substitution_effect[new_base_index]];
    
    ASSERT(pos_info.m_base_substitution_effect[new_base_index] != no_change_base_substitution, "Attempt to count base substitution that does not change base");
    
    string base_pair_key = BaseSubstitutionEffects::base_change_to_base_pair_change[base_change];	
    
    m_observed_base_pair_change_counts[base_change_type_key][base_pair_key] += inc;
    m_observed_base_pair_change_counts[base_change_type_key]["TOTAL"] += inc;
    
    m_observed_base_pair_change_counts["TOTAL"][base_pair_key] += inc;
    m_observed_base_pair_change_counts["TOTAL"]["TOTAL"] += inc;
    
  }
   
  
  BaseSubstitutionEffectPositionInfo 
  BaseSubstitutionEffects::position_info_1(
                                           cReferenceSequences& ref_seq_info, 
                                           string seq_id, 
                                           uint32_t pos_1 
                                           )
  {
    uint32_t pos_0 = pos_1 - 1;
    BaseSubstitutionEffectPositionInfo pos_info;
   
    pos_info.m_base_char = ref_seq_info[seq_id].get_sequence_1(pos_1);
    
    SequenceBaseSubstitutionEffects::iterator bse_it = m_bse[seq_id].begin() + pos_0 * 4;
    

    for (size_t i=0; i<4; i++) {
      BaseSubstitutionEffect this_bse = *bse_it;
      if (this_bse != no_change_base_substitution) {
        pos_info.m_base_type = snp_type_to_base_type[this_bse];
        break;
      }
      ++bse_it;
    }
    
    bse_it = m_bse[seq_id].begin() + pos_0 * 4;
    
    /*
    cout << BaseSubstitutionEffects::base_change_type_list[m_bse[seq_id][pos_0 * 4]] << endl;
    cout << BaseSubstitutionEffects::base_change_type_list[m_bse[seq_id][pos_0 * 4+1]] << endl;
    cout << BaseSubstitutionEffects::base_change_type_list[m_bse[seq_id][pos_0 * 4+2]] << endl;
    cout << BaseSubstitutionEffects::base_change_type_list[m_bse[seq_id][pos_0 * 4+3]] << endl;
    */
    
    pos_info.m_base_substitution_effect.insert(
                                              pos_info.m_base_substitution_effect.begin(),
                                              bse_it, bse_it+4
                                              );

    
    
    SequenceBaseCDSStrands::iterator bcs_it = m_bcs[seq_id].begin() + pos_0;
    pos_info.m_base_cds_strand = *bcs_it;
    
    return pos_info;
  }
  
  

  
  void MutationCountFile(
                         cReferenceSequences& ref_seq_info, 
                         vector<cGenomeDiff>& genome_diffs, 
                         string& output_file_name,
                         string& detailed_output_file_name,
                         bool base_substitution_statistics,
                         bool count_polymorphisms,
                         bool calculate_genome_size,
                         bool verbose
                         )
  {
      
    // Figure out the names of all "repeat" columns
    map<string,bool> mob_name_hash;
    map<string,bool> mob_mediated_name_hash;
    map<string,bool> con_mediated_name_hash;
    
    // create a copy so we don't alter the original list
    vector<cGenomeDiff> sorted_genome_diffs(genome_diffs);
    cGenomeDiff::sort_gd_list_by_treatment_population_time(sorted_genome_diffs);
        
    for (vector<cGenomeDiff>::iterator it=sorted_genome_diffs.begin(); it != sorted_genome_diffs.end(); ++it) {
      cGenomeDiff &gd = *it;
      
      diff_entry_list_t muts = gd.mutation_list();
      for (diff_entry_list_t::iterator it=muts.begin(); it != muts.end(); ++it) {
        cDiffEntry& mut = **it;
        if (mut._type == MOB) {
          mob_name_hash[mut["repeat_name"]] = true;
        }
        if (mut._type == CON) {
          if (mut.entry_exists("mediated"))
            con_mediated_name_hash[mut["mediated"]] = true;
        }
        if ( (mut._type == DEL) || (mut._type == AMP) || (mut._type == SUB) || (mut._type == INS)) {
          if (mut.entry_exists("mediated"))
            mob_mediated_name_hash[mut["mediated"]] = true;
        }
      }
    }
    vector<string> mob_name_list = map_keys_to_list<string,bool>(mob_name_hash);
    sort(mob_name_list.begin(), mob_name_list.end());
    vector<string> mob_mediated_name_list = map_keys_to_list<string,bool>(mob_mediated_name_hash);
    sort(mob_mediated_name_list.begin(), mob_mediated_name_list.end());
    vector<string> con_mediated_name_list = map_keys_to_list<string,bool>(con_mediated_name_hash);
    sort(con_mediated_name_list.begin(), con_mediated_name_list.end());
    
    
    // Handle metadata
    bool has_metadata_treatment(false);
    bool has_metadata_population(false);
    bool has_metadata_time(false);
    bool has_metadata_clone(false);
    
    // Handle other #=ITEM metadata things in the header in a general way
    set<string> other_metadata_headers;
    
    for (vector<cGenomeDiff>::iterator it=sorted_genome_diffs.begin(); it != sorted_genome_diffs.end(); ++it) {
      cGenomeDiff &gd = *it;
      
      if (gd.metadata.treatment.size()) has_metadata_treatment = true;
      if (gd.metadata.population.size()) has_metadata_population = true;
      if (gd.metadata.time != -1.0) has_metadata_time = true;
      if (gd.metadata.clone.size()) has_metadata_clone = true;

      for (map<string,string>::iterator itd=gd.metadata.breseq_data.begin(); itd != gd.metadata.breseq_data.end(); ++itd) {
        const string &key = to_upper(itd->first);
        const string &value = itd->second;
        other_metadata_headers.insert(key);
      }
    }
    
    vector<string> column_headers;
    column_headers.push_back("file");
    column_headers.push_back("sample");
    if (has_metadata_treatment) column_headers.push_back("treatment");
    if (has_metadata_population) column_headers.push_back("population");
    if (has_metadata_time) column_headers.push_back("time");
    if (has_metadata_clone) column_headers.push_back("clone");
    // Handle other #=ITEM metadata things in the header in a general way
    for (set<string>::iterator it=other_metadata_headers.begin(); it!=other_metadata_headers.end(); it++) {
      column_headers.push_back(*it);
    }
    column_headers.push_back("total");
    column_headers.push_back("base_substitution");
    column_headers.push_back("small_indel");
    column_headers.push_back("large_deletion");
    column_headers.push_back("large_insertion");
    column_headers.push_back("large_amplification");
    column_headers.push_back("large_substitution");
    column_headers.push_back("mobile_element_insertion");
    column_headers.push_back("gene_conversion");
    column_headers.push_back("inversion");
    if (calculate_genome_size) {
      column_headers.push_back("changed_bp");
      column_headers.push_back("deleted_bp");
      column_headers.push_back("inserted_bp");
      column_headers.push_back("final_size_bp");
    }
    column_headers.push_back("reference_bp");
    column_headers.push_back("called_reference_bp");
    
    vector<string> header_snp_types = prefix_each_in_vector(snp_types, "base_substitution.");
    column_headers.insert(column_headers.end(),header_snp_types.begin(), header_snp_types.end());
    
    vector<string> header_mob_name_list = prefix_each_in_vector(mob_name_list, "mobile_element.");
    column_headers.insert(column_headers.end(),header_mob_name_list.begin(), header_mob_name_list.end());
    
    vector<string> header_mob_mediated_name_list = prefix_each_in_vector(mob_mediated_name_list, "mobile_element.mediated.");
    column_headers.insert(column_headers.end(),header_mob_mediated_name_list.begin(), header_mob_mediated_name_list.end());
    
    vector<string> header_con_mediated_name_list = prefix_each_in_vector(con_mediated_name_list, "gene_conversion.mediated.");
    column_headers.insert(column_headers.end(),header_con_mediated_name_list.begin(), header_con_mediated_name_list.end());
    
    ofstream output_file(output_file_name.c_str());
    ASSERT(output_file.good(), "Error writing to file: " + output_file_name);
    
    BaseSubstitutionEffectCounts total_bsec;
    BaseSubstitutionEffects bse;
    
    // Populate information about the effects of every base substitution in the genome.
    if (base_substitution_statistics) {
      //uout("Calculating base substitution effects in reference sequences") << endl;
      bse.initialize_from_sequence(ref_seq_info);
      total_bsec.initialize_possible_totals(ref_seq_info, bse);
    }
    
    // Create the column headings for the detailed base substitution counts.
    if (base_substitution_statistics) {
     
     for (vector<string>::const_iterator snp_type = BaseSubstitutionEffectCounts::base_change_type_count_list.begin();
          snp_type != BaseSubstitutionEffectCounts::base_change_type_count_list.end(); ++snp_type) {
       for (vector<string>::const_iterator bp_change = BaseSubstitutionEffectCounts::base_pair_change_count_list.begin();
            bp_change != BaseSubstitutionEffectCounts::base_pair_change_count_list.end(); ++bp_change) {
         column_headers.push_back("POSSIBLE." + *snp_type + BaseSubstitutionEffects::separator + *bp_change);
       }
     }
     
      for (vector<string>::const_iterator snp_type = BaseSubstitutionEffectCounts::base_change_type_count_list.begin();
           snp_type != BaseSubstitutionEffectCounts::base_change_type_count_list.end(); ++snp_type) {
        for (vector<string>::const_iterator bp_change = BaseSubstitutionEffectCounts::base_pair_change_count_list.begin();
             bp_change != BaseSubstitutionEffectCounts::base_pair_change_count_list.end(); ++bp_change) {
          column_headers.push_back("OBSERVED." + *snp_type + BaseSubstitutionEffects::separator + *bp_change);
        }
      }
     
    } //if (base_substitution_statistics)
    
    
    bool detailed_output = detailed_output_file_name.size() != 0;
    
    // For detailed output
    vector<string> base_substitution_lines;
    vector<string> small_indel_lines;
    vector<string> large_deletion_lines;
    vector<string> large_insertion_lines;
    vector<string> large_amplification_lines;
    vector<string> large_substitution_lines;
    vector<string> mobile_element_insertion_lines;
    vector<string> gene_conversion_lines;
    vector<string> inversion_lines;
    
    output_file << join(column_headers, ",") << endl;
    
    for (vector<cGenomeDiff>::iterator it=sorted_genome_diffs.begin(); it != sorted_genome_diffs.end(); ++it) {
      cGenomeDiff &gd = *it;
      cout << "    Counting mutations " + gd.get_title() << endl << endl;
      
      vector<string> line_prefix_items;
      line_prefix_items.push_back(gd.get_file_name());
      line_prefix_items.push_back(gd.metadata.title);
      
      if (has_metadata_treatment) line_prefix_items.push_back(gd.metadata.treatment);
      if (has_metadata_population) line_prefix_items.push_back(gd.metadata.population);
      if (has_metadata_time)
        line_prefix_items.push_back((gd.metadata.time != -1.0) ? to_string<double>(gd.metadata.time) : "");
      if (has_metadata_clone) line_prefix_items.push_back(gd.metadata.clone);
      
      // Handle other #=ITEM metadata things in the header in a general way
      for (set<string>::iterator it=other_metadata_headers.begin(); it!=other_metadata_headers.end(); it++) {
        if (gd.metadata.breseq_data.find(*it) != gd.metadata.breseq_data.end()) {
          line_prefix_items.push_back(gd.metadata.breseq_data[*it]);
        } else {
          line_prefix_items.push_back("");
        }
      }
      
      
      string detailed_line_prefix = join(line_prefix_items, "\t");
      
      BaseSubstitutionEffectCounts this_bsec;
      // deep copy totals of entire sequence
      if (base_substitution_statistics)
        this_bsec = total_bsec;

      int32_t total_bp = ref_seq_info.get_total_length();
      
      // Complicated map storing a bunch of counts
      map<string,map<string,int32_t> > count;
      
      for(vector<string>::const_iterator snp_type = snp_types.begin(); snp_type != snp_types.end(); ++snp_type) {
        count["type"][*snp_type] = 0;
      }
      for(vector<string>::const_iterator mob_name = mob_name_list.begin(); mob_name != mob_name_list.end(); ++mob_name) {
        count["mob"][*mob_name] = 0;
      }
      
      count["base_substitution"][""] = 0;
      count["large_deletion"][""] = 0;
      count["small_indel"][""] = 0;
      count["large_insertion"][""] = 0;
      count["large_amplification"][""] = 0;
      count["large_substitution"][""] = 0;
      count["gene_conversion"][""] = 0;
      count["mobile_element_insertion"][""] = 0;
      count["inversion"][""] = 0;
      
      // BEGIN for each mutation
      uint32_t total_mutations(0);
      diff_entry_list_t mut_list = gd.mutation_list();
      for (diff_entry_list_t::iterator it=mut_list.begin(); it != mut_list.end(); ++it) {		
        
        cDiffEntry& mut = **it;
        
        // Don't count polymorphisms - could make it an option to partially count them
        if (!count_polymorphisms && mut.entry_exists(FREQUENCY) && from_string<double>(mut[FREQUENCY]) != 1.0) {
          if (verbose) cerr << "Skipping polymorphic: " << mut << endl;
          continue;
        }
          
        if (verbose) cerr << "Counting: " << mut << endl;
        total_mutations++;
        
        // Main count
        ASSERT(mut.entry_exists("mutation_category"), "mutation_category entry does not exist in mutation\n" + mut.as_string());
        count[mut["mutation_category"]][""]++;
        
        
        // Below we save some classes of mutations for additional output and
        // count lists of things sliced and diced in different ways.
        if (mut._type == SNP) {
          count["base_substitution"][""]++;
          if (base_substitution_statistics) {
            this_bsec.change_position_1_observed_totals(ref_seq_info, bse, mut[SEQ_ID], from_string<int32_t>(mut[POSITION]), mut[NEW_SEQ], +1);					
          }
          
          // @JEB 2017-12-14 new code to account for multiple SNP types if it overlaps multiple genes
          // Return the highest category of base substitution effect, according to the normal heirarchy

          vector<string> _snp_type_list = split(mut["snp_type"], cReferenceSequences::multiple_separator);
          BaseSubstitutionEffect bse(intergenic_base_substitution);
          for (vector<string>::iterator it=_snp_type_list.begin(); it != _snp_type_list.end(); it++) {
            bse = max(string_to_bse(*it), bse);
          }
          string _snp_type = bse_to_string(bse);
          
          count["type"][_snp_type]++;
          base_substitution_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
        
        } else if (mut._type == DEL) {
          
          if (mut.entry_exists("mediated"))
            count["mob_mediated"][mut["mediated"]]++;
          
          if ( mut["mutation_category"] == "large_deletion" ) {
            large_deletion_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          } else {
            small_indel_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          }
          
        } else if (mut._type == INS) {
          
          if (mut.entry_exists("mediated"))
            count["mob_mediated"][mut["mediated"]]++;
          
          if ( mut["mutation_category"] == "large_insertion" ) {
            large_insertion_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          } else {
            small_indel_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          }
          
        } else if (mut._type == SUB) {
          
          if (mut.entry_exists("mediated"))
            count["mob_mediated"][mut["mediated"]]++;
          
          
          if ( mut["mutation_category"] == "large_substitution" ) {
            large_substitution_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          } else {
            small_indel_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          }
          
        } else if (mut._type == CON) {
          
          gene_conversion_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
        
        } else if (mut._type == MOB) {
          int32_t rpos = -1;
          // This includes the MOB and any ins/del start/end adjustments, NOT duplication size
          string mob_region;
          string repeat_seq = ref_seq_info.repeat_family_sequence(mut["repeat_name"], 1, mut.entry_exists("mob_region") ? &mut["mob_region"] : NULL);
          
          int32_t this_length = repeat_seq.size();
          this_length += from_string<int32_t>(mut["duplication_size"]);
          
          count["mob"][mut["repeat_name"]]++;
          mobile_element_insertion_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          
        } else if (mut._type == AMP) {
          int32_t this_size = from_string<uint32_t>(mut[SIZE]) * (from_string<uint32_t>(mut["new_copy_number"]) - 1);
          
          if (mut.entry_exists("mediated"))
            count["mob_mediated"][mut["mediated"]]++;
          
          if ( mut["mutation_category"] == "large_amplification" ) {
            large_amplification_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          } else {
            small_indel_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          }
        
        } else if (mut._type == INV) {
          
          inversion_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          
        } else {
          ERROR("Could not count mutation:\n" + mut.as_string());
        }
      }
      
      // statistics for UN
      int32_t un_bp = 0;
      
      diff_entry_list_t un_list = gd.get_list(make_vector<gd_entry_type>(UN));
      for (diff_entry_list_t::iterator it=un_list.begin(); it!= un_list.end(); ++it) {
        cDiffEntry& un = **it;
        un_bp += from_string<int32_t>(un[END]) - from_string<int32_t>(un[START]) + 1;
        
        // subtract these positions from the possible observations of base pair statistics
        if (base_substitution_statistics)
          for (uint32_t pos_1 = from_string<uint32_t>(un[START]); pos_1 <= from_string<uint32_t>(un[END]); pos_1++)
            this_bsec.change_position_1_possible_totals(ref_seq_info, bse, un[SEQ_ID], pos_1, -1);					
        
      }
      // END for each mutation

      // USE APPLY TO CALCULATE SIZE CHANGES!!!
      // MUST BE DONE AFTER OTHER COUNTING BECAUSE IT CHANGES THE GENOME DIFF!!
      
      if (calculate_genome_size) {
        cReferenceSequences new_ref_seq_info = cReferenceSequences::deep_copy(ref_seq_info);
        gd.apply_to_sequences(ref_seq_info, new_ref_seq_info, false, 20);
      }
        
      int32_t called_bp = total_bp - un_bp;
      
      vector<string> this_columns = line_prefix_items;
      this_columns.push_back(to_string(total_mutations));
      this_columns.push_back(to_string(count["base_substitution"][""]));
      this_columns.push_back(to_string(count["small_indel"][""]));
      this_columns.push_back(to_string(count["large_deletion"][""]));
      this_columns.push_back(to_string(count["large_insertion"][""]));
      this_columns.push_back(to_string(count["large_amplification"][""]));
      this_columns.push_back(to_string(count["large_substitution"][""]));
      this_columns.push_back(to_string(count["mobile_element_insertion"][""]));
      this_columns.push_back(to_string(count["gene_conversion"][""]));
      this_columns.push_back(to_string(count["inversion"][""]));
      if (calculate_genome_size) {
        this_columns.push_back(gd.get_breseq_data("BASES-CHANGED"));
        this_columns.push_back(gd.get_breseq_data("BASES-DELETED"));
        this_columns.push_back(gd.get_breseq_data("BASES-INSERTED"));
        this_columns.push_back(gd.get_breseq_data("GENOME-SIZE-FINAL"));
      }
      this_columns.push_back(to_string(total_bp));
      this_columns.push_back(to_string(called_bp));
      
      vector<string> snp_type_counts = map_key_list_to_values_as_strings(count["type"], snp_types);
      this_columns.insert(this_columns.end(),snp_type_counts.begin(), snp_type_counts.end());
      
      vector<string> mob_type_counts = map_key_list_to_values_as_strings(count["mob"], mob_name_list);
      this_columns.insert(this_columns.end(),mob_type_counts.begin(), mob_type_counts.end());
      
      vector<string> mob_mediated_type_counts = map_key_list_to_values_as_strings(count["mob_mediated"], mob_mediated_name_list);
      this_columns.insert(this_columns.end(),mob_mediated_type_counts.begin(), mob_mediated_type_counts.end());
      
      vector<string> con_mediated_type_counts = map_key_list_to_values_as_strings(count["con_mediated"], con_mediated_name_list);
      this_columns.insert(this_columns.end(),con_mediated_type_counts.begin(), con_mediated_type_counts.end());
      
      
      if (base_substitution_statistics) {	
        
        for(vector<string>::iterator this_base_change_type = BaseSubstitutionEffectCounts::base_change_type_count_list.begin();
            this_base_change_type != BaseSubstitutionEffectCounts::base_change_type_count_list.end(); ++this_base_change_type) {
                    
          for(vector<string>::iterator this_base_pair_change = BaseSubstitutionEffectCounts::base_pair_change_count_list.begin();
              this_base_pair_change != BaseSubstitutionEffectCounts::base_pair_change_count_list.end(); ++this_base_pair_change) {
          
            this_columns.push_back( s( this_bsec.m_possible_base_pair_change_counts[*this_base_change_type][*this_base_pair_change] ) );
          }
          
        }
        
        for(vector<string>::iterator this_base_change_type = BaseSubstitutionEffectCounts::base_change_type_count_list.begin();
            this_base_change_type != BaseSubstitutionEffectCounts::base_change_type_count_list.end(); ++this_base_change_type) {
          
          for(vector<string>::iterator this_base_pair_change = BaseSubstitutionEffectCounts::base_pair_change_count_list.begin();
              this_base_pair_change != BaseSubstitutionEffectCounts::base_pair_change_count_list.end(); ++this_base_pair_change) {
            
            this_columns.push_back( s( this_bsec.m_observed_base_pair_change_counts[*this_base_change_type][*this_base_pair_change] ) );
          }          
        }
      } // if (base_substitution_statistics)
        
      output_file << join(this_columns, ",") << endl;
    }
    
    
    // Write the detailed output file
    if (detailed_output) {
      ofstream detailed_output_file(detailed_output_file_name.c_str());
      
      detailed_output_file << "BASE SUBSTITUTIONS: " << base_substitution_lines.size() << "\n";
      detailed_output_file << join(base_substitution_lines, "\n") << "\n\n";
      
      detailed_output_file << "SMALL INDELS: " << small_indel_lines.size() << "\n";
      detailed_output_file << join(small_indel_lines, "\n")<< "\n\n";
      
      detailed_output_file << "LARGE DELETIONS: " << large_deletion_lines.size() << "\n";
      detailed_output_file << join(large_deletion_lines, "\n")<< "\n\n";
      
      detailed_output_file << "LARGE INSERTIONS: " << large_insertion_lines.size() << "\n";
      detailed_output_file << join(large_insertion_lines, "\n")<< "\n\n";
      
      detailed_output_file << "LARGE AMPLIFICATIONS: " << large_amplification_lines.size() << "\n";
      detailed_output_file << join(large_amplification_lines, "\n")<< "\n\n";
      
      detailed_output_file << "LARGE SUBSTITUTIONS: " << large_substitution_lines.size() << "\n";
      detailed_output_file << join(large_substitution_lines, "\n")<< "\n\n";
      
      detailed_output_file << "MOBILE ELEMENT INSERTIONS: " << mobile_element_insertion_lines.size() << "\n";
      detailed_output_file << join(mobile_element_insertion_lines, "\n") << "\n\n";
      
      detailed_output_file << "GENE CONVERSIONS: " << gene_conversion_lines.size() << "\n";
      detailed_output_file << join(gene_conversion_lines, "\n")<< "\n";
      
      detailed_output_file << "INVERSIONS: " << inversion_lines.size() << "\n";
      detailed_output_file << join(inversion_lines, "\n")<< "\n";
    }
    
  }
  

} // namespace breseq
