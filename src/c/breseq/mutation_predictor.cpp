
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

#include "libbreseq/mutation_predictor.h"

#include "libbreseq/output.h"
#include "libbreseq/identify_mutations.h"
#include "libbreseq/resolve_alignments.h"

using namespace std;

namespace breseq {

  cReferenceSequences MutationPredictor::ref_seq_info;
  
	MutationPredictor::MutationPredictor(cReferenceSequences& _ref_seq_info)
	{
		 ref_seq_info = _ref_seq_info;
	}

	// Private methods

	cSequenceFeature* MutationPredictor::within_repeat(string seq_id, int32_t position)
	{
		cSequenceFeatureList& repeat_list = ref_seq_info[seq_id].m_repeats;
    cSequenceFeature* repeat= NULL;
    
    // by returning the last one we encounter that we are inside, 
    // we get the inner repeat in nested cases
    for(cSequenceFeatureList::iterator it = repeat_list.begin(); it != repeat_list.end(); it++) {
      cSequenceFeaturePtr& test_repeat = *it;
			if ((test_repeat->m_location.get_start_1() <= position) && (position <= test_repeat->m_location.get_end_1()))
				repeat = test_repeat.get();
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

	// look at SNPs and small indels predicted by read alignments.
	// be sure they are sorted by position
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
   */
  void MutationPredictor::prepare_junctions(Settings& settings, Summary& summary, cGenomeDiff& gd)
  {
    (void) settings;
    (void) summary;
    
    // For all that follows, we need information about repeat_regions overlapping the sides of junctions
    vector<gd_entry_type> jc_types = make_vector<gd_entry_type>(JC);
		diff_entry_list_t jc = gd.list(jc_types);
    
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
				cSequenceFeaturePtr is = ref_seq_info.find_closest_repeat_region_boundary(
                                                                                  n(j[side_key + "_position"]),
                                                                                  ref_seq_info[j[side_key + "_seq_id"]].m_repeats,
                                                                                  this_max_distance_to_repeat,
                                                                                  n(j[side_key + "_strand"])
                                                                                  );
				if (is.get() != NULL)
				{
					j["_" + side_key + "_is"] = "1";
					j["_" + side_key + "_is_start"] = s(is->m_location.get_start_1());
					j["_" + side_key + "_is_end"] = s(is->m_location.get_end_1());
          j["_" + side_key + "_is_name"] = (*is)["name"];
          j["_" + side_key + "_is_strand"] = s(is->m_location.get_strand());
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
    
    int32_t max_read_length = summary.sequence_conversion.max_read_length;
    
		// DEL prediction:
		// (1) there is a junction that exactly crosses the deletion boundary
		// (2) there is is no junction, but both ends of the deletion are in repeat sequences
		// (3) there is a junction between unique sequence and a repeat element
    
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
			mut
      ("seq_id", mc_item["seq_id"])
      ("position", mc_item["start"])
      ("size", s(n(mc_item["end"]) - n(mc_item["start"]) + 1));
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
      bool done = false;
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
          
          //it's possible that one or both sides are in repeat elements
          cSequenceFeature* r1_pointer = within_repeat(jc_item["side_1_seq_id"], side_1_position);
          cSequenceFeature* r2_pointer = within_repeat(jc_item["side_2_seq_id"], side_2_position);
          
          // one repeat cases where the end matches up exactly
          if (r1_pointer) 
          {
            // must match up to an end of the repeat
            if (side_1_position == static_cast<int32_t>(r1_pointer->m_location.get_start_1())
                || (side_1_position == static_cast<int32_t>(r1_pointer->m_location.get_end_1())))
            {
              mut["mediated"] = (*r1_pointer)["name"];
            }
          }
          // if it didn't match, then check possibility of a second repeat
          if (!mut.count("mediated") && r2_pointer) 
          {
            // must match up to an end of the repeat
            if ((side_2_position == static_cast<int32_t>(r2_pointer->m_location.get_start_1())
                 || (side_2_position) == static_cast<int32_t>(r2_pointer->m_location.get_end_1())))
            {
              mut["mediated"] = (*r2_pointer)["name"];
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
      
      
      cSequenceFeature* r1_pointer = within_repeat(mut["seq_id"], n(mut["position"]));
      cSequenceFeature* r2_pointer = within_repeat(mut["seq_id"], n(mut["position"]) + n(mut["size"]));
      
			///
			// (2) there is no junction, but both ends of the deletion are in different copies of the same repeat sequence
			///
      
			// Then we will adjust the coordinates to remove...
			if  (
           (r1_pointer != r2_pointer) 
           && (r1_pointer != NULL) 
           && (r2_pointer != NULL) 
           && ((*r1_pointer)["name"] == (*r2_pointer)["name"])
           )
			{
				cSequenceFeature& r1 = *r1_pointer, r2 = *r2_pointer;
        
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
				uint32_t slop_distance = max_read_length;
        
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
        if (n(mut["size"]) != 0) {
          // remember the name of the element
          mut["between"] = r1["name"];
          gd.add(mut);
          
          mc_it = mc.erase(mc_it); // iterator is now past the erased element
          mc_it--;                //We just removed the current jc, do not iterate.
          
          if (verbose)
            cout << "**** Ends of junction in copies of same repeat element ****\n";
          continue; // to next mc_item
        }
			}
      
			// Both sides were unique or redundant, nothing more we can do...
			if ( (r1_pointer == NULL && r2_pointer == NULL) || (r1_pointer != NULL && r2_pointer != NULL) )
				continue; // to next mc_item
      
			///
			// (3) there is a junction between unique sequence and a repeat element
			///
			cSequenceFeature& r = (r1_pointer != NULL) ? *r1_pointer : *r2_pointer;
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
					cout << "Check 2: " << r["name"] << " ne " << j["_" + j["_is_interval"] + "_is_name"] << endl;
				if (r["name"] != j["_" + j["_is_interval"] + "_is_name"])
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
        if (n(mut["size"]) == 0)
          continue;
        
				// OK, we're good!
				mut["mediated"] = r["name"];
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

  }
  
  void MutationPredictor::predictJCplusJCtoMOB(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc)
  {
    (void)summary;
    bool verbose = false;
    
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
      
			// Compile a list of the next possibilities within a certain length of bases
      vector<diff_entry_ptr_t> j2_list;
			vector<diff_entry_list_t::iterator> it_delete_list_2;
      
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
        
        it_delete_list_2.push_back(jc2_it);        
				j2_list.push_back(*jc2_it);
			}
      if (verbose)
        cout << "Size of J2 list: " << j2_list.size() << endl;
      
			//sort the $j2_list by reject reason and score
      
			sort(j2_list.begin(), j2_list.end(), sort_by_reject_score);
      
			// We need to go through all with the same coordinate (or within a certain coordinate stretch?)
			// because sometimes a failed junction will be in between the successful junctions
      for(size_t i=0; i<j2_list.size(); i++)
			{
        // we may have already used this j1
        if (jc1_erased) break;
        
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
        
				// Remove these predictions from the list
        // it's important to remove jc1_it LAST so that we don't accidentally invalidate
        // this iterator by erasing what it points to again with the it_delete_list_2 call
				jc.erase(it_delete_list_2[i]);
        jc1_it = jc.erase(jc1_it); // iterator is now past element erased
        jc1_erased = true; // this tells us to not increment the list

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
        if (!j1.entry_exists("user_defined")) {
          
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
        }
        
        bool j2_is_ambiguous = false;
        if (!j2.entry_exists("user_defined")) {
          
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
        
        if (j1.entry_exists(SIDE_1_READ_COUNT) && j1.entry_exists(SIDE_2_READ_COUNT) & j1.entry_exists(NEW_JUNCTION_READ_COUNT)
         && j2.entry_exists(SIDE_1_READ_COUNT) && j2.entry_exists(SIDE_2_READ_COUNT) & j2.entry_exists(NEW_JUNCTION_READ_COUNT)) {

        
          int32_t require_overlap = n(mut["duplication_size"]);
          
          
          if (verbose) cerr << "Before 1:" << endl << j1 << endl;
          assign_one_junction_read_counts(settings, summary, j1, require_overlap);
          j1["read_count_offset"] = mut["duplication_size"];
          if (verbose) cerr << "After 1:" << endl << j1 << endl;
          
          if (verbose) cerr << "Before 2:" << endl << j2 << endl;
          assign_one_junction_read_counts(settings, summary, j2, require_overlap);
          j2["read_count_offset"] = mut["duplication_size"];
          if (verbose) cerr << "After 2:" << endl << j2 << endl;
          
          // @JEB 01-04-13
          // Calculate a frequency for the mobile element insertion from the reads supporting the new and old junctions on each side
          
          if (settings.polymorphism_prediction) {
            
            
            double a1 = from_string<uint32_t>(j1[SIDE_1_READ_COUNT]);
            double b1 = from_string<uint32_t>(j1[SIDE_2_READ_COUNT]);
            double c1 = from_string<uint32_t>(j1[NEW_JUNCTION_READ_COUNT]);
            double d1 = 2.0;
            
            if (j1[SIDE_1_READ_COUNT] == "NA") {
              a1 = 0; //"NA" in read count sets value to 1 not 0
              d1--;
            }
            if (j1[SIDE_2_READ_COUNT] == "NA") {
              b1 = 0; //"NA" in read count sets value to 1 not 0
              d1--;
            }
            
            double a2 = from_string<uint32_t>(j2[SIDE_1_READ_COUNT]);
            double b2 = from_string<uint32_t>(j2[SIDE_2_READ_COUNT]);
            double c2 = from_string<uint32_t>(j2[NEW_JUNCTION_READ_COUNT]);
            double d2 = 2.0;
            
            if (j2[SIDE_1_READ_COUNT] == "NA") {
              a2 = 0; //"NA" in read count sets value to 1 not 0
              d2--;
            }
            if (j2[SIDE_2_READ_COUNT] == "NA") {
              b2 = 0; //"NA" in read count sets value to 1 not 0
              d2--;
            }
            
            double frequency;
            if (d1 && d2) {
              frequency = (c1 + c2) / (c1 + (a1 + b1)/d1 + c2 + (a2 + b2)/d2);
              mut[FREQUENCY] = formatted_double(frequency, settings.polymorphism_precision_places, true).to_string();
            } else if (d1) {
              frequency = (c1) / (c1 + (a1 + b1)/d1);
              mut[FREQUENCY] = formatted_double(frequency, settings.polymorphism_precision_places, true).to_string();
            } else if (d2) {
              frequency = (c2) / (c2 + (a2 + b2)/d2);
              mut[FREQUENCY] = formatted_double(frequency, settings.polymorphism_precision_places, true).to_string();
            } else {
              // Can't calculate a frequency if no sides of the junction fall in unique sequence
              mut[FREQUENCY] = "NA";          
            }
            
            // don't record mutations below the cutoff frequency
            ASSERT(mut.entry_exists(FREQUENCY), "FREQUENCY field does not exist for mutatation:\n" + mut.as_string());
            
            if (mut[FREQUENCY] != "NA") {
              if (frequency < settings.polymorphism_frequency_cutoff) {
                mut.add_reject_reason("POLYMORPHISM_FREQUENCY_CUTOFF");
                // @JEB 08-08-13 we might want to keep the mutation as rejected. This discards completely.
                break;
              }
              if (1-frequency < settings.polymorphism_frequency_cutoff) {
                mut[FREQUENCY] = "1";
              }
            }
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
        
				gd.add(mut);
				break; // next JC1
			}
      
      // because this is a list, we only advance if we didn't delete the current entry
      if (!jc1_erased) jc1_it++;
		}

  }
  
  void MutationPredictor::predictJCtoINSorSUBorDEL(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc)
  {
    (void)summary;
    (void)mc;
    bool verbose = false;
    
    // variables pulled from the settings
    int32_t avg_read_length = static_cast<int32_t>(summary.sequence_conversion.avg_read_length);
    
    // @JEB 03-09-2014 Changed this section to produce INS and DEL instead of AMP,
    //                 when the mutation size is less than or equal to threshold in settings. 
		///
    for(diff_entry_list_t::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) //JC
		{
			cDiffEntry& j = **jc_it;
      
      //cout << j["side_1_seq_id"] << " " << j["side_2_seq_id"] << endl;
      //cout << j["side_1_strand"] << " " << j["side_2_strand"] << endl;
      
      // must be on same sequence
      if ((j["side_1_seq_id"] != j["side_2_seq_id"]))
        continue;
      
      // must be in same orientation (implies strands are opposite)
      int32_t side_1_strand = n(j["side_1_strand"]);
      int32_t side_2_strand = n(j["side_2_strand"]);
      
      // to be safe about making predictions, neither can be ambiguous 
      if ( (n(j["side_1_redundant"]) == 1) || (n(j["side_2_redundant"]) == 1) )
        continue;
      
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
				j["circular_chromosome"] = "1";
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
      
      // Insertion or deletion must be smaller than read length to be predicted
      // by this evidence alone.
      if (abs(size) > avg_read_length) 
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
      
      // If we are in polymorphism mode, propagate the frequency from JC evidence to mutation
      if (settings.polymorphism_prediction) 
        mut[FREQUENCY] = j[FREQUENCY];
      
      // Finally, add it
      gd.add(mut);
		}

  }
  
  void MutationPredictor::predictRAtoSNPorDELorINSorSUB(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc)
  {
    (void)summary;
    (void)mc;
    (void)jc;
    bool verbose = false;
    
    // Pull settings variables    
    vector<gd_entry_type> ra_types = make_vector<gd_entry_type>(RA);
    diff_entry_list_t ra = gd.list(ra_types);
    
		///
		// Ignore RA that overlap DEL or MC unless they are user
		// They are due to low spurious coverage in deleted regions!
		///
    
    vector<gd_entry_type> del_types = make_vector<gd_entry_type>(DEL);
    diff_entry_list_t del = gd.list(del_types);
    
    // Don't add deleted flags if we are in targeted sequencing mode
    if (!settings.targeted_sequencing) {
      
      for(diff_entry_list_t::iterator ra_it = ra.begin(); ra_it != ra.end(); ra_it++) //RA
      {
        cDiffEntry& ra_item = **ra_it;
        
        bool next_ra = false;
        
        // fon't remove these, no matter what!
        if (ra_item.entry_exists("user_defined"))
          continue;
        
        for(diff_entry_list_t::iterator del_it = del.begin(); del_it != del.end(); del_it++) //DEL
        {
          cDiffEntry& del_item = **del_it;
          
          if (ra_item["seq_id"] != del_item["seq_id"])
            continue;
          
          // there might be a problem here with insert_position > 0
          if ( (n(ra_item["position"]) >= n(del_item["position"])) && (n(ra_item["position"]) <= n(del_item["position"]) + n(del_item["size"]) - 1) )
          {
            ra_item["deleted"] = "1";
            next_ra = true;
            break;
          }
        }
        if (next_ra) continue;
        
        for(diff_entry_list_t::iterator mc_it = mc.begin(); mc_it != mc.end(); mc_it++) //MC
        {
          cDiffEntry& mc_item = **mc_it;
          
          if (ra_item["seq_id"] != mc_item["seq_id"]) continue;
          
          if ( (n(ra_item["position"]) >= n(mc_item["start"])) && (n(ra_item["position"]) <= n(mc_item["end"])) )
          {
            ra_item["deleted"] = "1";
            break;
          }
        }
      }
    }
    
    // Don't use rejected evidence
    ra.remove_if(cDiffEntry::rejected_and_not_user_defined());
    
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
      
      // Sometimes a SNP might be called in a deleted area because the end was wrong,
			// but it was corrected using a junction. (This catches this case.)
			if ( (!item.entry_exists("user_defined")) && (item.entry_exists("reject") || item.entry_exists("deleted")) )
			  continue;
      
      // If we are predicting mixed bases and not polymorphisms, then don't create
      // mutations for mixed frequency predictions (leave them as unassigned RA evidence)
      // But do mark them as "mixed"
      if (settings.mixed_base_prediction && (item[FREQUENCY] != "1")) {
        item["mixed"] = "1";
        continue;
      }
      
			bool same = false;
			if (!first_time)
			{
				if ( ((mut["end"] == item["position"]) && (n(mut["insert_end"]) + 1 == n(item["insert_position"])))
						|| ((n(mut["end"]) + 1 == n(item["position"])) && (item["insert_position"] == "0")) )
					same = true;
        
        // This code is only safe if every mutation has a frequency
        if (settings.polymorphism_prediction) {
          if ( (item[FREQUENCY] != "1") || (mut[FREQUENCY] != "1") //don't join polymorphisms
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
        ("ref_seq", (item["ref_base"] != ".") ? item["ref_base"] : "")
        ("new_seq", (item["new_base"] != ".") ? item["new_base"] : "")
				;
        
        if (settings.polymorphism_prediction) {
          new_mut[FREQUENCY] = item.entry_exists(FREQUENCY) ? item[FREQUENCY] : "1";
        }
				mut = new_mut;
			}
			else
			{
				mut
        ("insert_end", item["insert_position"])
        ("end", item["position"])
				;
				if (item["ref_base"] != ".") mut["ref_seq"] += item["ref_base"];
				if (item["new_base"] != ".") mut["new_seq"] += item["new_base"];
				mut._evidence.push_back(item._id);
			}
		}
		//don't forget the last one
		if (!first_time) muts.push_back(mut);
    
		///
		// Finally, convert these items into the fields needed for the various types of mutations
		///
    
		for (uint32_t i = 0; i < muts.size(); i++) {
      
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
  
  // Normalizes INS/DEL mutations by shifting them to the highest coordinates possible
  // by repeat units.
  
  // Adds additional fields "repeat_length", "repeat_ref_copies", and "repeat_new_copies"
  // to any repeat that is above a certain threshold length in the original sequence (by default 6 bp).
  
  // There are some difficulties dealing with deletions that have been split into two
  // e.g. deletion of one G and then another G from a larger G repeat, for now we skip those!
  // Insertions don't have the same problem. We could correct the final repeat number and
  // so on, but it's not obvious how to do this: set it to the final state for all mutations)
  // or set it to the state after each mutation and make that the original state for the next one.
  
  void MutationPredictor::normalize_and_annotate_tandem_repeat_mutations(Settings& settings, Summary& summary, cGenomeDiff& gd)
  { 
    (void)settings;
    (void) summary;

    // Sort for reproducibility -- but note this sort will possibly be changed by the
    // position shifts that we make within this function!!
    gd.sort();
    
    uint32_t minimum_tandem_repeat_length = 5;
    uint32_t minimum_repeat_ref_copies = 2;
    uint32_t maximum_repeat_sequence_length_to_show = 40;

    
    // moving previous mutations back if they are put on top of existing;
    // Currently only implemented for DEL
    cDiffEntry* previous_mutation = NULL;
    int32_t previous_position;
    
    diff_entry_list_t test_muts = gd.mutation_list();
    
    // Add additional fields for INS or DEL mutations that are in tandem repeats
    for(diff_entry_list_t::iterator it = test_muts.begin(); it != test_muts.end(); it++) {
      cDiffEntry& mut = **it;
      if (mut._type == INS) {
        
        // This field only exists when manually added to split an insertion into separate events 
        // or in polymorphism mode when a mutation may be split into each inserted base.
        // This code may inappropriately shift the position in those cases.
        if (mut.entry_exists(INSERT_POSITION) || settings.polymorphism_prediction) continue;
        
        int32_t size = mut["new_seq"].size();
        int32_t position = from_string<int32_t>(mut["position"]);
        string mutation_sequence = mut["new_seq"];
        uint32_t repeat_unit_size;
        string repeat_unit_sequence;
        find_repeat_unit(mutation_sequence, repeat_unit_size, repeat_unit_sequence);
                
        
        normalizeINSposition(ref_seq_info[mut["seq_id"]], mut, repeat_unit_sequence);
        int32_t new_position = from_string<int32_t>(mut["position"]);
        string new_mutation_sequence = mut["new_seq"];
        
        // ---> stub for future development of checking versus previous mutation
        
        // repeat info may have changed, so reload
        position = new_position;
        mutation_sequence = new_mutation_sequence;
        find_repeat_unit(mutation_sequence, repeat_unit_size, repeat_unit_sequence);
        uint32_t num_repeat_units = size / repeat_unit_size;
        
        // Note shift to +1 to get to where the first unit of a repeat would be for INS
        uint32_t original_num_repeat_units = find_original_num_repeat_units(ref_seq_info[mut["seq_id"]], position+1, repeat_unit_sequence);

        // save normalized position even if we aren't a repeat
        mut["position"] = to_string<int32_t>(position);
        mut["new_seq"] = mutation_sequence;
        
        if (original_num_repeat_units * repeat_unit_size < minimum_tandem_repeat_length)
          continue;
        
        if (original_num_repeat_units + num_repeat_units < minimum_repeat_ref_copies)
          continue;
        
        mut.erase("repeat_seq");
        if (repeat_unit_size <= maximum_repeat_sequence_length_to_show)
          mut["repeat_seq"] = repeat_unit_sequence;
        mut["repeat_length"] = s(repeat_unit_size);
        mut["repeat_ref_copies"] = s(original_num_repeat_units);
        mut["repeat_new_copies"] = s(original_num_repeat_units + num_repeat_units);
      }
      
      else if (mut._type == DEL) {
        
        // If we are in polymorphism mode, don't shift the coordinates
        // because we can shift two adjacent (but separately) deleted bases
        // onto the same position!
        if (settings.polymorphism_prediction) continue;
        
        // If we are annotated as 'within', do not shift the coordinates,
        // because the reference sequence is actually different by the time we are applied
        
        // We are still potentially in danger of doing the wrong thing here,
        // because a mutation could be applied only after one with the 'before' tag, making the shift
        // incorrect. So, setting 'no_normalize' is an out that can be used.
        
        if (mut.entry_exists("within") || mut.entry_exists("no_normalize") ) continue;
        
        int32_t size = from_string<int32_t>(mut["size"]);
        
        int32_t position = from_string<int32_t>(mut["position"]);
        string mutation_sequence = ref_seq_info.get_sequence_1(mut["seq_id"], position, position + size - 1);
        
        uint32_t repeat_unit_size;
        string repeat_unit_sequence;
        find_repeat_unit(mutation_sequence, repeat_unit_size, repeat_unit_sequence);
        
        normalizeDELposition(ref_seq_info[mut["seq_id"]], mut, repeat_unit_sequence);
        
        // Normalize may actually change the sequence used for the repeat... so call again here.
        int32_t new_position = from_string<int32_t>(mut["position"]);
        string new_mutation_sequence = ref_seq_info.get_sequence_1(mut["seq_id"], position, position + size - 1);
        
        // check against previous mutations to see if we are within the same repeat
        if (previous_mutation && (previous_mutation->_type == DEL)) {

          // Put things back if we are moving this to the exact same base
          int32_t updated_prev_position = from_string<int32_t>((*previous_mutation)["position"]);
          if (new_position == updated_prev_position) {
            // revert the previous mutation
            (*previous_mutation)["position"] = to_string<int32_t>(previous_position);
            previous_position = position;
            previous_mutation = &mut;
            continue;
          }
        }
        
        // must set before jumping out or updating
        previous_position = position;
        previous_mutation = &mut;
        
        // repeat info may have changed, so reload
        position = new_position;
        mutation_sequence = new_mutation_sequence;

        find_repeat_unit(mutation_sequence, repeat_unit_size, repeat_unit_sequence);
        uint32_t num_repeat_units = size / repeat_unit_size;

        // must set before jumping out
        previous_mutation = &mut;
        
        // Note that we add the number of repeats that were deleted
        uint32_t original_num_repeat_units = find_original_num_repeat_units(ref_seq_info[mut["seq_id"]], position, repeat_unit_sequence) + num_repeat_units;
        
        if (original_num_repeat_units * repeat_unit_size < minimum_tandem_repeat_length)
          continue;
        
        if (original_num_repeat_units < minimum_repeat_ref_copies)
          continue;
        
        // We only count if there is at least one unit remaining...
        if (original_num_repeat_units - num_repeat_units == 0)
          continue;
        
        mut.erase("repeat_seq");
        if (repeat_unit_size <= maximum_repeat_sequence_length_to_show)
          mut["repeat_seq"] = repeat_unit_sequence;
        mut["repeat_length"] = s(repeat_unit_size);
        mut["repeat_ref_copies"] = s(original_num_repeat_units);
        mut["repeat_new_copies"] = s(original_num_repeat_units - num_repeat_units);
      }
      else { // not deletion or insertion
        previous_mutation = NULL;
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
      int32_t pos = from_string<uint32_t>(mut[POSITION]) - unit_size;
      // note shift back up of position by one is correct
      // because INS are after this position, but AMP start at this position
      mut["position"] = to_string<int32_t>(pos+1);
      mut["new_copy_number"] = mut["repeat_new_copies"];
      mut["size"] = mut["repeat_length"];
      
      // delete all the repeat information...
      mut.erase(NEW_SEQ);
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
		(void)settings; //TODO; unused?
    bool verbose = false; // for debugging

		///
		//  Preprocessing of JC evidence
    //  NB: This call is likely redundant in the normal pipeline, but preserved here so that 
    //  predict can be a stand-alone call that is not dependent on other processing.
		///
    
    prepare_junctions(settings, summary, gd);

    ///
    // Create master lists of evidence that will be culled as they are used
    // by functions below.
    ///
    
    // Do not use rejected junctions unless they are user-defined
    vector<gd_entry_type> jc_types = make_vector<gd_entry_type>(JC);
		diff_entry_list_t jc = gd.list(jc_types);
    jc.remove_if(cDiffEntry::rejected_and_not_user_defined());
    
    // Do not use rejected missing coverage evidence
    vector<gd_entry_type> mc_types = make_vector<gd_entry_type>(MC);
		diff_entry_list_t mc = gd.list(mc_types);
    mc.remove_if(cDiffEntry::field_exists("reject"));
    
		///
		// evidence MC + JC => DEL mutation
		///
    predictMCplusJCtoDEL(settings, summary, gd, jc, mc);
    
		///
		// evidence JC + JC = MOB mutation
		///
    
    predictJCplusJCtoMOB(settings, summary, gd, jc, mc);
		
		///
		// evidence JC => INS, SUB, DEL mutations
    ///
    
    predictJCtoINSorSUBorDEL(settings, summary, gd, jc, mc);
    
		///
		// evidence RA => SNP, DEL, INS, SUB
		///
    
    predictRAtoSNPorDELorINSorSUB(settings, summary, gd, jc, mc);
    
		
    ///
		// Read Alignments => SNP, DEL, INS, SUB
		///
    
    // Normalizing in polymorphism mode causes consistency issues when multiple INS or DEL are next to each other
    // (which is not possible in consensus mode).
    if (!settings.polymorphism_prediction) {
      normalize_and_annotate_tandem_repeat_mutations(settings, summary, gd);
    }
        
    ///
		// mutation INS => mutation AMP
		///
    normalize_INS_to_AMP(settings, summary, gd);
      
    ///////////////////////////////////////////////////////
    // Check to be sure the "frequency" field is present //
    // as appropriate in consensus/polymorphism mode.    //
    ///////////////////////////////////////////////////////
    
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
    
    uint32_t num_repeat_units = 0;
    
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
  
  void MutationPredictor::normalizeDELposition(cAnnotatedSequence& ref_seq, cDiffEntry& de, string& repeat_sequence)
  {
    // Don't shift mediated or between mutations
    if (de.entry_exists("mediated") || de.entry_exists("between")) return;
    
    int32_t test_position = from_string<int32_t>(de[POSITION]) + from_string<int32_t>(de[SIZE]);   
    // The offset of 1 ensures we are on the first base of a possible repeat unit (relative to the reference).
    
    // attempt to move forward by repeat units at a time from the current position
    while ( test_position + repeat_sequence.size() - 1 <= ref_seq.get_sequence_length()) {
      
      string test_sequence = ref_seq.get_sequence_1(test_position, test_position + repeat_sequence.size() - 1);
      
      if (test_sequence != repeat_sequence) break;
      test_position += repeat_sequence.size();
    }
    
    string new_position = to_string<int32_t>(test_position - from_string<int32_t>(de[SIZE]));
    
    if (new_position != de[POSITION]) {
      de["_original_aligned_position"] = de[POSITION]; // Save the original position for marking in alignments
      de[POSITION] = new_position;
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
  //          *
  // ABCABCABCABC
  //
  // For INS this is the position of the base before the insertion (same as the mutation position).
  // For DEL it is one position before the first deleted base (one before the mutation position).
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
      
      string test_sequence = ref_seq.get_sequence_1(test_position, test_position + repeat_sequence.size() - 1);
      
      if (test_sequence != repeat_sequence) break;
      
      num_repeat_units++;
    }
    
    return num_repeat_units;
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
      for(size_t i=0; i<seq.get_sequence_length(); ++i) {
        char base = seq.m_fasta_sequence.m_sequence[i];
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
        if (verbose) cout << f.SafeGet("name") << " " << f.get_start_1() << " " << f.get_end_1() << " " << f.get_strand() << endl;
        
        
        // By taking anything within a gene as NONCODING, it will also include **introns** in this category
        // instead of an alternative strategy which might be to give them their own category or call them intergenic.
        
        // Catches tRNA/rRNA/pseudogenes...
        if (f["type"] == "gene") {
          vector<cLocation> sub_locations = f.m_location.get_all_sub_locations();
          for(vector<cLocation>::iterator it3=sub_locations.begin(); it3!=sub_locations.end(); ++it3) {
            cLocation& loc = *it3;
            for (int32_t this_location_1=loc.get_start_1(); this_location_1<=loc.get_end_1(); this_location_1++) {
              int32_t this_location_0 = this_location_1-1;
              for (size_t b=0; b<BaseSubstitutionEffects::base_char_list.size(); b++)
                seq_bse[this_location_0*4+b] = max(seq_bse[this_location_0*4+b], f.m_pseudo ? pseudogene_base_substitution : noncoding_base_substitution);
            }
          }
        }
        
        // Remainder is only for coding sequences (and not pseudogenes)
        if ((f["type"] != "CDS") || f.m_pseudo)
          continue;
        
        // initialize gene structure
        cGeneFeature g(f);
        
        total_orfs++;
        
        // Piece together the gene - at each nucleotide make all three possible changes
        
        string this_codon = "NNN";
        vector<cLocation> sub_locations = g.m_location.get_all_sub_locations();
        size_t on_codon_pos_0 = 0; // The position within a codon... indexed to start at 0.
        vector<uint32_t> this_codon_locations_0(3, 0); // 0-indexed
        
        uint32_t total_nucleotide_length = 0; 
        int8_t strand = sub_locations.front().m_strand;
        
        for(vector<cLocation>::iterator it3=sub_locations.begin(); it3!=sub_locations.end(); ++it3) {
          cLocation& loc = *it3;
          total_nucleotide_length += loc.m_end - loc.m_start + 1;
        }
        
        // >>CODE_INDETERMINATE
        // If we have an indeterminate start or end, we may need to start out-of-frame
        // @JEB: this happens at the ends of sequences (usually contigs from a de novo assembly,
        //       when a gene runs into the end of a contig
        // We also need to adjust things so that codon #1 is either missing (indeterminate gene start)
        // or happens at the last codon present (indeterminate gene end)
        
        uint32_t total_amino_acid_length = total_nucleotide_length / 3;
        uint32_t on_codon_number_1 = 1; // The number of the amino acid / codon indexed to start at 1
        if (strand == -1) on_codon_number_1 = total_amino_acid_length;

        
        for(vector<cLocation>::iterator it3=sub_locations.begin(); it3!=sub_locations.end(); ++it3) {
          
          cLocation& loc = *it3;
          ASSERT(strand == loc.m_strand, "CDS has sublocations on different strands: " + g["name"]);
          
          for (int32_t pos_1=loc.m_start; pos_1<=loc.m_end; pos_1++) {
            int32_t pos_0 = pos_1 - 1;
            
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
            else
            {
              seq_bcs[pos_0] = (loc.get_strand() == +1 ? forward : reverse);
            }
            
            //// Handle codon synonymous/nonsynonymous changes
            this_codon_locations_0[on_codon_pos_0] = pos_0;
            this_codon[on_codon_pos_0] = seq.get_sequence_1(pos_1);
            
            on_codon_pos_0++;
            
            // The codon is filled, now make all mutations and assign to proper nucleotides
            if (on_codon_pos_0 == 3) {
              
              // change the codon to the right strand
              string original_codon = this_codon;
              if (strand == -1)
                original_codon = reverse_complement(original_codon);
              
              char original_amino_acid = cReferenceSequences::translate_codon(original_codon, g.translation_table, on_codon_number_1);
              
               
              for (int32_t test_codon_index=0; test_codon_index<3; test_codon_index++) {
                
                for (size_t b=0; b<BaseSubstitutionEffects::base_char_list.size(); b++) {
                  
                  char mut_base = BaseSubstitutionEffects::base_char_list[b];
                  
                  string test_codon = this_codon;
                  test_codon[test_codon_index] = mut_base;
                  
                  if (strand == -1)
                    test_codon = reverse_complement(test_codon);
                  
                  char mut_amino_acid = cReferenceSequences::translate_codon(test_codon, g.translation_table, on_codon_number_1);
                  
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
              
              on_codon_number_1 += strand;
              on_codon_pos_0 = 0;
            }
          }
        } // end sublocation loop
        
        // @JEB 01-24-2014 - will there be crashes with this? 
        // BUG REPORT: It fails when a gene is truncated at the end of a contig
        // TEST CASE: 
        if (on_codon_pos_0 != 0) {
          WARN("Number of base pairs in CDS not a multiple of 3: " + g["name"]);
        }
        // Currently these changes are not counted as anything when the codon overlaps the end of the sequence
        // that seems fine!
        
        /// code to debug consistency with Perl
        /*
        cout << ">" << g.name << endl;
        cout << codon_sequence << endl;
        cout << amino_acid_sequence << endl;
        */
        
      } // end feature loop
      
      if (verbose) {
        for(size_t i_0=0; i_0<seq.get_sequence_length(); ++i_0) {
          char base = seq.m_fasta_sequence.m_sequence[i_0];        
          cout << (i_0+1) << "\t" << base << "\t" << seq_bcs[i_0] << "\t" 
          << seq_bse[i_0*4+0] << "\t" << seq_bse[i_0*4+1] << "\t" << seq_bse[i_0*4+2] << "\t" << seq_bse[i_0*4+3] << endl;
        }
      }
      
      // For debugging consistency with Perl
      /*
      for(size_t pos_1=1; pos_1 <= seq.get_sequence_length(); pos_1++) {
        cout << pos_1 << " " << seq.get_sequence_1(pos_1);
        size_t pos_0 = pos_1-1;
        
        BaseSubstitutionEffectPositionInfo position_info = position_info_1(ref_seq_info, seq.m_seq_id, pos_1);
        BaseType bt = position_info.m_base_type;
        
        cout << " " << bt;
        
        //for (uint8_t base_index=0; base_index<4; ++base_index) {
        //  cout << " " << ((seq_bse[pos_0*4+base_index] == nonsynonymous_base_substitution) ? "1" : "0");
        //}
        
        cout << endl;
      }
      */
      
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
                         bool verbose
                         )
  {
    // Could be a parameter > this is a "large" mutation, <= this is an "small" mutation
    int32_t large_size_cutoff = 50;
      
    // Figure out the names of all "repeat" columns
    map<string,bool> mob_name_hash;
    map<string,bool> con_name_hash;
    
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
            mob_name_hash[mut["mediated"]] = true;
        }
        if (mut._type == DEL) {
          if (mut.entry_exists("mediated"))
            mob_name_hash[mut["mediated"]] = true;
        }
      }
    }
    vector<string> mob_name_list = map_keys_to_list<string,bool>(mob_name_hash);
    sort(mob_name_list.begin(), mob_name_list.end());
    vector<string> con_name_list = map_keys_to_list<string,bool>(con_name_hash);
    sort(con_name_list.begin(), con_name_list.end());
    
    vector<string> column_headers;
    column_headers.push_back("file");
    column_headers.push_back("sample");
    column_headers.push_back("treatment");
    column_headers.push_back("population");
    column_headers.push_back("time");
    column_headers.push_back("clone");
    column_headers.push_back("total");
    column_headers.push_back("base_substitution");
    column_headers.push_back("small_indel");
    column_headers.push_back("large_deletion");
    column_headers.push_back("large_insertion");
    column_headers.push_back("large_amplification");
    column_headers.push_back("large_substitution");
    column_headers.push_back("mobile_element_insertion");
    column_headers.push_back("gene_conversion");
    column_headers.push_back("deleted_bp");
    column_headers.push_back("inserted_bp");
    column_headers.push_back("repeat_inserted_bp");
    column_headers.push_back("called_bp");
    column_headers.push_back("total_bp");
    
    vector<string> header_snp_types = prefix_each_in_vector(snp_types, "base_substitution.");
    column_headers.insert(column_headers.end(),header_snp_types.begin(), header_snp_types.end());
    
    vector<string> header_mob_name_list = prefix_each_in_vector(mob_name_list, "mobile_element.");
    column_headers.insert(column_headers.end(),header_mob_name_list.begin(), header_mob_name_list.end());
    
    vector<string> header_con_name_list = prefix_each_in_vector(con_name_list, "gene_conversion.");
    column_headers.insert(column_headers.end(),header_con_name_list.begin(), header_con_name_list.end());
    
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

    
    output_file << join(column_headers, ",") << endl;
    
    for (vector<cGenomeDiff>::iterator it=sorted_genome_diffs.begin(); it != sorted_genome_diffs.end(); ++it) {
      cGenomeDiff &gd = *it;
      //uout("Counting mutations " + gd.metadata.run_name);
      
      vector<string> line_prefix_items = make_vector<string>(gd.get_file_name())
          (gd.metadata.title)(gd.metadata.treatment)(gd.metadata.population)
          ((gd.metadata.time != -1.0) ? to_string<double>(gd.metadata.time) : "")
          (gd.metadata.clone);
      
      string detailed_line_prefix = join(line_prefix_items, "\t");
      
      column_headers.push_back("file");
      column_headers.push_back("sample");
      column_headers.push_back("treatment");
      column_headers.push_back("population");
      column_headers.push_back("time");
      column_headers.push_back("clone");
      
      BaseSubstitutionEffectCounts this_bsec;
      // deep copy totals of entire sequence
      if (base_substitution_statistics)
        this_bsec = total_bsec;
      
      // Zero out counts
      int32_t total_deleted = 0;
      int32_t total_inserted = 0;
      int32_t total_repeat_inserted = 0;

      int32_t total_bp = ref_seq_info.total_length();
      
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
      
      
      // BEGIN for each mutation
      uint32_t total_mutations(0);
      diff_entry_list_t mut_list = gd.mutation_list();
      for (diff_entry_list_t::iterator it=mut_list.begin(); it != mut_list.end(); ++it) {		
        
        cDiffEntry& mut = **it;
        
        // Don't count mutations that were hidden by later deletions, but kept for phylogenetic inference.
        // @JEB: actually we do want to count these!
        //if (mut.is_marked_deleted()) {
        //  if (verbose) cerr << "Skipping deleted: " << mut << endl;
        //  continue;
        //}
        
        // Don't count polymorphisms - could make it an option to partially count them
        if (!count_polymorphisms && mut.entry_exists(FREQUENCY) && from_string<double>(mut[FREQUENCY]) != 1.0) {
          if (verbose) cerr << "Skipping polymorphic: " << mut << endl;
          continue;
        }
          
        if (verbose) cerr << "Counting: " << mut << endl;
        total_mutations++;
        
        if (mut._type == SNP) {
          count["base_substitution"][""]++;
          if (base_substitution_statistics) {
            this_bsec.change_position_1_observed_totals(ref_seq_info, bse, mut[SEQ_ID], from_string<int32_t>(mut[POSITION]), mut[NEW_SEQ], +1);					
          }
          count["type"][mut["snp_type"]]++;
          base_substitution_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
        }
        
        if (mut._type == DEL) {
          total_deleted += from_string<int32_t>(mut[SIZE]);
          
          if (mut.entry_exists("mediated"))
            count["mob"][mut["mediated"]]++;
          
          if (from_string<int32_t>(mut[SIZE]) > large_size_cutoff) {
            count["large_deletion"][""]++;
            large_deletion_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          } else {
            count["small_indel"][""]++;
            small_indel_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          }
        }
        
        if (mut._type == INS) {
          int32_t ins_size = mut[NEW_SEQ].size();
          
          total_inserted += ins_size;
          
          if (mut.entry_exists("mediated"))
            count["mob"][mut["mediated"]]++;
          
          if (ins_size > large_size_cutoff) {
            count["large_insertion"][""]++;
            large_insertion_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          } else {
            count["small_indel"][""]++;
            small_indel_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          }
        }
        
        if (mut._type == SUB) {
          int32_t old_size = from_string<int32_t>(mut[SIZE]);
          int32_t new_size = mut[NEW_SEQ].size();
          
          if (new_size - old_size > 0)
            total_inserted += new_size - old_size;
          else
            total_deleted += old_size - new_size;
          
          if (mut.entry_exists("mediated"))
            count["mob"][mut["mediated"]]++;
          
          if (abs(new_size - old_size) > large_size_cutoff) {
            count["large_substitution"][""]++;
            large_substitution_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          } else {
            count["small_indel"][""]++;
            small_indel_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          }
        }
        
        if (mut._type == CON) {
          // TODO: Need size change?
          
          if (mut.entry_exists("mediated"))
            count["con"][mut["mediated"]]++;
          count["gene_conversion"][""]++;
          gene_conversion_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
        }
        
        if (mut._type == MOB) {
          int32_t rpos = -1;
          // This includes the MOB and any ins/del start/end adjustments, NOT duplication size
          string mob_region;
          string repeat_seq = ref_seq_info.repeat_family_sequence(mut["repeat_name"], 1, mut.entry_exists("mob_region") ? &mut["mob_region"] : NULL);
          
          int32_t this_length = repeat_seq.size();
          this_length += from_string<int32_t>(mut["duplication_size"]);
          
          total_inserted += this_length;
          total_repeat_inserted += this_length;
          //print "Repeat $mut->{repeat_name} $this_length bp\n";
          
          count["mob"][mut["repeat_name"]]++;
          count["mobile_element_insertion"][""]++;
          mobile_element_insertion_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          
        } 
        
        if (mut._type == AMP) {
          int32_t this_size = mut.mutation_size_change(ref_seq_info);
          total_inserted += this_size;
          
          if (this_size > large_size_cutoff) {
            count["large_amplification"][""]++;
            large_amplification_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          } else {
            count["small_indel"][""]++;
            small_indel_lines.push_back(detailed_line_prefix + "\t" + mut.as_string());
          }
        }
      }
      
      // statistics for UN
      int32_t un_bp = 0;
      
      diff_entry_list_t un_list = gd.list(make_vector<gd_entry_type>(UN));
      for (diff_entry_list_t::iterator it=un_list.begin(); it!= un_list.end(); ++it) {
        cDiffEntry& un = **it;
        un_bp += from_string<int32_t>(un[END]) - from_string<int32_t>(un[START]) + 1;
        
        // subtract these positions from the possible observations of base pair statistics
        if (base_substitution_statistics)
          for (uint32_t pos_1 = from_string<uint32_t>(un[START]); pos_1 <= from_string<uint32_t>(un[END]); pos_1++)
            this_bsec.change_position_1_possible_totals(ref_seq_info, bse, un[SEQ_ID], pos_1, -1);					
        
      }
      // END for each mutation

      
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
      this_columns.push_back(to_string(total_deleted));
      this_columns.push_back(to_string(total_inserted));
      this_columns.push_back(to_string(total_repeat_inserted));
      this_columns.push_back(to_string(called_bp));
      this_columns.push_back(to_string(total_bp));
      
      vector<string> snp_type_counts = map_key_list_to_values_as_strings(count["type"], snp_types);
      this_columns.insert(this_columns.end(),snp_type_counts.begin(), snp_type_counts.end());
      
      vector<string> mob_type_counts = map_key_list_to_values_as_strings(count["mob"], mob_name_list);
      this_columns.insert(this_columns.end(),mob_type_counts.begin(), mob_type_counts.end());
      
      vector<string> con_type_counts = map_key_list_to_values_as_strings(count["con"], con_name_list);
      this_columns.insert(this_columns.end(),con_type_counts.begin(), con_type_counts.end());
      
      
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
      
      detailed_output_file << "BASE SUBSTITUTIONS:" << base_substitution_lines.size() << "\n";
      detailed_output_file << join(base_substitution_lines, "\n") << "\n\n";
      
      detailed_output_file << "SMALL INDELS:" << small_indel_lines.size() << "\n";
      detailed_output_file << join(small_indel_lines, "\n")<< "\n\n";
      
      detailed_output_file << "LARGE DELETIONS:" << large_deletion_lines.size() << "\n";
      detailed_output_file << join(large_deletion_lines, "\n")<< "\n\n";
      
      detailed_output_file << "LARGE INSERTIONS:" << large_insertion_lines.size() << "\n";
      detailed_output_file << join(large_insertion_lines, "\n")<< "\n\n";
      
      detailed_output_file << "LARGE AMPLIFICATIONS:" << large_amplification_lines.size() << "\n";
      detailed_output_file << join(large_amplification_lines, "\n")<< "\n\n";
      
      detailed_output_file << "LARGE SUBSTITUTIONS:" << large_substitution_lines.size() << "\n";
      detailed_output_file << join(large_substitution_lines, "\n")<< "\n\n";
      
      detailed_output_file << "MOBILE ELEMENT INSERTIONS:" << mobile_element_insertion_lines.size() << "\n";
      detailed_output_file << join(mobile_element_insertion_lines, "\n") << "\n\n";
      
      detailed_output_file << "GENE CONVERSIONS:" << gene_conversion_lines.size() << "\n";
      detailed_output_file << join(gene_conversion_lines, "\n")<< "\n";
    }
    
  }
  

} // namespace breseq
