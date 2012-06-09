
/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#include "libbreseq/mutation_predictor.h"

#include "libbreseq/output.h"


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

    return (n((*a)["pos_hash_score"]) > n((*b)["pos_hash_score"]));
  }

	bool MutationPredictor::sort_by_reject_score(const counted_ptr<cDiffEntry>& a, const counted_ptr<cDiffEntry>& b)
	{
		uint32_t a_reject_order = a->number_reject_reasons();
		uint32_t b_reject_order = b->number_reject_reasons();

		// sort by seq_id, position, fewer reject reasons, then score (highest to lowest)
    if (a_reject_order != b_reject_order) 
      return (a_reject_order < b_reject_order);
            
    return (n((*a)["pos_hash_score"]) > n((*b)["pos_hash_score"]));
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
	 Title   : predict
	 Usage   : $mp->predict();
	 Function: Predicts mutations from evidence in a GenomeDiff and adds them to it
	 Returns :

	*/
	void MutationPredictor::predict(Settings& settings, Summary& summary, cGenomeDiff& gd)
	{
		(void)settings; //TODO; unused?
    bool verbose = false; // for debugging

    int32_t avg_read_length = summary.sequence_conversion.avg_read_length;
    int32_t max_read_length = summary.sequence_conversion.max_read_length;
    
		///
		//  Preprocessing of JC evidence
		///
    
    // For all that follows, we need information about repeat_regions overlapping the sides of junctions
    vector<gd_entry_type> jc_types = make_vector<gd_entry_type>(JC);
		diff_entry_list_t jc = gd.list(jc_types);
    // Don't count rejected ones at all, this can be relaxed, but it makes MOB 
    // prediction much more complicated and prone to errors.
    
    // This sections just normalizes read counts to the average coverage of the correct sequence fragment
    for (diff_entry_list_t::iterator jc_it=jc.begin(); jc_it!=jc.end(); jc_it++) {
      cDiffEntry& j = **jc_it;
      
      double side_1_correction = (summary.sequence_conversion.avg_read_length - 1 - abs(from_string<double>(j["alignment_overlap"])) - from_string<double>(j["continuation_left"])) / (summary.sequence_conversion.avg_read_length - 1);
      
      if (j[SIDE_1_READ_COUNT] == "NA")
        j[SIDE_1_COVERAGE] = "NA";
      else
        j[SIDE_1_COVERAGE] = to_string<double>(from_string<double>(j[SIDE_1_READ_COUNT]) / summary.unique_coverage[j[SIDE_1_SEQ_ID]].average / side_1_correction);
       
      double side_2_correction = (summary.sequence_conversion.avg_read_length - 1 - abs(from_string<double>(j["alignment_overlap"])) - from_string<double>(j["continuation_right"])) / (summary.sequence_conversion.avg_read_length - 1);

      if (j[SIDE_2_READ_COUNT] == "NA")
        j[SIDE_2_COVERAGE] = "NA";
      else
        j[SIDE_2_COVERAGE] = to_string<double>(from_string<double>(j[SIDE_2_READ_COUNT]) / summary.unique_coverage[j[SIDE_2_SEQ_ID]].average);

      //corrects for overlap making it less likely for a read to span
      double overlap_correction = (summary.sequence_conversion.avg_read_length - 1 - abs(from_string<double>(j["alignment_overlap"])) - from_string<double>(j["continuation_left"]) - from_string<double>(j["continuation_right"])) / (summary.sequence_conversion.avg_read_length - 1);
      double new_junction_average_read_count = (summary.unique_coverage[j[SIDE_1_SEQ_ID]].average + summary.unique_coverage[j[SIDE_2_SEQ_ID]].average) / 2;
      
      j[NEW_JUNCTION_COVERAGE] = to_string<double>(from_string<double>(j[NEW_JUNCTION_READ_COUNT]) / new_junction_average_read_count / overlap_correction);
    }
      
    const int32_t max_distance_to_repeat = 50;

		for (diff_entry_list_t::iterator jc_it=jc.begin(); jc_it!=jc.end(); jc_it++)
		{
			cDiffEntry& j = **jc_it;
      
			j["_side_1_read_side"] = "-1";
			j["_side_2_read_side"] = "1";

			for (uint32_t side = 1; side <= 2; side++)
			{
				string side_key = "side_" + s(side);

				cSequenceFeaturePtr is = ref_seq_info.find_closest_repeat_region_boundary(
					n(j[side_key + "_position"]),
					ref_seq_info[j[side_key + "_seq_id"]].m_repeats,
					max_distance_to_repeat,
					n(j[side_key + "_strand"])
				);
				if (is.get() != NULL)
				{
					j["_" + side_key + "_is"] = "1";
					j["_" + side_key + "_is_start"] = s(is->m_location.get_start_1());
					j["_" + side_key + "_is_end"] = s(is->m_location.get_end_1());
          j["_" + side_key + "_is_name"] = (*is)["name"];
          j["_" + side_key + "_is_strand"] = s(is->m_location.get_strand());
				}
				
        uint32_t test_val = n(j[side_key + "_redundant"]);
				j[side_key + "_annotate_key"] = (j.entry_exists("_" + side_key + "_is_start") || n(j[side_key + "_redundant"])) ? "repeat" : "gene";
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
    
    // Now, remove rejected from the list after annotating them.
    jc.remove_if(cDiffEntry::field_exists("reject"));

		///
		// evidence MC + JC => DEL mutation
		///
    
    vector<gd_entry_type> mc_types = make_vector<gd_entry_type>(MC);
		diff_entry_list_t mc = gd.list(mc_types);
    mc.remove_if(cDiffEntry::field_exists("reject"));

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
      
      // We will be erasing inside the jc_it loop.  This is to keep
      // track of whether or not we should iterate to the next element.
      bool jc_it_iterator = true;

			for(diff_entry_list_t::iterator jc_it = jc.begin(); jc_it != jc.end(); ) //JC
			{
				cDiffEntry& jc_item = **jc_it;
        jc_it_iterator = true;

				if (jc_item["side_1_seq_id"] != mut["seq_id"] || jc_item["side_2_seq_id"] != mut["seq_id"])  {
          // Iterate it ONLY if we haven't erased something
          if(jc_it_iterator)jc_it++;  
          continue;  
        }
					
        
				// We always know that the lower coordinate part of the junction is first, hence these
				// assumptions about the strands hold.
				if (
					   (n(jc_item["side_1_position"]) == n(mut["position"])-1)
					&& (n(jc_item["side_1_strand"]) == -1)
					&& (n(jc_item["side_2_position"]) == n(mut["position"])+n(mut["size"]))
					&& (n(jc_item["side_2_strand"]) == +1)
           )
				{
          
          //it's possible that one or both sides are in repeat elements
          cSequenceFeature* r1_pointer = within_repeat(jc_item["side_1_seq_id"], n(jc_item["side_1_position"]));
          cSequenceFeature* r2_pointer = within_repeat(jc_item["side_2_seq_id"], n(jc_item["side_2_position"]));
          
          // one repeat cases where the end matches up exactly
          if (r1_pointer) 
          {
            // must match up to an end of the repeat
            if ((n(jc_item["side_1_position"]) == static_cast<int32_t>(r1_pointer->m_location.get_start_1()))
             || (n(jc_item["side_1_position"]) == static_cast<int32_t>(r1_pointer->m_location.get_end_1())))
            {
              mut["mediated"] = (*r1_pointer)["name"];
            }
          }
          else if (r2_pointer) 
          {
            // must match up to an end of the repeat
            if ((n(jc_item["side_2_position"]) == static_cast<int32_t>(r2_pointer->m_location.get_start_1()))
                || (n(jc_item["side_2_position"]) == static_cast<int32_t>(r2_pointer->m_location.get_end_1())))
            {
              mut["mediated"] = (*r2_pointer)["name"];
            }
          }    
                      
					mut._evidence.push_back(jc_item._id);
					jc_it = jc.erase(jc_it); // iterator is now past the erased element
          jc_it_iterator = false; //We just removed the current jc, do not iterate.
					gd.add(mut);
          if (verbose)
            cout << "**** Junction precisely matching deletion boundary found ****\n";
				}
        
        // Iterate it ONLY if we haven't erased something
        if(jc_it_iterator)jc_it++;
			}

//      -> //need to have within repeat also return the distance to the end flushness
      cSequenceFeature* r1_pointer = within_repeat(mut["seq_id"], n(mut["position"]));
      cSequenceFeature* r2_pointer = within_repeat(mut["seq_id"], n(mut["position"]) + n(mut["size"]));
      
			///
			// (2) there is is no junction, but both ends of the deletion are in different copies of the same repeat sequence
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

				// remember the name of the element
				mut["between"] = r1["name"];
				gd.add(mut);
        
        if (verbose)
          cout << "**** Ends of junction in copies of same repeat element ****\n";
				continue; // to next mc_item
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

        
				// need to adjust the non-unique coords
				if (redundant_deletion_side == -1)
				{
//         -> // Need to get actual copy that it is within here, instead of other copies.
          int32_t discrepancy_dist = r.get_end_1() - n(j[j["_is_interval"] + "_position"]);
          
					uint32_t move_dist = r.get_end_1() + 1 - n(mut["position"]);
					mut["position"] = s(n(mut["position"]) + move_dist);
					mut["size"] = s(n(mut["size"]) - move_dist);
				}
				else
				{
					int32_t move_dist = (n(mut["position"]) + n(mut["size"]) - 1) - (r.get_start_1() - 1);
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
        
        if (verbose)
          cout << "**** Junction with repeat element corresponding to deletion boundaries found ****\n";

        break; // done looking at jc_items
			}
		}


		///
		// evidence JC + JC = MOB mutation
		///

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

		jc.sort(MutationPredictor::sort_by_hybrid);

    for(diff_entry_list_t::iterator jc1_it = jc.begin(); jc1_it != jc.end(); jc1_it++) //JC1
		{
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
				jc1_it = jc.erase(jc1_it); // iterator is now past element erased
        jc1_it--;                  // and must be moved back because loop will move forward
				jc.erase(it_delete_list_2[i]);

				// Create the mutation, with evidence
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

				// get any unique junction sequence
				JunctionInfo j1i(j1["key"]);
				string j1_unique_read_sequence = j1i.unique_read_sequence;

				JunctionInfo j2i(j2["key"]);
				string j2_unique_read_sequence = j2i.unique_read_sequence;

				// _gap_left and _gap_right also refer to the top strand of the genome

				mut["_ins_start"] = "";
				mut["_ins_end"] = "";

				mut["_del_start"] = "0";
				mut["_del_end"] = "0";

        uint32_t start_1 = 0, end_1 = 0, pos_1 = 0;

        // sometimes the ends of the IS are not quite flush
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

				if (n(mut["_gap_left"]) >= 0)
				{
					if (verbose)
						cout << "J1 NF:" << j1_not_flush_seq << " U:" << j1_unique_read_sequence << endl;

					if (n(j1["_" + j1["_is_interval"] + "_read_side"]) != n(j1[j1["_is_interval"] + "_strand"]))
					{
						j1_not_flush_seq = reverse_complement(j1_not_flush_seq);
					}

					if (n(j1["_" + j1["_is_interval"] + "_read_side"]) == -1)
					{
						mut["_gap_left"] = j1_not_flush_seq + j1_unique_read_sequence;
					}
					else
					{
						mut["_gap_left"] = j1_unique_read_sequence + j1_not_flush_seq;
					}
					mut["_ins_start"] = mut["_gap_left"];
				}
				else if (n(mut["_gap_left"]) < 0)
				{
					mut["_del_start"] = s(abs(n(mut["_gap_left"])));
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


				if (n(mut["_gap_right"]) >= 0)
				{
					if (verbose)
						cout << "J2 NF:" << j2_not_flush_seq << " U:" << j2_unique_read_sequence << endl;

					if ( n(j2["_" + j2["_is_interval"] + "_read_side"]) * n(j2[j2["_is_interval"] + "_strand"]) == -1)
					{
						j2_not_flush_seq = reverse_complement(j2_not_flush_seq);
					}

					if (n(j2["_" + j2["_is_interval"] + "_read_side"]) == -1)
					{
						mut["_gap_right"] = j2_not_flush_seq + j2_unique_read_sequence;
					}
					else
					{
						mut["_gap_right"] = j2_unique_read_sequence + j2_not_flush_seq;
					}
					mut["_ins_end"] = mut["_gap_right"];
				}
				else if (n(mut["_gap_right"]) < 0)
				{
					mut["_del_end"] = s(abs(n(mut["_gap_right"])));
				}

				// At this point any added junction sequences are on the strand as you would see them in the alignment.
				// we may need to reverse complement and change sides.

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

				int32_t j1_is_overlap_length = (n(j1["_" + j1["_is_interval"] + "_read_side"]) == -1) ? n(j1["max_left"]) : n(j1["max_right"]);
				int32_t j2_is_overlap_length = (n(j2["_" + j2["_is_interval"] + "_read_side"]) == -1) ? n(j2["max_left"]) : n(j2["max_right"]);

				if (verbose) {
					cout << "J1 IS overlap length: " << j1_is_overlap_length << endl;
					cout << "J2 IS overlap length: " << j2_is_overlap_length << endl;
				}

        /*
        pos_1 = n(j1[j1["_is_interval"] + "_position"]);
        string j1_is_seq_matched = "";
				if (n(j1[j1["_is_interval"] + "_strand"]) == -1)
        {
          start_1 = pos_1 - j1_is_overlap_length - 1;
          end_1   = pos_1 - j1_not_flush_length;
					j1_is_seq_matched = ref_seq_info.get_sequence_1 (
            j1[j1["_is_interval"] + "_seq_id"],
            start_1,
            end_1
					);
					j1_is_seq_matched = reverse_complement(j1_is_seq_matched);
				}
				else
				{          
          start_1 = pos_1 + j1_not_flush_length;
          end_1   = pos_1 + j1_is_overlap_length - 1;
          j1_is_seq_matched = ref_seq_info.get_sequence_1 (
            j1[j1["_is_interval"] + "_seq_id"],
            start_1,
            end_1
					);
				}
        */
        /*

        pos_1 = n(j2[j2["_is_interval"] + "_position"]);
        string j2_is_seq_matched = "";
				if (n(j2[j2["_is_interval"] + "_strand"]) == -1)
				{
          start_1 = pos_1 - j2_is_overlap_length - 1;
          end_1   = pos_1 - j2_not_flush_length;
          j2_is_seq_matched = ref_seq_info.get_sequence_1 (
            j2[j2["_is_interval"] + "_seq_id"],
            start_1,
            end_1
          );
					j2_is_seq_matched = reverse_complement(j2_is_seq_matched);
				}
				else
				{
          start_1 = pos_1 + j2_not_flush_length;
          end_1   = pos_1 + j2_is_overlap_length - 1;
          j2_is_seq_matched = ref_seq_info.get_sequence_1 (
            j2[j2["_is_interval"] + "_seq_id"],
            start_1,
            end_1
					);
				}
         */

        // what are the actual sequences of this length at the end of the IS elements?
        start_1 = n(j1["_" + j1["_is_interval"] + "_is_start"]);
        end_1   = start_1 + j1_is_overlap_length - 1;
				string j1_left_is_sequence = ref_seq_info.get_sequence_1 (
          j1[j1["_is_interval"] + "_seq_id"],
          start_1,
          end_1
				);


        end_1   = n(j1["_" + j1["_is_interval"] + "_is_end"]);
        start_1 = end_1 - j1_is_overlap_length - 1;
        string j1_right_is_sequence = ref_seq_info.get_sequence_1 (
          j1[j1["_is_interval"] + "_seq_id"],
          start_1,
          end_1
				);
				j1_right_is_sequence = reverse_complement(j1_right_is_sequence);

				if (verbose) {
					cout << "J1 LEFT : " << j1_left_is_sequence << endl;
					cout << "J1 RIGHT: " << j1_right_is_sequence << endl;
				}

				// believe the direction if the sequences are different
				bool j1_is_ambiguous = (j1_left_is_sequence == j1_right_is_sequence);

        start_1 = n(j2["_" + j2["_is_interval"] + "_is_start"]);
        end_1   = start_1 +j2_is_overlap_length - 1;
        string j2_left_is_sequence = ref_seq_info.get_sequence_1 (
          j2[j2["_is_interval"] + "_seq_id"],
          start_1,
          end_1
				);

        end_1   = n(j2["_" + j2["_is_interval"] + "_is_end"]);
        start_1 = end_1 - j2_is_overlap_length - 1;
        string j2_right_is_sequence = ref_seq_info.get_sequence_1 (
          j2[j2["_is_interval"] + "_seq_id"],
          start_1,
          end_1
				);
				j2_right_is_sequence = reverse_complement(j2_right_is_sequence);

				// believe the direction if the sequences are different
				bool j2_is_ambiguous = (j2_left_is_sequence == j2_right_is_sequence);

				
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
				// $j1 may be the left side, rather than the right side of the insertion, if so...
				if (uc1_strand)
				{
					if (verbose) cout << "reverse right and left" << endl;

          // @JEB can't seem to get swap() to work here.          
          swap(mut["_ins_start"], mut["_ins_end"]);
          swap(mut["_del_start"], mut["_del_end"]);          
          /*
          diff_entry_value_t temp;
          
          temp = mut["_ins_start"];
          mut["_ins_start"] = mut["_ins_end"];
          mut["_ins_end"] = temp;
          
          temp = mut["_del_start"];
          mut["_del_start"] = mut["_del_end"];
          mut["_del_end"] = temp;
           */
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
        
				gd.add(mut);
				break; // next JC1
			}
		}

		///
		// evidence JC => INS, SUB, AMP mutations
		///
    for(diff_entry_list_t::iterator jc_it = jc.begin(); jc_it != jc.end(); jc_it++) //JC
		{
			cDiffEntry& j = **jc_it;

      //cout << j["side_1_seq_id"] << " " << j["side_2_seq_id"] << endl;
      //cout << j["side_1_strand"] << " " << j["side_2_strand"] << endl;
      
			if (
				// must be on same sequence
				   (j["side_1_seq_id"] != j["side_2_seq_id"])
				// must be in same orientation (implies strands are opposite)
				|| (n(j["side_1_strand"]) == n(j["side_2_strand"]))
			)
				continue;

			string seq_id = j["side_1_seq_id"];
      int32_t side_1_position = n(j["side_1_position"]);
      int32_t side_2_position = n(j["side_2_position"]);
      
			// We can assume that the lower coordinate will be first since this is NOT a deletion
			// (which would be handled above)
      // By this point any positive overlap should have been resolved.
			assert(n(j["overlap"]) <= 0);
      
			// mutation will always be after this position
			// Special case of circular chromosome
			if ( (side_1_position == 1) && ( side_2_position== ref_seq_info[ref_seq_info.seq_id_to_index(j["side_2_seq_id"])].m_length ) )
			{
				j["circular_chromosome"] = "1";
				continue;
			}
      
      // If we are predicting a very big insertion (longer than read length), 
      // it is likely spurious. Require other evidence to convert to a mutation.
      
      // side_2_position will always be greater than or equal to size_i_position
      ASSERT(side_2_position >= side_1_position, "Unexpected positions");
      int32_t size = side_2_position - side_1_position + 1;
      
      // This implies a deletion, and not counting the endpoint nucleotides (hence -1 instead of +1)
      // Note that this is a negative number!
      if (n(j["side_1_strand"]) == -1 )
        size = side_1_position - side_2_position + 1;
      
      // Insertion or deletion must be smaller than read length to be predicted
      // by this evidence alone.
      if (abs(size) > avg_read_length) 
        continue;

			// 'DEL' or 'AMP'
			if (!j.entry_exists("unique_read_sequence"))
			{        
				if (size < 0) // this is a deletion!
        {
          cDiffEntry mut;
          mut._type = DEL;
          mut
          ("seq_id", seq_id)
          ("position", s(side_1_position+1))
          ("size", s(-size))         // note adjustment due to +1 above
          ;
          mut._evidence = make_vector<string>(j._id);
          gd.add(mut);
        } 
        else // this is an amplification.		
        {		
          cDiffEntry mut;		
          mut._type = AMP;		
          mut		
          ("seq_id", seq_id)		
          ("position", s(side_1_position))		
          ("size", s(size))		
          ("new_copy_number", "2")		
          ;		
          mut._evidence = make_vector<string>(j._id);
          gd.add(mut);		
        }
      }
			// 'INS'
      //  INS predicted here are aligned with missing unique_read_sequence info.
      //  We need to grab it.
			else if (n(j["side_1_position"]) >= n(j["side_2_position"]))
			{
				string new_seq = j["unique_read_sequence"];        
        string ref_seq = ref_seq_info.get_sequence_1(seq_id, n(j["side_2_position"]), n(j["side_1_position"]));

				cDiffEntry mut;
        mut._type = INS;
				mut
					("seq_id", seq_id)
					("position", s(side_1_position))
					("new_seq", new_seq + ref_seq)
				;
        mut._evidence = make_vector<string>(j._id);

				gd.add(mut);
			}
      // 'INS'
      //  INS predicted here are aligned with missing unique_read_sequence info.
      //  Further the unique read sequence is aligned on the reverse strand.
      //  We need to grab it, proper like.
			else if ((n(j["side_2_strand"])  < 0) && (n(j["side_1_position"]) <= n(j["side_2_position"])))
			{
				string new_seq = reverse_complement(j["unique_read_sequence"]);
        string ref_seq = ref_seq_info.get_sequence_1(seq_id, n(j["side_1_position"]), n(j["side_2_position"]));
				cDiffEntry mut;
        mut._type = INS;
				mut
        ("seq_id", seq_id)
        ("position", j["side_2_position"])
        ("new_seq", new_seq + ref_seq)
				;
        mut._evidence = make_vector<string>(j._id);
        
				gd.add(mut);
			}
			// "INS" || "AMP"
      // @JEB Note: this kind of AMP only happens due to the way that reads are currently split.
      //            If that is resolved, then only the version above will be needed for AMP.
			else if (n(j["side_1_position"]) + 1 == n(j["side_2_position"]))
			{
        // Check to see if unique sequence matches sequence directly before
        size_t size = j["unique_read_sequence"].size();
        size_t prev_position = n(j["side_2_position"]) - size;
        string dup_check_seq = ref_seq_info.get_sequence_1(j["side_1_seq_id"], prev_position, prev_position + size - 1);
        
        if (j["unique_read_sequence"] == dup_check_seq)
        {
          cDiffEntry mut;
          mut._type = AMP;
          mut
          ("seq_id", seq_id)
          ("position", s(n(j["side_2_position"]) - j["unique_read_sequence"].size()))
          ("size", s(j["unique_read_sequence"].size()))
          ("new_copy_number", "2")
          ;
          mut._evidence = make_vector<string>(j._id);
          gd.add(mut);
        }
        else
        {
          cDiffEntry mut;
          mut._type = INS;
          mut
            ("seq_id", seq_id)
            ("position", s(side_1_position))
            ("new_seq", j["unique_read_sequence"])
          ;
          mut._evidence = make_vector<string>(j._id);
          gd.add(mut);
        }
			}
		}


		///
		// Read Alignments => SNP, DEL, INS, SUB
		///
    
    vector<gd_entry_type> ra_types = make_vector<gd_entry_type>(RA);
    diff_entry_list_t ra = gd.list(ra_types);

		///
		// Ignore RA that overlap DEL or MC
		// They are due to low spurious coverage in deleted regions!
		///

		{
      vector<gd_entry_type> del_types = make_vector<gd_entry_type>(DEL);
      vector<gd_entry_type> mc_types = make_vector<gd_entry_type>(MC);
			diff_entry_list_t del = gd.list(del_types);
      diff_entry_list_t mc = gd.list(mc_types);

      for(diff_entry_list_t::iterator ra_it = ra.begin(); ra_it != ra.end(); ra_it++) //RA
      {
        cDiffEntry& ra_item = **ra_it;

				bool next_ra = false;
        
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
    ra.remove_if(cDiffEntry::field_exists("reject"));

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
			if ( item.entry_exists("reject") || item.entry_exists("deleted"))
			  continue;

      // If we are predicting mixed bases and not polymorphisms, then don't create
      // mutations for mixed frequency predictions (leave them as unassigned RA evidence)
      if (settings.mixed_base_prediction && (item["frequency"] != "1"))
        continue;
      
			bool same = false;
			if (!first_time)
			{
				if ( ((mut["end"] == item["position"]) && (n(mut["insert_end"]) + 1 == n(item["insert_position"])))
						|| ((n(mut["end"]) + 1 == n(item["position"])) && (item["insert_position"] == "0")) )
					same = true;
				if ( (item["frequency"] != "1") || (mut["frequency"] != "1") //don't join polymorphisms
						|| (mut["seq_id"] != item["seq_id"]) )
					same = false;
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
					("frequency", item["frequency"])
				;
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

		for (uint32_t i = 0; i < muts.size(); i++)
		{
			mut = muts[i];
			// insertion and amplifiction
			if (mut["ref_seq"].size() == 0)
			{
        mut._type = INS;
				// unused fields
				mut.erase("ref_seq");
        
        // Check to see if unique sequence matches sequence directly before
        string dup_check_seq = ref_seq_info.get_circular_sequence_1(mut["seq_id"], n(mut["position"]) - (mut["new_seq"].size() - 1), mut["new_seq"].size());
        
        if((mut["new_seq"] == dup_check_seq) && (mut["new_seq"].size() > 1))
        {
          mut._type = AMP;
          mut["size"] = s(mut["new_seq"].size());          
          mut["new_copy_number"] = "2";
          mut.erase("new_seq");
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

			// we don't need these fields
			if (mut["frequency"] == "1") mut.erase("frequency");
			mut.erase("start");
			mut.erase("end");
			mut.erase("insert_start");
			mut.erase("insert_end");

			gd.add(mut);
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
  ("INTERGENIC")("NONCODING")("PSEUDOGENE")("SYNONYMOUS")("NONSYNONYMOUS")("NO_CHANGE")("TOTAL")  
  ;
 
   vector<string>  BaseSubstitutionEffects::base_type_list = make_vector<string>
  ("INTERGENIC")("NONCODING")("PSEUDOGENE")("PROTEIN")("TOTAL")
  ;  
  
  vector<string> BaseSubstitutionEffectCounts::base_pair_change_count_list = make_vector<string> 
  ("AT.GC")("AT.CG")("AT.TA")("CG.TA")("CG.AT")("CG.GC")("TOTAL")
  ;
  
  vector<string> BaseSubstitutionEffectCounts::base_change_type_count_list = make_vector<string> 
  ("INTERGENIC")("NONCODING")("PSEUDOGENE")("SYNONYMOUS")("NONSYNONYMOUS")("TOTAL")  
  ;
  
  map<BaseSubstitutionEffect,BaseType>  BaseSubstitutionEffects::snp_type_to_base_type = make_map<BaseSubstitutionEffect,BaseType>
  (intergenic_base_substitution,intergenic_base)
  (pseudogene_base_substitution,pseudogene_base)
  (noncoding_base_substitution,noncoding_base)
  (synonymous_base_substitution,protein_base)
  (nonsynonymous_base_substitution,protein_base)
  ;  
  
  
  void BaseSubstitutionEffects::initialize_from_sequence(cReferenceSequences& ref_seq_info) 
  {    
    bool count_synonymous_stop_codons = true;
    bool verbose = false;
    
    map<string,string> codon_synonymous_changes;
    map<string,string> codon_nonsynonymous_changes;
    map<string,string> codon_num_synonymous_changes;
    map<string,string> codon_position_mutation_synonymous;
    
    map<string,string> nonsynonymous_mutations;
    map<string,string> synonymous_mutations;
    
    uint32_t total_num_synonymous_changes = 0;
    uint32_t total_num_nonsynonymous_changes = 0;
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
        if (f[TYPE] == "gene") {
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
        if ((f[TYPE] != "CDS") || f.m_pseudo)
          continue;
        
        // initialize gene structure
        Gene g(f);
        
        total_orfs++;
        
        // Piece together the gene - at each nucleotide make all three possible changes
        
        string this_codon = "   ";
        vector<cLocation> sub_locations = g.m_location.get_all_sub_locations();
        size_t on_codon_index = 0;
        vector<uint32_t> this_codon_locations_0(0, 3); // 0-indexed
        
        uint32_t total_nucleotide_length = 0; 
        int8_t strand = sub_locations.front().m_strand;
        
        for(vector<cLocation>::iterator it3=sub_locations.begin(); it3!=sub_locations.end(); ++it3) {
          cLocation& loc = *it3;
          total_nucleotide_length += loc.m_end - loc.m_start + 1;
        }
        uint32_t total_amino_acid_length = total_nucleotide_length / 3;
        uint32_t on_codon_pos_1 = 1;
        if (strand == -1) on_codon_pos_1 = total_amino_acid_length;

        /// code to debug consistency with Perl
        /*
        string amino_acid_sequence;
        string codon_sequence;
        */
        
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
              if ((loc.get_strand() == 1) && (seq_bcs[pos_0] == reverse) )
                seq_bcs[pos_0] = conflict;
              if ((loc.get_strand() == -1) && (seq_bcs[pos_0] == forward) )
                seq_bcs[pos_0] = conflict;
            }
            else
            {
              seq_bcs[pos_0] = (loc.get_strand() == +1 ? forward : reverse);
            }
            
            //// Handle codon synonymous/nonsynonymous changes
            this_codon_locations_0[on_codon_index] = pos_0;
            this_codon[on_codon_index] = seq.get_sequence_1(pos_1);
            
            on_codon_index++;
            
            // The codon is filled, now make all mutations and assign to proper nucleotides
            if (on_codon_index == 3) {
              
              // change the codon to the right strand
              string original_codon = this_codon;
              if (strand == -1)
                original_codon = reverse_complement(original_codon);
              
              char original_amino_acid = cReferenceSequences::translate_codon(original_codon, g.translation_table, on_codon_pos_1);
              
              /*
              /// code to debug consistency with Perl
              if (strand == -1) {
                codon_sequence = original_codon + codon_sequence;
                amino_acid_sequence = original_amino_acid + amino_acid_sequence;
              } else {
                codon_sequence = codon_sequence + original_codon;
                amino_acid_sequence = amino_acid_sequence + original_amino_acid;
              }
              /// end debug code
              */
               
              for (int32_t test_codon_index=0; test_codon_index<3; test_codon_index++) {
                
                for (size_t b=0; b<BaseSubstitutionEffects::base_char_list.size(); b++) {
                  
                  char mut_base = BaseSubstitutionEffects::base_char_list[b];
                  
                  string test_codon = this_codon;
                  test_codon[test_codon_index] = mut_base;
                  
                  if (strand == -1)
                    test_codon = reverse_complement(test_codon);
                  
                  char mut_amino_acid = cReferenceSequences::translate_codon(test_codon, g.translation_table, on_codon_pos_1);
                  
                  if (mut_amino_acid == original_amino_acid)
                    seq_bse[this_codon_locations_0[test_codon_index]*4+b] = max(seq_bse[this_codon_locations_0[test_codon_index]*4+b], synonymous_base_substitution);
                  else
                    seq_bse[this_codon_locations_0[test_codon_index]*4+b] = max(seq_bse[this_codon_locations_0[test_codon_index]*4+b], nonsynonymous_base_substitution);
                  
                }
                
              }
              
              on_codon_pos_1 += strand;
              on_codon_index = 0;
            }
          }
        } // end sublocation loop
        ASSERT(on_codon_index == 0, "Number of base pairs in CDS not a multiple of 3: " + g["name"]);
        
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
                         bool base_substitution_statistics
                         )
  {
    // Could be a parameter > this is a "large" mutation, <= this is an "small" mutation
    int32_t large_size_cutoff = 50;
      
    // Figure out the names of all "repeat" columns
    map<string,bool> mob_name_hash;
    map<string,bool> con_name_hash;
    
    for (vector<cGenomeDiff>::iterator it=genome_diffs.begin(); it != genome_diffs.end(); ++it) {
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
    column_headers.push_back("sample");
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
    
    
    output_file << join(column_headers, ",") << endl;
    
    for (vector<cGenomeDiff>::iterator it=genome_diffs.begin(); it != genome_diffs.end(); ++it) {
      cGenomeDiff &gd = *it;
      //uout("Counting mutations " + gd.metadata.run_name);
      
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
      diff_entry_list_t mut_list = gd.mutation_list();
      for (diff_entry_list_t::iterator it=mut_list.begin(); it != mut_list.end(); ++it) {		
        
        cDiffEntry& mut = **it;
        
        // Don't count mutations that were hidden by later deletions, but kept for phylogenetic inference.
        if (mut.is_marked_deleted()) continue;
        
        if (mut._type == SNP) {
          count["base_substitution"][""]++;
          if (base_substitution_statistics) {
            this_bsec.change_position_1_observed_totals(ref_seq_info, bse, mut[SEQ_ID], from_string<int32_t>(mut[POSITION]), mut[NEW_SEQ], +1);					
          }
          count["type"][mut["snp_type"]]++;
        }
        
        if (mut._type == DEL) {
          total_deleted += from_string<int32_t>(mut[SIZE]);
          
          if (mut.entry_exists("mediated"))
            count["mob"][mut["mediated"]]++;
          
          if (from_string<int32_t>(mut[SIZE]) > large_size_cutoff)
            count["large_deletion"][""]++;
          else
            count["small_indel"][""]++;
        }
        
        if (mut._type == INS) {
          int32_t ins_size = mut[NEW_SEQ].size();
          
          total_inserted += ins_size;
          
          if (mut.entry_exists("mediated"))
            count["mob"][mut["mediated"]]++;
          
          if (ins_size > large_size_cutoff)
            count["large_insertion"][""]++;
          else
            count["small_indel"][""]++;
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
          
          if (abs(new_size - old_size) > large_size_cutoff)
            count["large_substitution"][""]++;
          else
            count["small_indel"][""]++;
        }
        
        if (mut._type == CON) {
          // TODO: Need size change?
          
          if (mut.entry_exists("mediated"))
            count["con"][mut["mediated"]]++;
          count["gene_conversion"][""]++;
        }
        
        if (mut._type == MOB) {
          int32_t rpos = -1;
          string repeat_seq = ref_seq_info.repeat_family_sequence(mut["repeat_name"], 1, rpos);
          
          int32_t this_length = repeat_seq.size();
          total_inserted += this_length;
          total_repeat_inserted += this_length;
          //print "Repeat $mut->{repeat_name} $this_length bp\n";
          
          count["mob"][mut["repeat_name"]]++;
          count["mobile_element_insertion"][""]++;
          
        } 
        
        if (mut._type == AMP) {
          int32_t this_size = mut.mutation_size_change(ref_seq_info);
          total_inserted += this_size;
          
          if (this_size > large_size_cutoff)
            count["large_amplification"][""]++;
          else
            count["small_indel"][""]++;
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
      
      vector<string> this_columns;
      
      this_columns.push_back(gd.metadata.run_name);
      this_columns.push_back(to_string(mut_list.size()));
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
  }
  

} // namespace breseq
