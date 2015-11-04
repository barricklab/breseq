/*****************************************************************************
 *
 * AUTHORS
 *
 *    Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 *    David B. Knoester
 *
 * LICENSE AND COPYRIGHT
 *
 *    Copyright (c) 2008-2010 Michigan State University
 *    Copyright (c) 2011-2012 The University of Texas at Austin
 *
 *    breseq is free software; you can redistribute it and/or modify it under the
 *    terms the GNU General Public License as published by the Free Software
 *    Foundation; either version 1, or (at your option) any later version.
 *
 *****************************************************************************/

#include "libbreseq/output.h"

#include "libbreseq/alignment_output.h"

using namespace std;

namespace breseq
{
  using namespace breseq::output;

bool verbose = false; //TODO Options
bool text = false; //TODO Options

char  alignment_output::s_internal_gap_character = '-';
char  alignment_output::s_end_gap_character = '-';
char  alignment_output::s_reference_match_character = '.';
char  alignment_output::s_reference_match_masked_character = ',';

uint8_t alignment_output::k_reserved_quality_junction_overlap = 253;
uint8_t alignment_output::k_reserved_quality_dont_highlight = 254;
uint8_t alignment_output::k_reserved_quality_max = 255;

  
alignment_output::Alignment_Output_Pileup::Alignment_Output_Pileup ( 
                                                                    const string& bam, 
                                                                    const string& fasta, 
                                                                    const bool show_ambiguously_mapped 
                                                                    )
        : pileup_base ( bam, fasta )
        , _show_ambiguously_mapped(show_ambiguously_mapped)
        , unique_start ( 0 )
        , unique_end ( 0 )
{
}

alignment_output::Alignment_Output_Pileup::~Alignment_Output_Pileup() {}

/*! Called for each position.*/
void alignment_output::Alignment_Output_Pileup::pileup_callback( const pileup& p )
{
  uint32_t& start_1( m_start_position_1 );
  uint32_t& end_1( m_end_position_1 );
  uint32_t reference_pos_1 = p.position_1();
  
  if ( verbose )
  {
    cout << "POSITION: " << reference_pos_1 << endl;
    cout << "ALIGNED:  " << p.size() << endl;
  }
  
  // Don't pay attention if no reads uniquely aligned here
  if ( ( reference_pos_1 < unique_start ) || ( reference_pos_1 > unique_end ) )
    return;
  
  // update beginning and end of each reference
  for (vector<Aligned_Reference>::iterator it=aligned_references.begin(); it!=aligned_references.end(); it++)
  {
    Aligned_Reference& aligned_reference = *it;
    if(aligned_reference.start == 0) aligned_reference.start = reference_pos_1;
    aligned_reference.end = reference_pos_1;
  }
  
  //## Cull the list to those we want to align
  //## Other alignments may overlap this position that DO NOT
  //## overlap the positions of interest. This removes them.
  
  vector<const pileup_wrapper*> alignments;
  for ( pileup::const_iterator itr_pileup = p.begin(); itr_pileup != p.end() ; itr_pileup++ )
  {
    const pileup_wrapper& a = *itr_pileup;
    if (aligned_reads.count(a.read_name())) 
    {
      alignments.push_back(&(*itr_pileup));
    }
  }
  
  int32_t temp_max_indel = 0;
  for ( vector<const pileup_wrapper*>::const_iterator itr_pileup = alignments.begin(); itr_pileup != alignments.end() ; itr_pileup ++ )
  {
    const pileup_wrapper& a = **itr_pileup;
    
    if (a.on_base_indel() > temp_max_indel )
    {
      temp_max_indel = a.on_base_indel();
    }
  }
  assert(temp_max_indel >= 0);
  int32_t max_indel = static_cast<uint32_t>(temp_max_indel);
  
  if ( verbose ) cout << "MAX INDEL: " << max_indel << endl;
  
  // ## Now add this position to the alignments
  
  /// This flag is set if the read was among those returned by the pileup
  /// if it is not set, we insert padding after the loop through the alignments!
  for (Aligned_Reads::iterator itr_read = aligned_reads.begin(); itr_read != aligned_reads.end(); itr_read++ )
  {
    Aligned_Read& aligned_read ( itr_read->second );
    aligned_read.updated = false;
  }
  
  // MAIN ALIGNMENT LOOP
  for ( vector<const pileup_wrapper*>::const_iterator itr_pileup = alignments.begin(); itr_pileup != alignments.end() ; itr_pileup ++ )
  {
    const pileup_wrapper& a = **itr_pileup;
    
    Aligned_Read& aligned_read = aligned_reads[a.read_name()];
    aligned_read.updated = true;
    
    // -1 for gap in this read, >0 for this many bases inserted in read
    int32_t indel = a.on_base_indel();      
    
    // Update the beginning and the end coords of what we are showing
    if (aligned_read.reference_start == 0) aligned_read.reference_start = reference_pos_1;
    aligned_read.reference_end = reference_pos_1;
    
    if ( aligned_read.start == 0 ) aligned_read.start = a.query_position_1();
    aligned_read.end = a.query_position_1();
    
    // ## READS: add aligned positions
    // FOR EACH INSERT COUNT
    for ( int32_t index = 0; index <= max_indel; index ++ )
    {
      // add gap or space
      if ( index > indel )
      {
        if ( verbose )
        {
          cout << "Adding gap: " << indel << endl;
        }
        
        // If this is a gap at the last position, then add space
        if (reference_pos_1 == a.reference_end_1())
        {
          aligned_read.aligned_bases += ' ';
          aligned_read.aligned_quals += char ( k_reserved_quality_dont_highlight );
        }
        else
        {
          aligned_read.aligned_bases += s_internal_gap_character;
          aligned_read.aligned_quals += char ( k_reserved_quality_max );
        }
      }
      // add the aligned base
      else
      {
        uint8_t quality = a.read_base_quality_0( a.query_position_0() + index );
        char base = a.read_base_char_0( a.query_position_0() + index );
        
        if ( !text )
        {
          uint8_t trim_left = a.trim_left();
          if ( ( trim_left != 0 ) && ( a.query_position_1() <= trim_left ) )
          {
            base = tolower ( base );           
          }
          
          uint8_t trim_right = a.trim_right();
          if ( ( trim_right != 0 ) && ( ( a.read_length() - a.query_position_0() ) <= trim_right ) )
          {
            base = tolower ( base );
          }
        }
        
        aligned_read.aligned_bases += base;
        aligned_read.aligned_quals += char ( quality );
      }
      
      if ( verbose )
      {
        //        cout << aligned_read.aligned_bases << " ";
        //        cout << aligned_read.seq_id << " ";
        //        cout << a.on_base_indel() << endl;
      }
    } // END FOR EACH INSERT COUNT
  }
  //END MAIN ALIGNMENT LOOP
  
  // ## READS: handle those with no base aligned to this position
  for ( Aligned_Reads::iterator itr_read = aligned_reads.begin();  itr_read != aligned_reads.end(); itr_read++ )
  {    
    Aligned_Read& aligned_read (itr_read->second );
    
    if ( aligned_read.updated ) continue;
    
    aligned_read.aligned_bases += repeat_char(' ', max_indel+1);
    aligned_read.aligned_quals += repeat_char(char(k_reserved_quality_dont_highlight), max_indel+1);
    if ( verbose ) cout << aligned_read.aligned_bases << " NOALIGN " << itr_read->first << endl;
  }
  
  
  // ##now handle the reference sequences
  char ref_base = p.reference_base_char_1(reference_pos_1);
  
  for ( uint32_t index = 0; index < aligned_references.size(); index ++ )
  {    
    Aligned_Reference& aligned_reference ( aligned_references[index] );    
    
    // Default is to show reference 
    char my_ref_base = ref_base;
    char my_ref_qual = char(k_reserved_quality_max); // highlighted as max quality
    
    // handle truncated starts
    if (aligned_reference.truncate_start) {
      // past part we want to show
      if (reference_pos_1 < aligned_reference.truncate_start - aligned_reference.overlap_not_assigned) {
        my_ref_base = s_end_gap_character; 
        my_ref_qual = char(k_reserved_quality_junction_overlap);
      }
      // in overlap where we don't want to highlight
      else if (reference_pos_1 < aligned_reference.truncate_start) {
        my_ref_base = tolower(my_ref_base);
        my_ref_qual = char(k_reserved_quality_junction_overlap);
      }
    }
    
    // handle truncated ends
    if (aligned_reference.truncate_end) {
      // past part we want to show
      if (reference_pos_1 > aligned_reference.truncate_end + aligned_reference.overlap_not_assigned) {
        my_ref_base = s_end_gap_character; 
        my_ref_qual = char(k_reserved_quality_junction_overlap);
      }
      // in overlap where we don't want to highlight
      else if (reference_pos_1 > aligned_reference.truncate_end) {
        my_ref_base = tolower(my_ref_base);
        my_ref_qual = char(k_reserved_quality_junction_overlap);
      }
    }
    
    aligned_reference.aligned_bases += my_ref_base;
    aligned_reference.aligned_quals += my_ref_qual;
    
    aligned_reference.aligned_bases += repeat_char(s_internal_gap_character, max_indel);
    aligned_reference.aligned_quals += repeat_char(char(k_reserved_quality_max), max_indel);
  }
  
  // ##also update any positions of interest for gaps
  for ( uint32_t insert_count = 0; insert_count <= static_cast<uint32_t>(max_indel); insert_count++ )
  {
    if (   ((m_insert_start <= insert_count) && (insert_count <= m_insert_end) && (reference_pos_1 == start_1) && (reference_pos_1 == end_1))
        || ((m_insert_start <= insert_count) && (reference_pos_1 == start_1) && (reference_pos_1 != end_1)) 
        || ((m_insert_end >= insert_count) && (reference_pos_1 == end_1) && (reference_pos_1 != start_1))
        || ((reference_pos_1 < end_1) && (reference_pos_1 > start_1)) )
    {
      aligned_annotation.aligned_bases += '|';
    }
    else
    {		
      aligned_annotation.aligned_bases += ' ';
    }
  }
}
//END create_alignment()

/*! Called for each read alignment.*/
void alignment_output::Alignment_Output_Pileup::fetch_callback ( const alignment_wrapper& a )
{
  // we only keep track of unique alignments 
  if ( !_show_ambiguously_mapped && a.is_redundant() ) return;
  
  // create a new empty structure and fill it  
  Aligned_Read aligned_read;
  aligned_read.seq_id = a.read_name();
  aligned_read.strand = a.strand();
  aligned_read.length = a.read_length();
  aligned_read.read_sequence = a.read_char_sequence();
  aligned_read.qual_sequence = a.read_base_quality_bam_string();
  aligned_read.is_redundant = a.is_redundant();
  aligned_read.read_name_style =  aligned_read.is_redundant ? "AM" : "NC";
  
  aligned_read.mapping_quality = a.mapping_quality();
  
  aligned_read.reference_start = a.reference_start_1();
  aligned_read.reference_end = a.reference_end_1();
  
  aligned_reads[aligned_read.seq_id] = aligned_read;
  
  if (verbose)
  {
    cout << a.read_name() << endl;
    cout << a.reference_start_1() << " " << a.reference_end_1() << endl;
  }
  
  if ( ( unique_start == 0 ) || ( unique_start > a.reference_start_1() ) )
  {
    unique_start = a.reference_start_1();
  }
  if ( ( unique_end == 0 ) || ( unique_end < a.reference_end_1() ) )
  {
    unique_end = a.reference_end_1();
  }
}

  
alignment_output::alignment_output ( string bam, 
                                    string fasta, 
                                    uint32_t maximum_to_align, 
                                    const uint32_t quality_score_cutoff,
                                    const int32_t junction_minimum_size_match,
                                    const bool mask_ref_matches,
                                    const bool show_ambiguously_mapped
                                    )
        : m_alignment_output_pileup ( bam, fasta, show_ambiguously_mapped )
        , m_aligned_reads ( m_alignment_output_pileup.aligned_reads )
        , m_quality_score_cutoff ( quality_score_cutoff )
        , m_maximum_to_align ( maximum_to_align )
        , m_junction_minimum_size_match ( junction_minimum_size_match )
        , m_mask_ref_matches (mask_ref_matches) 
        , m_show_ambiguously_mapped ( show_ambiguously_mapped )
        , m_is_junction(false)
        , m_is_junction_junction(false)
{
  // zero is the default if no arg provided, make a reasonable value
  if (m_maximum_to_align == 0)
    m_maximum_to_align = 200;
  
  // the default for m_quality_score_cutoff is zero, which makes sense
}

void alignment_output::create_alignment ( const string& region, cOutputEvidenceItem* output_evidence_item_ptr )
{
  // we need the target_id to properly fill out the reference sequence later
  uint32_t target_id, start_pos, end_pos, insert_start, insert_end;
  m_alignment_output_pileup.parse_region(region, target_id, start_pos, end_pos, insert_start, insert_end);
  
  // Check for special reference lines that are junctions...
  bool bDrawAnnotationLine = true;    // draw the | annotation line pointing to region?
  
  
  // Check for special output in reference and legend lines for junctions
  m_is_junction = false;
  m_is_junction_junction = false;

  if (output_evidence_item_ptr != NULL) {
    if ((*output_evidence_item_ptr)[PREFIX] == "JC") {
      m_is_junction = true;
      m_is_junction_junction = true;
    } else if ( ((*output_evidence_item_ptr)[PREFIX] == "JC_SIDE_1")
             || ((*output_evidence_item_ptr)[PREFIX] == "JC_SIDE_2") ) {
       m_is_junction = true;
    }
  }
  
  // For junctions JC evidence to set up split reference
  if (m_is_junction_junction) {
    cDiffEntry& jc_item = *(output_evidence_item_ptr->item);
    
    int32_t overlap = from_string<int32_t>(jc_item["overlap"]);
    int32_t alignment_overlap = from_string<int32_t>(jc_item["alignment_overlap"]);
    int32_t positive_alignment_overlap = max(alignment_overlap, 0);
    int32_t abs_overlap = max(abs(alignment_overlap), abs(overlap));
    
    // No need to draw annotation in this case.
    if(overlap >= 0)
      bDrawAnnotationLine = false;
    
    // ghost_* are the coordinates (and strands) to assign to the innermost highlighted base in each case 
    // (in reference coords)
    // truncate_* says to stop showing the reference after or before this position (in junction coords)
    // overlap basses extend past these coords, are not highlighted, and not numbered
    //
    // Translating from the coordinates of the junction candidate to how the overlap was resolved creates some difficulties here
        
    Aligned_Reference aligned_reference_1, aligned_reference_2;
 
    // This is the amount of a positive overlap assigned to each side
    aligned_reference_1.overlap_assigned = from_string<int32_t>(jc_item["side_1_overlap"]);
    aligned_reference_2.overlap_assigned = from_string<int32_t>(jc_item["side_2_overlap"]);
    
    aligned_reference_1.overlap_not_assigned = positive_alignment_overlap - aligned_reference_1.overlap_assigned;
    aligned_reference_2.overlap_not_assigned = positive_alignment_overlap - aligned_reference_2.overlap_assigned;
    
    aligned_reference_1.truncate_end = from_string<int32_t>(jc_item["flanking_left"]) + aligned_reference_1.overlap_assigned; 
    aligned_reference_2.truncate_start = from_string<int32_t>(jc_item["flanking_left"]) + 1 + abs_overlap - aligned_reference_2.overlap_assigned;
    
    aligned_reference_1.ghost_seq_id = jc_item[SIDE_1_SEQ_ID];    
    aligned_reference_2.ghost_seq_id = jc_item[SIDE_2_SEQ_ID]; 
    
    aligned_reference_1.ghost_strand = from_string<int16_t>(jc_item[SIDE_1_STRAND]);
    aligned_reference_2.ghost_strand = from_string<int16_t>(jc_item[SIDE_2_STRAND]);
    
    aligned_reference_1.ghost_end = from_string<int32_t>(jc_item[SIDE_1_POSITION]); 
    aligned_reference_2.ghost_start = from_string<int32_t>(jc_item[SIDE_2_POSITION]);     
    
    //cout << "Junction reference information" << endl;
    //cout << aligned_reference_1.truncate_start << "-" << aligned_reference_1.truncate_end << " overlap " << aligned_reference_1.overlap << endl;
    //cout << aligned_reference_2.truncate_start << "-" << aligned_reference_2.truncate_end << " overlap " << aligned_reference_2.overlap << endl;

    // common settings
    aligned_reference_1.target_id = target_id;
    aligned_reference_1.seq_id = m_alignment_output_pileup.target_name(target_id);
    aligned_reference_1.reference_length = m_alignment_output_pileup.target_length(target_id);
    aligned_reference_2.target_id = target_id;
    aligned_reference_2.seq_id = m_alignment_output_pileup.target_name(target_id);
    aligned_reference_2.reference_length = m_alignment_output_pileup.target_length(target_id);
    
    m_alignment_output_pileup.aligned_references.push_back(aligned_reference_1);
    m_alignment_output_pileup.aligned_references.push_back(aligned_reference_2);
  }
  else
  {
    ///Single Reference case
    Aligned_Reference aligned_reference;
    aligned_reference.target_id = target_id;
    aligned_reference.seq_id = m_alignment_output_pileup.target_name(target_id);
    aligned_reference.reference_length = m_alignment_output_pileup.target_length(target_id);
    aligned_reference.strand = +1;
    m_alignment_output_pileup.aligned_references.push_back(aligned_reference);
  }

  // Use fetch to determine the list of reads that will be aligned to this position
  m_alignment_output_pileup.set_print_progress(false);
  m_alignment_output_pileup.do_fetch ( region );

  if ( ( m_alignment_output_pileup.unique_start == 0 ) || ( m_alignment_output_pileup.unique_end == 0 ) )
  {
    m_error_message = "No reads uniquely aligned to region.";
    return;
  }

  // Reduce the number of reads to the maximum that we want to align
  Aligned_Reads& aligned_reads = m_alignment_output_pileup.aligned_reads;
  if ( (m_maximum_to_align != 0) && (aligned_reads.size() > m_maximum_to_align) ) {
    Aligned_Reads new_aligned_reads;
    m_error_message = "Only " + to_string(m_maximum_to_align) + " of " + to_string(aligned_reads.size()) + " total aligned reads displayed.";
    
    // create a list of the keys
    vector<string> key_list;
    for(Aligned_Reads::iterator ar = aligned_reads.begin(); ar != aligned_reads.end(); ar++) {
      key_list.push_back(ar->first);
    }
    
    double max = static_cast<double>(m_maximum_to_align);
    double num = static_cast<double>(key_list.size());
    // Pick the desired number of reads evenly across the reads returned
    for (size_t i=0; i<m_maximum_to_align; i++) {
      size_t keep_index = static_cast<size_t>(floor(static_cast<double>(i) / max * num));
      new_aligned_reads[key_list[keep_index]] = aligned_reads[key_list[keep_index]];
    }
    
    m_alignment_output_pileup.aligned_reads = new_aligned_reads;
  }  

  // Build the alignment with the pileup
  m_alignment_output_pileup.do_pileup( region );

  // do_pileup populates the following
  m_aligned_reads = m_alignment_output_pileup.aligned_reads;
  m_aligned_references = m_alignment_output_pileup.aligned_references;
  m_aligned_annotation = m_alignment_output_pileup.aligned_annotation;
  
  
  // This section updates the style of the name of the reads that don't cross
  // the ambiguous part of the junction (and therefore are not counted).
  // This applies to junction reads.
  
  
  // For junctions JC evidence to set up split reference
  if (output_evidence_item_ptr != NULL) {
    
    if ((*output_evidence_item_ptr)[PREFIX] == "JC") {
        
      uint32_t require_reference_start = m_alignment_output_pileup.aligned_references[0].truncate_end - m_alignment_output_pileup.aligned_references[0].overlap_assigned - m_junction_minimum_size_match + 1;
      uint32_t require_reference_end = m_alignment_output_pileup.aligned_references[1].truncate_start + m_alignment_output_pileup.aligned_references[1].overlap_assigned + m_junction_minimum_size_match - 1;
      
      for(Aligned_Reads::iterator ar = m_aligned_reads.begin(); ar != m_aligned_reads.end(); ar++) {
        
        Aligned_Read& a = ar->second;
        if ((a.reference_start > require_reference_start) ||  (a.reference_end < require_reference_end)) {
          a.read_name_style = "AO";
        }
      }
      
      // If there was an insertion in the junction, then the annotation line is too long by one 
      // character on each side so, let's fix that here.
      int32_t first_match = -1;
      int32_t last_match = -1;
      
      for(size_t i=0; i<m_aligned_annotation.aligned_bases.size(); i++) {
        if (m_aligned_annotation.aligned_bases[i] == '|') {
          last_match = i;
          if (first_match == -1) first_match = i;
        }
      }
      
      if (first_match != -1)
        m_aligned_annotation.aligned_bases[first_match] = ' ';
      if (last_match != -1)
        m_aligned_annotation.aligned_bases[last_match] = ' ';

    }
    
    else if (   ((*output_evidence_item_ptr)[PREFIX] == "JC_SIDE_1")
             || ((*output_evidence_item_ptr)[PREFIX] == "JC_SIDE_2") ) {   
      cDiffEntry& jc_item = *(output_evidence_item_ptr->item);

      // This stuff is all copied from the code for counting reads...
      int32_t start, end, non_dup_start, non_dup_end;
      
      uint32_t side_1_continuation = from_string<uint32_t>(jc_item["side_1_continuation"]);
      uint32_t side_2_continuation = from_string<uint32_t>(jc_item["side_2_continuation"]);
      
      
      int32_t extra_stranded_require_overlap = 0; 
      if (jc_item.entry_exists("read_count_offset"))
          extra_stranded_require_overlap = from_string<int32_t>(jc_item["read_count_offset"]);
      
      //if (output_evidence_item_ptr->parent_item->_type == MOB) {
      //extra_stranded_require_overlap = from_string<int32_t>((*(output_evidence_item_ptr->parent_item))["duplication_size"]);
      // }

      // This is for requiring a certain number of bases (at least one) to match past the normal point
      // where a read could be uniquely assigned to the junction (or a side)
      int32_t minimum_side_match_correction = m_junction_minimum_size_match - 1;
      
      int32_t alignment_overlap = from_string<int32_t>(jc_item[ALIGNMENT_OVERLAP]);
      int32_t non_negative_alignment_overlap = alignment_overlap;
      non_negative_alignment_overlap = max(0, non_negative_alignment_overlap);
      
      // New side 1
      if ( ((*output_evidence_item_ptr)[PREFIX] == "JC_SIDE_1") ) {
        int32_t side_1_strand = from_string<int32_t>(jc_item[SIDE_1_STRAND]);
        start = from_string<uint32_t>(jc_item[SIDE_1_POSITION]);
        int32_t overlap_correction = non_negative_alignment_overlap - from_string<int32_t>(jc_item[SIDE_1_OVERLAP]);
        
        
        if (side_1_strand == +1) {
          start = start - 1;       
          end = start + 1; 
          start -= overlap_correction;
          end += extra_stranded_require_overlap;
          start -= minimum_side_match_correction;
          end += minimum_side_match_correction;
          end += side_1_continuation;
          
          non_dup_start = start;
          non_dup_end = end - extra_stranded_require_overlap;
        } else {
          //start = start;
          end = start + 1;
          start -= extra_stranded_require_overlap;
          end += overlap_correction;
          start -= minimum_side_match_correction;
          end += minimum_side_match_correction;
          start -= side_1_continuation;
          
          non_dup_start = start + extra_stranded_require_overlap;
          non_dup_end = end;
        }
      }
      
      
      if ( ((*output_evidence_item_ptr)[PREFIX] == "JC_SIDE_2") ) {
        int32_t side_2_strand = from_string<int32_t>(jc_item[SIDE_2_STRAND]);
        start = from_string<uint32_t>(jc_item[SIDE_2_POSITION]);
        int32_t overlap_correction = non_negative_alignment_overlap - from_string<int32_t>(jc_item[SIDE_2_OVERLAP]);
        
        if (side_2_strand == +1) {
          start = start - 1;
          end = start + 1; 
          start -= overlap_correction;
          end += extra_stranded_require_overlap;
          start -= minimum_side_match_correction;
          end += minimum_side_match_correction;
          end += side_2_continuation;
          
          non_dup_start = start;
          non_dup_end = end - extra_stranded_require_overlap;
        } else {
          //start = start;
          end = start + 1;
          start -= extra_stranded_require_overlap;
          end += overlap_correction;
          start -= minimum_side_match_correction;
          end += minimum_side_match_correction;
          start -= side_2_continuation;
          
          non_dup_start = start + extra_stranded_require_overlap;
          non_dup_end = end;
        }
      }
      
      
      for(Aligned_Reads::iterator ar = m_aligned_reads.begin(); ar != m_aligned_reads.end(); ar++) {
        
        Aligned_Read& a = ar->second;
        
        if ((a.reference_start > static_cast<uint32_t>(start)) ||  (a.reference_end < static_cast<uint32_t>(end))) {
          a.read_name_style = "AO";
        }
        
        else if ((a.reference_start > static_cast<uint32_t>(non_dup_start)) ||  (a.reference_end < static_cast<uint32_t>(non_dup_end))) {
          a.read_name_style = "AP";
        }
          
        if (a.seq_id.length() >= 3) {
          if ( a.seq_id.substr(a.seq_id.length()-3, 3) == "-M1")
            a.read_name_style = "AO";
          if ( a.seq_id.substr(a.seq_id.length()-3, 3) == "-M2")
            a.read_name_style = "AO";
        }
        
      }      
    }
    
  }
  
  if(!bDrawAnnotationLine)
    m_aligned_annotation.aligned_bases = substitute(m_alignment_output_pileup.aligned_annotation.aligned_bases, "|" ," ");
   
  // now add the unaligned portions of each
  int32_t max_extend_left = 0;
  int32_t max_extend_right = 0;
  for ( Aligned_Reads::iterator itr_read = m_aligned_reads.begin(); itr_read != m_aligned_reads.end(); itr_read++ )
  {
    Aligned_Read& aligned_read = ( m_aligned_reads[ ( *itr_read ).first] );

    // read the padding off of each sequence
    int32_t left_padding_length = aligned_read.aligned_bases.find_first_not_of ( ' ' );
    int32_t right_padding_length = aligned_read.aligned_bases.length() - aligned_read.aligned_bases.find_last_not_of ( ' ' ) - 1;
         
    if ( verbose )
    {
      cout  << aligned_read.start << " " << aligned_read.end << " " << aligned_read.seq_id << endl;
    }
    
    int32_t extend_left = ( (int32_t)aligned_read.start-1) - left_padding_length;
    int32_t extend_right = ((int32_t)aligned_read.length - (int32_t)aligned_read.end ) - right_padding_length;
          
    if ( extend_left > max_extend_left )
      max_extend_left = extend_left;
    if ( extend_right > max_extend_right )
      max_extend_right = extend_right;
  }
  
  // Now add this much to every read 
  for ( Aligned_Reads::iterator itr_read = m_aligned_reads.begin(); itr_read != m_aligned_reads.end(); itr_read++ )
  {
    Aligned_Read& aligned_read ( itr_read->second );
    
    aligned_read.aligned_bases = repeat_char(' ', max_extend_left) + aligned_read.aligned_bases + repeat_char(' ', max_extend_right);
    aligned_read.aligned_quals = repeat_char(char(k_reserved_quality_max), max_extend_left) + aligned_read.aligned_quals + repeat_char(char(k_reserved_quality_max), max_extend_right);
  }
   
  if ( verbose )
  {
      cout << "Extend: " << max_extend_left;
      cout << " " << max_extend_right << endl;
  }

  // extend reference sequence as requested, be aware of ends of sequence
  // handle left side extending past end of reference sequence
  for ( uint32_t itr_ref = 0; itr_ref < m_aligned_references.size(); itr_ref++ )
  {    
      Aligned_Reference& aligned_reference = m_aligned_references[itr_ref];
    
      uint32_t ref_extend_left = max_extend_left;
      uint32_t ref_extend_right = max_extend_right;
    
      string ref_add_left="";
      string ref_add_right="";
      
      // handle attempts to extend past the beginning of the reference sequence
      if ( static_cast<int32_t>(aligned_reference.start) - max_extend_left < 1 )
      {
          ref_extend_left = (aligned_reference.start -1);
          ref_add_left += repeat_char(s_end_gap_character, max_extend_left - aligned_reference.start + 1);
      }
    
      // handle attempts to extend past end of reference sequence
      if ( ( aligned_reference.end + max_extend_right ) > aligned_reference.reference_length )
      {
          ref_extend_right = aligned_reference.reference_length - aligned_reference.end;
          ref_add_right += repeat_char(s_end_gap_character, aligned_reference.end + max_extend_right - aligned_reference.reference_length);
      }

      uint32_t pos = aligned_reference.start - ref_extend_left;
        
      // create the string to add to the left side
      while ( pos < aligned_reference.start )
      {
        
        char base = m_alignment_output_pileup.reference_base_char_1(target_id, pos);
        if ( ((aligned_reference.truncate_start != 0) && (pos < aligned_reference.truncate_start))
            || ((aligned_reference.truncate_end != 0) && (pos > aligned_reference.truncate_end)) )
        {
          base = s_end_gap_character;
        }
        ref_add_left += base;
        pos++;
      }

      // create the string to add to the right side
      pos = aligned_reference.end + 1;
      while ( pos <= aligned_reference.end + ref_extend_right )
      {
        char base = m_alignment_output_pileup.reference_base_char_1(target_id, pos);
        if ( ((aligned_reference.truncate_start != 0) && (pos < aligned_reference.truncate_start))
            || ((aligned_reference.truncate_end != 0) && (pos > aligned_reference.truncate_end) ))
        {
          base = s_end_gap_character;
        }
        ref_add_right += base;
        pos++;
      }
        
      // update positions
      aligned_reference.start -= ref_extend_left;
      aligned_reference.end += ref_extend_right;
 
      // update qual and sequence strings
      aligned_reference.aligned_bases = ref_add_left + aligned_reference.aligned_bases + ref_add_right;
      aligned_reference.aligned_quals = repeat_char(char(k_reserved_quality_max), max_extend_left) + aligned_reference.aligned_quals + repeat_char(char(k_reserved_quality_max), max_extend_right);
  }

  // extend annotation line
  m_aligned_annotation.aligned_bases = repeat_char(' ', max_extend_left) + m_aligned_annotation.aligned_bases + repeat_char(' ', max_extend_right);

  // now go in and replace the empty space adjacent to each with the rest of the read sequence
  for ( Aligned_Reads::iterator itr_read = m_aligned_reads.begin(); itr_read != m_aligned_reads.end(); itr_read++ )
  {
    Aligned_Read& aligned_read ( ( *itr_read ).second );

    // *
    // *  FIND LENGTH OF THE NON-SPACE PART OF READ
    //      my $left_pos =  length($1);
    //      my $right_pos = $left_pos + length($2);
    //

    /* EX: counting padding in reads.
     *  string x = "  TTTA  ";
     *              01234567
     *   int right = x.find_last_not_of(' ');
     *   int left = x.find_first_not_of(' ');
     *> left: 2
     *> right: 5
     */

    uint32_t left_pos = aligned_read.aligned_bases.find_first_not_of ( ' ' );
    uint32_t right_pos = aligned_read.aligned_bases.find_last_not_of ( ' ' );

    if ( aligned_read.start > 1 )
    {
      string add_seq = aligned_read.read_sequence.substr ( 0, ( aligned_read.start - 1 ) );
      aligned_read.aligned_bases.replace ( ( left_pos - add_seq.length() ), add_seq.length(), to_lower ( add_seq ) ); //start,size, common.h function

      add_seq = repeat_char(char(k_reserved_quality_dont_highlight), aligned_read.start -1);
      aligned_read.aligned_quals.replace ( ( left_pos - add_seq.length() ), add_seq.length(), add_seq );
    }
  
    if ( aligned_read.end < aligned_read.length )
    {
      string add_seq = aligned_read.read_sequence.substr ( aligned_read.end, ( aligned_read.length - aligned_read.end ) );
      aligned_read.aligned_bases.replace ( right_pos+1, add_seq.length(), to_lower ( add_seq ) );

      add_seq = repeat_char(char(k_reserved_quality_dont_highlight), aligned_read.length - aligned_read.end);
      aligned_read.aligned_quals.replace ( right_pos+1,add_seq.length(), add_seq );
    }
  }

  //  ## swap out the ghost seq_ids and coordinates
  for ( uint32_t itr_ref = 0; itr_ref < m_aligned_references.size(); itr_ref++ )
  {    
    Aligned_Reference& aligned_reference = m_aligned_references[itr_ref];
    
    if (aligned_reference.ghost_seq_id != "")
    {
      aligned_reference.seq_id = aligned_reference.ghost_seq_id;
    }
    
    if (aligned_reference.truncate_start != 0)
    {
      uint32_t len = aligned_reference.end - aligned_reference.truncate_start + 1;
      aligned_reference.start = aligned_reference.ghost_start;
      aligned_reference.end = aligned_reference.start + aligned_reference.ghost_strand * (len - 1);
      aligned_reference.strand = aligned_reference.ghost_strand; // due to the way strands are given in junctions
    }
    
    if (aligned_reference.truncate_end != 0)
    {
      uint32_t len = aligned_reference.truncate_end - aligned_reference.start + 1;
      aligned_reference.end = aligned_reference.ghost_end;
      aligned_reference.start = aligned_reference.end + aligned_reference.ghost_strand * (len - 1);
      aligned_reference.strand = -aligned_reference.ghost_strand;
    }
    
  }

  // Need to reverse the coords for some
  for ( Aligned_Reads::iterator itr_read = m_aligned_reads.begin(); itr_read != m_aligned_reads.end(); itr_read++ )
  {
    Aligned_Read& aligned_read ( ( *itr_read ).second );

    if ( aligned_read.strand == -1 )
    {
        aligned_read.start = aligned_read.length - aligned_read.start + 1;
        aligned_read.end = aligned_read.length - aligned_read.end + 1;
    }
  }
} //End create alignment

string alignment_output::html_alignment( const string& region, output::cOutputEvidenceItem* output_evidence_item_ptr )
{

  // this sets object values (not re-usable currently)
  create_alignment(region, output_evidence_item_ptr);
  
  string output;
  
  if (m_error_message.size() > 0)
  {
    output += m_error_message;
    output += "<br>";
  }
  
  // We must have aligned reads for the remainder of the code to be safe
  if (m_aligned_reads.size() == 0)
    return output;
  
  /// all built by create_alignment and stored as member, m_ , variables.
  set_quality_range(m_quality_score_cutoff);

  Sorted_Keys sorted_keys;
  for (Aligned_Reads::iterator itr_read = m_aligned_reads.begin(); itr_read != m_aligned_reads.end(); itr_read++)
  {
    Sorted_Key sorted_key;
    sorted_key.seq_id = itr_read->first;
    sorted_key.aligned_bases = itr_read->second.aligned_bases;
    sorted_keys.push_back(sorted_key);
  }
  std::sort(sorted_keys.begin(),sorted_keys.end(),alignment_output::sort_by_aligned_bases_length);
  
  output += "\n<style>\n";
  output += html_header_string();
  output += "</style>\n";
  output += start_table("style=\"background-color: rgb(255,255,255)\"");
  output += start_tr();
  output += start_td("style=\"font-size:10pt\""); 
  
  output += "\n\n<!-- Reference Begin -->";  
  for (uint32_t index = 0; index < m_aligned_references.size(); index++)  {
    output += html_alignment_line( m_aligned_references[index], NULL, true, false) + "<BR>";  }
  output += "\n<!-- Reference End -->\n";
  
  output += html_alignment_line( m_aligned_annotation, NULL, false, false ) + "<BR>";

  output += "\n\n<!-- Reads Begin -->";
  for (Sorted_Keys::iterator itr_key = sorted_keys.begin(); itr_key != sorted_keys.end(); itr_key ++)  {
    output += html_alignment_line( m_aligned_reads[itr_key->seq_id], &m_aligned_references, true, true) + "<BR>";  }
  output += "\n<!-- Reads End -->\n";
  
  output += html_alignment_line( m_aligned_annotation, NULL, false, false)  + "<BR>"; 
  
  output += "\n\n<!-- Reference Begin -->";
  for (uint32_t index = 0; index < m_aligned_references.size(); index++)  {
    output += html_alignment_line( m_aligned_references[index], NULL, true, false) + "<BR>";  }
  output += "\n<!-- Reference End -->\n\n";
    
  output += end_td() + end_tr() + end_table();
  
  output += "\n\n<!-- Legend Begin -->";
  output += html_legend();
  output += "\n<!-- Legend End -->\n\n";
  
  return output;
}

string alignment_output::text_alignment( const string& region, cOutputEvidenceItem* output_evidence_item_ptr )
{
  
  // this sets object values (not re-usable currently)
  create_alignment(region, output_evidence_item_ptr);
  
  string output = "";
  
  if (m_error_message.size() > 0)
  {
    output += m_error_message;
    output += "\n";
  }
  
  // We must have aligned reads for the remainder of the code to be safe
  if (m_aligned_reads.size() == 0)
    return output;
  
  /// all built by create_alignment and stored as member, m_ , variables.
  set_quality_range(m_quality_score_cutoff);
  
  // @JEB sorting is not working yet
  Sorted_Keys sorted_keys;
  for (Aligned_Reads::iterator itr_read = m_aligned_reads.begin(); itr_read != m_aligned_reads.end(); itr_read++)
  {
    Sorted_Key sorted_key;
    sorted_key.seq_id = itr_read->first;
    sorted_key.aligned_bases = itr_read->second.aligned_bases;
    sorted_keys.push_back(sorted_key);
  }
  std::sort(sorted_keys.begin(),sorted_keys.end(),alignment_output::sort_by_aligned_bases_length);
  
  
  for (uint32_t index = 0; index < m_aligned_references.size(); index++)  {
    output += text_alignment_line( m_aligned_references[index] , true) + "\n";  }
  
  output += text_alignment_line( m_aligned_annotation, false ) + "\n";
  
  for (Sorted_Keys::iterator itr_key = sorted_keys.begin(); itr_key != sorted_keys.end(); itr_key ++)  {
    output += text_alignment_line( m_aligned_reads[itr_key->seq_id], true) + "\n";  }
  
  output += text_alignment_line( m_aligned_annotation, false)  + "\n"; 
  
  for (uint32_t index = 0; index < m_aligned_references.size(); index++)  {
    output += text_alignment_line( m_aligned_references[index], true) + "\n";  }
  
  output += "\n";
  
  return output;
}  
  

// 'quality_score_cutoff' = below this value you get a special color -- for bad Illumina bases
void alignment_output::set_quality_range(const uint32_t quality_score_cutoff)
{
  // calculate a cumulative distribution of the bases we are showing
  vector<uint32_t> qc(255, 0);
  uint32_t total = 0;
  int32_t max_qual = INT_MIN;

  for ( Aligned_Reads::iterator itr_read = m_aligned_reads.begin(); itr_read != m_aligned_reads.end(); itr_read++ )
  {
    Aligned_Read& aligned_read(itr_read->second );
    for ( uint32_t index = 0; index < aligned_read.qual_sequence.length(); index++ )
    {
      uint8_t c = static_cast<uint8_t>(aligned_read.qual_sequence[index]);
      if ( c >= quality_score_cutoff )
      {
        qc[c]++;
        total++;
      }
      if (c > max_qual) max_qual = c;
    }
  }

  // problem
  if(max_qual == INT_MIN) return;
  
  vector<uint8_t> qual_to_color(max_qual+1,0);
  double cutoff_percentiles[] = {0, 0.03, 0.1, 0.3, 0.9, 1.0};
  uint32_t num_cutoff_percentiles = 6;

  uint32_t current_cutoff_level = 0;
  //##set up to this score to the zero level (which is a completely different color)
  uint32_t index;
  for ( index = 0; index <  quality_score_cutoff; index++ )
  {
      qual_to_color[index] = current_cutoff_level;
  }

  current_cutoff_level++;
  double cumq = 0;
  while ( index < qual_to_color.size() )
  {
    cumq += (double)qc[index] / (double)total;

    // this can increment by at most one per quality score
    if ( cumq > cutoff_percentiles[current_cutoff_level] )
    {
      current_cutoff_level++;
    }
    qual_to_color[index] = current_cutoff_level;
    index++;     
    if (cumq == 1.0) break;
  }
  //#last must be set to max
  qual_to_color.back() = num_cutoff_percentiles - 1;
  //#first must be set to min
  qual_to_color[quality_score_cutoff] = 1;

  //#if there are at least as many quality scores in existence as
  //#there are color levels to assign....
  if ( index > ( num_cutoff_percentiles - 1 ) )
  {
    //  #...redistribute such that there are no jumps in quality level
    uint32_t gap = 1;
    while ( gap )
    {
      gap = 0;
      uint32_t last = 0;

      for ( uint32_t index = 0; index < qual_to_color.size(); index++ )
      {
        if ( qual_to_color[index] > last + 1 )
        {
          qual_to_color[index]--;
          gap = 1;
        }
        last = qual_to_color[index];
      }
    }
  }

  //
  //##finally, this sets the cutoff levels
  uint32_t last = 0;
  vector<uint8_t> cutoff_levels;
  cutoff_levels.clear();
  cutoff_levels.push_back ( quality_score_cutoff );

  for ( uint32_t index = quality_score_cutoff; index < qual_to_color.size(); index++ )
  {
    if ( qual_to_color[index] > last )
    {
        cutoff_levels.push_back ( index );
        last = qual_to_color[index];
    }
  }

  m_quality_range.qual_to_color_index = qual_to_color;

  // DEBUG
  //for (size_t i=0; i<qual_to_color.size(); i++)
  //{
  //  cout << i << ": " << (uint32_t)m_quality_range.qual_to_color_index[i] << endl;
  //}
  
  m_quality_range.qual_cutoffs = cutoff_levels;
}


string alignment_output::html_header_string()
{
  map<char, vector<string> > base_color_hash;
  
  // Colors used when base in read mismatches reference base
  
  // Black
  base_color_hash['G'].push_back ( "rgb(255,255,0)" );
  base_color_hash['G'].push_back ( "rgb(230,230,230)" );
  base_color_hash['G'].push_back ( "rgb(210,210,210)" );
  base_color_hash['G'].push_back ( "rgb(140,140,140)" );
  base_color_hash['G'].push_back ( "rgb(70,70,70)" );
  base_color_hash['G'].push_back ( "rgb(0,0,0)" );

  // Blue
  base_color_hash['C'].push_back ( "rgb(255,255,0)" );
  base_color_hash['C'].push_back ( "rgb(160,160,255)" );
  base_color_hash['C'].push_back ( "rgb(120,120,255)" );
  base_color_hash['C'].push_back ( "rgb(60,60,255)" );
  base_color_hash['C'].push_back ( "rgb(20,20,255)" );
  base_color_hash['C'].push_back ( "rgb(0,0,150)" );

  // Red
  base_color_hash['A'].push_back ( "rgb(255,255,0)" );
  base_color_hash['A'].push_back ( "rgb(255,210,210)" );
  base_color_hash['A'].push_back ( "rgb(255,180,180)" );
  base_color_hash['A'].push_back ( "rgb(255,100,100)" );
  base_color_hash['A'].push_back ( "rgb(255,20,20)" );
  base_color_hash['A'].push_back ( "rgb(200,0,0)" );

  // Green
  base_color_hash['T'].push_back ( "rgb(255,255,0)" );
  base_color_hash['T'].push_back ( "rgb(210,255,210)" );
  base_color_hash['T'].push_back ( "rgb(180,255,180)" );
  base_color_hash['T'].push_back ( "rgb(100,255,100)" );
  base_color_hash['T'].push_back ( "rgb(20,255,20)" );
  base_color_hash['T'].push_back ( "rgb(0,200,0)" );

  /*
  base_color_hash['N'].push_back ( "rgb(128,0,128)" );
  base_color_hash['N'].push_back ( "rgb(128,0,128)" );
  base_color_hash['N'].push_back ( "rgb(128,0,128)" );
  base_color_hash['N'].push_back ( "rgb(128,0,128)" );
  base_color_hash['N'].push_back ( "rgb(128,0,128)" );
  base_color_hash['N'].push_back ( "rgb(128,0,128)" );
   */
  
  string header_style_string;
  
  // General colors -- used to read names, etc.
  header_style_string += ".NC {color: rgb(0,0,0); background-color: rgb(255,255,255)}\n";         // normal color
  header_style_string += ".AM {color: rgb(0,0,0); background-color: rgb(210,210,210)}\n";         // ambiguously mapped
  header_style_string += ".AO {color: rgb(0,0,0); background-color: rgb(255,205,204)}\n";         // ambiguous overlap (junctions)
  header_style_string += ".AP {color: rgb(0,0,0); background-color: rgb(255,229,204)}\n";         // ambiguous overlap (junctions) - but only due to target site duplication

  
  // Colors used by characters in alignments
  header_style_string += ".UN {color: rgb(120,120,120); background-color: rgb(255,255,255)}\n";   // unaligned base
  
  header_style_string += ".GAP {color: rgb(255,255,255); background-color: rgb(128,0,128)}\n";   // unaligned base
  
  // Colors used by same as reference bases
  header_style_string += ".R0 {color: rgb(170,170,170); background-color: rgb(255,255,175)}\n";
  header_style_string += ".R1 {color: rgb(170,170,170); background-color: rgb(240,240,240)}\n";
  header_style_string += ".R2 {color: rgb(135,135,135); background-color: rgb(240,240,240)}\n";
  header_style_string += ".R3 {color: rgb(90,90,90); background-color: rgb(240,240,240)}\n";
  header_style_string += ".R4 {color: rgb(45,45,45); background-color: rgb(240,240,240)}\n";
  header_style_string += ".R5 {color: rgb(0,0,0); background-color: rgb(240,240,240)}\n";
  
  for ( map<char, vector<string> >::iterator itr = base_color_hash.begin(); itr != base_color_hash.end(); itr++ )
  {
    for ( uint32_t index = 0; index < base_color_hash[ ( *itr ).first].size(); index ++ )
    {
      if ( index > 0 )
      {
          header_style_string += "." + to_string ( ( *itr ).first) + to_string(index) + "{color: rgb(255,255,255); background-color:"
                                 + ( *itr ).second[index] + "}\n";
      }
      else
      {
          header_style_string += "." + to_string ( ( *itr ).first) + to_string(index) + "{color: rgb(120,120,120); background-color:"
                                 + ( *itr ).second[index] + "}\n";
      }
    }
  }

  header_style_string += "CODE {font-family:Monaco, monospace;}\n";
  
  // Below and default CODE font seem to be broken on Chrome
  //header_style_string += "CODE {font-family:\"Lucida Console\", Monaco, monospace;}\n";
  //header_style_string += "CODE {font-family:\"Courier New\", Courier, monospace;}\n";
  
  no_color_index = base_color_hash['G'].size()-1;
  return header_style_string;
}



string alignment_output::html_alignment_line(const Alignment_Base& a, Aligned_References* r, const bool coords, const bool use_quality_range)
{
  string output;
  output += "<CODE>";

  if (!a.aligned_quals.empty())
  {
    if (a.aligned_bases.length() != a.aligned_quals.length())
    {
      cout << a.aligned_bases << endl;
      cout << a.aligned_quals << endl;
      cerr << "unequal aligned base and aligned quals" << endl;
    }
  }   
  string last_color = "";

  for(uint32_t index = 0; index < a.aligned_bases.length(); index++)
  {
    uint8_t q = k_reserved_quality_max;
    if(!a.aligned_quals.empty())
    {
      q = static_cast<uint8_t>(a.aligned_quals[index]); 
    }
    
    char b=a.aligned_bases[index];       
    string color = "";
    
    // Special quality scores 
    // k_reserved_quality_junction_overlap = In junction overlap
    // k_reserved_quality_dont_highlight = Don't highlight, color in black on white
    // k_reserved_quality_max = Color with maximum quality
    if ((q == k_reserved_quality_dont_highlight) || (q == k_reserved_quality_junction_overlap)) 
    {
      color = "NC";
    }
    
    else if ((b == s_internal_gap_character) || (b == s_end_gap_character))
    {
      color = "NC";
      if (r != NULL) { 
        bool same_as_ref = true;
        for (Aligned_References::iterator it = r->begin(); it != r->end(); it++) {
          if ((it->aligned_bases[index] != s_internal_gap_character) && (it->aligned_bases[index] != s_end_gap_character))
            same_as_ref = false;
        }
        if (!same_as_ref) color = "GAP";
      }
    }
    
    else if (!use_quality_range)
    {
      color = to_upper(to_string(b)) + to_string<uint32_t>(no_color_index);
    }
    else // use_quality_range
    {
      uint8_t color_num;
      if (q == k_reserved_quality_max) {
        color_num = 5; // Hard coded max value...
      } else { 
        color_num = m_quality_range.qual_to_color_index[q];
      }
      char b_upper = toupper(b);
      color = b_upper;

      // Code to change things to grey if they are the same as the reference
      if (r != NULL) {
        
        // Get the consensus reference base
        Aligned_References::iterator it = r->begin();        
        char r_upper = toupper(it->aligned_bases[index]);
        it++;
        for (; it != r->end(); it++) {
          char this_r_upper = toupper(it->aligned_bases[index]);
          if ( (r_upper == s_internal_gap_character) && (this_r_upper != s_internal_gap_character) )
            r_upper = this_r_upper;
        }
        
        // If we are the same as the reference or an 'N' base, then change base and color scheme
        if ( (r_upper == b_upper) || (b_upper == 'N') ) { 
          color = "R";
          // if the upper case character does not equal the lowercase, then this base was trimmed.
          if (m_mask_ref_matches)
            b = (b == b_upper) ? s_reference_match_character : s_reference_match_masked_character;
        }
      }
      
      color += to_string<uint32_t>(color_num);
    }
      
    ASSERT(!color.empty(), "Alignment Output encountered unassigned base color.");

    if (color != last_color)
    {
      if (!last_color.empty())
        output += "</font>"; 
      output += "<font class=\"" + color + "\">"; 
      last_color = color;
    }
    
    if (b == ' ') output += "&nbsp;";
    else if (b == '-') output += "&#8209;";
    else output += b;
  }
  if (!last_color.empty()) output += "</font>";

  if (a.show_strand) output += "&nbsp;&nbsp;" + html_alignment_strand(a.strand); 

  // write the seq_id and coords in non-breaking html
  if (!a.seq_id.empty())
  {
    // substitute breaking hyphens
    string seq_id = substitute(a.seq_id, "–", "&#8209;"); // en-dash
    seq_id = substitute(seq_id, "-", "&#8209;");          // hyphen
    
    if(coords)
    {
      seq_id += "/" + to_string<uint32_t>(a.start);
      seq_id += "&#8209;";
      seq_id += to_string<uint32_t>(a.end);   
    }
    
    output += "&nbsp;&nbsp;";
    
    output += "<font class=\"" + a.read_name_style + "\">";
    output += seq_id;
    output += "</font>"; 
  }   
  
  // write the mapping quality - it would be nice to justify this
  if (!m_is_junction_junction && (a.mapping_quality != -1)) {
    output += "&nbsp;(MQ=" + to_string(a.mapping_quality) + ")";
  }
  
  
  output += "</CODE>";
  
  return output;
}
  
string alignment_output::html_alignment_strand(const int8_t &strand)
{
  if (strand == -1) return "&lt;";
  if (strand == 1) return "&gt;";
    return ".";
}
  
string alignment_output::html_legend()
{
  string output; 
  output += "<br>";
  output += start_table("border=\"0\" cellspacing=\"1\" cellpadding=\"3\" style=\"font-size:10pt\""); 
  output += tr(th(ALIGN_LEFT, "Alignment Legend"));
  
  // create legend information
  output += start_tr();
  output += start_td();

  if (m_mask_ref_matches) {
    output += nonbreaking("Aligned base (shaded by quality score): ");
  } else {
    output += nonbreaking("Aligned base mismatch/match (shaded by quality score): ");
  }
  Alignment_Base  temp_a;
  Aligned_Reference temp_r;
  vector<Aligned_Reference> aligned_r;
  if (m_mask_ref_matches) {
    temp_r.aligned_bases = "TAGCGG"; // Different from below so colors are shown, except last bases
    temp_a.aligned_bases = "ATCGGg";
    temp_a.aligned_quals = repeat_char('\0', 6);
  } else {
    temp_r.aligned_bases = "TAGC/ATCG"; // Different from below so colors are shown, except last bases
    temp_a.aligned_bases = "ATCG/ATCG";
    temp_a.aligned_quals = repeat_char('\0', 4) + char(k_reserved_quality_dont_highlight) + repeat_char('\0', 4);
  }
  aligned_r.push_back(temp_r);
  temp_a.show_strand = false;
  output += html_alignment_line(temp_a, &aligned_r, false, true);
  
  for (uint8_t index = 1; index < m_quality_range.qual_cutoffs.size(); index++)
  {
    char c = m_quality_range.qual_cutoffs[index];   
    output += "<CODE>&nbsp;&lt;&nbsp;" + to_string<uint32_t>(c) + "&nbsp;&le;&nbsp;</CODE>";
    
    if (m_mask_ref_matches) {
      temp_a.aligned_quals = repeat_char(static_cast<char>(c), 6);
    } else {
      temp_a.aligned_quals = repeat_char(static_cast<char>(c), 4) + char(k_reserved_quality_dont_highlight) + repeat_char(static_cast<char>(c), 4);
    }
    
    temp_a.show_strand = false;
    
    output += html_alignment_line(temp_a, &aligned_r, false, true);
  }
  
  output += end_td();
  output += end_tr();
  
  output += start_tr();
  output += start_td();

  if (m_mask_ref_matches) {
    output += nonbreaking("Unaligned base: <code>atcg</code>");
    output += nonbreaking("    Matching base: <code>") + s_reference_match_character + "</code>";
    output += nonbreaking("    Masked matching base: <code>") + s_reference_match_masked_character + "</code>";
  } else {
    output += nonbreaking("Unaligned base: <code>atcg</code>");
    temp_a.aligned_bases = "atcg";
    temp_a.aligned_quals = repeat_char(static_cast<char>(k_reserved_quality_max), 4);
    aligned_r[0].aligned_bases = "ATCG"; 
    output += nonbreaking("    Masked matching base: ") + html_alignment_line(temp_a, &aligned_r, false, true);
  }
    
  aligned_r[0].aligned_bases = "-"; // Different from below so colors are shown, except last base
  temp_a.aligned_bases = "-";
  temp_a.aligned_quals = repeat_char('\0', 1);
  output += nonbreaking("    Alignment gap: ") + html_alignment_line(temp_a, &aligned_r, false, true);

  
  aligned_r[0].aligned_bases = "G"; // Different from below so colors are shown, except last base
  temp_a.aligned_bases = "-";
  temp_a.aligned_quals = repeat_char('\0', 1);
  output += nonbreaking("    Deleted base: ") + html_alignment_line(temp_a, &aligned_r, false, true);
  
  output += end_td();
  output += end_tr();
  
  if (m_show_ambiguously_mapped) {
    output += tr(td("<code><font class=\"AM\">&nbsp;read_name</font></code> Reads with multiple best matches in reference."));
    
  }
  
  if (m_is_junction) {
    
    output += tr(th(ALIGN_LEFT, "Reads not counted as support for junction"));
    output += tr(td("<font class=\"AO\">read_name</font></code> Not counted due to insufficient overlap past the breakpoint."));
    output += tr(td("<font class=\"AP\">read_name</font></code> Not counted due to not crossing MOB target site duplication."));
  }

  output += end_table();

  
  return output;
}
  
string alignment_output::text_alignment_line(const alignment_output::Alignment_Base& a, const bool coords)
{
  string output;
    
  for(uint32_t index = 0; index < a.aligned_bases.length(); index++)
  {
    string b=to_string(a.aligned_bases[index]);       
    output += b;
  }
  
  if (a.show_strand) output += "  " + text_alignment_strand(a.strand); 
  
  // write the seq_id and coords in non-breaking html
  if (!a.seq_id.empty())
  {
    string seq_id = a.seq_id;
    // substitute breaking hyphens
    if(coords)
    {
      seq_id += "/" + to_string<uint32_t>(a.start);
      seq_id += "-";
      seq_id += to_string<uint32_t>(a.end);   
    }
    
    output += "  ";
    output += seq_id;
  }   
  
  // write the mapping quality - it would be nice to justify this
  if (a.mapping_quality != -1) {
    output += " (MQ=" + to_string<>(a.mapping_quality) + ")";
  }
  
  return output;
}
  
string alignment_output::text_alignment_strand(const int8_t &strand)
{
  if (strand == -1) return "<";
  if (strand == 1) return ">";
  return ".";
}

} // namespace breseq
