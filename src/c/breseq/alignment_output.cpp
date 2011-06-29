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
 *    Copyright (c) 2011 The University of Texas at Austin
 *
 *    breseq is free software; you can redistribute it and/or modify it under the
 *    terms the GNU General Public License as published by the Free Software
 *    Foundation; either version 1, or (at your option) any later version.
 *
 *****************************************************************************/

#include "breseq/alignment_output.h"

using namespace std;

namespace breseq
{

bool verbose = true; //TODO Options
bool text = false; //TODO Options

alignment_output::Alignment_Output_Pileup::Alignment_Output_Pileup ( const string& bam, const string& fasta, const uint32_t maximum_to_align )
        : pileup_base ( bam, fasta )
        , unique_start ( 0 )
        , unique_end ( 0 )
        , total_reads ( 0 )
        , processed_reads ( 0 )
        , maximum_to_align ( maximum_to_align )
        , insert_start ( 0 ) //these should be initialized by create alignment
        , insert_end ( 0 ) //these should be initialized by create alignment
        , last_pos ( 0 )
{
}

alignment_output::Alignment_Output_Pileup::~Alignment_Output_Pileup() {}

alignment_output::alignment_output ( string bam, string fasta, uint32_t in_maximum_to_align, const uint32_t quality_score_cutoff )
        : m_alignment_output_pileup ( bam, fasta, in_maximum_to_align )
        ,m_aligned_reads ( m_alignment_output_pileup.aligned_reads )
{
}

void alignment_output::create_alignment ( const string& region )
{
  //  my $reference_length = $bam->length($seq_id);
  //  my $aligned_reads;

  //
  //  my @aligned_references = ( {} );
  //  ////## More complex information about multiple reference sequence lines was provided
  //  ////if (defined $options->{alignment_reference_info_list})
  //  ////{
  //  ////    @aligned_references = @{$options->{alignment_reference_info_list}};
  //  ////}
  //  ////foreach my $aligned_reference (@aligned_references)
  //  ////{
  //  ////    $aligned_reference->{seq_id} = $seq_id;
  //  ////    $aligned_reference->{strand} = 0;
  //  ////}
  
  // we need the target_id to properly fill out the reference sequence later
  uint32_t target_id, start_pos, end_pos;
  m_alignment_output_pileup.parse_region(region, target_id, start_pos, end_pos);
    
  
  ///Single Reference case
  Aligned_Reference aligned_reference;
  aligned_reference.target_id = target_id;
  aligned_reference.seq_id = m_alignment_output_pileup.target_name(target_id);
  m_alignment_output_pileup.aligned_references.push_back(aligned_reference);

  // @JEB TODO
  // Need to implement case where we are given a reference sequence list

  
  // Use fetch to determine the list of reads that will be aligned to this position
  m_alignment_output_pileup.do_fetch ( region );


  if ( ( m_alignment_output_pileup.unique_start == 0 ) || ( m_alignment_output_pileup.unique_end == 0 ) )
  {
    cout << "No unique start or end initialized" << endl;
    return;
  }

  if ( m_alignment_output_pileup.total_reads > m_alignment_output_pileup.maximum_to_align )
  {
    cout << "Reads exceeded maximum to display alignment. ";
    cout << m_alignment_output_pileup.total_reads << " reads. ";
    cout << "(Limit = " << m_alignment_output_pileup.maximum_to_align << ")" << endl;
    return;
  }

  // @JEB TODO
  // sampling a limited number of reads when there are too many needs to be implemented! (BELOW)
  //  my $message;
  //////  if ($self->{maximum_to_align} && ($total_reads > $self->{maximum_to_align}))
  //////  {
  //////      $message = "Only $self->{maximum_to_align} of $total_reads total aligned reads displayed.";
  //////      my $new_aligned_reads;
  //////      my @new_keys = shuffle(keys %$aligned_reads);
  //////      foreach (my $i=0; $i<$self->{maximum_to_align}; $i++)
  //////      {
  //////          $new_aligned_reads->{$new_keys[$i]} = $aligned_reads->{$new_keys[$i]};
  //////      }
  //
  //////      $aligned_reads = $new_aligned_reads;
  //////  }
  //

  // * Call
  // * do_pileup ();

  m_alignment_output_pileup.last_pos = 0;
  m_alignment_output_pileup.do_pileup ( region );

  //do_pileup populates the following
  m_aligned_reads = m_alignment_output_pileup.aligned_reads;
  m_aligned_references = m_alignment_output_pileup.aligned_references;
  m_aligned_annotation = m_alignment_output_pileup.aligned_annotation;

  // If nothing aligned return @JEB TODO: we need to return something
  // that indicates to the caller to write "no aligned reads"
  if ( m_aligned_reads.empty() )
  {
    cerr << "pileup function did not work as intended";
    return;
  }
 
  //  ##now add the unaligned portions of each

  int32_t max_extend_left = 0;
  int32_t max_extend_right = 0;
  for ( Aligned_Reads::iterator itr_read = m_aligned_reads.begin(); //TODO typedef map
          itr_read != m_aligned_reads.end(); itr_read++ )
    
  {
    Aligned_Read& aligned_read = ( m_aligned_reads[ ( *itr_read ).first] );

    // read the padding off of each sequence
    int32_t left_padding_length = aligned_read.aligned_bases.find_first_not_of ( ' ' );
    int32_t right_padding_length = aligned_read.aligned_bases.length() - aligned_read.aligned_bases.find_last_not_of ( ' ' );
         
    if ( verbose )
    {
      cout  << aligned_read.start << " " << aligned_read.end << " " << aligned_read.seq_id << endl;
    }
    
    int32_t extend_left = ( (int32_t)aligned_read.start - 1) - left_padding_length;
    int32_t extend_right = ((int32_t)aligned_read.length - (int32_t)aligned_read.end ) - right_padding_length;
          
    if ( extend_left > max_extend_left )
      max_extend_left = extend_left;
    if ( extend_right > max_extend_right )
      max_extend_right = extend_right;
  }
  
  // Now add this much to every read 
  for (Aligned_Reads::iterator itr_read = m_aligned_reads.begin(); //TODO typedef map
          itr_read != m_aligned_reads.end(); itr_read++ )
  {
    Aligned_Read& aligned_read ( itr_read->second );
    
    aligned_read.aligned_bases = repeat_char(' ', max_extend_left) + aligned_read.aligned_bases + repeat_char(' ', max_extend_right);
    aligned_read.aligned_quals = repeat_char(char(255), max_extend_left) + aligned_read.aligned_bases + repeat_char(char(255), max_extend_right);
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
          ref_add_left += repeat_char('-', max_extend_left - aligned_reference.start + 1);
      }
    
      // handle attempts to extend past end of reference sequence
      if ( ( aligned_reference.end + max_extend_right ) > aligned_reference.reference_length )
      {
          ref_extend_right = aligned_reference.reference_length - aligned_reference.end;
          ref_add_right += repeat_char('-', aligned_reference.end + max_extend_right - aligned_reference.reference_length);
      }

      uint32_t pos = aligned_reference.start - ref_extend_left;
        
      // create the string to add to the right side
      while ( pos < aligned_reference.start )
      {
        
        char base = m_alignment_output_pileup.reference_base_char_1(target_id, pos);
        if ( (aligned_reference.truncate_start != 0) && (pos < aligned_reference.truncate_start) 
            || (aligned_reference.truncate_end != 0) && (pos > aligned_reference.truncate_end) )
        {
          base = '.';
        }
        ref_add_left += base;
        pos++;
      }

      // create the string to add to the right side
      pos = aligned_reference.end + 1;
      while ( pos <= aligned_reference.end + ref_extend_right )
      {
        char base = m_alignment_output_pileup.reference_base_char_1(target_id, pos);
        if ( (aligned_reference.truncate_start != 0) && (pos < aligned_reference.truncate_start) 
          || (aligned_reference.truncate_end != 0) && (pos > aligned_reference.truncate_end) )
        {
          base = '.';
        }
        ref_add_right += base;
        pos++;
      }
    
      for (uint32_t i= aligned_reference.start-ref_extend_left; i<=aligned_reference.start-1; i++)
        ref_add_left += m_alignment_output_pileup.reference_base_char_1(target_id, i);
    
      for (uint32_t i= aligned_reference.end+1; i<=aligned_reference.end + ref_extend_right; i++)
        ref_add_right += m_alignment_output_pileup.reference_base_char_1(target_id, i);
    
      // update positions
      aligned_reference.start -= ref_extend_left;
      aligned_reference.end += ref_extend_right;
 
      // update qual and sequence strings
      aligned_reference.aligned_bases = ref_add_left + aligned_reference.aligned_bases + ref_add_right;
      aligned_reference.aligned_quals = repeat_char(char(255), max_extend_left) + aligned_reference.aligned_quals + repeat_char(char(255), max_extend_right);
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

      add_seq = repeat_char(char(254), aligned_read.start -1);
      aligned_read.aligned_quals.replace ( ( left_pos - add_seq.length() ), add_seq.length(), add_seq );
    }
  
    if ( aligned_read.end < aligned_read.length )
    {
      string add_seq = aligned_read.read_sequence.substr ( aligned_read.end, ( aligned_read.length - aligned_read.end ) );
      aligned_read.aligned_bases.replace ( right_pos, add_seq.length(), to_lower ( add_seq ) );

      add_seq = repeat_char(char(254), aligned_read.length - aligned_read.end);
      aligned_read.aligned_quals.replace ( right_pos,add_seq.length(), add_seq );
    }
  }

  // @JEB TODO
  //  ## swap out the ghost seq_ids
  //////  foreach my $ar (@aligned_references)
  //////  {
  //      $ar->{seq_id} = $ar->{ghost_seq_id} if (defined $ar->{ghost_seq_id});
  //
  //      if (defined $ar->{truncate_start})
  //      {
  //          $ar->{start} = $ar->{truncate_start};
  //          my $len = $ar->{end} - $ar->{start} + 1;
  //          $ar->{start} = $ar->{ghost_start};
  //          $ar->{end} = $ar->{start} + $ar->{ghost_strand} * ($len - 1);
  //      }
  //      if (defined $ar->{truncate_end})
  //      {
  //          $ar->{end} = $ar->{truncate_end};
  //          my $len = $ar->{end} - $ar->{start} + 1;
  //          $ar->{end} = $ar->{ghost_end};
  //          $ar->{start} = $ar->{end} + $ar->{ghost_strand} * ($len - 1);
  //      }
  //////  }


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

string alignment_output::html_alignment ( const string& region )
{

  // this sets object values (not re-usable currently)
  create_alignment(region);
  
  string output = "";
  //return p . "No reads uniquely align to region." if (!defined $alignment_info);
  //$output .= p . "$alignment_info->{message}" if ($alignment_info->{message});
  //return $output if (!defined $alignment_info->{aligned_reads});
  //
  //my $aligned_reads = $alignment_info->{aligned_reads};
  //my @aligned_references = @{$alignment_info->{aligned_references}};
  //my $aligned_annotation = $alignment_info->{aligned_annotation};
  //my $quality_range = $self->set_quality_range($aligned_reads, $options);
  //
  /// all built by create_alignment and stored as member, m_ , variables.
    
  set_quality_range();

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
  
  
  output += "<style>"+create_header_string()+"</style>";
  output += "<table style=\"background-color: rgb(255,255,255)\">";
  output += "<tr><td style=\"font-size:10pt\">";    
  
  for (uint index = 0; index < m_aligned_references.size(); index++)
  {
    output += html_alignment_line( m_aligned_references[index] , true ,false) + "<BR>";
  }
  output += html_alignment_line( m_aligned_annotation, false, false ) + "<BR>";

  for (Sorted_Keys::iterator itr_key = sorted_keys.begin(); itr_key != sorted_keys.end(); itr_key ++)
  {
    output += html_alignment_line( m_aligned_reads[itr_key->seq_id], true, true) + "<BR>";
  }
  output += html_alignment_line( m_aligned_annotation, false, false)  + "<BR>"; 
  
  for (uint index = 0; index < m_aligned_references.size(); index++)
  {
    output += html_alignment_line( m_aligned_references[index], true, false) + "<BR>";
  }
  output += "<BR>";

  // create legend information
  output += "<CODE>Base quality scores:&nbsp</CODE>";
  Alignment_Base temp_a;
  temp_a.aligned_bases = "ATCG";
  temp_a.aligned_quals = repeat_char('\0', 4);
  temp_a.show_strand = false;

  output += html_alignment_line(temp_a, false, true);
  for (uint8_t index = 1; index < m_quality_range.qual_cutoffs.size(); index++)
  {
    char c = m_quality_range.qual_cutoffs[index];   
    output += "<CODE>&nbsp;&lt;&nbsp;" + to_string<uint32_t>(c) + "&nbsp;&le;&nbsp;</CODE>";

    temp_a.aligned_bases = "ATCG";
    temp_a.aligned_quals = repeat_char(static_cast<char>(c), 4);
    temp_a.show_strand = false;

    output += html_alignment_line(temp_a, false, true);
  }
  
  output += "</TABLE></TR></TD>";
  return output;
}

/*! Called for each position.*/
void alignment_output::Alignment_Output_Pileup::pileup_callback ( const pileup& p )
{
  
  //  #create the alignment via "pileup"
  //  my $pileup_function = sub { 
    ///EX: region = REL606-5:1-15; 1 is start, 15 is end
    uint32_t& start_1( m_start_position_1 );
    uint32_t& end_1( m_end_position_1 );
    
    
    // my ($seq_id,$pos,$pileup) = @_;
    uint32_t pos_1 = p.position_1();
    
    // print "POSITION: $pos\n" if ($verbose);
    if ( verbose ) //TODO Options
        cout << "POSITION: " << pos_1 << endl;
    
    // return if ($pos < $unique_start);
        // return if ($pos > $unique_end);
        if ( ( pos_1 < unique_start ) || ( pos_1 > unique_end ) )
          return;
        //
        // foreach my $aligned_reference (@aligned_references)
        // {
        //   $aligned_reference->{start} = $pos if (!defined $aligned_reference->{start});
        //   $aligned_reference->{end} = $pos;
        // }
        //
        
    ///For Multiple Reference case
//     for ( uint32_t index_aligned_reference = 0;
//          index_aligned_reference < aligned_references.size();
//          index_aligned_reference++ )
//          {
//            //@JEB: 0 counts as undefined since we use 1-indexed coords
//            if (aligned_references[index_aligned_reference].start == 0) 
//            {
//              aligned_references[index_aligned_reference].start = pos_1;  
//            }
//              aligned_references[index_aligned_reference].end =pos_1;
//          }
//TODO Must be better way to get this information then calling it in pileup...
        
   ///Single Reference case
   Aligned_Reference& aligned_reference (aligned_references[0]);
   aligned_reference.reference_name = p.target_name();
   aligned_reference.reference_length = p.target_length();
   if(aligned_reference.start == 0)
     aligned_reference.start = pos_1;
   aligned_reference.end = pos_1;
   
   //TEST with multiple reference sequences

    ////      ## Cull the list to those we want to align
    ////      ## Other alignments may overlap this position that DO NOT
    ////      ## overlap the positions of interest. This removes them.
    ////      @$pileup = grep { defined $aligned_reads->{$_->alignment->display_name} } @$pileup;
    //
    //      ## Find the maximum indel count
    //      ## We will add gaps to reads with an indel count lower than this.
    //      my $max_indel = 0;
    //      my $alignment_spans_position;
    int temp_max_indel = 0;
    vector<string> alignment_spans_position;

    // ALIGNMENT: for my $p (@$pileup)
    for ( pileup::const_iterator itr_pileup = p.begin();
            itr_pileup != p.end() ; itr_pileup ++ )
    {
        //## alignment spans the position unless this is the first base...
        //$alignment_spans_position->{$p->alignment->display_name} = 1 if ($p->qpos > 0);
        ///KNOWN qpos is zero indexed
        if (itr_pileup->query_position_0() > 0)
        {
            alignment_spans_position.push_back ( itr_pileup->read_name() );
        }
          
        //$max_indel = $p->indel if ($p->indel > $max_indel);
        if (itr_pileup->indel() > temp_max_indel )
        {
          temp_max_indel = itr_pileup->indel();
        }
    }
    assert(temp_max_indel >= 0);
    uint32_t max_indel = static_cast<uint32_t>(temp_max_indel);
  
    if ( verbose ) //TODO Options
        cout << "MAX INDEL: " << max_indel << endl;

    // ## Reference only positions, with no aligned reads
    // ## are never called, so we keep track of the last position to add them
    //     if (defined $last_pos && ($last_pos < $pos - 1))
    //     {
    //         $last_pos++;
    //         while ($last_pos < $pos)
    //         {
    // ## READS: add gaps to all
    //             foreach my $key (keys %$aligned_reads)
    //             {
    //                 my $aligned_read = $aligned_reads-> {$key};
    //                 $aligned_read-> {aligned_bases} .= ($alignment_spans_position-> {$key}) ? '.' : ' ';
    //                 $aligned_read-> {aligned_quals} .= chr(255);
    //             }
    //
    // ## REFERENCE SEQUENCES: add actual bases
    //             my $ref_bases = $bam->segment($seq_id,$last_pos,$last_pos)->dna;
    //             foreach my $aligned_reference (@aligned_references)
    //             {
    //                 my $my_ref_bases = $ref_bases;
    //                 ////                 $my_ref_bases = '.' if ($aligned_reference->{truncate_start} && ($last_pos < $aligned_reference->{truncate_start}))
    //                 ////                                  || ($aligned_reference->{truncate_end} && ($last_pos > $aligned_reference->{truncate_end}));
    //                 $aligned_reference-> {aligned_bases} .= $my_ref_bases;
    //                 $aligned_reference-> {aligned_quals} .= chr(255);
    //             }
    //
    // ## ANNOTATIONS: add gaps
    //             $aligned_annotation-> {aligned_bases} .=
    //                 ( (($insert_start == 0) && ($last_pos == $start)) || (($last_pos > $start) && ($last_pos <= $end)) )
    //                 ? '|' : ' ';
    //
    //             $last_pos++;
    //         }
    //     }
    //     $last_pos = $pos;


    //     //NOT NEEDED IN C++ PORT?
    //     if ( ( last_pos != 0 ) && ( last_pos < pos ) )
    //     {
    //         last_pos++;
    //         while ( last_pos < p.position_1() )
    //         {
    //             //READS: add gaps to all
    //             for ( Aligned_Reads::iterator itr_reads = aligned_reads.begin(); 
    //                     itr_reads != aligned_reads.end(); itr_reads++ )
    //             {
    //                 aligned_reads[ ( *itr_reads ).first].aligned_bases += '.';
    //                 aligned_reads[itr_reads->first].aligned_quals += char ( 255 );
    //             }
    //             //REFERENCE SEQUENCES: add actual bases
    //             for ( uint32_t index_aligned_reference;
    //                     index_aligned_reference < aligned_references.size();
    //                     index_aligned_reference++ )
    //             {
    //                 aligned_references[index_aligned_reference].aligned_bases += p.reference_base_char_1 ( pos );
    //                 aligned_references[index_aligned_reference].aligned_quals += char ( 255 );
    //             }//TEST with multiple reference sequences
    //             //ANNOTATIONS: add gaps
    //
    //
    //
    //             if ( ( insert_start == 0 ) && ( last_pos == start )
    //                     || ( last_pos > start ) && ( last_pos <= end ) )
    //             {
    //                 aligned_annotation.aligned_bases += '|';
    //             }
    //             else
    //             {
    //                 aligned_annotation.aligned_bases += ' ';
    //             }
    //             last_pos++;
    //         }
    //     }
    //     last_pos = pos;

    // ## END adding reference only positions.
    //
    // ## Now add this position to the alignments
    //     my $updated;
		/// Reset "updated flag" @JEB
    /// This flag is set if the read was among those returned by the pileup
    /// if it is not set, we insert padding after the loop through the alignments!
    for (Aligned_Reads::iterator itr_read = aligned_reads.begin(); itr_read != aligned_reads.end(); itr_read++ )
    {
        Aligned_Read& aligned_read ( itr_read->second );
        aligned_read.updated = false;
    }
    // ALIGNMENT:
    //     for my $p (@$pileup)
    //     {
    //BEGIN FIRST ALIGNMENT LABEL LOOP
    for ( pileup::const_iterator itr_pileup = p.begin(); itr_pileup != p.end(); itr_pileup++ )
    {
        // my $a = $p->alignment;
        Aligned_Read& aligned_read = aligned_reads[itr_pileup->read_name()];
      
        //$updated-> {$a->display_name} = 1;
        aligned_read.updated = true;
        //
        // ##this setup gives expected behavior for indels!
        //my $indel = $p->indel;
        // $indel = 0 if ($indel < 0);
        // $indel = -1 if ($p->is_del);
        bool indel = itr_pileup->indel();

        // ## Which read are we on?
      
        // my $aligned_read = $aligned_reads-> {$a->display_name};
        // $aligned_read-> {strand} = ($a->reversed) ? -1 : +1;
        aligned_read.strand = itr_pileup->strand();
      
        // $aligned_read-> {reference_start} = $pos if (!defined $aligned_read-> {reference_start});
        if (aligned_read.reference_start == 0) aligned_read.reference_start = pos_1;
      
        // $aligned_read-> {reference_end} = $pos;
        aligned_read.reference_end = pos_1;
      
        // $aligned_read-> {start} = $p->qpos+1 if (!defined $aligned_read-> {start});
        if ( aligned_read.start == 0 )
        {
          aligned_read.start = itr_pileup->query_position_1();
        }
          
        //         $aligned_read-> {end} = $p->qpos+1;
				aligned_read.end = itr_pileup->query_position_1(); 

        // ## READS: add aligned positions
        // for (my $i=0; $i<=$max_indel; $i++)
				for ( uint32_t index = 0; index <= max_indel; index ++ )
        {
          // if ($i > $indel)
					if ( index > indel )
					{
            //print $p->indel . "\n" if ($verbose);
        		if ( verbose )
          	{
              cout << indel << endl;
            }
            //$aligned_read-> {aligned_bases} .= '.';
            aligned_read.aligned_bases += ',';
            
            // $aligned_read-> {aligned_quals} .= chr(255);
            aligned_read.aligned_quals += char ( 255 );
					}
          // else
					else
          {
            // my $quality = $a->qscore->[$p->qpos+$i];
        		uint8_t quality = itr_pileup->read_base_quality_0( itr_pileup->query_position_0() +index );
            // my $base  = substr($a->qseq, $p->qpos+$i,1);
        		char base = itr_pileup->read_char_sequence().substr ( ( *itr_pileup ).query_position_0() +index,1 ) [0];

            // if (!$options-> {text})
						if ( !text ) //TODO Options
            {
              // my $trim_left = $a->aux_get('XL');
        			uint8_t trim_left = ( *itr_pileup ).trim_left();
              // $base = "\L$base" if ((defined $trim_left) && ($p->qpos+1 <= $trim_left));
							if ( ( trim_left > 0 ) && ( itr_pileup->query_position_1() <= trim_left ) ) 
              {
                base = tolower ( base );           
              }
              
              // my $trim_right = $a->aux_get('XR');
						  uint8_t trim_right = ( *itr_pileup ).trim_right();
              // $base = "\L$base" if ((defined $trim_right) && ($a->l_qseq-$p->qpos <= $trim_right));
							if ( ( trim_right !=0 ) && ( ( itr_pileup->read_length() - itr_pileup->query_position_0() ) <= trim_right ) )
              {
              	base = tolower ( base );
              }
						}
            // $aligned_read-> {aligned_bases} .= $base;
        	  aligned_read.aligned_bases += base;
            // $aligned_read-> {aligned_quals} .= chr($quality);
						aligned_read.aligned_quals += char ( quality );
					}
          // print $aligned_read-> {aligned_bases} . " " .$aligned_read-> {seq_id} . " " . $p->indel . "\n" if ($verbose);
          if ( verbose )
          {
          	cout << aligned_read.aligned_bases << " ";
            cout << aligned_read.seq_id << " ";
            cout << ( *itr_pileup ).indel() << endl;
          }

        }
    }
    //END FIRST ALGINMENT LABEL LOOP
   
    //
    // ## READS: handle those with no
    //     foreach my $key (keys %$aligned_reads)
    //     {
		for ( Aligned_Reads::iterator itr_read = aligned_reads.begin(); 
            itr_read != aligned_reads.end(); itr_read++ )
    {    
      // my $aligned_read = $aligned_reads-> {$key};
			Aligned_Read& aligned_read (itr_read->second );
      
      // next if ($updated-> {$key});
    	if ( aligned_read.updated ) continue;
      
      // $aligned_read-> {aligned_bases} .= ' ' x ($max_indel+1);
			for(uint32_t i = 0; i < (max_indel + 1); i++)
      {
				aligned_read.aligned_bases += ' '; //BUG dont .append(' ', max_indel+1)
      }
      
      // $aligned_read-> {aligned_quals} .= chr(254) x ($max_indel+1);
    	for(uint32_t i = 0; i < (max_indel + 1); i++)
      {
				aligned_read.aligned_quals += char(254); //BUG dont .append(char(254), max_indel+1)
      }
      // print $aligned_read-> {aligned_bases} . " NOALIGN " . $key . " " . "\n" if ($verbose);
     	if ( verbose )
        cout << aligned_read.aligned_bases << " NOALIGN " << itr_read->first << endl;
    }
    //
    // ##now handle the reference sequence
    // my $ref_base = $bam->segment($seq_id,$pos,$pos)->dna;
    char ref_base = p.reference_base_char_1(pos_1);

    // foreach my $aligned_reference (@aligned_references)
		for ( uint32_t index = 0; index < aligned_references.size(); index ++ )
    {
      // my $my_ref_base = $ref_base;
    	char my_ref_base = ref_base;
      //// $my_ref_base = '.' if ($aligned_reference-> {truncate_start} && ($last_pos < $aligned_reference-> {truncate_start}))
      //// || ($aligned_reference-> {truncate_end} && ($last_pos > $aligned_reference-> {truncate_end}));
      
			Aligned_Reference& aligned_reference ( aligned_references[index] );	
      
      // $aligned_reference-> {aligned_bases} .= $my_ref_base;
			aligned_reference.aligned_bases += my_ref_base;
      
      // $aligned_reference-> {aligned_bases} .= '.' x ($max_indel) if ($max_indel > 0);
			if ( max_indel >  0 )
      	aligned_reference.aligned_bases.append ( '.',max_indel );
        // $aligned_reference-> {aligned_quals} .= chr(255) x ($max_indel+1);
      for(uint32_t i = 0 ; i < (max_indel + 1); i++)
				aligned_reference.aligned_quals += char(255); //BUG dont .append ( char ( 255 ), max_indel +1 );
   
				///TODO assigned here for single reference case?
				///aligned_reference.base = my_ref_base;
        ///aligned_reference.reference_length = p.target_length();
        ///aligned_reference.reference_name = p.target_name();
		}
    //
  
    // @JEB these 
    ///For basic case these are 0, defined in the user defined Region
  
    // ##also update any positions of interest for gaps
  
    // for (my $insert_count=0; $insert_count<=$max_indel; $insert_count++)
	  for ( uint32_t insert_count = 0; insert_count <= max_indel; insert_count++ )
    {
      // if (   (($insert_start <= $insert_count) && ($insert_count <= $insert_end) && ($pos == $start) && ($pos == $end))
    	if ( (( insert_start <= insert_count ) && ( insert_count <= insert_end ) && ( pos_1 == start_1 ) && ( pos_1 == end_1 ) )
          // || (($insert_start <= $insert_count) && ($pos == $start) && ($pos != $end))
					|| ( ( insert_start <= insert_count ) && ( pos_1 == start_1 ) && ( pos_1 != end_1 ) )
    	// || (($pos < $end) && ($pos > $start))
			 		|| ( ( pos_1 < end_1 ) && ( pos_1 > start_1 ) )
    	// || (($insert_end <= $insert_count) && ($pos == $end) && ($pos != $start)) )
					|| ( ( insert_end <= insert_count ) && ( pos_1 == end_1 ) && ( pos_1 != start_1 ) ) )
			{
        
        // $aligned_annotation-> {aligned_bases} .= '|';
    		aligned_annotation.aligned_bases += '|';
			}
    	// else
	  	else
			{		
        //$aligned_annotation-> {aligned_bases} .= ' ';
    		aligned_annotation.aligned_bases += ' ';
     	}
    }

}
//END create_alignment()

/*! Called for each read alignment.*/
void alignment_output::Alignment_Output_Pileup::fetch_callback ( const alignment& a )
{
  // we only keep track of unique alignments 
  if ( a.is_redundant() ) return;
  
  total_reads++;

  if ( total_reads > maximum_to_align ) return;

  // create a new empty structure and fill it
  aligned_reads[a.read_name()] = Aligned_Read();
  Aligned_Read& aligned_read = aligned_reads[a.read_name()];
  
  aligned_read.seq_id = a.read_name();
  aligned_read.length = a.read_length();
  aligned_read.read_sequence = a.read_char_sequence();
  aligned_read.qual_sequence = string ( ( char* ) a.read_base_quality_sequence() );

  aligned_reads[aligned_read.seq_id] = aligned_read;

  if ( ( unique_start == 0 ) || ( unique_start > a.reference_start_1() ) )
  {
      unique_start = a.reference_start_1();
  }
  if ( ( unique_end == 0 ) || ( unique_end < a.reference_end_1() ) )
  {
      unique_end = a.reference_end_1();
  }
}

// 'quality_score_cutoff' = below this value you get a special color -- for bad Illumina bases
void alignment_output::set_quality_range(const uint32_t quality_score_cutoff)
{
    // calculate a cumulative distribution of the bases we are showing
    vector<uint32_t> qc(255, 0);
    uint32_t total = 0;

    for ( Aligned_Reads::iterator itr_read = m_aligned_reads.begin(); itr_read != m_aligned_reads.end(); itr_read++ )
    {
      Aligned_Read& aligned_read(itr_read->second );
      for ( uint32_t index = 0; index < aligned_read.qual_sequence.length(); index ++ ) //TODO check if 1 indexed or not
      {
        uint8_t c = static_cast<uint8_t>(aligned_read.qual_sequence[index]);
        if ( c >= quality_score_cutoff )
        {
          qc[c]++;
          total++;
        }
      }
    }
    
    map<uint,uint8_t> qual_to_color;
    uint32_t cutoff_percentiles[] = {0, 0.03, 0.1, 0.3, 0.9, 1.0};
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
    while ( index < qc.size() )
    {
        cumq += (double)qc[index] / (double)total;

        // this can increment by at most one per quality score
        if ( cumq > cutoff_percentiles[current_cutoff_level] )
        {
          current_cutoff_level++;
        }
        qual_to_color[index] = current_cutoff_level;
        index++;       
    }
    //#last must be set to max
    qual_to_color[index - 1] = num_cutoff_percentiles - 1;
    //#first must be set to min
    qual_to_color[quality_score_cutoff] = 1;

    //#if there are at least as many quality scores in existence as
    //#there are color levels to assign....

    if ( qual_to_color.size() > ( num_cutoff_percentiles - 1 ) )
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
            qual_to_color[index - 1]++;
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

    for ( uint32_t index = quality_score_cutoff;
            index < qual_to_color.size(); index++ )
    {
        if ( qual_to_color[index] > last )
        {
            cutoff_levels.push_back ( index );
            last = qual_to_color[index];
        }
    }

    m_quality_range.qual_to_color_index = qual_to_color;
    m_quality_range.qual_cutoffs = cutoff_levels;
}


string alignment_output::create_header_string()
{
    map<char, vector<string> > base_color_hash;
    base_color_hash['G'].push_back ( "rgb(255,255,0)" );
    base_color_hash['G'].push_back ( "rgb(230,230,230)" );
    base_color_hash['G'].push_back ( "rgb(210,210,210)" );
    base_color_hash['G'].push_back ( "rgb(140,140,140)" );
    base_color_hash['G'].push_back ( "rgb(70,70,70)" );
    base_color_hash['G'].push_back ( "rgb(0,0,0)" );

    base_color_hash['C'].push_back ( "rgb(255,255,0)" );
    base_color_hash['C'].push_back ( "rgb(160,160,255)" );
    base_color_hash['C'].push_back ( "rgb(120,120,255)" );
    base_color_hash['C'].push_back ( "rgb(60,60,255)" );
    base_color_hash['C'].push_back ( "rgb(0,0,255)" );
    base_color_hash['C'].push_back ( "rgb(0,0,150)" );

    base_color_hash['A'].push_back ( "rgb(255,255,0)" );
    base_color_hash['A'].push_back ( "rgb(255,210,210)" );
    base_color_hash['A'].push_back ( "rgb(255,180,180)" );
    base_color_hash['A'].push_back ( "rgb(255,100,100)" );
    base_color_hash['A'].push_back ( "rgb(255,20,20)" );
    base_color_hash['A'].push_back ( "rgb(200,0,0)" );

    base_color_hash['T'].push_back ( "rgb(255,255,0)" );
    base_color_hash['T'].push_back ( "rgb(210,255,210)" );
    base_color_hash['T'].push_back ( "rgb(180,255,180)" );
    base_color_hash['T'].push_back ( "rgb(100,255,100)" );
    base_color_hash['T'].push_back ( "rgb(20,255,20)" );
    base_color_hash['T'].push_back ( "rgb(0,200,0)" );

    base_color_hash['N'].push_back ( "rgb(128,0,128)" );
    base_color_hash['N'].push_back ( "rgb(128,0,128)" );
    base_color_hash['N'].push_back ( "rgb(128,0,128)" );
    base_color_hash['N'].push_back ( "rgb(128,0,128)" );
    base_color_hash['N'].push_back ( "rgb(128,0,128)" );
    base_color_hash['N'].push_back ( "rgb(128,0,128)" );

    string header_style_string;
    header_style_string += ".NC {color: rgb(0,0,0); background-color: rgb(255,255,255)}\n";
    header_style_string += ".UN {color: rgb(120,120,120); background-color: rgb(255,255,255)}\n";

  
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
                header_style_string += "."+to_string ( ( *itr ).first) + to_string(index) + "{color: rgb(120,120,120); background-color:"
                                       + ( *itr ).second[index] + "}\n";
            }
        }
    }
    no_color_index = base_color_hash['G'].size()-1;
    return header_style_string;
}



string alignment_output::html_alignment_line(const alignment_output::Alignment_Base& a, const bool coords, const bool use_quality_range)
{
  string output;
  output += "<CODE>";

  if (!a.aligned_quals.empty())
  {
    if (a.aligned_bases.length() != a.aligned_quals.length())
    {
      cout << a.aligned_bases << endl;;
      cout << a.aligned_quals << endl;
      cerr << "unequal aligned base and aligned quals" <<endl;
    }
  }   
  string last_color = "";

  for(uint32_t index = 0; index < a.aligned_bases.length(); index++)
  {
    uint8_t q = 255;
    if(!a.aligned_quals.empty())
    {
      q = static_cast<uint8_t>(a.aligned_quals[index]); 
    }
    
    string b=to_string(a.aligned_bases[index]);       
    string color = "";
  
    // no base quality provided -- assume BEST if not space
    if (q != 254)
    {
      if ( (q == 255) || (!use_quality_range) )
      {
        if (b == " ")
          color = "NC";   
        else
          color = to_upper(to_string(b)) + to_string<uint32_t>(no_color_index);
      }
      else if ( (use_quality_range) &&
         ( (b.find(".") == string::npos) || (b.find("-") == string::npos) ) )
      {
        uint8_t color_num = m_quality_range.qual_to_color_index[q];
        color = to_upper(b) + to_string<uint32_t>(color_num);
      }
    }     
  
    if (b == " ") b = "&nbsp";
    if (color.empty())
    {
      color = "UN";
    } 

    if (color != last_color)
    {
      if (!last_color.empty())
        output += "</font>"; 
      output += "<font class=\"" + color + "\">"; 
      last_color = color;
    }
    output += b;
  }
  if (!last_color.empty()) output += "</font>";

  if (a.show_strand) output += "&nbsp;&nbsp;" + html_alignment_strand(a.strand); 

  // write the seq_id and coords in non-breaking html
  if (!a.seq_id.empty())
  {
    string seq_id = substitute(a.seq_id, "-", "&#8209;");
        
    if(coords)
    {
      seq_id += "/" + to_string<uint32_t>(a.start);
      seq_id += "&#8209;";
      seq_id += to_string<uint32_t>(a.end);   
    }
    output += "&nbsp;&nbsp;" + seq_id;
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

} // namespace breseq




