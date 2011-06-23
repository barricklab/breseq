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
#include "breseq/common.h"
#include "breseq/pileup.h"


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

alignment_output::alignment_output ( string bam, string fasta, uint32_t maximum_to_align )
        : m_alignment_output_pileup ( bam, fasta, maximum_to_align )
        ,m_aligned_reads ( m_alignment_output_pileup.aligned_reads )
{
}

void alignment_output::create_alignment ( const string bam, const string fasta, const string region )
{
    //sub create_alignment
    //{
    //  my ($self, $bam_path, $fasta_path, $region, $options) = @_;
    //  my $verbose = $options->{'verbose'};
    //
    //  my ($seq_id, $start, $end, $insert_start, $insert_end) = Breseq::Shared::region_to_coords($region, $bam);
    //  $region = "$seq_id:$start-$end";
    //  print "$bam_path  $fasta_path  $region\n" if ($verbose);
    //
    //  my $reference_length = $bam->length($seq_id);
    //  my $aligned_reads;
    //  my $aligned_annotation;
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
    //
    //  ## Need to ignore positions on the left and right that are only overlapped by
    //  ## Redundant reads. Currently Pileup will start and end on those coords.
    //  ## Since it doesn't know about redundancy marking.
    //
    //  my $unique_start;
    //  my $unique_end;
    //
    //  my $total_reads = 0;
    //  my $processed_reads = 0;
    //
    //*Call
    //* do fetch_callback
    //  $bam->fetch($region, $fetch_function);

    //
    //  ## There are not uniquely aligned reads...
    //  if (!defined $unique_start || !defined $unique_end)
    //  {
    //      return;
    //  }

    //
    //  ## If there are WAY too many reads, such that a pileup might take forever, bail...
    //  if ($self->{maximum_to_make_alignment} && ($total_reads > $self->{maximum_to_make_alignment}))
    //  {
    //      return {
    //          message => "Reads exceeded maximum to display alignment. $total_reads reads. (Limit = $self->{maximum_to_make_alignment})",
    //      };
    //  }
    //

    m_alignment_output_pileup.do_fetch ( region );
    //TEST that unique_end, unique_start and aligned_reads map compare to perl bam2aln output


    if ( ( m_alignment_output_pileup.unique_start == 0 ) || ( m_alignment_output_pileup.unique_end == 0 ) )
    {
        cout << "No unique start or end initialized" << endl;
        return;
    }

    if ( m_alignment_output_pileup.total_reads > m_alignment_output_pileup.maximum_to_align )
    {
        cout << "Reads exceeded maximum to display alignment. ";
        cout << m_alignment_output_pileup.total_reads << " reads. ";
        cout << "(Limit = " << m_alignment_output_pileup.maximum_to_align;
        cout << ")" << endl;
        return;
    }
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
    // *
    //
    //  $bam->pileup($region, $pileup_function);

    m_alignment_output_pileup.last_pos = 0;
    m_alignment_output_pileup.do_pileup ( region );
    //do_pileup is needed to create the following

    m_aligned_reads = m_alignment_output_pileup.aligned_reads;
    m_aligned_references = m_alignment_output_pileup.aligned_references;
    m_aligned_annotation = m_alignment_output_pileup.aligned_annotation;
    //NOTE "READ-4421-M1" at [74]
    //
    //
    //    #### IF NOTHING ALIGNED, RETURN undef
    //  return undef if (scalar keys %$aligned_reads == 0);
    //


    if ( m_aligned_reads.empty() )
    {
        cerr << "No aligned reads created in pileup function";
        return;
    }
    //  ##now add the unaligned portions of each
    //  my $max_extend_left = 0;
    //  my $max_extend_right = 0;


    int32_t max_extend_left = 0;
    int32_t max_extend_right = 0;

    //  foreach my $key (keys %$aligned_reads)
    //  {
    //      my $aligned_read = $aligned_reads->{$key};
    //
    //      $aligned_read->{aligned_bases} =~ m/^(\s*)\S+(\s*)$/;
    //
    //#     print "\"$aligned_read->{aligned_bases}\"\n";
    //
    //      print "$aligned_read->{start} $aligned_read->{end}\n" if ($verbose);
    //

    for ( Aligned_Reads::iterator itr_reads = m_aligned_reads.begin(); //TODO typedef map
            itr_reads != m_aligned_reads.end(); itr_reads++ )
    {
        Aligned_Read& aligned_read = ( m_aligned_reads[ ( *itr_reads ).first] );
        //$aligned_read->{aligned_bases} =~ m/^(\s*)\S+(\s*)$/;
        size_t left_padding_length = aligned_read.aligned_bases.find_first_not_of ( ' ' );
        size_t right_padding_length = ( aligned_read.aligned_bases.length() - 1 )
                                        - aligned_read.aligned_bases.find_last_not_of ( ' ' );

      
        if ( verbose )
        {
            cout  << aligned_read.start << " " << aligned_read.end << " " << aligned_read.seq_id << endl;
        }
        //      my $extend_left =  ($aligned_read->{start}-1) - length($1);
        //      my $extend_right = ($aligned_read->{length}-$aligned_read->{end}) - length($2);
        //
        //      $max_extend_left = $extend_left if ($extend_left > $max_extend_left);
        //      $max_extend_right = $extend_right if ($extend_right > $max_extend_right);
        //  }
        //

        int32_t extend_left = ( aligned_read.start - 1 ) - left_padding_length;
        int32_t extend_right = ( aligned_read.length - aligned_read.end ) - right_padding_length;

        if ( extend_left > max_extend_left )
            max_extend_left = extend_left;
        if ( extend_right > max_extend_right )
            max_extend_right = extend_right;
    }
    //
    //  #now add this much to every one
    //  foreach my $key (keys %$aligned_reads)
    //  {
    //      my $aligned_read = $aligned_reads->{$key};
    //      $aligned_read->{aligned_bases} = (' ' x $max_extend_left) . $aligned_read->{aligned_bases} . (' ' x $max_extend_right);
    //      $aligned_read->{aligned_quals} = (chr(255) x $max_extend_left) . $aligned_read->{aligned_quals} . (chr(255) x $max_extend_right);
    //  }
    //
    //  print "Extend: $max_extend_left $max_extend_right\n" if ($verbose);
    for ( map<string, Aligned_Read>::iterator itr_reads = m_aligned_reads.begin(); //TODO typedef map
            itr_reads != m_aligned_reads.end(); itr_reads++ )
    {
        Aligned_Read& aligned_read ( ( *itr_reads ).second );

        aligned_read.aligned_bases.insert ( 0, max_extend_left, ' ' );
        aligned_read.aligned_bases.append ( ' ', max_extend_right );

        aligned_read.aligned_quals.insert ( 0, max_extend_left, char ( 255 ) );
        aligned_read.aligned_quals.append ( char ( 255 ), max_extend_right );
    }

    if ( verbose )
    {
        cout << "Extend: " << max_extend_left;
        cout << " " << max_extend_right << endl;
    }


    //  ###extend reference sequence as requested, be aware of ends of sequence
    //  ## handle left side extending past end
    //  foreach my $aligned_reference (@aligned_references)
    //  {
    for ( uint32_t itr_aligned_reference = 0; itr_aligned_reference <
            m_aligned_references.size(); itr_aligned_reference++ )
    {
        Aligned_Reference& aligned_reference
        = ( m_aligned_references[itr_aligned_reference] );
        //      my $ref_extend_left = $max_extend_left;
        //      my $ref_add_left = '';
        uint32_t ref_extend_left = max_extend_left;
        string ref_add_left;
        ref_add_left.clear();
        //         if ($aligned_reference-> {start}-$max_extend_left < 1)
        //         {
        //             $ref_extend_left = $aligned_reference-> {start} - 1;
        //             $ref_add_left .= '-' x (1 - ($aligned_reference-> {start}-$max_extend_left));
        //         }

        if ( ( aligned_reference.start - max_extend_left ) < 1 )
        {
            ref_extend_left = aligned_reference.start -1;
            ref_add_left.append ( '-', 1 - ( aligned_reference.start
                                             - max_extend_left ) );
        }
        //       ## handle right side extending past end
        //       my $ref_extend_right = $max_extend_right;
        //       my $ref_add_right = '';
        uint32_t ref_extend_right = max_extend_right;
        string ref_add_right;
        ref_add_right.clear();

        //       if ($aligned_reference-> {end}+$max_extend_right > $reference_length)
        //       {
        //           $ref_extend_right = $reference_length - $aligned_reference-> {end};
        //           $ref_add_right .= '-' x (($aligned_reference-> {end}+$max_extend_right) - $reference_length);
        //       }
        //
        if ( ( aligned_reference.end +max_extend_right ) > aligned_reference.reference_length )
        {
            ref_extend_right = aligned_reference.reference_length - aligned_reference.end;
            ref_add_right.append ( '-', ( aligned_reference.end +max_extend_right )
                                   - aligned_reference.reference_length );
        }
        //       my $pos;
        //       $pos = $aligned_reference-> {start}-$ref_extend_left;
        uint32_t pos = aligned_reference.start - ref_extend_left;
        //       while ($pos < $aligned_reference-> {start})
        //       {
        //           my $base = $bam->segment($aligned_reference-> {seq_id},$pos,$pos)->dna;
        //           $base = '.' if ($aligned_reference-> {truncate_start} && ($pos < $aligned_reference-> {truncate_start}))
        //                       || ($aligned_reference-> {truncate_end} && ($pos > $aligned_reference-> {truncate_end}));
        //           $ref_add_left .= $base;
        //           $pos++;
        //       }
        while ( pos < aligned_reference.start )
        {
            char base = aligned_reference.base;
            //TODO implement trucate start and end
            ref_add_left += base;
            pos++;
        }
        //
        //       $pos = $aligned_reference-> {end}+1;
        //       while ($pos <= $aligned_reference-> {end}+$ref_extend_right)
        //       {
        //           my $base = $bam->segment($aligned_reference-> {seq_id},$pos,$pos)->dna;
        //           $base = '.' if ($aligned_reference-> {truncate_start} && ($pos < $aligned_reference-> {truncate_start}))
        //                       || ($aligned_reference-> {truncate_end} && ($pos > $aligned_reference-> {truncate_end}));
        //           $ref_add_right .= $base;
        //           $pos++;
        //       }
        pos = aligned_reference.end + 1;
        while ( pos <= aligned_reference.end + ref_extend_right )
        {
            char base = aligned_reference.base;
            ref_add_right += base;
            pos++;
        }
        //
        //       #   $ref_add_left .= $bam->segment($aligned_reference->{seq_id},$aligned_reference->{start}-$ref_extend_left,$aligned_reference->{start}-1)->dna;
        //       #   $ref_add_right .= $bam->segment($aligned_reference->{seq_id},$aligned_reference->{end}+1,$aligned_reference->{end}+$ref_extend_right)->dna;
        //
        //       $aligned_reference-> {start} -= $ref_extend_left;
        //       $aligned_reference-> {end} += $ref_extend_right;
        //       $aligned_reference-> {aligned_bases} = $ref_add_left . $aligned_reference-> {aligned_bases} . $ref_add_right;
        //       $aligned_reference-> {aligned_quals} = (chr(255) x $max_extend_left) . $aligned_reference-> {aligned_quals} . (chr(255) x $max_extend_right);
        //       }
        aligned_reference.start -= ref_extend_left;
        aligned_reference.end += ref_extend_right;

        aligned_reference.aligned_bases.insert ( 0, ref_add_left );
        aligned_reference.aligned_bases.append ( ref_add_right );

        aligned_reference.aligned_quals.insert ( 0, max_extend_left, char ( 255 ) );
        aligned_reference.aligned_quals.append ( char ( 255 ), max_extend_right );
    }
    //
    //  #extend annotation line
    //  $aligned_annotation->{aligned_bases} = (' ' x $max_extend_left) . $aligned_annotation->{aligned_bases} . (' ' x $max_extend_right);
    //
    m_aligned_annotation.aligned_bases.insert ( 0, max_extend_left, ' ' );
    m_aligned_annotation.aligned_bases.append ( ' ', max_extend_left );
    //  #now go in and replace the empty space adjacent to each with the rest of the read sequence
    //  foreach my $key (keys %$aligned_reads)
    //  {
    for ( map<string, Aligned_Read>::iterator itr_reads = m_aligned_reads.begin(); //TODO typedef map
            itr_reads != m_aligned_reads.end(); itr_reads++ )
    {
        //      my $aligned_read = $aligned_reads->{$key};
        //
        //      $aligned_read->{aligned_bases} =~ m/^(\s*)(\S+)/;
        Aligned_Read& aligned_read ( ( *itr_reads ).second );

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

        //      if ($aligned_read->{start} > 1)
        //      {
        if ( aligned_read.start > 1 )
        {
            //          my $add_seq = substr $aligned_read->{read_sequence}, 0, $aligned_read->{start}-1;
            //          substr($aligned_read->{aligned_bases}, $left_pos-length($add_seq), length($add_seq)) = "\L$add_seq";
            //
            string add_seq = aligned_read.read_sequence.substr ( 0, ( aligned_read.start - 1 ) );
            aligned_read.aligned_bases.replace ( ( left_pos - add_seq.length() ), add_seq.length(), to_lower ( add_seq ) ); //start,size, common.h function

            //#to color according to quals
            //#         $add_seq = substr $aligned_read->{qual_sequence}, 0, $aligned_read->{start}-1;
            //
            //          $add_seq = chr(254) x ($aligned_read->{start}-1);
            //          substr($aligned_read->{aligned_quals}, $left_pos-length($add_seq), length($add_seq)) = $add_seq;
            //      }
            add_seq.clear();
            add_seq.append ( char ( 254 ), aligned_read.start - 1 );
            aligned_read.aligned_quals.replace ( ( left_pos - add_seq.length() ), add_seq.length(), add_seq );
        }
        //      if ($aligned_read->{end} < $aligned_read->{length})
        //      {
        //          my $add_seq = substr $aligned_read->{read_sequence}, $aligned_read->{end}, $aligned_read->{length}-$aligned_read->{end};
        //          substr($aligned_read->{aligned_bases}, $right_pos, length($add_seq)) = "\L$add_seq";
        //
        if ( aligned_read.end < aligned_read.length )
        {
            string add_seq = aligned_read.read_sequence.substr ( aligned_read.end, ( aligned_read.length - aligned_read.end ) );
            aligned_read.aligned_bases.replace ( right_pos, add_seq.length(), to_lower ( add_seq ) );


            //#to color according to quals
            //#         $add_seq = substr $aligned_read->{qual_sequence}, $aligned_read->{end}, $aligned_read->{length}-$aligned_read->{end};
            //
            //          $add_seq = chr(254) x ($aligned_read->{length}-$aligned_read->{end});
            //          substr($aligned_read->{aligned_quals}, $right_pos, length($add_seq)) = $add_seq;
            //      }
            add_seq.clear();
            add_seq.append ( char ( 254 ), ( aligned_read.length - aligned_read.end ) );
            aligned_read.aligned_quals.replace ( right_pos,add_seq.length(), add_seq );
        }
        //  }
        //
    }
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
    //
    //  ## Need to reverse the coords for some
    //  foreach my $key (keys %$aligned_reads)
    //  {
    for ( map<string, Aligned_Read>::iterator itr_reads = m_aligned_reads.begin(); //TODO typedef map
            itr_reads != m_aligned_reads.end(); itr_reads++ )
    {
        //      my $aligned_read = $aligned_reads->{$key};

        Aligned_Read& aligned_read ( ( *itr_reads ).second );

        //      if ($aligned_read->{strand} == -1)
        //      {
        //          ($aligned_read->{start}, $aligned_read->{end}) = ($aligned_read->{length} - $aligned_read->{start} + 1, $aligned_read->{length} - $aligned_read->{end} + 1);
        //      }
        //  }

        if ( aligned_read.strand == -1 )
        {
            aligned_read.start = aligned_read.length - aligned_read.start + 1;
            aligned_read.end = aligned_read.length - aligned_read.end + 1;
        }

    }
    //
    //  return {
    //      aligned_reads => $aligned_reads,
    //      aligned_references => \@aligned_references,
    //      aligned_annotation => $aligned_annotation,
    //      message => $message,
    //  };
    //
} //End create alignment

string alignment_output::html_alignment ( const string& region )
{
    //sub html_alignment
    //{
    //my $verbose = 0;
    //my ($self, $bam_path, $fasta_path, $region, $options) = @_;
    //
    //my $alignment_info = $self->create_alignment($bam_path, $fasta_path, $region, $options);
    //
    //my $output = '';
    //return p . "No reads uniquely align to region." if (!defined $alignment_info);
    //$output .= p . "$alignment_info->{message}" if ($alignment_info->{message});
    //return $output if (!defined $alignment_info->{aligned_reads});
    //
    //my $aligned_reads = $alignment_info->{aligned_reads};
    //my @aligned_references = @{$alignment_info->{aligned_references}};
    //my $aligned_annotation = $alignment_info->{aligned_annotation};
    //my $quality_range = $self->set_quality_range($aligned_reads, $options);
    //
    string output;
    Aligned_Reads& aligned_reads ( m_aligned_reads );
    Aligned_References& aligned_references ( m_aligned_references );
    Aligned_Annotation& aligned_annotation ( m_aligned_annotation );
    Quality_Range& quality_range ( m_quality_range );
    //my @sorted_keys = sort { -($aligned_reads->{$a}->{aligned_bases} cmp $aligned_reads->{$b}->{aligned_bases}) } keys %$aligned_reads;
    Sorted_Keys sorted_keys;
    for (Aligned_Reads::iterator itr_read = aligned_reads.begin();
            itr_read != aligned_reads.end(); itr_read++)
    {
        Sorted_Key sorted_key;
        sorted_key.seq_id = itr_read->second.seq_id;
        sorted_key.aligned_bases_length = itr_read->second.aligned_bases.length();
        sorted_keys.push_back(sorted_key);
    }

    sort(sorted_keys.begin(),sorted_keys.end(),sort_by_aligned_bases_length);


    //$output .= style($self->{header_style_string});
    //$output .= start_table({-style=>"background-color: rgb(255,255,255)"}) . start_Tr() . start_td({-style=>"font-size:10pt"});
    //
    output += "<style>"+create_header_string()+"</style>";
    output += "<table style=""background-color: rgb(255,255,255)"">";
    output += "<tr><td style=""font-size:10pt""><code><font class=""-5"">";

    //foreach my $aligned_reference (@aligned_references)
    //{
    //  $output .= $self->_html_alignment_line($aligned_reference, 1) . br;
    //}
    //$output .= $self->_html_alignment_line($aligned_annotation, 0) . br;
    //
    //     foreach my $key (@sorted_keys)
    //     {
    //         $output .= $self->_html_alignment_line($aligned_reads->{$key}, 0, $quality_range) . br;
    //     }
    //     $output .= $self->_html_alignment_line($aligned_annotation, 0) . br;
    //     foreach my $aligned_reference (@aligned_references)
    //     {
    //         $output .= $self->_html_alignment_line($aligned_reference, 1) . br;
    //     }
    //     $output .= br;
    //
    //     ## create legend information
    //
    //     $output .= start_code . "Base quality scores:&nbsp;" . end_code;
    //     $output .= $self->_html_alignment_line({aligned_bases => 'ATCG', => aligned_quals => pack('CCCC',0,0,0,0)}, 0,  $quality_range);
    //     for (my $i=((defined $options->{quality_score_cutoff}) ? 1 : 2); $i<scalar @{$quality_range->{qual_cutoffs}}; $i++)
    //     {
    //         my $c = $quality_range->{qual_cutoffs}->[$i];
    //         $output .= start_code . "&nbsp;&lt;&nbsp;$c&nbsp;&le;&nbsp;" . end_code;
    //         $output .= $self->_html_alignment_line({aligned_bases => 'ATCG', => aligned_quals => pack('CCCC',$c,$c,$c,$c)}, 0,  $quality_range);
    //     }
    //
    //     $output .= end_table() . end_Tr() . end_td();
    //
    //     return $output;
    // }
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
    for ( uint32_t index_aligned_reference;
            index_aligned_reference < aligned_references.size();
            index_aligned_reference++ )
    {
      // @JEB: 0 counts as undefined since we use 1-indexed coords
      if (aligned_references[index_aligned_reference].start == 0) 
      {
        aligned_references[index_aligned_reference].start = pos_1;  
      }
      aligned_references[index_aligned_reference].end =pos_1;
    }//TEST with multiple reference sequences

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
        //print "MAX INDEL: $max_indel\n" if ($verbose);
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
    //             for ( map<string, Aligned_Read>::const_iterator itr_reads = aligned_reads.begin(); //TODO typedef map
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
    for (Aligned_Reads::iterator itr_reads = aligned_reads.begin(); itr_reads != aligned_reads.end(); itr_reads++ )
    {
        Aligned_Read& aligned_read ( itr_reads->second );
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
				aligned_read.end = itr_pileup->query_position_1(); //BUG correctly initialized here for all but READ 4421-M1

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
							if ( ( trim_left > 0 ) && ( itr_pileup->query_position_1() <= trim_left ) ) //BUG FIX trim_left >0
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
		for ( Aligned_Reads::iterator itr_read = aligned_reads.begin(); //TODO typedef map
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
		for ( uint32_t itr = 0; itr < aligned_references.size(); itr ++ )
    {
      // my $my_ref_base = $ref_base;
    	char my_ref_base = ref_base;
      //// $my_ref_base = '.' if ($aligned_reference-> {truncate_start} && ($last_pos < $aligned_reference-> {truncate_start}))
      //// || ($aligned_reference-> {truncate_end} && ($last_pos > $aligned_reference-> {truncate_end}));
      
			Aligned_Reference& aligned_reference ( aligned_references[itr] );	
      
      // $aligned_reference-> {aligned_bases} .= $my_ref_base;
			aligned_reference.aligned_bases += my_ref_base;
      
      // $aligned_reference-> {aligned_bases} .= '.' x ($max_indel) if ($max_indel > 0);
			if ( max_indel >  0 )
      	aligned_reference.aligned_bases.append ( '.',max_indel );
        // $aligned_reference-> {aligned_quals} .= chr(255) x ($max_indel+1);
				aligned_reference.aligned_quals.append ( char ( 255 ), max_indel +1 );
   
				//TODO Need these for other uses, initialize here?
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


    // ## Retrieve all unique alignments overlapping position with "fetch"
    // ## This lets us know how many slots we need to reserve for alignments.
    // my $fetch_function = sub {
    //     my ($a) = @_;
    // #print $a->display_name,' ',$a->cigar_str,"\n";
    //     my $redundancy = $a->aux_get('X1');
    //
    //     if ((!defined $redundancy) || ($redundancy == 1))
    //     {
    //         $total_reads++;
    //
    //         return undef if ($self-> {maximum_to_make_alignment} && ($total_reads > $self-> {maximum_to_make_alignment}));
    //
    //         my $aligned_read;
    //         $aligned_read-> {seq_id} = $a->display_name;
    //         $aligned_read-> {length} = $a->l_qseq;
    //         $aligned_read-> {read_sequence} = $a->qseq;
    //         $aligned_read-> {qual_sequence} = $a->_qscore;
    //
    // ## save in the hash, creating a spot for each read we will be aligning
    //         $aligned_reads-> {$a->display_name} = $aligned_read;
    //
    // ## keep track of the earliest and latest coords we see in UNIQUE alignments
    //         $unique_start = $a->start if (!defined $unique_start || $unique_start > $a->start);
    //         $unique_end   = $a->end   if (!defined $unique_end   || $unique_end   < $a->end  );
    //     }
    // };
    if ( !a.is_redundant() )
    {
        total_reads++;

        if ( total_reads > maximum_to_align )
            return;

        Aligned_Read aligned_read;
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
}


void alignment_output::set_quality_range()
{
    // sub set_quality_range
    // {
    //     my ($self, $aligned_reads, $options) = @_;

    //
    //     my @qc;
    //     my $total = 0;
    //     my $quality_score_cutoff = 0;
    //     $quality_score_cutoff = $options->{quality_score_cutoff} if (defined $options->{quality_score_cutoff});
    //
    map<uint32_t, uint32_t> qc;
    uint32_t total = 0;
    uint32_t quality_score_cutoff = 0;

    // foreach my $key (keys %$aligned_reads)
    // {
    for ( Aligned_Reads::iterator itr_reads = m_aligned_reads.begin(); //TODO typedef map
            itr_reads != m_aligned_reads.end(); itr_reads++ )
    {
        // my $aligned_read = $aligned_reads->{$key};
        Aligned_Read& aligned_read = ( ( *itr_reads ).second );
        //  foreach my $c (split (//, $aligned_read->{qual_sequence}))
        //  {
        for ( uint32_t index = 0; index < aligned_read.qual_sequence.length(); index ++ ) //TODO check if 1 indexed or not
        {
            uint32_t c = uint32_t ( aligned_read.qual_sequence[index] );
            //    if (ord($c) >= $quality_score_cutoff)
            //    {
            if ( c >= quality_score_cutoff )
            {
                //      $qc[ord($c)]++;
                //      $total++;
                if ( qc.find ( c ) !=qc.end() )
                {
                    qc[c]++;
                    total++;
                }
                else
                {
                    qc[c] = 1;
                }
                //    }
            }
            //  }
        }
        //}
    }
    //
    //my @qual_to_color;
    //my @cutoff_percentiles = (0, 0.03, 0.1, 0.3, 0.9, 1.0);
    //my $current_cutoff_level = 0;
    //
    vector<uint32_t> qual_to_color;
    uint32_t cutoff_values[] = {0, 0.03, 0.1, 0.3, 0.9, 1.0};
    vector<uint32_t> cutoff_percentiles ( cutoff_values, cutoff_values
                                          + sizeof ( cutoff_values ) / sizeof ( uint32_t ) );
    uint32_t current_cutoff_level = 0;


    //##set up to this score to the zero level (which is a completely different color)
    //my $i;
    //for ($i=0; $i<$quality_score_cutoff; $i++)
    //{
    //$qual_to_color[$i] = $current_cutoff_level;
    //}
    //     $current_cutoff_level++;
    uint32_t index;
    for ( index = 0; index <  quality_score_cutoff; index++ )
    {
        qual_to_color[index] = current_cutoff_level;
    }
    current_cutoff_level++;
    //
    //my $cumq = 0;
    //while ($i < scalar @qc)
    //{
    uint32_t cumq = 0;
    while ( qc.find ( index ) != qc.end() )
    {
        //  $cumq += $qc[$i] / $total if (defined $qc[$i]);
        //TODO ASK passes if defined by entering while loop?
        cumq+= ( qc[index] / total );
        //
        // #this can increment by at most one per quality score
        //  if ($cumq > $cutoff_percentiles[$current_cutoff_level])
        //  {
        if ( cumq > cutoff_percentiles[current_cutoff_level] )
        {
            //    $current_cutoff_level++;
            //  }
            current_cutoff_level++;
        }
        // $qual_to_color[$i] = $current_cutoff_level;
        qual_to_color[index] = current_cutoff_level;
        // } continue {
        //         $i++
        //     }
        index++; //TODO ASK
    }
    //#last must be set to max
    //$qual_to_color[$i-1] = scalar(@cutoff_percentiles)-1;
    qual_to_color[index - 1] = cutoff_percentiles.size() - 1;
    //#first must be set to min
    //$qual_to_color[$quality_score_cutoff] = 1;
    qual_to_color[quality_score_cutoff] = 1;
    //
    //#if there are at least as many quality scores in existence as
    //#there are color levels to assign....
    //if ((scalar(@qual_to_color) > scalar(@cutoff_percentiles)-1))
    //{
    if ( qual_to_color.size() > ( cutoff_percentiles.size() - 1 ) )
    {
        //  #...redistribute such that there are no jumps in quality level
        //  my $gap = 1;
        uint32_t gap = 1;
        //  while ($gap)
        //  {
        while ( gap )
        {
            //    $gap = 0;
            gap = 0;
            //    my $last = 0;
            uint32_t last = 0;
            //    for (my $i=0; $i<scalar @qual_to_color; $i++)
            //    {
            for ( uint32_t index = 0; index < qual_to_color.size(); index++ )
            {
                //      if ($qual_to_color[$i] > $last + 1)
                //      {
                if ( qual_to_color[index] > last + 1 )
                {
                    //        $qual_to_color[$i-1]++;
                    qual_to_color[index - 1]++;
                    //        $gap = 1;
                    gap = 1;

                    //      }
                }
                //      $last = $qual_to_color[$i];
                last = qual_to_color[index];
                //    }
            }
            //  }
        }
        //}
    }
    //
    //##finally, this sets the cutoff levels
    //my $last = 0;
    uint32_t last = 0;
    //my @cutoff_levels = ($quality_score_cutoff);
    vector<uint32_t> cutoff_levels;
    cutoff_levels.clear();
    cutoff_levels.push_back ( quality_score_cutoff );
    //for (my $i=$quality_score_cutoff; $i<scalar @qual_to_color; $i++)
    //{
    for ( uint32_t index = quality_score_cutoff;
            index < qual_to_color.size(); index++ )
    {
        //  if ($qual_to_color[$i] > $last)
        //  {
        if ( qual_to_color[index] > last )
        {
            //    push @cutoff_levels, $i;
            cutoff_levels.push_back ( index );
            //    $last = $qual_to_color[$i];
            last = qual_to_color[index];
            //  }
        }
        //}
    }
    //
    //  return {qual_to_color_index => \@qual_to_color, qual_cutoffs => \@cutoff_levels };
    //
    m_quality_range.qual_to_color_index = qual_to_color;
    m_quality_range.qaul_cutoffs = cutoff_levels;
    // return 1
}


string alignment_output::create_header_string()
{
    // our $base_colors_hash = {
    //     'G' => [ "rgb(255,255,0)", "rgb(230,230,230)", "rgb(210,210,210)", "rgb(140,140,140)", "rgb(70,70,70)",  "rgb(0,0,0)"     ],
    //     'C' => [ "rgb(255,255,0)", "rgb(160,160,255)", "rgb(120,120,255)", "rgb(60,60,255)",   "rgb(0,0,255)",   "rgb(0,0,150)"   ],
    //     'A' => [ "rgb(255,255,0)", "rgb(255,210,210)", "rgb(255,180,180)", "rgb(255,100,100)", "rgb(255,20,20)", "rgb(200,0,0)"   ],
    //     'T' => [ "rgb(255,255,0)", "rgb(210,255,210)", "rgb(180,255,180)", "rgb(100,255,100)", "rgb(20,255,20)", "rgb(0,200,0)"   ],
    //     'N' => [ "rgb(128,0,128)", "rgb(128,0,128)",   "rgb(128,0,128)",   "rgb(128,0,128)",   "rgb(128,0,128)", "rgb(128,0,128)" ],
    // };
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

    //$self-> {header_style_string} = '';
    string header_style_string;
    header_style_string.clear();
    //$self-> {header_style_string} .= "\.NC {color: rgb(0,0,0); background-color: rgb(255,255,255)}\n"; #no color
    //$self-> {header_style_string} .= "\.UN {color: rgb(120,120,120); background-color: rgb(255,255,255)}\n"; #unaligned
    header_style_string += "\\.NC {color: rgb(0,0,0); background-color: rgb(255,255,255)}\n";
    header_style_string += "\\.UN {color: rgb(120,120,120); background-color: rgb(255,255,255)}\n";
    //foreach my $key (keys %$base_colors_hash)
    //{
    for ( map<char, vector<string> >::iterator itr = base_color_hash.begin();
            itr != base_color_hash.end(); itr++ )
    {
        //  for (my $i=0; $i<scalar @ {$base_colors_hash->{$key}}; $i++)
        //  {
        for ( uint32_t index = 0; index < base_color_hash[ ( *itr ).first].size(); index ++ )
        {
            if ( index > 0 )
            {
                //    if ($i>0)
                //    {

                //      $self-> {header_style_string} .= "\.$key$i \{color: rgb(255,255,255); background-color: $base_colors_hash->{$key}->[$i]\}\n";
                header_style_string += "\\." + to_string ( ( *itr ).first ) + "\\{color: rgb(255,255,255); background-color:"
                                       + ( *itr ).second[index] + "\\ \n";
                //    }
            }
            //    else
            //    {
            else
            {
                //      $self-> {header_style_string} .= "\.$key$i \{color: rgb(120,120,120); background-color: $base_colors_hash->{$key}->[$i]\}\n";
                header_style_string += "\\."+to_string ( ( *itr ).first ) + "\\{color: rgb(120,120,120); background-color:"
                                       + ( *itr ).second[index] + "\\ \n";
                //    }
            }
            //   }
        }
        //}
    }
    //
    //     $self-> {no_color_index} = scalar(@ {$base_colors_hash->{'G'}}) - 1;
    return header_style_string;
}

} // namespace breseq


