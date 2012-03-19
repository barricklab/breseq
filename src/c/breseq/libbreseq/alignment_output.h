/*****************************************************************************
 * 
 * AUTHORS
 * 
 *  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 *  David B. Knoester
 * 
 * LICENSE AND COPYRIGHT
 * 
 *  Copyright (c) 2008-2010 Michigan State University
 *  Copyright (c) 2011 The University of Texas at Austin
 * 
 *  breseq is free software; you can redistribute it and/or modify it under the
 *  terms the GNU General Public License as published by the Free Software
 *  Foundation; either version 1, or (at your option) any later version.
 * 
 *****************************************************************************/

#ifndef _BRESEQ_ALIGNMENT_OUTPUT_H_
#define _BRESEQ_ALIGNMENT_OUTPUT_H_

#include "common.h"
#include "pileup_base.h"
#include "pileup.h"
#include "alignment.h"
#include "candidate_junctions.h"

using namespace std;

namespace breseq
{
  
  /*! This class is a FACTORY for generating HTML alignments
   */
  class alignment_output
  {
  public:
    //! returns more information about aligned reads given a sequence id string.
    struct Alignment_Base
    {
      Alignment_Base() 
        : seq_id("")
        , start(0)
        , end(0)
        , aligned_bases("")
        , aligned_quals("")
        , strand(0)
        , is_redundant(false)
        , show_strand(true)
        {}
      
      string seq_id;
      uint32_t start;
      uint32_t end;
      string aligned_bases;
      string aligned_quals;
      int8_t strand;
      bool is_redundant;
      bool show_strand;
      
    };
    //! returns more information about aligned reads given a sequence id string.
    struct Aligned_Read : Alignment_Base
    {
      Aligned_Read() 
        : length(0)
        , read_sequence("")
        , qual_sequence("")
        , reference_start(0) 
        , reference_end(0) 
        , updated(false)
      { }
            
      uint32_t length;
      string read_sequence;
      string qual_sequence;
      uint32_t reference_start;
      uint32_t reference_end;
      bool updated; //whether the read was updated at this pileup iteration already
    };
    //! returns more information about an aligned references
    struct Aligned_Reference : Alignment_Base
    {
      Aligned_Reference() 
        : reference_length(0)
        , reference_name("")  
        , truncate_start(0)
        , truncate_end(0)
        , target_id(-1)
        , ghost_strand(0)
        , ghost_start(0)
        , ghost_end(0)
        , ghost_seq_id("")
        , overlap(0)
         {}
         
      uint32_t reference_length;
      string reference_name;
      uint32_t truncate_start, truncate_end;
      uint32_t target_id;
      int8_t ghost_strand;
      uint32_t ghost_start, ghost_end;
      string ghost_seq_id;
      uint32_t overlap;
    };
    
    struct Aligned_Annotation: Alignment_Base
    { 
      Aligned_Annotation() 
      { show_strand = false; }
    };
    
    struct Sorted_Key
    {
      Sorted_Key() 
        : seq_id("")
       , aligned_bases ("")
      { }
      
      Sorted_Key(const Sorted_Key& copy)
      {
        seq_id = copy.seq_id;
        aligned_bases = copy.aligned_bases;
      }


      string seq_id;
      string aligned_bases;
    };
    
    //!Helper struct for set_quality_range
    typedef struct
    {
      vector<uint8_t> qual_to_color_index;
      vector<uint8_t> qual_cutoffs;
    }Quality_Range;
    
    typedef map<string, Aligned_Read> Aligned_Reads;
    typedef vector<Aligned_Reference> Aligned_References;
    typedef vector<Sorted_Key> Sorted_Keys;
    
  private:
    //! Builds Aligned_Reads, Aligned_References and Aligned_Annotation
    class Alignment_Output_Pileup : public pileup_base
    {
    public:
      //! Constructor.
      Alignment_Output_Pileup ( 
                               const string& bam, 
                               const string& fasta, 
                               const bool show_ambiguously_mapped 
                               );
      //! Destructor.
      virtual ~Alignment_Output_Pileup();
      //! Called for each genome position.
      virtual void pileup_callback ( const pileup& aligned_reference );
      //! Called for each aligned read.
      virtual void fetch_callback ( const alignment_wrapper& a );
      
      bool _show_ambiguously_mapped;
      Aligned_Reads aligned_reads;
      Aligned_References aligned_references;
      Aligned_Annotation aligned_annotation;
      uint32_t unique_start; //used in create alignment and passed to fetch
      uint32_t unique_end; //used in create alignment and passed to fetch
      uint32_t total_reads;
      uint32_t processed_reads;
      
      uint32_t max_indel;
      char base;
    };
    
    Alignment_Output_Pileup m_alignment_output_pileup;
    Aligned_Reads m_aligned_reads;
    Aligned_References m_aligned_references;
    Aligned_Annotation m_aligned_annotation;
    Quality_Range m_quality_range;
    uint32_t m_quality_score_cutoff;
    string m_error_message;
    uint32_t m_maximum_to_align;
    bool m_show_ambiguously_mapped;

    
  public:
    //! Constructor.
    alignment_output ( 
                      string bam, 
                      string fasta, 
                      uint32_t in_maximum_to_align = 0, 
                      const uint32_t quality_score_cutoff = 0,
                      const bool show_ambiguously_mapped = false
                      );
    //! Output an HTML alignment.
    string html_alignment ( const string& region, const string& corrected="" );
    void create_alignment ( const string& region, const string& corrected="" );
    void set_quality_range(const uint32_t quality_score_cutoff = 0);
  private:
    uint32_t no_color_index;
    string create_header_string();
    string html_alignment_line(const Alignment_Base& a, const bool coords, const bool use_quality_range);
    string html_alignment_strand(const int8_t &strand);
    
    static bool sort_by_aligned_bases_length ( const Sorted_Key& a, const Sorted_Key& b )
    {
      return ( a.aligned_bases.compare(b.aligned_bases) > 0 );
    }
  };
  
  
}//end namespace breseq
#endif
