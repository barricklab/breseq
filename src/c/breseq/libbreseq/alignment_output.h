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
 *  Copyright (c) 2011-2022 The University of Texas at Austin
 * 
 *  breseq is free software; you can redistribute it and/or modify it under the
 *  terms the GNU General Public License as published by the Free Software
 *  Foundation; either version 1, or (at your option) any later version.
 * 
 *****************************************************************************/

#ifndef _BRESEQ_ALIGNMENT_OUTPUT_H_
#define _BRESEQ_ALIGNMENT_OUTPUT_H_

#include "common.h"
#include "pileup.h"
#include "candidate_junctions.h"


namespace breseq
{
  
  // Pre-declaration
  namespace output {
    class cOutputEvidenceItem;  
  }
  using namespace output;
  
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
        , strand(0)
        , start_is_clipped(false)
        , end_is_clipped(false)
        , aligned_bases("")
        , aligned_quals("")
        , is_redundant(false)
        , show_strand(true)
        , mapping_quality(-1)
        , read_name_style("NC")
        {}
      
      string seq_id;
      uint32_t start;
      uint32_t end;
      int8_t strand;
      bool start_is_clipped;
      bool end_is_clipped;
      string aligned_bases;
      string aligned_quals;
      bool is_redundant;
      bool show_strand;
      int32_t mapping_quality;
      string read_name_style;
      
      static const string gap_string;
      
      void truncate(size_t start_pos, size_t end_pos) {
        size_t len = end_pos - start_pos + 1;
        
        string clipped_start = aligned_bases.substr(0, start_pos);
        string clipped_end = aligned_bases.substr(end_pos+1, aligned_bases.length()-end_pos);
        
        // Debug code
        /*
        // Count up how many bases and whether any bases were removed
        size_t actual_length = 0;
        for (size_t i=0; i<aligned_bases.length(); i++) {
          
          if (gap_string.find(aligned_bases[i]) == std::string::npos) {
            actual_length++;
          }
        }
        
        string unclipped = aligned_bases.substr(start_pos, len);
        
        cout << actual_length << endl;
        cout << start << endl;
        cout << end << endl;
        cout << static_cast<int32_t>(strand) << endl;
        cout << aligned_bases.length();
        cout << clipped_start << endl;
        cout << endl;
        cout << clipped_end << endl;
        cout << endl;
        cout << unclipped << endl;
        cout << endl;
        */
        // Count up how many bases and whether any bases were removed
        for (size_t i=0; i<clipped_start.length(); i++) {
          
          if (gap_string.find(clipped_start[i]) == std::string::npos) {
            start_is_clipped = true;
            start += strand;
          }
        }
        
        for (size_t i=0; i<clipped_end.length(); i++) {
          if (gap_string.find(clipped_end[i]) == std::string::npos) {
            end_is_clipped = true;
            end -= strand;
          }
        }
        
        aligned_bases = aligned_bases.substr(start_pos, len);

        // May be empty
        if (aligned_quals.length()) {
          
          aligned_quals = aligned_quals.substr(start_pos, len);
        }
        
      }
      
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
        , overlap_assigned(0)
        , overlap_not_assigned(0)
         {}
         
      uint32_t reference_length;
      string reference_name;
      uint32_t truncate_start, truncate_end;
      uint32_t target_id;
      int8_t ghost_strand;
      uint32_t ghost_start, ghost_end;
      string ghost_seq_id;
      uint32_t overlap_assigned, overlap_not_assigned;
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
                               const bool show_ambiguously_mapped,
                               const uint32_t minimum_mapping_quality = 0
                               );
      //! Destructor.
      virtual ~Alignment_Output_Pileup();
      //! Called for each genome position.
      virtual void pileup_callback ( const pileup& aligned_reference );
      //! Called for each aligned read.
      virtual void fetch_callback ( const alignment_wrapper& a );
      
      bool _show_ambiguously_mapped;
      uint32_t _minimum_mapping_quality;
      
      Aligned_Reads aligned_reads;
      Aligned_References aligned_references;
      Aligned_Annotation aligned_annotation;
      uint32_t unique_start; //used in create alignment and passed to fetch
      uint32_t unique_end;   //used in create alignment and passed to fetch
      
      uint32_t max_indel;
      char base;
    };

    // Internal data and classes
    Alignment_Output_Pileup m_alignment_output_pileup;
    Aligned_Reads m_aligned_reads;
    Aligned_References m_aligned_references;
    Aligned_Annotation m_aligned_annotation;
    
    // Display options
    Quality_Range m_quality_range;
    uint32_t m_quality_score_cutoff;
    string m_error_message;
    uint32_t m_maximum_flanking_columns;
    uint32_t m_maximum_to_align;
    int32_t m_junction_minimum_size_match;
    bool m_mask_ref_matches;    // Show matches to reference as '.' rather than base.
    bool m_show_ambiguously_mapped;
    uint32_t m_minimum_mapping_quality; // Show reads as redundantly matched if their mapping quality is below this
    bool m_is_junction, m_is_junction_junction; // Former for JC, latter for JC+junction part of it
    
    // Characters substituted for gaps and reference matches
    static char s_internal_gap_character;             // Default = "-"
    static char s_end_gap_character;                  // Default = "-"
    static char s_reference_match_character;          // Default = "."
    static char s_reference_match_masked_character;   // Default = ","
    
    static uint8_t k_reserved_quality_start;
    static uint8_t k_reserved_quality_junction_overlap;
    static uint8_t k_reserved_quality_dont_highlight;
    static uint8_t k_reserved_quality_max;
    
  public:
    //! Constructor.
    alignment_output ( 
                      string bam, 
                      string fasta,
                      uint32_t in_maximum_flanking_columns = 0,
                      uint32_t in_maximum_to_align = 0,
                      const uint32_t quality_score_cutoff = 0,
                      const int32_t junction_minimum_size_match = 1,
                      const bool mask_ref_matches = false,
                      const bool show_ambiguously_mapped = false,
                      const uint32_t minimum_mapping_quality = 0
                      );
    //! Output an HTML alignment.
    void create_alignment ( const string& region, cOutputEvidenceItem * output_evidence_item = NULL );
    string html_alignment ( const string& region, cOutputEvidenceItem * output_evidence_item = NULL );
    string text_alignment ( const string& region, cOutputEvidenceItem * output_evidence_item = NULL );
    string json_alignment ( const string& region, cOutputEvidenceItem * output_evidence_item = NULL );
    void set_quality_range(const uint32_t quality_score_cutoff = 0);
  private:
    uint32_t no_color_index;
    string html_header_string();
    string html_alignment_line(const Alignment_Base& a, Aligned_References* r, const bool coords, const bool use_quality_range, const bool show_coord_clipping);
    string html_alignment_strand(const int8_t &strand);
    string html_legend();
    string text_alignment_line(const Alignment_Base& a, const bool coors, const bool show_coord_clipping);
    string text_alignment_strand(const int8_t &strand);
    json json_alignment_line(const Alignment_Base& a);

    
    static bool sort_by_aligned_bases_length ( const Sorted_Key& a, const Sorted_Key& b )
    {
      return ( a.aligned_bases.compare(b.aligned_bases) > 0 );
    }
  };
  
}//end namespace breseq
#endif
