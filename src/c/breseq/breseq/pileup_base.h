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

#ifndef _BRESEQ_PILEUP_BASE_H_
#define _BRESEQ_PILEUP_BASE_H_

#include "common.h"

using namespace std;
namespace breseq {

// pre-decs
class pileup;
class alignment;
int first_level_pileup_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data);
int first_level_fetch_callback(bam1_t *b, void *data);

//! Helper struct to manage a single reference sequence.
struct reference_sequence {
  public:
    reference_sequence(const string& fasta_filename, const string& target);
    ~reference_sequence();

    faidx_t* m_ref; //!< FAI file handle.
    char* m_seq; //!< Reference sequence (ascii).
    int m_len; //<! Length of reference sequence.
// 
  private:
    // not allowed:
    reference_sequence(const reference_sequence& that);
    reference_sequence& operator=(const reference_sequence& that);
};


/*! Class to assist in developing pileup-related functionality.
 */
class pileup_base {
  public:
    //! Type for a list of reference sequences.
    typedef vector<reference_sequence*> refseq_list_t;

    //! Constructor.
    pileup_base(const string& bam, const string& fasta);

    //! Destructor.
    virtual ~pileup_base();

    //! Retrieve the reference sequence for the given target and fai index.
    char* get_refseq(uint32_t target) const;

    //! Retrieve the name of the given target.
    const char* target_name(uint32_t target) const {
        return m_bam->header->target_name[target];
    }

    const uint32_t num_targets() const {
      return m_bam->header->n_targets;
    }
  
    //! Retrieve the length of the given target.
    const uint32_t target_length(uint32_t target) const {
      return m_refs[target]->m_len;
    }

    char reference_base_char_1(uint32_t target, uint32_t pos1) const  {
      return get_refseq(target)[pos1-1];
    } ;
  
    char reference_base_char_0(uint32_t target, uint32_t pos0) const  {
      return get_refseq(target)[pos0];
    } ;

    // handle this reference sequence position during pileup?
    bool handle_position(uint32_t pos_1);

    //! Do the pileup;  (Callback for each position, with information about read alignments there.)
    //  Note that we call for all positions that have been skipped (with zero alignments)
    //  which is unlike the default SAM behaviour.
    void do_pileup();

    //! Do the pileup, but only on specified region.
    void do_pileup(const string& region, bool clip = false, uint32_t downsample = 0);

    //! Do the fetch, (Callback for each read alignment to region.)
    void do_fetch(const string& region);

    //! Pileup callback.
    virtual void pileup_callback(const pileup& p) {
        assert(false);
    };

    //! Fetch callback.
    virtual void fetch_callback(const alignment& a) {
        assert(false);
    };
  
    //! Called before pileup starts a target.
    virtual void at_target_start(uint32_t tid) { }
  
    //! Called after the pileup completed a target.
    virtual void at_target_end(const uint32_t tid) { }
  
    //! Pass through to BAM.
    //! @JEB make this transparently handle insert_start and insert_end REL606:13.1-16.0 using split_on_any
    void parse_region(const string& region, uint32_t& target_id, uint32_t& start_pos_1, uint32_t& end_pos_1)
    {
      int temp_target_id, temp_start_pos, temp_end_pos;
      bam_parse_region(m_bam_header, region.c_str(), &temp_target_id, &temp_start_pos, &temp_end_pos); 
      target_id = static_cast<uint32_t>(temp_target_id);
      start_pos_1 = static_cast<uint32_t>(temp_start_pos)+1; // bam_parse_region returns zero indexed start
      end_pos_1 = static_cast<uint32_t>(temp_end_pos);       // bam_parse_region returns one indexed end
    }
  
    void parse_region(const string& region, uint32_t& target_id, uint32_t& start_pos_1, uint32_t& end_pos_1, uint32_t& insert_start, uint32_t& insert_end)
    {
      insert_start = 0;
      insert_end = 0;
            
      size_t start = 0;
      size_t end = 0;
      
      end = region.find_first_of(':', start);
      assert(end != string::npos);
      string target_name = region.substr(start, end - start);
      start = end+1;
      
      string start_pos_1_string;
      string insert_start_string("0");
      end = region.find_first_of('.', start);
      if (end == string::npos) {
        end = region.find_first_of('-', start);
        assert(end != string::npos);
        start_pos_1_string = region.substr(start, end - start);
        start = end+1;
      }
      else
      {
        start_pos_1_string = region.substr(start, end - start);
        start = end+1;
        end = region.find_first_of('-', start);
        assert(end != string::npos);
        insert_start_string = region.substr(start, end - start);
        start = end+1;
      }

      string end_pos_1_string;
      string insert_end_string("0");
      end = region.find_first_of('.', start);
      if (end == string::npos) {
        end = string::npos;
        end_pos_1_string = region.substr(start, end - start);
        start = end+1;
      }
      else
      {
        end_pos_1_string = region.substr(start, end - start);
        start = end+1;
        end = string::npos;
        insert_end_string = region.substr(start, end - start);
        start = end+1;
      }

      string new_region = target_name + ":" + start_pos_1_string + "-" + end_pos_1_string;
      parse_region(new_region, target_id, start_pos_1, end_pos_1);
      
      insert_start = from_string<uint32_t>(insert_start_string);
      insert_end = from_string<uint32_t>(insert_end_string);
    }
  
  inline void set_print_progress(bool print_progress) { m_print_progress = print_progress; }

  protected:
    friend int first_level_pileup_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data);
    friend int first_level_fetch_callback(const bam1_t *b, void *data);

    samfile_t* m_bam; //!< BAM file handle.
    bam_header_t* m_bam_header;
    bam_index_t* m_bam_index;
    bamFile m_bam_file;

    uint32_t m_last_position_1;        // last position handled by pileup
    uint32_t m_start_position_1;       // requested start, 0 = whole fragment
    uint32_t m_end_position_1;         // requested end,   0 = whole fragment
    uint32_t m_insert_start;           // requested insert start (#bases after reference base inserted in read)
    uint32_t m_insert_end;             // requested insert end (#bases after reference base inserted in read)
    uint32_t m_clip_start_position_1;  // clip columns handled starting here, 0 = off
    uint32_t m_clip_end_position_1;    // clip columns handled ending here,   0 = off
    uint32_t m_downsample;
    
    refseq_list_t m_refs; //!< Reference sequences.
    uint32_t m_last_tid; //!< The "last target" for which the first-level-callback was called. -1 = none
  
    bool m_print_progress;
};

} // breseq

#endif
