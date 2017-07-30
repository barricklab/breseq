/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2017 The University of Texas at Austin

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
class alignment_wrapper;
int first_level_pileup_callback(uint32_t tid, uint32_t pos, int32_t n, const bam_pileup1_t *pile, void *data);
int first_level_fetch_callback(bam1_t *b, void *data);

//! Helper struct to manage a single reference sequence.
struct reference_sequence {
  public:
    reference_sequence(faidx_t* m_ref, const string& fasta_filename, const string& target);
    ~reference_sequence();

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
      ASSERT(target < num_targets(), "Requested target_id [" + to_string(target) + "] was not found in FASTA file [" + m_fasta_file_name + "]");
      return m_bam->header->target_name[target];
    }

    uint32_t num_targets() const {
      return m_bam->header->n_targets;
    }
  
    //! Retrieve the length of the given target.
    uint32_t target_length(uint32_t target) const {
      return m_refs[target]->m_len;
    }

    char reference_base_char_1(uint32_t tid, uint32_t pos1) const  {
      ASSERT(pos1 <= target_length(tid), "1-indexed position out of bounds! targed_id: " + to_string(tid) 
             + " position: " + to_string(pos1) + " [1," + to_string(target_length(tid)) + "].");
      ASSERT(pos1 >= 1, "1-indexed position out of bounds! targed_id: " + to_string(tid) 
             + " position: " + to_string(1) + " [1," + to_string(target_length(tid)) + "].");
      return get_refseq(tid)[pos1-1];
    } ;
  
    char reference_base_char_0(uint32_t tid, uint32_t pos0) const  {
      ASSERT(pos0 < target_length(tid), "0-indexed position out of bounds! targed_id: " + to_string(tid) 
             + " position: " + to_string(pos0) + " [0," + to_string(target_length(tid)-1) + "].");
      //ASSERT(pos0 >= 0, "0-indexed position out of bounds! targed_id: " + to_string(tid) 
      //       + " position: " + to_string(0) + " [0," + to_string(target_length(tid)-1) + "].");
      return get_refseq(tid)[pos0];
    } ;

    // handle this reference sequence position during pileup?
    bool handle_position(uint32_t pos_1);

    //! Do the pileup;  (Callback for each position, with information about read alignments there.)
    //  Note that we call for all positions that have been skipped (with zero alignments)
    //  which is unlike the default SAM behaviour.
    void do_pileup();

    //! Do the pileup, but only on specified region.
    void do_pileup(const string& region, bool clip = false, uint32_t downsample = 0);

    //! Do the pileup, but only on specified seq_ids.
    void do_pileup(const set<string>& seq_ids);
  
    //! Do the fetch, (Callback for each read alignment to region.)
    void do_fetch(const string& region);

    //! Pileup callback.
    virtual void pileup_callback(const pileup& p) {
      (void)p;
      ASSERT(false, "pileup_callback not defined for class");
    };


    //! Fetch callback.
    virtual void fetch_callback(const alignment_wrapper& a) {
      (void)a;
      ASSERT(false, "fetch_callback not defined for class");
    };
  
    //! Called before pileup starts a target.
    virtual void at_target_start(const uint32_t tid) { (void)tid; }

    virtual void at_target_start_first_level_callback(const uint32_t tid) { 
      if (m_print_progress) {
        cerr << "  REFERENCE: " << m_bam->header->target_name[tid] << endl;
        cerr << "  LENGTH: " << m_bam->header->target_len[tid] << endl;
      }
      at_target_start(tid);
    }
  
    //! Called after the pileup completed a target.
    virtual void at_target_end(const uint32_t tid) { (void)tid; }
  
    virtual void at_target_end_first_level_callback(const uint32_t tid) { 
      at_target_end(tid);
    }
  
    vector<string> valid_seq_ids() const
    {
      vector<string> ret_val;
      for(uint32_t i=0; i< num_targets(); i++) {
        ret_val.push_back(target_name(i));
      }
      return ret_val;
    }
  
    int32_t seq_id_to_target_id(const string& seq_id) const
    {
      for(uint32_t i=0; i< num_targets(); i++) {
        if (target_name(i) == seq_id) {
          return i;
        }
      }
      return -1;
    }
  
  
    //! Special parsing for seq_id:start.insert_start-end.insert_end form.
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

      // Save target id
      int32_t temp_target_id = seq_id_to_target_id(target_name);
      // Target was not found.
      ASSERT(temp_target_id != -1, "Target seq id was not found for region [" + target_name + "] using FASTA file [" + m_fasta_file_name + "].\n" + "Valid seq ids: " + join(this->valid_seq_ids(), ", ") );
      target_id = static_cast<uint32_t>(temp_target_id);
      
      //
      start_pos_1 = from_string<uint32_t>(start_pos_1_string);
      end_pos_1 = from_string<uint32_t>(end_pos_1_string);
      
      insert_start = from_string<uint32_t>(insert_start_string);
      insert_end = from_string<uint32_t>(insert_end_string);
    }
  
    void parse_region(const string& region, uint32_t& target_id, uint32_t& start_pos_1, uint32_t& end_pos_1)
    {
      uint32_t temp_insert_start, temp_insert_end;
      parse_region(region, target_id, start_pos_1, end_pos_1, temp_insert_start, temp_insert_end);
    }
  
  inline void set_print_progress(bool print_progress) { m_print_progress = print_progress; }

  protected:
    friend int first_level_pileup_callback(uint32_t tid, uint32_t pos, int32_t n, const bam_pileup1_t *pile, void *data);
    friend int first_level_fetch_callback(const bam1_t *b, void *data);

    samfile_t* m_bam; //!< BAM file handle.
    bam_header_t* m_bam_header;
    bam_index_t* m_bam_index;
    bamFile m_bam_file;
    faidx_t* m_faidx;
    string m_bam_file_name;
    string m_fasta_file_name;


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
