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


#include "pileup_base.h"
#include "pileup.h"

using namespace std;

namespace breseq {
  
/*! Constructor.
 
 Open the FASTA and read the target reference sequence.
 */
  
// @JEB : remove reference_sequence object and 
// only load reference sequences as requested from a single fai index
// stored in the pileup object. This is currently very inefficient.
reference_sequence::reference_sequence(faidx_t* m_ref, const string& fasta_filename, const string& target)
: m_seq(0), m_len(0) {

  (void)fasta_filename;
	m_seq = fai_fetch(m_ref, target.c_str(), &m_len);  
  // need to close the file
  
	assert(m_seq);
	assert(m_len > 0);
}


/*! Destructor.
 */
reference_sequence::~reference_sequence() {
	if(m_seq) { free(m_seq); }
}


/*! Constructor for single-BAM, >=0 FASTA.
 
 \todo Change this to lazily load & cache reference sequences as they are needed,
 instead of loading all sequences at once.
 */
pileup_base::pileup_base(const string& bam, const string& fasta)
: m_bam(0), m_bam_header(0), m_bam_index(0), m_faidx(0), m_bam_file_name(bam), m_fasta_file_name(fasta),
  m_last_position_1(0), m_start_position_1(0), m_end_position_1(0), m_clip_start_position_1(0), m_clip_end_position_1(0),
  m_downsample(0), m_last_tid(static_cast<uint32_t>(-1)), m_print_progress(false)
{
  m_bam = hts_open(bam.c_str(), "rb");
  ASSERT(m_bam, "Could not load bam file: " + bam + "\nCheck that file exists and is in right format!");

  m_bam_index = sam_index_load(m_bam, bam.c_str());
  ASSERT(m_bam_index, "Could not load bam index file: " + bam + "\nCheck to be sure " + bam + ".bai exists!");

  m_bam_header = sam_hdr_read(m_bam);
  assert(m_bam_header);

	// load all the reference sequences:
  m_faidx = fai_load(fasta.c_str());
  assert(m_faidx);

  for(int i=0; i<m_bam_header->n_targets; ++i) {
		reference_sequence* refseq = new reference_sequence(m_faidx, fasta, m_bam_header->target_name[i]);
		assert(static_cast<unsigned int>(refseq->m_len) == m_bam_header->target_len[i]);
		m_refs.push_back(refseq);
	}
}


/*! Destructor.
 */
pileup_base::~pileup_base() {
  hts_close(m_bam);
  hts_idx_destroy(m_bam_index);
  sam_hdr_destroy(m_bam_header);
  fai_destroy(m_faidx);
  
  // Clean up vector of pointers
  while(!m_refs.empty())
  {
	  if (m_refs.back() != NULL)
		  delete m_refs.back();
	  m_refs.pop_back();
  }
}


/*! Retrieve the reference sequence for the given target and fasta index.
 */
char* pileup_base::get_refseq(uint32_t target) const {
	assert(static_cast<uint32_t>(target) < m_refs.size());
	return m_refs[target]->m_seq;
}

/*! Guards for whether to handle a position during pileup
 */
bool pileup_base::handle_position(uint32_t pos_1) {

  // Print progress (1-indexed position)
  if(m_print_progress) {
    uint32_t target_total = target_length(m_last_tid);
    if ((pos_1 % 10000 == 0) || (pos_1 == target_total)) {
      ostringstream progress_message;
      progress_message << "    POSITION:" << setw(12) << right << pos_1 << "/" << target_total;
      print_progress_line(progress_message.str());
    }
  }
  
  if ( m_clip_start_position_1 && (pos_1 < m_clip_start_position_1) ) return false;
  if ( m_clip_end_position_1   && (pos_1 > m_clip_end_position_1  ) ) return false;
  if ( m_downsample && ( (pos_1 + m_start_position_1) % m_downsample != 0 ) ) return false;

	return true;
}


/*! First-level callback, used to re-route the callback from samtools to the virtual
 function defined in breseq::pileup_base.
 */
int first_level_pileup_callback(uint32_t tid, uint32_t pos, int32_t n, const bam_pileup1_t *pile, void *data) {
  bool debug = false;
  
  if (debug) cout << "pileup: " << tid << ", pos: " << pos << "n" << endl;
  
	pileup_base* pb = reinterpret_cast<pileup_base*>(data);
  
  /////////////////////////
  // We've changed targets (and may have skipped a target entirely)

  // handle special case at the beginning
  if (pb->m_last_tid == UNDEFINED_UINT32) {
    pb->m_last_tid = 0;
    pb->at_target_start_first_level_callback(pb->m_last_tid); 
  }

  
  while(pb->m_last_tid != tid) {
    
    // 1) Finish all positions in the last target 
    for (uint32_t on_pos_1 = pb->m_last_position_1+1; on_pos_1 <= pb->target_length(pb->m_last_tid); on_pos_1++) {
      
      // Missing positions
      if ( pb->handle_position(on_pos_1) ) {
        pileup p(pb->m_last_tid,on_pos_1,0,NULL,*pb);
        pb->pileup_callback(p);
      }
      
      pb->m_last_position_1 = on_pos_1;
    }
    
    // 2) Finish the last target
    if (tid > 0) pb->at_target_end_first_level_callback(pb->m_last_tid);
    
    // 3) Move to next target and reset position to zero
    pb->m_last_tid++;
    pb->m_last_position_1 = 0;
    
    // Start a new target if necessary
    if (pb->m_last_position_1 == 0) {
      pb->at_target_start_first_level_callback(pb->m_last_tid); 
    }
  }
  
  uint32_t this_pos_1 = pos+1;

  //////////////////////////
  // Handle missing positions between last handled and current on same target id
  for (uint32_t on_pos_1 = pb->m_last_position_1+1; on_pos_1 <= this_pos_1-1; on_pos_1++) {
    
    // Missing positions
    if ( pb->handle_position(on_pos_1) ) {
      pileup p(tid,on_pos_1,0,NULL,*pb);
      pb->pileup_callback(p);
    }
    pb->m_last_position_1 = on_pos_1;
  }


  //////////////////////////
  // Handle current position
        
  if ( pb->handle_position(this_pos_1) ) {
    pileup p(tid,this_pos_1,n,pile,*pb);
    pb->pileup_callback(p);
  } 
  pb->m_last_position_1 = this_pos_1;

	return 0;
}

/*! First-level callback, used to re-route the callback from samtools to the virtual
 function defined in breseq::pileup_base.
 */
int first_level_fetch_callback(const bam1_t *b, void *data)
{
 	pileup_base* pb = reinterpret_cast<pileup_base*>(data);
  alignment_wrapper a(b);
  pb->fetch_callback(a);
  
  return 0;
}


/*! Helper struct and read function for whole-file pileup via bam_plp_init.
 */
struct PileupReadData {
  samFile   *fp;
  bam_hdr_t *hdr;
};

static int pileup_read_func(void *data, bam1_t *b) {
  PileupReadData *d = reinterpret_cast<PileupReadData*>(data);
  return sam_read1(d->fp, d->hdr, b);
}

/*! Run the pileup.
 */
void pileup_base::do_pileup() {

  // Start the first target...
  m_last_position_1 = 0;
  m_downsample = 0;
  m_start_position_1 = 0;
  m_end_position_1 = 0;
  m_clip_start_position_1 = 0;
  m_clip_end_position_1 = 0;
  m_last_tid = UNDEFINED_UINT32;

  // Rewind to the start of alignments so we read the whole file
  // (hts_set_opt or re-open would be alternatives; sam_read1 reads sequentially)
  PileupReadData rd = {m_bam, m_bam_header};
  bam_plp_t plp = bam_plp_init(pileup_read_func, &rd);
  // Set maxcnt to 1e9 to match the bundled modified samtools behaviour
  // (original bundled code patched sam.c: iter->maxcnt = 1000000000)
  bam_plp_set_maxcnt(plp, 1000000000);

  int tid, pos, n;
  const bam_pileup1_t *pile;
  while ((pile = bam_plp_auto(plp, &tid, &pos, &n)) != NULL) {
    first_level_pileup_callback((uint32_t)tid, (uint32_t)pos, n, pile, this);
  }
  bam_plp_destroy(plp);

  // Handle positions after the last one handled by the pileup.
  // => This includes target id's that might not have been handled at all...
  for (uint32_t t = m_last_tid; t < num_targets(); t++) {

    // We need to start this target
    if (m_last_position_1 == 0) {
      at_target_start_first_level_callback(m_last_tid);
    }

    for (uint32_t on_pos_1 = m_last_position_1+1; on_pos_1 <= m_bam_header->target_len[t]; on_pos_1++) {
      if ( handle_position(on_pos_1) ) {
        pileup p(t, on_pos_1, 0, NULL, *this);
        pileup_callback(p);
        m_last_position_1 = on_pos_1;
      }
    }

    // We finished this target
    at_target_end_first_level_callback(m_last_tid);
    m_last_tid++;
    m_last_position_1 = 0;
  }
}


/*! Helper struct and read function for region-based pileup via bam_plp_init.
 */
struct IterReadData {
  samFile    *fp;
  bam_hdr_t  *hdr;
  hts_itr_t  *iter;
};

static int iter_read_func(void *data, bam1_t *b) {
  IterReadData *d = reinterpret_cast<IterReadData*>(data);
  return sam_itr_next(d->fp, d->iter, b);
}

/*! Run the pileup.

    !clip means that we handle all alignment positions that reads overlapping this position overlap
    clip means we stop and end at the indicated alignment columns, even if reads extend past them
 */
void pileup_base::do_pileup(const string& region, bool clip, uint32_t downsample) {

  uint32_t target_id, start_pos_1, end_pos_1, insert_start, insert_end;
  parse_region(region.c_str(), target_id, start_pos_1, end_pos_1, insert_start, insert_end);

  m_downsample = downsample; // init
  m_start_position_1 = start_pos_1;
  m_end_position_1 = end_pos_1;
  m_insert_start = insert_start;
  m_insert_end = insert_end;

  m_last_tid = target_id;
  m_last_position_1 = 0; // reset
  m_clip_start_position_1 = 0; // reset
  m_clip_end_position_1 = 0; // reset

  at_target_start_first_level_callback(target_id);

  if (clip) {
    m_clip_start_position_1 = start_pos_1;
    m_clip_end_position_1 = end_pos_1;
    assert(m_clip_start_position_1 > 0); // prevent underflow of unsigned
  }

  // sam_itr_queryi uses 0-indexed half-open [start, end)
  hts_itr_t *iter = sam_itr_queryi(m_bam_index, (int)target_id,
                                    (int)(start_pos_1 - 1), (int)end_pos_1);
  ASSERT(iter, "Could not create iterator for region: " + region);

  IterReadData ird = {m_bam, m_bam_header, iter};
  bam_plp_t plp = bam_plp_init(iter_read_func, &ird);
  bam_plp_set_maxcnt(plp, 1000000000);

  int tid2, pos2, n;
  const bam_pileup1_t *pile;
  while ((pile = bam_plp_auto(plp, &tid2, &pos2, &n)) != NULL) {
    first_level_pileup_callback((uint32_t)tid2, (uint32_t)pos2, n, pile, this);
  }
  bam_plp_destroy(plp);
  hts_itr_destroy(iter);

  // handle positions after the last one handled by the pileup
  // unlike the case above, we only need to finish this target
  for (uint32_t on_pos_1 = m_last_position_1+1; on_pos_1 <= m_end_position_1; on_pos_1++) {
    if ( handle_position(on_pos_1) ) {
      pileup p(m_last_tid, on_pos_1, 0, NULL, *this);
      pileup_callback(p);
      m_last_position_1 = on_pos_1;
    }
  }
  at_target_end_first_level_callback(target_id);
}
  
/*! Version of pileup that only processes certain seq_ids
 */   
  
void pileup_base::do_pileup(const set<string>& seq_ids) {
  
  for(set<string>::iterator it=seq_ids.begin(); it!=seq_ids.end(); it++) {
        
    // We need to find the target id from the name, this is a rather inefficient way to do it.
    const string& this_seq_id = *it;
    uint32_t this_tid=0;
    bool found = false;
    while( !found && (this_tid < num_targets()) ) {
      if (target_name(this_tid) == this_seq_id) {
        found = true;
      }
      else {
        this_tid++;
      }
    }
    ASSERT(found, "Could not find seq_id: " + this_seq_id); 
    string region_str = this_seq_id + ":1-" + to_string(target_length(this_tid)); 
                                                        
    this->do_pileup(region_str);
  }
}

  
void pileup_base::do_fetch(const string& region) {

  uint32_t target_id, start_pos_1, end_pos_1, insert_start, insert_end;
  parse_region(region.c_str(), target_id, start_pos_1, end_pos_1, insert_start, insert_end);
  m_insert_start = insert_start;
  m_insert_end = insert_end;

  hts_itr_t *iter = sam_itr_queryi(m_bam_index, (int)target_id,
                                    (int)(start_pos_1 - 1), (int)end_pos_1);
  ASSERT(iter, "Could not create iterator for region: " + region);

  bam1_t *b = bam_init1();
  while (sam_itr_next(m_bam, iter, b) >= 0) {
    alignment_wrapper a(b);
    fetch_callback(a);
  }
  bam_destroy1(b);
  hts_itr_destroy(iter);
}

} //end namespace breseq
