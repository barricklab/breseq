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


#include "libbreseq/pileup_base.h"
#include "libbreseq/pileup.h"

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

  (void)fasta_filename; //TODO: unused??
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
: m_bam(0), m_bam_header(0), m_bam_index(0), m_bam_file(0), m_faidx(0), m_bam_file_name(bam), m_fasta_file_name(fasta), 
  m_last_position_1(0), m_start_position_1(0), m_end_position_1(0), m_clip_start_position_1(0), m_clip_end_position_1(0), 
  m_downsample(0), m_last_tid(static_cast<uint32_t>(-1)), m_print_progress(false)
{
	m_bam = samopen(bam.c_str(), "rb", 0);
  ASSERT(m_bam, "Could not load bam file: " + bam + "\nCheck that file exists and is in right format!");
  
  //m_bam_file = m_bam->x.bam;
  
  m_bam_file = bam_open(bam.c_str(), "r");
  ASSERT(m_bam_file, "Could not load bam file: " + bam + "\nCheck that file exists and is in right format!");
  
  m_bam_index = bam_index_load(bam.c_str());
  ASSERT(m_bam_index, "Could not load bam index file: " + bam + "\nCheck to be sure " + bam + ".bai exists!");

  m_bam_header = bam_header_read(m_bam_file);
  assert(m_bam_header);
  
	// load all the reference sequences:
  m_faidx = fai_load(fasta.c_str());
  assert(m_faidx);
  
  for(int i=0; i<m_bam->header->n_targets; ++i) {
		reference_sequence* refseq = new reference_sequence(m_faidx, fasta, m_bam->header->target_name[i]);
		assert(static_cast<unsigned int>(refseq->m_len) == m_bam->header->target_len[i]);
		m_refs.push_back(refseq);
	}
}


/*! Destructor.
 */
pileup_base::~pileup_base() {
	samclose(m_bam);
  bam_close(m_bam_file);
  bam_header_destroy(m_bam_header);
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
  if(m_print_progress && (pos_1 % 10000 == 0) ) {
    cerr << "    POSITION:" << pos_1 << endl;
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
  
  // Do the samtools pileup
	sampileup(m_bam, BAM_DEF_MASK, first_level_pileup_callback, this);
  
  //Handle positions after the last one handled by the pileup.
  // => This includes target id's that might not have been handled at all...

  for (uint32_t tid = m_last_tid; tid < num_targets(); tid++) {
    
    // We need to start this target
    if (m_last_position_1 == 0) {
      at_target_start_first_level_callback(m_last_tid);
    }
    
    for (uint32_t on_pos_1 = m_last_position_1+1; on_pos_1 <= m_bam->header->target_len[tid]; on_pos_1++) {
      if ( handle_position(on_pos_1) ) {
        pileup p(tid,on_pos_1,0,NULL,*this);
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


// Note: This code snippet copied from Perl Module XS Bio::Samtools
// Original Author: Lincoln Stein

int add_pileup_line (const bam1_t *b, void *data) {
  bam_plbuf_t *pileup = (bam_plbuf_t*) data;
  bam_plbuf_push(b,pileup);
  return 0;
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
  
  bam_plbuf_t        *pileup_buff;
  pileup_buff = bam_plbuf_init(first_level_pileup_callback,this);
  bam_fetch(m_bam_file,m_bam_index,target_id,start_pos_1-1,end_pos_1,(void*)pileup_buff,add_pileup_line);
  // bam_fetch expected 0 indexed start_pos and 1 indexed end_pos
  
  bam_plbuf_push(NULL,pileup_buff); // This clears out the clipped right regions... call before at_end!
  bam_plbuf_destroy(pileup_buff);
  
  // handle positions after the last one handled by the pileup
  // unlike the case above, we only need to finish this target
  
  for (uint32_t on_pos_1 = m_last_position_1+1; on_pos_1 <= m_end_position_1; on_pos_1++) {
    if ( handle_position(on_pos_1) ) {
      pileup p(m_last_tid,on_pos_1,0,NULL,*this);
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
  
  // should throw if target not found!
  bam_fetch(m_bam_file,m_bam_index,target_id,start_pos_1-1,end_pos_1,this,first_level_fetch_callback);
  // bam_fetch expected 0 indexed start_pos and 1 indexed end_pos
}

} //end namespace breseq
