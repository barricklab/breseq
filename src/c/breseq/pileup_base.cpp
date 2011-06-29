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


#include "breseq/pileup_base.h"
#include "breseq/pileup.h"

using namespace std;

namespace breseq {
  
/*! Constructor.
 
 Open the FASTA and read the target reference sequence.
 */
reference_sequence::reference_sequence(const std::string& fasta_filename, const std::string& target)
: m_ref(0), m_seq(0), m_len(0) {
	m_ref = fai_load(fasta_filename.c_str());
	assert(m_ref);

	m_seq = fai_fetch(m_ref, target.c_str(), &m_len);
	assert(m_seq);
	assert(m_len > 0);
}


/*! Destructor.
 */
reference_sequence::~reference_sequence() {
	if(m_ref) { fai_destroy(m_ref);	}
	if(m_seq) { free(m_seq); }
}


/*! Constructor for single-BAM, >=0 FASTA.
 
 \todo Change this to lazily load & cache reference sequences as they are needed,
 instead of loading all sequences at once.
 */
pileup_base::pileup_base(const std::string& bam, const std::string& fasta)
: m_bam(0), m_bam_header(0), m_bam_index(0), m_bam_file(0), m_last_position_1(0), 
m_start_position_1(0), m_end_position_1(0), m_clip_start_position_1(0), m_clip_end_position_1(0), m_downsample(0),
m_last_tid(static_cast<uint32_t>(-1))
{
	using namespace std;
	m_bam = samopen(bam.c_str(), "rb", 0);
	assert(m_bam);

	// load all the reference sequences:
	for(int i=0; i<m_bam->header->n_targets; ++i) {
		cerr << "  REFERENCE: " << m_bam->header->target_name[i] << endl;
		cerr << "  LENGTH: " << m_bam->header->target_len[i] << endl;
		reference_sequence* refseq = new reference_sequence(fasta, m_bam->header->target_name[i]);
		assert(static_cast<unsigned int>(refseq->m_len) == m_bam->header->target_len[i]);
		m_refs.push_back(refseq);
	}
  
  m_bam_file = bam_open(bam.c_str(), "r");
  assert(m_bam_file);
  
  m_bam_index = bam_index_load(bam.c_str());
  assert(m_bam_index);

  m_bam_header = bam_header_read(m_bam_file);
  assert(m_bam_header);
}


/*! Destructor.
 */
pileup_base::~pileup_base() {
	samclose(m_bam);
  bam_close(m_bam_file);
  bam_header_destroy(m_bam_header);

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
	assert(static_cast<std::size_t>(target) < m_refs.size());
	return m_refs[target]->m_seq;
}

/*! Guards for whether to handle a position during pileup
 */
bool pileup_base::handle_position(uint32_t pos_1) {

  // Print progress (1-indexed position)
  if((pos_1 % 10000) == 0) {
    std::cerr << "    POSITION:" << pos_1 << std::endl;
  }
  
  if ( m_clip_start_position_1 && (pos_1 < m_clip_start_position_1) ) return false;
  if ( m_clip_end_position_1   && (pos_1 > m_clip_end_position_1  ) ) return false;
  if ( m_downsample && ( (pos_1 + m_start_position_1) % m_downsample != 0 ) ) return false;

	return true;
}


/*! First-level callback, used to re-route the callback from samtools to the virtual
 function defined in breseq::pileup_base.
 */
int first_level_pileup_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data) {
  
	pileup_base* pb = reinterpret_cast<pileup_base*>(data);
  
  /////////////////////////
  // We've changed targets (and may have skipped a target entirely)

  // handle special case at the beginning
  if (pb->m_last_tid == UNDEFINED) {
    pb->m_last_tid = 0;
    pb->at_target_start(pb->m_last_tid); 
  }

  
  while(pb->m_last_tid != tid) {
    
    // 1) Finish all positions in the last target 
    for (uint32_t on_pos_1 = pb->m_last_position_1+1; on_pos_1 <= pb->target_length(pb->m_last_tid); on_pos_1++) {
      
      // Missing positions
      if ( pb->handle_position(on_pos_1) ) {
        pileup p(tid,on_pos_1,0,NULL,*pb);
        pb->pileup_callback(p);
      }
      
      pb->m_last_position_1 = on_pos_1;
    }
    
    // 2) Finish the last target
    if (tid > 0) pb->at_target_end(pb->m_last_tid);
    
    // 3) Move to next target and reset position to zero
    pb->m_last_tid++;
    pb->m_last_position_1 = 0;
    
    // Start a new target if necessary
    if (pb->m_last_position_1 == 0) {
      pb->at_target_start(pb->m_last_tid); 
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
  alignment a(b);
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
  m_last_tid = UNDEFINED;
  
  // Do the samtools pileup
	sampileup(m_bam, BAM_DEF_MASK, first_level_pileup_callback, this);
  
  //Handle positions after the last one handled by the pileup.
  // => This includes target id's that might not have been handled at all...

  
  for (uint32_t tid = m_last_tid; tid < num_targets(); tid++) {
  
    // We need to start this target
    if (m_last_position_1 == 0) {
      at_target_start(m_last_tid);
    }
    
    for (uint32_t on_pos_1 = m_last_position_1+1; on_pos_1 <= m_bam->header->target_len[tid]; on_pos_1++) {
      if ( handle_position(on_pos_1) ) {
        pileup p(tid,on_pos_1,0,NULL,*this);
        pileup_callback(p);        
        m_last_position_1 = on_pos_1;
      }
    }
    
    // We finished this target
    at_target_end(m_last_tid);
    m_last_position_1 = 0;
  }
  // @JEB strictly speaking, should also handle any remaining fragments in list....
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
  
  uint32_t target_id, start_pos_1, end_pos_1;
  parse_region(region.c_str(), target_id, start_pos_1, end_pos_1); 
  // start_pos is one less than input??
  
  // should throw if target not found!

  m_downsample = downsample; // init
  m_start_position_1 = start_pos_1;
  m_end_position_1 = end_pos_1;
  
  m_last_tid = UNDEFINED;
  m_last_position_1 = 0; // reset
  m_clip_start_position_1 = 0; // reset
  m_clip_end_position_1 = 0; // reset

  if (clip) {
    m_clip_start_position_1 = start_pos_1;
    m_clip_end_position_1 = end_pos_1;
    assert(m_clip_start_position_1 > 0); // prevent underflow of unsigned

    //@JEB test removal
    //m_last_position_1 = m_clip_start_position_1-1; // So that nothing is done at start leading up to requested position
  }
    
  bam_plbuf_t        *pileup_buff;
  pileup_buff = bam_plbuf_init(first_level_pileup_callback,this);
  bam_fetch(m_bam_file,m_bam_index,target_id,start_pos_1-1,end_pos_1,(void*)pileup_buff,add_pileup_line);
  // bam_fetch expected 1 indexed start_pos and 0 indexed end_pos

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
  at_target_end(m_last_tid);
}

  
void pileup_base::do_fetch(const string& region) {
  
  uint32_t target_id;
  uint32_t start_pos_1;
  uint32_t end_pos_1;
  parse_region(region.c_str(), target_id, start_pos_1, end_pos_1);
  
  // should throw if target not found!
  //	cout << this->unique_start << endl;
  bam_fetch(m_bam_file,m_bam_index,target_id,start_pos_1-1,end_pos_1,this,first_level_fetch_callback);
  // bam_fetch expected 1 indexed start_pos and 0 indexed end_pos
	//cout << this->unique_start << endl;
}

} //end namespace breseq
