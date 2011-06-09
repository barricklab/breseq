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

#include <assert.h>
#include <iostream>
#include <sstream>

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
m_start_position_1(0), m_end_position_1(0), m_clip_start_position_1(0), m_clip_end_position_1(0), m_downsample(0)
{
	using namespace std;
	m_bam = samopen(bam.c_str(), "rb", 0);
	assert(m_bam);

	// load all the reference sequences:
	for(int i=0; i<m_bam->header->n_targets; ++i) {
		cerr << "  REFERENCE: " << m_bam->header->target_name[i] << endl;
		cerr << "  LENGTH: " << m_bam->header->target_len[i] << endl;
		boost::shared_ptr<reference_sequence> refseq(new reference_sequence(fasta, m_bam->header->target_name[i]));
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
}


/*! Retrieve the reference sequence for the given target and fasta index.
 */
char* pileup_base::get_refseq(uint32_t target) const {
	assert(static_cast<std::size_t>(target) < m_refs.size());
	return m_refs[target]->m_seq;
}

/*! Guards for whether to handle a position during pileup
 */
bool pileup_base::handle_position(uint32_t pos) {

  if ( m_clip_start_position_1 && (pos < m_clip_start_position_1) ) return false;
  if ( m_clip_end_position_1   && (pos > m_clip_end_position_1  ) ) return false;
  if ( m_downsample && ( (pos + m_start_position_1) % m_downsample != 0 ) ) return false;

	return true;
}


/*! First-level callback, used to re-route the callback from samtools to the virtual
 function defined in breseq::pileup_base.
 */
int first_level_pileup_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data) {
  
	pileup_base* pb = reinterpret_cast<pileup_base*>(data);
		
	// if _last_tid is initialized, and is different than tid, then we've changed targets.  
	// call at_end() for the previous target:
	if(pb->m_last_tid && (*pb->m_last_tid != tid)) {
    pb->m_last_position_1 = 0;
		pb->at_end(*pb->m_last_tid, pb->m_bam->header->target_len[*pb->m_last_tid]);
	}

	// update _last_tid to the current tag (this is effectively a lag):
	pb->m_last_tid = tid;
  
	uint32_t this_pos_1 = pos+1;
  
  for (uint32_t on_pos_1 = pb->m_last_position_1+1; on_pos_1 <= this_pos_1; on_pos_1++) {
    
    // Print progress (1-indexed position)
    if((on_pos_1 % 10000) == 0) {
      std::cerr << "    POSITION:" << on_pos_1 << std::endl;
    }
    
    // Check for clipping and downsampling
    // Send empty pileup for skipped positions
    if ( pb->handle_position(on_pos_1) ) {
      
      if (on_pos_1 == this_pos_1) {
        pileup p(tid,on_pos_1-1,n,pile,*pb);
        pb->pileup_callback(p);
      } else {
        pileup p(tid,on_pos_1-1,0,NULL,*pb);
        pb->pileup_callback(p);
      }
    }
    
    // always update last position because we had an alignment...at_end expects this behavior
    pb->m_last_position_1 = on_pos_1;
    
  }

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
  m_last_position_1 = 0; // reset
  m_downsample = 0;
  m_start_position_1 = 0;
  m_end_position_1 = 0;
  m_clip_start_position_1 = 0;
  m_clip_end_position_1 = 0;
	sampileup(m_bam, BAM_DEF_MASK, first_level_pileup_callback, this);
	at_end(*m_last_tid, m_bam->header->target_len[*m_last_tid]);
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
void pileup_base::do_pileup(std::string region, bool clip, uint32_t downsample) {

  int target_id, start_pos, end_pos;
  bam_parse_region(m_bam_header, region.c_str(), &target_id, &start_pos, &end_pos); 
  // start_pos is one less than input??
  start_pos++;
  
  // should throw if target not found!

  m_downsample = downsample; // init

  m_start_position_1 = start_pos;
  m_end_position_1 = end_pos;

  m_last_position_1 = 0; // reset
  m_clip_start_position_1 = 0; // reset
  m_clip_end_position_1 = 0; // reset

  if (clip) {
    m_clip_start_position_1 = start_pos;
    m_clip_end_position_1 = end_pos;
    assert(m_clip_start_position_1 > 0); // prevent underflow of unsigned
    m_last_position_1 = m_clip_start_position_1-1; // So that nothing is done at start leading up to requested position
  }
  
  bam_plbuf_t        *pileup;
  pileup = bam_plbuf_init(first_level_pileup_callback,this);
  bam_fetch(m_bam_file,m_bam_index,target_id,start_pos,end_pos,(void*)pileup,add_pileup_line);
  bam_plbuf_push(NULL,pileup); // This clears out the clipped right regions... call before at_end!
  bam_plbuf_destroy(pileup);
  
  // Call at end with the last position we handled
  at_end(target_id, m_last_position_1);
}

  
void pileup_base::do_fetch(std::string region) {
  
  int target_id, start_pos, end_pos;
  bam_parse_region(m_bam_header, region.c_str(), &target_id, &start_pos, &end_pos); 
  
  // should throw if target not found!
//	cout << this->unique_start << endl;
  bam_fetch(m_bam_file,m_bam_index,target_id,start_pos,end_pos,this,first_level_fetch_callback);
	//cout << this->unique_start << endl;
}


} //end namespace breseq
