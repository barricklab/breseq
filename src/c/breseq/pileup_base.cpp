/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2010 Michigan State University

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#include <assert.h>
#include <iostream>
#include <sstream>

#include "breseq/pileup_base.h"
#include "breseq/pileup.h"

/*! Constructor.
 
 Open the FASTA and read the target reference sequence.
 */
breseq::reference_sequence::reference_sequence(const std::string& fasta_filename, const std::string& target)
: _ref(0), _seq(0), _len(0) {
	_ref = fai_load(fasta_filename.c_str());
	assert(_ref);

	_seq = fai_fetch(_ref, target.c_str(), &_len);
	assert(_seq);
	assert(_len > 0);
}


/*! Destructor.
 */
breseq::reference_sequence::~reference_sequence() {
	if(_ref) { fai_destroy(_ref);	}
	if(_seq) { free(_seq); }
}


/*! Constructor for single-BAM, >=0 FASTA.
 
 \todo Change this to lazily load & cache reference sequences as they are needed,
 instead of loading all sequences at once.
 */
breseq::pileup_base::pileup_base(const std::string& bam, const std::string& fasta)
: _bam(0) {
	using namespace std;
	_bam = samopen(bam.c_str(), "rb", 0);
	assert(_bam);

	// load all the reference sequences:
	for(int i=0; i<_bam->header->n_targets; ++i) {
		cerr << "  REFERENCE: " << _bam->header->target_name[i] << endl;
		cerr << "  LENGTH: " << _bam->header->target_len[i] << endl;
		boost::shared_ptr<reference_sequence> refseq(new reference_sequence(fasta, _bam->header->target_name[i]));
		assert(static_cast<unsigned int>(refseq->_len) == _bam->header->target_len[i]);
		_refs.push_back(refseq);
	}
}


/*! Destructor.
 */
breseq::pileup_base::~pileup_base() {
	samclose(_bam);
}


/*! Retrieve the reference sequence for the given target and fasta index.
 */
char* breseq::pileup_base::get_refseq(uint32_t target) const {
	assert(static_cast<std::size_t>(target) < _refs.size());
	return _refs[target]->_seq;
}


/*! First-level callback, used to re-route the callback from samtools to the virtual
 function defined in breseq::pileup_base.
 */
int breseq::first_level_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data) {
	pileup_base* pb = reinterpret_cast<pileup_base*>(data);
		
	// if _last_tid is initialized, and is different than tid, then we've changed targets.  
	// call at_end() for the previous target:
	if(pb->_last_tid && (*pb->_last_tid != tid)) {
		pb->at_end(*pb->_last_tid, pb->_bam->header->target_len[*pb->_last_tid]);
	}

	// update _last_tid to the current tag (this is effectively a lag):
	pb->_last_tid = tid;
	
  // Print progress (1-indexed position)
	if(((pos+1) % 10000) == 0) {
		std::cerr << "    POSITION:" << (pos+1) << std::endl;
	}
	
	pileup p(tid,pos,n,pile,*pb);
	pb->callback(p);
	return 0;
}


/*! Run the pileup.
 */
void breseq::pileup_base::do_pileup() {
	sampileup(_bam, BAM_DEF_MASK, first_level_callback, this);
	at_end(*_last_tid, _bam->header->target_len[*_last_tid]);
}
