#include <assert.h>
#include <iostream>
#include <sstream>
#include "pileup.h"


breseq::reference_sequence::reference_sequence(const std::string& fasta_filename, const std::string& region_name)
: _ref(0), _seq(0), _len(0) {
	_ref = fai_load(fasta_filename.c_str());
	assert(_ref);

	_seq = fai_fetch(_ref, region_name.c_str(), &_len);	
	assert(_seq);
	assert(_len > 0);
}


breseq::reference_sequence::~reference_sequence() {
	if(_ref) { fai_destroy(_ref);	}
	if(_seq) { free(_seq); }
}


/*! Constructor for single-BAM, single-FASTA error counting.
 */
breseq::pileup_base::pileup_base(const std::string& bam, const std::vector<std::string>& fastas)
: _bam(0) {
	using namespace std;
	_bam = samopen(bam.c_str(), "rb", 0);
	assert(_bam);
	assert(_bam->header->n_targets == static_cast<int>(fastas.size()));
	
	// load the reference sequence for each target:
	for(int i=0; i<_bam->header->n_targets; ++i) {
		cerr << "  REFERENCE: " << _bam->header->target_name[i] << endl;
		cerr << "  LENGTH: " << _bam->header->target_len[i] << endl;
		boost::shared_ptr<reference_sequence> refseq(new reference_sequence(fastas[i], _bam->header->target_name[i]));
		_ref_seqs.push_back(refseq);
	}
}


/*! Destructor.
 */
breseq::pileup_base::~pileup_base() {
	samclose(_bam);
	// reference sequence destructors called automatically.
}


/*! First-level callback, used to re-route the callback from samtools to the virtual
 function defined in breseq::pileup.
 */
int breseq::first_level_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data) {
	breseq::pileup_base* p = reinterpret_cast<breseq::pileup_base*>(data);
	
	if((pos % 10000) == 0) {
		std::cerr << "    POSITION:" << pos << std::endl;
	}
	
	assert(tid < p->_ref_seqs.size());
	assert(pos < static_cast<uint32_t>(p->_ref_seqs[tid]->_len));
	p->callback(p->_ref_seqs[tid]->_seq,tid,pos,n,pile);
	return 0;
}


/*! Run the pileup.
 */
void breseq::pileup_base::pileup() {
	sampileup(_bam, -1, first_level_callback, this);
}
