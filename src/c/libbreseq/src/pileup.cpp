#include <assert.h>
#include <iostream>
#include <sstream>

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
char* breseq::pileup_base::get_refseq(int target) {
	assert(static_cast<std::size_t>(target) < _refs.size());
	return _refs[target]->_seq;
}


/*! First-level callback, used to re-route the callback from samtools to the virtual
 function defined in breseq::pileup.
 */
int first_level_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data) {
	breseq::pileup_base* p = reinterpret_cast<breseq::pileup_base*>(data);
	
	if((pos % 10000) == 0) {
		std::cerr << "    POSITION:" << pos << std::endl;
	}
	
	p->callback(tid,pos,n,pile);
	return 0;
}

/*! Run the pileup.
 */
void breseq::pileup_base::pileup() {
	sampileup(_bam, BAM_DEF_MASK, first_level_callback, this);
}
