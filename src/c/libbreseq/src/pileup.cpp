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
    cerr.flush();
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

/*! Get the query start or end from the cigar string of an alignment
*/

int32_t breseq::pileup_base::query_start(const bam1_t*& a) {
  
  // traverse the cigar array
  uint32_t* cigar = bam1_cigar(a); // cigar array for this alignment
  int32_t pos = 1;

  for(uint32_t j=0; j<=a->core.n_cigar; j++) {
    uint32_t op = cigar[j] & BAM_CIGAR_MASK;
    uint32_t len = cigar[j] >> BAM_CIGAR_SHIFT;

    // if we encounter padding, or a gap in reference then we are done
    if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
        break;
    }
    
    pos += len;
  }
  
  return pos;
}

int32_t breseq::pileup_base::query_end(const bam1_t*& a) {
  
  // traverse the cigar array
  uint32_t* cigar = bam1_cigar(a); // cigar array for this alignment
  int32_t pos = bam_cigar2qlen(&a->core, cigar); // total length of the query

  for(uint32_t j=(a->core.n_cigar-1); j>0; --j) {
    uint32_t op = cigar[j] & BAM_CIGAR_MASK;
    uint32_t len = cigar[j] >> BAM_CIGAR_SHIFT;
    
    // if we encounter padding, or a gap in reference then we are done
    if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
      break;
    }
    pos -= len;
  }

  return pos;
}

/*! First-level callback, used to re-route the callback from samtools to the virtual
 function defined in breseq::pileup.
 */
int first_level_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data) {
	breseq::pileup_base* p = reinterpret_cast<breseq::pileup_base*>(data);
	
  // print position, correcting for 0-indexing versus 1-indexing
	if(((pos+1) % 10000) == 0) {
		std::cerr << "    POSITION:" << (pos+1) << std::endl;
    std::cerr.flush();
	}
	
	p->callback(tid,pos,n,pile);
	return 0;
}

/*! Run the pileup.
 */
void breseq::pileup_base::pileup() {
	sampileup(_bam, BAM_DEF_MASK, first_level_callback, this);
}
