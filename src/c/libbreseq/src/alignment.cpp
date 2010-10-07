#include "breseq/alignment.h"


/*! Constructor.
 */
breseq::alignment::alignment(const bam_pileup1_t* p)
: _p(p)
, _a(p->b) {
}


/*! Does this alignment have any redundancies?
 */
bool breseq::alignment::is_redundant() const {
	return (redundancy() > 1);
}


/*! Number of redundancies at this alignment.
 */
int32_t breseq::alignment::redundancy() const {
	return bam_aux2i(bam_aux_get(_a,"X1"));
}


/*! Calculate the length of this query out to the last non-clip, non-skip.
 */
int32_t breseq::alignment::query_length() const {
	uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment
	int32_t qlen = bam_cigar2qlen(&_a->core, cigar); // total length of the query
	return qlen;
}


/*! Retrieve the index of the read file that contained this alignment
 */
int32_t breseq::alignment::fastq_file_index() const {
	return bam_aux2i(bam_aux_get(_a,"X2"));
}

//! Has this alignment been trimmed?
bool breseq::alignment::is_trimmed() const {
	// is our query position in the left-side trimmed region?
	uint8_t *auxl = bam_aux_get(_a,"XL");
	if(auxl) {
		if((query_position()+1) <= bam_aux2i(auxl)) {
			return true;
		}
	}
	
	// is our query position in the right-side trimmed region?
	uint8_t *auxr = bam_aux_get(_a,"XR");
	if(auxr) {
		if((query_length()-(query_position()+1)) <= bam_aux2i(auxr)) {
			return true;
		}
	}
	
	return false;
}

//			my $trimmed = 0;
//			my $trim_left = $a->aux_get('XL');  
//			my $trim_right = $a->aux_get('XR');
//			$trimmed = 1 if ((defined $trim_left) && ($p->qpos+1 <= $trim_left));
//			$trimmed = 1 if ((defined $trim_right) && ($a->query->length-$p->qpos <= $trim_right));
//			
//##also trimmed if up to next position and this is an indel.
//			if ($indel != 0)
//			{
//#remember that qpos is 0-indexed
//				$trimmed = 1 if ((defined $trim_left) && ($p->qpos+1 <= $trim_left)); #so this is position <1
//				$trimmed = 1 if ((defined $trim_right) && ($a->query->length-($p->qpos+1)+1 <= $trim_right)); #this is position+1
//			}


/*! Retrieve the start and end coordinates of the aligned part of the read.
 */
std::pair<int32_t,int32_t> breseq::alignment::query_bounds() const {
  uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment
	int32_t start=1, end=bam_cigar2qlen(&_a->core,cigar);
	
	// start:
  for(uint32_t i=0; i<=_a->core.n_cigar; i++) {
    uint32_t op = cigar[i] & BAM_CIGAR_MASK;
    uint32_t len = cigar[i] >> BAM_CIGAR_SHIFT;
    // if we encounter padding, or a gap in reference then we are done
    if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
			break;
    }
    start += len;
  }
	
	// end:
  for(uint32_t i=(_a->core.n_cigar-1); i>0; --i) {
    uint32_t op = cigar[i] & BAM_CIGAR_MASK;
    uint32_t len = cigar[i] >> BAM_CIGAR_SHIFT;    
    // if we encounter padding, or a gap in reference then we are done
    if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
      break;
    }
    end -= len;
  }

	return std::make_pair(start,end);
}


/*! Get the query start or end from the cigar string of an alignment
 */
int32_t breseq::alignment::query_start() const {
  // traverse the cigar array
  uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment
  int32_t pos = 1;
	
  for(uint32_t j=0; j<=_a->core.n_cigar; j++) {
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


int32_t breseq::alignment::query_end() const {
  // traverse the cigar array
  uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment
  int32_t pos = bam_cigar2qlen(&_a->core, cigar); // total length of the query
	
  for(uint32_t j=(_a->core.n_cigar-1); j>0; --j) {
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
