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

#include "breseq/alignment.h"

namespace breseq {

/*! Constructor.
 */
alignment::alignment(const bam_pileup1_t* p)
: _p(p)
, _a(p->b) {
}
  
/*! Constructor.
 only alignment!
 */
alignment::alignment(const bam1_t* a)
: _a(a)
{
}

/*! Does this alignment have any redundancies?
 */
bool alignment::is_redundant() const {
	return (redundancy() > 1);
}


/*! Number of redundancies at this alignment.
 */
uint32_t alignment::redundancy() const {
	return bam_aux2i(bam_aux_get(_a,"X1"));
}


/*! Calculate the length of this query out to the last non-clip, non-skip.
 */
uint32_t alignment::query_length() const {
	uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment
	uint32_t qlen = bam_cigar2qlen(&_a->core, cigar); // total length of the query
	return qlen;
}


/*! Retrieve the index of the read file that contained this alignment
 */
uint32_t alignment::fastq_file_index() const {
	return bam_aux2i(bam_aux_get(_a,"X2"));
}

//! Has this alignment been trimmed?
bool alignment::is_trimmed() const {
	// is our query position in the left-side trimmed region?
	uint8_t *auxl = bam_aux_get(_a,"XL");
	if(auxl) {
		if((query_position_1()) <= (uint32_t)bam_aux2i(auxl)) {
			return true;
		}
	}
	
	// is our query position in the right-side trimmed region?
	uint8_t *auxr = bam_aux_get(_a,"XR");
	if(auxr) {
		if((query_length()-(query_position_1())) <= (uint32_t)bam_aux2i(auxr)) {
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


std::pair<uint32_t,uint32_t> alignment::query_bounds_0() const {
  pair<int32_t,int32_t> qb = query_bounds_1();
  qb.first--;
  qb.second--;
  return qb;
}
void alignment::query_bounds_0(int32_t& start, int32_t& end) const {
	pair<int32_t,int32_t> qb = query_bounds_0();
	start = (int32_t)qb.first;
	end = (int32_t)qb.second;
}

/*! Retrieve the start and end coordinates of the aligned part of the read.
 */
std::pair<uint32_t,uint32_t> alignment::query_bounds_1() const {
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
void alignment::query_bounds_1(int32_t& start, int32_t& end) const {
	pair<uint32_t,uint32_t> qb = query_bounds_1();
	start = (int32_t)qb.first;
	end = (int32_t)qb.second;
}

/*! Get the query start or end from the cigar string of an alignment
 */
uint32_t alignment::query_start_1() const {
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


uint32_t alignment::query_end_1() const {
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
  
uint32_t alignment::reference_end_0() const {

	//TODO: Can you just do "return bam_calend(&_a->core, bam1_cigar(_a));"?
  uint32_t pos = reference_start_0();
  
  uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment
	
  for(uint32_t j=(_a->core.n_cigar-1); j>0; --j) {
    uint32_t op = cigar[j] & BAM_CIGAR_MASK;
    uint32_t len = cigar[j] >> BAM_CIGAR_SHIFT;
    
    // if we encounter padding, or a gap in reference then we are done
    if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
      pos += len;
    }
  }

  return pos;
}



uint32_t alignment::base_repeat_0(uint32_t q_pos_0) const {

  uint8_t this_base_bam = query_base_bam_0(q_pos_0);
  uint32_t base_repeat = 0;
  if (!reversed()) {
    while ( q_pos_0 < query_end_0()) {
      q_pos_0++;
      if (this_base_bam != query_base_bam_0(q_pos_0)) break;
      base_repeat++;
    }  
  } else {
    while (q_pos_0 > 0) {
      q_pos_0--;
      if (this_base_bam != query_base_bam_0(q_pos_0)) break;
      base_repeat++;
    }    
  }
  
  return base_repeat;
}


vector<alignment> alignment::tam_next_read_alignments(tamFile tam, bam_header_t* header, alignment* last_alignment, bool paired)
{
	int num_to_slurp = (paired) ? 2 : 1;
	string last_read_name;
	vector<alignment> al_ref;
	if (last_alignment != NULL)
	{
		last_read_name = last_alignment->query_name();
		al_ref.push_back(*last_alignment);
	}

	int num_slurped = 0;
	while (true)
	{
		bam1_t* last_alignment_bam;
		int bytes = sam_read1(tam, header, last_alignment_bam);
		last_alignment = new alignment(last_alignment_bam);

		//returns bytes == -1 if EOF reached
		if (bytes < 0)
		{
			last_alignment = NULL;
			return al_ref;
		}

		string read_name = last_alignment->query_name();

		if ( (last_read_name.size() > 0) && (read_name != last_read_name) && (++num_slurped == num_to_slurp) )
			break;

		if (last_read_name.size() == 0)
			last_read_name = read_name;

		al_ref.push_back(*last_alignment);
	}

	return al_ref;
}

void alignment::tam_write_read_alignments(ofstream& fh, bam_header_t* header, int32_t fastq_file_index, vector<alignment> al, vector<Trim>* trims)
{
	for (int32_t i = 0; i < al.size(); i++)
	{
		alignment a = al[i];

		stringstream aux_tags_ss;
		aux_tags_ss << "AS:i:" << a.aux_get("AS") << "\t" << "X1:i:" << al.size() << "\t" << "X2:i:" << fastq_file_index;

		if (trims != NULL && trims->size() > i)
		{
			Trim trim = (*trims)[i];
			aux_tags_ss << "\t" << "XL:i:" << trim.L << "\t" << "XR:i:" << trim.R;
		}

		string aux_tags = aux_tags_ss.str();

		string* qscore = (string*)(a.quality_scores());
		string quality_score_string = (qscore == NULL) ? "" : *qscore;
		for (int32_t j = 0; j < quality_score_string.size(); j++)
			quality_score_string[j] = quality_score_string[j] + 33;

		uint32_t* cigar_list = a.cigar_array();
		stringstream cigar_string_ss;

		for (int32_t j = 0; j <= a.cigar_array_length(); j++) //foreach my $c (@$cigar_list)
		{
			uint32_t op = cigar_list[i] & BAM_CIGAR_MASK;
			uint32_t len = cigar_list[i] >> BAM_CIGAR_SHIFT;
			cigar_string_ss << len << op; //$cigar_string += $c->[1] + $c->[0];
		}
		string cigar_string = cigar_string_ss.str();

		vector<string> ll;
		ll.push_back(a.query_name());
		ll.push_back(boost::lexical_cast<string>(fix_flags(a.flag())));
		ll.push_back(header->target_name[a.reference_target_id()]);
		ll.push_back(boost::lexical_cast<string>(a.reference_start_0()));
		ll.push_back(boost::lexical_cast<string>(a.quality()));
		ll.push_back(cigar_string);

		//something strange in new version... such that mate_start sometimes
		//returns 1 even though there is no mate
		if (a.flag() & BAM_FPROPER_PAIR != 0)
		{
			ll.push_back("*");
			ll.push_back(0);
			ll.push_back(0);
		}
		else
		{
			ll.push_back("=");
			ll.push_back(boost::lexical_cast<string>(a.mate_start_1()));
			ll.push_back(boost::lexical_cast<string>(a.isize()));
		}

		ll.push_back(boost::lexical_cast<string>(*(a.query_bam_sequence())));
		ll.push_back(quality_score_string);
		ll.push_back(aux_tags);

		fh << join(ll, "\t") << endl;
	}
}

}

