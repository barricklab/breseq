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


using namespace std;

namespace breseq {

const char alignment::op_to_char[10] = "MIDNSHP=X";
  
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
: _p(NULL), _a(a) 
{
}
  
alignment::alignment()
: _p(NULL), _a(NULL) 
{
}
  
alignment::~alignment()
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

 Replaced with cigar_query_length()
uint32_t alignment::query_length() const {
	uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment
	uint32_t qlen = bam_cigar2qlen(&_a->core, cigar); // total length of the query
	return qlen;
}
*/

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
		if((read_length()-(query_position_1())) <= (uint32_t)bam_aux2i(auxr)) {
			return true;
		}
	}
	
	return false;
}
//! Return number of locations on left of sequence to be trimmed
uint32_t alignment::trim_left() const {
  uint8_t *auxl = bam_aux_get(_a,"XL");
  if(auxl)
  {
    return (uint32_t)bam_aux2i(auxl); 
  }
  else
  {
    return 0;
  }
}
//! Return number of locations on right of sequence to be trimmed
uint32_t alignment::trim_right() const {
  uint8_t *auxr = bam_aux_get(_a,"XR");
  if(auxr)
  {
    return (uint32_t)bam_aux2i(auxr); 
  }
  else 
  {
    return 0;
  }
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
  pair<uint32_t,uint32_t> qb = query_bounds_1();
  qb.first--;
  qb.second--;
  return qb;
}
void alignment::query_bounds_0(uint32_t& start, uint32_t& end) const {
	pair<uint32_t,uint32_t> qb = query_bounds_0();
	start = qb.first;
	end = qb.second;
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
void alignment::query_bounds_1(uint32_t& start, uint32_t& end) const {
	pair<uint32_t,uint32_t> qb = query_bounds_1();
	start = (int32_t)qb.first;
	end = (int32_t)qb.second;
}
  
/*! Retrieve the start and end coordinates of the aligned part of the read.
    switch start and end if on opposite reference strand
 */
std::pair<uint32_t,uint32_t> alignment::query_stranded_bounds_1() const {
  
  pair<uint32_t,uint32_t> qb = query_bounds_1();
  uint32_t start = (int32_t)qb.first;
  uint32_t end = (int32_t)qb.second;
  
  if (reversed()) {
    return std::make_pair(read_length() - end + 1, read_length() - start + 1);
  }
  
  return std::make_pair(start,end);
}
void alignment::query_stranded_bounds_1(uint32_t& start, uint32_t& end) const {
  pair<uint32_t,uint32_t> qb = query_stranded_bounds_1();
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
  
  for(uint32_t j=(_a->core.n_cigar-1); j>0; j--) {
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
  
/*uint32_t alignment::reference_end_0() const {

  uint32_t pos = reference_start_0();
  
  uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment
	
  for(uint32_t j=0; j<_a->core.n_cigar; j++) {
    uint32_t op = cigar[j] & BAM_CIGAR_MASK;
    uint32_t len = cigar[j] >> BAM_CIGAR_SHIFT;
    
    // if we encounter padding, or a gap in reference then we are done
    if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
      pos += len;
    }
  }
  pos -= 1; // to get inclusive coords

  return pos;
}*/



uint32_t alignment::base_repeat_0(uint32_t q_pos_0) const {

  uint8_t this_base_bam = read_base_bam_0(q_pos_0);
  uint32_t base_repeat = 0;
  if (!reversed()) {
    while ( q_pos_0 < query_end_0()) {
      q_pos_0++;
      if (this_base_bam != read_base_bam_0(q_pos_0)) break;
      base_repeat++;
    }  
  } else {
    while (q_pos_0 > 0) {
      q_pos_0--;
      if (this_base_bam != read_base_bam_0(q_pos_0)) break;
      base_repeat++;
    }    
  }
  
  return base_repeat;
}
  
tam_file::tam_file(const string& tam_file_name, const string& fasta_file_name, ios_base::openmode mode) 
  : bam_header(NULL), input_tam(NULL)
{
  if (mode == ios_base::in)
  {
    open_read(tam_file_name, fasta_file_name);
  }
  else
  {
    open_write(tam_file_name, fasta_file_name);
  }
}
  
tam_file::~tam_file()
{
  sam_close(input_tam);
  if (bam_header) bam_header_destroy(bam_header);
  free_loaded_alignments();
}
  
void tam_file::free_loaded_alignments()
{
  for(vector<bam1_t*>::iterator it=loaded_alignments.begin(); it<loaded_alignments.end(); it++)
  {
    bam_destroy1(*it);
  }
  loaded_alignments.clear();
}

  
void tam_file::open_read(const string& tam_file_name, const string& fasta_file_name)
{
  string faidx_file_name(fasta_file_name);
  faidx_file_name += ".fai";
  
  input_tam = sam_open(tam_file_name.c_str()); // or die("Could not open reference same file: $reference_sam_file_name");
  assert(input_tam);
  bam_header = sam_header_read2(faidx_file_name.c_str()); // or die("Error reading reference fasta index file: $reference_faidx_file_name");
}

void tam_file::open_write(const string& tam_file_name, const string& fasta_file_name)
{
  string faidx_file_name(fasta_file_name);
  faidx_file_name += ".fai";

  output_tam.open(tam_file_name.c_str(), ios_base::out);
  assert(output_tam.is_open());
  bam_header = sam_header_read2(faidx_file_name.c_str());
}

bool tam_file::read_alignments(alignment_list& alignments, bool paired)
{
  alignments.clear();
  
	int num_to_slurp = (paired) ? 2 : 1;
	string last_read_name;
	if (loaded_alignments.size() > 0)
	{
    alignment last_alignment(loaded_alignments.back());
		last_read_name = last_alignment.read_name();
    
    loaded_alignments.resize(loaded_alignments.size()-1);
		alignments.push_back(last_alignment);
	}

	int num_slurped = 0;
	while (true)
	{
		bam1_t* last_alignment_bam = new bam1_t;
    last_alignment_bam = bam_init1();
    int bytes = sam_read1(input_tam,bam_header,last_alignment_bam);
    if (bytes < 0) break;
    
    loaded_alignments.push_back(last_alignment_bam);
    
    alignment last_alignment(last_alignment_bam);

		string read_name = last_alignment.read_name();
      
		if (last_read_name.size() == 0)
    {
			last_read_name = read_name;
    }
    else
    {
      if (read_name != last_read_name) break;
    }
    
    alignments.push_back(last_alignment);
  }

  return (alignments.size() > 0);
}
  
void tam_file::write_alignments(int32_t fastq_file_index, alignment_list& alignments, vector<Trim>* trims)
{
	for (uint32_t i = 0; i < alignments.size(); i++)
	{
		alignment& a = alignments[i];

		stringstream aux_tags_ss;
    
		aux_tags_ss << "AS:i:" << a.aux_get_i("AS") << "\t" << "X1:i:" << alignments.size() << "\t" << "X2:i:" << fastq_file_index;

		if (trims != NULL && trims->size() > i)
		{
			Trim trim = (*trims)[i];
			aux_tags_ss << "\t" << "XL:i:" << trim.L << "\t" << "XR:i:" << trim.R;
		}

		string aux_tags = aux_tags_ss.str();

		uint8_t* qscore = a.read_base_quality_sequence();
    stringstream quality_score_ss;

		for (uint32_t j = 0; j < a.read_length(); j++)
    {
			quality_score_ss << static_cast<char>(*qscore + 33);
      qscore++;
    }
    string quality_score_string = quality_score_ss.str();
    
		string cigar_string = a.cigar_string();

		vector<string> ll;
		ll.push_back(a.read_name());
    //cerr << a.read_name() << endl;
		ll.push_back(to_string(fix_flags(a.flag())));
		ll.push_back(bam_header->target_name[a.reference_target_id()]);
		ll.push_back(to_string(a.reference_start_1()));
		ll.push_back(to_string<uint32_t>(a.quality()));
		ll.push_back(cigar_string);

		//part of a pair?
		if ((a.flag() & BAM_FPROPER_PAIR) == 0)
		{
			ll.push_back("*");
			ll.push_back("0");
			ll.push_back("0");
		}
		else
		{
			ll.push_back("=");
			ll.push_back(to_string<int32_t>(a.mate_start_1()));
			ll.push_back(to_string<int32_t>(a.isize()));
		}

    
		ll.push_back(a.read_char_sequence());
		ll.push_back(quality_score_string);
		ll.push_back(aux_tags);

		output_tam << join(ll, "\t") << endl;
	}
}
  
void tam_file::write_split_alignment(uint32_t min_indel_split_len, const alignment& a)
{


}

}
