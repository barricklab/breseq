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

#include "libbreseq/alignment.h"

#include "libbreseq/fastq.h"
#include "libbreseq/reference_sequence.h"

using namespace std;

namespace breseq {

const char alignment_wrapper::op_to_char[10] = "MIDNSHP=X";

/*! Does this alignment have any redundancies?
 */
bool alignment_wrapper::is_redundant() const {
	return (redundancy() > 1);
}


/*! Number of redundancies at this alignment.
 */
uint32_t alignment_wrapper::redundancy() const {
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
uint32_t alignment_wrapper::fastq_file_index() const {
	return bam_aux2i(bam_aux_get(_a,"X2"));
}

//! Return number of locations on left of sequence to be trimmed
uint32_t alignment_wrapper::trim_left() const {
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
uint32_t alignment_wrapper::trim_right() const {
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


pair<uint32_t,uint32_t> alignment_wrapper::query_bounds_0(uint32_t min_qual) const {
  pair<uint32_t,uint32_t> qb = query_bounds_1(min_qual);
  if (qb.first != UNDEFINED_UINT32) qb.first--;
  if (qb.second != UNDEFINED_UINT32) qb.second--;
  return qb;
}
void alignment_wrapper::query_bounds_0(uint32_t& start, uint32_t& end, uint32_t min_qual) const {
	pair<uint32_t,uint32_t> qb = query_bounds_0(min_qual);
	start = qb.first;
	end = qb.second;
}

/*! Retrieve the start and end coordinates of the aligned part of the read.
 *  
 *  returns {-1,-1} if the entire read has lower than min_qual
 */
pair<uint32_t,uint32_t> alignment_wrapper::query_bounds_1(uint32_t min_qual) const {
  uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment
	uint32_t start=1, end=bam_cigar2qlen(&_a->core,cigar);
	
	// start:
  uint32_t i;
  uint32_t op;
  uint32_t len;
  for(i=0; i<=_a->core.n_cigar; i++) {
    op = cigar[i] & BAM_CIGAR_MASK;
    len = cigar[i] >> BAM_CIGAR_SHIFT;
    // if we encounter padding, or a gap in reference then we are done
    if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
			break;
    }
    start += len;
  }
  
  //move past low quality bases
  if (min_qual) {
    
    for(; i<=_a->core.n_cigar; i++) {
      
      if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP) && (op != BAM_CINS)) {
        for(uint32_t j=0; j<len; j++) {
          if (read_base_quality_1(start+j) > min_qual) {
            start = start + j;
            goto finish_start;
          }
        }
      }
      start += len;
      op = cigar[i] & BAM_CIGAR_MASK;
      len = cigar[i] >> BAM_CIGAR_SHIFT;
    }
    finish_start: ;
    
    if (start == read_length())
      return make_pair(UNDEFINED_UINT32,UNDEFINED_UINT32);
  }
  
	// end:
  for(i=_a->core.n_cigar-1; i>0; --i) {
    op = cigar[i] & BAM_CIGAR_MASK;
    len = cigar[i] >> BAM_CIGAR_SHIFT;    
    // if we encounter padding, or a gap in reference then we are done
    if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP) ) {
      break;
    }
    end -= len;
  }
  
  //move past low quality bases
  if (min_qual) {

    for(; i>0; --i) {
      
      if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP) && (op != BAM_CINS)) {
        for(uint32_t j=0; j<len; j++) {
          if (read_base_quality_1(end-j) > min_qual) {
            end = end-j;
            goto finish_end;
          }
        }
      }
      end -= len;
      op = cigar[i] & BAM_CIGAR_MASK;
      len = cigar[i] >> BAM_CIGAR_SHIFT;
    }
  finish_end: ;
  }

	return make_pair(start,end);
}
void alignment_wrapper::query_bounds_1(uint32_t& start, uint32_t& end, uint32_t min_qual) const {
	pair<uint32_t,uint32_t> qb = query_bounds_1(min_qual);
	start = qb.first;
	end = qb.second;
}
  
/*! Retrieve the start and end coordinates of the aligned part of the read.
    switch start and end if on opposite reference strand
 */
std::pair<uint32_t,uint32_t> alignment_wrapper::query_stranded_bounds_1() const {
  
  pair<uint32_t,uint32_t> qb = query_bounds_1();
  uint32_t start = (int32_t)qb.first;
  uint32_t end = (int32_t)qb.second;
  
  if (reversed()) {
    return std::make_pair(read_length() - end + 1, read_length() - start + 1);
  }
  
  return std::make_pair(start,end);
}
void alignment_wrapper::query_stranded_bounds_1(uint32_t& start, uint32_t& end) const {
  pair<uint32_t,uint32_t> qb = query_stranded_bounds_1();
  start = (int32_t)qb.first;
  end = (int32_t)qb.second;
}

/*! Get the query start or end from the cigar string of an alignment
 */
uint32_t alignment_wrapper::query_start_1() const {
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


uint32_t alignment_wrapper::query_end_1() const {
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
  
uint32_t alignment_wrapper::reference_start_0(uint32_t min_qual) const
{
  uint32_t start = _a->core.pos;
  
  if (min_qual) {
    
    uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment
    
    // start:
    uint32_t i;
    uint32_t op;
    uint32_t len;
    for(i=0; i<=_a->core.n_cigar; i++) {
      op = cigar[i] & BAM_CIGAR_MASK;
      len = cigar[i] >> BAM_CIGAR_SHIFT;
      // if we encounter padding, or a gap in reference then we are done
      if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
        break;
      }
    }
    
    //move past low quality bases
    for(; i<=_a->core.n_cigar; i++) {
      
      if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP) && (op != BAM_CDEL)) {
        for(uint32_t j=0; j<len; j++) {
          if (read_base_quality_1(start+j) > min_qual) {
            start = start + j;
            goto finish_start;
          }
        }
      }
      start += len;
      op = cigar[i] & BAM_CIGAR_MASK;
      len = cigar[i] >> BAM_CIGAR_SHIFT;
    }
  finish_start:
    
    // whole thing was below quality
    if (i>_a->core.n_cigar)
      return UNDEFINED_UINT32;
  }
  
  return start; 
}
  
uint32_t alignment_wrapper::reference_start_1(uint32_t min_qual) const
{
  uint32_t start = reference_start_0(min_qual);
  if (start != UNDEFINED_UINT32) start++;
  return start;
}

uint32_t alignment_wrapper::reference_end_0(uint32_t min_qual) const
{ 
  uint32_t end = reference_end_1(min_qual);
  if (end != UNDEFINED_UINT32) end--;
  return end; 
}
uint32_t alignment_wrapper::reference_end_1(uint32_t min_qual ) const
{ 
  uint32_t end = bam_calend(&_a->core, bam1_cigar(_a));
  
  if (min_qual) {
    uint32_t query_pos_1 = read_length();
    uint32_t* cigar = bam1_cigar(_a); // cigar array for this alignment

    uint32_t i;
    uint32_t op;
    uint32_t len;
    for(i=_a->core.n_cigar-1; i>0; --i) {
      op = cigar[i] & BAM_CIGAR_MASK;
      len = cigar[i] >> BAM_CIGAR_SHIFT;    
      // if we encounter padding, or a gap in reference then we are done
      if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP)) {
        break;
      }
      query_pos_1 -=len;
    }
  
  //move past low quality bases
    
    for(; i>0; --i) {
      
      if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP) && (op != BAM_CDEL)) {
        for(uint32_t j=0; j<len; j++) {
          if (read_base_quality_1(end-j) > min_qual) {
            end = end-j;
            goto finish_end;
          }
        }
        end -= len;
      }

      op = cigar[i] & BAM_CIGAR_MASK;
      len = cigar[i] >> BAM_CIGAR_SHIFT;
    }
  finish_end:

    
    // whole thing was below quality
    if (i == 0)
      return UNDEFINED_UINT32;
  }

  
  
  return end; 
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



uint32_t alignment_wrapper::base_repeat_0(uint32_t q_pos_0) const {

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
  
void alignment_wrapper::num_matches_from_end(const cReferenceSequences& ref_seq_info, bool dir, int32_t overlap, int32_t& qry_mismatch_pos, int32_t& ref_mismatch_pos)
{
  bool verbose = false;
  
  // Read
  // start, end, length  of the matching sequence (in top strand coordinates)
  uint32_t q_start, q_end;
  this->query_bounds_1(q_start, q_end);
  uint32_t q_length = q_end - q_start + 1;
  
  bool reversed = this->reversed();
  
  // sequence of the matching part of the query (top genomic strand)
  string q_str = read_char_sequence().substr(q_start - 1, q_length);
  
  // Reference
  // start, end, length of match in reference
  uint32_t r_start, r_end;
  this->reference_bounds_1(r_start, r_end);
  uint32_t r_length = r_end - r_start + 1;
  
  // sequence of match in reference (top genomic strand)
  string r_str = ref_seq_info.get_sequence_1(this->reference_target_id(), r_start, r_end);
  
  if (verbose)
  {
    cout << "====> Num Matches from End" << endl;
    cout << this->read_name() << endl;
    cout << "direction: " << dir << endl;
    cout << "Read sequence: " << this->read_char_sequence() << endl;
    cout << "Read Match coords: " << q_start << "-" << q_end << " " << reversed << endl;
    cout << "Read Match sequence: " << q_str << endl;
    cout << "Reference Match coords: " << r_start << "-" << r_end << endl;
    cout << "Reference Match sequence: " << r_str << endl;
  }
  
  if (reversed) dir = !dir;
  if (!dir)
  {
    q_str = reverse_string(q_str);
    r_str = reverse_string(r_str);
  }
  
  vector<pair<char,uint16_t> > cigar_pair_array = this->cigar_pair_array();
  
  char op_0 = cigar_pair_array.front().first;
  char op_last = cigar_pair_array.back().first;
  
  //remove soft padding
  if (op_0 == BAM_CSOFT_CLIP) cigar_pair_array.erase(cigar_pair_array.begin());
  if (op_last == BAM_CSOFT_CLIP) cigar_pair_array.erase(cigar_pair_array.end()-1);
  if (!dir) reverse(cigar_pair_array.begin(), cigar_pair_array.end());
  
  qry_mismatch_pos = -1;
  ref_mismatch_pos = -1;
  
  uint32_t positive_overlap = (overlap > 0) ? static_cast<uint32_t>(overlap) : 0;
  uint32_t r_pos = 0;
  uint32_t q_pos = 0;
  uint32_t len_0 = cigar_pair_array.front().second;
  while ((q_pos < positive_overlap) && (r_pos < r_length) && (q_pos < q_length))
  {
    // get rid of previous items...
    if (len_0 == 0) cigar_pair_array.erase(cigar_pair_array.begin());
    
    // subtract one from the length
    len_0 = cigar_pair_array.front().second;
    len_0--;
    cigar_pair_array.front().second = len_0;
    
    // handle indels
    op_0 = cigar_pair_array.front().first;
    if (op_0 == BAM_CINS)
    {
      q_pos++;
      ref_mismatch_pos = r_pos;
      qry_mismatch_pos = q_pos;
    }
    else if (op_0 == BAM_CDEL)
    {
      r_pos++;
      ref_mismatch_pos = r_pos;
      qry_mismatch_pos = q_pos;
    }
    else
    {
      if (q_str[q_pos] != r_str[r_pos])
      {
        ref_mismatch_pos = r_pos;
        qry_mismatch_pos = q_pos;
      }
      r_pos++;
      q_pos++;
    }
    
    if (verbose)
      cout << r_pos << " " << q_pos << endl;
  }
  
  // make 1 indexed...
  if (qry_mismatch_pos >= 0) qry_mismatch_pos++;
  if (ref_mismatch_pos >= 0) ref_mismatch_pos++;
  
  if (verbose)
  {
    if (qry_mismatch_pos >= 0)
      cout << "Query Mismatch At = " << qry_mismatch_pos << "Reference Mismatch At = " << ref_mismatch_pos << " | Overlap = " << overlap << endl;
    cout << "<====" << endl;
  }
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
}
  
void tam_file::open_read(const string& tam_file_name, const string& fasta_file_name)
{
  string faidx_file_name(fasta_file_name);
  faidx_file_name += ".fai";
  
  // Alternately, we could automatically generate the index.
  ASSERT(file_exists(faidx_file_name.c_str()), "FAI file for FASTA file does not exist: " + faidx_file_name 
         + "\nTry running the command:\nsamtools faidx " + fasta_file_name);
  
  input_tam = sam_open(tam_file_name.c_str());
  if (!input_tam) {
    cerr << "Could not open tam file: " << tam_file_name << endl;
    assert(input_tam);
  }
  
  // Read header from SAM file and discard (really should use it)
  bam_header_t* bam_header_no_care = sam_header_read(input_tam);
  if (bam_header_no_care) bam_header_destroy(bam_header_no_care);
  
  // but we keep a header
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
  (void)paired;
  alignments.clear();

	//uint32_t num_to_slurp = (paired) ? 2 : 1; //TODO:unused?
	string last_read_name = "";
	if (last_alignment.get() != NULL)
	{
		last_read_name = last_alignment->read_name();    
		alignments.push_back(last_alignment);
    last_alignment = counted_ptr<bam_alignment>(NULL);
	}

	//int num_slurped = 0;//TODO:unused?
	while (true)
	{
    bam_alignment* this_alignment_bam = new bam_alignment(); 
    
    counted_ptr<bam_alignment> this_alignment(this_alignment_bam);
    
    int32_t bytes = sam_read1(input_tam,bam_header, this_alignment_bam);
    if (bytes < 0) break;
    
    last_alignment = this_alignment;

		string read_name = this_alignment->read_name();
      
		if (last_read_name.size() == 0)
    {
			last_read_name = read_name;
    }
    else
    {
      if (read_name != last_read_name) break; 
    }
    alignments.push_back(this_alignment);
  }
  return (!alignments.empty());
}
 
void tam_file::write_alignments(
                                int32_t fastq_file_index, 
                                const alignment_list& alignments, 
                                vector<Trims>* trims, 
                                const cReferenceSequences* ref_seq_info_ptr,
                                bool shift_gaps
                                )
{
  (void)shift_gaps;
  
  uint32_t i=-1;
  for (alignment_list::const_iterator it=alignments.begin(); it != alignments.end(); it++)
	{
    i++;
		bam_alignment& a = *(it->get());

		stringstream aux_tags_ss;
    
		aux_tags_ss << "AS:i:" << a.aux_get_i("AS") << "\t" << "X1:i:" << alignments.size() << "\t" << "X2:i:" << fastq_file_index;

		if ((trims != NULL) && (trims->size() > i))
		{
			Trims trim = (*trims)[i];
			aux_tags_ss << "\t" << "XL:i:" << trim.L << "\t" << "XR:i:" << trim.R;
		}

		string aux_tags = aux_tags_ss.str();

		string quality_score_string = a.read_base_quality_char_string();
		string cigar_string;
    
    // @JEB experimental!!
    // Fix the cigar string by shifting gaps if asked for!
    if (0)
    //if (ref_seq_info_ptr && shift_gaps)
    {
      cigar_string =  shifted_cigar_string(a, *ref_seq_info_ptr);
    }
    else
    {
      cigar_string = a.cigar_string();
    }
    

		vector<string> ll;
		ll.push_back(a.read_name());
		//if (verbose) cerr << a.read_name() << endl;
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

void tam_file::write_split_alignment(uint32_t min_indel_split_len, const alignment_wrapper& a)
{
  // Debug
  //if (a.read_name() == "3:234041") {
  //  cout << "stop here << endl";
  //}
  
	uint32_t q_length = a.read_length();

  string qual_string = a.read_base_quality_char_string();
	string qseq_string = a.read_char_sequence();

	uint32_t rpos = a.reference_start_1();
	uint32_t qpos = a.query_start_1();

	vector<pair<char,uint16_t> > cigar_list = a.cigar_pair_array();
	char op;
	uint32_t len;
	uint32_t i = 0;
	while (i < cigar_list.size())
	{
		uint32_t rstart = rpos;
		uint32_t qstart = qpos;

		string cigar_string = "";
		while (i < cigar_list.size()) //CIGAR
		{
			op = cigar_list[i].first;
			len = cigar_list[i].second;

			//#print "$cigar_string\n";
			//#print Dumper($c);
      
      if (op == 'I')
			{
				if (len >= min_indel_split_len) break;
				qpos += len;
			}
			else if (op == 'D')
			{
				if (len >= min_indel_split_len) break;
				rpos += len;
			}
			else if (op == 'M')
			{
				qpos += len;
				rpos += len;
			}
			
      // @JEB we must skip over soft padding at the ends, it will be added back in one chunk
      if (op != 'S')
      {
        cigar_string += to_string(len) + op;
      }
      i++;
    }
    
		if (qpos > qstart)
		{
			//add padding to the sides of the match
			uint32_t left_padding = qstart - 1;
			uint32_t right_padding = q_length - qpos + 1;

			cigar_string = ( left_padding > 0 ? to_string(left_padding) + "S" : "" ) + cigar_string + ( right_padding > 0 ? to_string(right_padding) + "S" : "" );

			//#print "Q: $qstart to " . ($qpos-1) . "\n";

      vector<string> ll = make_vector<string>
				(a.read_name())
				(to_string(a.flag()))
				(to_string(this->bam_header->target_name[a.reference_target_id()]))
				(to_string(rstart))
				(to_string(a.quality()))
				(cigar_string)
				("*")("0")("0") //mate info
				(qseq_string)
				(qual_string)
				//#$aux_tags
			;
			output_tam << join(ll, "\t") << endl;
		}
    
    // If the inserted region matches to the next part of the match, then we want to adjust where that begins
    if ( (op == 'I') && (qpos - len >= 1) )
    {
      string previous_string = qseq_string.substr(qpos - len - 1, len);
      string insert_string = qseq_string.substr(qpos - 1, len);
      if (previous_string == insert_string)
      {
        //cout << a.read_name() << endl;
        rpos -= len;
        i++;
        ASSERT(cigar_list[i].first == 'M', "Expected match in cigar string when resolving repeat insert.");
        cigar_list[i].second += len;
      }
    }

		//move up to the next match position
		while (i < cigar_list.size()) //CIGAR
		{
			op = cigar_list[i].first;
			len = cigar_list[i].second;

			if (op == 'I')
			{
				qpos += len;
			}
			else if (op == 'D')
			{
				rpos += len;
			}
			//else if (op == 'S')
			//{
			//	qpos += len;
			//	i++;
			//}
			else if (op == 'M')
			{
				break;
			}
      i++;
		}
	}

}

// Project a read alignment from a candidate junction to the reference sequence
//  and write out the result in a TAM file.
// This is probably the most complicated function in all of breseq.
// Abandon all hope, ye who enter here.
/*
	$a, 					# CJ: SAM alignment object for the read to the candidate junction
	$fastq_file_index, 		# which fastq file this read came from
	$seq_id, 				# REFERENCE: sequence id
	$reference_pos, 		# REFERENCE: position of this junction side
	$reference_strand, 		# REFERENCE: strand of this junction side (-1 or +1)
	$reference_overlap, 	# REFERENCE: amount of overlap in the reference coords on this side
	$junction_side, 		# CJ: side of the junction (0 or 1) that we are writing
	$junction_flanking, 	# CJ: number of bases before overlap in the candidate junction sequence
	$junction_overlap, 		# CJ: amount of overlap in the candidate junction sequence that we aligned to
	$trim					# CJ: list with two items, indicating what the trim on each end is
*/
void tam_file::write_moved_alignment(
                                     const alignment_wrapper& a, 
                                     const string& rname,
                                     uint32_t fastq_file_index, 
                                     const string& seq_id, 
                                     int32_t reference_pos, 
                                     int32_t reference_strand, 
                                     int32_t reference_overlap, 
                                     uint32_t junction_side, 
                                     int32_t junction_flanking, 
                                     int32_t junction_overlap, 
                                     const Trims* trim
                                     )
{
	bool verbose = false;

	if (verbose)
	{
		cerr << "qname                 = " << a.read_name() << endl;
		cerr << "rname                 = " << rname << endl;
		cerr << "seq_id                = " << seq_id << endl;
		cerr << "reference_pos	       = " << reference_pos << endl;
		cerr << "reference_strand      = " << reference_strand << endl;
		cerr << "reference_overlap     = " << reference_overlap << endl;
		cerr << "junction_side         = " << junction_side << endl;
		cerr << "junction_flanking     = " << junction_flanking << endl;
		cerr << "junction_overlap      = " << junction_overlap << endl;
		cerr << "alignment->start      = " << a.reference_start_1() << "  alignment->end = " << a.reference_end_1() << endl;
		//cerr << "trim<<<<<<" << endl;
		//cerr << Dumper($trim);
	}

	// Which strand of the read are we on? Controls whether CIGAR is reversed
	int8_t read_strand = ((junction_side == 1) ? -1 : +1) * reference_strand;

	if (verbose)
		cerr << "read strand = " << static_cast<int32_t>(read_strand) << endl;

	uint32_t a_read_start, a_read_end;
	a.query_bounds_1(a_read_start, a_read_end);
	if (verbose)
		cerr << "a_read_start = " << a_read_start << ", a_read_end = " << a_read_end << endl;


	// Setup all of the original read information
	uint8_t* qscore = a.read_base_quality_bam_sequence(); // quality score array
	uint32_t q_length = a.read_length();

	vector<uint8_t> qual_scores;
	stringstream qseq_ss;
	for (uint32_t i = 0; i < q_length; i++)
	{
		qual_scores.push_back(*qscore);
		qscore++;
	}
	string seq = a.read_char_sequence();

	vector<pair<char,uint16_t> > cigar_list = a.cigar_pair_array();
	uint16_t flags = a.flag();

	// Remove soft padding from CIGAR (since it does not correspond to the
	// aligned positions we are dealing with. It is added back later.
	uint32_t left_padding = 0;
	uint32_t right_padding = 0;
	if (cigar_list[0].first == 'S')
	{
		left_padding = cigar_list.front().second;
		cigar_list.erase(cigar_list.begin());
	}
	if (cigar_list.back().first == 'S')
	{
		right_padding = cigar_list.back().second;
		cigar_list.pop_back();
	}

	// If we are operating on the opposite read strand,
	// Reverse complement sequence, reverse quals, and toggle strand bit
	if (read_strand == -1)
	{
		seq = reverse_complement(seq);
		reverse(qual_scores.begin(), qual_scores.end());
		flags ^= 16; //bitwise XOR to flip strand
	}

	// this isn't allowed!
	assert(reference_overlap >= 0);

	// $junction_pos gives the position in the CJS
	// where we want to split the read and only count one side
	// For side 1, we go up to this coordinate
	// Side 2 begins after this coordinate
	uint32_t junction_pos = junction_flanking;
	if (junction_side == 1)
	{
		// Offset position to include additional sequence on this side
		junction_pos += reference_overlap;
	}
	else if (junction_side == 2)
	{
		// Offset past the common part of the alignment overlap:
		//   for negative values, this is a gap
		//   for positive values, this is the common sequence
		junction_pos += abs(junction_overlap);

		// Offset backwards for any REMAINING positive overlap.
		junction_pos -= reference_overlap;
	}
	if (verbose) cerr << "junction_pos = " << junction_pos << endl;

	////
	// split the CIGAR list into two halves and keep track of their length in the read
	///

	//if (verbose) cerr << "Original CIGAR:" << endl << Dumper($cigar_list);

	// We want to determine how much of the read matches on each side
	// of this position, use the CIGAR to correct for indels:
	// At the same time, split the CIGAR

	vector<pair<char,uint16_t> > side_1_cigar_list, side_2_cigar_list;

	uint32_t test_read_pos = a_read_start;
	uint32_t test_junction_pos = a.reference_start_1();

	// it's possible that due to overlap we are already past the position we want
	if (test_junction_pos > junction_pos)
	{
		test_junction_pos = junction_pos;
	}
	else
	{
		// Remove CIGAR operations until we have enough length for side 1
		while (cigar_list.size() > 0)
		{
			pair<char,uint16_t> c = cigar_list[0];
			cigar_list.erase(cigar_list.begin());
			char op = c.first;
			uint16_t n = c.second;
			if (op == 'I') //insertion in read relative to reference
			{
				test_read_pos += n;
			}
			else
			{
				// If we shot past the correct position, backtrack
				int32_t overshoot = test_junction_pos + n - junction_pos - 1;
				if (overshoot > 0)
				{
					// reduce $n so that overshoot is removed
					n -= overshoot;
					// push back the reduced match length onto the CIGAR
					// so that it can become part of the side 2 match
					pair<char,uint16_t> new_element = make_pair<char,uint16_t>(c.first, overshoot);
					cigar_list.insert(cigar_list.begin(), new_element);
				}
				// After $n may have been adjusted add it to both positions
				test_junction_pos += n;
				if (op != 'D') //if not deletion in read relative to reference
					test_read_pos += n;
			}

			pair<char,uint16_t> new_element = make_pair<char,uint16_t>(c.first, n);
			side_1_cigar_list.push_back(new_element);
			if (test_junction_pos >= junction_pos) break;
		}
	}

  side_2_cigar_list = cigar_list;
	// Use the remaining CIGAR operations to construct side 2
  // @JEB below causes memory error
  //	copy(cigar_list.begin(), cigar_list.end(), side_2_cigar_list.begin());

	if (verbose) cout << "test_read_pos = " << test_read_pos << endl;
	if (verbose) cout << "test_junction_pos = " << test_junction_pos << endl;

	// Determine the matched length on each side of the junction
	//  In the read:
	int32_t total_read_match_length = a_read_end - a_read_start + 1;
	int32_t side_1_read_match_length = test_read_pos - a_read_start;
	int32_t side_2_read_match_length = total_read_match_length - side_1_read_match_length;
	int32_t read_match_length = (junction_side == 1) ? side_1_read_match_length : side_2_read_match_length;
	if (verbose) {
		cerr << "total_read_match_length = " << total_read_match_length << endl;
		cerr << "side_1_read_match_length = " << side_1_read_match_length << endl;
		cerr << "side_2_read_match_length = " << side_2_read_match_length << endl;
		cerr << "read_match_length = " << read_match_length << endl;
	}

	//  In the candidate junction:
	int32_t total_junction_match_length = a.reference_end_1() - a.reference_start_1() + 1;
	int32_t side_1_junction_match_length = test_junction_pos - a.reference_start_1();
	if (side_1_junction_match_length < 0) side_1_junction_match_length = 0;
	int32_t side_2_junction_match_length = total_junction_match_length - side_1_junction_match_length;
  
	int32_t junction_match_length = (junction_side == 1) ? side_1_junction_match_length : side_2_junction_match_length;
	if (verbose) {
		cerr << "total_junction_match_length = " << total_junction_match_length << endl;
		cerr << "side_1_junction_match_length = " << side_1_junction_match_length << endl;
		cerr << "side_2_junction_match_length = " << side_2_junction_match_length << endl;
		cerr << "junction_match_length = " << junction_match_length << endl;
	}

	// we could still be short of the junction, which means we will
	// have to offset the reference coordinate of this piece of the match
	// both of these compute positive numbers for how short we are
	int32_t short_of_junction;
	if (junction_side == 1)
	{
		short_of_junction =  junction_pos - (a.reference_start_1() + total_junction_match_length - 1);
		//we end short of the junction if < 0, so we have to offset the reference position by this
	}
	// or started matching after the junction
	else //if (junction_side == 2)
	{
		short_of_junction =  a.reference_start_1() - junction_pos - 1;
	}
	if (short_of_junction < 0) short_of_junction = 0;

	// get the right side of the junction
	cigar_list = (junction_side == 1) ? side_1_cigar_list : side_2_cigar_list;

	if (verbose) {
		cerr << "Short of junction = " << short_of_junction << endl;
		
		// Lots of debug output to make sure the CIGAR list is proper...
		//cerr << "CIGAR for each side:<< endl . Dumper(\@side_1_cigar_list, \@side_2_cigar_list) if ($verbose);
		//cerr << "CIGAR for this junction side:<< endl .Dumper($cigar_list) if ($verbose);

		// Add original padding to one end and padding to the other side representing
		// the piece that was not used (is aligned to the other side of the junction)
		cerr << "Left Padding = " << left_padding << ", Right Padding = " << right_padding << endl;
	}

	//additional padding on the end that is blocked
	if (junction_side == 2)
		left_padding += side_1_read_match_length;
	else //if (junction_side == 1)
		right_padding += side_2_read_match_length;

	if (verbose)
		cerr << "Adjusted Left Padding = " << left_padding << ", Adjusted Right Padding = " << right_padding << endl;

	if (left_padding > 0)
	{
		pair<char,uint16_t> new_element = make_pair<char,uint16_t>('S', left_padding);
		cigar_list.insert(cigar_list.begin(), new_element);
	}
	if (right_padding > 0)
	{
		pair<char,uint16_t> new_element = make_pair<char,uint16_t>('S', right_padding);
		cigar_list.push_back(new_element);
	}

	if (read_strand == -1)
		reverse(cigar_list.begin(), cigar_list.end());

	//if (verbose) cerr << "Final CIGAR:" << endl << Dumper($cigar_list);

	/// It's possible there may not actually be any match on this side
	///  in cases of overlap. We must bail or get negative values
	///  in the resulting CIGAR string.
	//#	return if (($junction_side == 1) && ($side_1_ref_match_length < 0));
	//#	return if (($junction_side == 2) && ($side_2_ref_match_length < 0));

	// Determine the reference coordinate we will write out for this junction side match.
	// Recall:
	//  strand == 1 means this is the lowest coordinate of that junction side sequence
	//  strand == 0 means this is the highest coordinate
	int32_t reference_match_start = (reference_strand == 1) ? reference_pos + short_of_junction : reference_pos - (junction_match_length - 1) - short_of_junction;

	////
	//// Convert the CIGAR list back to a CIGAR string
	////
	//// at the same time check to make sure the length
	//// is corect and that there are no negative nums
	stringstream cigar_string_ss;
	uint32_t cigar_length = 0;
	for (uint32_t i = 0; i < cigar_list.size(); i++) //CIGAR
	{
		char op = cigar_list[i].first;
		uint32_t len = static_cast<uint32_t>(cigar_list[i].second);

		assert(len > 0);
		cigar_string_ss << len << op;
		if (op != 'D') cigar_length += len;
	}
	string cigar_string = cigar_string_ss.str();
  
	////
	//// Assemble the quality score string
	////
	stringstream quality_score_ss;
	for (uint32_t i = 0; i < qual_scores.size(); i++)
    {
		quality_score_ss << static_cast<char>(qual_scores[i] + 33);
	}
	string quality_score_string = quality_score_ss.str();

	////
	//// Setup custom aux tags
	////

	stringstream aux_tags_ss;
	aux_tags_ss << "AS:i:" << a.aux_get_i("AS") << "\t" << "X1:i:" << 1 << "\t" << "X2:i:" << fastq_file_index;

	//this flag indicates this is a junction match and which side of the match is in the middle of the read across the junction
	int32_t within_side = (reference_strand == 1) ? junction_side : (junction_side + 1) % 2;
	aux_tags_ss << "\t" << "XJ:i:" << within_side;

	//handle putting the trims in the right places
	//need to be aware if read is trimmed out of existence??
	if (trim != NULL)
	{
		string trim_left = (junction_side == 1) ? to_string(trim->L+left_padding) : "0";
		string trim_right = (junction_side == 1) ? "0" : to_string(trim->R+right_padding);
		if (read_strand == -1) swap(trim_left, trim_right);
		aux_tags_ss << "\t" << "XL:i:" << trim_left << "\t" << "XR:i:" << trim_right;
	}

	string aux_tags = aux_tags_ss.str();

	////
	//// Create the TAM line and print
	////

  vector<string> ll = make_vector<string>
		(a.read_name() + "-M" + to_string(junction_side))
		(to_string(fix_flags(a.flag())))
		(seq_id)
		(to_string(reference_match_start))
		(to_string(a.quality()))
		(cigar_string)
		(a.proper_pair() ? "=" : "*")
		(to_string(a.mate_start_1()))
		(to_string(a.isize()))
		(seq)
		(quality_score_string)
		(aux_tags)
	;
	string l = join(ll, "\t") + "\n";
	if (verbose) cout << l;

	assert(cigar_length == q_length);
	output_tam << l;
}

}
