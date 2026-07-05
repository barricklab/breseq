/*****************************************************************************

 AUTHORS

   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors

 LICENSE AND COPYRIGHT

   Copyright (c) 2008-2010 Michigan State University
   Copyright (c) 2011-2025 The University of Texas at Austin
   Copyright (c) 2025-     Michigan State University

   breseq is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2, or (at your option) any later version.

   SPDX-License-Identifier: GPL-2.0-or-later

*****************************************************************************/

#include "alignment.h"

#include "fastq.h"
#include "reference_sequence.h"

using namespace std;

namespace breseq {

const char alignment_wrapper::op_to_char[10] = "MIDNSHP=X";

string alignment_wrapper::read_char_stranded_sequence_1(uint32_t start_1, uint32_t end_1) const { 
  ASSERT(start_1 <= end_1, "Start (" + to_string(start_1) + ") must be less than end (" + to_string(end_1) + ").");
  
  string s = read_char_sequence();
  if (reversed()) s = reverse_complement(s);
  s = s.substr(start_1-1, end_1-start_1+1);
  return s;
}

  
/*! Does this alignment have any redundancies?
 */
bool alignment_wrapper::is_redundant() const {
	return (redundancy() > 1);
}


/*! Number of redundancies at this alignment.
    Defaults to 1 when custome breseq tag is missing
 */
uint32_t alignment_wrapper::redundancy() const {
      
  uint8_t *aux = bam_aux_get(_a, "X1");
  if (aux == NULL) return 1;
  return bam_aux2i(aux);
}



/*! Calculate the length of this query out to the last non-clip, non-skip.

 Replaced with cigar_query_length()
uint32_t alignment::query_length() const {
	uint32_t* cigar = bam_get_cigar(_a); // cigar array for this alignment
	uint32_t qlen = bam_cigar2qlen(_a->core.n_cigar, cigar); // total length of the query
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
  uint32_t* cigar = bam_get_cigar(_a); // cigar array for this alignment
	uint32_t start=1, end=bam_cigar2qlen(_a->core.n_cigar, cigar);
	
  /*
  if (this->read_name()=="1:369") {
    cout << "debug" << end;
  }
  */
  
	// start:
  uint32_t i;
  uint32_t op;
  uint32_t len;
  for(i=0; i<=_a->core.n_cigar; i++) {
    op = cigar[i] & BAM_CIGAR_MASK;
    len = cigar[i] >> BAM_CIGAR_SHIFT;
    
    // Skip any operations that don't involve aligning bases to the reference
    if( (op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP) ) {
			break;
    }
    
    // Only move forward our position in read for soft clipping
    if(op == BAM_CSOFT_CLIP) {
      start += len;
    }
  }
  
  //move past low quality bases
  if (min_qual) {
    
    for(; i<=_a->core.n_cigar; i++) {
      
      // Skip operations that don't involve read bases being aligned
      if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP) && (op != BAM_CINS)) {
        for(uint32_t j=0; j<len; j++) {
          if (read_base_quality_1(start) > min_qual) {
            goto finish_start;
          }
          start++;
        }
      }
      
      // Only move forward our position in read for soft clipping
      if (op == BAM_CSOFT_CLIP) {
        start += len;
      }
      
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
    // skip any operations that don't involve aligning bases to the reference
    if( (op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP) ) {
      break;
    }
    
    // Only move forward our position in read for soft clipping
    if(op == BAM_CSOFT_CLIP) {
      end -= len;
    }
  }
  
  //move past low quality bases
  if (min_qual) {

    for(; i>0; --i) {
      
      // Skip operations that don't involve read bases being aligned
      if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP) && (op != BAM_CINS)) {
        for(uint32_t j=0; j<len; j++) {
          if (read_base_quality_1(end) > min_qual) {
            goto finish_end;
          }
          end--;
        }
      }
      // Only move backward our position in read for soft clipping
      if (op == BAM_CSOFT_CLIP) {
        end -= len;
      }
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
  uint32_t* cigar = bam_get_cigar(_a); // cigar array for this alignment
  int32_t pos = 1;
  
  for(uint32_t j=0; j<=_a->core.n_cigar; j++) {
    uint32_t op = cigar[j] & BAM_CIGAR_MASK;
    uint32_t len = cigar[j] >> BAM_CIGAR_SHIFT;
		
    // if we encounter padding, or a gap in reference then we are done
    if(op != BAM_CSOFT_CLIP) {
			break;
    }
    
    pos += len;
  }
  
  return pos;
}


uint32_t alignment_wrapper::query_end_1() const {
  // traverse the cigar array
  uint32_t* cigar = bam_get_cigar(_a); // cigar array for this alignment
  int32_t pos = bam_cigar2qlen(_a->core.n_cigar, cigar); // total length of the query includes soft clipped nucleotides
  
  uint32_t j;
  for(j=(_a->core.n_cigar-1); j>0; j--) {
    uint32_t op = cigar[j] & BAM_CIGAR_MASK;
    uint32_t len = cigar[j] >> BAM_CIGAR_SHIFT;
    
    // if we encounter non-padding padding, or a gap in reference then we are done
    if(op != BAM_CSOFT_CLIP) {
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
    
    uint32_t* cigar = bam_get_cigar(_a); // cigar array for this alignment
    
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
          if (read_base_quality_1(start) > min_qual) {
            goto finish_start;
          }
          start++;
        }
      }
      // only move reference position for ref-skip or deletion in read releative to reference 
      if ((op == BAM_CREF_SKIP) || (op == BAM_CDEL)) {
        start += len;
      }
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
  uint32_t end = bam_endpos(_a);
  
  if (min_qual) {
    uint32_t* cigar = bam_get_cigar(_a); // cigar array for this alignment

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
    }
  
  //move past low quality bases
    
    for(; i>0; --i) {
      
      if((op != BAM_CSOFT_CLIP) && (op != BAM_CHARD_CLIP) && (op != BAM_CREF_SKIP) && (op != BAM_CDEL)) {
        for(uint32_t j=0; j<len; j++) {
          if (read_base_quality_1(end) > min_qual) {
            goto finish_end;
          }
          end--;
        }
      }
      
      if ((op == BAM_CREF_SKIP) || (op == BAM_CDEL)) {
        end -=len;
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
  
  uint32_t* cigar = bam_get_cigar(_a); // cigar array for this alignment
	
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
  
  // Get the cigar array (reverse if direction is reversed)
  vector<pair<char,uint16_t> > cigar_pair_array = this->cigar_pair_char_op_array();
  if (!dir) reverse(cigar_pair_array.begin(), cigar_pair_array.end());

  char op_0 = cigar_pair_array.front().first;
  char op_last = cigar_pair_array.back().first;
  
  //remove soft padding
  if (op_0 == 'S') {
    cigar_pair_array.erase(cigar_pair_array.begin());
  }
  if (op_last == 'S') {
    cigar_pair_array.erase(cigar_pair_array.end()-1);
  }
  
  qry_mismatch_pos = -1;  // 1-indexed
  ref_mismatch_pos = -1;  // 1-indexed
  
  bool last_encountered_indel = false;
  uint32_t positive_overlap = (overlap > 0) ? static_cast<uint32_t>(overlap) : 0;
  uint32_t r_pos = 0; // 1-indexed
  uint32_t q_pos = 0; // 1-indexed
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
    if (op_0 == 'I')
    {
      last_encountered_indel = true;
      q_pos++;
      ref_mismatch_pos = r_pos;
      qry_mismatch_pos = q_pos;
    }
    else if (op_0 == 'D')
    {
      last_encountered_indel = true;
      r_pos++;
      ref_mismatch_pos = r_pos;
      qry_mismatch_pos = q_pos;
    }
    else
    {
      r_pos++;
      q_pos++;
      if (q_str[q_pos-1] != r_str[r_pos-1])
      {
        last_encountered_indel = false;
        ref_mismatch_pos = r_pos;
        qry_mismatch_pos = q_pos;
      }

    }
    
    if (verbose)
      cout << r_pos << " " << r_str[r_pos-1] << " " << q_pos << " " << q_str[q_pos-1] << endl;
  }
  
  if (verbose)
  {
    if (qry_mismatch_pos >= 0)
      cout << "Query Mismatch At = " << qry_mismatch_pos << " Reference Mismatch At = " << ref_mismatch_pos << endl;
  }
  
  // Go back and add exact matches, which may extend after correcting for an indel!
  // @JEB - 2013-10-12 - this code is not perfect: 
  //  (1) It will not let the new match give any more than the current overlap
  //  (2) It will not correct for indels that are outside of the current overlap
  if (last_encountered_indel) {
    if (verbose)
      cout << "Indel re-correction of overlap" << endl;
    
    r_pos = ref_mismatch_pos;
    q_pos = qry_mismatch_pos;  
      
    while ((r_pos != 0) && (q_pos != 0))
    {
      if (verbose)
        cout << r_pos << " " << r_str[r_pos-1] << " " << q_pos << " " << q_str[q_pos-1] << endl;

      if ( r_str[r_pos-1] != q_str[q_pos-1] ) break;
       
      qry_mismatch_pos--;
      ref_mismatch_pos--;   
      
      r_pos--;
      q_pos--;
    }
  }
  
  if (verbose)
  {
    if (qry_mismatch_pos >= 0)
      cout << "Query Mismatch At = " << qry_mismatch_pos << " Reference Mismatch At = " << ref_mismatch_pos << " | Overlap = " << overlap << endl;
    cout << "<====" << endl;
  }
}
  
/*! Build a BAM header from a FASTA index file.
 *  Replaces the deprecated sam_header_read2() function removed in modern htslib.
 *  Accepts either the FASTA file itself or the .fai index path.
 */
bam_hdr_t* make_bam_header_from_faidx(const string& fasta_file_name) {
  // Strip .fai suffix to get base FASTA path for fai_load
  string fasta = fasta_file_name;
  if (fasta.size() > 4 && fasta.substr(fasta.size() - 4) == ".fai")
    fasta = fasta.substr(0, fasta.size() - 4);

  faidx_t *fai(NULL);
  if (file_exists(fasta.c_str())) {
    fai = fai_load(fasta.c_str());
  }

  // There is an instance where we call this with a non-existent FASTA file and
  // it should not throw an error in resolve_alignments
  string hdr_text = "@HD\tVN:1.6\tSO:unsorted\n";
  if (fai) {
    int n = faidx_nseq(fai);
    for (int i = 0; i < n; i++) {
      const char *name = faidx_iseq(fai, i);
      int len = faidx_seq_len(fai, name);  // faidx_seq_len returns int (portable across htslib versions)
      hdr_text += string("@SQ\tSN:") + name + "\tLN:" + to_string(len) + "\n";
    }
    fai_destroy(fai);
  }
  
  bam_hdr_t *hdr = sam_hdr_parse(hdr_text.size(), hdr_text.c_str());
  ASSERT(hdr, "Could not parse BAM header from FASTA index: " + fasta);
  return hdr;
}

bam_file::bam_file(const string& bam_file_name, const string& fasta_file_name, ios_base::openmode mode)
: bam_header(NULL), m_bam_file(NULL)
{
  if (mode == ios_base::in)
  {
    open_read(bam_file_name, fasta_file_name);
  }
  else
  {
    open_write(bam_file_name, fasta_file_name);
  }
}

bam_file::~bam_file()
{
  if (m_bam_file) { hts_close(m_bam_file); m_bam_file=NULL; }
  if (bam_header) sam_hdr_destroy(bam_header);
}

void bam_file::open_read(const string& bam_file_name, const string& fasta_file_name)
{
  string faidx_file_name(fasta_file_name);
  faidx_file_name += ".fai";

  // Alternately, we could automatically generate the index.
  ASSERT(file_exists(faidx_file_name.c_str()), "FAI file for FASTA file does not exist: " + faidx_file_name
         + "\nTry running the command:\nsamtools faidx " + fasta_file_name);

  m_bam_file = hts_open(bam_file_name.c_str(), "rb");
  ASSERT(m_bam_file, "Could not open bam file: " + bam_file_name);

  bam_header = sam_hdr_read(m_bam_file);
}

void bam_file::open_write(const string& bam_file_name, const string& fasta_file_name)
{
  m_bam_file = hts_open(bam_file_name.c_str(), "wb");
  ASSERT(m_bam_file, "Could not open bam file: " + bam_file_name);
  bam_header = make_bam_header_from_faidx(fasta_file_name);
  (void) sam_hdr_write(m_bam_file, bam_header);
}

bool bam_file::read_alignments(alignment_list& alignments, bool paired)
{
  (void)paired;
  alignments.clear();
  
  string last_read_name = "";
  if (m_last_alignment.get() != NULL)
  {
    last_read_name = m_last_alignment->read_name();
    alignments.push_back(m_last_alignment);
    m_last_alignment = counted_ptr<bam_alignment>(NULL);
  }
  
  while (true)
  {
    bam_alignment* this_alignment_bam = new bam_alignment();
    
    counted_ptr<bam_alignment> this_alignment(this_alignment_bam);
    
    int32_t bytes = sam_read1(m_bam_file, bam_header, this_alignment_bam);
    if (bytes < 0) break;
    
    m_last_alignment = this_alignment;
    
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

bool soft_clip_alignment_ends(
                              vector<pair<char,uint16_t> >& cigar_list,
                              uint32_t& reference_start_1,
                              const string& read_seq_top_strand,
                              const cReferenceSequences& ref_seq_info,
                              const string& seq_id,
                              bool protect_left,
                              bool protect_right,
                              uint32_t left_trim_reads,
                              uint32_t right_trim_reads
                              )
{
  // Minimum fraction of a chewed-back terminal run that must be mismatches for it to be
  // soft-clipped. Tunable; higher = more conservative (clips only denser mismatch runs).
  // At 0.33 the walk always examines (at least) the first 3 non-trimmed bases at each end.
  const double kMinEndMismatchFraction = 0.33;

  if (cigar_list.empty()) return false;

  // An end already carrying a soft-clip is left alone (this also protects the junction/middle
  // side of -M1/-M2 split reads, which always carries padding 'S').
  if (cigar_list.front().first == 'S') protect_left = true;
  if (cigar_list.back().first == 'S') protect_right = true;
  if (protect_left && protect_right) return false;

  // Reference span (M/D consume reference) needed to fetch the aligned reference region.
  uint32_t ref_span = 0;
  for (size_t i = 0; i < cigar_list.size(); i++) {
    char op = cigar_list[i].first;
    if (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X') ref_span += cigar_list[i].second;
  }
  if (ref_span == 0) return false;
  const string ref_seq = ref_seq_info[seq_id].get_sequence_1(reference_start_1, reference_start_1 + ref_span - 1);

  // Expand the (non-soft-clip) core of the CIGAR into per-column records. Each column advances
  // the read and/or the reference and is flagged as a mismatch (a substituted M base, or any
  // inserted/deleted base -- each indel base counts as one mismatch column).
  struct Col { char op; bool mismatch; bool read_adv; bool ref_adv; uint32_t read_pos; };
  vector<Col> cols;
  uint32_t read_i = 0; // 0-based into read_seq_top_strand
  uint32_t ref_i = 0;  // 0-based into ref_seq
  for (size_t c = 0; c < cigar_list.size(); c++) {
    char op = cigar_list[c].first;
    uint16_t n = cigar_list[c].second;
    if (op == 'S' || op == 'H') { if (op == 'S') read_i += n; continue; }
    for (uint16_t k = 0; k < n; k++) {
      Col col; col.op = op; col.read_pos = read_i; // read base this column consumes (D: boundary)
      if (op == 'I') {
        col.mismatch = true;  col.read_adv = true;  col.ref_adv = false; read_i++;
      } else if (op == 'D' || op == 'N') {
        col.mismatch = true;  col.read_adv = false; col.ref_adv = true;  ref_i++;
      } else { // M / = / X
        bool differs = (op == 'X') ? true
                     : (op == '=') ? false
                     : (read_seq_top_strand[read_i] != ref_seq[ref_i]);
        col.mismatch = differs; col.read_adv = true; col.ref_adv = true; read_i++; ref_i++;
      }
      cols.push_back(col);
    }
  }
  if (cols.empty()) return false;

  // Walk inward from one end, returning how many columns to chew back into a soft-clip.
  // A base within the pre-clip trimmed tip that merely *matches* the reference is "unknown" -- in a
  // low-complexity / indel-adjacent region a match does not prove correct alignment -- so it is
  // skipped (counts toward neither the examined total n nor the mismatch total m), though a clip may
  // still extend physically through it. A trimmed base that *mismatches* is a real mismatch and
  // counts normally. Examine the next counted column only if a mismatch there could keep the
  // mismatch fraction >= threshold ((m+1)/(n+1) >= T); record the deepest mismatch whose fraction
  // (m/n) is >= T -- returned as a physical column depth.
  struct Walk {
    static uint32_t clip_columns(const vector<Col>& cols, bool from_left, double T,
                                 uint32_t trim_reads, uint32_t read_length) {
      uint32_t n = 0, m = 0, clip = 0;
      for (uint32_t p = 1; p <= cols.size(); p++) {
        size_t idx = from_left ? (p - 1) : (cols.size() - p);
        bool trimmed = from_left ? (cols[idx].read_pos < trim_reads)
                                 : (cols[idx].read_pos + trim_reads >= read_length);
        if (trimmed && !cols[idx].mismatch) continue; // trimmed match = unknown; clip may pass through
        if (static_cast<double>(m + 1) < T * (n + 1)) break; // even a mismatch here can't reach T
        n++;
        if (cols[idx].mismatch) {
          m++;
          if (static_cast<double>(m) >= T * n) clip = p;
        }
      }
      return clip;
    }
  };

  const uint32_t read_length = static_cast<uint32_t>(read_seq_top_strand.size());
  uint32_t left_clip = protect_left  ? 0 : Walk::clip_columns(cols, true,  kMinEndMismatchFraction, left_trim_reads, read_length);
  uint32_t right_clip = protect_right ? 0 : Walk::clip_columns(cols, false, kMinEndMismatchFraction, right_trim_reads, read_length);

  // A soft-clip must not be left adjacent to a deletion (leading/trailing 'D' is invalid), so
  // absorb any deletion columns exposed at the new boundary into the clip.
  while (left_clip > 0 && left_clip < cols.size() && cols[left_clip].op == 'D') left_clip++;
  while (right_clip > 0 && right_clip < cols.size() && cols[cols.size() - 1 - right_clip].op == 'D') right_clip++;

  if (left_clip == 0 && right_clip == 0) return false;

  // Degenerate case: the two clips would consume the whole alignment -- leave it unchanged.
  if (left_clip + right_clip >= cols.size()) return false;

  // Tally read/reference bases removed at each end.
  uint32_t lead_S = (cigar_list.front().first == 'S') ? cigar_list.front().second : 0;
  uint32_t trail_S = (cigar_list.back().first == 'S') ? cigar_list.back().second : 0;
  uint32_t left_read = 0, left_ref = 0, right_read = 0, right_ref = 0;
  for (uint32_t j = 0; j < left_clip; j++) { if (cols[j].read_adv) left_read++; if (cols[j].ref_adv) left_ref++; }
  for (uint32_t j = 0; j < right_clip; j++) { size_t idx = cols.size() - 1 - j; if (cols[idx].read_adv) right_read++; if (cols[idx].ref_adv) right_ref++; }

  // Rebuild the CIGAR: soft-clip + surviving middle columns (run-length encoded) + soft-clip.
  vector<pair<char,uint16_t> > out;
  uint32_t new_lead_S = lead_S + left_read;
  uint32_t new_trail_S = trail_S + right_read;
  if (new_lead_S > 0) out.push_back(make_pair('S', static_cast<uint16_t>(new_lead_S)));
  for (uint32_t j = left_clip; j < cols.size() - right_clip; j++) {
    char op = cols[j].op;
    if (!out.empty() && out.back().first == op && op != 'S') out.back().second += 1;
    else out.push_back(make_pair(op, static_cast<uint16_t>(1)));
  }
  if (new_trail_S > 0) out.push_back(make_pair('S', static_cast<uint16_t>(new_trail_S)));

  cigar_list = out;
  reference_start_1 += left_ref; // clipping the left end shifts the alignment start rightward
  return true;
}

void bam_file::write_alignments(
                                int32_t fastq_file_index,
                                const alignment_list& alignments,
                                const SequenceTrimsList* trims_list,
                                const cReferenceSequences* ref_seq_info_ptr,
                                bool shift_gaps,
                                bool soft_clip_ends
                                )
{

  uint32_t i=-1;
  for (alignment_list::const_iterator it=alignments.begin(); it != alignments.end(); it++) {

    i++;
    bam_alignment& a = *(it->get());

    // --- CIGAR + POS (done first so a soft-clipped alignment's new coordinates are available
    //     when we (re)compute the XL/XR trims below) ---
    string cigar_string;
    uint32_t reference_start_1 = a.reference_start_1();
    bool soft_clipped = false;
    vector<pair<char,uint16_t> > cigar_pair; // populated only on the soft-clip path

    if (soft_clip_ends && ref_seq_info_ptr) {
      // Gap-shift (matching shifted_cigar_string) then chew back mis-mapped ends into soft-clips.
      cigar_pair = a.cigar_pair_char_op_array();
      string read_seq = a.read_char_sequence();
      if (shift_gaps) {
        const string ref_region = (*ref_seq_info_ptr)[a.reference_target_id()].get_sequence_1(a.reference_start_1(), a.reference_end_1());
        shift_indels_in_cigar_array(cigar_pair, ref_region, read_seq);
      }
      // Pre-clip trims: their trimmed tip bases are "unknown" to the end-clipping heuristic.
      Trims pre_trim; pre_trim.L = 0; pre_trim.R = 0;
      if (trims_list != NULL) pre_trim = get_alignment_trims(a, *trims_list);
      soft_clipped = soft_clip_alignment_ends(cigar_pair, reference_start_1, read_seq, *ref_seq_info_ptr,
                                              bam_header->target_name[a.reference_target_id()], false, false,
                                              pre_trim.L, pre_trim.R);
      cigar_string = alignment_wrapper::cigar_op_array_to_cigar_string(cigar_pair);
    } else if (ref_seq_info_ptr && shift_gaps) {
      cigar_string = shifted_cigar_string(a, *ref_seq_info_ptr);
    } else {
      cigar_string = a.cigar_string();
    }

    stringstream aux_tags_ss;

    uint32_t as;
    bool AS_found = a.aux_get_i("AS", as);
    ASSERT(AS_found, "Could not find required tag AS for alignment.");

    aux_tags_ss << "AS:i:" << as << "\t" << "X1:i:" << alignments.size() << "\t" << "X2:i:" << fastq_file_index;

    if (trims_list != NULL) {
      Trims trim;
      if (soft_clipped) {
        // A soft-clipped read ends at a different reference position and has a larger soft-clip
        // offset, so compute its trims from the new coordinates (mirrors get_alignment_trims():
        // genomic trim at the mapped end + the bases lying outside the match, i.e. the soft-clip).
        uint32_t tid = a.reference_target_id();
        uint32_t ref_span = 0;
        for (size_t c = 0; c < cigar_pair.size(); c++) {
          char op = cigar_pair[c].first;
          if (op=='M'||op=='D'||op=='N'||op=='='||op=='X') ref_span += cigar_pair[c].second;
        }
        uint32_t new_ref_start_0 = reference_start_1 - 1;
        uint32_t new_ref_end_0 = new_ref_start_0 + ref_span - 1;
        uint32_t lead_S = (cigar_pair.front().first == 'S') ? cigar_pair.front().second : 0;
        uint32_t trail_S = (cigar_pair.back().first == 'S') ? cigar_pair.back().second : 0;
        trim.L = (*trims_list)[tid].left_trim_0(new_ref_start_0) + lead_S;
        trim.R = (*trims_list)[tid].right_trim_0(new_ref_end_0) + trail_S;
      } else {
        trim = get_alignment_trims(a, *trims_list);
      }
      aux_tags_ss << "\t" << "XL:i:" << trim.L << "\t" << "XR:i:" << trim.R;
    }

    string xp;
    if (a.aux_get_Z("XP", xp)) {
      aux_tags_ss << "\t" << "XP:Z:" << xp;
    }

    string aux_tags = aux_tags_ss.str();

    string quality_score_string = a.read_base_quality_char_string();
    if (quality_score_string[0] == ' ') {
      quality_score_string = alignments.read_base_quality_char_string;
      if (alignments.read_base_quality_char_string_reversed ^ a.reversed())
        quality_score_string = reverse_string(quality_score_string);
    }
    ASSERT(quality_score_string.size() > 0, "Attempt to write read with no quality scores: " + a.read_name());

    vector<string> ll;
    ll.push_back(a.read_name());
    ll.push_back(to_string(fix_flags(a.flag())));
    ll.push_back(bam_header->target_name[a.reference_target_id()]);
    ll.push_back(to_string(reference_start_1));
    ll.push_back(to_string<uint32_t>(a.mapping_quality()));
    ll.push_back(cigar_string);

    if ((a.flag() & BAM_FPAIRED) == 0) {
      ll.push_back("*");
      ll.push_back("0");
      ll.push_back("0");
    } else if (a.mate_reference_target_id() == a.reference_target_id()) {
      ll.push_back("=");
      ll.push_back(to_string<int32_t>(a.mate_start_1()));
      ll.push_back(to_string<int32_t>(a.insert_size()));
    } else {
      ll.push_back(bam_header->target_name[a.mate_reference_target_id()]);
      ll.push_back(to_string<int32_t>(a.mate_start_1()));
      ll.push_back("0"); // TLEN is undefined when mates are on different reference sequences
    }

    ll.push_back(a.read_char_sequence());
    ll.push_back(quality_score_string);
    ll.push_back(aux_tags);

    {
      string sam_line = join(ll, "\t");
      kstring_t ks = KS_INITIALIZE;
      kputs(sam_line.c_str(), &ks);
      bam1_t* b = bam_init1();
      ASSERT(sam_parse1(&ks, bam_header, b) >= 0, "sam_parse1 failed: " + sam_line);
      ASSERT(sam_write1(m_bam_file, bam_header, b) >= 0, "sam_write1 failed");
      bam_destroy1(b);
      ks_free(&ks);
    }
  }
}

void bam_file::write_moved_alignment(
                                     const alignment_wrapper& a,
                                     const string& junction_reference_name,
                                     uint32_t fastq_file_index,
                                     const string& seq_id,
                                     int32_t reference_pos,
                                     int32_t reference_strand,
                                     int32_t reference_overlap,
                                     uint32_t junction_side,
                                     int32_t junction_flanking,
                                     int32_t junction_overlap,
                                     const alignment_list& alignments,
                                     const SequenceTrimsList* trims_list,
                                     const cReferenceSequences* ref_seq_info_ptr,
                                     bool shift_gaps
                                     )
{
	bool verbose = false;

	// Which strand of the read are we on? Controls whether CIGAR is reversed
	int8_t read_strand = ((junction_side == 1) ? -1 : +1) * reference_strand;

	uint32_t a_read_start, a_read_end;
	a.query_bounds_1(a_read_start, a_read_end);

	// Setup all of the original read information
	uint32_t q_length = a.read_length();

  vector<char> qual_scores;
  {
    string quality_score_string = a.read_base_quality_char_string();
    if (quality_score_string[0] == ' ') {
      quality_score_string = alignments.read_base_quality_char_string;
      if (alignments.read_base_quality_char_string_reversed ^ a.reversed())
        quality_score_string = reverse_string(quality_score_string);
    }
    ASSERT(quality_score_string.size() > 0, "Attempt to write read with no quality scores: " + a.read_name());

    for (uint32_t i = 0; i < q_length; i++)
    {
      qual_scores.push_back(quality_score_string[i]);
    }
  }

	string seq = a.read_char_sequence();

	vector<pair<char,uint16_t> > cigar_list = a.cigar_pair_char_op_array();

	// Remove soft padding from CIGAR (since it does not correspond to the
	// aligned positions we are dealing with. It is added back later.
	uint16_t left_padding = 0;
	uint16_t right_padding = 0;
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

	////
	// split the CIGAR list into two halves and keep track of their length in the read
	///

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
			uint8_t op = c.first;
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
					pair<uint8_t,uint16_t> new_element(op, overshoot);
					cigar_list.insert(cigar_list.begin(), new_element);
				}
				// After $n may have been adjusted add it to both positions
				test_junction_pos += n;
				if (op != 'D') //if not deletion in read relative to reference
					test_read_pos += n;
			}

			pair<uint8_t,uint16_t> new_element(op, n);
			side_1_cigar_list.push_back(new_element);
			if (test_junction_pos > junction_pos) break;
		}
	}

	// Use the remaining CIGAR operations to construct side 2
  side_2_cigar_list = cigar_list;

	// Determine the matched length on each side of the junction
	//  In the read:
	int32_t total_read_match_length = a_read_end - a_read_start + 1;
	int32_t side_1_read_match_length = test_read_pos - a_read_start;
	int32_t side_2_read_match_length = total_read_match_length - side_1_read_match_length;
	int32_t read_match_length = (junction_side == 1) ? side_1_read_match_length : side_2_read_match_length;

	//  In the candidate junction:
	int32_t total_junction_match_length = a.reference_end_1() - a.reference_start_1() + 1;
	int32_t side_1_junction_match_length = test_junction_pos - a.reference_start_1();
	if (side_1_junction_match_length < 0) side_1_junction_match_length = 0;
	int32_t side_2_junction_match_length = total_junction_match_length - side_1_junction_match_length;

	int32_t junction_match_length = (junction_side == 1) ? side_1_junction_match_length : side_2_junction_match_length;

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
	else
	{
		short_of_junction =  a.reference_start_1() - junction_pos - 1;
	}
	if (short_of_junction < 0) short_of_junction = 0;

	// get the right side of the junction
	cigar_list = (junction_side == 1) ? side_1_cigar_list : side_2_cigar_list;

	// Shift homopolymer-length indels in this half's CIGAR as far right as the real genome
	// (seq_id, NOT the candidate junction sequence a is aligned to) allows -- same normalization
	// shifted_cigar_string() does for the resolved reference BAM, applied here before padding/
	// strand-reversal while cigar_list is still bam-native (real M/I/D content only, matching
	// the reference-forward slice below).
	if (ref_seq_info_ptr && shift_gaps && (junction_match_length > 0))
	{
		int32_t reference_match_start_for_shift = (reference_strand == 1) ? reference_pos + short_of_junction : reference_pos - (junction_match_length - 1) - short_of_junction;
		uint32_t read_start_for_shift = (junction_side == 1) ? a_read_start : (a_read_start + side_1_read_match_length);
		const string junction_side_ref_seq = (*ref_seq_info_ptr)[seq_id].get_sequence_1(reference_match_start_for_shift, reference_match_start_for_shift + junction_match_length - 1);
		const string junction_side_read_seq = a.read_char_sequence().substr(read_start_for_shift - 1, read_match_length);
		shift_indels_in_cigar_array(cigar_list, junction_side_ref_seq, junction_side_read_seq);
	}

	//additional padding on the end that is blocked
	if (junction_side == 2)
		left_padding += side_1_read_match_length;
	else
		right_padding += side_2_read_match_length;

  // It's possible that we will be putting soft padding up to an inserted or deleted base
  // in the cigar string on the end that we are padding from splitting the read - correct here
  // Example without correction: 438S1I43M1S (side 2 of match)
  if (cigar_list.size() >= 1) {
    if (junction_side == 2)
    {
      pair<char,uint16_t> c = cigar_list.front();
      char op = c.first;
      uint16_t n = c.second;
      if (op == 'I') {
        left_padding += n;
        cigar_list.erase(cigar_list.begin());
      }
    }
    else
    {
      pair<char,uint16_t> c = cigar_list.back();
      char op = c.first;
      uint16_t n = c.second;
      if (op == 'I') {
        right_padding += n;
        cigar_list.pop_back();
      }
    }
  }

  // Add the padding to the cigar list
	if (left_padding > 0)
	{
		pair<char,uint16_t> new_element = make_pair('S', left_padding);
		cigar_list.insert(cigar_list.begin(), new_element);
  }
	if (right_padding > 0)
	{
		pair<char,uint16_t> new_element = make_pair('S', right_padding);
		cigar_list.push_back(new_element);
	}

	if (read_strand == -1)
		reverse(cigar_list.begin(), cigar_list.end());

	// Determine the reference coordinate we will write out for this junction side match.
	// Recall:
	//  strand == 1 means this is the lowest coordinate of that junction side sequence
	//  strand == 0 means this is the highest coordinate
	int32_t reference_match_start = (reference_strand == 1) ? reference_pos + short_of_junction : reference_pos - (junction_match_length - 1) - short_of_junction;

	// Chew back mis-mapped bases at the OUTER end of this split read into soft-clipping, the
	// same as regular reference alignments. The junction (split) side is the middle of the
	// original read across the junction and must never be clipped; protect it explicitly (it
	// also always carries padding 'S', which soft_clip_alignment_ends() leaves alone anyway).
	if (ref_seq_info_ptr) {
		bool junction_on_left = ((junction_side == 2) == (read_strand == 1));
		// Pre-clip trim on the outer (non-junction) end -- its trimmed tip bases are "unknown" to
		// the end-clipping heuristic. The junction side is protected, so its trim value is unused.
		uint32_t left_trim = 0, right_trim = 0;
		if ((trims_list != NULL) && (reference_match_start >= 1)) {
			uint32_t tid = ref_seq_info_ptr->seq_id_to_index(seq_id);
			uint32_t seq_len = ref_seq_info_ptr->get_sequence_length(seq_id);
			uint32_t pre_ref_span = 0;
			for (size_t c = 0; c < cigar_list.size(); c++) {
				char op = cigar_list[c].first;
				if (op=='M'||op=='D'||op=='N'||op=='='||op=='X') pre_ref_span += cigar_list[c].second;
			}
			uint32_t pre_start_0 = static_cast<uint32_t>(reference_match_start) - 1;
			uint32_t pre_end_0 = pre_start_0 + pre_ref_span - 1;
			// This runs before the all-soft-padded early return, so a degenerate side can have
			// out-of-range (e.g. circular-wrapped) coordinates; fall back to 0 (unused) trim then.
			if (!junction_on_left && pre_start_0 < seq_len) left_trim  = (*trims_list)[tid].left_trim_0(pre_start_0);
			if ( junction_on_left && pre_end_0  < seq_len) right_trim = (*trims_list)[tid].right_trim_0(pre_end_0);
		}
		uint32_t rms = static_cast<uint32_t>(reference_match_start);
		soft_clip_alignment_ends(cigar_list, rms, seq, *ref_seq_info_ptr, seq_id,
		                         /*protect_left=*/junction_on_left, /*protect_right=*/!junction_on_left,
		                         left_trim, right_trim);
		reference_match_start = static_cast<int32_t>(rms);
	}

	////
	//// Convert the CIGAR list back to a CIGAR string
	////
	//// at the same time check to make sure the length
	//// is correct and that there are no negative nums
	stringstream cigar_string_ss;
	uint32_t cigar_length = 0;
  bool all_soft_padded = true;
	for (uint32_t i = 0; i < cigar_list.size(); i++) //CIGAR
	{
		char op = cigar_list[i].first;
		uint32_t len = static_cast<uint32_t>(cigar_list[i].second);

		assert(len > 0);
		cigar_string_ss << len << op;
		if (op != 'D') cigar_length += len;
    if (op != 'S') all_soft_padded = false;
	}
	string cigar_string = cigar_string_ss.str();

  ////
  //// Don't write if this side turns out to be all soft-padded!  Issue #146
  ////
  if (all_soft_padded) return;

	////
	//// Assemble the quality score string
	////
	stringstream quality_score_ss;
	for (uint32_t i = 0; i < qual_scores.size(); i++)
  {
		quality_score_ss << static_cast<char>(qual_scores[i]);
	}
	string quality_score_string = quality_score_ss.str();

	////
	//// Setup custom aux tags
	////

	stringstream aux_tags_ss;

  uint32_t as;
  bool AS_found = a.aux_get_i("AS", as);
  ASSERT(AS_found, "Could not find required tag AS for alignment.");

	aux_tags_ss << "AS:i:" << as << "\t" << "X1:i:" << 1 << "\t" << "X2:i:" << fastq_file_index;

	//this flag indicates this is a junction match and which side of the match is in the middle of the read across the junction
	int32_t within_side = (reference_strand == 1) ? junction_side : (junction_side + 1) % 2;
	aux_tags_ss << "\t" << "XJ:i:" << within_side;

	// XL/XR trims for this split read. The junction (middle-of-read) side is never trimmed --
	// it is a trustworthy boundary, not a read end -- so it gets 0. The outer genomic side is
	// trimmed from the real-genome trims_list at this read's final (post soft-clip) coordinates,
	// exactly as write_alignments() does for a whole read (genomic trim + soft-clip offset).
	if ((trims_list != NULL) && (ref_seq_info_ptr != NULL) && (reference_match_start >= 1))
	{
		bool junction_on_left = ((junction_side == 2) == (read_strand == 1));
		uint32_t tid = ref_seq_info_ptr->seq_id_to_index(seq_id);
		uint32_t ref_span = 0;
		for (size_t c = 0; c < cigar_list.size(); c++) {
			char op = cigar_list[c].first;
			if (op=='M'||op=='D'||op=='N'||op=='='||op=='X') ref_span += cigar_list[c].second;
		}
		uint32_t new_ref_start_0 = static_cast<uint32_t>(reference_match_start) - 1;
		uint32_t new_ref_end_0 = new_ref_start_0 + ref_span - 1;
		uint32_t lead_S = (cigar_list.front().first == 'S') ? cigar_list.front().second : 0;
		uint32_t trail_S = (cigar_list.back().first == 'S') ? cigar_list.back().second : 0;
		uint32_t xl = junction_on_left ? 0 : ((*trims_list)[tid].left_trim_0(new_ref_start_0) + lead_S);
		uint32_t xr = junction_on_left ? ((*trims_list)[tid].right_trim_0(new_ref_end_0) + trail_S) : 0;
		aux_tags_ss << "\t" << "XL:i:" << xl << "\t" << "XR:i:" << xr;
	}

	string aux_tags = aux_tags_ss.str();

	////
	//// Create the BAM line and write it
	////

  vector<string> ll = make_vector<string>
		(a.read_name() + "-M" + to_string(junction_side))
		(to_string(fix_flags(a.flag())))
		(seq_id)
		(to_string(reference_match_start))
    ("255")
		(cigar_string)
		(a.proper_pair() ? "=" : "*")
		(to_string(a.mate_start_1()))
		(to_string(a.insert_size()))
		(seq)
		(quality_score_string)
		(aux_tags)
	;
	string l = join(ll, "\t") + "\n";
	if (verbose) cout << l;

	assert(cigar_length == q_length);
  {
    kstring_t ks = KS_INITIALIZE;
    kputs(l.c_str(), &ks);
    bam1_t* b = bam_init1();
    ASSERT(sam_parse1(&ks, bam_header, b) >= 0, "sam_parse1 failed: " + l);
    ASSERT(sam_write1(m_bam_file, bam_header, b) >= 0, "sam_write1 failed");
    bam_destroy1(b);
    ks_free(&ks);
  }
}

void bam_file::write_split_alignment(uint32_t min_indel_split_len, const alignment_wrapper& a, const alignment_list& alignments, const cReferenceSequences& ref_seq_info)
{
  bool debug = false;

  uint32_t q_length = a.read_length();
  string qseq_string = a.read_char_sequence();

  string quality_score_string = a.read_base_quality_char_string();
  if (quality_score_string[0] == ' ') {
    quality_score_string = alignments.read_base_quality_char_string;
    if (alignments.read_base_quality_char_string_reversed ^ a.reversed())
      quality_score_string = reverse_string(quality_score_string);
  }
  ASSERT(quality_score_string.size() > 0, "Attempt to write read with no quality scores: " + a.read_name());

  uint32_t rpos = a.reference_start_1();
  uint32_t qpos = a.query_start_1();

  vector<pair<char,uint16_t> > cigar_list = a.cigar_pair_char_op_array();
  vector<pair<char,uint16_t> > split_cigar_list;

  char op;
  uint32_t len;
  uint32_t i = 0;
  while (i < cigar_list.size())
  {
    uint32_t rstart = rpos;
    uint32_t qstart = qpos;

    string cigar_string = "";
    while (i < cigar_list.size())
    {
      op = cigar_list[i].first;
      len = cigar_list[i].second;

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
        split_cigar_list.push_back( make_pair(op, len) );
      }
      i++;
    }

    if (qpos > qstart)
    {
      // Save old values and use these if we move things
      uint32_t ins_updated_qpos = qpos;
      uint32_t ins_updated_rpos = rpos;

      // Extend the match if the inserted bases match the reference genome!
      // including past the inserted bases if necessary
      if (op == 'I') {
        if (debug) cout << "Testing insertion.." << endl;
        uint32_t tid = a.reference_target_id();
        while ( (ins_updated_qpos < a.query_end_1() )
               && (qseq_string[ins_updated_qpos - 1] == ref_seq_info[tid].get_sequence_1(ins_updated_rpos))  )  {
          ins_updated_qpos++;
          ins_updated_rpos++;
          ASSERT(split_cigar_list[split_cigar_list.size()-1].first =='M', "Attempt to modify nonmatch cigar operation.");
          split_cigar_list[split_cigar_list.size()-1].second +=1;
        }
      }

      //add padding to the sides of the match
      uint32_t left_padding = qstart - 1;
      uint32_t right_padding = q_length - ins_updated_qpos + 1;

      cigar_string = alignment_wrapper::cigar_op_array_to_cigar_string(split_cigar_list);
      split_cigar_list.clear();

      cigar_string = ( left_padding > 0 ? to_string(left_padding) + "S" : "" ) + cigar_string + ( right_padding > 0 ? to_string(right_padding) + "S" : "" );

      vector<string> ll = make_vector<string>
        (a.read_name())
        (to_string(a.flag()))
        (to_string(this->bam_header->target_name[a.reference_target_id()]))
        (to_string(rstart))
        (to_string(a.mapping_quality()))
        (cigar_string)
        ("*")("0")("0")
        (qseq_string)
        (quality_score_string)
      ;

      {
        string sam_line = join(ll, "\t");
        kstring_t ks = KS_INITIALIZE;
        kputs(sam_line.c_str(), &ks);
        bam1_t* b = bam_init1();
        ASSERT(sam_parse1(&ks, bam_header, b) >= 0, "sam_parse1 failed: " + sam_line);
        ASSERT(sam_write1(m_bam_file, bam_header, b) >= 0, "sam_write1 failed");
        bam_destroy1(b);
        ks_free(&ks);
      }
    }

    // If the inserted region matches to the next part of the match, then we want to adjust where that begins
    // (written as qpos > len, not qpos - len >= 1, since qpos/len are unsigned and the insertion can be
    // as long as or longer than the read prefix before it -- qpos - len would then wrap around instead of
    // going negative, defeating the check and sending a huge position into substr() below)
    if ( (op == 'I') && (qpos > len) )
    {
      string previous_string = qseq_string.substr(qpos - len - 1, len);
      string insert_string = qseq_string.substr(qpos - 1, len);
      if (previous_string == insert_string)
      {
        rpos -= len;
        i++;
        ASSERT(cigar_list[i].first == 'M', "Expected match in cigar string when resolving repeat insert.");
        cigar_list[i].second += len;
      }
    }

    //move up to the next match position
    while (i < cigar_list.size())
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
      else if (op == 'M')
      {
        break;
      }
      i++;
    }
  }
}


}
