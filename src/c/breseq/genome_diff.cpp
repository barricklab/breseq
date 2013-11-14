/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2012 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#include "libbreseq/genome_diff.h"
#include "libbreseq/reference_sequence.h"
#include "libbreseq/candidate_junctions.h"
#include "libbreseq/output.h"

namespace breseq {

// Common keywords used for diff entries:
const char* TYPE="type";
const char* ID="id";
const char* PID="parent";
const char* SEQ_ID="seq_id";
const char* START="start";
const char* END="end";
const char* START_RANGE="start_range";
const char* END_RANGE="end_range";
const char* LEFT_OUTSIDE_COV="left_outside_cov";
const char* LEFT_INSIDE_COV="left_inside_cov";
const char* RIGHT_INSIDE_COV="right_inside_cov";
const char* RIGHT_OUTSIDE_COV="right_outside_cov";
const char* POSITION="position";
const char* INSERT_POSITION="insert_position";
const char* QUALITY="quality";
const char* POLYMORPHISM_QUALITY="polymorphism_quality";
const char* GENOTYPE_QUALITY="genotype_quality";
const char* REF_BASE="ref_base";
const char* NEW_BASE="new_base";
const char* FREQUENCY="frequency";
const char* REJECT="reject";
const char* REF_COV="ref_cov";
const char* NEW_COV="new_cov";
const char* TOT_COV="tot_cov";
const char* ERROR="error";
  
//For JC
const char* SIDE_1_SEQ_ID = "side_1_seq_id";
const char* SIDE_1_POSITION = "side_1_position";
const char* SIDE_1_STRAND = "side_1_strand";
const char* SIDE_1_REDUNDANT = "side_1_redundant";
  
const char* SIDE_2_SEQ_ID = "side_2_seq_id";
const char* SIDE_2_POSITION = "side_2_position";
const char* SIDE_2_STRAND = "side_2_strand";
const char* SIDE_2_REDUNDANT = "side_2_redundant";

const char* SIDE_1_READ_COUNT="side_1_read_count";
const char* SIDE_2_READ_COUNT="side_2_read_count";
const char* NEW_JUNCTION_READ_COUNT="new_junction_read_count";
const char* NEW_JUNCTION_FREQUENCY = "new_junction_frequency";

const char* SIDE_1_COVERAGE = "side_1_coverage";
const char* SIDE_2_COVERAGE = "side_2_coverage";
const char* NEW_JUNCTION_COVERAGE = "new_junction_coverage";

  
map<gd_entry_type, vector<string> > line_specification = make_map<gd_entry_type, vector<string> >
//! seq_id and positions are already parameters in cDiffEntry
//## mutations
(SNP,make_vector<string> ("seq_id")("position")("new_seq"))
(SUB,make_vector<string> ("seq_id")("position")("size")("new_seq"))
(DEL,make_vector<string> ("seq_id")("position")("size"))
(INS,make_vector<string> ("seq_id")("position")("new_seq"))
(MOB,make_vector<string> ("seq_id")("position")("repeat_name")("strand")("duplication_size"))
(INV,make_vector<string> ("seq_id")("position")("size"))
(AMP,make_vector<string> ("seq_id")("position")("size")("new_copy_number"))
(CON,make_vector<string> ("seq_id")("position")("size")("region"))

//## evidence
(RA,make_vector<string> ("seq_id")("position")("insert_position")("ref_base")("new_base"))
(MC,make_vector<string> ("seq_id")("start")("end")("start_range")("end_range"))
(JC,make_vector<string> ("side_1_seq_id")("side_1_position")("side_1_strand")("side_2_seq_id")("side_2_position")("side_2_strand")("overlap"))
(CN,make_vector<string> ("seq_id")("start")("end")("copy_number"))
(UN,make_vector<string> ("seq_id")("start")("end"))

//## validation
(CURA,make_vector<string> ("expert"))
(FPOS,make_vector<string> ("expert"))
(PHYL,make_vector<string> ("gd"))
(TSEQ,make_vector<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
(PFLP,make_vector<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
(RFLP,make_vector<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
(PFGE,make_vector<string> ("seq_id")("enzyme"))
(NOTE,make_vector<string> ("note"))
(MASK,make_vector<string> ("seq_id")("position")("size")) 

  
; // end line specifications


const vector<string>gd_entry_type_lookup_table =
  make_vector<string>("UNKNOWN")("SNP")("SUB")("DEL")("INS")("MOB")("AMP")("INV")("CON")("RA")("MC")("JC")("CN")("UN")("CURA")("FPOS")("PHYL")("TSEQ")("PFLP")("RFLP")("PFGE")("NOTE")("MASK");


// Field order.
static const char* s_field_order[] = { 
  SEQ_ID,
  POSITION,
  INSERT_POSITION,
  REF_BASE,
  NEW_BASE,
  START,
  END,
  START_RANGE,
  END_RANGE,
0}; // required trailing null.


/*! Constructor.
 */
  
cDiffEntry::cDiffEntry()
: _type(UNKNOWN)
, _id("")
, _evidence()
{
}  
  
cDiffEntry::cDiffEntry(const gd_entry_type type)
: _type(type)
, _id("")
, _evidence()
{
}

cDiffEntry::cDiffEntry(const string &line, uint32_t line_number, cFileParseErrors* file_parse_errors)
  : _type(UNKNOWN)
  ,_id("")
  ,_evidence()
{
  cDiffEntry& de = *this;
  vector<string> tokens = split(line, "\t");

  if (tokens.size() < 3) {
    if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Could not determine type, id, or parent_id.", true);
    return;
  }

  uint32_t COLUMN = 0;
  //Type.
  string type = tokens[COLUMN];
  RemoveLeadingTrailingWhitespace(type);
  de._type = cDiffEntry::type_to_enum(type);

  if (de._type == UNKNOWN) {
    if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Unknown type for entry.", true);
    return;
  }
  ++COLUMN;

  //Id.
  de._id = tokens[COLUMN];
  RemoveLeadingTrailingWhitespace(de._id);
  ++COLUMN;

  //Evidence.
  if (de.is_mutation()) {
    de._evidence = split(tokens[COLUMN], ",");
    for (uint32_t i = 0; i < de._evidence.size(); ++i) {
      RemoveLeadingTrailingWhitespace(de._evidence[i]);
    }
  }
  ++COLUMN;

  //Specs.
  const vector<string>& specs = line_specification[de._type];

  if (tokens.size() < specs.size() ) {
    if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Expected " + breseq::to_string(specs.size()) + " tab-delimited fixed columns for entry", true);
    return;
  }
  
  for (uint32_t i = 0; i < specs.size(); ++i) {
    if (COLUMN < tokens.size()) {
      de[specs[i]] = tokens[COLUMN];
      RemoveLeadingTrailingWhitespace(de[specs[i]]);
      ++COLUMN;
    } else {
      
      if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Line specification [" + specs[i] + "] has no value.", true);
      return;
    }
  }

  //Fields.
  while(COLUMN < tokens.size() && tokens[tokens.size() - 1].size()) {
    cKeyValuePair kvp(tokens[COLUMN], '=');
    if (kvp.valid()) {
      string key   = kvp.get_key();
      string value = kvp.get_value();
      RemoveLeadingTrailingWhitespace(key);
      RemoveLeadingTrailingWhitespace(value);
      de[key] = value;
    } else {
      if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Field " + kvp + " is not a key=value pair. Ignoring this key.", false);
    }
    ++COLUMN;
  }

  return;
}
  
gd_entry_type cDiffEntry::type_to_enum(string type) {
  for (size_t i = 0; i < gd_entry_type_lookup_table.size(); i++) {
    if (type == gd_entry_type_lookup_table[i]) {
      return (gd_entry_type) i;
    }
  }
  return UNKNOWN;
}
  
//Comparing IDs here will currently break cGenomeDiff::merge and cGenomeDiff::subtract
bool cDiffEntry::operator==(const cDiffEntry& de)
{
  //! Case: Easy if not same type
  if (this->_type != de._type) {
    return false;
  }
  //! Case: Same type, but are fields that are common equal?
  else {
    // Get common keys
    const vector<diff_entry_key_t>& specs = line_specification[this->_type];
    for(vector<diff_entry_key_t>::const_iterator it = specs.begin(); it != specs.end(); it++) {
      const diff_entry_key_t& spec(*it);
      if ((*this)[spec] != de.find(spec)->second)
        return false;
    }
    return true;
  }
}

bool cDiffEntry::is_mutation() const
{
  return gd_entry_type_lookup_table[_type].size() == 3;
}

bool cDiffEntry::is_evidence() const
{
  return gd_entry_type_lookup_table[_type].size() == 2;
}

bool cDiffEntry::is_validation() const
{
  return gd_entry_type_lookup_table[_type].size() == 4;
}
 
  
uint32_t cDiffEntry::get_start()
{
  switch (this->_type) {
    case SNP:
    case SUB:
    case DEL:
    case INS:
    case MOB:
    case INV:
    case AMP:
    case CON:
      return from_string<uint32_t>((*this)[POSITION]);
    case UN:
      return from_string<uint32_t>((*this)[START]);
    case MASK:
    default:
      ERROR("cDiffEntry::get_start not implemented for type: " + gd_entry_type_lookup_table[this->_type]);
  }
  return 0;
}
  
uint32_t cDiffEntry::get_end()
{
  switch (this->_type) {
    case SNP:
      return from_string<uint32_t>((*this)[POSITION]);
    case SUB:
    case DEL:
    case INV:
    case AMP:
    case CON:
      return from_string<uint32_t>((*this)[POSITION]) + from_string<uint32_t>((*this)[SIZE]);
    case INS:
      return from_string<uint32_t>((*this)[POSITION]) + (*this)[NEW_SEQ].length();
    case MOB:
      return from_string<uint32_t>((*this)[POSITION]) + from_string<uint32_t>((*this)["duplication_size"]);
    case UN:
      return from_string<uint32_t>((*this)[END]);
    case MASK:
    default:
      ERROR("cDiffEntry::get_end not implemented for type: " + gd_entry_type_lookup_table[this->_type]);
  }
  return 0;
}
  
int32_t cDiffEntry::mutation_size_change(cReferenceSequences& ref_seq_info)
{  
  switch (this->_type) {
    case SNP: return 0; break;
    case SUB: return (*this)[NEW_SEQ].length() - from_string<uint32_t>((*this)[SIZE]); break;
    case INS: return (*this)[NEW_SEQ].length(); break;
    case DEL: return -(from_string<uint32_t>((*this)[SIZE])); break;
    case AMP: return from_string<uint32_t>((*this)[SIZE]) * (from_string<uint32_t>((*this)["new_copy_number"]) - 1); break;
    case MASK: return 0; break;
      
    case MOB:
    {
      int32_t size = from_string<int32_t>((*this)["repeat_size"]) + from_string<int32_t>((*this)["duplication_size"]);
      if (this->entry_exists("del_start"))
        size -= from_string<uint32_t>((*this)["del_start"]);
      if (this->entry_exists("del_end"))
        size -= from_string<uint32_t>((*this)["del_end"]);
      if (this->entry_exists("ins_start"))
        size += (*this)["ins_start"].length();
      if (this->entry_exists("ins_end"))
        size += (*this)["ins_end"].length();
      return size;
      break;
    }
      
    case CON:
    {
      uint32_t replace_target_id, replace_start, replace_end;
      ref_seq_info.parse_region((*this)["region"], replace_target_id, replace_start, replace_end);  
      int32_t size = from_string<uint32_t>((*this)[SIZE]);
      
      int32_t new_size = (replace_end > replace_start) ? replace_end - replace_start + 1 : replace_start - replace_end + 1;
      return  new_size - size;
      break;
    }
    default:
      ASSERT(false, "Unable to calculate mutation size change.");
      return UNDEFINED_INT32;
  }
}
  
/*! Marshal this diff entry into an ordered list of fields.
 */
void cDiffEntry::marshal(vector<string>& s) const {
  s.push_back(gd_entry_type_lookup_table[_type]);
  s.push_back(_id);
  
  s.push_back(join(_evidence, ","));
  
  // deep copy all fields:
  cDiffEntry cp= *this;
  
  // marshal specified fields in-order, removing them from the copy after they've 
  // been printed:
  
  vector<string>& f = line_specification[_type];
  
  for (vector<string>::iterator it=f.begin(); it != f.end(); it++)
  {
    diff_entry_map_t::iterator iter=cp.find(*it);
    
    ASSERT(iter != cp.end(), "Did not find required field '" + *it + "' to write in entry id " + _id + " of type '" + gd_entry_type_lookup_table[_type] + "'.");
    
    s.push_back(iter->second);
    cp.erase(iter);
    
  }
  
  // marshal whatever's left, unless it's an empty field or _begins with an underscore
  for(diff_entry_map_t::iterator i=cp.begin(); i!=cp.end(); ++i) {
    
    assert(i->first.size());
    if (i->first.substr(0,1) == "_") continue;
    
    if (i->second.size() > 0) 
    {
      s.push_back(i->first + "=" + i->second);
    }
  }
}

// Created the line to be printed
string cDiffEntry::to_string(void) const
{
  vector<string> fields;
  marshal(fields);
  return join(fields, "\t");
}
 
size_t cDiffEntry::number_reject_reasons()
{
  if(this->entry_exists(REJECT))
  {
    return this->get_reject_reasons().size();
  }
  return 0;
}

vector<string> cDiffEntry::get_reject_reasons()
{
  vector<string> return_value;
  if (this->entry_exists(REJECT)) {
    return_value = split((*this)[REJECT], ",");
  } 
  return return_value;
}

/*! Add reject reason to diff entry.
 */
void cDiffEntry::add_reject_reason(const string &reason) {

  if (this->find(REJECT) == this->end()) {
    (*this)[REJECT] = reason;
  }
  // exists already, make comma separated list
  else {
    string reject = (*this)[REJECT];
    reject += ",";
    reject +=reason; 
  }
}
  
  
//! Remove all information except required fields
cDiffEntry cDiffEntry::to_spec() const
{
  cDiffEntry de(_type);
  
  const vector<diff_entry_key_t>& specs = line_specification[_type];
  
  for(vector<diff_entry_key_t>::const_iterator it = specs.begin(); it != specs.end(); it++) {
    const diff_entry_key_t& spec(*it);
    de[spec] = this->get(spec);
  }
  
  return de;
}
  
void cDiffEntry::normalize_to_sequence(const cAnnotatedSequence &sequence, bool verbose)
{
  //! Step: Initialize parameters.
  //Only diff entries applicable to de sequence can be considered.
  assert(this->entry_exists("seq_id") || (*this)["seq_id"] == sequence.m_seq_id);
  assert(this->entry_exists("position"));
  
  //Sequences should be viewed as having index + one offset.
  const uint32_t pos_1 = strtoul((*this)["position"].c_str(), NULL, 0);
  if (!pos_1) {
    this->insert(pair<string, string>("comment_out", "True"));
    return;
  }
  
  /*! Step: Type specific normalizations. For some, the initial parameters given
   can be altered to be normalized, for others the parameters can't be altered
   and the mutation is not valid for the given reference sequence. */
  switch (this->_type)
  {
    case DEL:
    {
      /*     
       +++EXAMPLE DEL: pos_1 = 1, size = 3
       
       START:
       
       v pos_2
       Ref: ACTACTATCACACTAATACAATA
       ^ pos_1
       
       seq_1 = ACT
       seq_2 = ACT
       NOT VALID, ACT == ACT
       
       NEXT(pos_1++):
       
       v pos_2
       Ref: ACTACTATCACACTAATACAATA
       ^ pos_1
       
       seq_1 = CTA
       seq_2 = CTA
       NOT VALID, CTA == CTA
       
       NEXT(pos_1++):
       
       v pos_2
       Ref: ACTACTATCACACTAATACAATA
       ^ pos_1
       
       seq_1 = TAC
       seq_2 = TAT
       ^ first mismatch, normalize pos_1 to here.
       
       THEREFOR: pos_1 = 8. 
       */
      
      assert(this->entry_exists("size"));
      
      //! Step: Initialize parameters.
      typedef string::const_iterator base_itr_t;
      typedef pair<base_itr_t, base_itr_t> base_pair_t;
      const size_t n = strtoul((*this)["size"].c_str(), NULL, 0);
      
      if (!n) {
        this->insert(pair<string, string>("comment_out", "True"));
        return;
      }
      
      /*! Step: Attempt to normalize the start position by iterating through new
       start positions. */
      uint32_t norm_pos_1 = pos_1;
      for(int i = 0;;++i) {
        
        uint32_t seq1_pos_1 = pos_1 + i;
        uint32_t seq1_end_1 = seq1_pos_1 + n - 1;
        
        uint32_t seq2_pos_1 = seq1_end_1 + 1;
        uint32_t seq2_end_1 = seq2_pos_1 + n - 1;
        
        const string &seq1 =
        sequence.get_circular_sequence_1(seq1_pos_1, n);
        assert(seq1.size() == n);
        
        const string &seq2 =
        sequence.get_circular_sequence_1(seq2_pos_1, n);
        assert(seq2.size() == n);
        
        
        const base_pair_t &base_pair =
        mismatch(seq1.begin(), seq1.end(), seq2.begin());
        
        bool valid = base_pair.first != seq1.end();
        
        
        if (valid) {
          norm_pos_1 += i + (base_pair.first - seq1.begin());
          sprintf((*this)["position"], "%u", norm_pos_1);
          if (verbose) {
            sprintf((*this)["norm_pos"], "%u_to_%u", pos_1, norm_pos_1);
          }
        }
        
        if (verbose) {
          if (valid && norm_pos_1 == pos_1) {
            cerr << "VALID POSITION: " << pos_1 << endl;
          }
          else if (valid && norm_pos_1 != pos_1) {
            cerr << "NORMALIZED POSITON: " << pos_1 << " to " << norm_pos_1 << ' ' << endl;
          } else {
            cerr << "INVALID POSITION: " << pos_1 << endl;
          }
          
          cerr << "[Mutation]: " << this->to_spec().to_string() << endl;
          cerr << "Sequence 1 [" << seq1_pos_1 << '-' << seq1_end_1 << "]: " << seq1 << endl;
          
          cerr << "Sequence 2 [" << seq2_pos_1 << '-' << seq2_end_1 << "]: " << seq2 << endl;
          
          cerr << "Deleted Sequence [" << norm_pos_1 << '-' << norm_pos_1 + n - 1 << "]: "; 
          cerr << sequence.get_circular_sequence_1(norm_pos_1, n) << endl;
          
          string prev_seq = sequence.get_circular_sequence_1(norm_pos_1 - 10, 20);
          string left_seq = sequence.get_circular_sequence_1(norm_pos_1 - 10, 10);
          string right_seq = sequence.get_circular_sequence_1(norm_pos_1 + n , 10);
          
          cerr << "Previous Sequence:  " << prev_seq << endl;
          cerr << "Resulting Junction: " << left_seq << '\t' << right_seq << endl;
          
          cerr << endl;
        }
        
        if (valid) {
          break;
        }
        
      }
      
    } break;
      
    case INS:
    {
      /*
       +++EXAMPLE INS: pos_1 = 1, new_seq = CAA
       
       START:
       
       v pos_2
       Ref: ACTACTATCACACTAATACAATA
       ^ pos_1
       
       seq_1 = CAA
       seq_2 = CTA
       ^ match, could be predicted as a SNP, NOT VALID
       
       NEXT(pos_1++):
       Ref: ACTACTATCACACTAATACAATA
       ^ pos_1
       seq_1 = CAA
       seq_2 = TAC
       ^ mismatch, won't be predicted as a SNP, normalize pos_1 to here.
       
       THEREFOR: pos_1 = 2.
       
       
       +++EXAMPLE INS -> AMP: pos_1 = 1, new_seq = CTA
       
       START:
       Ref: ACTACTATCACACTAATACAATA
       ^ pos_1
       seq_1 = CTA
       seq_2 = CTA
       NOT VALID, CTA == CTA, possible AMP.
       
       NEXT(pos_1 += 3):
       Ref: ACTACTATCACACTAATACAATA
       ^ pos_1
       seq_1 = CTA
       seq_2 = CTA
       NOT VALID, CTA == CTA, possible AMP.
       
       NEXT(pos_1 += 3):
       Ref: ACTACTATCACACTAATACAATA
       ^ pos_1
       seq_1 = CTA
       seq_2 = TCA
       VALID, Passes SNP check and will be evaluated as an AMP.
       
       THEREFOR: INS->AMP, pos_1 = 1, size = 3, new_copy_number = 2, orig_type = INS.
       
       */
      assert(this->entry_exists("new_seq"));
      
      //! Step: Initialize parameters.
      const string first = (*this)["new_seq"];
      const size_t n = first.size();
      assert(n);
      
      /*! Step: Attempt to normalize the start position by iterating through new
       start positions. */
      uint32_t i = pos_1;
      
      bool bAmp, bSnp;
      bAmp = true;
      bSnp = true;
      
      while(bAmp || bSnp)
      {
        bAmp = true;
        bSnp = true;
        
        for(;;i += n)
        {
          const string second = sequence.get_circular_sequence_1(i + 1, n);
          assert(second.size());
          
          if (first != second)  {
            bAmp = false;
            break;  }
        }
        
        for(;;i += 1)
        {
          const string second = sequence.get_circular_sequence_1(i + 1, n);
          assert(second.size());
          
          if (first[0] != second[0])  {
            bSnp = false;
            break;  }
        }
      }
      
      //! Step: Determine if the start pos needs to be normilized.
      bool is_new_position = pos_1 != i;
      if (is_new_position) {
        sprintf((*this)["position"], "%u", i);
        sprintf((*this)["norm_pos"], "%u_to_%u", pos_1, i);
      }
      
      //! Step: Determine if this is actually an AMP.
      if ((sequence.get_circular_sequence_1(i - (n - 1), n) == first)  && (n > 1)) {
        this->_type = AMP;
        sprintf((*this)["size"], "%u", n);
        sprintf((*this)["new_copy_number"], "2");
        sprintf((*this)["orig_type"], "INS");
      }
    } break;
      
    case SNP:
    {
      assert(this->entry_exists("new_seq"));
      assert((*this)["new_seq"].size() == 1);
      
      //! Step: Initialize parameters.
      const base_char first = (*this)["new_seq"][0];
      
      //! Step: Compare bases.
      const base_char second = sequence.get_circular_sequence_1(pos_1, 1)[0];
      
      //! Step: If bases are not equal then it's not a valid snp.
      bool is_base_not_equal = (first != second);
      if (!is_base_not_equal) {
        sprintf((*this)["norm"], "is_not_valid");
      }
    } break;
      
    case SUB:
    {
      assert(this->entry_exists("size"));
      assert(this->entry_exists("new_seq"));
      
      //! Step: Initialize parameters.
      const size_t n = strtoul((*this)["size"].c_str(), NULL, 0);
      assert(n);
      const string first = (*this)["new_seq"];
      assert(first.size());
      
      const string second = sequence.get_circular_sequence_1(pos_1, n);
      const string third = sequence.get_circular_sequence_1(pos_1 + n, 1);
      
      if (first.find_first_of(second) != string::npos &&
          first.find_first_of(third) != string::npos ) {
        sprintf((*this)["norm"], "is_not_valid");
      }
    } break;
      
    case MOB:
    {
      assert(this->entry_exists("repeat_name"));
      assert(this->entry_exists("strand"));
      assert(this->entry_exists("duplication_size"));
    } break;
      
    case INV:
    {
      assert(this->entry_exists("size"));
      const size_t n =  strtoul((*this)["size"].c_str(), NULL, 0);
      assert(n);
    } break;
      
    case AMP:
    {
      assert(this->entry_exists("size"));
      assert(this->entry_exists("new_copy_number"));
    } break;
      
    case CON:
    {
      assert(this->entry_exists("size"));
      assert(this->entry_exists("region"));
    } break;
      
      
    default:
      break;
      
  }//End switch.
  
  return;
}


  

/*! Constructor.
 */
cGenomeDiff::cGenomeDiff(const string& filename)
 : _default_filename(filename)
 , _unique_id_counter(0) 
{
 read(filename);  
}

/*! Merge Constructor.
 */
cGenomeDiff::cGenomeDiff(cGenomeDiff& merge1, cGenomeDiff& merge2, bool unique, bool new_id, bool verbose)
 : _unique_id_counter(0)
{
  this->merge(merge1, unique, new_id, verbose);
  this->merge(merge2, unique, new_id, verbose);
}
  
diff_entry_ptr_t cGenomeDiff::parent(const cDiffEntry& evidence)
{
  for(diff_entry_list_t::iterator itr_test_item = _entry_list.begin();
      itr_test_item != _entry_list.end(); itr_test_item ++) { 
    
    cDiffEntry& test_item = **itr_test_item;
    
    for(vector<string>::iterator itr = test_item._evidence.begin();
        itr != test_item._evidence.end(); itr ++) { 
      string& test_evidence_id = (*itr);
      if(test_evidence_id == evidence._id) {      
        return diff_entry_ptr_t(*itr_test_item);
      }
    }
  }
  return diff_entry_ptr_t(NULL);
}
  
/*! Read a genome diff(.gd) from the given file to class member
 _entry_list
 */
cFileParseErrors cGenomeDiff::read(const string& filename, bool suppress_errors) {
  ifstream in(filename.c_str());
  ASSERT(in.good(), "Could not open file for reading: " + filename);
  uint32_t line_number = 1;
  cFileParseErrors parse_errors(filename);
  
  //! Step: Handle header parameters.
  //Example header:
  //#=GENOME_DIFF 1.0
  //#=AUTHOR    Barrick JE
  //#=REFSEQ    Genbank:NC_012967.1
  //#=READSEQ   SRA:SRR066595
  
  /*#=GENOME_DIFF version must be initialized for this file to be recognized
   as a genome diff file. */
  metadata.version = "";
  while(in.peek() == '#') {
    in.get();
    if (in.peek() != '=') {
      in.unget();
      break;
    } else {
      in.unget();
    }
    
    string whole_line = "";
    string second_half = "";
    
    getline(in, whole_line);
    line_number++;
    
    vector<string> split_line = split_on_whitespace(whole_line);
    
    if(split_line.size() > 1)second_half = split_line[1];
    for(size_t j = 2; j < split_line.size(); j++)  {
      second_half += " " + split_line[j];  }    
    
    //Evaluate key.
    if (split_line[0] == "#=GENOME_DIFF" && split_line.size() > 1) {
      metadata.version = second_half;
    }
    else if (split_line[0] == "#=AUTHOR") { 
      metadata.author = second_half;
    }
    else if (split_line[0] == "#=REFSEQ" && split_line.size() > 1) {
      metadata.ref_seqs.resize(metadata.ref_seqs.size() + 1);
      metadata.ref_seqs.back() = second_half;
    }
    else if (split_line[0] == "#=READSEQ" && split_line.size() > 1) {
      metadata.read_seqs.resize(metadata.read_seqs.size() + 1);
      metadata.read_seqs.back() = second_half;
    }
    else if (split_line[0] == "#=TITLE" && split_line.size() > 1) {
      metadata.run_name = second_half;
      replace(metadata.run_name.begin(), metadata.run_name.end(), ' ', '_');
    }
    else if (split_line[0].substr(0, 2) == "#=" && split_line.size() > 1) {                                                                     
      string key = split_line[0].substr(2, split_line[0].size());
      this->add_breseq_data(key, second_half);
      continue;
    } else {
      //Warn if unknown header lines are encountered.
      parse_errors.add_line_error(line_number, whole_line, "Header line not recognized and will be ignored.", false);
    }
  }
  
  /*Error if #=GENOME_DIFF is not found. Genome diff files are required to have
   this header line. */
  
  if (metadata.version.empty()) {
    
    parse_errors.add_line_error(0,"", "No #=GENOME_DIFF XX header line in this file.", true);
    if (!suppress_errors) {
      parse_errors.print_errors();
      exit(0);
    } else { 
      return parse_errors;
    }
  }
  
  /*If the run_name/title is not set by a #=TITLE tag in the header info then
   default to the name of the file. */
  if (metadata.run_name.empty()) {
    string *run_name = &metadata.run_name;
    *run_name = filename.substr(0, filename.find(".gd"));
    if (run_name->find('/') !=  string::npos) {
      run_name->erase(0, run_name->find_last_of('/') + 1);
    }
  }
  
  //! Step: Handle the diff entries.
  while (in.good()) {
    string line = "";
    std::getline(in, line);
    line_number++;
    //Warn if commented out or a possibly blank line is encountered.
    if (line.empty()) {
      continue;
    } else if (line[0] == '#') {
      //printf("Discarding Genome Diff comment file:%s line:\n%s\n", filename.c_str(), line.c_str());
      continue;
    } else if (line.find_first_not_of(' ') == string::npos) {
      continue;
    }
    const cDiffEntry de(line, line_number, &parse_errors);
    if (de._type != UNKNOWN) add(de);
    
  }
  
  if (!suppress_errors) {
    parse_errors.print_errors();
    if (parse_errors.fatal() )
    {
      ERROR("Parse errors in loading GenomeDiff File. Not safe to continue");
      exit(1);
    }
  }
  
  return parse_errors;
}

/*! Write this genome diff to a file.
 NOTES:
 1) If you want a diff entry to be commented out(prefix with '#') add the key
 "comment_out" to the diff entry.
 */
void cGenomeDiff::write(const string& filename) {
  string basename = cString(filename).get_base_name();
  string dir = cString(filename).remove_ending(basename);
  if (dir.size()) {
    create_path(dir);
  }
  ofstream os(filename.c_str());
  
  
  //! Step: Header lines.
  /*Always write version tag. It's how we identify it as a genome diff file
   in cGenomeDiff::read(). */
  fprintf(os, "#=GENOME_DIFF\t%s\n", metadata.version.c_str());
  
  fprintf(os, "#=AUTHOR\t%s\n", metadata.author.c_str());
  
  for (size_t i = 0; i < metadata.ref_seqs.size(); i++) {
    fprintf(os, "#=REFSEQ\t%s\n", metadata.ref_seqs[i].c_str());
  }
  
  for (size_t i = 0; i < metadata.read_seqs.size(); i++) {
    fprintf(os, "#=READSEQ\t%s\n", metadata.read_seqs[i].c_str());
  }
  if (!metadata.breseq_data.empty()) {
    for (map<key_t,string>::iterator it = metadata.breseq_data.begin();
         it != metadata.breseq_data.end(); it ++) {
      fprintf(os, "#=%s\t%s\n", it->first.c_str(), it->second.c_str());
    }
  }
  
  // sort
  this->sort();
  
  for(diff_entry_list_t::iterator it=_entry_list.begin(); it!=_entry_list.end(); ++it) {
    if (!(*it)->entry_exists("comment_out")) {
      fprintf(os, "%s\n", (**it).to_string().c_str());
    } else {
      (*it)->erase("comment_out");
      fprintf(os, "#%s\n", (**it).to_string().c_str());
    }
  }
  os.close();
}
  
//! Call to assure that every entry in a cGenomeDiff
//  has the same SEQ_ID as the ref_seq_info.
//  This will also try and check entry positions and
//  sizes against the length of the ref_seq_info.
bool cGenomeDiff::is_valid(cReferenceSequences& ref_seq_info, bool verbose)
{  
  //Go through every cDiffEntry for this cGenomeDiff
  for (diff_entry_list_t::iterator entry = _entry_list.begin(); entry != _entry_list.end(); entry++)
  {    
    //Grab the information based on type for this entry
    const list_t spec = line_specification[(**(entry))._type];
    
    //For every entry of this type go through each of their fields
    for(size_t i = 0; i < spec.size(); i++)
    {
      //Does the current field contain a SEQ_ID?
      if(spec[i].find(SEQ_ID) != string::npos)
      {
        bool no_match = true;
        vector<string> failure;
        
        for (vector<cAnnotatedSequence>::iterator it_as = ref_seq_info.begin(); it_as < ref_seq_info.end(); it_as++)
        {
          //If the current entry has any SEQ_ID fields, does that SEQ_ID
          // field match any from ref_seq_info?
          if((**(entry))[spec[i]] == it_as->m_seq_id)
            no_match = false;
          
          failure.push_back(it_as->m_seq_id);
        }
        
        //If we couldn't find any matching SEQ_IDs, perform notification and failure.
        if(no_match)
        {
          cout << "LOADED REFERENCE\t" << *failure.begin() << endl;
          for(vector<string>::iterator it_fail = ++failure.begin(); it_fail < failure.end(); it_fail++)cout << "\t\t\t\t\t" << *it_fail << endl;
          cout << "LOADED GENOMEDIFF\t" << (**(entry))[spec[i]] << " ID=" << (**(entry))._id << endl;
          return false;
        }
      }
      
      string temp_seq_id = "";
      
      //Find the previous SEQ_ID; the one associated with this position
      // To do this, we iterate backwards from the current field.
      for(int32_t x = 0 + i; x >= 0; --x)
      {
        //Does the current entry have a SEQ_ID field?
        if(spec[x].find(SEQ_ID) != string::npos)
        {
          temp_seq_id = (**(entry))[spec[x]];
          break;
        }
      }
      
      //Does the current entry have any POSITION fields?
      // WARNING!
      // Position check might fail if the entry has multiple
      // position fields as well as a size field.
      if(spec[i].find(POSITION) != string::npos)
      {        
        int32_t temp_pos = from_string<int32_t>((**(entry))[spec[i]]);
        
        //Grab the information based on type for this entry
        const list_t temp_spec = line_specification[(**(entry))._type];
        
        //For every entry of this type go through each of their fields
        for(size_t u = 0 + i; u < temp_spec.size(); u++)
        {
          //Does the current entry have a SIZE field?
          if(spec[u] == SIZE)
          {
            //Add the size field to the position, size should include the current
            // position, so subtract 1.
            temp_pos += from_string<int32_t>((**(entry))[temp_spec[u]]) - 1;
          }
        }        
        
        for(vector<cAnnotatedSequence>::iterator it_as = ref_seq_info.begin(); it_as < ref_seq_info.end(); it_as++)
        {
          //Does the position in this field, plus any field that might indicate
          //size, exceed the length of the relevant sequence in ref_seq_info?
          if(it_as->m_seq_id == temp_seq_id && temp_pos > it_as->m_length)
          {
            cout << "LOADED REFERENCE\t" << it_as->m_seq_id << "\tLENGTH:\t" << it_as->m_length << endl;
            cout << "LOADED GENOMEDIFF\t" << temp_seq_id << "\tID=" << (**(entry))._id << "\t" << temp_pos << " POTENTIAL POSITION" << endl;
            return false;
          }
        }
      }
      
      //Does the current entry have a START field?
      if(spec[i] == START)
      {        
        for(vector<cAnnotatedSequence>::iterator it_as = ref_seq_info.begin(); it_as < ref_seq_info.end(); it_as++)
        {
          //Does this entry have a START larger than ref_seq_info length?
          if(it_as->m_seq_id == temp_seq_id && from_string<int32_t>((**(entry))[spec[i]]) > it_as->m_length)
          {
            cout << "LOADED REFERENCE\t" << it_as->m_seq_id << "\tLENGTH:\t" << it_as->m_length << endl;
            cout << "LOADED GENOMEDIFF\t" << temp_seq_id << "\tID=" << (**(entry))._id << "\t" << (**(entry))[spec[i]] << " START" << endl;
            return false;
          }
        }
      }
      
      //Does the current entry have a END field?
      if(spec[i] == END)
      {        
        for(vector<cAnnotatedSequence>::iterator it_as = ref_seq_info.begin(); it_as < ref_seq_info.end(); it_as++)
        {
          //Does this entry have a END larger than ref_seq_info length?
          if(it_as->m_seq_id == temp_seq_id && from_string<int32_t>((**(entry))[spec[i]]) > it_as->m_length)
          {
            cout << "LOADED REFERENCE\t" << it_as->m_seq_id << "\tLENGTH:\t" << it_as->m_length << endl;
            cout << "LOADED GENOMEDIFF\t" << temp_seq_id << "\tID=" << (**(entry))._id << "\t" << (**(entry))[spec[i]] << " END" << endl;
            return false;
          }
        }        
      }
    }
  }
  
  //Notify the user that the files match.
  if(verbose)
    cout << "** LOADED FILES MATCH **" << endl;
  
  return true;
}

/*! Find the next unused unique id.
 */
uint32_t cGenomeDiff::new_unique_id()
{ 
  uint32_t assigned_id = ++_unique_id_counter;
  
  while (unique_id_used.count(assigned_id))
  {
    assigned_id++;
  }
  return assigned_id;
}

  
/*! Add evidence to this genome diff.
 */
diff_entry_ptr_t cGenomeDiff::add(const cDiffEntry& item, bool lowest_unique) {
  
  ASSERT(item._type != UNKNOWN, "Tried to add item of type UNKNOWN to Genome Diff file.");
  
  // allocating counted_ptr takes care of disposal
  cDiffEntry* diff_entry_copy = new cDiffEntry(item);
  diff_entry_ptr_t added_item(diff_entry_copy);
  _entry_list.push_back(added_item);
  
  uint32_t new_id = 0;    
  
  if ((added_item->_id.size() == 0) || lowest_unique || unique_id_used.count(from_string<uint32_t>(added_item->_id)))  {
    new_id = new_unique_id();
    added_item->_id = to_string(new_id);
  }
  else  {
    new_id = from_string<uint32_t>(added_item->_id);
  }  
  
  unique_id_used[new_id] = true;
  
  return added_item;
}
  
void cGenomeDiff::remove(cGenomeDiff::group group)
{
  this->sort();
  diff_entry_list_t::iterator it1 = _entry_list.begin();
  diff_entry_list_t::iterator it2 = _entry_list.begin();
  
  //Mutations.
  while (it2 != _entry_list.end()) {
    if (!(**it2).is_mutation()) break;
    ++it2;
  }
  if (group == cGenomeDiff::MUTATIONS) {
    _entry_list.erase(it1, it2);
    return;
  }
  
  //Evidence.
  it1 = it2; 
  while (it2 != _entry_list.end()) {
    if (!(**it2).is_evidence()) break;
    ++it2;
  }
  if (group == cGenomeDiff::EVIDENCE) {
    _entry_list.erase(it1, it2);
    return;
  }
  
  //Validation.
  it1 = it2; 
  while (it2 != _entry_list.end()) {
    if (!(**it2).is_validation()) break;
    ++it2;
  }
  if (group == cGenomeDiff::VALIDATION) {
    _entry_list.erase(it1, it2);
    return;
  }
  
  return;
}

/*! Given an id return the entry if it exists. NULL otherwise.
 */ 

diff_entry_ptr_t cGenomeDiff::find_by_id(string _id)
{
  for (diff_entry_list_t::iterator itr_diff_entry = _entry_list.begin();
       itr_diff_entry != _entry_list.end(); itr_diff_entry++)
  {
    if ( (*itr_diff_entry)->_id == _id)
      return *itr_diff_entry;
  }
  return diff_entry_ptr_t(NULL);
}


  
  
/*! Given a list of types, search and return the cDiffEntry's within diff_entry_list_t whose 
 * _type parameter matches one of those within input types. 
 */ 
diff_entry_list_t cGenomeDiff::list(const vector<gd_entry_type>& types)
{
  // default is to have to types
  if (types.size() == 0)
    return _entry_list;
  
  diff_entry_list_t return_list;
  
  for (diff_entry_list_t::iterator itr_diff_entry = _entry_list.begin();
       itr_diff_entry != _entry_list.end(); itr_diff_entry++)
  {
    for (vector<gd_entry_type>::const_iterator requested_type = types.begin();
         requested_type != types.end(); requested_type++)
    {
      if((*itr_diff_entry)->_type == *requested_type)
        return_list.push_back(*itr_diff_entry);
    }
  }
  
  return return_list;
}


diff_entry_list_t cGenomeDiff::show_list(const vector<gd_entry_type>& types)
{
  diff_entry_list_t ret_list = list(types);
  ret_list.remove_if(cDiffEntry::fields_exist(make_vector<diff_entry_key_t>("deleted")));
  ret_list.remove_if(cDiffEntry::fields_exist(make_vector<diff_entry_key_t>("no_show")));
  return ret_list;
}
  
  
/*! Return all cDiffEntrys within _entry_list whose _id matches one
 * of those within input's item._evidence
 */ 
diff_entry_list_t cGenomeDiff::mutation_evidence_list(const cDiffEntry& item)
{
  diff_entry_list_t return_list;
  
  //return diff_entries with matching evidence
  for (vector<string>::const_iterator itr_i = item._evidence.begin(); itr_i != item._evidence.end(); itr_i ++) 
  {  
    const string& evidence = *itr_i;
    
    for (diff_entry_list_t::iterator itr_j = _entry_list.begin(); itr_j != _entry_list.end(); itr_j ++)
    {  
      cDiffEntry& entry = **itr_j;
      
      if (entry._id == evidence)
        return_list.push_back(*itr_j);
    }   
  }
  return return_list;
}
  
// return items with types that are 3 characters long
diff_entry_list_t cGenomeDiff::mutation_list()
{
  this->sort();
  diff_entry_list_t::iterator it = _entry_list.begin();
  while (it != _entry_list.end()) {
    if (!(**it).is_mutation()) {
      break;
    } else {
      ++it;
    }
  }
  
  return diff_entry_list_t(_entry_list.begin(), it);
}

// return items with types that are 2 characters long
diff_entry_list_t cGenomeDiff::evidence_list()
{
  diff_entry_list_t mut_list;
  
  for(diff_entry_list_t::iterator itr = _entry_list.begin();
      itr != _entry_list.end(); itr ++) {
    cDiffEntry& item = **itr;
    if(item.is_evidence()) {
      mut_list.push_back(*itr);
    }
  }
  
  return mut_list;
}

// return items with types that are 4 characters long
diff_entry_list_t cGenomeDiff::validation_list()
{
  this->sort();
  diff_entry_list_t::reverse_iterator it = _entry_list.rend();
  while (it != _entry_list.rbegin()) {
    if (!(**it).is_validation()) {
      --it;
      break;
    } else {
      ++it; 
    }
  }
  return diff_entry_list_t(it.base(), _entry_list.end());
}

//! Removes all GD entries that aren't used as evidence.
void cGenomeDiff::filter_not_used_as_evidence(bool verbose)
{
  // Yes, I know the bool is useless.
  map<string,bool> used_evidence;
  
  diff_entry_list_t muts = this->mutation_list();
  //Iterate through all mutations
  for (diff_entry_list_t::iterator it = muts.begin(); it != muts.end(); it++)
  {    
    //For every piece of evidence this entry has
    for(uint32_t iter = 0; iter < (**it)._evidence.size(); iter++)
    {
      //Each piece of evidence will only get assigned to the map once.
      used_evidence[(**it)._evidence[iter]] = true;
    }
  }
  
  // We will be erasing inside the it loop.  This is to keep
  // track of whether or not we should iterate to the next element.
  bool it_iterate = true;
  
  //Iterate through all entries
  for (diff_entry_list_t::iterator it = _entry_list.begin(); it != _entry_list.end(); )
  {
    bool it_iterate = true;
    
    //Is this ID in our map of used_evidence?
    if(!(used_evidence.count((**it)._id)) && (**it).is_evidence())
    {
      //Inform the user
      if(verbose){cout << "NOT USED AS EVIDENCE: " << (**it)._id << endl;}
      
      //Remove this entry from the list.
      it = _entry_list.erase(it);
      
      //We just removed the current feauture, do not iterate.
      it_iterate = false;
    }
    
    // Iterate it ONLY if we haven't erased something.
    if(it_iterate)it++;
  }
}

/*-----------------------------------------------------------------------------
 * returns entries NOT used as evidence by other entries. 
 *
 *-----------------------------------------------------------------------------*/
diff_entry_list_t cGenomeDiff::filter_used_as_evidence(const diff_entry_list_t& input)
{
  // first we make a map with everything used as evidence by any entry in the entire genome diff
  map<string,bool> used_as_evidence;
  for (diff_entry_list_t::const_iterator it = _entry_list.begin(); it != _entry_list.end(); it++)
  {
    const diff_entry_ptr_t& de = *it;
    for (vector<string>::const_iterator ev_it = de->_evidence.begin(); ev_it != de->_evidence.end(); ev_it++) 
    {  
      used_as_evidence[*ev_it] = true;
    }   
  }
  
  // then construct a list of all items in input with ids not in this map
  diff_entry_list_t return_list;
  for (diff_entry_list_t::const_iterator it = input.begin(); it != input.end(); it++)
  {
    const diff_entry_ptr_t& de = *it;
    if ( !used_as_evidence.count(de->_id) )
      return_list.push_back(de);
  }
  
  return return_list;
}

// Helper function for mutation_unknown and mutation_deleted
bool cGenomeDiff::mutation_in_entry_of_type(cDiffEntry mut, const gd_entry_type type)
{
  uint32_t start = mut.get_start();
  uint32_t end = mut.get_end();
  
  diff_entry_list_t check_list = list(make_vector<gd_entry_type>(type));
  
  for (diff_entry_list_t::iterator itr = check_list.begin(); itr != check_list.end(); itr++) {
    
    cDiffEntry de(**itr);
    if ( (start >= de.get_start()) && (end <= de.get_end()) ) {
      return true;
    }
  }
  return false;
}

//! Subtract mutations using gd_ref as reference.
void cGenomeDiff::set_subtract(cGenomeDiff& gd_ref, bool verbose)
{
  // We will be erasing inside the it loop.  This is to keep
  // track of whether or not we should iterate to the next element.
  bool it_iterate = true;
  
  //Iterate through all the entries
  for (diff_entry_list_t::iterator it = _entry_list.begin(); it != _entry_list.end(); )
  {
    it_iterate = true;
    //The current entry we're looking at
    cDiffEntry& entry = **it;
    
    //if (verbose) cout << entry << endl;
    
    //Is the entry a mutation?
    if(entry.is_mutation())
    {
      //Iterate through all the entries we're checking against.
      for (diff_entry_list_t::iterator it_ref = gd_ref._entry_list.begin(); it_ref != gd_ref._entry_list.end(); it_ref++)
      {
        //The current entry we're looking at
        cDiffEntry& entry_ref = **it_ref;
        
        if(!entry_ref.is_mutation())
          continue;
        
        //if (verbose) cout << "  " << entry_ref << endl;
        
        //Does the current entry match any of the reference entries?
        if(entry == entry_ref)
        {
          //Notify the user of the action.
          if(verbose){cout << "REMOVE\t" << to_string(entry._type) << "\t" << entry._id << endl;}
          it = _entry_list.erase(it);
          //We just removed the current feauture, do not iterate.
          it_iterate = false;
          break; // Done comparing to this mutaion.
        }
      }
    }
    
    // Iterate it ONLY if we haven't erased something.
    if(it_iterate)it++;
  }
}
  
void cGenomeDiff::set_intersect(cGenomeDiff &gd, bool verbose)
{
  (void)verbose; //unused
  
  set<cDiffEntry> seen;
  diff_entry_list_t muts = gd.mutation_list();
  
  for (diff_entry_list_t::iterator it = muts.begin();
       it != muts.end(); ++it) {
    seen.insert(**it);
  }
  
  this->sort();
  set<string> ids;
  diff_entry_list_t::iterator it = _entry_list.begin();
  //Handle mutations, store evidence id of deleted mutations
  //to later delete.
  while (it != _entry_list.end()) {
    if (!(**it).is_mutation()) break;
    
    if (!seen.count(**it)) {
      for (uint32_t i = 0; i < (**it)._evidence.size(); ++i) {
        ids.insert((**it)._evidence[i]); 
      }
      it = _entry_list.erase(it);
    } else {
      ++it;
    }
  }
  //Delete evidence that matches it.
  while (it != _entry_list.end()) {
    if (ids.count((**it)._id)) {
      it = _entry_list.erase(it);
    } else {
      ++it;
    }
  }
}
  
void cGenomeDiff::set_union(cGenomeDiff& gd, bool verbose)
{
  (void)verbose; //unused
  this->fast_merge(gd);
  this->remove(cGenomeDiff::EVIDENCE);
  this->remove(cGenomeDiff::VALIDATION);
  this->unique(); 
}

// Keeps only one copy when equivalent entries are encountered
void cGenomeDiff::unique()
{
  bool (*comp_fn) (const diff_entry_ptr_t&, const diff_entry_ptr_t&) = diff_entry_ptr_sort;
  typedef set<diff_entry_ptr_t, bool(*)(const diff_entry_ptr_t&, const diff_entry_ptr_t&)> diff_entry_set_t;
  //Filter entries.
  diff_entry_set_t seen(comp_fn);
  //Store pointers to mutations.
  map<string, vector<diff_entry_ptr_t> > keep_ids;
  //Store ids of evidence to erase.
  set<string> erase_ids;
  
  this->sort();
  diff_entry_list_t::iterator it = _entry_list.begin();
  while (it != _entry_list.end()) {
    if (!(**it).is_mutation()) break;
    
    const vector<string>& ids = (**it)._evidence;
    uint32_t n = ids.size();
    
    //Is mutation unique?
    //Case: true.
    if (seen.insert(*it).second) { 
      for (uint32_t i = 0; i < n; ++i) {
        keep_ids[ids[i]].push_back(*it);
      }
      ++it;
    } 
    //Case: false.
    else { 
      for (uint32_t i = 0; i < n; ++i) {
        erase_ids.insert(ids[i]);
      }
      it = _entry_list.erase(it);
    }
  }
  
  seen.clear(); //Re-use to filter the evidence.
  while (it != _entry_list.end()) {
    //Keep this evidence?
    //Case: unkown.
    if (keep_ids.count((**it)._id) && erase_ids.count((**it)._id)) {
      //Is evidence unique?
      //Case: true.
      if (seen.insert(*it).second) {
        ++it;
      } 
      //Case: false.
      else {
        it = _entry_list.erase(it);
        //Update mutations that may of been using this id.
        for (uint32_t i = 0; i < keep_ids[(**it)._id].size(); ++it) {
          vector<string>* evidence = &(*keep_ids[(**it)._id][i])._evidence;
          vector<string>::iterator jt = std::remove(evidence->begin(), evidence->end(), (**it)._id);
          //evidence->erase(jt);
        }
      }
    } 
    //Case: false.
    else if (!keep_ids.count((**it)._id) && erase_ids.count((**it)._id)) {
      it = _entry_list.erase(it);
    } 
    //Case: false.
    else if (!keep_ids.count((**it)._id) && !erase_ids.count((**it)._id)) {
      stringstream ss;
      ss << "\tRemoving [entry]:\t" << **it << endl;
      ss << "\tfrom [file]:\t" << this->_default_filename << endl;
      ss << "\tbecause no mutation referenced it's ID." << endl; 
      WARN(ss.str());
      it = _entry_list.erase(it);
    }
    
    //Case: true.
    else {
      ++it;
    }
  }
  
  return;
}



//! Merge GenomeDiff information using gd_new as potential new info.
//  Evidence IDs that are not unique are given new IDs.  Mutations
//  that refer to this evidence have their evidence updated as well.
//
//  bool unique:  TRUE will NOT merge entries that match existing entries.
//                FALSE WILL merge entries that match existing entries.
//
//  bool new_id:  TRUE will give assign all new entries with the lowest available ID.
//                FALSE will allow all entries to retain their IDs if they are unique.
void cGenomeDiff::merge(cGenomeDiff& gd_new, bool unique, bool new_id, bool verbose)
{
  uint32_t old_unique_ids = unique_id_used.size();
  
  //Iterate through all the potential new entries
  for (diff_entry_list_t::iterator it_new = gd_new._entry_list.begin(); it_new != gd_new._entry_list.end(); it_new++)
  {
    //The current potential new entry we're looking at
    cDiffEntry& entry_new = **it_new;
    bool new_entry = true;
    
    //Iterate through all the entries in the current list.
    for (diff_entry_list_t::iterator it_cur = _entry_list.begin(); it_cur != _entry_list.end() && unique; it_cur++)
    {
      //The current entry we're looking at
      cDiffEntry& entry = **it_cur;
      
      //Does the new entry match the current entry?
      if(entry == entry_new)
      {
        //Existing matching entry found, this is not new
        new_entry = false;
        break;
      }
    }
    
    //We definately have a new entry
    if(new_entry)
    {      
      //Add the new entry to the existing list
      add(entry_new, new_id);
      
      //Notify user of new entry
      if(verbose)cout << "\tNEW ENTRY\t" << entry_new._id << " >> " << _entry_list.back()->_id << "\t" << gd_new._default_filename << endl;
    }
  }
  
  //Iterate through all the entries in the new list.
  //This is where we update the evidence IDs for mutations.
  for (diff_entry_list_t::iterator it = _entry_list.begin(); it != _entry_list.end(); it++)
  {
    // @JEB: optimization: we don't need to do this for evidence items.
    if ( (*it)->is_evidence() ) continue;
    
    //Is this one of the new entries?
    if(from_string<uint32_t>((**it)._id) > old_unique_ids && (**it).is_mutation())
    {                
      //For every piece of evidence this entry has
      for(int32_t iter = 0; iter < (int32_t)(**it)._evidence.size(); iter++)
      {
        bool found_match = false;
        
        //Iterate through all the potential new entries
        for (diff_entry_list_t::iterator it_new = gd_new._entry_list.begin(); it_new != gd_new._entry_list.end(); it_new++)
        {            
          //Does this evidence ID match an ID in the old list?
          if((**it)._evidence[iter] == (**it_new)._id && !found_match)
          {
            //Iterate through all the current entries
            for (diff_entry_list_t::iterator it_cur =_entry_list.begin(); it_cur != _entry_list.end(); it_cur++)
            {
              //Does the new entry match the current entry?
              if((**it_cur) == (**it_new))
              {
                //Notify user of the update
                if(verbose)cout << "\tID: " << (**it)._id  << "\tEVIDENCE\t" << (**it)._evidence[iter] << " >> " << (**it_cur)._id << endl;
                
                //Change the evidence ID to it's new ID in the new updated list
                (**it)._evidence[iter] = (**it_cur)._id;
                found_match = true;  
                break;
              }
            }
          }
        }
        
        // If we've gone through all the lists, and can't find the evidence
        // then remove the evidence entry completely, as it matches to nothing.
        if(!found_match)
        {
          //Notify user of the update
          if(verbose)cout << "\tID: " << (**it)._id  << "\tEVIDENCE  \t" << (**it)._evidence[iter] << " >> REMOVED" << endl;
          
          (**it)._evidence.erase((**it)._evidence.begin() + iter--);
        }        
      }
    }
  }
  
  //Notify user of the update
  if(verbose)cout << "\tMERGE DONE - " << gd_new._default_filename << endl;
  
}


void cGenomeDiff::fast_merge(const cGenomeDiff& gd)
{  
  diff_entry_list_t gd_list = gd.list();
  for(diff_entry_list_t::const_iterator it=gd_list.begin(); it!= gd_list.end(); it++) {
    this->add(**it);
  }
}
  
// This function is for simplifying the numbering
void cGenomeDiff::reassign_unique_ids()
{
  this->sort();
  
  //Handle mutations.
  _unique_id_counter = 0;
  
  map<string, vector<diff_entry_ptr_t> > id_table;
  diff_entry_list_t::iterator it = _entry_list.begin();
  while (it != _entry_list.end()) {
    if (!(**it).is_mutation()) break;
    (**it)._id = to_string(++_unique_id_counter);
    
    for (uint32_t i = 0; i < (**it)._evidence.size(); ++i) {
      id_table[(**it)._evidence[i]].push_back(*it);
    }
    (**it)._evidence.clear();
    
    ++it;
  }
  
  //Handle the rest.
  while (it != _entry_list.end()) {
    string new_id = to_string(++_unique_id_counter);
    
    if (id_table.count((**it)._id)) {
      for (uint32_t i = 0; i < id_table[(**it)._id].size(); ++i) {
        id_table[(**it)._id][i]->_evidence.push_back(new_id);
      }
    }
    
    (**it)._id = new_id;
    ++it;
  }
}

  
// All fields must be assigned in this table and be required fields of the gd entries.
map<gd_entry_type, cGenomeDiff::sort_fields_item> diff_entry_sort_fields = make_map<gd_entry_type, cGenomeDiff::sort_fields_item>
  (SNP, cGenomeDiff::sort_fields_item(1, SEQ_ID, POSITION))
  (SUB, cGenomeDiff::sort_fields_item(1, SEQ_ID, POSITION))
  (DEL, cGenomeDiff::sort_fields_item(1, SEQ_ID, POSITION))
  (INS, cGenomeDiff::sort_fields_item(1, SEQ_ID, POSITION))
  (MOB, cGenomeDiff::sort_fields_item(1, SEQ_ID, POSITION))
  (AMP, cGenomeDiff::sort_fields_item(1, SEQ_ID, POSITION))
  (INV, cGenomeDiff::sort_fields_item(1, SEQ_ID, POSITION))
  (CON, cGenomeDiff::sort_fields_item(1, SEQ_ID, POSITION))
  (RA,  cGenomeDiff::sort_fields_item(2, SEQ_ID, POSITION))
  (MC,  cGenomeDiff::sort_fields_item(3, SEQ_ID, START))
  (JC,  cGenomeDiff::sort_fields_item(4, "side_1_seq_id", "side_1_position"))
  (CN,  cGenomeDiff::sort_fields_item(5, SEQ_ID, START))
  (UN,  cGenomeDiff::sort_fields_item(6, SEQ_ID, START))
  (CURA, cGenomeDiff::sort_fields_item(7, "expert", "expert"))
  (FPOS, cGenomeDiff::sort_fields_item(7, "expert", "expert"))
  (PHYL, cGenomeDiff::sort_fields_item(7, "gd", "gd"))
  (TSEQ, cGenomeDiff::sort_fields_item(7, "seq_id", "primer_1_start"))
  (PFLP, cGenomeDiff::sort_fields_item(7, "seq_id", "primer_1_start"))
  (RFLP, cGenomeDiff::sort_fields_item(7, "seq_id", "primer_1_start"))
  (PFGE, cGenomeDiff::sort_fields_item(7, "seq_id", "enzyme"))
  (NOTE, cGenomeDiff::sort_fields_item(7, "note", "note"))
  (MASK, cGenomeDiff::sort_fields_item(7, SEQ_ID, POSITION))
;

map<gd_entry_type, uint8_t> sort_order = make_map<gd_entry_type, uint8_t>
  (SNP, 2)
  (SUB, 4)
  (DEL, 1)
  (INS, 3)
  (MOB, 5)
  (AMP, 6)
  (INV, 7)
  (CON, 8)
  (RA,  9)
  (MC,  10)
  (JC,  11)
  (CN,  12)  
  (UN,  13)
  (CURA, 14)
  (FPOS, 15)
  (PHYL, 16)
  (TSEQ, 17)
  (PFLP, 18)
  (RFLP, 19)
  (PFGE, 20)
  (NOTE, 20)
  (MASK, 20)
;


/*! Write this genome diff to a file.
 */
bool cGenomeDiff::diff_entry_ptr_sort(const diff_entry_ptr_t& a, const diff_entry_ptr_t& b) {

  gd_entry_type a_type = a->_type;
  gd_entry_type b_type = b->_type;

  sort_fields_item a_sort_fields = diff_entry_sort_fields[a_type];
  sort_fields_item b_sort_fields = diff_entry_sort_fields[b_type];
  
  
  if (a_sort_fields._f1 < b_sort_fields._f1) {
    return true;
  } else if (a_sort_fields._f1 > b_sort_fields._f1) {
    return false;
  }
  
  string a_sort_field_2 = (*a)[a_sort_fields._f2];
  string b_sort_field_2 = (*b)[b_sort_fields._f2];
  
  if (a_sort_field_2 < b_sort_field_2) {
    return true;
  } else if (a_sort_field_2 > b_sort_field_2) {
    return false;
  }  

  uint32_t a_sort_field_3 = from_string<uint32_t>((*a)[a_sort_fields._f3]);
  uint32_t b_sort_field_3 = from_string<uint32_t>((*b)[b_sort_fields._f3]);
  
  if (a_sort_field_3 < b_sort_field_3) {
    return true;
  } else if (a_sort_field_3 > b_sort_field_3) {
    return false;
  }  
  
  uint8_t a_sort_order = sort_order[a_type];
  uint8_t b_sort_order = sort_order[b_type];

  if (a_sort_order < b_sort_order) {
    return true;
  } else if (a_sort_order > b_sort_order) {
    return false;
  } 
  
  // last sort by id
  uint32_t a_sort_id = from_string(a->_id);
  uint32_t b_sort_id = from_string(b->_id);

  if (a_sort_id < b_sort_id) {
    return true;
  } else if (a_sort_id > b_sort_id) {
    return false;
  } 
  
  return false;
}

  
  

//! Call to generate random mutations.
void cGenomeDiff::random_mutations(string exclusion_file,
                                   string type,
                                   uint32_t n_muts,
                                   uint32_t buffer,
                                   cAnnotatedSequence& ref,
                                   bool verbose)
{
  //Parse input option into mutation type.
  //Also determine minimum size and maximum size if provided by user. 
  vector<string> type_options = split_on_any(type, ":-,");
  string mut_type = type_options[0];
  uint32_t min_size = 1, max_size = 1;
  uint32_t min_copy_number = 2, max_copy_number = 2;
  switch(type_options.size())
  {
    case 1:  {
    }  break;
      
    case 2:  {
      min_size = un(type_options[1]);
      max_size = un(type_options[1]);
    }  break;
      
    case 3:  {
      min_size = un(type_options[1]);
      max_size = un(type_options[2]);
    }  break;
    
    case 4:  {
      min_size = un(type_options[1]);
      max_size = un(type_options[2]);
      min_copy_number = un(type_options[3]);
      max_copy_number = un(type_options[3]);
    }  break;
    
    case 5:  {
      min_size = un(type_options[1]);
      max_size = un(type_options[2]);
      min_copy_number = un(type_options[3]);
      max_copy_number = un(type_options[4]);
    }  break;
      
    default:      
      ERROR("CANNOT PARSE: " + type);
  }

  // These are repeat regions that we want to avoid when placing mutations
  cFlaggedRegions repeat_match_regions;
  if (exclusion_file.size()) {
    repeat_match_regions.read_mummer(exclusion_file, ref);
  }


  //Container for checking that future simulated mutations are a buffered distance apart.
  cFlaggedRegions used_mutation_regions;

  const uint32_t max_attempts = 1000;
  uint32_t n_attempts = max_attempts;

  srand(cSimFastqSequence::SEED_VALUE);
  buffer +=1;

  /* TYPICAL WORKFLOW:
   *    1) Chose a random mutation size.
   *
   *    2) Chose a random position.
   *        a) Create and normalize a mutation at that position.
   *        b) Check that the normalized position and size will not overlap an excluded region plus buffer.
   *        c) Check that the normalized position and size not be within buffer distance of 
   *          another mutation.
   *        d) If necessary, check for any mutation specific changes after normalization.
   *
   *    3) Repeat #2 until a valid position is found or the max number of attempts for that 
   *      mutation size is reached.
   *        a) If the max number of attempts is reached go back to #1.
   *        b) If valid, add mutation to GD and to used_mutation_regions.
   *
   *    4) Repeat #1 until N valid mutations are found or the max number of attempts for simulating
   *      N mutations is reached.
   *        a) If the max number of attempts is reached then no more mutations can be simulated.
   *
   */

  if (mut_type == "SNP" || mut_type == "INS" || mut_type == "DEL") {      
    
    while (n_muts && n_attempts) {

      cDiffEntry new_item;
      new_item._type = cDiffEntry::type_to_enum(mut_type);
      new_item["seq_id"] = ref.m_seq_id;

      uint32_t pos_1 = 0;
      uint32_t size = mut_type == "SNP" ? 0 : rand() % (max_size - min_size + 1) + min_size;
      uint32_t n_size_attempts = max_attempts;

      while (n_size_attempts) {
        pos_1 = (rand() % (ref.get_sequence_size() - buffer - size)) + buffer;
        new_item["position"] = to_string(pos_1);

        if (mut_type == "SNP") {
          new_item["new_seq"] = cSimFastqSequence::get_random_error_base(ref.get_sequence_1(pos_1));        
        } 
        else 
        if (mut_type == "INS") {
          string* new_seq = &new_item["new_seq"];
          new_seq->resize(size);
          generate(new_seq->begin(), new_seq->end(), cSimFastqSequence::get_random_insertion_base);
        }
        else
        if (mut_type == "DEL") {
          new_item["size"] = s(size);
        }

        new_item.normalize_to_sequence(ref, verbose);

        pos_1 = un(new_item["position"]);

        bool is_excluded = repeat_match_regions.is_flagged(pos_1 - buffer, pos_1 + size + buffer);
        bool is_near_mutation = used_mutation_regions.is_flagged(pos_1 - buffer, pos_1 + size + buffer);

        if (is_excluded || is_near_mutation) {
          --n_size_attempts;
        } else {
          break;
        }
      }

      if (n_size_attempts == 0) {
        --n_attempts;
      } else {
        --n_muts, n_attempts = max_attempts;
        this->add(new_item);
        used_mutation_regions.flag_region(pos_1, pos_1 + size);
      }

    }
  }
  else if (mut_type == "RMD") {
    //Get all available IS elements.
    cSequenceFeatureList repeats = ref.m_repeats;
    ASSERT(repeats.size(), "No repeat_regions / ISX elements in reference sequence.");
    CHECK(n_muts <= repeats.size(), "Too many deletions requested, creating a potential maximum of " + s(repeats.size()));


    //Unflag repeat-match regions which affect an IS element and then store them.
    cFlaggedRegions IS_element_regions;
    for (cSequenceFeatureList::iterator it = repeats.begin(); it != repeats.end(); ++it) {
        uint32_t start_1 = (*it)->get_start_1(), end_1   = (*it)->get_end_1();
        IS_element_regions.flag_region(start_1, end_1 + 1); //Want end_1 + 1 because a deletion will occur there.



        cFlaggedRegions::regions_t regions = repeat_match_regions.regions(start_1, end_1);

        if (regions.size()) {
            uint32_t lower = regions.begin()->first; 
            uint32_t upper = regions.rbegin()->second;
            
            cerr << "\tRemoving repeat-match excluded region: " << lower << "-" << upper << endl;
            cerr << "\t\tFor IS element: " << (**it)["name"] << "\t" << (*it)->get_start_1() << "-" << (*it)->get_end_1() << endl;
            cerr << endl;

            repeat_match_regions.remove(regions);

          }
    }

    while (n_muts && repeats.size() && n_attempts) {
      vector<cDiffEntry> valid_items;
      cSequenceFeatureList::iterator it;

      //Randomly choose deletion size.
      uint32_t pos_1 = 0;
      uint32_t size = rand() % (max_size - min_size + 1) + min_size;
      uint32_t n_size_attempts = max_attempts;
      while (n_size_attempts) {
        //Randomly choose a repeat_region.
        it = repeats.begin();
        advance(it, rand() % repeats.size());

        //Collect potentially valid left_side and right_side positions.
        vector<int32_t> valid_pos_1;
        uint32_t start_1 = (*it)->get_start_1(), end_1   = (*it)->get_end_1();
        valid_pos_1.push_back(start_1 - size);
        valid_pos_1.push_back(end_1 + 1);

        //Evaluate that once normalized as DELs, the potentially valid mutations are not within an excluded region.
        for (vector<int32_t>::iterator jt = valid_pos_1.begin(); jt != valid_pos_1.end(); ++jt) {
          pos_1 = *jt ; 

          cDiffEntry temp_item;        
          temp_item._type = DEL;
          temp_item["seq_id"]   = ref.m_seq_id;
          temp_item["position"] = s(pos_1);        
          temp_item["size"]     = s(size);        

          bool overlaps_multiple_IS = IS_element_regions.regions(pos_1, pos_1 + size).size() > 1;
          bool is_excluded = repeat_match_regions.is_flagged(pos_1 - buffer, pos_1 + size + buffer);
          bool is_within_buffer = used_mutation_regions.is_flagged(pos_1 - buffer, pos_1 + size + buffer);

          if (!is_excluded && !is_within_buffer && !overlaps_multiple_IS) {
            valid_items.push_back(temp_item);
          }
        }

        if (valid_items.empty()) {
          --n_size_attempts;
        } else { 
          break;
        }
      }

      if (n_size_attempts == 0) {
        --n_attempts;
      } else {
        --n_muts, n_attempts = max_attempts;
        cDiffEntry new_item = valid_items[rand() % valid_items.size()];
        new_item["mediated"] = (**it)["name"];
        this->add(new_item);
        pos_1 = un(new_item["position"]);
        size = un(new_item["size"]);
        used_mutation_regions.flag_region(pos_1, pos_1 + size);
        repeat_match_regions.flag_region((*it)->get_start_1(), (*it)->get_end_1());

        if (verbose) {
          cerr << "[ISX]: " + (**it)["name"] << "\t[start_1]: " << (*it)->get_start_1() << "\t[end_1]: " << (*it)->get_end_1() << endl;
          cerr << "\t[DEL]: " << new_item << endl;
          cerr << endl;
        }
        repeats.erase(it);

      }

    }
    //Validate mutations.

    //Check that neighboring mutations are a buffered distance from each other.
    diff_entry_list_t muts = this->mutation_list();
    diff_entry_list_t::iterator it = muts.begin();
    diff_entry_list_t::iterator jt = it;
    for (++jt; jt != muts.end(); ++it, ++jt) {
      //Left mutation.
      uint32_t lpos_1 = un((**it)["position"]);
      uint32_t lsize = un((**it)["size"]);

      //Right mutation.
      uint32_t rpos_1 = un((**jt)["position"]);

      ASSERT(static_cast<uint32_t>(abs(static_cast<int32_t>(lpos_1 + lsize - rpos_1))) >=  buffer,
          "Mutation: " + (*it)->to_spec().to_string() + "\n" +
          "\tand\n" +
          "Mutation: " + (*it)->to_spec().to_string() +"\n" +
          "are less then the buffered distance away from each other.");
    }

    //Check that each mutation only effects one IS element.
    for (it = muts.begin(); it != muts.end(); ++it) {
      uint32_t pos_1 = un((**it)["position"]);
      uint32_t size = un((**it)["size"]);
      uint32_t n_IS_elements = IS_element_regions.regions(pos_1, pos_1 + size).size();
      ASSERT(n_IS_elements == 1,
          "Mutation overlaps more then one IS element: " + (*it)->to_spec().to_string());
    }

    //Display IS elements that could not be used.
    for (cSequenceFeatureList::iterator it = repeats.begin(); it != repeats.end(); ++it) {
      WARN("Un-used IS element: " + (**it)["name"]+ '\t' + s((*it)->get_start_1()) + '-' + s((*it)->get_end_1()));
    }

  }
  else
  if (mut_type == "AMP") {
    while (n_muts && n_attempts) {
      uint32_t pos_1 = 0;
      uint32_t size = rand() % (max_size - min_size + 1) + min_size;
      
      cDiffEntry new_item;        
      new_item._type = AMP;
      new_item["seq_id"] = ref.m_seq_id;

      uint32_t n_size_attempts = max_attempts;

      while (n_size_attempts) {
        pos_1 = (rand() % (ref.get_sequence_size() - buffer - size)) + buffer;

        new_item["position"] = to_string(pos_1);
        new_item["size"] = s(size);
        new_item["new_copy_number"] = s((rand() % (max_copy_number - min_copy_number + 1)) + min_copy_number);

        new_item.normalize_to_sequence(ref);
        new_item.erase("norm_pos");
        uint32_t norm_pos_1 = un(new_item["position"]);

        bool is_excluded      = repeat_match_regions.is_flagged(pos_1 - buffer, pos_1 + size + buffer);
        bool is_near_mutation = used_mutation_regions.is_flagged(pos_1 - buffer, pos_1 + size + buffer);
        bool is_new_pos_1     = pos_1 != norm_pos_1;
        bool is_new_type      = to_string(new_item._type) != mut_type;

        if (is_excluded  || is_near_mutation || is_new_pos_1 || is_new_type) {
          --n_size_attempts;
        } else {
          break;
        }

      }

      if (n_size_attempts == 0) {
        --n_attempts;
      } else {
        --n_muts, n_attempts = max_attempts;
        new_item.erase("norm");
        this->add(new_item);
        used_mutation_regions.flag_region(pos_1, pos_1 + (size * un(new_item["new_copy_number"])));
        if (verbose) {
          cerr << "\t" << new_item << endl;
        }

      }

    }

  }
  else 
  if (mut_type == "MOB") {
    //@JEB: Update.. it does. But this comment from Matt is too fun to remove.
    //ERROR("THE BRESEQENSTEIN MUTATION GENERATOR DOES NOT YET HANDLE MOBS\nESPECIALLY IF THEY HAVE TORCHES AND PITCHFORKS");

    cSequenceFeatureList repeats = ref.m_repeats;
    ASSERT(repeats.size(), "No repeat_regions / ISX elements in reference sequence.");

    while (n_muts && n_attempts) {
      uint32_t pos_1 = 0;
      uint32_t size = rand() % (max_size - min_size + 1) + min_size;
      int8_t strand = rand() % 2 ? 1 : -1;
      cSequenceFeatureList:: iterator it = repeats.begin();
      advance(it, rand() % repeats.size());

      uint32_t n_size_attempts = max_attempts;
      while (n_size_attempts) {
        pos_1 = (rand() % (ref.get_sequence_size() - buffer)) + buffer;
        string temp_seq = ref.get_sequence_1(pos_1, pos_1 + size - 1);
        cDiffEntry temp_item;        
        temp_item._type = INS;
        temp_item["seq_id"] = ref.m_seq_id;
        temp_item["position"] = s(pos_1);
        temp_item["new_seq"] = temp_seq;
        
        temp_item.normalize_to_sequence(ref);
        uint32_t norm_pos_1   = un(temp_item["position"]);

        bool is_excluded      = repeat_match_regions.is_flagged(norm_pos_1 - buffer, norm_pos_1 + size + buffer);
        bool is_near_mutation = used_mutation_regions.is_flagged(norm_pos_1 - buffer, norm_pos_1 + size + buffer);
        bool is_new_pos_1     = pos_1 != norm_pos_1; 
        bool is_new_seq       = temp_seq != temp_item["new_seq"];
        bool is_not_INS       = temp_item._type != INS;

        if (is_excluded || is_near_mutation || is_new_pos_1 || is_new_seq || is_not_INS) {
          --n_size_attempts;
        } else {
          break;
        }

      }

      if (n_size_attempts == 0) {
        --n_attempts;
      } else {
        --n_muts, n_attempts = max_attempts;
        cDiffEntry new_item;        
        new_item._type = MOB;
        new_item["seq_id"] = ref.m_seq_id;
        new_item["position"] = s(pos_1);
        new_item["repeat_name"] = (*it)->SafeGet("name");
        new_item["strand"] = s(strand);        
        new_item["duplication_size"] = s(size);

        new_item.normalize_to_sequence(ref);
        new_item.erase("norm_pos");

        this->add(new_item);
        used_mutation_regions.flag_region(pos_1, pos_1 + size);

        if (verbose) {
          cerr << "\t" << new_item << endl;
        }

      }

    }
  }
  else { 
    ERROR("MUTATION TYPE NOT HANDLED: " + mut_type);
  }
  
  CHECK(max_attempts == n_attempts, "Forced to halt mutation generation.\nAttempted " +
        s(max_attempts - n_attempts) + " times to generate another mutation.\n" + 
        "It's likely that it's no longer possible to add new mutations.");
  
}
  
void cGenomeDiff::shift_positions(cDiffEntry &current_mut, cReferenceSequences& ref_seq_info, bool verbose)
{  
  int32_t delta = current_mut.mutation_size_change(ref_seq_info);
  if (verbose)
    cout << "Shift size: " << delta << endl;
  if (delta == UNDEFINED_INT32)
    WARN("Size change not defined for mutation.");
  
  uint32_t offset = from_string<uint32_t>(current_mut[POSITION]);
  
  // Only offset if past the duplicated part of a MOB (allows putting a mutation in the first copy)
  if (current_mut._type == MOB) {
    offset += from_string<uint32_t>(current_mut["duplication_size"]);
  }
  
  diff_entry_list_t muts = this->mutation_list();
  for (diff_entry_list_t::iterator itr = muts.begin(); itr != muts.end(); itr++) {
    cDiffEntry& mut = **itr;
    
    if (mut._type == INV)
      ERROR("shift_positions() cannot handle inversions yet!");
  

    // Check to see whether this is listed as being nested within the current item
    cDiffEntry* nested_within(NULL);
    int32_t nested_copy = -1;
    if (mut.entry_exists("nested_id")) {
      
      if (current_mut._id == mut["nested_id"]) {
        
        nested_within = &current_mut;
        
        ASSERT(current_mut._type == AMP, "Attempt to use \"nested\" to put mutation in non-AMP mutation.\n" + current_mut.to_string() + "\nreferenced from\n" + mut.to_string());
        ASSERT(current_mut[SEQ_ID] == mut[SEQ_ID], "Nested mutation not in same reference sequence?");
        
        if (mut.entry_exists("nested_copy")) {
          nested_copy = from_string<int32_t>(mut["nested_copy"]);
          int32_t new_copy_number = from_string<int32_t>(current_mut["new_copy_number"]);
          ASSERT((nested_copy > 0) && (nested_copy <= new_copy_number), "Copy number of mutation nested within is not in a valid range for the number of copies that exist. Requested copy: " + to_string(nested_copy) + " Total copies: " + to_string(new_copy_number));
        } else {
          ERROR("Expected \"nested_copy=\" additional field to be present in item with \"nexted_id\":\n" + mut.to_string())
        }
      }
    }
    
    
    // Special case, we are nested within current_mut
    uint32_t position = from_string<uint32_t>(mut[POSITION]);
    if (nested_within) {
      
      uint64_t special_delta = from_string<uint64_t>(current_mut[SIZE]) * (nested_copy-1);
      mut[POSITION] = to_string(position + special_delta);
      
    // Normal behavior -- offset mutations later in same reference sequence
    } else {
      if (current_mut[SEQ_ID] == mut[SEQ_ID] && position > offset)
        mut[POSITION] = to_string(position + delta);
    }

    
  }
}

// The return sequence replaces the target site duplication
string cGenomeDiff::mob_replace_sequence(cReferenceSequences& ref_seq_info, cDiffEntry& mut)
{

  ASSERT(mut._type == MOB, "Attempt to get mob_replace_sequence for non-MOB diff entry.");
  
  string new_seq_string = "";
  string rep_string = "";
  string duplicate_sequence = "";
  int32_t iDelStart = 0;
  int32_t iDelEnd = 0;
  int32_t iInsStart = 0;
  int32_t iInsEnd = 0;        
  int32_t iDupLen = 0;
  int32_t uRPos = -1;
  
  //@JEB - these checks are really unnecessary, since these values always exist for a MOB
  
  //Size to delete from start of repeat string.
  if(mut.entry_exists("del_start")) 
    iDelStart = from_string<int32_t>(mut["del_start"]);
  if(mut.entry_exists("del_end"))   
    iDelEnd = from_string<int32_t>(mut["del_end"]);
  ASSERT((iDelStart >= 0) && (iDelEnd >= 0), (to_string(mut._type) + " " + mut._id) + " - NEGATIVE DELETION");
  
  if(mut.entry_exists("repeat_pos"))
    uRPos = from_string<int32_t>(mut["repeat_pos"]);        
  if(mut.entry_exists("duplication_size"))
    iDupLen = from_string<int32_t>(mut["duplication_size"]);
  
  // @JEB: correct here to look for where the repeat is in the original ref_seq_info???
  // This saves us from possibly looking at a shifted location...
  rep_string = ref_seq_info.repeat_family_sequence(mut["repeat_name"], from_string<int16_t>(mut["strand"]), uRPos);
  
  mut["repeat_size"] = to_string(rep_string.length()); // saving this for shifting
  
  // This is the string we're going to pass to be inserted.
  // It will eventually contain the repeat string, insertions
  // and the duplicate_sequence.
  new_seq_string = rep_string;
  
  // Do we have deletes?  Go ahead and delete them from the repeat.
  if(iDelStart)
    new_seq_string.replace(0,iDelStart,"");
  if(iDelEnd)
    new_seq_string.resize(new_seq_string.size() - iDelEnd);
  
  // If there are any inserts, put them in front of or behind the repeat sequence.
  if(mut.entry_exists("ins_start")) {
    new_seq_string = mut["ins_start"] + rep_string;
  }
  
  if(mut.entry_exists("ins_end"))   {
    new_seq_string += mut["ins_end"];
  }
  
  return new_seq_string;
}
  
// This function will use the current GD and apply it to the new_ref_seq_info.
// When calling this function make SURE that you load ref_seq_info and
// new_ref_seq_info separately.
//
void cGenomeDiff::apply_to_sequences(cReferenceSequences& ref_seq_info, cReferenceSequences& new_ref_seq_info, bool verbose)
{    
  uint32_t count_SNP = 0, count_SUB = 0, count_INS = 0, count_DEL = 0, count_AMP = 0, count_INV = 0, count_MOB = 0, count_CON = 0, count_MASK = 0;
  
  // Handle all mutation types, plus MASK four-letter type.
  diff_entry_list_t mutation_list = this->mutation_list();
  diff_entry_list_t mask_list = this->list(make_vector<gd_entry_type>(MASK));
  
  mutation_list.insert(mutation_list.end(), mask_list.begin(), mask_list.end());
  
  for (diff_entry_list_t::iterator itr_mut = mutation_list.begin(); itr_mut != mutation_list.end(); itr_mut++)
  {
    cDiffEntry& mut(**itr_mut);
    uint32_t position = from_string<uint32_t>(mut[POSITION]);
    
    switch (mut._type) 
    {
      case SNP :
      {
        count_SNP++;
        
        new_ref_seq_info.replace_sequence_1(mut[SEQ_ID], position, position, mut[NEW_SEQ], (to_string(mut._type) + " " + mut._id), verbose);
        
        if (verbose)
        {
          cout << "SNP: 1 bp => " << mut[NEW_SEQ] << " at position " << position << endl;
          cout << "  shift: none" << endl;
        }
      } break;
        
      case SUB:
      {
        count_SUB++;
        
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);
        
        new_ref_seq_info.replace_sequence_1(mut[SEQ_ID], position, position + size - 1, mut[NEW_SEQ], (to_string(mut._type) + " " + mut._id), verbose);
        
        if (verbose)
        {
          cout << "SUB: " << size << " => " << mut[NEW_SEQ] << endl;
          cout << "   shift +" << mut[NEW_SEQ].length() << " bp at position " <<  position << endl;
        }
      } break;
        
      case INS:
      {          
        count_INS++;
        
        new_ref_seq_info.insert_sequence_1(mut[SEQ_ID], position, mut[NEW_SEQ], (to_string(mut._type) + " " + mut._id), verbose);
        
        if (verbose)
        {
          cout << "INS: 0 => " << mut[NEW_SEQ] << endl;
          cout << "   shift +" << mut[NEW_SEQ].length() << " bp at position " <<  position << endl;
        }
      } break;
        
      case DEL:
      {
        count_DEL++;
        
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);
        
        new_ref_seq_info.replace_sequence_1(mut[SEQ_ID], position, position + size -1, "", (to_string(mut._type) + " " + mut._id), verbose);
        
        if (verbose)
        {
          cout << "DEL: " << size << " => 0" << endl;
          cout << "   shift -" << size << " bp at position " <<  position << endl;
        }
      } break;
        
      case MASK:
      {
        count_MASK++;
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);
        string mask_string(size, 'N');
        new_ref_seq_info.replace_sequence_1(mut[SEQ_ID], position, position + size - 1, mask_string, (to_string(mut._type) + " " + mut._id), verbose);
        if (verbose)
        {
          cout << "MASK: " << size << " => N" << endl;
          cout << "   shift: none" << endl;
        }
        
      } break;
        
      case AMP:
      {
        count_AMP++;
        
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);
        
        //Build duplicate sequence
        string dup;
        for (uint32_t i = 1; i < from_string<uint32_t>(mut["new_copy_number"]); i++)
          dup.append(new_ref_seq_info.get_sequence_1(mut[SEQ_ID], position, position+size-1));
        ASSERT(!dup.empty(), "Duplicate sequence is empty. You may have specified an AMP with a new copy number of 1.");
        
        new_ref_seq_info.insert_sequence_1(mut[SEQ_ID], position-1, dup, (to_string(mut._type) + " " + mut._id), verbose);
        
        if (verbose)
        {
          cout << "AMP: 0" << " => " << dup << endl;
          cout << "   shift +" << dup.length() << " bp at position " << position << endl;
        }        
      } break;
        
      case INV:
      {
        count_INV++;
        
        WARN("INV: mutation type not handled yet");
      } break;
        
      case MOB:
      {
        ASSERT(mut["strand"] != "?", "Unknown MOB strand (?)\n" + mut.to_string());
        
        count_MOB++;
        string new_seq_string = "";
        string rep_string = "";
        string duplicate_sequence = "";
        int32_t iDelStart = 0;
        int32_t iDelEnd = 0;
        int32_t iInsStart = 0;
        int32_t iInsEnd = 0;        
        int32_t iDupLen = 0;
        int32_t uRPos = -1;
        
        if (verbose) cout << "APPLYING MOB" << endl;
        if (verbose) cout << mut << endl;
        
        //Size to delete from start of repeat string.
        if(mut.entry_exists("del_start")) iDelStart = from_string<int32_t>(mut["del_start"]);
        if(mut.entry_exists("del_end"))   iDelEnd = from_string<int32_t>(mut["del_end"]);
        ASSERT((iDelStart >= 0) && (iDelEnd >= 0), (to_string(mut._type) + " " + mut._id) + " - NEGATIVE DELETION");
        
        if(mut.entry_exists("repeat_pos"))uRPos = from_string<int32_t>(mut["repeat_pos"]);        
        if(mut.entry_exists("duplication_size"))iDupLen = from_string<int32_t>(mut["duplication_size"]);
        
        // @JEB: correct here to look for where the repeat is in the original ref_seq_info???
        // This saves us from possibly looking at a shifted location...
        cSequenceFeature chosen_feature; // return value
        string seq_id;                   // return value
        rep_string = ref_seq_info.repeat_family_sequence(mut["repeat_name"], from_string<int16_t>(mut["strand"]), uRPos, &seq_id, &chosen_feature);
        
        if (verbose) cout << ">Sequence of repeat chosen for " << mut["repeat_name"] << endl;
        if (verbose) cout << rep_string << endl;
        
        mut["repeat_size"] = to_string(rep_string.length()); // saving this for shifting
        
        // This is the string we're going to pass to be inserted.
        // It will eventually contain the repeat string, insertions
        // and the duplicate_sequence.
        new_seq_string = rep_string;
        
        // Do we have deletes?  Go ahead and delete them from the repeat.
        if(iDelStart)new_seq_string.replace(0,iDelStart,"");
        if(iDelEnd)new_seq_string.resize(new_seq_string.size() - iDelEnd);
        
        // The position of a MOB is the first position that is duplicated
        // Inserting at the position means we have to copy the duplication
        // in FRONT OF the repeat sequence
        if(iDupLen > 0) {
          duplicate_sequence = new_ref_seq_info.get_sequence_1(mut[SEQ_ID], position, position + (from_string<uint32_t>(mut["duplication_size"]) - 1));
        }
        else if (iDupLen < 0) {
          // A negative duplication length indicates that this many bases were deleted from the original genome starting at the specified base
          new_ref_seq_info.replace_sequence_1(mut[SEQ_ID], position, position + abs(iDupLen)-1, "");
        }
        
        // If there are any inserts, put them in front of or behind the
        // repeat sequence.
        if(mut.entry_exists("ins_start")) {new_seq_string = mut["ins_start"] + rep_string;iInsStart = mut["ins_start"].length();}
        if(mut.entry_exists("ins_end"))   {new_seq_string += mut["ins_end"];iInsEnd = mut["ins_end"].length();}
        
        new_seq_string = mob_replace_sequence(ref_seq_info, mut);
        
        // Add on the duplicated sequence.  This happens AFTER
        // we have inserted any insertions.
        new_seq_string = duplicate_sequence + new_seq_string;        
        
        if (verbose) cout << ">Sequence that will be inserted after position: " << position-1 << endl;
        if (verbose) cout << new_seq_string << endl;
        
        // Insert our newly minted sequence.
        new_ref_seq_info.insert_sequence_1(mut[SEQ_ID], position-1, new_seq_string, (to_string(mut._type) + " " + mut._id), verbose);
        
        // We've repeated the sequence, now it's time to repeat all the features
        // inside of and including the repeat region.
        new_ref_seq_info.repeat_feature_1(mut[SEQ_ID], position+iInsStart+duplicate_sequence.size(), iDelStart, iDelEnd, ref_seq_info, mut["repeat_name"], from_string<int16_t>(mut["strand"]) * chosen_feature.get_strand(), uRPos, verbose);
        
        if (verbose)
        {
          cout << "MOB: " << ( (iDupLen > 0 ) ? "0" : to_string(abs(iDupLen)) ) << " => " << new_seq_string << endl;
          cout << "   shift +" << new_seq_string.length() - iDelStart - iDelEnd << " bp at position " << position << endl;
        }          
      } break;
        
      case CON:
      {
        count_CON++;
        
        uint32_t size = from_string<uint32_t>(mut[SIZE]);
        
        uint32_t replace_target_id, replace_start, replace_end;
        new_ref_seq_info.parse_region(mut["region"], replace_target_id, replace_start, replace_end);
        ASSERT(replace_start != replace_end, "Cannot process CON mutation with end == start. ID:" + to_string(mut._id));
        
        int8_t strand = (replace_start < replace_end) ?  +1 : -1;
        
        if (strand == -1) {
          swap(replace_start, replace_end);
        }
        
        // @JEB: correct here to look for where the replacing_sequence is in the original ref_seq_info.
        // This saves us from possible looking at a shifted location...
        string replacing_sequence = ref_seq_info[replace_target_id].get_sequence_1(replace_start, replace_end);
        
        if (strand == -1) {
          replacing_sequence = reverse_complement(replacing_sequence);
        }
        
        string displaced_sequence = new_ref_seq_info.get_sequence_1(mut[SEQ_ID], position, position + size - 1);
        
        new_ref_seq_info.replace_sequence_1(mut[SEQ_ID], position, position + size - 1, replacing_sequence, (to_string(mut._type) + " " + mut._id), verbose); 
        
        if (verbose)
        {
          cout << "CON: " << displaced_sequence << " => " << replacing_sequence << endl;
          int32_t size_change = static_cast<int32_t>(replacing_sequence.length()) - static_cast<int32_t>(displaced_sequence.length());
          cout << "   shift " << ((size_change >= 0) ? "+" : "") << size_change << " bp at position " << position << endl;
        }       
      } break;
        
      default:
        ASSERT(false, "Can't handle mutation type: " + to_string(mut._type));
    }
    
    this->shift_positions(mut, new_ref_seq_info, verbose);
  }
  
  if (verbose)
  {    
    cout << "MUTATION COUNT" << endl;
    cout << "\tSNP: " << count_SNP << endl;
    cout << "\tSUB: " << count_SUB << endl;
    cout << "\tINS: " << count_INS << endl;
    cout << "\tDEL: " << count_DEL << endl;
    cout << "\tAMP: " << count_AMP << endl;
    cout << "\tINV: " << count_INV << endl;
    cout << "\tMOB: " << count_MOB << endl;
    cout << "\tCON: " << count_CON << endl;
    cout << "\tMASK: " << count_MASK << endl;
  }
  
  //Cleanup.  If any of the sequences are of zero length, remove them.
  for (vector<cAnnotatedSequence>::iterator it_as = new_ref_seq_info.begin(); it_as < new_ref_seq_info.end(); it_as++) {
    if(!it_as->m_length){new_ref_seq_info.erase(it_as);it_as--;}
  }
}
  
void cGenomeDiff::normalize_to_sequence(cReferenceSequences &ref)
{
  diff_entry_list_t muts = this->mutation_list();
  for (diff_entry_list_t::iterator it = muts.begin();
       it != muts.end(); it++) {
    cDiffEntry &mut = **it;
    mut.normalize_to_sequence(ref[mut["seq_id"]]);
  }
}

cGenomeDiff cGenomeDiff::validate(cGenomeDiff& ctrl, cGenomeDiff& test, bool verbose)
{
  bool (*comp_fn) (const diff_entry_ptr_t&, const diff_entry_ptr_t&) = diff_entry_ptr_sort;
  typedef set<diff_entry_ptr_t, bool(*)(const diff_entry_ptr_t&, const diff_entry_ptr_t&)> diff_entry_set_t;
  
  //Isolate control and test mutations into sets for quick lookup.
  diff_entry_list_t muts = ctrl.mutation_list();
  diff_entry_set_t ctrl_muts(comp_fn);
  copy(muts.begin(), muts.end(), inserter(ctrl_muts, ctrl_muts.begin()));
  
  muts = test.mutation_list();
  diff_entry_set_t test_muts(comp_fn);
  copy(muts.begin(), muts.end(), inserter(test_muts, test_muts.begin()));
  
  if (verbose) {
    printf("\tComparing %u control mutations versus %u test mutations.\n\n",
           static_cast<unsigned int>(ctrl_muts.size()), static_cast<unsigned int>(test_muts.size()));
  }
  
  /* Combine alike mutations, however we may lose information like supporting evidence IDs
   * which we will search ctrl_muts and test_muts for.
   */
  diff_entry_set_t uniq_muts(comp_fn);
  std::set_union(
                 ctrl_muts.begin(),
                 ctrl_muts.end(),
                 test_muts.begin(),
                 test_muts.end(),
                 inserter(uniq_muts, uniq_muts.begin())
                 );
  
  uint32_t n_tp = 0, n_fn = 0, n_fp = 0;
  
  cGenomeDiff ret_val;
  ret_val.metadata = test.metadata;
  
  for (diff_entry_set_t::iterator it = uniq_muts.begin(); it != uniq_muts.end(); ++it) {
    assert(ctrl_muts.count(*it) || test_muts.count(*it));
    
    string key = "";
    diff_entry_list_t evidence;
    diff_entry_set_t::iterator it_ctrl = ctrl_muts.find(*it);
    diff_entry_set_t::iterator it_test = test_muts.find(*it);
    bool in_ctrl = (it_ctrl != ctrl_muts.end());
    bool in_test = (it_test != test_muts.end());
    if (in_ctrl && in_test) {
      key = "TP";
      ++n_tp;
      evidence = test.mutation_evidence_list(**it_test);
      if (evidence.empty()) {
        evidence = ctrl.mutation_evidence_list(**it_ctrl);
      }
    }
    else if (in_ctrl && !in_test) {
      key = "FN";
      ++n_fn;
      evidence = ctrl.mutation_evidence_list(**it_ctrl);
    }
    else if (!in_ctrl && in_test) {
      key = "FP";
      ++n_fp;
      evidence = test.mutation_evidence_list(**it_test);
    } 
    
    if (verbose) {
      string temp = "";
      if (key == "TP") temp = "[True  Positive]:\t";
      if (key == "FN") temp = "[False Negative]:\t";
      if (key == "FP") temp = "[False Positive]:\t";
      if (key == "")   temp = "[ERROR]         :\t";
      cout << "\t\t"<< temp << **it << endl;
    }
    
    (**it)["compare"] = key;
    (**it)._evidence.clear();
    
    //Add supporting evidence and assign new IDs to the mutation.
    diff_entry_ptr_t mut = ret_val.add(**it);
    if (evidence.size()) {
      for (diff_entry_list_t::iterator jt = evidence.begin(); jt != evidence.end(); ++jt) {
        (**jt)["compare"] = key;
        mut->_evidence.push_back(ret_val.add(**jt)->_id);
      }
    } else {
      mut->_evidence.push_back(".");
    }
    
  }
  
  if (verbose) {
    printf("\tUpdating mutation's evidence.\n"); 
  }
  
  //Add TP|FN|FP header info.
  string value = "";
  sprintf(value, "%u|%u|%u", n_tp, n_fn, n_fp);
  ret_val.add_breseq_data("TP|FN|FP", value);
  
  printf("\t#=TP|FN|FP	%s \t for %s versus %s \n",
         value.c_str(),
         ctrl._default_filename.c_str(),
         test._default_filename.c_str());
  
  return ret_val;
}

  
// Returns true if two junction sequences are considered equivalent
//  -- subsequence or subsequence of reverse complement
bool equivalent_junction_sequences(string s1, string s2) {
  if (s1.size() < s2.size()) {
    
    if (s2.find(s1) != string::npos ) {
      return true;
    }
    if (s2.find(reverse_complement(s1)) != string::npos) {
      return true;
    }
  }
  else {
    if (s1.find(s2) != string::npos ) {
      return true;
    }
    if (s1.find(reverse_complement(s2)) != string::npos) {
      return true; 
    }
  }
  return false;
}
  
typedef map<string, diff_entry_ptr_t> jc_data_t;  
  
// Helper function that returns keys of all equivalent junctions
vector<string> equivalent_junction_keys(jc_data_t& jcs, string& key)
{
  vector<string> matching_keys;
  
  for (jc_data_t::iterator it = jcs.begin(); it != jcs.end(); it++) {
    
    if (equivalent_junction_sequences(key, it->first))
      matching_keys.push_back(it->first);
  }
  
  return matching_keys;
}

cGenomeDiff cGenomeDiff::validate_evidence(cReferenceSequences& sequence, 
                                           uint32_t buffer,
                                           uint32_t shorten_length,
                                           cGenomeDiff& ctrl,
                                           cGenomeDiff& test,
                                           bool verbose) {
  
  //TODO currently only compares JC evidence.
  cGenomeDiff ret_val;
  ret_val.metadata = test.metadata;
  
  
  // START JC Evidence Block
  //
  // Conceptually we build up two sets containing equivalent junction sequences.
  // One is from the control set and one is from the test set.
  // 
  // * The control set must have no equivalent junctions in it.
  //
  // * The test set isn't penalized for predicting an equivalent junction many times
  //   Subsequent appearances of the same junction are not forwarded.
  
  // This keeps track of all unique sequences.
  // We assume that each sequence will only occur ONCE in each genome diff file.
  // Otherwise you can have more True Positives than the total that should be possible if you
  // predict the same junction twice...
  
  jc_data_t ctrl_jc;
  
  diff_entry_list_t ctrl_list = ctrl.list(make_vector<gd_entry_type>(JC));
  ctrl_list.remove_if(cDiffEntry::field_exists("circular_chromosome")); 
  
  ////////////////////////////
  //      CONTROL list      //
  ////////////////////////////
  
  for (diff_entry_list_t::iterator it = ctrl_list.begin(); it != ctrl_list.end(); ++it) {
    cDiffEntry& jc = **it;
    
    string jc_segment = CandidateJunctions::construct_junction_sequence(sequence, jc, buffer, true);
    jc["segment"] = jc_segment;
    ASSERT(jc_segment.size(), "Could not locate JC sequence for: " + (*it)->to_string());
    
    // This properly deals with reverse-complements and subsequences
    vector<string> equivalent_segments = equivalent_junction_keys(ctrl_jc, jc_segment);
    if ( equivalent_segments.size() > 0 ) {
      ERROR("Duplicate junction sequence in control data set for entry:\n" + jc.to_string() + "\n" + ctrl_jc[jc_segment]->to_string());
    }
    
    // We have to shorten the control segments to allow for this situation to be equivalent...
    // GATCTAGTCATGCTAC
    //  ATCTAGTCATGCTACG
    
    jc_segment = jc_segment.substr(shorten_length, jc_segment.size()-2*shorten_length);
    
    ctrl_jc[jc_segment] = *it;
  }
  
  
  /////////////////////////
  //      TEST list      //
  /////////////////////////
  
  // For recovering full information about a junction from the set
  jc_data_t test_jc;
  
  diff_entry_list_t test_list = test.list(make_vector<gd_entry_type>(JC));
  test_list.remove_if(cDiffEntry::field_exists("circular_chromosome")); 
  
  int32_t i = 0;
  for (diff_entry_list_t::iterator it = test_list.begin(); it != test_list.end(); ++it) {
    
    cDiffEntry& jc = **it;
    
    string jc_segment = CandidateJunctions::construct_junction_sequence(sequence, jc, buffer, true);
    jc["segment"] = jc_segment;
    ASSERT(jc_segment.size(), "Could not locate JC sequence for: " + (*it)->to_string());
    
    vector<string> equivalent_segments = equivalent_junction_keys(test_jc, jc_segment);
    
    for (vector<string>::iterator its=equivalent_segments.begin(); its != equivalent_segments.end(); its++) {
      
      diff_entry_ptr_t prev_jc = test_jc[*its];
      
      WARN("Duplicate junction sequence in test data set for entry:\n" + jc.to_string() + "\n" + ctrl_jc[jc_segment]->to_string());
      
      if (verbose) {
        cerr << "*** Merged two junctions:" << endl;
        cerr << jc << endl;
        cerr << "AND:" << endl;
        cerr << *(prev_jc) << endl;
      }
      
      // Save the max score with the new item
      if (n(jc["score"]) < n((*prev_jc)["score"])) {
        jc["score"] = (*prev_jc)["score"];
      }
      
      if (verbose) {
        cerr << "Result:" << endl;
        cerr << jc << endl;
      }
      
      test_jc.erase(*its);
    }
    
    test_jc[jc_segment] = *it;
  }
  
  //////////////////////////////////////
  //     Assignment of TP, FP, FN     //
  //////////////////////////////////////
  
  uint32_t n_tp = 0, n_fn = 0, n_fp = 0;
  
  // Create the new item, which must either be a true-positive or a false-positive!
  
  for (jc_data_t::iterator it = test_jc.begin(); it != test_jc.end(); it++) {
    
    string test_junction_seq = it->first;
    cDiffEntry& jc = *(it->second);
    
    cDiffEntry new_item = jc.to_spec();
    new_item["segment"] = test_junction_seq;
    
    if (jc.count("score")) {
      new_item["score"] = jc["score"];
    } else if (jc.count("neg_log10_pos_hash_p_value")) {
      new_item["score"] = jc["neg_log10_pos_hash_p_value"];
    } else {
      new_item["score"] = "9999999";
    }
    
    vector<string> equivalent_segments = equivalent_junction_keys(ctrl_jc, test_junction_seq);
    
    // One item found = TP
    if (equivalent_segments.size() == 1) {    
      new_item["compare"] = "TP";
      ++n_tp;
      // delete from control set
      ctrl_jc.erase(equivalent_segments[0]);
      // Zero items found = FP
    } else if (equivalent_segments.size() == 0) {
      new_item["compare"] = "FP";
      ++n_fp;
    }
    // More than one item found = ERROR! (Should be ruled out above.)
    else {
      ERROR("More than one match in control found for junction:\n" + jc.to_string());
    }
    
    ret_val.add(new_item);
  }
  
  // Now iterate through remaining control items, which are false-negatives
  for (jc_data_t::iterator it = ctrl_jc.begin(); it != ctrl_jc.end(); it++) {
    string test_junction_seq = it->first;
    cDiffEntry& jc = *(it->second);
    
    cDiffEntry new_item = jc.to_spec();
    new_item["segment"] = test_junction_seq;
    new_item["compare"] = "FN";
    new_item["score"] = "0"; // Control items don't have a score...
    ++n_fn;
    ret_val.add(new_item);
  }
  
  //Add TP|FN|FP header info.
  string value = "";
  sprintf(value, "%u|%u|%u", n_tp, n_fn, n_fp);
  ret_val.add_breseq_data("TP|FN|FP", value);
  
  printf("\t#=TP|FN|FP	%s \t for %s versus %s \n",
         value.c_str(),
         ctrl._default_filename.c_str(),
         test._default_filename.c_str());
  
  return ret_val;
}

// @JEB: This is only partially implemented -- it adds JC for mutations involving these
void cGenomeDiff::mutations_to_evidence(cReferenceSequences &ref_seq, bool remove_mutations)
{
  diff_entry_list_t muts = mutation_list();
  this->remove(EVIDENCE);
  this->remove(VALIDATION);
  for(diff_entry_list_t::iterator it=muts.begin(); it!=muts.end(); ++it) {
    cDiffEntry& de = **it;
    if (de._type == MOB) {
      // Return the actual copy of the sequence feature found...
      cSequenceFeature repeat_feature;
      string repeat_seq_id;
      int32_t repeat_position = -1;
      ref_seq.repeat_family_sequence(de["repeat_name"], +1, repeat_position, &repeat_seq_id, &repeat_feature);
      
      cDiffEntry de1;
      de1._type = JC;
      de1["side_1_seq_id"] = de["seq_id"];
      de1["side_1_position"] = s(n(de["position"]) + n(de["duplication_size"]) - 1);
      de1["side_1_strand"] = "-1";
      de1["side_2_seq_id"] = repeat_seq_id;
      int32_t strand = repeat_feature.m_location.m_strand * n(de["strand"]);
      de1["side_2_position"] = (strand > 0) ? s(repeat_feature.m_location.m_start) : s(repeat_feature.m_location.m_end);
      de1["side_2_strand"] = s(strand);
      de1["overlap"] = "0";
      (*it)->_evidence.push_back(this->add(de1)->_id);
      
      cDiffEntry de2;
      de2._type = JC;
      de2["overlap"] = "0";
      de2["side_1_seq_id"] = de["seq_id"];
      de2["side_1_position"] = s(n(de["position"]));
      de2["side_1_strand"] = "1";
      de2["side_2_seq_id"] = repeat_seq_id;
      de2["side_2_position"] = (strand > 0) ? s(repeat_feature.m_location.m_end) : s(repeat_feature.m_location.m_start);
      de2["side_2_strand"] = s(-strand);
      (*it)->_evidence.push_back(this->add(de2)->_id);
      
    } else if (de._type == DEL) {
      
      cDiffEntry de1;
      de1._type = JC;
      de1["overlap"] = "0";
      de1["side_1_seq_id"] = de["seq_id"];
      de1["side_1_position"] = s(n(de["position"])-1);
      de1["side_1_strand"] = "-1";
      de1["side_2_seq_id"] = de["seq_id"];
      de1["side_2_position"] = s(n(de["position"])+n(de["size"]));
      de1["side_2_strand"] = "1";
      (*it)->_evidence.push_back(this->add(de1)->_id);
      
    }
  }
  if (remove_mutations) {
    this->remove(MUTATIONS);
  }
}
  

void cGenomeDiff::write_jc_score_table(cGenomeDiff& compare, string table_file_path, bool verbose) {
  //assert(compare.metadata.breseq_data.count("TP|FN|FP"));

  typedef map<float, map<string, uint32_t> > table_t;
  table_t table;
  
  diff_entry_list_t jc = compare.list(make_vector<gd_entry_type>(JC));
  double max_score = 0;
  uint32_t total_gold_standard_predictions = 0;
  for (diff_entry_list_t::iterator it = jc.begin(); it != jc.end(); ++it) {
    if (!(*it)->count("score")) {
      cerr << "No score value for: " + (*it)->to_string() << endl;
      continue;
    }
    double score = from_string<double>((**it)["score"]);
    //score = roundp<10>(score);
    max_score = max(max_score, score);

    assert((*it)->count("compare"));
    string compare = (**it)["compare"];

    if ( (compare == "TP") || (compare == "FN") ) {
      total_gold_standard_predictions++; 
    }
    if (compare == "FN") continue;
    
    if (!table.count(score)) {
      table[score]["TP"] = 0;
      table[score]["FN"] = 0;
      table[score]["FP"] = 0;
    }
  
    table[score][compare]++;
  }


  ofstream out(table_file_path.c_str());
  ASSERT(out, "Could not write to file: " + table_file_path);
  
  out << "score" << '\t' << "TP" << '\t' << "FN" << '\t' << "FP" << endl;
  if (verbose) {
    cerr << "\t\tscore" << '\t' << "TP" << '\t' << "FN" << '\t' << "FP" << endl;
  }
  uint32_t n_tp = 0, n_fn = 0, n_fp = 0;
  
  //TODO @GRC for sorting tophat scoring parameters, remove if/when needed.
  bool is_tophat_scoring = compare.metadata.author == "tophat" ? true : false;
  if (is_tophat_scoring) {
    for (table_t::reverse_iterator it = table.rbegin(); it != table.rend(); it++) {
      double i = it->first;
      //double i = 0; i <= max_score; i += .1f) {
      
      if (table.count(i)) {
        n_tp += table[i]["TP"];
        // The total number of TP and FN at a given score is equal to the total minus the number of TP up to that point.
        n_fn =  total_gold_standard_predictions - n_tp;
        n_fp += table[i]["FP"];

        out << i << '\t' << n_tp << '\t' << n_fn << '\t' << n_fp << endl;
        if (verbose) {
          cerr << "\t\t" << i << '\t' << n_tp << '\t' << n_fn << '\t' << n_fp << endl;
        }

      }

    }
  } else {
    for (table_t::iterator it = table.begin(); it != table.end(); it++) {
      double i = it->first;
      //double i = 0; i <= max_score; i += .1f) {
      
      if (table.count(i)) {
        n_tp += table[i]["TP"];
        // The total number of TP and FN at a given score is equal to the total minus the number of TP up to that point.
        n_fn =  total_gold_standard_predictions - n_tp;
        n_fp += table[i]["FP"];

        out << i << '\t' << n_tp << '\t' << n_fn << '\t' << n_fp << endl;
        if (verbose) {
          cerr << "\t\t" << i << '\t' << n_tp << '\t' << n_fn << '\t' << n_fp << endl;
        }

      }

    }
  }
  out.close();
}

  
// Convert GD file to VCF file
//
// VCF Format https://github.com/amarcket/vcf_spec
// This function was written by consulting v4.2
// (which gives the example below)
//
// ##fileformat=VCFv4.2
// ##fileDate=20090805
// ##source=myImputationProgramV3.1
// ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
// ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
// ##phasing=partial
// ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
// ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
// ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
// ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
// ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
// ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
// ##FILTER=<ID=q10,Description="Quality below 10">
// ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
// ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
// ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
// #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
//   20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
//   20 17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3 0/0:41:3
//   20 1110696 rs6040355 A G,T 67 PASS NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2 2/2:35:4
//   20 1230237 . T . 47 PASS NS=3;DP=13;AA=T GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
//   20 1234567 microsat1 GTC G,GTCT 50 PASS NS=3;DP=9;AA=G GT:GQ:DP 0/1:35:4 0/2:17:2 1/1:40:3
//
// IMPORTANT: Although it doesn't seem to specify -- ALL columns must be separated by TABS, not spaces
  
void cGenomeDiff::write_vcf(const string &vcffile, cReferenceSequences& ref_seq_info)
{  
  ofstream output( vcffile.c_str() );
  output << "##fileformat=VCFv4.1" << endl;
  output << "##fileDate" << endl;
  output << "##source=breseq_GD2VCF_converterter" << endl;
  //output << "##reference=" << endl;
  
  
  // Write contig information.
  
  for(cReferenceSequences::iterator it=ref_seq_info.begin(); it!=ref_seq_info.end(); it++) {
    output << "##contig=<ID=" << it->m_seq_id << ",length=" << it->get_sequence_length() << ">" << endl;
  }
  
  output << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">" << endl;
  output << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;

  
  // Write header line
  //output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample" << endl;
  output << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << endl;

  diff_entry_list_t muts = mutation_list();  
  for (diff_entry_list_t::iterator it=muts.begin(); it!=muts.end(); it++) {

    cDiffEntry& mut = **it;
    uint32_t pos = from_string<uint32_t>(mut[POSITION]);
    
    string chrom = mut[SEQ_ID];
    string id = ".";
    string ref;
    string alt;
    string qual = ".";
    string filter = ".";
    string info;
    string format = "GT";
    string sample = "1/1";
    
    // For now, only allele frequency info field
    double freq = 1.0;
    if (mut.count(FREQUENCY))freq = from_string<double>(mut[FREQUENCY]);
    string AF = formatted_double(freq, 4).to_string(); // allele frequency
    info = "AF=" + AF;
    
    switch (mut._type) 
    {
      case SNP:
      { 
        alt = mut[NEW_SEQ];
        ref = ref_seq_info.get_sequence_1(mut[SEQ_ID], pos, pos);
        
        // Carry forward qulity from related RA evidence
        diff_entry_list_t ev = mutation_evidence_list(mut);
        if (ev.size() == 1) 
          qual = ev.front()->count("genotype_quality") ? (*ev.front())["genotype_quality"] : ".";

      } break;
        
      case SUB:
      {
        // We know Ref seq is not just a "." in this case.
        
        alt = mut[NEW_SEQ];
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);
        ref = ref_seq_info.get_sequence_1(mut[SEQ_ID], pos, pos + size - 1);

        // Carry forward qulity from related RA evidence
        // @JEB should do something with multiple evidence
        diff_entry_list_t ev = mutation_evidence_list(mut);
        if (ev.size() == 1) 
          qual = ev.front()->count("genotype_quality") ? (*ev.front())["genotype_quality"] : ".";
        
      } break;
        
      case INS:
      {          
        // Ref has to be something - so take base insertion was after
        // Complication: It could be inserted at the beginning of the sequence.
        // In this case, we take the base after.
        
        bool before_base = (pos != 0);
        ref = before_base ? ref_seq_info.get_sequence_1(mut[SEQ_ID], pos, pos) : ref_seq_info.get_sequence_1(mut[SEQ_ID], pos+1, pos+1);

        // Correct position of first base shown, if necessary
        if (before_base) pos--;
        
        alt = before_base ? ref + mut[NEW_SEQ] : mut[NEW_SEQ] + ref;
        
        // Carry forward qulity from related RA evidence
        // @JEB should do something with multiple evidence
        diff_entry_list_t ev = mutation_evidence_list(mut);
        if (ev.size() == 1) 
          qual = ev.front()->count("genotype_quality") ? (*ev.front())["genotype_quality"] : ".";
        
      } break;
        
      case DEL:
      {
        // Complication: Alt has to be something (not .)
        // So, take first base before deletion or after if it is the first base
        
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);
        ref = ref_seq_info.get_sequence_1(mut[SEQ_ID], pos, pos + size - 1);
        
        bool before_base = (pos != 1);
        alt = before_base ? ref_seq_info.get_sequence_1(mut[SEQ_ID], pos-1, pos-1) : ref_seq_info.get_sequence_1(mut[SEQ_ID], pos, pos);
        ref = before_base ? alt + ref : ref + alt;

        // Correct position of first base shown, if necessary
        if (before_base) pos--;
        
        // Carry forward qulity from related RA evidence
        // @JEB should do something with multiple evidence
        diff_entry_list_t ev = mutation_evidence_list(mut);
        if (ev.size() == 1) 
          qual = ev.front()->count("genotype_quality") ? (*ev.front())["genotype_quality"] : ".";
        
      } break;
        
      case AMP:
      {        
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);
        
        //Build duplicate sequence
        string alt;
        for (uint32_t i = 0; i < from_string<uint32_t>(mut["new_copy_number"]); i++)
          alt.append(ref_seq_info.get_sequence_1(mut[SEQ_ID], pos, pos+size-1));
        ASSERT(!alt.empty(), "Duplicate sequence is empty. You may have specified an AMP with a new copy number of 1.");
        
        ref = ref_seq_info.get_sequence_1(mut[SEQ_ID], pos, pos+size-1);

        // @JEB should do something with multiple evidence
        diff_entry_list_t ev = mutation_evidence_list(mut);
        if (ev.size() == 1) qual = ev.front()->count("genotype_quality") ? (*ev.front())["genotype_quality"] : ".";
        
      } break;
        
      case INV:
      {        
        WARN("INV: mutation type not handled yet");
      } break;
        
      case MOB:
      {
        // Ref has to be something - so take first base before deletion
        // Or after if there is none before
        
        bool before_base = (pos != 1);
        ref = before_base ? ref_seq_info.get_sequence_1(mut[SEQ_ID], pos-1, pos) : ref_seq_info.get_sequence_1(mut[SEQ_ID], pos, pos + 1);
        //ref = before_base ? alt + ref : ref + alt;
        
        // This includes the IS and all adjacent duplicated or deleted nucleotides  
        string new_seq_string = mob_replace_sequence(ref_seq_info, mut);
        
        // The position of a MOB is the first position that is duplicated
        // Inserting at the position means we have to copy the duplication
        // in FRONT OF the repeat sequence
        
        string duplicate_sequence;
        int32_t iDupLen = from_string<int32_t>(mut["duplication_size"]);
        if (iDupLen < 0) {
          cerr << "Warning: MOB with negative target site insertion not handled. Ignoring:" << endl << mut << endl;
          break;
        }
        if (iDupLen > 0) {
          duplicate_sequence = ref_seq_info.get_sequence_1(mut[SEQ_ID], pos, pos + iDupLen - 1);
        }
        
        // Add on the duplicated sequence.  This happens AFTER
        // we have inserted any insertions.
        new_seq_string = duplicate_sequence + new_seq_string;        
        
        if (before_base) 
          alt = ref + duplicate_sequence + new_seq_string;
        else
          alt = duplicate_sequence + new_seq_string + ref;
        
        // Correct position of first base shown, if necessary
        if (before_base) pos--;
        
      } break;
        
      case CON:
      {        
        uint32_t size = from_string<uint32_t>(mut[SIZE]);
        
        uint32_t replace_target_id, replace_start, replace_end;
        ref_seq_info.parse_region(mut["region"], replace_target_id, replace_start, replace_end);
        ASSERT(replace_start != replace_end, "Cannot process CON mutation with end == start. ID:" + to_string(mut._id));
        
        int8_t strand = (replace_start < replace_end) ?  +1 : -1;
        
        if (strand == -1)
          swap(replace_start, replace_end);
        
        // @JEB: correct here to look for where the replacing_sequence is in the original ref_seq_info.
        // This saves us from possible looking at a shifted location...
        alt = ref_seq_info[replace_target_id].get_sequence_1(replace_start, replace_end);
        
        if (strand == -1)
          alt = reverse_complement(alt);
        
        ref = ref_seq_info.get_sequence_1(mut[SEQ_ID], pos, pos + size - 1);
        
      } break;
        
      default:
        WARN("Can't handle mutation type: " + to_string(mut._type) + "\nIt will be skipped.");
    }
    
    //output << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << qual << "\t" << filter << "\t" << info << "\t" << format << "\t" << sample << endl;
    
        output << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << qual << "\t" << filter << "\t" << info << endl;
  }
  
  output.close();
}
  
// Convert GD file to GVF file
void cGenomeDiff::write_gvf(const string &gvffile, cReferenceSequences& ref_seq_info, bool snv_only)
{  

  diff_entry_list_t diff_entry_list = this->list();
  diff_entry_list_t::iterator it = diff_entry_list.begin();
  
  // Stores the features
  vector< vector<string> > features;
  vector< vector<string> > featuresGVF;
  
  // We only write entries for mutations
  diff_entry_list_t mut_list = this->mutation_list();
  
  // Processes the features
  // gvf[0]: ID of reference
  // gvf[1]: Source
  // gvf[2]: Type
  // gvf[3]: Start
  // gvf[4]: End
  // gvf[5]: Score
  // gvf[6]: Strand
  // gvf[7]: Phase
  // gvf[8]: Attributes    
  
  for (diff_entry_list_t::iterator it=mut_list.begin(); it != mut_list.end(); it++)
  {
    cDiffEntry& de = **it;
    
    vector<string> gvf(9,"");
    
    for( int j=5; j<8; j++ ){
      gvf[j] = ".";
    }
    
    // Common to all entries
    // SeqID
    gvf[0] = de[SEQ_ID];
    // Source
    gvf[1] = "breseq";
    // Position
    gvf[3] = de[POSITION];
    
    if( de._type == SNP )
    {
      // Type
      gvf[2] = "SNV";
      // End
      gvf[4] = de[POSITION];
      // Strand
      gvf[6] = "+";
      // Attributes - Reference base
      gvf[8].append(";Reference_seq=").append( ref_seq_info.get_sequence_1(de[SEQ_ID], from_string(de[POSITION]),  from_string(de[POSITION])) );
      // Attributes - New base
      gvf[8].append("Variant_seq=").append( de[NEW_BASE] );
      
      diff_entry_list_t ev_list = mutation_evidence_list(de);
      ASSERT(ev_list.size() == 1, "Did not find RA evidence supporting SNP\n" + to_string(de))
      cDiffEntry& ev = *(ev_list.front());
      
      // Score
      gvf[5] = ev[QUALITY];
        
      // Attributes - Total Reads 
      vector<string> covs = split( ev[TOT_COV], "/" );
      uint32_t cov = from_string<uint32_t>(covs[0]) + from_string<uint32_t>(covs[1]);
      gvf[8] = gvf[8].append(";Total_reads=").append(to_string(cov));
      
      // Attributes - Variant Reads 
      vector<string> variant_covs = split( ev[NEW_COV], "/" );
      uint32_t variant_cov = from_string<uint32_t>(variant_covs[0]) + from_string<uint32_t>(variant_covs[1]);
      gvf[8] = gvf[8].append(";Variant_reads=").append(to_string(variant_cov));
        
      // Attributes - Frequency 
      gvf[8].append(";Variant_freq=").append( ev[FREQUENCY] );
      

      if (de.entry_exists("snp_type")) {
        if (de["snp_type"] == "nonsynonymous") {
          gvf[8].append(";Variant_effect=non_synonymous_codon");
        }
        else if (de["snp_type"] == "synonymous") {
          gvf[8].append(";Variant_effect=synonymous_codon");
        }
        else if (de["snp_type"] == "intergenic") {
          gvf[8].append(";Variant_effect=intergenic_variant");
        }
        else if (de["snp_type"] == "RNA") {
          gvf[8].append(";Variant_effect=nc_transcript_variant");
        }
        else if (de["snp_type"] == "pseudogene") {
          gvf[8].append(";Variant_effect=nc_transcript_variant");
        }
      }
    } // END of SNP
    
    else if( de._type == SUB ){
      gvf[2] = "indel";
      // End
      gvf[4] = to_string(from_string(de[POSITION]) + from_string(de[SIZE])); 
      gvf[8].append("Reference_seq=").append( ref_seq_info.get_sequence_1(de[SEQ_ID], from_string(de[POSITION]), from_string(de[POSITION]) + from_string(de[SIZE]) - 1));
      gvf[8].append(";Variant_seq=").append( de[NEW_SEQ] );
    }
    
    else if( de._type == DEL ){
      gvf[2] = "deletion";
      gvf[4] = gvf[3];
      gvf[8].append("Reference_seq=").append( ref_seq_info.get_sequence_1(de[SEQ_ID], from_string(de[POSITION]), from_string(de[POSITION]) + from_string(de[SIZE]) - 1) );
      gvf[8].append(";Variant_seq=").append( "." );
    }
    
    else if( de._type == INS ){
      gvf[2] = "insertion";
      gvf[4] = gvf[3];
      gvf[8].append("Reference_seq=").append( "." );
      gvf[8].append(";Variant_seq=").append( de[NEW_SEQ] );
    }
    
    else if( de._type == MOB ){
      gvf[2] = "mobile_element_insertion";
      gvf[4] = gvf[3];
      //Strand
      if( from_string(de["strand"]) > 0 )
        gvf[6] = "+";
      else
        gvf[6] = "-";
      gvf[8].append("Reference_seq=").append( "." );
      gvf[8].append(";Variant_seq=").append( mob_replace_sequence(ref_seq_info, de) );
    }
    
    else if( de._type == AMP )
    {
       gvf[2] = "copy_number_gain";
       stringstream ss;
       gvf[4] = gvf[3];
       gvf[8].append("Reference_seq=").append( "." );
       gvf[8].append(";Variant_seq=").append( ref_seq_info.get_sequence_1(de[SEQ_ID], from_string(de[POSITION]), from_string(de[POSITION]) + from_string(de[SIZE]) - 1) );
    }
    else if( de._type == INV ){
      gvf[2] = "inversion";
      gvf[4] = to_string(from_string(de[POSITION]) + from_string(de[SIZE]) - 1);
    }
    else if( de._type == CON ){
      gvf[2] = "substitution";
      gvf[4] = gvf[3];
      
      uint32_t tid, start_pos, end_pos;
      ref_seq_info.parse_region(de["region"], tid, start_pos, end_pos);
      
      gvf[8].append("Reference_seq=").append( ref_seq_info.get_sequence_1(de[SEQ_ID], from_string(de[POSITION]), from_string(de[POSITION]) + from_string(de[SIZE]) - 1) );
      gvf[8].append(";Variant_seq=").append( ref_seq_info.get_sequence_1(tid, start_pos, end_pos ));
    }
    
    // ID attribute
    if( gvf[8].compare( "" ) == 0 || ( gvf[8].size()>8 && !gvf[8].substr(0,3).compare("ID=") == 0) ){
      string s = "";
      s.append("ID=").append(gvf[0]).append(":").append(gvf[1]).append(":");
      s.append(gvf[2]).append(":").append(gvf[3]).append(";");
      s.append(gvf[8]);
      gvf[8] = s;
    }
    
    if (!snv_only || (de._type == SNP))
      featuresGVF.push_back(gvf);
  }
  
  // Write results to file
  ofstream output( gvffile.c_str() );
  output << "##gff-version 3" << endl;
  output << "##gvf-version 1.0" << endl;
  output << "" << endl;
  output << "##source-method Source=breseq;Type=SNV;Dbxref=http://barricklab.org/breseq;Comment=Mapping and variant calling with breseq;" << endl;
  output << "" << endl;
  for( size_t i=0; i<featuresGVF.size(); i++ ){
    for( size_t j=0; j<featuresGVF[i].size(); j++ ){
      output << featuresGVF[i][j] << "\t";
    }
    output << "\n";
  }
  output.close();
  
}

  
void cGenomeDiff::read_vcf(const string &file_name)
{
  //VCF Column order.
  enum {CHROM = 0, POS, ID, REF, ALT, QUAL, FILTER, INFO};
  
  ifstream in(file_name.c_str());
  string line = "";
  
  while (!in.eof()) {
    getline(in, line);
    if (!in.good() || line.empty()) break;
    //Discard header lines for now.
    if (line[0] == '#') continue;
    const vector<string> &tokens = split(line, "\t");
    
    cDiffEntry de;
    const cString &rseq = tokens[REF], &aseq = tokens[ALT];
    if (aseq.contains(',')) {
      //! TODO: Classify these mutations, GATK doesn't predict these(?)
      //        de._type = SUB;
      //        de[SEQ_ID]   = tokens[CHROM];
      //        de[POSITION] = tokens[POS];
      //        de[SIZE]     = to_string(tokens[ALT].size());
      //        de[NEW_SEQ]  = tokens[ALT];
      WARN("Can't classify line: " + line);
      continue;
    }
    if (rseq.size() > aseq.size()) {
      de._type = DEL;
      de[SEQ_ID]   = tokens[CHROM];
      de[POSITION] = to_string(from_string<uint32_t>(tokens[POS]) + aseq.size());
      de[SIZE]     = to_string(rseq.size() - aseq.size());
    }
    else if (rseq.size() < aseq.size()) {
      de._type = INS;
      de[SEQ_ID]   = tokens[CHROM];
      de[POSITION] = tokens[POS];
      de[NEW_SEQ]  = cString(aseq).remove_starting(rseq);
    }
    else if (rseq.size() == 1) {
      de._type = SNP;
      de[SEQ_ID]   = tokens[CHROM];
      de[POSITION] = tokens[POS];
      de[NEW_SEQ]  = aseq;
    }
    else {
      WARN("Can't classify line: " + line);
      continue;
    }
    
    const vector<string> &info_tokens = split(tokens[INFO], ";");
    //size_t n = info_tokens.size();
    for (size_t i = 0, n = info_tokens.size(); i < n; ++i) {
      cKeyValuePair kvp(info_tokens[i], '=');
      de[kvp.get_key()] = kvp.get_value();
    }
    
    add(de);
  }
  
}
  
    
void write_phylip( const vector<string>& gd_file_names, 
               const vector<string>& reference_file_names,
               const string& output_phylip_file_name)
{
  
  (void) gd_file_names;
  (void) reference_file_names;
  (void) output_phylip_file_name;
  
  ERROR("Not implemented");
  /*
   
   my @mutations = $composite_gd->mutation_list;
   
   if ($ignore_ambiguous)
   {
   MUT: for (my $i=0; $i<scalar @mutations; $i++)
   {
   my $mut = $mutations[$i];
   foreach my $gd_name (@$gd_name_list_ref)
   {
   my $f = $mut->{"frequency_$gd_name"};
   if ($f eq '?')
   {
   $mut->{"frequency_$gd_name"} = ($mut->{type} eq 'SNP') ? $mut->{_ref_seq} : 'A';
   }
   }
   } 
   }
   
   open KEYFILE, ">$phylip_input_file.key" or die "Could not open file: $phylip_input_file.key";
   MUT: for (my $i=0; $i<scalar @mutations; $i++)
   {
   my $mut = $mutations[$i];
   print KEYFILE +join("\t", $i, $mut->{type}, $mut->{position}) . "\n"; 
   }
   close KEYFILE;
   
   open PHYLIP, ">$phylip_input_file" or die "Could not open file: $phylip_input_file";
   #first line is number of species, then number of sites
   print PHYLIP +(scalar @$gd_name_list_ref) . " " . +(scalar @mutations) . "\n";
   
   foreach my $gd_name (@$gd_name_list_ref)
   {			
   my $mut_str;
   foreach my $mut (@mutations)
   {
   my $f = $mut->{"frequency_$gd_name"};
   my $b = '?';
   if ($f eq '?')
   {
   if ($phylip_for_mix)
   {
   $b = '?'; #ancestor
   }
   else
   {
   $b = 'N';
   }
   }
   elsif (($f == 1) || ($f eq 'H'))
   {
   if ($phylip_for_mix)
   {
   $b = '1';
   }
   else
   {
   $b = ($mut->{type} eq 'SNP') ? $mut->{new_seq} : 'T';
   }
   }
   elsif ($f == 0)
   {
   if ($phylip_for_mix)
   {
   $b = '0';
   }
   else
   {
   $b = ($mut->{type} eq 'SNP') ? $mut->{_ref_seq} : 'A';
   }
   }
   else
   {
   die "Unknown frequency: $f\n";
   }
   $mut_str .= $b;				
   }
   
   printf PHYLIP "%-10.10s%s\n", $gd_name, $mut_str;
   }
   
   close PHYLIP;
   */
}
 
// Unlike other conversion functions, takes a list of gd files
void cGenomeDiff::GD2Circos(const vector<string> &gd_file_names, 
                            const vector<string> &reference_file_names,
                            const string &circos_directory,
                            double distance_scale,
                            double feature_scale){
  
  cGenomeDiff combined_gd;
  
  int32_t number_of_mutations = 0;
  
  for (size_t i = 0; i < gd_file_names.size(); i++){
    cGenomeDiff single_gd(gd_file_names[i]);
    combined_gd.merge(single_gd, false);
    number_of_mutations += single_gd.mutation_list().size();
  }
  
  cReferenceSequences ref;
  ref.LoadFiles(reference_file_names);
  ref.annotate_mutations(combined_gd, true);
  const vector<string> seq_ids(ref.seq_ids());
  
  string make_me;
  
  Settings settings("");
  
  
  create_path(circos_directory);
  create_path(circos_directory + "/data");
  create_path(circos_directory + "/etc");
  
  //copy run script
  copy_file(settings.program_data_path + "/run_circos.sh", circos_directory + "/run_circos.sh");
  
  //filling circos_dir/etc
  
  vector<string> conf_names = make_vector<string>
  ("ideogram.conf")
  ("karyotype.and.layout.conf")
  ("indels.conf")
  ("mobs.conf")
  ("mutations.conf")
  ("combined_circos.conf")
  ;
  
  for (size_t i = 0; i < conf_names.size(); i++){
    copy_file(settings.program_data_path + "/" + conf_names[i], circos_directory + "/etc/" + conf_names[i]);
  }
  
  //modifying circos_dir/etc with scale values
  stringstream command;
  double distance_value = 0.4 * distance_scale;
  double space_value = 0.25 * distance_value;
  double feature_value = 5 * feature_scale;
  
  map<string, string> replacement_map = make_map<string, string>
  (" = inner_distance_value_1", " = " + to_string(distance_value, 10) + "r")
  (" = inner_distance_value_2", " = " + to_string(distance_value + .01, 10) + "r")
  (" = inner_ticks", " = " + to_string(distance_value, 10) + "r")
  (" = syn_axis_value_1", " = " + to_string(distance_value + (2 * space_value), 10) + "r")
  (" = syn_axis_value_2", " = " + to_string(distance_value + (2 * space_value) + .01, 10) + "r")
  (" = nonsyn_axis_value_1", " = " + to_string(distance_value + (3 * space_value), 10) + "r")
  (" = nonsyn_axis_value_2", " = " + to_string(distance_value + (3 * space_value) + .01, 10) + "r")
  (" = npi_axis_value_1", " = " + to_string(distance_value + (4 * space_value), 10) + "r")
  (" = npi_axis_value_2", " = " + to_string(distance_value + (4 * space_value) + .01, 10) + "r")
  (" = mob_1_axis_value_1", " = " + to_string(distance_value + (5 * space_value), 10) + "r")
  (" = mob_1_axis_value_2", " = " + to_string(distance_value + (5 * space_value) + .01, 10) + "r")
  (" = mob_2_axis_value_1", " = " + to_string(distance_value + (6 * space_value), 10) + "r")
  (" = mob_2_axis_value_2", " = " + to_string(distance_value + (6 * space_value) + .01, 10) + "r")
  (" = outer_axis_value_1", " = " + to_string(distance_value + (7 * space_value), 10) + "r")
  (" = outer_axis_value_2", " = " + to_string(distance_value + (7 * space_value) + .01, 10) + "r")
  (" = outer_ticks", " = " + to_string(distance_value + (7 * space_value), 10) + "r")
  (" = indel_distance", " = " + to_string(distance_value + (7 * space_value) + .01, 10) + "r")
  (" = mob_distance", " = " + to_string(distance_value + (5.5 * space_value), 10) + "r - " + to_string(feature_value * 2.5) + "p")
  //(" = mob_distance", " = " + to_string(distance_value + (5 * space_value), 10) + "r + " + to_string((400 * .7 * space_value * .5) - (feature_value * 2.5), 10) + "p")
  (" = syn_distance", " = " + to_string(distance_value + (2 * space_value), 10) + "r - " + to_string(feature_value * 1.5, 10) + "p")
  (" = nonsyn_distance", " = " + to_string(distance_value + (3 * space_value), 10) + "r - " + to_string(feature_value * 1.5, 10) + "p")
  (" = npi_distance", " = " + to_string(distance_value + (4 * space_value), 10) + "r - " + to_string(feature_value * 1.5, 10) + "p")
  (" = indel_value", " = " + to_string(feature_value * 4, 10) + "p")
  (" = mob_value", " = " + to_string(feature_value * 5, 10) + "p")
  (" = snp_value", " = " + to_string(feature_value * 3, 10) + "p")
  (" = ind_syn_axis_value_1", " = " + to_string(1 + (1 * space_value), 10) + "r")
  (" = ind_syn_axis_value_2", " = " + to_string(1 + (1 * space_value) + .01, 10) + "r")
  (" = ind_nonsyn_axis_value_1", " = " + to_string(1 + (2 * space_value), 10) + "r")
  (" = ind_nonsyn_axis_value_2", " = " + to_string(1 + (2 * space_value) + .01, 10) + "r")
  (" = ind_npi_axis_value_1", " = " + to_string(1 + (3 * space_value), 10) + "r")
  (" = ind_npi_axis_value_2", " = " + to_string(1 + (3 * space_value) + .01, 10) + "r")
  (" = ind_syn_distance", " = " + to_string(1 + (1 * space_value), 10) + "r - " + to_string(feature_value * 1.5, 10) + "p")
  (" = ind_nonsyn_distance", " = " + to_string(1 + (2 * space_value), 10) + "r - " + to_string(feature_value * 1.5, 10) + "p")
  (" = ind_npi_distance", " = " + to_string(1 + (3 * space_value), 10) + "r - " + to_string(feature_value * 1.5, 10) + "p")
  (" = space_value_in_pixels", " = " + to_string((300 * space_value), 10) + "p")
  (" = ind_scale", " = " + to_string((.7 * distance_scale), 10) + "r")
  (" = label_size_value", " = " + to_string((16 * distance_scale), 10) + "p")
  (" = label_offset_value", " = " + to_string((-7.5 * distance_scale), 10) + "r")
  ;
  
  for (size_t i = 0; i < conf_names.size(); i++){
    replace_file_contents_using_map(circos_directory + "/etc/" + conf_names[i],
                                    circos_directory + "/etc/circos_temp.conf",
                                    replacement_map);
    copy_file(circos_directory + "/etc/circos_temp.conf",
              circos_directory + "/etc/" + conf_names[i]);
  }
  
  //filling circos_dir/data
  
  ofstream karyotype_file;
  ofstream empty_file;
  
  make_me = circos_directory + "/data/karyotype.txt";
  karyotype_file.open(make_me.c_str());
  make_me = circos_directory + "/data/empty_data.txt";
  empty_file.open(make_me.c_str());
  
  //keeps track of current position when examining sequence sizes of genomes
  int32_t current_position = 0;
  
  int32_t half_ref_length;
  half_ref_length = int32_t(ref.total_length() / 2) + 1;
  
  for (uint32_t i = 0; i < seq_ids.size(); i++){
    uint32_t seq_size;
    seq_size = ref[seq_ids[i]].get_sequence_size();
    
    karyotype_file << "chr - " << seq_ids[i] << " 1 1 " <<
    seq_size << " black" << endl;
    empty_file << seq_ids[i] << " 1 2 1" << endl;
    
    //if seq_size goes past halfway point of total length of genomes,
    //add this sequence and its bounds to the left_side vector.
    current_position += seq_size;
  }
  
  karyotype_file.close();
  empty_file.close();
  
  //minimum tile size width for indel graph
  const int32_t MIN_WIDTH = static_cast<int32_t>(floor(static_cast<double>(ref.total_length()) * 0.000));
  const int32_t MIN_DISPLAY_LENGTH = 51; //inclusive
  
  ofstream indel_file;
  ofstream mob_file;
  
  ofstream synonymous_mutation_file;
  ofstream nonsynonymous_mutation_file;
  ofstream npi_mutation_file;
  
  make_me = circos_directory + "/data/indels_data.txt";
  indel_file.open(make_me.c_str());
  make_me = circos_directory + "/data/mobs_data.txt";
  mob_file.open(make_me.c_str());
  
  make_me = circos_directory + "/data/syn_data.txt";
  synonymous_mutation_file.open(make_me.c_str());
  make_me = circos_directory + "/data/nonsyn_data.txt";
  nonsynonymous_mutation_file.open(make_me.c_str());
  make_me = circos_directory + "/data/npi_data.txt";
  npi_mutation_file.open(make_me.c_str());
  
  map <string, string> mob_colors;
  
  //colors for mobs
  const char* c_colors[] = {"vvdred", "vvdgreen", "vvdblue", "vvdorange", "vvdpurple",
    "dred", "dgreen", "dblue", "dorange", "dpurple",
    "red", "green", "blue", "orange", "purple",
    "lred", "lgreen", "lblue", "lorange", "lpurple",
    "vvlred", "vvlgreen", "vvlblue",  "vvlorange", "vvlpurple"};
  
  map<string,bool> pre_assigned_mob_colors = make_map<string,bool>
  ("IS1"  , true )("IS186", true )("IS3"  , true  ) 
  ("IS150", true )("IS911", true )("IS4"  , true  )
  ("IS2"  , true )("IS30" , true )("IS600", true  )
  ;
  
  vector<string> colors(c_colors, c_colors + 25);
  string color;
  int32_t next_color = 0;
  
  //reference sequence MOBs
  for(size_t i = 0; i < ref.size(); i++){
    cAnnotatedSequence& ref_seq = ref[i];
    
    for(cSequenceFeatureList::iterator it = ref_seq.m_repeats.begin(); it != ref_seq.m_repeats.end(); it++){
      cSequenceFeature& seq_feature = **it;
      int32_t middle = int32_t(seq_feature.m_location.get_start_1() + seq_feature.m_location.get_end_1()) / 2;
      
      string color;
      
      // Color assignment -- prefer preassigned, then grab next from list
      // and assign that color permanently to copies of this repeat
      if (pre_assigned_mob_colors.count(seq_feature["name"])){
        color = seq_feature["name"];
        mob_colors[seq_feature["name"]] = color;
      }
      else if (mob_colors.count(seq_feature["name"]) == 0){
        color = colors[next_color];
        mob_colors[seq_feature["name"]] = color;
        next_color++;
      }
      else{
        color = mob_colors[seq_feature["name"]];
      }
      
      mob_file << ref_seq.m_seq_id << " " <<
      middle << " " <<
      middle << " " <<
      "i" << ((seq_feature.m_location.m_strand == 1)? "right" : "left" ) << " " <<
      "color=" << color << endl;
    }
  }
  
  
  
  diff_entry_list_t gd_data = combined_gd.mutation_list();
  
  for (diff_entry_list_t::iterator it = gd_data.begin(); it != gd_data.end(); it++){
    
    cDiffEntry diff = **it;
    
    int32_t width;
    string direction;
    
    
    if (diff._type == INS){
      width = diff["new_seq"].size();
      if (width < MIN_DISPLAY_LENGTH){
        continue;
      }
      if (width < MIN_WIDTH){
        width = MIN_WIDTH;
      }
      
      indel_file << diff["seq_id"] << " " <<
      diff["position"] << " " <<
      from_string<int32_t>(diff["position"]) + width << " " <<
      "color=green" << endl;
    }
    else if (diff._type == AMP){
      width = from_string<int32_t>(diff["size"]);
      if (width < MIN_DISPLAY_LENGTH){
        continue;
      }
      if (width < MIN_WIDTH){
        width = MIN_WIDTH;
      }
      
      indel_file << diff["seq_id"] << " " <<
      diff["position"] << " " <<
      from_string<int32_t>(diff["position"]) + width << " " <<
      "color=green" << endl;
    }
    else if (diff._type == DEL){
      width = from_string<int32_t>(diff["size"]);
      if (width < MIN_DISPLAY_LENGTH){
        continue;
      }
      if (width < MIN_WIDTH){
        width = MIN_WIDTH;
      }
      indel_file << diff["seq_id"] << " " <<
      diff["position"] << " " <<
      from_string<int32_t>(diff["position"]) + width << " " <<
      "color=red" << endl;
      
      // Show new IS for IS mediated deletions
      if (diff.count("mediated")) {
        
        if (mob_colors.count(diff["mediated"])) { 
          // either the end or the beginning is in an IS element
          
          int32_t max_distance_to_repeat_1 = 0;
          int32_t max_distance_to_repeat_2 = 0;
          cSequenceFeaturePtr feat1 = cReferenceSequences::find_closest_repeat_region_boundary(n(diff["position"]) - 1, ref[diff["seq_id"]].m_repeats, max_distance_to_repeat_1,-1);
          cSequenceFeaturePtr feat2 = cReferenceSequences::find_closest_repeat_region_boundary(n(diff["position"]) + n(diff["size"]) + 1 - 1, ref[diff["seq_id"]].m_repeats, max_distance_to_repeat_2,1);
          
          if (!feat1.get() && !feat2.get()) {
            cerr << diff << endl;
            ASSERT(false,"Could not find mediating repeat.");
          }  
          
          cSequenceFeature& seq_feature = feat1.get() ? *(feat1.get()) : *(feat2.get());
          
          cAnnotatedSequence ref_seq = ref[diff["seq_id"]];
          
          int32_t middle = feat1.get() ? n(diff["position"]) + n(diff["size"]) - 1 : n(diff["position"]);
          
          string color;
          
          // Color assignment -- prefer preassigned, then grab next from list
          // and assign that color permanently to copies of this repeat
          if (pre_assigned_mob_colors.count(seq_feature["name"])){
            color = seq_feature["name"];
          }
          else if (mob_colors.count(seq_feature["name"]) == 0){
            color = colors[next_color];
            mob_colors[seq_feature["name"]] = color;
            next_color++;
          }
          else{
            color = mob_colors[seq_feature["name"]];
          }
          
          mob_file << ref_seq.m_seq_id << " " <<
          middle << " " <<
          middle << " " <<
          "o" << ((seq_feature.m_location.m_strand == 1)? "right" : "left" ) << " " <<
          "color=" << color << endl;
        }
      }
    }
    else if(diff._type == SNP || diff._type == SUB){
      if (diff["snp_type"] == "synonymous"){
        synonymous_mutation_file << diff["seq_id"] << " " <<
        diff["position"] << " " <<
        diff["position"] << endl;
      }
      else if (diff["snp_type"] == "nonsynonymous"){
        nonsynonymous_mutation_file << diff["seq_id"] << " " <<
        diff["position"] << " " <<
        diff["position"] << endl;
      }
      else{
        npi_mutation_file << diff["seq_id"] << " " <<
        diff["position"] << " " <<
        diff["position"] << endl;
      }
      
    }
    
    else if(diff._type == MOB){
      
      // Color assignment -- prefer preassigned, then grab next from list
      // and assign that color permanently to copies of this repeat
      if (pre_assigned_mob_colors.count(diff["repeat_name"])){
        color = diff["repeat_name"];
      }
      else if (mob_colors.count(diff["repeat_name"]) == 0){
        color = colors[next_color];
        mob_colors[diff["repeat_name"]] = color;
        next_color++;
      }
      else{
        color = mob_colors[diff["repeat_name"]];
      }
      
      mob_file << diff["seq_id"] << " " <<
      diff["position"] << " " <<
      diff["position"] << " " <<
      "o" << ((n(diff["strand"]) == 1)? "right" : "left" ) << " " <<
      "color=" << color << endl;
    }
  }
  
  
  indel_file.close();
  mob_file.close();
  synonymous_mutation_file.close();
  nonsynonymous_mutation_file.close();
  npi_mutation_file.close();
  
  char current_dir[1024];
  ASSERT(getcwd(current_dir, sizeof(current_dir)), "Linux function call getcwd() has failed");
  command.str("");
  command << "cd " << circos_directory << "; bash run_circos.sh;";// cd " << current_dir;
  SYSTEM(command.str());
  
}  
  
//@JEB deprecated 08-16-2013
/*
 void MIRA2GD(const string &mira_file_name, const string &gd_file_name){
 //this was made on accident :(
 ifstream mira_file;
 
 mira_file.open(mira_file_name.c_str());
 
 if (!mira_file){
 cerr << "Could not open " << mira_file_name << endl;
 return;
 }
 
 string line;
 
 cGenomeDiff gd;
 
 while (getline(mira_file, line)){
 vector <string> items = split_on_whitespace(line);
 
 if (items[0] == "#"){
 continue;
 }
 
 
 cDiffEntry mut(UNKNOWN);
 
 if (items[9] == "basechange"){
 mut._type = SNP;
 mut["seq_id"] = items[0].substr(0, items[0].length() - 3);
 mut["position"] = items[4];
 mut["new_seq"] = items[10].substr(items[10].length() - 1);
 }
 else if (items[9] == "insertion"){
 mut._type = INS;
 mut["seq_id"] = items[0].substr(0, items[0].length() - 3);
 mut["position"] = items[4].substr(items[4].find(":") + 1);
 mut["new_seq"] = items[10].substr(items[10].length() - 1);
 }
 else if (items[9] == "deletion"){
 mut._type = DEL;
 mut["seq_id"] = items[0].substr(0, items[0].length() - 3);
 mut["position"] = items[4];
 mut["size"] = "1";
 }
 if (mut._type != UNKNOWN) {
 gd.add(mut);
 }
 }
 gd.write(gd_file_name);
 mira_file.close();
 }
 */
  
cFlaggedRegions& cFlaggedRegions::read_mummer(string file_path, cAnnotatedSequence& ref_seq) {
  /*
   * mummer -maxmatch -b -c -l 36 REL606.fna REL606.fna > exclude
   * 
   * Input file format:
   > SeqID
   1544819    3595894       266
   3398753    4156036      200
   2712737    2713011        38
   ......
   > SeqID Reverse
   Start1     Start2    Length
   1544819    3595894       266
   3398753    4156036       200
   2712737    2713011        38
   *
   * Columns are: Start1     Start2    Length
   * In the first block, the match goes forward from both coords
   * In the second block, the match goes up forward the first coord and downward from the second coord
   */
  
  ifstream in(file_path.c_str());
  assert(in);
  
  //Remove header.
  in.ignore(1000, '\n');
  
  bool on_reverse_block = false;
  
  string line = "";
  while (getline(in, line)) {
    const vector<string>& tokens = split_on_whitespace(line);
    uint32_t first  = from_string<uint32_t>(tokens[0]);
    uint32_t second = from_string<uint32_t>(tokens[1]);
    uint32_t size   = from_string<uint32_t>(tokens[2]);
    
    if ((tokens[0] == ">") && (tokens[2] == "Reverse")) {
      on_reverse_block = true;
      //cerr << "Found reverse block" << endl;
      continue;
    }
    
    if (!on_reverse_block) {
      this->flag_region(first,  first  + size - 1);
      this->flag_region(second, second + size - 1);
      string seq1 = ref_seq.get_sequence_1(first,  first  + size - 1);
      string seq2 = ref_seq.get_sequence_1(second, second + size - 1);
      ASSERT(seq1 == seq2, "Problem with line in MUMmer output file. Not repeat:\n" + line);
    } else {
      this->flag_region(first,  first  + size - 1);
      this->flag_region(second - size + 1, second);
      string seq1 = ref_seq.get_sequence_1(first,  first + size - 1);
      string seq2 = ref_seq.get_sequence_1(second - size + 1, second);
      seq2 = reverse_complement(seq2);
      ASSERT(seq1 == seq2, "Problem with line in MUMmer output file. Not repeat:\n" + line);
    }
    
  }
  
  //this->print();
  
  return *this;
}

cFlaggedRegions& cFlaggedRegions::read_nucmer_tab_coords(string file_path) {
  /*
   
   @JEB 08-12-2013 this method does not find all exact matches.... DO NOT USE!!
   
   options.addUsage("An exclusion file should be generated using these MUMmer commands:");
   options.addUsage("  nucmer --maxmatch -p reference reference.fasta reference.fasta");
   options.addUsage("  delta-filter -i 100 -l 50 reference.delta > reference.filtered.delta");
   options.addUsage("  show-coords -T reference.filtered.delta > reference.exclude.coords");
   options.addUsage("The resulting reference.exclude.coords file can be used with the --exclude option.");
   
   * Input file format:
   /Users/jbarrick/tmp/sv/REL606.fna /Users/jbarrick/tmp/sv/REL606.fna
   NUCMER
   [S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[FRM]	[TAGS]
   1	4629812	1	4629812	4629812	4629812	100.00	1	1	REL606	REL606
   15386	16730	430844	429500	1345	1345	100.00	1	-1	REL606	REL606
   
   *
   * Generated by MUMmer: nucmer --maxmatch of a sequence on itself, 
   * see help for gdtools simulate-mutations
   */
  
  ifstream in(file_path.c_str());
  assert(in);
  
  //Remove header.
  in.ignore(1000, '\n');
  in.ignore(1000, '\n');
  in.ignore(1000, '\n');
  in.ignore(1000, '\n');
  in.ignore(1000, '\n');  // this is the first coord line which is sequence  to itself
  
  
  string line = "";
  while (getline(in, line)) {
    const vector<string>& tokens = split_on_whitespace(line);
    uint32_t first;  
    uint32_t second;
    first = from_string<uint32_t>(tokens[0]);
    second = from_string<uint32_t>(tokens[1]);
    if (first > second) swap(first, second);
    this->flag_region(first, second);
    
    first = from_string<uint32_t>(tokens[2]);
    second = from_string<uint32_t>(tokens[3]);
    if (first > second) swap(first, second);
    this->flag_region(first, second);
  }
  
  return *this;
}

void cFlaggedRegions::write(string file_path) {
  ofstream out(file_path.c_str());
  assert(out);
  
  out << "Start1" << '\t' << "End1" << endl;
  
  for (regions_t::iterator it = m_regions.begin(); it != m_regions.end(); ++it) {
    out << it->first << '\t' << it->second << endl;
  }
  out.close();
  
  return;
}

void cFlaggedRegions::print() {
  for (regions_t::iterator it = m_regions.begin(); it != m_regions.end(); ++it) {
    cout << it->first << '\t' << it->second << endl;
  }
  return;
}

bool cFlaggedRegions::overlaps(uint32_t pos_1, region_t region) {
  return (region.first <= pos_1 && pos_1 <= region.second); 
}

bool cFlaggedRegions::overlaps(region_t region_1, region_t region_2) {
  region_t min = region_1 < region_2 ? region_1 : region_2;
  region_t max = region_1 > region_2 ? region_1 : region_2;
  
  return ((min.second >= max.first) && (max.second >= min.first)); 
}

cFlaggedRegions& cFlaggedRegions::flag_region(uint32_t start_1, uint32_t end_1) {
  end_1 = end_1 == 0 ? start_1 : end_1;
  assert(start_1 <= end_1);
  
  regions_t regions = this->regions(start_1, end_1);
  
  //Remove and merge overlapping regions if needed.
  if (regions.size()) { 
    this->remove(regions);
    
    region_t front = *regions.begin();
    start_1 = min(start_1, front.first);
    
    region_t back = *regions.end();
    end_1 = max(end_1, back.second);
    
  }
  
  m_regions.insert(make_pair(start_1, end_1 + 1));
  
  return *this;
}

/* @JEB: Not tested
 cFlaggedRegions& cFlaggedRegions::unflag_region(uint32_t start_1, uint32_t end_1) {
 end_1 = end_1 == 0 ? start_1 : end_1;
 ASSERT(start_1 <= end_1, "[start_1]: " + s(start_1) + " is greater than [end_1]: " +s(end_1));
 
 regions_t regions = this->regions(start_1, end_1);
 
 this->remove(regions);
 
 //Add front segment if a region was partially overlapping.
 region_t front = *regions.begin();
 if (front.first < start_1) {
 this->flag_region(front.first, start_1 - 1);
 }
 
 //Add back segment if a region was partially overlapping.
 region_t back = *regions.rbegin();
 if (end_1 < back.second) {
 this->flag_region(end_1 + 1, back.second );
 }
 
 return *this;
 }
 */


bool cFlaggedRegions::is_flagged(uint32_t start_1, uint32_t end_1) {
  end_1 = end_1 == 0 ? start_1 : end_1;
  ASSERT(start_1 <= end_1, "[start_1]: " + s(start_1) + " is greater than [end_1]: " +s(end_1));
  
  return this->regions(start_1, end_1).size();
}

cFlaggedRegions& cFlaggedRegions::remove(regions_t regions) {
  for (regions_t::iterator it = regions.begin(); it != regions.end(); ++it) {
    m_regions.erase(*it);
  }
  
  return *this;
}

cFlaggedRegions::regions_t cFlaggedRegions::regions(uint32_t start_1, uint32_t end_1) {
  end_1 = (end_1 == 0) ? start_1 : end_1;
  ASSERT(start_1 <= end_1, "[start_1]: " + s(start_1) + " is greater than [end_1]: " +s(end_1));
  
  if ((start_1 == 0) && (end_1 == 0)) {
    return m_regions;
  }
  
  //Check that the lower bounds does no overlap the lower region.
  regions_t::iterator it_lower = m_regions.upper_bound(make_pair(start_1, 0));
  if (it_lower != m_regions.begin() && !this->overlaps(start_1, *--it_lower)) {
    ++it_lower;
  }
  
  regions_t::iterator it_upper = m_regions.upper_bound(make_pair(end_1 + 1, 0));
  
  return regions_t(it_lower, it_upper);
}



}//namespace bresesq


