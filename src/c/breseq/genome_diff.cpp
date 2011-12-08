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

#include "libbreseq/genome_diff.h"
#include "libbreseq/annotated_sequence.h"
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
const char* REF_BASE="ref_base";
const char* NEW_BASE="new_base";
const char* FREQUENCY="frequency";
const char* REJECT="reject";
const char* REF_COV="ref_cov";
const char* NEW_COV="new_cov";
const char* TOT_COV="tot_cov";
const char* ERROR="error";


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
(UN,make_vector<string> ("seq_id")("start")("end"))

//## validation
(CURA,make_vector<string> ("expert"))
(FPOS,make_vector<string> ("expert"))
(PHYL,make_vector<string> ("gd"))
(TSEQ,make_vector<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
(PFLP,make_vector<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
(RFLP,make_vector<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
(PFGE,make_vector<string> ("seq_id")("enzyme"))
  
; // end line specifications


const vector<string>gd_entry_type_lookup_table =
make_vector<string>("TYPE_UNKOWN")("SNP")("SUB")("DEL")("INS")("MOB")("AMP")("INV")("CON")("RA")("MC")("JC")("UN")("CURA")("FPOS")("PHYL")("TSEQ")("PFLP")("RFLP")("PFGE");


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
cDiffEntry::cDiffEntry(const gd_entry_type type)
: _type(type)
, _id("")
, _evidence()
{
}

cDiffEntry::cDiffEntry()
: _type(TYPE_UNKOWN)
, _id("")
, _evidence()
{
}

cDiffEntry::cDiffEntry(const string &line)
  : _type(TYPE_UNKOWN)
  ,_id("")
  ,_evidence()
{
  cDiffEntry *de = this;

  stringstream ss_line(line);
  typedef vector<string> string_vector_t;

  //! Step: Get the type. (Convert from string to enumeration.)
  string type = "";
  getline(ss_line, type, '\t');

  //Find the string value in the lookup table and cast index to enumeration.
  const size_t lookup_table_size = gd_entry_type_lookup_table.size();
  for (size_t i = 0; i < lookup_table_size; i++) {
    if (type == gd_entry_type_lookup_table[i]) {
      de->_type = (gd_entry_type) i;
      break;
    }
  }

  /*Check that the type is set. If there is an error confirm that gd_entry_type
   and the lookup table are identical in size and order and contain this type.*/
  if (de->_type == TYPE_UNKOWN) {
    string message = "";
    sprintf(message,
            "cDiffEntry::cDiffEntry(%s): Could not determine type.",
            line.c_str()
            );
    ERROR(message);
  }

  //! Step: Get the id.
  getline(ss_line, de->_id, '\t');

  /*! Step: If this diff entry is a mutation we need to get it's supporting
  evidence. */
  if (de->is_mutation()) {
    string evidence = "";
    getline(ss_line, evidence, '\t');

    string_vector_t evidence_vector = split(evidence, ",");

    const size_t evidence_vector_size = evidence_vector.size();
    de->_evidence.resize(evidence_vector_size);

    for (size_t i = 0; i < evidence_vector_size; i++) {
      de->_evidence[i] = evidence_vector[i];
    }

  } else {
    /*Evidence and validation diff entry types do not have evidence parameters
    so we discard this segment from the stream and continue on. */
    ss_line.ignore(254,'\t');
  }

  /*! Step: Get line specifications and assign them. */
  const string_vector_t &diff_entry_specs = line_specification[de->_type];
  const size_t diff_entry_specs_size = diff_entry_specs.size();
  for (size_t i = 0; i < diff_entry_specs_size; i++) {
     const diff_entry_key_t &key = diff_entry_specs[i];
     diff_entry_value_t value = "";
     getline(ss_line, value, '\t');
     de->insert(pair<diff_entry_key_t, diff_entry_value_t>(key, value));
  }

  /*! Step: Get the rest of the fields.*/
  while (ss_line.good()) {
    string key_value_pair = "";
    getline(ss_line, key_value_pair, '\t');
    if (key_value_pair.empty()) {
      //Not enough of a reason to exit, could be from user input.
      string message = "";
      sprintf(message,
              "cDiffEntry::cDiffEntry(%s): Empty field found.",
              line.c_str()
              );
      WARN(message);
      break;
    }

    size_t equal_sign_pos = 0;
    const size_t key_value_pair_size = key_value_pair.size();
    for (size_t i = 0; i < key_value_pair_size; i++) {
      if (key_value_pair[i] == '=') {
        equal_sign_pos = i;
        break;
      }
    }

    //Confirm this is a key=value field.
    if (!equal_sign_pos) {
      string message = "";
      sprintf(message,
              "cDiffEntry::cDiffEntry(%s): Field %s is not a key=value pair.",
              line.c_str(), key_value_pair.c_str()
              );
      ERROR(message);
    }

    const string &key   = key_value_pair.substr(0, equal_sign_pos);
    const string &value = key_value_pair.substr(equal_sign_pos + 1,
                                                key_value_pair_size);

    de->insert(pair<string, string>(key, value));
  }
  //All done!
  return;

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


/*! Marshal this diff entry into an ordered list of fields.
 */
void cDiffEntry::marshal(field_list_t& s) const {
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
  
vector<string> cDiffEntry::get_reject_reasons()
{
  vector<string> return_value;
  if (this->entry_exists("reject")) {
    return_value = split((*this)["reject"], ",");
  } 
  return return_value;
}

size_t cDiffEntry::number_reject_reasons()
{
  if(this->entry_exists(REJECT))
  {
    return this->get_reject_reasons().size();
  }
  return 0;
}
  
/*! Add reject reason to diff entry.
 */
void add_reject_reason(cDiffEntry& de, const string &reason) {

  if (de.find(REJECT) == de.end()) {
      de[REJECT] = reason;
  }
  // exists already, make comma separated list
  else {
    string reject = de[REJECT];
    reject += ",";
    reject +=reason; 
  }

}

uint32_t number_reject_reasons(cDiffEntry& de) {

	if (!de.entry_exists(REJECT) || de[REJECT].size() == 0)
		return 0;

	uint32_t reason_count = 1;
	for (uint32_t i = 0; i < de[REJECT].size(); i++)
		if (de[REJECT][i] == ',') reason_count++;
	return reason_count;
}

string cDiffEntry::to_string(void) const
{
  field_list_t fields;
  marshal(fields);

  return join(fields, "\t");
}

/*! Output operator for a diff entry.
 */
ostream& operator<<(ostream& out, const cDiffEntry& de) {
	field_list_t fields;
	de.marshal(fields);
	out << join(fields, "\t");
	return out;
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
void cGenomeDiff::add(const cDiffEntry& item, bool lowest_unique) {

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
}

//! Subtract mutations using gd_ref as reference.
void cGenomeDiff::subtract(cGenomeDiff& gd_ref, bool verbose)
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
  
cGenomeDiff cGenomeDiff::fast_merge(const cGenomeDiff& gd1, const cGenomeDiff& gd2)
{
  cGenomeDiff gd = gd1;
  
  diff_entry_list_t gd_list = gd2.list();
  for(diff_entry_list_t::const_iterator it=gd_list.begin(); it!= gd_list.end(); it++) {
    gd.add(**it);
  }
  return gd;
}

/*! Read a genome diff(.gd) from the given file to class member
  _entry_list
 */

void cGenomeDiff::read(const string& filename) {
  ifstream in(filename.c_str());
  ASSERT(in.good(), "Could not open file for reading: " + filename);

  //! Step: Handle header parameters.
  //Example header:
  //#=GENOME_DIFF 1.0
  //#=AUTHOR    Barrick JE
  //#=REFSEQ    Genbank:NC_012967.1
  //#=READSEQ   SRA:SRR066595

  /*#=GENOME_DIFF version must be initialized for this file to be recognized
   as a genome diff file. */
  metadata.version = "";
  while(in.peek() == 35) { //35 is ASCII of '#'.
    string key = "";
    getline(in, key, ' ');

    //Discard whitespace between key and value.
    while (in.peek() == 32) { //32 is ASCII of ' '
      in.ignore();
    }

    //Evaluate key.
    if (key == "#=GENOME_DIFF") {
      getline(in, metadata.version);
    }
    else if (key == "#=AUTHOR") {
      getline(in, metadata.author);
    }
    else if (key == "#=REFSEQ") {
      metadata.ref_seqs.resize(metadata.ref_seqs.size() + 1);
      getline(in, metadata.ref_seqs.back());
    }
    else if (key == "#=READSEQ") {
      metadata.read_seqs.resize(metadata.read_seqs.size() + 1);
      getline(in, metadata.read_seqs.back());
    }
    else if (key.substr(0, 2) == "#=") {
      key.erase(0,2);
      getline(in, metadata.breseq_data[key]);
    }
    else {
      //Warn if unkown header lines are encountered.
      string value = "";
      getline(in, value);
      string message = "";
      sprintf(message,
              "cGenomeDiff:read(%s): Header line: %s %s is not recognized.",
              filename.c_str(), key.c_str(), value.c_str()
             );
      WARN(message);
    }
  }

  /*Error if #=GENOME_DIFF is not found. Genome diff files are required to have
   this header line. */
  if (metadata.version.empty()) {
    string message = "";
    sprintf(message,
            "cGenomeDiff:read(%s): No #=GENOME_DIFF XX header line in this file.",
            filename.c_str()
            );
    ERROR(message);
  }

  //! Step: Handle the diff entries.
  while (in.good()) {
    string line = "";
    getline(in, line);
    //Warn if commented out or a possibly blank line is encountered.
    if (line.empty() || line[0] == '#' || line[0] == ' ') {
      size_t line_size = line.size();
      //Handle blank lines.
      bool line_is_blank = true;
      for (size_t i = 0; i < line_size; i++) {
        if (line[i] != ' ') {
          line_is_blank = false;
        }
      }
      if (line_is_blank) {
        continue;
      }
      //Handle commented lines.
      string message = "";
      sprintf(message,
              "cGenomeDiff:read(%s): Discarding line : %s",
              filename.c_str(), line.c_str()
             );
      WARN(message);
    } else {
      add(cDiffEntry(line));
    }
  }

  //All done!
  return;
}


map<gd_entry_type, sort_fields_item> diff_entry_sort_fields = make_map<gd_entry_type, sort_fields_item>
  (SNP, sort_fields_item(1, SEQ_ID, POSITION))
  (SUB, sort_fields_item(1, SEQ_ID, POSITION))
  (DEL, sort_fields_item(1, SEQ_ID, POSITION))
  (INS, sort_fields_item(1, SEQ_ID, POSITION))
  (MOB, sort_fields_item(1, SEQ_ID, POSITION))
  (AMP, sort_fields_item(1, SEQ_ID, POSITION))
  (INV, sort_fields_item(1, SEQ_ID, POSITION))
  (CON, sort_fields_item(1, SEQ_ID, POSITION))
  (RA,  sort_fields_item(2, SEQ_ID, POSITION))
  (MC,  sort_fields_item(2, SEQ_ID, START))
  (JC,  sort_fields_item(2, "side_1_seq_id", "side_1_position"))
  (UN,  sort_fields_item(3, SEQ_ID, START))
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
  (UN,  12)
;


/*! Write this genome diff to a file.
 */
bool diff_entry_ptr_sort(const diff_entry_ptr_t& a, const diff_entry_ptr_t& b) {

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

bool diff_entry_sort(const cDiffEntry &_a, const cDiffEntry &_b) {
  cDiffEntry a =_a;
  cDiffEntry b = _b;

  gd_entry_type a_type = a._type;
  gd_entry_type b_type = b._type;

  sort_fields_item a_sort_fields = diff_entry_sort_fields[a_type];
  sort_fields_item b_sort_fields = diff_entry_sort_fields[b_type];


  if (a_sort_fields._f1 < b_sort_fields._f1) {
    return true;
  } else if (a_sort_fields._f1 > b_sort_fields._f1) {
    return false;
  }

  string a_sort_field_2 = a[a_sort_fields._f2];
  string b_sort_field_2 = b[b_sort_fields._f2];

  if (a_sort_field_2 < b_sort_field_2) {
    return true;
  } else if (a_sort_field_2 > b_sort_field_2) {
    return false;
  }

  uint32_t a_sort_field_3 = from_string<uint32_t>(a[a_sort_fields._f3]);
  uint32_t b_sort_field_3 = from_string<uint32_t>(b[b_sort_fields._f3]);

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
  uint32_t a_sort_id = from_string(a._id);
  uint32_t b_sort_id = from_string(b._id);

  if (a_sort_id < b_sort_id) {
    return true;
  } else if (a_sort_id > b_sort_id) {
    return false;
  }

  return false;
}

/*! Write this genome diff to a file.
 */
void cGenomeDiff::write(const string& filename) {
  ofstream os(filename.c_str());

  //! Step: Header lines.
  /*Always write version tag. It's how we identify it as a genome diff file
   in cGenomeDiff::read(). */
  fprintf(os, "#=GENOME_DIFF %s\n", metadata.version.c_str());

  /*If adding additionaly header lines is requested (usually for the final
    genome diff file in breseq's pipeline) then add them here. */
  if (metadata.write_header) {
    fprintf(os, "#=AUTHOR    %s\n", metadata.author.c_str());

    for (size_t i = 0; i < metadata.ref_seqs.size(); i++) {
      fprintf(os, "#=REFSEQ    %s\n", metadata.ref_seqs[i].c_str());
    }

    for (size_t i = 0; i < metadata.read_seqs.size(); i++) {
      fprintf(os, "#=READSEQ   %s\n", metadata.read_seqs[i].c_str());
    }
    if (!metadata.breseq_data.empty()) {
      for (map<key_t,string>::iterator it = metadata.breseq_data.begin();
           it != metadata.breseq_data.end(); it ++) {
        fprintf(os, "#=%s %s\n", it->first.c_str(), it->second.c_str());
      }
    }
  }
  
  // sort
  _entry_list.sort(diff_entry_ptr_sort);
  
  for(diff_entry_list_t::iterator i=_entry_list.begin(); i!=_entry_list.end(); ++i) {
    fprintf(os, "%s\n", (**i).to_string().c_str());
	}
  os.close();
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


void cGenomeDiff::add_breseq_data(const key_t &key, const string& value)
{
  this->metadata.breseq_data.insert(pair<string,string>(key, value));
}


// Convert GD file to GVF file
void GDtoGVF( const string &gdfile, const string &gvffile, bool snv_only){
  
  cGenomeDiff gd_object(gdfile);
  diff_entry_list_t diff_entry_list = gd_object.list();
  diff_entry_list_t::iterator it = diff_entry_list.begin();
  
  // Stores the features
  vector< vector<string> > features;
  vector< vector<string> > featuresGVF;
  // Keeps track of the index of the entry associated with a particular evidence number
  map< string, int > eDict;

  //Read input into array
  ifstream gd( gdfile.c_str() );
  string line;
  getline( gd, line);

  while( !gd.eof() )
  {
    // split line on tabs
    char * cstr = new char [line.size()+1];
    strcpy (cstr, line.c_str());

    if( cstr[0] == '#' ){ getline(gd,line); continue; }
    vector<string> feature;
    char * pch;
    pch = strtok(cstr,"\t");

    while (pch != NULL)
    {
        feature.push_back(pch);
        pch = strtok (NULL, "\t");
        
    }  
    features.push_back(feature);

    // If it is evidence, note its index
    if( feature[0].size() == 2 ){
        eDict[ feature[1]] = (int)features.size()-1;
    }

    delete[] cstr;
    getline(gd,line);
      
  }
  gd.close();
    
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

  for( size_t i=0; i<features.size() ; i++ )
  {
    if (it == diff_entry_list.end()) break;
    cDiffEntry& de = **it;
    it++;
    
    vector<string> gvf(9,"");

    for( int j=5; j<8; j++ ){
      gvf[j] = ".";
    }

    if( features[i].size() <= 4 || features[i][0].size() == 2 ){
      continue;
    }

    // SeqID
    gvf[0] = features[i][3];
    // Source
    gvf[1] = "breseq";
    // Type
    gvf[3] = features[i][4];

    if( features[i][0].compare( "SNP") == 0 )
    {
      gvf[2] = "SNV";
      stringstream ss;
      ss << atoi( gvf[3].c_str() );
      ss >> gvf[4];
      gvf[8].append("Variant_seq=").append( features[i][5] );
      
      //Look for evidence information
      vector<string> evidenceNums = split( features[i][2], "," );
      vector<string> evidence = features[ eDict[ evidenceNums[0] ] ];
      
      gvf[6] = "+";
      
      if (evidence.size() > 0) {
        gvf[8].append(";Reference_seq=").append( evidence[5] );
      }
      
      for( size_t j=0; j<evidence.size(); j++ ){
        string s = evidence[j];
        if( s.size()>8 && s.substr(0,8).compare("quality=") == 0){
          gvf[5] = s.substr(8,s.size()-8);
        }
        if( s.size()>8 && s.substr(0,8).compare("new_cov=") == 0){
          s = s.substr(8,s.size()-8);
          vector<string> covs = split( s, "/" );
          uint32_t cov = from_string<uint32_t>(covs[0]) + from_string<uint32_t>(covs[1]);
          gvf[8] = gvf[8].append(";Variant_reads=").append(to_string(cov));
        }
      
        if( s.size()>8 && s.substr(0,8).compare("tot_cov=") == 0){
          s = s.substr(8,s.size()-8);
          vector<string> covs = split( s, "/" );
          uint32_t cov = from_string<uint32_t>(covs[0]) + from_string<uint32_t>(covs[1]);
          gvf[8] = gvf[8].append(";Total_reads=").append(to_string(cov));
        }
      }
      
      for( size_t j=0; j<features[i].size(); j++ ){
        if( features[i][j].size()>10 && features[i][j].substr(0,10).compare("frequency=") == 0){
          gvf[8].append(";Variant_freq=").append( features[i][j].substr(10,features[i][j].size()-10 ) );
        }
      }
      
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
          ERROR("Mutation annotated as pseudogene encountered. Not clear how to annotate.");
        }
      }
    }

    else if( features[i][0].compare("SUB") == 0 ){
      //Look for evidence information
      vector<string> evidenceNums = split( features[i][2], "," );
      string s = "";
      for( size_t j=0; j<evidenceNums.size(); j++ ){
        vector<string> e = features[ eDict[ evidenceNums[j] ] ];
        s.append( e[5] );
      }
      gvf[8].append("Reference_seq=").append(s);
      gvf[2] = "substitution";
      stringstream ss;
      ss << atoi( features[i][4].c_str() ) + atoi( features[i][5].c_str() );
      ss >> gvf[4]; 
      gvf[8].append(";Variant_seq=").append( features[i][6] );
    }

    else if( features[i][0].compare("DEL") == 0){
      gvf[2] = "deletion";
      stringstream ss; int length = atoi(features[i][4].c_str()) + atoi(features[i][5].c_str());
      ss << length;
      ss >> gvf[4];
    }

    else if( features[i][0].compare("INS") == 0 ){
      gvf[2] = "insertion";
      gvf[4] = gvf[3];
      gvf[8] = gvf[8].append("Variant_seq=").append( features[i][5] );
    }

    else if( features[i][0].compare("MOB") == 0 ){
      gvf[2] = "transposable_element";
      gvf[4] = gvf[3];
      gvf[8] = string("ID=").append( features[i][5] );
      //Strand
      if( atoi( features[i][6].c_str() ) > 0 ){
          gvf[6] = "+";
      }
      else{ gvf[6] = "-"; }
    }

    else if( features[i][0].compare("AMP") == 0 )
    {
      WARN("Input file contains AMP features that are currently not converted to GVF");
      /* @JEB TODO This code seg faults. Broken 2011-11-17.
      int x, y; 
      gvf[2] = "insertion";
      stringstream ss;
      ss << gvf[3] << gvf[5]; ss >> x; ss >> y; x += y; ss << x;
      ss >> gvf[4];
      gvf[3] = gvf[4];
      gvf[8].append( "Variant_seq=" );
      for( int i=0; i < (atoi(features[i][6].c_str())-1)*atoi(features[i][5].c_str()); i++ ){
        gvf[8].append("N");
      }
      */
    }
    else if( features[i][0].compare("INV") == 0 ){
      gvf[2] = "inversion";
      gvf[3] = features[i][4];
      stringstream ss;
      ss << atoi( features[i][4].c_str() ) + atoi( features[i][5].c_str() );
      ss >> gvf[4]; 
    }
    else if( features[i][0].compare("CON") == 0 ){
      gvf[2] = "substitution";
      gvf[4] = gvf[3];
      stringstream ss;
      ss << atoi( features[i][4].c_str() ) + atoi( features[i][5].c_str() );
      ss >> gvf[4];
      gvf[8].append(";Variant_seq=").append( features[i][6] );
    }

    // ID attribute
    if( gvf[8].compare( "" ) == 0 || ( gvf[8].size()>8 && !gvf[8].substr(0,3).compare("ID=") == 0) ){
      string s = "";
      s.append("ID=").append(gvf[0]).append(":").append(gvf[1]).append(":");
      s.append(gvf[2]).append(":").append(gvf[3]).append(";");
      s.append(gvf[8]);
      gvf[8] = s;
    }

    if (!snv_only || (gvf[2] == "SNV"))
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

void VCFtoGD( const string& vcffile, const string& gdfile ){
    // Stores the features
    vector< vector<string> > featuresVCF;
    vector< vector<string> > featuresGD;
    vector< vector<string> > evidences;
    // Keeps track of the index of the entry associated with a particular evidence number
    
    //Read input into array
    ifstream vcf( vcffile.c_str() );
    string line;
    getline( vcf, line);
    
    while( !vcf.eof() ){
        // split line on tabs
        char * cstr = new char [line.size()+1];
        strcpy (cstr, line.c_str());
        
        if( cstr[0] == '#' ){ getline(vcf,line); continue; }
        vector<string> feature;
        char * pch;
        pch = strtok(cstr,"\t");
        
        while (pch != NULL)
        {
            feature.push_back(pch);
            pch = strtok (NULL, "\t");
            
        }  
        featuresVCF.push_back(feature);
        
        delete[] cstr;
        getline(vcf,line);
        
    }
    vcf.close();
    
    for( size_t i=0; i<featuresVCF.size(); i++ ){
        vector<string> gd(9,"");
        vector<string> ev(9,"");
        
        // SeqID
        if( featuresVCF[i][3].size() > featuresVCF[i][4].size() ){
            gd[0] = "DEL";
        }
        else if( featuresVCF[i][3].size() < featuresVCF[i][4].size() ){
            gd[0] = "INS";
        }
        else if( featuresVCF[i][3].size() == 1 ){
            gd[0] = "SNP";
        }
        else{
            gd[0] = "SUB";
        }
        
        stringstream ss; ss << i+1;
        ss >> gd[1];
        ss << i+featuresVCF.size()+1;
        ss >> gd[2];
        gd[3] = featuresVCF[i][0];
        gd[4] = featuresVCF[i][1];
        gd[5] = featuresVCF[i][4];
        
        ev[0] = "RA";
        ss << i+featuresVCF.size();
        ss >> ev[1];
        ev[2] = ".";
        ev[3] = featuresVCF[i][0];
        ev[4] = featuresVCF[i][1];
        ev[5] = "0";
        ev[6] = featuresVCF[i][3];
        ev[7] = featuresVCF[i][4];
        ev[8] = string("quality=").append( featuresVCF[i][5] );
        
        featuresGD.push_back(gd);
        evidences.push_back(ev);
        
    }
    
    // Write results to file
    ofstream output( gdfile.c_str() );
    for( size_t i=0; i<featuresGD.size(); i++ ){
        for( size_t j=0; j<featuresGD[i].size(); j++ ){
            output << featuresGD[i][j] << "\t";
        }
        output << "\n";
    }
    output.close();
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


// return items with types that are 3 characters long
diff_entry_list_t cGenomeDiff::mutation_list()
{
  diff_entry_list_t mut_list;

   for(diff_entry_list_t::iterator itr = _entry_list.begin();
       itr != _entry_list.end(); itr ++) {
     cDiffEntry& item = **itr;
     if(item.is_mutation()) {
       mut_list.push_back(*itr);
     }
   }

  return mut_list;
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

diff_entry_ptr_t cGenomeDiff::parent(const cDiffEntry& item)
{
  for(diff_entry_list_t::iterator itr_test_item = _entry_list.begin();
      itr_test_item != _entry_list.end(); itr_test_item ++) { 
    cDiffEntry& test_item = **itr_test_item;
    for(vector<string>::iterator itr = test_item._evidence.begin();
        itr != test_item._evidence.end(); itr ++) { 
      string& test_evidence_id = (*itr);
      if(test_evidence_id == item._id)      
        return diff_entry_ptr_t(*itr_test_item);
    }
  }
  return diff_entry_ptr_t(NULL);
}

bool cGenomeDiff::mutation_unknown(cDiffEntry mut)
{
  
  if (mut._type == SNP) {
    return interval_un(from_string<uint32_t>(mut[POSITION]),
                       from_string<uint32_t>(mut[POSITION]));
  }

  //Should be updated to new unknown that includes linkage
  if (mut._type == INS) {
    return interval_un(from_string<uint32_t>(mut[POSITION]),
                       from_string<uint32_t>(mut[POSITION])) + 1;
  }
  if (mut._type == DEL) {

//#doesn't work b/c evidence list may not be correct here
//		## only call unknowns if all support is RA
//#		my $only_ra_evidence = 1;
//#		foreach my $ev ($self->mutation_evidence_list($mut))
//#		{
//#			print Dumper($ev);
//#			$only_ra_evidence &&= $ev->{type} eq 'RA';
//#		}
//#		print Dumper($mut);
//#		print "Only RA evidence? $only_ra_evidence\n";		
//#		return 0 if (!$only_ra_evidence);
		return interval_un(from_string<uint32_t>(mut[POSITION]),
                       from_string<uint32_t>(mut[POSITION]) + from_string<uint32_t>(mut["size"]) - 1);
	}
	

  if (mut._type == SUB) {
    return interval_un(from_string<uint32_t>(mut[POSITION]),
                       from_string<uint32_t>(mut[POSITION]) - 1);
  }
  
  return false;
}

void cGenomeDiff::add_reject_reasons(cDiffEntry item, const string& reject)
{
  if (item.entry_exists(REJECT))
    item[REJECT] += ",";
  else 
    item[REJECT] += reject;
}

bool 
cGenomeDiff::interval_un(const uint32_t& start,const uint32_t& end)
{
  diff_entry_list_t un_list = list(make_vector<gd_entry_type>(UN));

  for (diff_entry_list_t::iterator itr = un_list.begin();
       itr != un_list.end(); itr++) {
    cDiffEntry un(**itr);

    if (start >= from_string<uint32_t>(un[START]) &&
        end <= from_string<uint32_t>(un[END])) {
      return true;
    }
  }
  return false;
}

cReferenceSequences cGenomeDiff::apply_to_sequences(cReferenceSequences& ref_seq_info, bool verbose)
{
  // copy the reference sequence info
  cReferenceSequences new_ref_seq_info(ref_seq_info);
    
  uint32_t count_SNP = 0, count_SUB = 0, count_INS = 0, count_DEL = 0, count_AMP = 0, count_INV = 0, count_MOB = 0, count_CON = 0;

  diff_entry_list_t mutation_list = this->mutation_list();
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

      case AMP:
      {
        count_AMP++;
          
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);

        //Build duplicate sequence
        string dup;
        for (uint32_t i = 1; i < from_string<uint32_t>(mut["new_copy_number"]); i++)
          dup.append(new_ref_seq_info.get_sequence_1(mut[SEQ_ID], position, position+size-1));
        ASSERT(!dup.empty(), "Duplicate sequence is empy.");
          
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
        //ASSERT(!mut.entry_exists("ins_start") && !mut.entry_exists("ins_end") && !mut.entry_exists("del_start") && !mut.entry_exists("del_end"),
        //       "MOB: does not handle ins_start, ins_end, del_start, del_end yet.");
        ASSERT(mut["strand"] != "?", "Unknown repeat strand");
        
        count_MOB++;
        string new_seq_string = "";
        int32_t iDelStart = 0;
        int32_t iDelEnd = 0;
        int32_t iInsStart = 0;
        int32_t iInsEnd = 0;
        
        if(mut.entry_exists("del_start")) iDelStart = from_string<uint32_t>(mut["del_start"]);
        if(mut.entry_exists("del_end"))   iDelEnd = from_string<uint32_t>(mut["del_end"]);
        ASSERT((iDelStart >= 0) && (iDelEnd >= 0), (to_string(mut._type) + " " + mut._id) + " - NEGATIVE DELETION");

        // @JEB: correct here to look for where the repeat is in the original ref_seq_info.
        // This saves us from possibly looking at a shifted location...
        string rep_string = ref_seq_info.repeat_family_sequence(mut["repeat_name"], from_string<int16_t>(mut["strand"]));
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
        string duplicate_sequence = new_ref_seq_info.get_sequence_1(mut[SEQ_ID], position, position + (from_string<uint32_t>(mut["duplication_size"]) - 1));
        
        // If there are any inserts, put them in front of or behind the
        // repeat sequence.
        if(mut.entry_exists("ins_start")) {new_seq_string = mut["ins_start"] + rep_string;iInsStart = mut["ins_start"].length();}
        if(mut.entry_exists("ins_end"))   {new_seq_string += mut["ins_end"];iInsEnd = mut["ins_end"].length();}
        
        // Add on the duplicated sequence.  This happens AFTER
        // we have inserted any insertions.
        new_seq_string = duplicate_sequence + new_seq_string;        
        
        // Insert our newly minted sequence.
        new_ref_seq_info.insert_sequence_1(mut[SEQ_ID], position-1, new_seq_string, (to_string(mut._type) + " " + mut._id), verbose);
        
        // We've repeated the sequence, now it's time to repeat all the features
        // inside of and including the repeat region.
        new_ref_seq_info.repeat_feature_1(mut[SEQ_ID], position+iInsStart+duplicate_sequence.size(), iDelStart, iDelEnd, ref_seq_info, mut["repeat_name"], from_string<int16_t>(mut["strand"]), verbose);
          
        if (verbose)
        {
          cout << "MOB: 0" << " => " << new_seq_string << endl;
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
        
        Strand strand = (replace_start < replace_end) ?  POS_STRAND : NEG_STRAND;
        
        if (strand == NEG_STRAND) {
          swap(replace_start, replace_end);
        }
        
        // @JEB: correct here to look for where the replacing_sequence is in the original ref_seq_info.
        // This saves us from possible looking at a shifted location...
        string replacing_sequence = ref_seq_info[replace_target_id].get_sequence_1(replace_start, replace_end);

        if (strand == NEG_STRAND) {
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
  }
  
  //Cleanup.  If any of the sequences are of zero length, remove them.
  for (vector<cAnnotatedSequence>::iterator it_as = new_ref_seq_info.begin(); it_as < new_ref_seq_info.end(); it_as++) {
    if(!it_as->m_length){new_ref_seq_info.erase(it_as);it_as--;}
  }
  
  return new_ref_seq_info;
}


void cGenomeDiff::shift_positions(cDiffEntry &item, cReferenceSequences& ref_seq_info, bool verbose)
{  
  int32_t delta = item.mutation_size_change(ref_seq_info);
  if (verbose)
    cout << "Shift size: " << delta << endl;
  if (delta == UNDEFINED_INT32)
    WARN("Size change not defined for mutation.");

  uint32_t offset = from_string<uint32_t>(item[POSITION]);
  bool inversion = false;

  diff_entry_list_t muts = this->mutation_list();
  for (diff_entry_list_t::iterator itr = muts.begin(); itr != muts.end(); itr++) {
    cDiffEntry& mut = **itr;

    if (inversion) {
      WARN("shift_positions cannot handle inversions yet!");
    } else {
      uint32_t position = from_string<uint32_t>(mut[POSITION]);
      if (item[SEQ_ID] == mut[SEQ_ID] && position > offset)
        mut[POSITION] = to_string(position + delta);
    }

  }
}

//>! Return the
  
int32_t cDiffEntry::mutation_size_change(cReferenceSequences& ref_seq_info)
{  
  switch (this->_type) {
    case SNP: return 0; break;
    case SUB: return (*this)[NEW_SEQ].length() - from_string<uint32_t>((*this)[SIZE]); break;
    case INS: return (*this)[NEW_SEQ].length(); break;
    case DEL: return -(from_string<uint32_t>((*this)[SIZE])); break;
    case AMP: return from_string<uint32_t>((*this)[SIZE]) * (from_string<uint32_t>((*this)["new_copy_number"]) - 1); break;

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

cGenomeDiff cGenomeDiff::compare_genome_diff_files(const cGenomeDiff &control, const cGenomeDiff &test)
{
  cGenomeDiff control_gd = control;
  cGenomeDiff test_gd = test;
  cGenomeDiff new_gd; //ret val

  const diff_entry_list_t &control_mutations = control_gd.mutation_list();
  const diff_entry_list_t &test_mutations = test_gd.mutation_list();

  typedef set<cDiffEntry> diff_entry_set_t;
  diff_entry_set_t control_diff_mutations_set;
  diff_entry_set_t test_diff_mutations_set;

  //! Step: Strip counted_ptr<> for STL use.
  for (diff_entry_list_t::const_iterator it = control_mutations.begin();
       it != control_mutations.end(); it++) {
    control_diff_mutations_set.insert(**it);
  }
  for (diff_entry_list_t::const_iterator it = test_mutations.begin();
       it != test_mutations.end(); it++) {
    test_diff_mutations_set.insert(**it);
  }

  //! Step: True positive mutations
  diff_entry_set_t true_positive_mutations_set;
  set_intersection(control_diff_mutations_set.begin(), control_diff_mutations_set.end(),
                   test_diff_mutations_set.begin(), test_diff_mutations_set.end(),
                   inserter(true_positive_mutations_set, true_positive_mutations_set.begin()));

  //Remove true_positives from previous sets
  for (diff_entry_set_t::const_iterator it = true_positive_mutations_set.begin();
       it != true_positive_mutations_set.end(); it++) {
    control_diff_mutations_set.erase(*it);
    test_diff_mutations_set.erase(*it);
  }

  //Find difference and then determine if false-positive or false-negative
  diff_entry_set_t difference_mutations_set;
  set_difference(control_diff_mutations_set.begin(), control_diff_mutations_set.end(),
                   test_diff_mutations_set.begin(), test_diff_mutations_set.end(),
                   inserter(difference_mutations_set, difference_mutations_set.begin()));

  //! Step: False Negative
  diff_entry_set_t false_negative_mutations_set;
  set_intersection(difference_mutations_set.begin(), difference_mutations_set.end(),
                   control_diff_mutations_set.begin(), control_diff_mutations_set.end(),
                   inserter(false_negative_mutations_set, false_negative_mutations_set.begin()));
  //Remove false negatives from test set
  for (diff_entry_set_t::const_iterator it = true_positive_mutations_set.begin();
       it != true_positive_mutations_set.end(); it++) {
    test_diff_mutations_set.erase(*it);
  }

  //! Step: False positive
  diff_entry_set_t false_positive_mutations_set;
  set_intersection(difference_mutations_set.begin(), difference_mutations_set.end(),
                   test_diff_mutations_set.begin(), test_diff_mutations_set.end(),
                   inserter(false_positive_mutations_set, false_positive_mutations_set.begin()));

//! Step: Display results;
  const size_t &num_true_positive = true_positive_mutations_set.size();
  const size_t &num_false_negative = false_negative_mutations_set.size();
  const size_t &num_false_positive = false_positive_mutations_set.size();

  printf("#=TP|FN|FP	%i|%i|%i \t for %s versus %s \n",
         static_cast<int>(num_true_positive), static_cast<int>(num_false_negative), static_cast<int>(num_false_positive),
         control_gd._default_filename.c_str(), test_gd._default_filename.c_str());

  //! Step: Add validation=<TP/FP/FN> tags to new_gd.
  //True positive
  for (diff_entry_set_t::iterator it = true_positive_mutations_set.begin();
       it != true_positive_mutations_set.end(); it++) {
    cDiffEntry de = *it;
    de.insert(pair<string,string>("validation", "TP"));
    new_gd.add(de);
  }

  //False negative
  for (diff_entry_set_t::iterator it = false_negative_mutations_set.begin();
       it != false_negative_mutations_set.end(); it++) {
    cDiffEntry de = *it;
    de.insert(pair<string,string>("validation", "FN"));
    new_gd.add(de);
  }

   //False positive
  for (diff_entry_set_t::iterator it = false_positive_mutations_set.begin();
       it != false_positive_mutations_set.end(); it++) {
    cDiffEntry de = *it;
    de.insert(pair<string,string>("validation", "FP"));
    new_gd.add(de);
  }

  return new_gd;
}


}//namespace bresesq



