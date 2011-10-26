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
//! seq_id and positions are already parameters in diff_entry
//## mutations
(SNP,make_list<string> ("seq_id")("position")("new_seq"))
(SUB,make_list<string> ("seq_id")("position")("size")("new_seq"))
(DEL,make_list<string> ("seq_id")("position")("size"))
(INS,make_list<string> ("seq_id")("position")("new_seq"))
(MOB,make_list<string> ("seq_id")("position")("repeat_name")("strand")("duplication_size"))
(DEL,make_list<string> ("seq_id")("position")("size"))
(INV,make_list<string> ("seq_id")("position")("size"))
(AMP,make_list<string> ("seq_id")("position")("size")("new_copy_number"))
(CON,make_list<string> ("seq_id")("position")("size")("region"))

//## evidence
(RA,make_list<string> ("seq_id")("position")("insert_position")("ref_base")("new_base"))
(MC,make_list<string> ("seq_id")("start")("end")("start_range")("end_range"))
(JC,make_list<string> ("side_1_seq_id")("side_1_position")("side_1_strand")("side_2_seq_id")("side_2_position")("side_2_strand")("overlap"))
(UN,make_list<string> ("seq_id")("start")("end"))

//## validation
(CURA,make_list<string> ("expert"))
(FPOS,make_list<string> ("expert"))
(PHYL,make_list<string> ("gd"))
(TSEQ,make_list<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
(PFLP,make_list<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
(RFLP,make_list<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
(PFGE,make_list<string> ("seq_id")("enzyme"))
  
; // end line specifications


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
diff_entry::diff_entry(const gd_entry_type type)
: _type(type)
, _id("")
, _evidence()
{
}

diff_entry::diff_entry()
: _type(TYPE_UNKOWN)
, _id("")
, _evidence()
{
}

bool diff_entry::is_mutation() const
{
  const size_t size = sizeof(gd_entry_mutation_types) / sizeof(gd_entry_mutation_types[0]);
  return count(gd_entry_mutation_types, gd_entry_mutation_types + size, this->_type);
}

bool diff_entry::is_evidence() const
{
  const size_t size = sizeof(gd_entry_evidence_types) / sizeof(gd_entry_evidence_types[0]);
  return count(gd_entry_evidence_types, gd_entry_evidence_types + size, this->_type);
}

bool diff_entry::is_validation() const
{
  const size_t size = sizeof(gd_entry_validation_types) / sizeof(gd_entry_validation_types[0]);
  return count(gd_entry_validation_types, gd_entry_validation_types + size, this->_type);
}
  
  
//Comparing IDs here will currently break genome_diff::merge and genome_diff::subract
bool diff_entry::operator==(const diff_entry& de)
{
  diff_entry item = de;
  //! Case: Easy if not same type
  if (this->_type != item._type) {
    return false;
  }
  //! Case: Same type, but are fields that are common equal?
  else {
    // Get common keys
    const vector<key_t>& specs = line_specification[this->_type];
    for(vector<key_t>::const_iterator it = specs.begin();
        it != specs.end(); it++) {
      const key_t& spec(*it);
      if (_fields[spec] != item._fields[spec]) {
        return false;
      } else {
        continue;
      }
    }
    return true;
  }
}


/*! Marshal this diff entry into an ordered list of fields.
 */
void diff_entry::marshal(field_list_t& s) {
  s.push_back(to_string(_type));
	s.push_back(_id);
  
  s.push_back(join(_evidence, ","));

  
	// copy all fields:
	map_t cp=_fields;

	// marshal specified fields in-order, removing them from the copy after they've 
	// been printed:
  
  vector<string>& f = line_specification[_type];

  for (vector<string>::iterator it=f.begin(); it != f.end(); it++)
  {
		map_t::iterator iter=cp.find(*it);
    
    ASSERT(iter != cp.end(), "Did not find required field '" + *it + "' to write in entry id " + _id + " of type '" + to_string(_type) + "'.");
    
    s.push_back(iter->second);
    cp.erase(iter);
  
	}
	
	// marshal whatever's left, unless it's an empty field or _begins with an underscore
	for(map_t::iterator i=cp.begin(); i!=cp.end(); ++i) {
    
    assert(i->first.size());
    if (i->first.substr(0,1) == "_") continue;
    
    if (i->second.size() > 0) 
    {
      s.push_back(i->first + "=" + i->second);
    }
  }
}
  
vector<string> diff_entry::get_reject_reasons()
{
  vector<string> return_value;
  if (this->entry_exists("reject")) {
    return_value = split((*this)["reject"], ",");
  } 
  return return_value;
}

size_t diff_entry::number_reject_reasons()
{
  if(this->entry_exists(REJECT))
  {
    return this->get_reject_reasons().size();
  }
  return 0;
}
  
/*! Add reject reason to diff entry.
 */
void add_reject_reason(diff_entry& de, const string &reason) {

  if (de._fields.find(REJECT) == de._fields.end()) {
      de[REJECT] = reason;
  }
  // exists already, make comma separated list
  else {
    string reject = de[REJECT];
    reject += ",";
    reject +=reason; 
  }

}

uint32_t number_reject_reasons(diff_entry& de) {

	if (!de.entry_exists(REJECT) || de[REJECT].size() == 0)
		return 0;

	uint32_t reason_count = 1;
	for (uint32_t i = 0; i < de[REJECT].size(); i++)
		if (de[REJECT][i] == ',') reason_count++;
	return reason_count;
}


/*! Output operator for a diff entry.
 */
ostream& operator<<(ostream& out, diff_entry& de) {
	field_list_t fields;
	de.marshal(fields);
	out << join(fields, "\t");
	return out;
}
/*! Constructor.
 */
genome_diff::genome_diff(const string& filename)
 : _default_filename(filename)
 , _unique_id_counter(0) 
{
 read(filename);  
}

/*! Merge Constructor.
 */
genome_diff::genome_diff(genome_diff& merge1, genome_diff& merge2)
 : _unique_id_counter(0)
{
  // calling add() makes sure numbers are assigned appropriately
  for(diff_entry_list::iterator it=merge1._entry_list.begin(); it != merge1._entry_list.end(); it++)
  {
    this->add(*(it->get()));
  }
 
  for(diff_entry_list::iterator it=merge2._entry_list.begin(); it != merge2._entry_list.end(); it++)
  {
    this->add(*(it->get()));
  }  
}

  
uint32_t genome_diff::new_unique_id() 
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
void genome_diff::add(const diff_entry& item) {

  // allocating counted_ptr takes care of disposal
  diff_entry* diff_entry_copy = new diff_entry(item); 
  diff_entry_ptr added_item(diff_entry_copy);
  _entry_list.push_back(added_item);
    
  if ((added_item->_id.size() == 0) || unique_id_used.count(from_string<uint32_t>(added_item->_id)))
  {
    uint32_t new_id = new_unique_id();
    added_item->_id = to_string(new_id); 
    unique_id_used[new_id] = true;
  }
  else
  {
    uint32_t new_id = from_string<uint32_t>(added_item->_id);
    unique_id_used[new_id] = true;
  }
  
// # sub add
// # {
// #   my ($self, $item) = @_;
// #   my @missing_required_columns = ();
// # 
// #   ## no ID, give it a new one (need to re-assign id's later...)
// #   if ( !defined $item->{id} )
// #   {
// #     $item->{id} = $self->new_unique_id;
// #   }
// #   elsif ( $self->used_unique_id( $item->{id}) && !(($item->{id} eq '.') || ($item->{id} eq '+') || ($item->{id} eq '?')))
// #   {
// #     $self->warn("Ignoring attempt to add item with an existing id: $item->{id}");
// #     return;
// #   }
// # 
// #   sub check_required_field
// #   {
// #     my ($item, $field, $missing_ref) = @_;
// #     push @$missing_ref, $field if (!defined $item->{$field});
// #   }
// #   
// #   ## check to be sure the item has required fields, or auto-populate them
// #   $item->{type} = '' if (!defined $item->{type});
// #   
// #   my $spec = $line_specification->{$item->{type}};
// #   if (!defined $spec)
// #   {
// #     $self->warn("Type \'$item->{type}\' is not recognized. Ignoring item.");
// #     return;
// #   }
// #   
// #   ## check for required fields
// #   foreach my $key (@$spec)
// #   {
// #     check_required_field($item, $key, \@missing_required_columns);
// #   }
// # 
// #   if (scalar @missing_required_columns > 0)
// #   {
// #     $self->warn("GenomeDiff::Ignoring item of type \'$item->{type}\' that is missing required field(s):" . join (',', @missing_required_columns));
// #     return;
// #   }
// # 
// #   ## these are all required columns
// #   $item->{SORT_1} = $tag_sort_fields->{$item->{type}}->[0];
// #   $item->{SORT_2} = $item->{$tag_sort_fields->{$item->{type}}->[1]};
// #   $item->{SORT_3} = $item->{$tag_sort_fields->{$item->{type}}->[2]};
// #   
// #   push @{$self->{list}}, $item;
// #   $self->mark_unique_id($item->{id});
// # }
  
  
  
}

  //! Subtract mutations using gd_ref as reference.
  void genome_diff::subtract(genome_diff& gd_ref, bool verbose)
  {
    //Iterate through all the entries
    for (diff_entry_list::iterator it = _entry_list.begin(); it != _entry_list.end(); it++)
    {
      //The current entry we're looking at
      diff_entry& entry = **it;
      
      //Is the entry a mutation?
      if(entry.is_mutation())
      {
        //Iterate through all the entries we're checking against.
        for (diff_entry_list::iterator it_ref = gd_ref._entry_list.begin(); it_ref != gd_ref._entry_list.end(); it_ref++)
        {
          //The current entry we're looking at
          diff_entry& entry_ref = **it_ref;
          
          //Does the current entry match any of the reference entries?
          if(entry == entry_ref)
          {
            //Notify the user of the action.
            //This will currently output the enumerated number representing the _type.
            if(verbose){cout << "REMOVE\t" << entry._type << "\t" << entry._id << endl;}
            _entry_list.erase(it);
            it--;
          }            
        }
      }        
    }
  }
  
  //! Merge GenomeDiff information using gd_new as potential new info.
  void genome_diff::merge(genome_diff& gd_new, bool verbose)
  {
    uint32_t old_unique_ids = unique_id_used.size();
    
    //Iterate through all the potential new entries
    for (diff_entry_list::iterator it_new = gd_new._entry_list.begin(); it_new != gd_new._entry_list.end(); it_new++)
    {
      //The current potential new entry we're looking at
      diff_entry& entry_new = **it_new;
      bool new_entry = true;
      
      //Iterate through all the entries in the current list.
      for (diff_entry_list::iterator it_cur = _entry_list.begin(); it_cur != _entry_list.end(); it_cur++)
      {
        //The current entry we're looking at
        diff_entry& entry = **it_cur;
        
        //Does the new entry match the current entry?
        if(entry == entry_new)
        {
          //Existing matching entry found, this is not new
          bool new_entry = false;
          break;
        }
      }
      
      //We definately have a new entry
      if(new_entry)
      {
        //Notify user of new entry
        if(verbose)cout << "NEW ENTRY\t" << entry_new._id << "\t" << gd_new._default_filename << endl;
        
        //Add the new entry to the existing list
        add(entry_new);        
      }
    }
    
    //Iterate through all the entries in the new list.
    //This is where we update the evidence IDs for mutations.
    for (diff_entry_list::iterator it = _entry_list.begin(); it != _entry_list.end(); it++)
    {
      //Is this one of the new entries?
      if(from_string<uint32_t>((**it)._id) > old_unique_ids)
      {                
        //For every piece of evidence this entry has
        for(uint32_t iter = 0; iter < (**it)._evidence.size(); iter++)
        {
          bool found_match = false;
          
          //Iterate through all the potential new entries
          for (diff_entry_list::iterator it_new = gd_new._entry_list.begin(); it_new != gd_new._entry_list.end(); it_new++)
          {            
            //Does this evidence ID match an ID in the old list?
            if((**it)._evidence[iter] == (**it_new)._id && !found_match)
            {
              //Iterate through all the current entries
              for (diff_entry_list::iterator it_cur =_entry_list.begin(); it_cur != _entry_list.end(); it_cur++)
              {
                //Does the new entry match the current entry?
                if((**it_cur) == (**it_new))
                {
                  //Change the evidence ID to it's new ID in the new updated list
                  (**it)._evidence[iter] = (**it_cur)._id;
                  found_match = true;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }


/*! Read a genome diff(.gd) from the given file to class member
  _entry_list
 */

void genome_diff::read(const string& filename) {
  ifstream IN(filename.c_str());
  if(!IN.good())
    cerr << "Could not open file for reading: " << filename << endl;

  ::list<string> lines;
  string line;
  
  if(IN.is_open()) {
    while(IN.good()) {
      getline(IN,line);
      lines.push_back(line);
    }
    IN.close();
  }
  //Can't find version in first line
  if(lines.front().find("#=GENOME_DIFF") == string::npos)
  {
    cerr << "Could not match version line in file" << endl;
    assert(false);
  }
  
  //strip paths and .gd off of filename to obtain run_id
  //EX: a filename of 
  //RS0001_Woods2011/RJW1129.gd
  //or
  //RJW1129.gd
  //becomes RJW1129
  size_t run_id_end = filename.find(".gd");
  if(size_t run_id_begin = filename.find_last_of("\\") == string::npos) {
    metadata.run_id = filename.substr(0, run_id_end);
  } else {
    metadata.run_id = filename.substr(run_id_begin, run_id_end);
  }
 
  // Check for and read in metadata
  //#=GENOME_DIFF 1.0
  //#=AUTHOR    Barrick JE 
  //#=REFSEQ    Genbank:NC_012967.1
  //#=READSEQ   SRA:SRR066595
  while (lines.front().find_first_of("#") == 0)
  {
    string line = lines.front();
  
    //#=GENOME_DIFF 1.0
    if (line.find("#=GENOME_DIFF") != string::npos) {
      vector<string> version_split = split(line, " ");
      metadata.version = version_split.back();
    }

    //#=AUTHOR    Barrick JE 
    if (line.find("#=AUTHOR") != string::npos) {
      size_t author_begins = line.find_first_not_of(" ", 8);
      size_t author_ends = line.length();
      metadata.author = line.substr(author_begins, author_ends);
    }

    //#=REFSEQ    Genbank:NC_012967.1
    if (line.find("#=REFSEQ") != string::npos) {
      vector<string> ref_seq_split = split(line, ":");
      metadata.ref_seq = ref_seq_split.back(); 
    }

    //#=READSEQ   SRA:SRR066595 
    if (line.find("#=READSEQ") != string::npos) {
      vector<string> read_seq_split = split(line, ":");
      metadata.read_seq.push_back(read_seq_split.back());
    }

    lines.pop_front();
  }
  
  
  while(!lines.empty())
  {    
    string line = lines.front();
    //Check that line isn't empty
    if (line.empty()) {
      lines.erase(lines.begin());
      continue;
    }
    //Remove lines containing evidence
//    string type = split(line, "\t").front();
//   if(type.length() > 3) {
//      lines.erase(lines.begin());
//      continue;
//    }
      
    // # $self->add($self->_line_to_item($l));
    add(_line_to_item(line));
    
    lines.pop_front();
  }

  // match common fields - type id pid seqid
		
  // match type-specific fields

 	
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
bool diff_entry_sort(const diff_entry_ptr& a, const diff_entry_ptr& b) {

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



/*! Write this genome diff to a file.
 */
void genome_diff::write(const string& filename) {
	ofstream ofs(filename.c_str());
  ofs << "#=GENOME_DIFF 1.0" << endl;
  //! Place any extra data gathered in breseq pipeline in header
  if (!this->metadata.breseq_data.empty()) {
    for (map<key_t,string>::iterator it = metadata.breseq_data.begin();
         it != metadata.breseq_data.end(); it ++) {
      const key_t& key = it->first;
      const string& value = it->second;
      ofs << "#" << key << "=" << value << endl;
    }
  }
  
  // sort
  _entry_list.sort(diff_entry_sort);
  
	for(diff_entry_list::iterator i=_entry_list.begin(); i!=_entry_list.end(); ++i) {
		ofs << (**i) << endl;
	}
	ofs.close();
}
  
//! Call to assure that that this GenomeDiff seq_id
//  matches the passed in cReferenceSequences seq_id
//  returns true if it matches
//  returns false if it doesn't match
bool genome_diff::is_valid_seq_id(cReferenceSequences& ref_seq_info, bool verbose)
{  
  if((**(_entry_list.begin()))[SEQ_ID] != ref_seq_info.begin()->m_seq_id)
  {
    cout << "LOADED REFERENCE\t" << ref_seq_info.begin()->m_seq_id << endl;
    cout << "LOADED GENOMEDIFF\t" << (**(_entry_list.begin()))[SEQ_ID] << endl;
    cout << "LOADED FILES\t\tDIFFER" << endl;
    return false;
  }
  
  if(verbose)
  {
    cout << "LOADED REFERENCE\t" << ref_seq_info.begin()->m_seq_id << endl;
    cout << "LOADED GENOMEDIFF\t" << (**(_entry_list.begin()))[SEQ_ID] << endl;
    cout << "LOADED FILES\t\tMATCH" << endl;
  } 
  
  return true;
}


void genome_diff::add_breseq_data(const key_t &key, const string& value)
{
  this->metadata.breseq_data.insert(pair<string,string>(key, value));
}


// Convert GD file to GVF file
void GDtoGVF( const string &gdfile, const string &gvffile, bool snv_only){
    
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

  for( size_t i=0; i<features.size(); i++ )
  {
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
      
      gvf[8].append(";Reference_seq=").append( evidence[5] );
      
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

/*! Given a list of types, search and return the diff_entry's within diff_entry_list whose 
 * _type parameter matches one of those within input types. 
 */ 
diff_entry_list genome_diff::list(const vector<gd_entry_type>& types)
{
  // default is to have to types
  if (types.size() == 0)
    return _entry_list;
  
  diff_entry_list return_list;
  
  for (diff_entry_list::iterator itr_diff_entry = _entry_list.begin(); 
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
  
diff_entry_list genome_diff::show_list(const vector<gd_entry_type>& types)
{
  diff_entry_list ret_list = list(types);
  ret_list.remove_if(diff_entry::fields_exist(make_list<diff_entry::key_t>("deleted")));
  ret_list.remove_if(diff_entry::fields_exist(make_list<diff_entry::key_t>("no_show")));  
  return ret_list;
}

/*-----------------------------------------------------------------------------
 * returns entries NOT used as evidence by other entries. 
 *
 *-----------------------------------------------------------------------------*/
diff_entry_list genome_diff::filter_used_as_evidence(const diff_entry_list& input)
{
  // first we make a map with everything used as evidence by any entry in the entire genome diff
  map<string,bool> used_as_evidence;
  for (diff_entry_list::const_iterator it = _entry_list.begin(); it != _entry_list.end(); it++) 
  {
    const diff_entry_ptr& de = *it;
    for (vector<string>::const_iterator ev_it = de->_evidence.begin(); ev_it != de->_evidence.end(); ev_it++) 
    {  
      used_as_evidence[*ev_it] = true;
    }   
  }
  
  // then construct a list of all items in input with ids not in this map
  diff_entry_list return_list;
  for (diff_entry_list::const_iterator it = input.begin(); it != input.end(); it++) 
  {
    const diff_entry_ptr& de = *it;
    if ( !used_as_evidence.count(de->_id) )
      return_list.push_back(de);
  }

  return return_list;
}


diff_entry genome_diff::_line_to_item(const string& line)
{
  list_t line_list = split(line, "\t");
  diff_entry item;
  item._type = to_type(shift<string>(line_list));
  item._id = shift<string>(line_list);
  string evidence_string = shift<string>(line_list);
  item._evidence = split(evidence_string, ",");

  const list_t spec = line_specification[item._type];

  // make sure it is a recognized type
  ASSERT(!spec.empty(), "Type '" + to_string(item._type) + "' is not recognized for line:\n" + line );

  for(size_t i = 0; i < spec.size(); i++)
  {
    string key = spec[i];
    string next = shift<string>(line_list);

    if(next.empty())
    {
     WARN("Number of required items is less than expected for type " + to_string(item._type));
      assert(false);
    }

    ASSERT(next.find("=") == string::npos,
            "Unexpected key=value pair '" + next + "' encountered for required item '" + key 
            + "' in type '" + to_string(item._type) + "' line:\n" + line);
    
    item[key] = next;
  }


  for(list_t::iterator itr = line_list.begin();
      itr != line_list.end(); itr ++)
  {
    string key_value_pair(*itr); 
    if(key_value_pair.empty()) continue;
    //assert(regex_m("=",key_value_pair));
    assert(key_value_pair.find("=") != string::npos); 
    vector<string> matched = split(key_value_pair,"=");

    if(matched.empty())
    {
      cerr << "Not a key value pair" << key_value_pair <<  endl << line;
      assert(false);
      continue;
    }
    assert(matched.size()==2);
    string item_key = matched[0];
    string item_value = matched[1]; 
    item[item_key] = item_value;
  }
  
 return item;
}

// return items with types that are 3 characters long
diff_entry_list genome_diff::mutation_list()
{
  diff_entry_list mut_list;

   for(diff_entry_list::iterator itr = _entry_list.begin();
       itr != _entry_list.end(); itr ++) {
     diff_entry& item = **itr;
     if(item.is_mutation()) {
       mut_list.push_back(*itr);
     }
   }

  return mut_list;
}

// return items with types that are 2 characters long
diff_entry_list genome_diff::evidence_list()
{
  diff_entry_list mut_list;
  
  for(diff_entry_list::iterator itr = _entry_list.begin();
      itr != _entry_list.end(); itr ++) {
    diff_entry& item = **itr;
    if(item.is_evidence()) {
      mut_list.push_back(*itr);
    }
  }
  
  return mut_list;
}


/*! Return all diff_entrys within _entry_list whose _id matches one
 * of those within input's item._evidence
 */ 
diff_entry_list genome_diff::mutation_evidence_list(const diff_entry& item)
{
  diff_entry_list return_list;
  
  //return diff_entries with matching evidence
  for (vector<string>::const_iterator itr_i = item._evidence.begin(); itr_i != item._evidence.end(); itr_i ++) 
  {  
    const string& evidence = *itr_i;
    
    for (diff_entry_list::iterator itr_j = _entry_list.begin(); itr_j != _entry_list.end(); itr_j ++) 
    {  
      diff_entry& entry = **itr_j;
    
      if (entry._id == evidence)
        return_list.push_back(*itr_j);
    }   
  }
  return return_list;
}

diff_entry_ptr genome_diff::parent(const diff_entry& item)
{
  for(diff_entry_list::iterator itr_test_item = _entry_list.begin();
      itr_test_item != _entry_list.end(); itr_test_item ++) { 
    diff_entry& test_item = **itr_test_item;
    for(vector<string>::iterator itr = test_item._evidence.begin();
        itr != test_item._evidence.end(); itr ++) { 
      string& test_evidence_id = (*itr);
      if(test_evidence_id == item._id)      
        return diff_entry_ptr(*itr_test_item);
    }
  }
  return diff_entry_ptr(NULL);
}

bool genome_diff::mutation_unknown(diff_entry mut)
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

void genome_diff::add_reject_reasons(diff_entry item, const string& reject)
{
  if (item.entry_exists(REJECT))
    item[REJECT] += ",";
  else 
    item[REJECT] += reject;
}

bool 
genome_diff::interval_un(const uint32_t& start,const uint32_t& end)
{
  diff_entry_list un_list = list(make_list<gd_entry_type>(UN));

  for (diff_entry_list::iterator itr = un_list.begin();
       itr != un_list.end(); itr++) {
    diff_entry un(**itr);

    if (start >= from_string<uint32_t>(un[START]) &&
        end <= from_string<uint32_t>(un[END])) {
      return true;
    }
  }
  return false;
}

cReferenceSequences genome_diff::apply_to_sequences(cReferenceSequences& ref_seq_info, bool verbose)
{
  // copy the reference sequence info
  cReferenceSequences new_ref_seq_info(ref_seq_info);
  
  ASSERT(is_valid_seq_id(new_ref_seq_info, verbose), "Loaded files do not match");
    
  uint32_t count_SNP = 0, count_SUB = 0, count_INS = 0, count_DEL = 0, count_AMP = 0, count_INV = 0, count_MOB = 0, count_CON = 0;

  diff_entry_list mutation_list = this->mutation_list();
  for (diff_entry_list::iterator itr_mut = mutation_list.begin(); itr_mut != mutation_list.end(); itr_mut++) 
  {
    diff_entry& mut(**itr_mut);
    uint32_t position = from_string<uint32_t>(mut[POSITION]);
        
    switch (mut._type) 
    {
      case SNP :
      {
        count_SNP++;
        string type_temp = "SNP"; 
        
        if (verbose)
        {
          cout << "SNP: 1 bp => " << mut[NEW_SEQ] << " at position " << position << endl;
          cout << "  shift: none" << endl;
        }
        
        new_ref_seq_info.replace_sequence_1(mut[SEQ_ID], position, position, mut[NEW_SEQ], (type_temp + " " + mut._id), verbose);
      } break;

      case SUB:
      {
        count_SUB++;
        string type_temp = "SUB";
          
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);
          
        if (verbose)
        {
          cout << "SUB: " << size << " => " << mut[NEW_SEQ] << endl;
          cout << "   shift +" << mut[NEW_SEQ].length() << " bp at position " <<  position << endl;
        }

        new_ref_seq_info.replace_sequence_1(mut[SEQ_ID], position, position + size - 1, mut[NEW_SEQ], (type_temp + " " + mut._id), verbose);
      } break;

      case INS:
      {          
        count_INS++;
        string type_temp = "INS";
          
        if (verbose)
        {
          cout << "INS: 0 => " << mut[NEW_SEQ] << endl;
          cout << "   shift +" << mut[NEW_SEQ].length() << " bp at position " <<  position << endl;
        }
          
        new_ref_seq_info.insert_sequence_1(mut[SEQ_ID], position, mut[NEW_SEQ], (type_temp + " " + mut._id), verbose);
      } break;

      case DEL:
      {
        count_DEL++;
        string type_temp = "DEL";
          
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);
          
        if (verbose)
        {
          cout << "DEL: " << size << " => 0" << endl;
          cout << "   shift -" << size << " bp at position " <<  position << endl;
        }
          
        new_ref_seq_info.replace_sequence_1(mut[SEQ_ID], position, position + size -1, "", (type_temp + " " + mut._id), verbose);
      } break;

      case AMP:
      {
        count_AMP++;
        string type_temp = "AMP";
          
        const uint32_t& size = from_string<uint32_t>(mut[SIZE]);

        //Build duplicate sequence
        string dup;
        for (uint32_t i = 1; i < from_string<uint32_t>(mut["new_copy_number"]); i++)
          dup.append(new_ref_seq_info.get_sequence_1(mut[SEQ_ID], position, position+size-1));
        ASSERT(!dup.empty(), "Duplicate sequence is empy.");
          
        if (verbose)
        {
          cout << "AMP: 0" << " => " << dup << endl;
          cout << "   shift +" << dup.length() << " bp at position " << position << endl;
        }

        new_ref_seq_info.insert_sequence_1(mut[SEQ_ID], position-1, dup, (type_temp + " " + mut._id), verbose);
      } break;
        
      case INV:
      {
        count_INV++;
          
        WARN("INV: mutation type not handled yet");
      } break;

      case MOB:
      {
        count_MOB++;
        string type_temp = "MOB";
          
        ASSERT(!mut.entry_exists("ins_start") && !mut.entry_exists("ins_end") && !mut.entry_exists("del_start") && !mut.entry_exists("del_end"),
               "MOB: does not handle ins_start, ins_end, del_start, del_end yet.");
        ASSERT(mut["strand"] != "?", "Unknown repeat strand");

        // @JEB: correct here to look for where the repeat is in the original ref_seq_info.
        // This saves us from possible looking at a shifted location...
        string seq_string = ref_seq_info.repeat_family_sequence(mut["repeat_name"], from_string<int16_t>(mut["strand"]));
        mut["repeat_size"] = to_string(seq_string.length()); // saving this for shifting

        string duplicate_sequence = new_ref_seq_info.get_sequence_1(mut[SEQ_ID], position, position + from_string<uint32_t>(mut["duplication_size"]) - 1);
          
        if (verbose)
        {
          cout << "MOB: 0" << " => " << duplicate_sequence + seq_string << endl;
          cout << "   shift +" << (duplicate_sequence.length() + seq_string.length())  << " bp at position " << position << endl;
        }
          
        new_ref_seq_info.insert_sequence_1(mut[SEQ_ID], position-1, duplicate_sequence + seq_string, (type_temp + " " + mut._id), verbose);
      } break;

      case CON:
      {
        count_CON++;
        string type_temp = "CON";
          
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
          
        if (verbose)
        {
          cout << "CON: " << displaced_sequence << " => " << replacing_sequence << endl;
          int32_t size_change = static_cast<int32_t>(replacing_sequence.length()) - static_cast<int32_t>(displaced_sequence.length());
          cout << "   shift " << ((size_change >= 0) ? "+" : "") << size_change << " bp at position " << position << endl;
        }

        new_ref_seq_info.replace_sequence_1(mut[SEQ_ID], position, position + size - 1, replacing_sequence, (type_temp + " " + mut._id), verbose);        
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
  
  return new_ref_seq_info;
}


void genome_diff::shift_positions(diff_entry &item, cReferenceSequences& ref_seq_info, bool verbose)
{
  ASSERT(is_valid_seq_id(ref_seq_info), "Loaded files do not match");
  
  int32_t delta = item.mutation_size_change(ref_seq_info);
  if (verbose)
    cout << "Shift size: " << delta << endl;
  if (delta == UNDEFINED_INT32)
    WARN("Size change not defined for mutation.");

  uint32_t offset = from_string<uint32_t>(item[POSITION]);
  bool inversion = false;

  diff_entry_list muts = this->mutation_list();
  for (diff_entry_list::iterator itr = muts.begin(); itr != muts.end(); itr++) {
    diff_entry& mut = **itr;

    if (inversion) {
      WARN("shift_positions cannot handle inversions yet!");
    } else {
      int32_t position = from_string<int32_t>(mut[POSITION]);
      if (from_string<uint32_t>(mut[POSITION]) > offset)
        mut[POSITION] = to_string(position + delta);
    }

  }
}

//>! Return the
  
int32_t diff_entry::mutation_size_change(cReferenceSequences& ref_seq_info)
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


}//namespace bresesq



