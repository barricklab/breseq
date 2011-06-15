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

#include "breseq/genome_diff.h"
#include <boost/algorithm/string.hpp>
#include <boost/assign/list_of.hpp>

using namespace breseq;

// Common keywords used for diff entries:
const char* breseq::TYPE="type";
const char* breseq::ID="id";
const char* breseq::PID="parent";
const char* breseq::SEQ_ID="seq_id";
const char* breseq::START="start";
const char* breseq::END="end";
const char* breseq::START_RANGE="start_range";
const char* breseq::END_RANGE="end_range";
const char* breseq::LEFT_OUTSIDE_COV="left_outside_cov";
const char* breseq::LEFT_INSIDE_COV="left_inside_cov";
const char* breseq::RIGHT_INSIDE_COV="right_inside_cov";
const char* breseq::RIGHT_OUTSIDE_COV="right_outside_cov";
const char* breseq::POSITION="position";
const char* breseq::INSERT_POSITION="insert_position";
const char* breseq::QUALITY="quality";
const char* breseq::POLYMORPHISM_QUALITY="polymorphism_quality";
const char* breseq::REF_BASE="ref_base";
const char* breseq::NEW_BASE="new_base";
const char* breseq::FREQUENCY="frequency";
const char* breseq::REJECT="reject";
const char* breseq::REF_COV="ref_cov";
const char* breseq::NEW_COV="new_cov";
const char* breseq::TOT_COV="tot_cov";
const char* breseq::ERROR="error";

// Types of diff entries:
const char* breseq::SNP="SNP";
const char* breseq::SUB="SUB";
const char* breseq::DEL="DEL";
const char* breseq::INS="INS";
const char* breseq::MOB="MOB";
const char* breseq::AMP="AMP";
const char* breseq::INV="INV";
const char* breseq::CON="CON";

const char* breseq::RA="RA";
const char* breseq::MC="MC";
const char* breseq::JC="JC";
const char* breseq::UN="UN";

//our $line_specification = {
//## mutations
//'SNP' => ['seq_id', 'position', 'ref_seq', 'new_seq'],
//'SUB' => ['seq_id', 'position', 'ref_seq', 'new_seq'],
//'DEL' => ['seq_id', 'position', 'size'],
//'INS' => ['seq_id', 'position', 'new_seq'],
//'MOB' => ['seq_id', 'position', 'repeat_name', 'strand', 'duplication_size', 'gap_left', 'gap_right'],
//#	'MOB' => ['seq_id', 'position', 'repeat_name', 'strand', 'duplication_size', 'del_left', 'del_right', 'ins_left', 'ins_right'],	
//'DUP' => ['seq_id', 'position', 'size'],
//'INV' => ['seq_id', 'position', 'size'],
//
//## evidence
//'RA' => ['seq_id', 'position', 'insert_position', 'ref_base', 'new_base'],
//'MC' => ['seq_id', 'start', 'end'],
//'JC' => ['side_1_seq_id', 'side_1_position', 'side_1_strand', 'side_2_seq_id', 'side_2_position', 'side_2_strand', 'overlap'],
//'UN' => ['seq_id', 'start', 'end'],
//};
//
//our $tag_sort_fields = {
//'SNP' => [1, 'seq_id', 'position'],
//'SUB' => [1, 'seq_id', 'position'],
//'DEL' => [1, 'seq_id', 'position'],
//'INS' => [1, 'seq_id', 'position'],
//'MOB' => [1, 'seq_id', 'position'],
//'DUP' => [1, 'seq_id', 'position'],
//'INV' => [1, 'seq_id', 'position'],
//'RA' => [2, 'seq_id', 'position'],
//'MC' => [2, 'seq_id', 'start'],
//'JC' => [2, 'side_1_seq_id', 'side_1_position'],
//'UN' => [3, 'seq_id', 'start'],
//};

// Field order.
static const char* s_field_order[] = { 
breseq::SEQ_ID,
breseq::POSITION,
breseq::INSERT_POSITION,
breseq::REF_BASE,
breseq::NEW_BASE,
breseq::START,
breseq::END,
breseq::START_RANGE,
breseq::END_RANGE,
0}; // required trailing null.


/*! Constructor.
 */
breseq::diff_entry::diff_entry(const string& t, const string& id, const string& parents)
: _type(t)
, _id(id)
, _parents(parents) {
}


/*! Marshal this diff entry into an ordered list of fields.
 */
void breseq::diff_entry::marshal(field_list_t& s) {
	s.push_back(_type);
	s.push_back(_id);
	s.push_back(_parents);
	
	// copy all fields:
	map_t cp=_fields;

	// marshal specified fields in-order, removing them from the copy after they've 
	// been printed:
	const char* f=s_field_order[0];
	for(size_t i=0; ; ++i) {
		f = s_field_order[i];
		if(f == 0) { break; }
	
		map_t::iterator iter=cp.find(f);
		if(iter != cp.end()) {
			s.push_back(iter->second);
			cp.erase(iter);
		}
	}
	
	// marshal whatever's left:
	for(map_t::iterator i=cp.begin(); i!=cp.end(); ++i) {
		s.push_back(i->first + "=" + i->second);
	}
}

/*! Add reject reason to diff entry.
 */
void breseq::add_reject_reason(diff_entry& de, const string &reason) {

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


/*! Output operator for a diff entry.
 */
ostream& breseq::operator<<(ostream& out, breseq::diff_entry& de) {
	field_list_t fields;
	de.marshal(fields);
	out << join(fields, "\t");
	return out;
}


/*!
 */
breseq::ra::ra(const string& id, const string& parents) : diff_entry(RA, id, parents) {
}


/*!
 */
breseq::mc::mc(const string& id, const string& parents) : diff_entry(MC, id, parents) {
}


/*!
 */
breseq::un::un(const string& id, const string& parents) : diff_entry(UN, id, parents) {
}


/*! Add evidence to this genome diff.
 */
void breseq::genome_diff::add(const diff_entry& v) {
	boost::shared_ptr<diff_entry> de(v.clone());
	_entry_list.push_back(de); 
}


/*! Read a genome diff from the given file.
 */
void breseq::genome_diff::read(const string& filename) {
	using namespace std;
	ifstream ifs(filename.c_str());
	
	while(!ifs.eof()) {
		string line;
		getline(ifs, line); // read a line upto eof or '\n'.
		
		// strip all characters trailing a '#', unless it's escaped.
		
		// if !line.empty
		
		// match common fields - type id pid seqid
		
		// match type-specific fields
		
	}
}


map<string, sort_fields_item> diff_entry_sort_fields = make_map<string, sort_fields_item>
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

map<string, uint8_t> sort_order = make_map<string, uint8_t>
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
bool breseq::diff_entry_sort(boost::shared_ptr<diff_entry> a, boost::shared_ptr<diff_entry> b) {

  string a_type = a->_type;
  string b_type = b->_type;

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
  uint32_t b_sort_field_3 = from_string<uint32_t>((*b)[b_sort_fields._f3]);;
  
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
  
  return false;
}



/*! Write this genome diff to a file.
 */
void breseq::genome_diff::write(const string& filename) {
	ofstream ofs(filename.c_str());
	ofs << "#=GENOME_DIFF 1.0" << endl;
  
  // sort
  sort(_entry_list.begin(), _entry_list.end(), diff_entry_sort);
  
	for(entry_list_t::iterator i=_entry_list.begin(); i!=_entry_list.end(); ++i) {
		ofs << (**i) << endl;
	}
	ofs.close();
}
