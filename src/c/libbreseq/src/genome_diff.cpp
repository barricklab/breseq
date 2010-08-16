#include "breseq/genome_diff.h"
#include <boost/algorithm/string.hpp>
#include <fstream>

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
const char* breseq::REF_BASE="ref_base";
const char* breseq::NEW_BASE="new_base";
const char* breseq::FREQUENCY="frequency";
const char* breseq::REJECT="reject";
const char* breseq::REF_COV="ref_cov";
const char* breseq::NEW_COV="new_cov";
const char* breseq::TOT_COV="tot_cov";

// Types of diff entries:
const char* breseq::RA="RA";
const char* breseq::MC="MC";
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
0}; // required trailing null.

// These fields require a modification for 0- vs. 1-indexing before IO.
static const char* s_requires_index_mod[] = {
breseq::POSITION,
breseq::START,
breseq::END,
0}; // required trailing null.


/*! Constructor.
 */
breseq::diff_entry::diff_entry(const std::string& t, const std::string& id, const std::string& parents)
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
	
	// @dk: this is a weird one - position-related fields need to be altered due to
	// differences in indexing base (0- vs. 1-indexing).
	const char* f=s_requires_index_mod[0];
	for(std::size_t i=0; ; ++i) {
		f = s_requires_index_mod[i];
		if(f == 0) { break; }
		
		map_t::iterator iter=cp.find(f);
		if(iter != cp.end()) {
			iter->second = boost::get<uint32_t>(iter->second) + 1;
		}
	}	
	

	// marshal specified fields in-order, removing them from the copy after they've 
	// been printed:
	f=s_field_order[0];
	for(std::size_t i=0; ; ++i) {
		f = s_field_order[i];
		if(f == 0) { break; }
	
		map_t::iterator iter=cp.find(f);
		if(iter != cp.end()) {
			string_visitor v;
			boost::apply_visitor(v,iter->second);
			s.push_back(v.str());
			cp.erase(iter);
		}
	}
	
	// marshal whatever's left:
	for(map_t::iterator i=cp.begin(); i!=cp.end(); ++i) {
		string_visitor v;
		boost::apply_visitor(v,i->second);
		s.push_back(i->first + "=" + v.str());
	}
}


/*! Output operator for a diff entry.
 */
std::ostream& breseq::operator<<(std::ostream& out, breseq::diff_entry& de) {
	field_list_t fields;
	de.marshal(fields);
	out << boost::join(fields, "\t");
	return out;
}


/*!
 */
breseq::ra::ra(const std::string& id, const std::string& parents) : diff_entry(RA, id, parents) {
}


/*!
 */
breseq::mc::mc(const std::string& id, const std::string& parents) : diff_entry(MC, id, parents) {
}


/*!
 */
breseq::un::un(const std::string& id, const std::string& parents) : diff_entry(UN, id, parents) {
}


/*! Add evidence to this genome diff.
 */
void breseq::genome_diff::add(const diff_entry& v) {
	boost::shared_ptr<diff_entry> de(v.clone());
	_entry_list.push_back(de); 
}


/*! Read a genome diff from the given file.
 */
void breseq::genome_diff::read(const std::string& filename) {
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


/*! Write this genome diff to a file.
 */
void breseq::genome_diff::write(const std::string& filename) {
	std::ofstream ofs(filename.c_str());
	ofs << "#=GENOME DIFF 1.0" << std::endl;
	for(entry_list_t::iterator i=_entry_list.begin(); i!=_entry_list.end(); ++i) {
		ofs << (**i) << std::endl;
	}	
}
