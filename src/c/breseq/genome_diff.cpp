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
#include <list>

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


// Types of diff entries:
const char* SNP="SNP";
const char* SUB="SUB";
const char* DEL="DEL";
const char* INS="INS";
const char* MOB="MOB";
const char* AMP="AMP";
const char* INV="INV";
const char* CON="CON";

const char* RA="RA";
const char* MC="MC";
const char* JC="JC";
const char* UN="UN";

//our $line_specification = {
map<string, vector<string> > line_specification =make_map<string, vector<string> > 
//! seq_id and positions are already parameters in diff_entry
//## mutations
("SNP",make_list<string> ("seq_id")("position")("new_seq"))
("SUB",make_list<string> ("seq_id")("position")("size")("new_seq"))
("DEL",make_list<string> ("seq_id")("position")("size"))
("INS",make_list<string> ("seq_id")("position")("new_seq"))
("MOB",make_list<string> ("seq_id")("position")("repeat_name")("strand")("duplication_size")("gap_left")("gap_right"))
("DEL",make_list<string> ("seq_id")("position")("size"))
("INV",make_list<string> ("seq_id")("position")("size"))
("AMP",make_list<string> ("seq_id")("position")("size")("new_copy_number"))
//## evidence
("RA",make_list<string> ("seq_id")("position")("insert_position")("ref_base")("new_base"))
("MC",make_list<string> ("seq_id")("start")("end")("start_range")("end_range"))
("JC",make_list<string> ("side_1_seq_id")("side_1_position")("side_1_strand")("side_2_seq_id")("side_2_position")("side_2_strand")("overlap"))
("UN",make_list<string> ("seq_id")("start")("end"))
;


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
diff_entry::diff_entry(const string& type)
: _type(type)
, _id("")
, _evidence()
{
}

diff_entry::diff_entry()
: _type("")
, _id("")
, _evidence()
{
}



/*! Marshal this diff entry into an ordered list of fields.
 */
void diff_entry::marshal(field_list_t& s) {
	s.push_back(_type);
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
		if(iter != cp.end()) {
			s.push_back(iter->second);
			cp.erase(iter);
		}
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
  counted_ptr<diff_entry> added_item(diff_entry_copy);
  _entry_list.push_back(added_item);
    
  if (added_item->_id.size() == 0)
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
bool diff_entry_sort(const diff_entry_ptr& a, const diff_entry_ptr& b) {

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
  
  // sort
  _entry_list.sort(diff_entry_sort);
  
	for(diff_entry_list::iterator i=_entry_list.begin(); i!=_entry_list.end(); ++i) {
		ofs << (**i) << endl;
	}
	ofs.close();
}

// Convert GD file to GVF file
void GDtoGVF( const string &file, const string &gdfile ){
    
    // Stores the features
    vector< vector<string> > features;
    vector< vector<string> > featuresGVF;
    // Keeps track of the index of the entry associated with a particular evidence number
    map< string, int > eDict;
    
    //Read input into array
    ifstream gd( gdfile.c_str() );
    string line;
    getline( gd, line);
    
    while( !gd.eof() ){
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
    // f[0]: Type
    // f[1]: ID
    // f[2]: Evidences
    // f[3]: seq_id (ID of reference)
    // f[4]: Position
    // f[5]: Varies
    // f[6]: Varies
    // f[7]: Varies
    // f[8]: Varies
    
    for( size_t i=0; i<features.size(); i++ ){
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
        
        if( features[i][0].compare( "SNP") == 0 ){
            
            gvf[2] = "SNV";
            stringstream ss;
            ss << atoi( gvf[3].c_str() ) + 1;
            ss >> gvf[4];
            gvf[8].append("Variant_seq=").append( features[i][5] );
            
            //Look for evidence information
            vector<string> evidenceNums = split( features[i][2], "," );
            vector<string> evidence = features[ eDict[ evidenceNums[0] ] ];
            
            gvf[8].append(";Reference_seq=").append( evidence[5] );
            
            for( size_t j=0; j<evidence.size(); j++ ){
                string s = evidence[j];
                if( s.size()>8 && s.substr(0,8).compare("quality=") == 0){
                    gvf[5] = s.substr(8,s.size()-8);
                }
                if( s.size()>8 && s.substr(0,8).compare("new_cov=") == 0){
                    s = s.substr(8,s.size()-8);
                    vector<string> covs = split( s, "/" );
                    gvf[8] = gvf[8].append(";Variant_reads=").append(covs[0]).append(";Total_reads=").append(covs[1]);
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
        
        else if( features[i][0].compare("AMP") == 0 ){
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
        
        featuresGVF.push_back(gvf);
    }
    
    // Write results to file
    ofstream output( gdfile.c_str() );
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
        ss << i+1+featuresVCF.size();
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
diff_entry_list genome_diff::list(vector<string> types)
{
  diff_entry_list return_list;
  
  for (diff_entry_list::iterator itr_diff_entry = _entry_list.begin(); 
       itr_diff_entry != _entry_list.end(); itr_diff_entry++)
    {
      for (vector<string>::iterator requested_type = types.begin();
           requested_type != types.end(); requested_type++)
      {
        if((*itr_diff_entry)->_type == *requested_type)
          return_list.push_back(*itr_diff_entry);
      }
    }
  
  return return_list;
}

/*-----------------------------------------------------------------------------
 * returns mutations not used as evidence. 
 *
 * Poorly implemented right now, look into future refractor into returning
 * list<counted_ptr<diff_entry> >
 *-----------------------------------------------------------------------------*/
::list<counted_ptr<diff_entry> > genome_diff::filter_used_as_evidence(const diff_entry_list& input)
{
  ::list<counted_ptr<diff_entry> > return_list(input.begin(),input.end());
// Can't change first loop to ::list because you will get segmentation faults after 
// the .remove() method. Check if .end() is called only during initialization or 
// during each itr_outer != input.end() ?
  for (diff_entry_list::const_iterator itr_outer = input.begin();
       itr_outer != input.end(); itr_outer++) {
    string& input_id = (**itr_outer)._id;

    for (diff_entry_list::const_iterator itr_inner = _entry_list.begin();
         itr_inner != _entry_list.end(); itr_inner ++) {  
      diff_entry& test_item = **itr_inner;
      
        if (count(test_item._evidence.begin(), 
            test_item._evidence.end(),
            input_id) > 0) {
          return_list.remove(*itr_outer);
          break;
        }    
    
    }    
  }
  return return_list;
// # }
}


diff_entry genome_diff::_line_to_item(const string& line)
{
  list_t line_list = split(line, "\t");
  diff_entry item;
  item._type = shift<string>(line_list); 
  item._id = shift<string>(line_list);
  string evidence_string = shift<string>(line_list);
  item._evidence = split(evidence_string, ",");

  const list_t spec = line_specification[item._type];  

  if(spec.empty())
  {
    cerr << "Type " << item._type << "is not recognized for line # :"<< endl << line;
    assert(false);
  }


  for(size_t i = 0; i < spec.size(); i++)
  {
    string key = spec[i];
    string next = shift<string>(line_list);

    if(next.empty())
    {
      cerr << "Number of required items is less than expected for type " << item._type << " line:" << endl << line; 
      assert(false);
    }

    //if(regex_m("=",next))
    if (next.find("=") != string::npos)
    {
      cerr << "Unexpected key=value pair \'$next\' encountered for required item" << key << " in type " << item._type << " line:" << endl << line; 
      assert(false);
    }
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
  
/// Dealing with JC is inconvenient here// TODO @GRC Confirm not needed?
///##############################
// #   ### We do some extra convenience processing for junctions...
// #   if ($item->{type} eq 'JC')
// #   {
// #     foreach my $side_key ('side_1', 'side_2')
// #     {
// #       foreach my $key ('seq_id', 'position', 'strand')
// #       {
// #         $item->{"_$side_key"}->{$key} = $item->{"$side_key\_$key"};
// #       }
// #       $item->{"_$side_key"}->{type} = 'NA';
// #     }
// #   }
///###############################
  
 return item;
}

diff_entry_list genome_diff::mutation_list()
{
  diff_entry_list mut_list;

   for(diff_entry_list::iterator itr = _entry_list.begin();
       itr != _entry_list.end(); itr ++) {
     diff_entry& item = **itr;
     size_t type_length = item._type.length();
     if(type_length == 2) {
       mut_list.push_back(*itr);
     }
   }

  return mut_list;
}


/*! Return all diff_entrys within _entry_list whose _type matches one
 * of those within input's item._evidence
 */ 
::list<counted_ptr<diff_entry> > genome_diff::mutation_evidence_list(const diff_entry& item)
{
  ::list<counted_ptr<diff_entry> > return_list;
  vector<string> evidence_list = item._evidence;
  //return diff_entrys with matching evidence
  for (vector<string>::iterator itr_i = evidence_list.begin();
       itr_i != evidence_list.end(); itr_i ++) {  
    string& evidence = *itr_i;
    
    for (diff_entry_list::iterator itr_j = _entry_list.begin();
         itr_j != _entry_list.end(); itr_j ++) {  
      diff_entry& entry = **itr_j;
    
      if (entry._id == evidence) {
        return_list.push_back(*itr_j);
      }
    
    }   
  }
  return return_list;
}

// @JEB: we need to have this return a counted_ptr<diff_entry> (which can be NULL)
// we should create a typedef for counted_ptr<diff_entry> = diff_entry_ptr
diff_entry_ptr genome_diff::parent(const diff_entry& item)
{
  for(diff_entry_list::iterator itr_test_item = _entry_list.begin();
      itr_test_item != _entry_list.end(); itr_test_item ++) { 
    diff_entry& test_item = **itr_test_item;
    for(vector<string>::iterator itr = test_item._evidence.begin();
        itr != test_item._evidence.end(); itr ++) { 
      string& test_evidence_id = (*itr);
      if(test_evidence_id == item._id)      
        return counted_ptr<diff_entry>(*itr_test_item);
    }
  }
  return counted_ptr<diff_entry>(NULL);
}

bool genome_diff::mutation_unknown(diff_entry mut)
{
//sub mutation_unknown
//{
//	my ($self, $mut) = @_;
//
//	if ($mut->{type} eq 'SNP')
//	{
//		return $self->interval_un($mut->{position}, $mut->{position});
//	}
  if (mut._type == SNP) {
    return interval_un(from_string<uint32_t>(mut[POSITION]),
                       from_string<uint32_t>(mut[POSITION]));
  }
//	
//	## should be updated to new unknown that includes linkage
//	if ($mut->{type} eq 'INS')
//	{
//		return $self->interval_un($mut->{position}, $mut->{position}+1);
//	}
  //Should be updated to new unknown that includes linkage
  if (mut._type == INS) {
    return interval_un(from_string<uint32_t>(mut[POSITION]),
                       from_string<uint32_t>(mut[POSITION])) + 1;
  }
//
//	if ($mut->{type} eq 'DEL')
//	{
//
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
//		return $self->interval_un($mut->{position}, $mut->{position}+$mut->{size}-1);
//	}
//	
//	if ($mut->{type} eq 'SUB')
//	{
//		return $self->interval_un($mut->{position}, $mut->{position}+$mut->{size}-1);
//	}
  if (mut._type == SUB) {
    return interval_un(from_string<uint32_t>(mut[POSITION]),
                       from_string<uint32_t>(mut[POSITION]) - 1);
  }
//	
//	return 0;
  return false;
//}
}
//

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
  diff_entry_list un_list = list(make_list<string>(UN));

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

}//namespace bresesq



