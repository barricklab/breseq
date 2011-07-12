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
map<string, vector<string> > line_specification =make_map<string, vector<string> > 
//! seq_id and positions are already parameters in diff_entry
//## mutations
//'SNP' => ['seq_id', 'position', 'ref_seq', 'new_seq'],
("SNP",make_list<string> ("seq_id")("position")("ref_seq")) ///BUG 1 too many parameters, new_seq?
//'SUB' => ['seq_id', 'position', 'ref_seq', 'new_seq'],
("SUB",make_list<string> ("seq_id")("position")("ref_seq")("new_seq"))
//'DEL' => ['seq_id', 'position', 'size'],
("DEL",make_list<string> ("seq_id")("position")("size"))
//'INS' => ['seq_id', 'position', 'new_seq'],
("INS",make_list<string> ("seq_id")("position")("new_seq"))
//'MOB' => ['seq_id', 'position', 'repeat_name', 'strand', 'duplication_size', 'gap_left', 'gap_right'],
("MOB",make_list<string> ("seq_id")("position")("repeat_name")("strand")("duplication_size")("gap_left")("gap_right"))
//'DUP' => ['seq_id', 'position', 'size'],
("DEL",make_list<string> ("seq_id")("position")("size"))
//'INV' => ['seq_id', 'position', 'size'],
("INV",make_list<string> ("seq_id")("position")("size"))
//
//## evidence  ("")
//'RA' => ['seq_id', 'position', 'insert_position', 'ref_base', 'new_base'],
("RA",make_list<string> ("seq_id")("position")("insert_position")("ref_base")("new_base"))
//'MC' => ['seq_id', 'start', 'end'],
("MC",make_list<string> ("seq_id")("start")("end"))
//'JC' => ['side_1_seq_id', 'side_1_position', 'side_1_strand', 'side_2_seq_id', 'side_2_position', 'side_2_strand', 'overlap'],
("JC",make_list<string> ("side_1_seq_id")("side_1_position")("side_1_strand")("side_2_seq_id")("side_2_position")("side_2_strand")("overlap"))
//'UN' => ['seq_id', 'start', 'end'],
("UN",make_list<string> ("seq_id")("start")("end"))
//};
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
breseq::diff_entry::diff_entry(const string& type, const string& id, vector<string> positions)
: _type(type)
, _id(id)
, _evidence(positions) {
}

breseq::diff_entry::diff_entry()
: _type("")
, _id("")
, _evidence(make_list<string>("")) {
}


/*! Marshal this diff entry into an ordered list of fields.
 */
void breseq::diff_entry::marshal(field_list_t& s) {
	s.push_back(_type);
	s.push_back(_id);
  
  for(uint8_t i = 0; i < _evidence.size(); i++)
	  s.push_back(_evidence[i]);
	
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
/*! Constructor.
 */
breseq::genome_diff::genome_diff(const string& filename)
 : _default_filename(filename)
 , _current_id(0) 
{
 read(filename);  
}

/*! Add evidence to this genome diff.
 */
void breseq::genome_diff::add(const diff_entry& item) {
	_entry_list.push_back(item);
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


/*! Read a genome diff(.gd) from the given file to build the vector lines
 */

void breseq::genome_diff::read(const string& filename) {
  typedef vector<string> Lines;
  typedef string Line;
  
  ifstream IN(filename.c_str());
  if(!IN.good())
    cerr << "Could not open file for reading: " << filename << endl;
  
 Lines lines;
  char const line_delim = '\n';
  for(Line line; getline(IN, line, line_delim);)
    lines.push_back(line);
// #   ## read version from first line
// #   my $l = shift @lines;
  Line l = shift(lines);
// #   ($l =~ m/#=GENOME_DIFF\s+(\d+)/) or ($l =~ m/#=GENOMEDIFF\s+(\d+)/)  
// #   or $self->throw("Could not match version line in file $self->{file_name}.");
  if(!regex_m("#=GENOME_DIFF",l)
      && !regex_m("#=GENOMEDIFF",l))
    cerr << "Could not match version line in file" << endl;
  Lines version_split = split(l," ");
  assert(version_split.size() == 2);
// #   $self->{version} = $1;
  version = version_split.back();
// #   ## read header information
// #   
// #   ## read data
// #   while ($l = shift @lines)
// #   {
  while(!lines.empty())
  {
    l=shift(lines);
    ///TODO metadata 
    
    // # $self->add($self->_line_to_item($l));
    diff_entry temp = _line_to_item(l);
  }
  

  

// #   while ($l = shift @lines)
// #   {
// #     if ($l =~ m/^\s*#=(\S+)\s+(.+)/) {
// #       #metadata line key value pair
// #       push @{$self->{metadata}->{$1}}, $2;
// #       push @{$self->{metadata_lines}}, $l;
// #       
// #     } elsif ($l =~ m/^\s+#=(.+)/) {
// #     } else {
// #       $self->add($self->_line_to_item($l));
// #     }
// #   }
// #   close IN;
  IN.close();

  
  
  
  
		
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
bool breseq::diff_entry_sort(diff_entry a, diff_entry b) {

  string a_type = a._type;
  string b_type = b._type;

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
		ofs << (*i) << endl;
	}
	ofs.close();
}

// Convert GD file to GVF file
void breseq::GDtoGVF( const string &file, const string &gdfile ){
    
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
    
    for( int i=0; i<features.size(); i++ ){
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
            
            for( int j=0; j<evidence.size(); j++ ){
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
            
            for( int j=0; j<features[i].size(); j++ ){
                if( features[i][j].size()>10 && features[i][j].substr(0,10).compare("frequency=") == 0){
                    gvf[8].append(";Variant_freq=").append( features[i][j].substr(10,features[i][j].size()-10 ) );
                }
            }
            
        }
        
        else if( features[i][0].compare("SUB") == 0 ){
            //Look for evidence information
            vector<string> evidenceNums = split( features[i][2], "," );
            string s = "";
            for( int j=0; j<evidenceNums.size(); j++ ){
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
    for( int i=0; i<featuresGVF.size(); i++ ){
        for( int j=0; j<featuresGVF[i].size(); j++ ){
            output << featuresGVF[i][j] << "\t";
        }
        output << "\n";
    }
    output.close();
    
}

void breseq::VCFtoGD( const string& vcffile, const string& gdfile ){
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
    
    for( int i=0; i<featuresVCF.size(); i++ ){
        vector<string> gd(9,"");
        vector<string> ev(9,"");
        
        // SeqID
        gd[0] = "SNP";
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
    for( int i=0; i<featuresGD.size(); i++ ){
        for( int j=0; j<featuresGD[i].size(); j++ ){
            output << featuresGD[i][j] << "\t";
        }
        output << "\n";
    }
    output.close();
}
namespace breseq {

/*! Given a list of types, search and return the diff_entry's within entry_list_t whose 
 * _type parameter matches one of those within input types. 
 */ 
genome_diff::entry_list_t genome_diff::list(vector<string> types)
{
// # sub list
// # {
// #   my ($self, @types) = @_;
// #     
// #   ## return ALL
// #   if (scalar @types == 0)
// #   {
  if(types.size())
  {
// #     return @{$self->{list}};
    return _entry_list;
// #   }
  }
// #   ## return only requested types
// #   else
// #   {
  else
  {
// #     my %requested;
    
// #     foreach my $type (@types)
// #     {
// #       $requested{$type} = 1;
// #     }
// #     my @return_list = grep { $requested{$_->{type}} } @{$self->{list}};
// #     return @return_list;
// #   }
    entry_list_t return_list;
    for (entry_list_t::iterator itr_diff_entry = _entry_list.begin(); 
       itr_diff_entry != _entry_list.end() ;itr_diff_entry ++)
    {
      for (make_list<string>::iterator requested_type = types.begin();
           requested_type != types.end(); requested_type++)
      {
        if(itr_diff_entry->_type == *requested_type)
          return_list.push_back(*itr_diff_entry);
      }
    }
    return return_list;
  }
// #   
// #   return undef;
  ///TODO return undef
// # }
}

genome_diff::entry_list_t genome_diff::filter_used_as_evidence(entry_list_t list)
{
// # sub filter_used_as_evidence
// # {
// #   my ($self, @list) = @_;
// #   
// #   IN: for (my $i=0; $i<scalar @list; $i++)
// #   {
// #     my $in_item = $list[$i];
// #     foreach my $test_item ($self->list)
// #     {
// #       foreach my $test_evidence_id (@{$test_item->{evidence}})
// #       {
// #         if ($test_evidence_id ==  $in_item->{id})
// #         {
// #           splice @list, $i, 1;
// #           $i--;
// #           next IN;
// #         }
// #       }
// #     }
// #   }
// #   return @list;
return list;
// # }
}


diff_entry genome_diff::_line_to_item(string line)
{

// # sub _line_to_item
// # {
// #   my ($self, $line) = @_;
// #   my @line_list = split /\t/, $line;
  list_t line_list = split(line, "\t");
// #   
// #   ##remove items at the end that are empty 
// #   while ( scalar(@line_list) && ($line_list[-1] =~ m/^\s+$/) )
// #   {
// #     pop @line_list;
// #   }
  ///^^^NOT Applicable?
// #   
// #   my $item = {};
  diff_entry item;
// #   $item->{type} = shift @line_list;
  item._type = shift(line_list); 
// #   $item->{id} = shift @line_list;
  item._id = shift(line_list);
// #   my $evidence_string = shift @line_list;
  string evidence_string = shift(line_list);
// #   @{$item->{evidence}} = split /,/, $evidence_string;
  item._evidence = split(evidence_string, ",");
// #     
// #   my $spec = $line_specification->{$item->{type}};
  const list_t spec = line_specification[item._type];  
// #   if (!defined $spec)
// #   {
  if(spec.empty())
  {
// #     $self->warn("Type \'$item->{type}\' is not recognized for line:\n$line");
    cerr << "Type " << item._type << "is not recognized for line # :"<< endl << line;
// #     return undef;
    ///TODO undef
// #   }
  }
// #   
// #   ######## Temporary transition code for 'DUP' => AMP
// #   if ($item->{type} eq 'DUP')
// #   {
  if(item._type == "DUP")
  {
// #     $item->{type} = 'AMP';
    item._type = "AMP";
// #     $item->{new_copy_number} = 2;
    item["new_copy_number"] = 2;
// #   }
  }
// #   
// #   ######## Temporary transition code for 'MOB'
// #   if ($item->{type} eq 'MOB')
// #   {
  if(item._type == "MOB")
  {
// #     my @spec_items = grep {!($_ =~ m/=/)} @line_list;
    list_t spec_items = grep(false, "=", line_list);    
// #     
// #     if (scalar(@spec_items) == scalar(@$spec) + 2)
// #     {
    if(spec_items.size() == spec.size() + 2)
    {
// #       my $gap_left = $line_list[5];
      string gap_left = line_list[5];
// #       if ($gap_left =~ m/^-/)
// #       {
      if(regex_m("-", gap_left))
      {
// #         $item->{del_start} = abs($gap_left);
        item["del_start"] = to_string(abs(atoi(gap_left.c_str()))); ///TODO GROSS
// #       }
      }
// #       else
// #       {
      else
      {
// #         $item->{ins_start} = $gap_left;
        item["ins_start"] = gap_left;
// #       }
      }
///TODO TEST if abs(atoi("-7")) equals 7
// # 
// #       my $gap_right = $line_list[6];
      string gap_right = line_list[6];
// #       if ($gap_right =~ m/^-/)
// #       {
      if(regex_m("-",gap_right))
      {
// #         $item->{del_end} = abs($gap_left);
        item["del_end"] = abs(atoi(gap_left.c_str())); ///TODO GROSS / really gap_left?
// #       }
      }
// #       else
// #       {
      else
      {
// #         $item->{ins_end} = $gap_left;
        item["ins_end"] = gap_left;
// #       }
      }
// #       
// #       ##remove these items
// #       splice @line_list, 5, 2;
      splice(line_list, 5, 2); 
// #     }
    }
// #   }
  }
// #   
// #   

// #   foreach my $key (@$spec)
// #   {
  for(uint8_t i = 0; i < spec.size(); i++)
  {
    string key = spec[i];
// #     my $next = shift @line_list;
   string next = shift(line_list);
// #     if (!defined $next)
// #     {
    if(next.empty())
    {
// #       $self->warn("Number of required items is less than expected for type \'$item->{type}\' line:\n$line");
      cerr << "Number of required items is less than expected for type " << item._type << " line:" << endl << line; 
// #       return undef;
      ///TODO return undef
      
// #     }
    }
// #     if ($next =~ m/=/) 
// #     {
    if(regex_m("=",next))
    {
// #       $self->warn("Unexpected key=value pair \'$next\' encountered for required item \'$key\' in type \'$item->{type}\' line:\n$line");
      cerr << "Unexpected key=value pair \'$next\' encountered for required item" << key << " in type " << item._type << " line:" << endl << line; 
// # #     return undef;
      ///TODO return undef
// #     }
    }
// #     
// #     $item->{$key} = $next;
  item[key] = next;
// #   }
  }
  
// # 
// #   ## Remainder of the line is comprised of non-required key value pairs
// #   foreach my $key_value_pair (@line_list)
// #   {
  for(list_t::iterator itr = line_list.begin();
      itr != line_list.end(); itr ++)
  {
       string key_value_pair(*itr); 
// #     next if (!$key_value_pair);
       if(key_value_pair.empty()) continue;
// #     next if ($key_value_pair =~ m/^\s*$/);
       if(!regex_m("=",key_value_pair)) continue;//TODO ASK
// #     my $matched = ($key_value_pair =~ m/^(.+)=(.+)$/);
         vector<string> matched = split("=", key_value_pair);
// #     if (!$matched)
// #     {
       if(matched.empty())
       {
// #       $self->warn("Not a key value pair \'$key_value_pair\' line:\n$line");
         cerr << "Not a key value pair" << key_value_pair <<  endl << line;
// #       next;
         continue;
// #     }   
       }
// #     
// #     my ($item_key, $item_value) = ($1, $2); +
         string item_key = matched[0], item_value = matched[1]; 
         // #     $item->{$item_key} = $item_value;
         item[item_key] = item_value;
// #   }
  }
/// Dealing with JC is inconvenient here
///#############################3
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
// # 
// #   return $item;
 return item;
// # } 
}




}//namespace bresesq
