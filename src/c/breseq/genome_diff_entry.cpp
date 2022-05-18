/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2008-2010 Michigan State University
 Copyright (c) 2011-2022 The University of Texas at Austin
 
 breseq is free software; you can redistribute it and/or modify it under the
 terms the GNU General Public License as published by the Free Software
 Foundation; either version 1, or (at your option) any later version.
 
 *****************************************************************************/

#include "libbreseq/genome_diff_entry.h"
#include "libbreseq/reference_sequence.h"

namespace breseq {
  // Common keywords used for diff entries:
  
  const string cDiffEntry::unprintable_key_prefix("_");
  
  // Shared
  const char* SEQ_ID="seq_id";
  const char* START="start";
  const char* END="end";
  const char* STRAND="strand";
  const char* POSITION="position";
  const char* INSERT_POSITION="insert_position";
  const char* PHYLOGENY_ID="phylogeny_id";
  const char* FREQUENCY="frequency";
  const char* REJECT="reject";
  const char* USER_DEFINED="user_defined";
  const char* POLYMORPHISM_REJECT="polymorphism_reject";
  const char* CONSENSUS_REJECT="consensus_reject";
  const char* MEDIATED="mediated";
  const char* BETWEEN="between";
  
  // For APPLY
  const char* WITHIN="within";
  const char* BEFORE="before";
  const char* APPLY_SIZE_ADJUST = "apply_size_adjust";
  
  // For DEL
  const char* SIZE="size";
  
  // For INS
  const char* NEW_SEQ="new_seq";
  
  // For MOB
  const char* REPEAT_NAME = "repeat_name";
  const char* DUPLICATION_SIZE = "duplication_size";
  const char* INS_START = "ins_start";
  const char* INS_END = "ins_end";
  const char* DEL_START = "del_start";
  const char* DEL_END = "del_end";
  const char* MOB_REGION = "mob_region";
  
  // For INS/DEL
  const char* REPEAT_SEQ = "repeat_seq";
  const char* REPEAT_LENGTH = "repeat_length";
  const char* REPEAT_REF_COPIES = "repeat_ref_copies";
  const char* REPEAT_NEW_COPIES = "repeat_new_copies";
  
  // For AMP
  const char* NEW_COPY_NUMBER = "new_copy_number";
  const char* MEDIATED_STRAND = "mediated_strand";
  
  // For CON/INT
  const char* REGION = "region";
  const char* REPLACE_SIZE = "replace_size";
  
  //For RA
  // old + new required field
  const char* REF_BASE="ref_base";
  // old fields to maintain RA definition
  const char* NEW_BASE="new_base";
  const char* REF_COV="ref_cov";
  const char* NEW_COV="new_cov";
  
  // new fields
  const char* MAJOR_BASE="major_base";
  const char* MINOR_BASE="minor_base";
  const char* MAJOR_COV="major_cov";
  const char* MINOR_COV="minor_cov";
  const char* TOTAL_COV="total_cov";
  const char* PREDICTION = "prediction";
  const char* SCORE="score";
  const char* CONSENSUS_SCORE="consensus_score";
  const char* POLYMORPHISM_SCORE="polymorphism_score";
  const char* POLYMORPHISM_FREQUENCY="polymorphism_frequency";
  const char* MAJOR_FREQUENCY="major_frequency";            // frequency of major allele
  const char* POLYMORPHISM_EXISTS="polymorphism_exists";    // internal flag for running R script
  
  //For MC
  const char* START_RANGE="start_range";
  const char* END_RANGE="end_range";
  const char* LEFT_OUTSIDE_COV = "left_outside_cov";
  const char* LEFT_INSIDE_COV = "left_inside_cov";
  const char* RIGHT_INSIDE_COV = "right_inside_cov";
  const char* RIGHT_OUTSIDE_COV = "right_outside_cov";
  
  //For JC
  const char* SIDE_1_SEQ_ID = "side_1_seq_id";
  const char* SIDE_1_POSITION = "side_1_position";
  const char* SIDE_1_STRAND = "side_1_strand";
  const char* SIDE_1_REDUNDANT = "side_1_redundant";
  
  const char* SIDE_2_SEQ_ID = "side_2_seq_id";
  const char* SIDE_2_POSITION = "side_2_position";
  const char* SIDE_2_STRAND = "side_2_strand";
  const char* SIDE_2_REDUNDANT = "side_2_redundant";
  
  const char* OVERLAP = "overlap";
  const char* UNIQUE_READ_SEQUENCE = "unique_read_sequence";
  
  const char* SIDE_1_READ_COUNT="side_1_read_count";
  const char* SIDE_2_READ_COUNT="side_2_read_count";
  const char* NEW_JUNCTION_READ_COUNT="new_junction_read_count";
  const char* NEW_JUNCTION_FREQUENCY = "new_junction_frequency";
  
  const char* SIDE_1_COVERAGE = "side_1_coverage";
  const char* SIDE_2_COVERAGE = "side_2_coverage";
  const char* NEW_JUNCTION_COVERAGE = "new_junction_coverage";
  
  //For CN
  const char* COPY_NUMBER = "copy_number";
  
  const int32_t k_num_line_specification_common_prefix_columns = 3;
  
  const map<gd_entry_type, vector<string> > line_specification = make_map<gd_entry_type, vector<string> >
  //! seq_id and positions are already parameters in cDiffEntry
  //## mutations
  (SNP,make_vector<string> (SEQ_ID)(POSITION)(NEW_SEQ))
  (SUB,make_vector<string> (SEQ_ID)(POSITION)(SIZE)(NEW_SEQ))
  (DEL,make_vector<string> (SEQ_ID)(POSITION)(SIZE))
  (INS,make_vector<string> (SEQ_ID)(POSITION)(NEW_SEQ))
  (MOB,make_vector<string> (SEQ_ID)(POSITION)(REPEAT_NAME)(STRAND)(DUPLICATION_SIZE))
  (AMP,make_vector<string> (SEQ_ID)(POSITION)(SIZE)(NEW_COPY_NUMBER))
  (INV,make_vector<string> (SEQ_ID)(POSITION)(SIZE))
  (CON,make_vector<string> (SEQ_ID)(POSITION)(REPLACE_SIZE)(REGION))
  (INT,make_vector<string> (SEQ_ID)(POSITION)(REPLACE_SIZE)(REGION))
  
  //## evidence
  (RA,make_vector<string> (SEQ_ID)(POSITION)(INSERT_POSITION)(REF_BASE)(NEW_BASE))
  (MC,make_vector<string> (SEQ_ID)(START)(END)(START_RANGE)(END_RANGE))
  (JC,make_vector<string> (SIDE_1_SEQ_ID)(SIDE_1_POSITION)(SIDE_1_STRAND)(SIDE_2_SEQ_ID)(SIDE_2_POSITION)(SIDE_2_STRAND)(OVERLAP))
  (CN,make_vector<string> (SEQ_ID)(START)(END)(COPY_NUMBER))
  (UN,make_vector<string> (SEQ_ID)(START)(END))
  
  //## validation
  (CURA,make_vector<string> ("expert"))
  (FPOS,make_vector<string> ("expert"))
  (PHYL,make_vector<string> ("gd"))
  (TSEQ,make_vector<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
  (PFLP,make_vector<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
  (RFLP,make_vector<string> ("seq_id")("primer_1_start")("primer_1_end")("primer_2_start")("primer_2_end"))
  (PFGE,make_vector<string> ("seq_id")("enzyme"))
  (NOTE,make_vector<string> ("note"))
  (MASK,make_vector<string> (SEQ_ID)(POSITION)(SIZE))
  
  ; // end line specifications
  
  // These specs include addition fields used when determined equal mutations and sorting!
  // IMPORTANT: They include fields that may not always be defined.
  // NOTE: 'unique' gets added to ALL specs for determining uniqueness
  const map<gd_entry_type, vector<string> > extended_line_specification = make_map<gd_entry_type, vector<string> >
  (INS,make_vector<string> (SEQ_ID)(POSITION)(INSERT_POSITION)(NEW_SEQ))
  (MOB,make_vector<string> (SEQ_ID)(POSITION)(REPEAT_NAME)(STRAND)(DUPLICATION_SIZE)(INS_START)(INS_END)(DEL_START)(DEL_END)(MOB_REGION))
  (AMP,make_vector<string> (SEQ_ID)(POSITION)(SIZE)(NEW_COPY_NUMBER)(MEDIATED)(MEDIATED_STRAND)(MOB_REGION))
  (RA,make_vector<string>  (SEQ_ID)(POSITION)(INSERT_POSITION)(REF_BASE)(NEW_BASE))
  (JC,make_vector<string>  (SIDE_1_SEQ_ID)(SIDE_1_POSITION)(SIDE_1_STRAND)(SIDE_2_SEQ_ID)(SIDE_2_POSITION)(SIDE_2_STRAND)(OVERLAP)(UNIQUE_READ_SEQUENCE))
  ;
  
  enum diff_entry_field_variable_t {
    kDiffEntryFieldVariableType_BaseSequence,
    kDiffEntryFieldVariableType_PositiveInteger,
    kDiffEntryFieldVariableType_NonNegativeInteger,
    kDiffEntryFieldVariableType_Integer,
    kDiffEntryFieldVariableType_Strand, // must be -1 or +1
    kDiffEntryFieldVariableType_PositiveInteger_ReverseSort,
  };
  
  map<string, diff_entry_field_variable_t > diff_entry_field_variable_types = make_map<string, diff_entry_field_variable_t>
  (POSITION, kDiffEntryFieldVariableType_PositiveInteger)
  (START, kDiffEntryFieldVariableType_PositiveInteger)
  (END, kDiffEntryFieldVariableType_PositiveInteger)
  (SIZE, kDiffEntryFieldVariableType_PositiveInteger_ReverseSort)
  (STRAND, kDiffEntryFieldVariableType_Strand)
  (NEW_SEQ, kDiffEntryFieldVariableType_BaseSequence)
  (NEW_COPY_NUMBER, kDiffEntryFieldVariableType_PositiveInteger)
  (DUPLICATION_SIZE, kDiffEntryFieldVariableType_Integer)
  (REPLACE_SIZE, kDiffEntryFieldVariableType_NonNegativeInteger)
  (DEL_START, kDiffEntryFieldVariableType_PositiveInteger)
  (DEL_END, kDiffEntryFieldVariableType_PositiveInteger)
  (INS_START, kDiffEntryFieldVariableType_BaseSequence)
  (INS_END, kDiffEntryFieldVariableType_BaseSequence)
  (INSERT_POSITION, kDiffEntryFieldVariableType_Integer)
  (MEDIATED_STRAND, kDiffEntryFieldVariableType_Strand)
  //JC item
  (SIDE_1_POSITION, kDiffEntryFieldVariableType_PositiveInteger)
  (SIDE_1_STRAND, kDiffEntryFieldVariableType_Strand)
  (SIDE_2_POSITION, kDiffEntryFieldVariableType_PositiveInteger)
  (SIDE_2_STRAND, kDiffEntryFieldVariableType_Strand)
  (OVERLAP, kDiffEntryFieldVariableType_Integer)
  (UNIQUE_READ_SEQUENCE, kDiffEntryFieldVariableType_BaseSequence)
  ;
  
  const vector<string>gd_entry_type_lookup_table =
  make_vector<string>("UNKNOWN")("SNP")("SUB")("DEL")("INS")("MOB")("AMP")("INV")("CON")("INT")("RA")("MC")("JC")("CN")("UN")("CURA")("FPOS")("PHYL")("TSEQ")("PFLP")("RFLP")("PFGE")("NOTE")("MASK");
  
  // Used when determining what fields need to be updated if ids are renumbered
  // accounts for key=mutation_id:copy_index notation.
  const vector<string> gd_keys_with_ids = make_vector<string>("before")("within");
    
  ////
  // Begin sorting variables
  ////
  // All fields must be assigned in this table and be required fields of the gd entries.
  map<gd_entry_type, cDiffEntry::sort_fields_item> diff_entry_sort_fields = make_map<gd_entry_type, cDiffEntry::sort_fields_item>
  (SNP,  cDiffEntry::sort_fields_item(1, SEQ_ID, POSITION))
  (SUB,  cDiffEntry::sort_fields_item(1, SEQ_ID, POSITION))
  (DEL,  cDiffEntry::sort_fields_item(1, SEQ_ID, POSITION))
  (INS,  cDiffEntry::sort_fields_item(1, SEQ_ID, POSITION))
  (MOB,  cDiffEntry::sort_fields_item(1, SEQ_ID, POSITION))
  (AMP,  cDiffEntry::sort_fields_item(1, SEQ_ID, POSITION))
  (INV,  cDiffEntry::sort_fields_item(1, SEQ_ID, POSITION))
  (CON,  cDiffEntry::sort_fields_item(1, SEQ_ID, POSITION))
  (INT,  cDiffEntry::sort_fields_item(1, SEQ_ID, POSITION))
  (NOTE, cDiffEntry::sort_fields_item(2, "note", "note"))
  (RA,   cDiffEntry::sort_fields_item(3, SEQ_ID, POSITION))
  (MC,   cDiffEntry::sort_fields_item(4, SEQ_ID, START))
  (JC,   cDiffEntry::sort_fields_item(5, SIDE_1_SEQ_ID, SIDE_1_POSITION))
  (CN,   cDiffEntry::sort_fields_item(6, SEQ_ID, START))
  (UN,   cDiffEntry::sort_fields_item(7, SEQ_ID, START))
  (CURA, cDiffEntry::sort_fields_item(8, "expert", "expert"))
  (FPOS, cDiffEntry::sort_fields_item(8, "expert", "expert"))
  (PHYL, cDiffEntry::sort_fields_item(8, "gd", "gd"))
  (TSEQ, cDiffEntry::sort_fields_item(8, "seq_id", "primer_1_start"))
  (PFLP, cDiffEntry::sort_fields_item(8, "seq_id", "primer_1_start"))
  (RFLP, cDiffEntry::sort_fields_item(8, "seq_id", "primer_1_start"))
  (PFGE, cDiffEntry::sort_fields_item(8, "seq_id", "enzyme"))
  (MASK, cDiffEntry::sort_fields_item(8, SEQ_ID, POSITION))
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
  (INT, 9)
  (RA,  10)
  (MC,  11)
  (JC,  12)
  (CN,  13)
  (UN,  14)
  (CURA, 15)
  (FPOS, 16)
  (PHYL, 17)
  (TSEQ, 18)
  (PFLP, 19)
  (RFLP, 20)
  (PFGE, 21)
  (NOTE, 21)
  (MASK, 21)
  ;
  ////
  // End sorting variables
  ////
  
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
    // this is a hidden field
    de["_line_number"] = to_string<uint32_t>(line_number);
    vector<string> tokens = split(line, "\t");
    
    if (tokens.size() < k_num_line_specification_common_prefix_columns) {
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
    de._evidence = split(tokens[COLUMN], ",");
    for (uint32_t i = 0; i < de._evidence.size(); ++i) {
      RemoveLeadingTrailingWhitespace(de._evidence[i]);
    }
    ++COLUMN;
    
    //Specs.
    const vector<string>& specs = line_specification.find(de._type)->second;
    
    if (tokens.size() - COLUMN < specs.size() ) {
      if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Expected " + to_string(specs.size() + COLUMN) + " tab-delimited columns for entry", true);
      return;
    }
    
    for (uint32_t i = 0; i < specs.size(); ++i) {
      de[specs[i]] = tokens[COLUMN];
      RemoveLeadingTrailingWhitespace(de[specs[i]]);
      ++COLUMN;
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
        
        
        // Certain keys are only allowed for specific entries
        if (key == MEDIATED) {
          if ( (de._type != DEL) && (de._type != AMP) && (de._type != SUB) ) {
            if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Key 'mediated' is only allowed for entries of type DEL, AMP, or SUB.", true);
          }
        }
        
        if (key == MEDIATED_STRAND) {
          if ( (de._type != AMP) ) {
            if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Key 'mediated_strand' is only allowed for entries of type AMP.", true);
          }
        }
        
        if (key == BETWEEN) {
          if ( (de._type != DEL) && (de._type != AMP) && (de._type != CON) ) {
            if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Key 'between' is only allowed for entries of type MOB, AMP, or CON.", true);
          }
        }
        
        // Deprecated keys
        
        if (key == "nested_within") {
          if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Key 'nested_within' is DEPRECATED and will be ignored. Use 'within' instead.", false);
        } else if (key == "nested_copy") {
          if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Key 'nested_copy' IS DEPRECATED and will be ignored. Use 'within' instead.", false);
        } else if (key == "after") {
          if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Key 'after' is DEPRECATED and will be ignored. Use 'within' or 'before' instead.", false);
        }
        
      } else {
        if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Field " + kvp + " is not a key=value pair. Ignoring this key.", false);
      }
      ++COLUMN;
    }
    
    // Certain keys must occur together (if one is there, the other had better be there)
    if (de._type == AMP) {
      if (de.entry_exists(MEDIATED) != de.entry_exists(MEDIATED_STRAND)) {
        if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Only one key of 'mediated' and 'mediated_strand' is supplied for this AMP. Both must be present to describe the mutation. Did you mean to use 'between' instead?", true);
      }
    }
    
    // Certain keys are only valid with certain entries
    if ( de.entry_exists(APPLY_SIZE_ADJUST) ) {
      if ( (de._type != AMP) && (de._type != DEL) && (de._type != SUB) && (de._type != CON) && (de._type != INT)  && (de._type != INV) ) {
        if (file_parse_errors) file_parse_errors->add_line_error(line_number, line, "Key 'apply_size_adjust' is only allowed for AMP, CON, INT, DEL, INV, and SUB mutations.", true);
      }
    }
    
    // We need to implicitly treat INS with no insert_position as having 1 for that
    //   but also flag so that it will be removed when we write
    if ( (de._type == INS) && (!de.entry_exists(INSERT_POSITION))) {
      de[INSERT_POSITION] = "1";
      de["_dont_print_insert_position"] = "1";
    }
    
    return;
  }
  
  bool cDiffEntry::valid_id(string& test_id)  {
    
    // Record that we used this id, if it is currently valid
    int32_t test_id_int;
    bool is_int = is_integer(test_id, test_id_int);
    return (is_int && (test_id_int >= 1));
  }
  
  gd_entry_type cDiffEntry::type_to_enum(string type) {
    for (size_t i = 0; i < gd_entry_type_lookup_table.size(); i++) {
      if (type == gd_entry_type_lookup_table[i]) {
        return (gd_entry_type) i;
      }
    }
    return UNKNOWN;
  }
  
  // Checks all specification fields for expected types and adds errors to parse_errors
  void cDiffEntry::valid_field_variable_types(cFileParseErrors& parse_errors) {
    
    const vector<diff_entry_key_t>& spec = extended_line_specification.count(this->_type)
    ? extended_line_specification.find(this->_type)->second : line_specification.find(this->_type)->second;
    
    uint32_t field_count = k_num_line_specification_common_prefix_columns; // skipping fixed fields
    for(map<string,string>::iterator it=this->begin(); it!=this->end(); ++it) {
      field_count++;
      string key = it->first;
      if (diff_entry_field_variable_types.count(key) == 0) continue;
      
      diff_entry_field_variable_t variable_type = diff_entry_field_variable_types[key];
      string value = it->second;
      
      if (variable_type == kDiffEntryFieldVariableType_BaseSequence) {
        if (!is_base_sequence(value))
          parse_errors.add_line_error(from_string<uint32_t>((*this)["_line_number"]), this->as_string(), "Expected base sequence for field containing only characters 'ATCGN' for field #" + to_string<uint32_t>(field_count) + " (" + key +  "): [" + key + "] instead of [" + value + "]." , true);
        continue;
      }
      
      /// Other types ---> are all numerical
      int32_t ret_val;
      bool integral = is_integer(value, ret_val);
      
      if (!integral) {
        parse_errors.add_line_error(from_string<uint32_t>((*this)["_line_number"]), this->as_string(), "Expected integral value for field " + to_string<uint32_t>(field_count) + ": [" + key + "] instead of [" + value + "]." , true);
        continue;
      }
      
      
      switch(variable_type)
      {
        case kDiffEntryFieldVariableType_PositiveInteger:
        case kDiffEntryFieldVariableType_PositiveInteger_ReverseSort:
          
          if (ret_val <= 0) {
            parse_errors.add_line_error(from_string<uint32_t>((*this)["_line_number"]), this->as_string(), "Expected positive integral value for field " + to_string<uint32_t>(field_count) + ": [" + key + "] instead of [" + value + "]."  , true);
          }
          break;
          
        case kDiffEntryFieldVariableType_NonNegativeInteger:
          if (ret_val < 0) {
            parse_errors.add_line_error(from_string<uint32_t>((*this)["_line_number"]), this->as_string(), "Expected zero or positive integral value for field " + to_string<uint32_t>(field_count) + ": [" + key + "] instead of [" + value + "]."  , true);
          }
          break;
          
        case kDiffEntryFieldVariableType_Strand:
          if ((ret_val != -1) && (ret_val != 1)) {
            parse_errors.add_line_error(from_string<uint32_t>((*this)["_line_number"]), this->as_string(), "Expected strand value (-1/1) for field " + to_string<uint32_t>(field_count) + ": [" + key + "] instead of [" + value + "]." , true);
          }
          break;
          
          // already tested
        case kDiffEntryFieldVariableType_BaseSequence:
        case kDiffEntryFieldVariableType_Integer:

          break;
      }
    }
  }
  
  //Static comparison function used for comparison operators
  int32_t cDiffEntry::compare(const cDiffEntry& a, const cDiffEntry& b)
  {
    gd_entry_type a_type = a._type;
    gd_entry_type b_type = b._type;
    
    //////////////////////////////////////////////////////////////////
    // First we sort according to output order
    
    cDiffEntry::sort_fields_item a_sort_fields = diff_entry_sort_fields[a_type];
    cDiffEntry::sort_fields_item b_sort_fields = diff_entry_sort_fields[b_type];
    
    if (a_sort_fields._f1 < b_sort_fields._f1) {
      return -1;
    } else if (a_sort_fields._f1 > b_sort_fields._f1) {
      return +1;
    }
    
    string a_sort_field_2 = a.find(a_sort_fields._f2)->second;
    string b_sort_field_2 = b.find(b_sort_fields._f2)->second;
    
    if (a_sort_field_2 < b_sort_field_2) {
      return -1;
    } else if (a_sort_field_2 > b_sort_field_2) {
      return +1;
    }
    
    uint32_t a_sort_field_3 = from_string<uint32_t>(a.find(a_sort_fields._f3)->second);
    uint32_t b_sort_field_3 = from_string<uint32_t>(b.find(b_sort_fields._f3)->second);
    
    if (a_sort_field_3 < b_sort_field_3) {
      return -1;
    } else if (a_sort_field_3 > b_sort_field_3) {
      return +1;
    }
    
    // Prefer certain mutation types before others
    uint8_t a_sort_order = sort_order[a_type];
    uint8_t b_sort_order = sort_order[b_type];
    
    if (a_sort_order < b_sort_order) {
      return -1;
    } else if (a_sort_order > b_sort_order) {
      return +1;
    }
    
    //////////////////////////////////////////////////////////////////
    // Then we sort for all fields defined in the line specification
    
    // Get full line spec
    vector<diff_entry_key_t> specs = extended_line_specification.count(a_type)
    ? extended_line_specification.find(a_type)->second : line_specification.find(a_type)->second;
    
    // always add these for uniqueness testing
    specs.push_back("phylogeny_id");
    specs.push_back("unique");
    specs.push_back("population_id");
    
    for(vector<diff_entry_key_t>::const_iterator it = specs.begin(); it != specs.end(); it++) {
      const diff_entry_key_t& spec(*it);
      
      bool a_exists = a.entry_exists(spec);
      bool b_exists = b.entry_exists(spec);
      
      // look for other breaks
      if (!a_exists && !b_exists) continue;
      if (!b_exists) return +1;
      if (!a_exists) return -1;
      
      // Perform the proper type of comparison
      // Default is a string if not provided...
      
      if (!diff_entry_field_variable_types.count(spec) ||
          (diff_entry_field_variable_types[spec] == kDiffEntryFieldVariableType_BaseSequence) ) {
        
        const string& a_val = a.find(spec)->second;
        const string& b_val = b.find(spec)->second;
        
        if (a_val < b_val)
          return -1;
        else if (a_val > b_val)
          return +1;
        
      } else {
        
        switch(diff_entry_field_variable_types[spec]) {
          case kDiffEntryFieldVariableType_PositiveInteger:
          {
            uint32_t a_val = from_string<uint32_t>(a.find(spec)->second);
            uint32_t b_val = from_string<uint32_t>(b.find(spec)->second);
            
            if (a_val < b_val)
              return -1;
            else if (a_val > b_val)
              return +1;
            break;
          }
            
          case kDiffEntryFieldVariableType_PositiveInteger_ReverseSort:
          {
            uint32_t a_val = from_string<uint32_t>(a.find(spec)->second);
            uint32_t b_val = from_string<uint32_t>(b.find(spec)->second);
            
            if (a_val < b_val)
              return +1;
            else if (a_val > b_val)
              return -1;
            break;
          }
            
          case kDiffEntryFieldVariableType_Integer:
          case kDiffEntryFieldVariableType_NonNegativeInteger:
          case kDiffEntryFieldVariableType_Strand:
          {
            int32_t a_val = from_string<int32_t>(a.find(spec)->second);
            int32_t b_val = from_string<int32_t>(b.find(spec)->second);
            
            if (a_val < b_val)
              return -1;
            else if (a_val > b_val)
              return +1;
            break;
          }
            // handled above
          case kDiffEntryFieldVariableType_BaseSequence:
            break;
        }
      }
    }
    
    //////////////////////////////////////////////////////////////////
    // Return zero if the entries are equal (according to major specs)
    return 0;
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
  
  cReferenceCoordinate cDiffEntry::get_reference_coordinate_start() const
  {
    switch (this->_type) {
      case SNP:
      case SUB:
      case DEL:
      case INV:
      case AMP:
      case CON:
      case INT:
      case MASK:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)));
      case INS:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)), this->entry_exists(INSERT_POSITION) ? from_string<uint32_t>(this->at(INSERT_POSITION)) : 1);
      case MOB:
      {
        if (from_string<uint32_t>(this->at("duplication_size")) == 0) {
          return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)), 1);
        } else {
          return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)));
        }
      }
      case RA:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)), from_string<uint32_t>(this->at(INSERT_POSITION)));
      case MC:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(START)));
      case UN:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(START)));
      default:
        ERROR("cDiffEntry::get_reference_coordinate_start not implemented for type: " + gd_entry_type_lookup_table[this->_type]);
    }
    return 0;
  }
  
  cReferenceCoordinate cDiffEntry::get_reference_coordinate_end() const
  {
    switch (this->_type) {
      case SNP:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)));
      case SUB:
      case DEL:
      case INV:
      case AMP:
      case MASK:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)) + from_string<uint32_t>(this->at(SIZE)) - 1);
      case CON:
      case INT:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)) + from_string<uint32_t>(this->at(REPLACE_SIZE)) - 1);
      case INS:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)), this->entry_exists(INSERT_POSITION) ? from_string<uint32_t>(this->at(INSERT_POSITION)) : 1);
      case MOB:
      {
        if (from_string<int32_t>(this->at("duplication_size")) == 0) {
          return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)) + from_string<int32_t>(this->at("duplication_size")), 1);
        } else {
          return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)) + abs(from_string<int32_t>(this->at("duplication_size"))) - 1);
        }
      }
      case RA:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(POSITION)), from_string<uint32_t>(this->at(INSERT_POSITION)));
      case MC:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(END)));
      case UN:
        return cReferenceCoordinate(from_string<uint32_t>(this->at(END)));
        
      default:
        ERROR("cDiffEntry::get_reference_coordinate_end not implemented for type: " + gd_entry_type_lookup_table[this->_type]);
    }
    return 0;
  }
  
  bool cDiffEntry::located_within(const cDiffEntry &within) const
  {
    return(
           (this->get(SEQ_ID) == within.get(SEQ_ID))
           && (this->get_reference_coordinate_start() >= within.get_reference_coordinate_start())
           && (this->get_reference_coordinate_end() <= within.get_reference_coordinate_end())
           );
  }
  
  int32_t cDiffEntry::mutation_size_change(cReferenceSequences& ref_seq_info) const
  {
    int32_t size_change(0);
    
    switch (this->_type) {
        
      case SNP:
        size_change = 0;
        break;
        
      case SUB:
        size_change = - from_string<uint32_t>(this->get(SIZE));
        if (this->entry_exists(APPLY_SIZE_ADJUST)) {
          size_change -= from_string<int32_t>(this->get(APPLY_SIZE_ADJUST));
          ASSERT(size_change < 0, "Mutation has zero or negative size after adding APPLY_SIZE_ADJUST:\n" + this->as_string() );
        }
        size_change += this->get(NEW_SEQ).length();
        break;
        
      case INS:
        size_change = this->get(NEW_SEQ).length();
        break;
        
      case DEL:
        size_change = -(from_string<uint32_t>(this->get(SIZE)));
        if (this->entry_exists(APPLY_SIZE_ADJUST)) {
          size_change -= from_string<int32_t>(this->get(APPLY_SIZE_ADJUST));
          ASSERT(size_change < 0, "Mutation has zero or negative size after adding APPLY_SIZE_ADJUST:\n" + this->as_string() );
        }
        break;
        
      case INV:
        size_change = 0;
        break;
        
      case MASK:
        size_change = 0;
        break;
        
      case AMP:
      {
        size_change = from_string<uint32_t>(this->get(SIZE)) * (from_string<uint32_t>(this->get("new_copy_number")) - 1);
        
        // Special case of mediated AMP
        if (this->entry_exists(MEDIATED)) {
          
          cFeatureLocation repeat_feature_picked;
          string seq_id_picked;
          string mediated_string;
          string mob_region;
          
          if (this->entry_exists(MOB_REGION)) {
            mob_region = this->get(MOB_REGION);
          }
          
          mediated_string = ref_seq_info.repeat_family_sequence(this->get(MEDIATED), from_string<int16_t>(this->get(MEDIATED_STRAND)), mob_region.size() ? &mob_region : NULL, &seq_id_picked, &repeat_feature_picked);
          
          size_change += mediated_string.size() * (from_string<uint32_t>(this->get("new_copy_number")) - 1);
        }
        
        if (this->entry_exists(APPLY_SIZE_ADJUST)) {
          size_change += from_string<int32_t>(this->get(APPLY_SIZE_ADJUST)) * (from_string<uint32_t>(this->get("new_copy_number")) - 1);
          ASSERT(size_change > 0, "Mutation has zero or negative size after adding APPLY_SIZE_ADJUST:\n" + this->as_string() );
        }
        break;
      }
        
      case MOB:
      {
        // @JEB: Important: repeat_size is not a normal attribute and must be set before calling this function
        //       Also: this size includes the target site duplication
        ASSERT(this->entry_exists("repeat_size"), "Repeat size field does not exist for entry:\n" + this->as_string());
        size_change = from_string<int32_t>(this->get("repeat_size")) + from_string<int32_t>(this->get("duplication_size"));
        if (this->entry_exists("del_start"))
          size_change -= from_string<uint32_t>(this->get("del_start"));
        if (this->entry_exists("del_end"))
          size_change -= from_string<uint32_t>(this->get("del_end"));
        if (this->entry_exists("ins_start"))
          size_change += this->get("ins_start").length();
        if (this->entry_exists("ins_end"))
          size_change += this->get("ins_end").length();
        break;
      }
      
      case CON:
      case INT:
      {
        uint32_t replace_target_id, replace_start, replace_end;
        ref_seq_info.parse_region(this->get("region"), replace_target_id, replace_start, replace_end);
        size_change = from_string<uint32_t>(this->get(REPLACE_SIZE));
        
        if (this->entry_exists(APPLY_SIZE_ADJUST)) {
          size_change += from_string<int32_t>(this->get(APPLY_SIZE_ADJUST));
          ASSERT(size_change > 0, "Mutation has zero or negative size after adding APPLY_SIZE_ADJUST:\n" + this->as_string() );
        }
        
        int32_t new_size = (replace_end > replace_start) ? replace_end - replace_start + 1 : replace_start - replace_end + 1;
        size_change =  new_size - size_change;
        break;
      }
        
      default:
        ASSERT(false, "Unable to calculate mutation size change for this type of entry:\n" + this->as_string());
        return UNDEFINED_INT32;
    }
    
    return size_change;
  }
  
  // shift_offset normally means the position of the mutation doing the shifting
  // shift_offset of -1 means we are within the current mutation
  //   => we don't change its size, but we may shift its position in a special way for AMP/MOB/INS
  //   inset_pos is the index after the current position // = 0 for actually at this position
  //
  // shift_size is the change in size within the specified interval
  //
  // shift_replace size is the size of the new sequence being inserted
  //
  //   For an INS mutation =>
  //      shift_offset = position
  //      insert_pos = 1 (meaning after the base referred to in shift_offset)
  //      shift_size = number of bases inserted
  //
  
  void cDiffEntry::mutation_shift_position(const string& seq_id, const cReferenceCoordinate& shift_start, const cReferenceCoordinate& shift_end, int32_t shift_size)
  {
    
    ASSERT(this->is_mutation(), "Attempt to shift position of non-mutation cDiffEntry");
    
    // negative shift_size means deletion; positive shift_size means insertion
    if (shift_size == 0) return;
    if (seq_id != (*this)[SEQ_ID]) return;
    
    // anything that has a 'size' potentially needs to be adjusted if the shift position and size overlaps
    cReferenceCoordinate original_start = this->get_reference_coordinate_start();
    cReferenceCoordinate original_end = this->get_reference_coordinate_end();
    cReferenceCoordinate final_start = original_start;
    cReferenceCoordinate final_size = original_end - original_start + cReferenceCoordinate(1);
    
    // Special case where the mutation is 'within' this mutation
    if (shift_start.get_position() < 0) {
      
      if ((shift_start <= original_start) && (shift_end >= original_end)) {
        // change overlaps both sides
        // ERROR("Deletion completely encompasses other mutation");
        final_start = shift_start;
        final_size = 0;
      } else if ((original_start >= shift_start) && (original_start <= shift_end)) {
        // change overlaps left side
        final_start = shift_start;
        final_size = original_end - original_start + 1 - (shift_end - original_start);
      } else if ((original_end >= shift_start) && (original_end <= shift_end)) {
        // change overlaps right side
        final_start = original_start;
        final_size = shift_start - original_start;
      } else if ((shift_start > original_start) && (shift_end < original_end)) {
        // change is contained within
        final_start = original_start;
        final_size = original_end - original_start + 1 + shift_size; // original size minus size of deletion
      } else if (original_start >= shift_start) {
        final_start = original_start + shift_size;
      }
      
      // Normal case where we shift the mutation
    } else {
      
      if ((shift_start >= original_start) && (shift_end <= original_end)) {
        final_size = original_end - original_start + 1 + shift_size;
      } else if (original_start >= shift_start) {
        final_start = original_start + shift_size;
      }
    }
    
    (*this)[POSITION] = to_string(final_start.get_position());
    if (this->entry_exists(SIZE)) {
      (*this)[SIZE] = to_string(final_size.get_position());
    }
  }
  
  
  // Reverse-complements without changing position
  void cDiffEntry::mutation_reverse_complement() {
    
    switch (this->_type) {
        
        // Nothing special for these cases
      case DEL:
      case INV:
      case AMP:
        break;
        
      case SNP:
      case SUB:
      case INS:
        (*this)[NEW_SEQ] = reverse_complement((*this)[NEW_SEQ]);
        break;
        
      case CON:
      case INT:
      {
        // flip coordinates of region
        string seq_id;
        uint32_t start_pos_1, end_pos_1;
        cReferenceSequences::parse_region((*this)[REGION], seq_id, start_pos_1, end_pos_1);
        (*this)[REGION] = seq_id + ":" + s(end_pos_1) + "-" + s(start_pos_1);
        break;
      }
      case MOB:
      {
        (*this)[STRAND] = s(-n((*this)[STRAND]));
        string del_start = this->entry_exists(DEL_START) ? (*this)[DEL_START] : "";
        string del_end = this->entry_exists(DEL_END) ? (*this)[DEL_END] : "";
        string ins_start = this->entry_exists(INS_START) ? reverse_complement((*this)[INS_START]) : "";
        string ins_end = this->entry_exists(INS_END) ? reverse_complement((*this)[INS_END]) : "";
        this->erase(DEL_START);
        this->erase(DEL_END);
        this->erase(INS_START);
        this->erase(INS_END);
        
        // flips start and end
        if (!del_end.empty()) (*this)[DEL_START] = del_end;
        if (!del_start.empty()) (*this)[DEL_END] = del_start;
        if (!ins_end.empty()) (*this)[INS_START] = ins_end;
        if (!ins_start.empty()) (*this)[INS_END] = ins_start;
        break;
      }
        
      default:
        ERROR("Attempt to reverse complement cDiffEntry of unsupported type:" + to_string(this->_type));
        
    }
  }
  
  // Updates positions and reverse-complements
  void cDiffEntry::mutation_invert_position_sequence(cDiffEntry& inverting_mutation) {
    ASSERT(this->is_mutation(), "Attempt to invert position of non-mutation cDiffEntry");
    
    // are we on the right sequence fragment?
    if (inverting_mutation[SEQ_ID] != (*this)[SEQ_ID]) return;
    
    int32_t start_inversion = from_string<int32_t>(inverting_mutation[POSITION]);
    int32_t end_inversion = start_inversion + from_string<int32_t>(inverting_mutation[SIZE]) - 1;
    
    int32_t position = from_string<int32_t>((*this)[POSITION]);
    
    // Flip things that are totally contained
    // Error if something overlaps the edges
    // shift_replace_size has the size of the inversion (shift_size = 0)
    
    // Size should be how many reference bases the mutation affects!
    // It will be one for a SNP, etc., insertion, but we won't use it.
    int32_t ref_size = this->get_reference_coordinate_end().get_position() - this->get_reference_coordinate_start().get_position() + 1;
    
    // Does not overlap
    if (position + ref_size < start_inversion)
      return;
    else if (position > end_inversion)
      return;
    
    // Contained within, invert coordinates
    else if ((position >= start_inversion) && (position + ref_size - 1 <= end_inversion)) {
      
      //Reverse complement the mutation
      this->mutation_reverse_complement();
      
      if (this->_type==INS) {
        
        int32_t insert_pos;
        if (this->entry_exists(INSERT_POSITION))
          insert_pos = from_string<int32_t>((*this)[INSERT_POSITION]);
        
        // case of INS, because it refers to between positions, we need to substract an extra one from the position
        // We also need to flip insert_pos order for any mutations that are strung together with multiple insert counts
        // (we do this by taking the negative). We must re-sort after applying an INV for this reason!!
        
        (*this)[POSITION] = to_string<int32_t>(end_inversion - (position - start_inversion) - 1);
        if (this->entry_exists(INSERT_POSITION))
          (*this)[INSERT_POSITION] = to_string<int32_t>(-insert_pos);
        
      } else { // other cases
        
        (*this)[POSITION] = to_string<int32_t>(end_inversion - (position + (ref_size - 1) - start_inversion));
        
      }
      return;
    }
    
    // Overlaps end
    else {
      WARN("This mutation:\n" + this->as_string() + "\nextends across an endpoint of inversion:\n" + inverting_mutation.as_string() + "\nWhen applying the inversion, its sequence will not be reverse complemented, and its coordinates will not be shifted.");
    }
    
  }
  
  
  /* Used for marking 'adjacent' mutations in APPLY/NORMALIZE and could be used for COUNT if merged into that code
   */
  
  bool cDiffEntry::is_small_mutation(uint32_t large_size_cutoff)
  {
    cDiffEntry& mut = *this;
    
    if (mut._type == SNP) {
      return true;
    }
    
    if (mut._type == DEL) {
      if (from_string<uint32_t>(mut[SIZE]) <= large_size_cutoff) {
        return true;
      }
      
      if (mut.entry_exists(REPEAT_NEW_COPIES)) {
        return true;
      }
    }
    
    if (mut._type == INS) {
      
      if (mut[NEW_SEQ].size() <= large_size_cutoff) {
        return true;
      }
      
      if (mut.entry_exists(REPEAT_NEW_COPIES)) {
        return true;
      }
    }
    
    if (mut._type == SUB) {
      int32_t old_size = from_string<int32_t>(mut[SIZE]);
      int32_t new_size = mut[NEW_SEQ].size();
      
      if (static_cast<uint32_t>(abs(new_size - old_size)) <= large_size_cutoff) {
        return true;
      }
    }
    
    return false;
  }
  
  /*! Add 'mediate', 'adjacent', and between tags to a mutation
   */
  
  void cDiffEntry::annotate_repeat_hotspots(cReferenceSequences& new_ref_seq_info, int32_t slop_distance, int32_t size_cutoff_AMP_becomes_INS_DEL_mutation, bool remove_old_tags, bool warn_after_mode)
  {
    cDiffEntry& mut = *this;
    
    if (remove_old_tags) {
      mut.erase("adjacent");
      mut.erase("between");
      mut.erase("mediated");
    }
    
    
    map<string,string> nearby_tags;
    if (mut.entry_exists("adjacent")) nearby_tags["adjacent"] = mut["adjacent"];
    if (mut.entry_exists("between")) nearby_tags["between"] = mut["between"];
    if (mut.entry_exists("mediated")) nearby_tags["mediated"] = mut["mediated"];
    
    cAnnotatedSequence& this_seq = new_ref_seq_info[mut[SEQ_ID]];
    
    int64_t mut_start_1 = mut.get_reference_coordinate_start().get_position();
    int64_t mut_end_1 = mut.get_reference_coordinate_end().get_position();
    
    string both_close_key = "ignore";
    string one_close_key = "ignore";
    
    // Only these types can be 'within' and 'between' and have a SIZE attribute
    if ((mut._type == DEL) || (mut._type == AMP)) {
      int64_t size = from_string<int64_t>(mut[SIZE]);
      
      // short ones aren't mediated, just adjacent
      if (size > size_cutoff_AMP_becomes_INS_DEL_mutation) {
        both_close_key = "between";
        one_close_key = "mediated";
        
        if (mut._type == AMP) {
          one_close_key = "ignore"; // this ignores it
        }
      }
    }
    
    if (this->is_small_mutation()) {
      both_close_key = "adjacent";
      one_close_key = "adjacent";
    }
    
    if (mut._type == MOB) {
      one_close_key = "ignore"; // this ignores it
      both_close_key = "ignore";
    }
    
    // We make no assumptions about the directions of relevant IS elements in between/mediated here.
    int32_t tmp_slop_distance = slop_distance;
    cFeatureLocation* start_repeat = cReferenceSequences::find_closest_repeat_region_boundary(mut_start_1, this_seq.m_repeats, tmp_slop_distance, -1, true);
    if (start_repeat == NULL) {
      tmp_slop_distance = slop_distance;
      start_repeat = cReferenceSequences::find_closest_repeat_region_boundary(mut_start_1, this_seq.m_repeats, tmp_slop_distance, 1, true);
    }
    
    tmp_slop_distance = slop_distance;
    cFeatureLocation* end_repeat = cReferenceSequences::find_closest_repeat_region_boundary(mut_end_1, this_seq.m_repeats, tmp_slop_distance, 1, true);
    if (end_repeat == NULL) {
      tmp_slop_distance = slop_distance;
      end_repeat = cReferenceSequences::find_closest_repeat_region_boundary(mut_end_1, this_seq.m_repeats, tmp_slop_distance, -1, true);
    }
    
    if ((start_repeat != NULL) && (end_repeat != NULL)) {
      if ( (*(start_repeat->get_feature()))["name"] != (*(end_repeat->get_feature()))["name"]) {
        // different names is an odd case - WARN and don't assign anything
        WARN("Mutation has boundaries near two different repeat families, saving only the first one:\n" + mut.as_string());
        nearby_tags[one_close_key] = (*(start_repeat->get_feature()))["name"];
      } else {
        nearby_tags[both_close_key] = (*(start_repeat->get_feature()))["name"];
      }
      
    } else if (start_repeat != NULL) {
      nearby_tags[one_close_key] = (*(start_repeat->get_feature()))["name"];
    } else if (end_repeat != NULL) {
      nearby_tags[one_close_key] = (*(end_repeat->get_feature()))["name"];
    }
    
    if (!warn_after_mode) {
      // Finally reassign (this assumes none are removed)
      if (nearby_tags.count("adjacent")) mut["adjacent"] = nearby_tags["adjacent"];
      if (nearby_tags.count("between")) mut["between"] = nearby_tags["between"];
      if (nearby_tags.count("mediated")) mut["mediated"] = nearby_tags["mediated"];
    } else {
      // In this mode we just show a warning
      
      // Be careful to not set these accidentally to NULL by looking them up
      // when they don't exist
      bool was_adjacent = mut.entry_exists("adjacent");
      bool was_between = mut.entry_exists("between");
      bool was_mediated = mut.entry_exists("mediated");
      
      if (mut["adjacent"] != nearby_tags["adjacent"]) {
        WARN("Possible 'adjacent' tag should be added with value (" + nearby_tags["adjacent"] + ") due to later mutation:\n" + mut.as_string());
      }
      
      // Why these checks? Deletion that is between two IS looks like mediated after it occurs...
      if (!mut.entry_exists("between") && !mut.entry_exists("mediated")) {
        
        if ((mut["between"] != nearby_tags["between"])) {
          WARN("Possible 'between' tag should be added due to later mutation:\n" + mut.as_string());
        }
        if (mut["mediated"] != nearby_tags["mediated"]) {
          WARN("Possible 'mediated' tag should be added due to later mutation:\n" + mut.as_string());
        }
        
      }
      
      if (!was_adjacent) mut.erase("adjacent");
      if (!was_between) mut.erase("between");
      if (!was_mediated) mut.erase("mediated");
    }
  }
  
  
  /*! Marshal this diff entry into an ordered list of fields.
   */
  void cDiffEntry::marshal(vector<string>& s, bool include_unprintable_fields) const {
    s.push_back(gd_entry_type_lookup_table[_type]);
    s.push_back(_id);
    
    // Use a dot "." if no evidence is provided
    string evidence_string = join(_evidence, ",");
    s.push_back(evidence_string.size() ? evidence_string : ".");
    
    // deep copy all fields:
    cDiffEntry cp= *this;
    
    // marshal specified fields in-order, removing them from the copy after they've
    // been printed:
    
    const vector<string>& f = line_specification.find(_type)->second;
    
    for (vector<string>::const_iterator it=f.begin(); it != f.end(); it++)
    {
      diff_entry_map_t::iterator iter=cp.find(*it);
      
      ASSERT(iter != cp.end(), "Did not find required field '" + *it + "' to write in entry id " + _id + " of type '" + gd_entry_type_lookup_table[_type] + "'.");
      
      s.push_back(iter->second);
      cp.erase(iter);
      
    }
    
    // Remove insert_position is it was added
    if ( (cp._type == INS) && (cp.entry_exists("_dont_print_insert_position" )) ) {
      cp.erase(INSERT_POSITION);
    }
    
    // marshal whatever's left, unless it _begins with an underscore or is empty (a placeholder)
    for(diff_entry_map_t::iterator i=cp.begin(); i!=cp.end(); ++i) {
      
      assert(i->first.size());
      if (!include_unprintable_fields && is_unprintable_key(i->first)) continue;
      if (i->second.empty()) continue;
      
      // Be sure the entry is non-empty! Would rather have this as a check.
      /*
       ASSERT(i->second.size(), "Attempt to write genome diff entry with blank value in key=value pair where [key]=" + i->first + "\n" + join(s,"\t"));
       */
      
      s.push_back(i->first + "=" + i->second);
    }
  }
  
  // Created the line to be printed
  string cDiffEntry::as_string(bool include_unprintable_fields) const
  {
    vector<string> fields;
    marshal(fields, include_unprintable_fields);
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
  
  vector<string> cDiffEntry::get_reject_reasons(const string& field)
  {
    vector<string> return_value;
    if (this->entry_exists(field)) {
      return_value = split((*this)[field], ",");
    }
    return return_value;
  }
  
  /*! Add reject reason to diff entry.
   */
  void cDiffEntry::add_reject_reason(const string &reason, const string& field) {
    
    if (!this->entry_exists(field)) {
      (*this)[field] = reason;
    }
    // exists already, make comma separated list
    else {
      vector<string> current_reject_reasons = split((*this)[field], ",");
      bool reject_exists = false;
      for (vector<string>::iterator it=current_reject_reasons.begin(); it != current_reject_reasons.end(); it++) {
        if (*it == reason) reject_exists = true;
      }
      if (!reject_exists) {
        current_reject_reasons.push_back(reason);
        (*this)[field] = join(current_reject_reasons, ",");
      }
    }
  }
  
  /*! Add reject reason to diff entry.
   */
  void cDiffEntry::remove_reject_reason(const string &reason, const string& field) {
    
    if (!this->entry_exists(field)) return;
    
    vector<string> current_reject_reasons = split((*this)[field], ",");
    vector<string> new_reject_reasons;
    
    for (vector<string>::iterator it=current_reject_reasons.begin(); it != current_reject_reasons.end(); it++) {
      if (*it != reason) {
        new_reject_reasons.push_back(*it);
      }
    }
    
    if (current_reject_reasons.size() > 0) {
      (*this)[field] = join(current_reject_reasons, ",");
    }
  }
  
  
  //! Remove all information except required fields
  cDiffEntry cDiffEntry::to_spec() const
  {
    cDiffEntry de(_type);
    
    // Get full line spec
    const vector<diff_entry_key_t>& specs = extended_line_specification.count(_type)
    ? extended_line_specification.find(_type)->second : line_specification.find(_type)->second;
    
    for(vector<diff_entry_key_t>::const_iterator it = specs.begin(); it != specs.end(); it++) {
      const diff_entry_key_t& spec(*it);
      
      if (this->count(spec)) {
        de[spec] = this->get(spec);
      }
    }
    
    return de;
  }
  
  void cDiffEntry::strip_to_spec()
  {
    cDiffEntry de(_type);
    de = this->to_spec();
    *this = de;
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
          sequence.get_sequence_1_start_size(seq1_pos_1, n);
          assert(seq1.size() == n);
          
          const string &seq2 =
          sequence.get_sequence_1_start_size(seq2_pos_1, n);
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
            
            cerr << "[Mutation]: " << this->to_spec().as_string() << endl;
            cerr << "Sequence 1 [" << seq1_pos_1 << '-' << seq1_end_1 << "]: " << seq1 << endl;
            
            cerr << "Sequence 2 [" << seq2_pos_1 << '-' << seq2_end_1 << "]: " << seq2 << endl;
            
            cerr << "Deleted Sequence [" << norm_pos_1 << '-' << norm_pos_1 + n - 1 << "]: ";
            cerr << sequence.get_sequence_1_start_size(norm_pos_1, n) << endl;
            
            string prev_seq = sequence.get_sequence_1_start_size(norm_pos_1 - 10, 20);
            string left_seq = sequence.get_sequence_1_start_size(norm_pos_1 - 10, 10);
            string right_seq = sequence.get_sequence_1_start_size(norm_pos_1 + n , 10);
            
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
            const string second = sequence.get_sequence_1_start_size(i + 1, n);
            assert(second.size());
            
            if (first != second)  {
              bAmp = false;
              break;  }
          }
          
          for(;;i += 1)
          {
            const string second = sequence.get_sequence_1_start_size(i + 1, n);
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
        if ((sequence.get_sequence_1_start_size(i - (n - 1), n) == first)  && (n > 1)) {
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
        const base_char second = sequence.get_sequence_1_start_size(pos_1, 1)[0];
        
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
        
        const string second = sequence.get_sequence_1_start_size(pos_1, n);
        const string third = sequence.get_sequence_1_start_size(pos_1 + n, 1);
        
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
      case INT:
      {
        assert(this->entry_exists(REPLACE_SIZE));
        assert(this->entry_exists(REGION));
      } break;
        
        
      default:
        break;
        
    }//End switch.
    
    return;
  }
  
}
