/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2011-20122010 Michigan State University
 
 breseq is free software; you can redistribute it and/or modify it under the  
 terms the GNU General Public License as published by the Free Software 
 Foundation; either version 1, or (at your option) any later version.
 
 *****************************************************************************/

#include "libbreseq/reference_sequence.h"

#include "libbreseq/error_count.h"

using namespace std;

namespace breseq {
  
  const vector<string> snp_types = make_vector<string>("nonsynonymous")("synonymous")("noncoding")("pseudogene")("intergenic");
  
  const string BULL_DUMMY_SEQ_ID = "__BULL_DUMMY_SEQ_ID__";
  
  // Replace Sequence with Input
  // Every position between and including start_1 and end_1 will be replaced with replacement_seq.
  // This function will shift everything else.
  // mut_type is used to append the type of mutation to the feature notes.
  // verbose outputs messages to console.
  // 
  // Successfully checks all three feature lists.  m_features, m_genes, and m_repeats.
  void cAnnotatedSequence::replace_sequence_1(int32_t start_1, int32_t end_1, const string &replacement_seq, string mut_type, bool verbose)
  {
    ASSERT(start_1 <= end_1, "start (" + to_string(start_1) + ") not less than or equal to end (" + to_string(end_1) + ")");
    m_fasta_sequence.m_sequence.replace(start_1-1, end_1-start_1+1, replacement_seq);
    
    //Temporary variable for the amount to shift the start and end positions
    // Example:
    //   If we're replacing positions 5 and 6, we need to potentially shift  
    //   position 7 down by 2 even though the difference between 5 and 6 is 1.
    //   This is why we add 1.  We also need to take the new replacement
    //   length into account.
    int32_t shift = (end_1-start_1 + 1) - replacement_seq.length();
    
    //Modify the length of the sequence
    m_length -= shift;
    
    //Notify what mutation is being worked on
    if(verbose)cout << "** " << mut_type << " **" << endl;
    
    // We will be erasing inside the it loop.  This is to keep
    // track of whether or not we should iterate to the next element.
    bool it_iterate = true;
    
    //Iterate through all the features
    for (list<cSequenceFeaturePtr>::iterator it = m_features.begin(); it != m_features.end(); )
    {
      it_iterate = true;
      
      //The current feature we're looking at
      cSequenceFeature& feat = **it;
      
      //Does the feature start and end inside of the replacement?
      if(feat.m_location.get_start_1() >= start_1 && feat.m_location.get_end_1() <= end_1)
      {
        //We'll also be checking other lists to see if they
        //contain a copy of the current feature
        list<cSequenceFeaturePtr>::iterator gene_it;
        list<cSequenceFeaturePtr>::iterator repeat_it;
        
        //Is the feature in any extra lists?
        //If so, annihilate it.
        if(find_feature(m_genes, feat, gene_it)){m_genes.erase(gene_it);}
        if(find_feature(m_repeats, feat, repeat_it)){m_repeats.erase(repeat_it);}
        
        //Notify the user of the action
        if(verbose){cout << "REMOVED\t" << feat["type"] << "\t" << feat.m_gff_attributes["ID"] << " " << feat.m_gff_attributes["Name"] << endl;}
        
        //Remove the current feature
        it = m_features.erase(it);
        
        //We just removed the current feauture, do not iterate.
        it_iterate = false;
      }
      
      //Does the feature end after the replacement starts?
      else if(feat.m_location.get_end_1() > start_1 && feat.m_location.get_end_1() < end_1)
      {
        //Temporary variable for the new end position
        uint32_t end_temp = start_1 - 1;
        
        //Notify the user of the action
        if(verbose){cout << "MODIFY\t" << feat["type"]<< "\t" << feat.m_gff_attributes["ID"] << " " << feat.m_gff_attributes["Name"] << endl;}
        
        //Modify the end of the feature
        feat.m_location.set_end_1(end_temp);
        
        //Mark it as pseudo
        feat.flag_pseudo(verbose);
        
        //Modify the notes for this feature
        //If this feature is already pseudo or a region, do nothing.
        //if(feat["type"] != "region" && feat["type"] != "source" && !feat.m_gff_attributes.count("Mutation from " + mut_type))feat.m_gff_attributes["Note"].push_back("Mutation from " + mut_type);
      }
      
      //Everything that starts after the replacement starts needs to be shifted          
      else if(feat.m_location.get_start_1() > start_1)
      {
        //Does the feature start before the replacement ends?
        if(feat.m_location.get_start_1() < end_1)
        {
          //Temporary variable for the new start postion
          uint32_t start_temp = end_1 + 1;
          
          //Notify the user of the action
          if(verbose){cout << "MODIFY\t" << feat["type"] << "\t" << feat.m_gff_attributes["ID"] << " " << feat.m_gff_attributes["Name"] << endl;}
          
          //Modify the start of the feature
          feat.m_location.set_start_1(start_temp);
          
          //Mark it as pseudo
          feat.flag_pseudo(verbose);
          
          //Modify the notes for this feature
          //if(feat["type"] != "region" && feat["type"] != "source" && !feat.m_gff_attributes.count("Mutation from " + mut_type))feat.m_gff_attributes["Note"].push_back("Mutation from " + mut_type);
        }             
        
        //Is there any reason to shift?
        if(shift)
        {          
          //Modify the both the start and end of the feature
          feat.m_location.set_start_1(feat.m_location.get_start_1() - shift);
          feat.m_location.set_end_1(feat.m_location.get_end_1() - shift);
          
          //Notify the user of the action
          if(verbose){cout << "SHIFT\t" << feat["type"] << "\t" << feat.m_gff_attributes["ID"] << " " << feat.m_gff_attributes["Name"] << endl;}
        }
      }
      
      //Any feature the encompases the replaced sequence needs to be resized
      else if(feat.m_location.get_start_1() <= start_1 && feat.m_location.get_end_1() >= end_1)
      {                          
        //Is there anything to modify?
        if(shift)
        {
          //Temporary variable for the new end position
          uint32_t end_temp = feat.m_location.get_end_1() - shift;
          
          //Notify the user of the action
          if(verbose){cout << "MODIFY\t" << feat["type"] << "\t" << feat.m_gff_attributes["ID"] << " " << feat.m_gff_attributes["Name"] << endl;}
          
          //Modify the just the end of the feature
          feat.m_location.set_end_1(end_temp);
          
          //Mark it as pseudo
          feat.flag_pseudo(verbose);
          
          //Modify the notes for this feature
          //if(feat["type"] != "region" && feat["type"] != "source" && !feat.m_gff_attributes.count("Mutation from " + mut_type))feat.m_gff_attributes["Note"].push_back("Mutation from " + mut_type);
        }
      }
      
      // Iterate it ONLY if we haven't erased something.
      if(it_iterate)it++;
    }
  }
  
  // Inserts AFTER the input position
  // Place insertion_seq at first position after pos_1.
  // This function will shift everything else.
  // mut_type is used to append the type of mutation to the feature notes.
  // verbose outputs messages to console.
  // 
  // Successfully checks all three feature lists.  m_features, m_genes, and m_repeats.
  void cAnnotatedSequence::insert_sequence_1(int32_t pos_1, const string &insertion_seq, string mut_type, bool verbose)
  {
    m_fasta_sequence.m_sequence.insert(pos_1, insertion_seq);
    
    //Variable for insertion length, only want to call the
    //function once
    uint32_t insert_length = insertion_seq.length();
    
    //Modify the length of the sequence
    m_length += insert_length;
    
    //Notify what mutation is being worked on
    if(verbose)cout << "** " << mut_type << " **" << endl;
    
    //Iterate through all the features
    for (list<cSequenceFeaturePtr>::iterator it = m_features.begin(); it != m_features.end(); it++)
    {
      //The current feature we're looking at
      cSequenceFeature& feat = **it;
         
      //Does the feature end after the insertion?
      //If it ends on the same postion, do nothing
      //because we're adding it AFTER.
      if(feat.m_location.get_end_1() > pos_1)
      {
        //Does the feature start after the insertion?
        //Starting on the same postion will mean we
        //do NOT alter the start postion.
        if(feat.m_location.get_start_1() > pos_1)
        {                
          //Notify the user of the upcomming action
          if(verbose){cout << "SHIFT\t" << feat["type"] << "\t" << feat.m_gff_attributes["ID"] << " " << feat.m_gff_attributes["Name"] << endl;}
          
          //Shift the entire feature down the line
          feat.m_location.set_start_1(feat.m_location.get_start_1() + insert_length);
          feat.m_location.set_end_1(feat.m_location.get_end_1() + insert_length);
        }
        else //If we can't move the start, only move the end.  This is a modification of the feature
        {
          //Temporary variable for the new end position
          uint32_t end_temp = feat.m_location.get_end_1() + insert_length;
          
          //Notify the user of the action
          if(verbose){cout << "MODIFY\t" << feat["type"] << "\t" << feat.m_gff_attributes["ID"] << " " << feat.m_gff_attributes["Name"] << endl;}
          
          //Modify the end position
          feat.m_location.set_end_1(feat.m_location.get_end_1() + insert_length);
          
          //Mark it as pseudo
          feat.flag_pseudo(verbose);
          
          //Modify the notes for this feature
          //if(feat["type"] != "region" && feat["type"] != "source" && !feat.m_gff_attributes.count("Mutation from " + mut_type))feat.m_gff_attributes["Note"].push_back("Mutation from " + mut_type);
        }            
      }
    }
  }
  
  // Repeat Feature at Position
  // Look through ref_seq_info for repeat_name.
  // Using the strand, insert the repeat feature at supplied position.
  // Also repeat any features inside the repeat.
  // verbose outputs messages to console.
  void cAnnotatedSequence::repeat_feature_1(int32_t pos, int32_t start_del, int32_t end_del, cReferenceSequences& ref_seq_info, const string &repeat_name, int8_t strand, int32_t region_pos, bool verbose)
  {
    // Go through all the reference sequences in ref_seq_info
    for (vector<cAnnotatedSequence>::iterator itr_seq = ref_seq_info.begin(); itr_seq != ref_seq_info.end(); itr_seq++)
    {
      cAnnotatedSequence& this_seq = *itr_seq;
      cSequenceFeatureList& repeats = this_seq.m_repeats;
      
      // This sorting should be completely unnecessary. To be
      // safe, we sort here only because the following code
      // assumes that these lists are in order from least to greatest.
      this_seq.m_features.sort();
      this_seq.m_genes.sort();
      this_seq.m_repeats.sort();
      
      // Go through the repeats of the sequence
      for (cSequenceFeatureList::iterator itr_rep = repeats.begin(); itr_rep != repeats.end(); itr_rep++)
      {
        cSequenceFeature& rep = **itr_rep;
        
        string rep_fam_seq = ref_seq_info.repeat_family_sequence(repeat_name, strand, region_pos);
        
        // Is this the repeat we're looking for?
        if (rep.SafeGet("name") == repeat_name && region_pos == rep.m_location.get_start_1())
        {      
          cSequenceFeatureList& features = this_seq.m_features;
          cSequenceFeatureList feat_list_new;        
                 
          // Go through ALL the features of the (this_seq) sequence
          for (cSequenceFeatureList::iterator itr_feat = features.begin(); itr_feat != features.end(); itr_feat++)
          {
            cSequenceFeature& feat = **itr_feat;
            
            // Does this feature start and end inside of the repeat?
            if(feat.get_start_1() >= rep.get_start_1() && feat.get_end_1() <= rep.get_end_1() && feat["type"] != "region" && feat["type"] != "source")
            {
              // Create a brand new feature, that we can love and cuddle.
              // This is where we copy the feature and all the attributes.
              cSequenceFeaturePtr fp(new cSequenceFeature(feat));              
              cSequenceFeature& feat_new = *fp;              
              feat_new.m_gff_attributes = feat.m_gff_attributes;              
              
              // Depending on the strand of the mutation,
              // we juggle the starts and ends here.
              if(strand < 0)
              {                
                feat_new.m_location.set_start_1( pos + (rep.get_end_1() - feat.get_end_1()) - start_del);  //Set start to position plus the difference betwen the repeat's end, and this feature's end
                if(start_del && feat_new.get_start_1() < pos)  {  //If the feature start got shifted below the repeat start, set the start to the pos and flag it as pseudo
                  feat_new.m_location.set_start_1( pos);
                  feat_new.flag_pseudo();  }
                feat_new.m_location.set_end_1( feat_new.get_start_1() + (feat.get_end_1() - feat.get_start_1()));  //Set end to start plus feature length
                if(end_del && (pos + (rep.get_end_1() - rep.get_start_1())) < (feat_new.get_end_1() + end_del))  {  //Do we have an end_del? Does this feature end in an area that got deleted?
                  feat_new.m_location.set_end_1( (pos + (rep.get_end_1() - rep.get_start_1())));
                  feat_new.flag_pseudo();  }
                feat_new.m_location.set_strand( -feat.get_strand());
                feat_list_new.push_front(fp);
              }
              else
              {
                feat_new.m_location.set_start_1( pos + (feat.get_start_1() - rep.get_start_1()) - start_del);  //Set start to position plus the difference between the repeat feature's start, and this feature's start
                if(start_del && feat_new.get_start_1() < pos)  {  //If the feature start got shifted below the repeat start, set the start to the pos and flag it as pseudo
                  feat_new.m_location.set_start_1( pos);
                  feat_new.flag_pseudo();  }
                feat_new.m_location.set_end_1( feat_new.get_start_1() + (feat.get_end_1() - feat.get_start_1()));  //Set end to start plus feature length
                if(end_del && (pos + (rep.get_end_1() - rep.get_start_1())) < (feat_new.get_end_1() + end_del))  {  //Do we have an end_del? Does this feature end in an area that got deleted?
                  feat_new.m_location.set_end_1( (pos + (rep.get_end_1() - rep.get_start_1())));
                  feat_new.flag_pseudo();  }
                feat_list_new.push_back(fp);
              }
              
              // This is where we let the user know that "very imporant" things are going on
              if(verbose){cout << "REPEAT\t" << feat["type"] << "\t" << feat.m_gff_attributes["ID"] << " " << feat.m_gff_attributes["Name"] << endl;}
            }
          }
          
          // Add these to all the correct feature lists, and in the right order.
          for (cSequenceFeatureList::reverse_iterator itr_feat = feat_list_new.rbegin(); itr_feat != feat_list_new.rend(); itr_feat++)
          {
            cSequenceFeature& feat = **itr_feat;
            
            // Load certain information into the main hash, so breseq knows to use it.
            // It's unlikely this will be necessary, but you can never be too cautious
            if (feat.m_gff_attributes.count("Note"))
              feat["product"] = join(feat.m_gff_attributes["Note"], ",");
            
            if (feat.m_gff_attributes.count("Alias"))
              feat["accession"] = join(feat.m_gff_attributes["Alias"], ",");
            
            if (feat.m_gff_attributes.count("Name"))
              feat["name"] = join(feat.m_gff_attributes["Name"], ",");
            
            this->feature_push_back(*itr_feat);
          }
          
          // Sort, because while I did manage to implement a way to insert
          // all of these new features at the correct positions and in the
          // right order without these functions, it was a stupid function
          m_features.sort();
          m_genes.sort();
          m_repeats.sort();
          
          return;
        }
      }
    }
    
    ASSERT(false, "Unknown Repeat\n\tType:\t" + repeat_name + "\n\tPos:\t" + to_string(region_pos));
  }
  
  // Load a complete collection of files and verify that sufficient information was loaded
  void cReferenceSequences::LoadFiles(const vector<string>& file_names)
  {
    list<string> sorted_unique_file_names(file_names.begin(), file_names.end());

    sorted_unique_file_names.unique();//Removes non-unique file names

    for(list<string>::const_iterator it = sorted_unique_file_names.begin(); it != sorted_unique_file_names.end(); it++)
    {
      this->LoadFile(*it);
    }
    
    this->Verify();
    this->m_initialized = true;
  }
  
  void cReferenceSequences::LoadFile(const string& file_name)
  {
    ifstream in(file_name.c_str());
    ASSERT(in.good(), "Could not open reference file: " +file_name);
    
    str_uint old_load = m_seq_id_loaded;
    map<string,pair<uint32_t,uint32_t> > old_seq_info;
    for(vector<cAnnotatedSequence>::iterator its = this->begin(); its != this->end(); its++)
    {      
      old_seq_info[(*its).m_seq_id].first = (*its).m_features.size();
      old_seq_info[(*its).m_seq_id].second = get_sequence_length((*its).m_seq_id);
    }

    //! Step 1: Determine what file format it is from first line
    string first_line;
    getline(in, first_line);
    in.close();

    vector<string> split_first_line = split_on_whitespace(first_line);
    FileType file_type = UNKNOWN;
    // Fasta?
    if (first_line[0] == '>') {
      file_type = FASTA;
    }
    // GenBank?
    else if (first_line.find("LOCUS") != string::npos) {
      file_type = GENBANK;
    }
    // GFF?
    else if (first_line.find("##gff-version 3") != string::npos) {
      file_type = GFF3;
    }
    // BULL? Lines have three "words"
    else if (split_first_line.size() == 3) {
      file_type = BULL;
    }
    // The file was not identified and will most likely cause Breseq to fail.
    else if(file_type == UNKNOWN) {
      ERROR("Could not determine format of reference file: " + file_name);
    }
    //! Step 2: Load appropriate file
    switch (file_type) {
      case GENBANK:
      {
        ifstream gbk_in(file_name.c_str());
        ASSERT(in.good(), "Could not open GenBank file: " +file_name);
        string line;
        bool eof_found = false;
        while (!gbk_in.eof() && getline(gbk_in,line))
        {
          if(GetWord(line) == "//")  {
            eof_found = true;  }
        }
        
        ASSERT(eof_found, file_name + "\nTHIS GENBANK FILE APPEARS TO BE INCOMPLETE.\nMAKE SURE YOU DOWNLOADED THE ENTIRE FILE.\nFILE NEEDS TO END WITH '//'.");
        
        ReadGenBank(file_name);
      }break;
        
      case FASTA:
      {
        ReadFASTA(file_name);
      }break;
        
      case BULL:
      {
        ReadBull(file_name);
      }break;
        
      case GFF3:
      {
        ReadGFF(file_name);
      }break;
        
      default:
        WARN("Could not load the reference file: " +file_name);
    }
      
    // To upper the most recently loaded sequence.
    to_upper((this->back()).m_fasta_sequence.m_sequence);
    
    //Here we check to see we haven't loaded some of the same information again.
    for(str_uint::iterator i = old_load.begin(); i != old_load.end(); i++)
    {      
      //This SEQ_ID has been loaded before
      if(m_seq_id_loaded[i->first] != i->second)
      {
        //Have features and the sequence BOTH been loaded?
        ASSERT(!(old_seq_info[i->first].first && old_seq_info[i->first].second), file_name + "\nANOTHER FILE(S) HAS ALREADY LOADED ALL INFORMATION FOR: " + i->first);
        
        //Have we only loaded features?
        if(old_seq_info[i->first].first)
        {
          ASSERT(old_seq_info[i->first].first == (*this)[i->first].m_features.size(), file_name + "\nANOTHER FILE HAS ALREADY LOADED FEATURES FOR: " + i->first);
        }
        //Have we only loaded the sequence?
        else if(old_seq_info[i->first].second)
        {
          switch (file_type) {
              
              //We've already loaded the sequence and we just tried to load a GENBANK?
            case GENBANK:{
              ifstream gbk_in(file_name.c_str());
              ASSERT(in.good(), "Could not open GenBank file: " +file_name);
              string line;
              while (!gbk_in.eof() && getline(gbk_in,line))
              {
                if(line.find("ORIGIN") != string::npos)
                {
                  getline(gbk_in,line);
                  if(GetWord(line) != "//")  {
                    ERROR(file_name + "\nANOTHER FILE HAS ALREADY LOADED THE SEQUENCE FOR: " + i->first);  }
                }
              }
            }break;
              
              //We've already loaded the sequence and we just tried to load a FASTA?
            case FASTA:{
              ERROR(file_name + "\nANOTHER FILE HAS ALREADY LOADED THE SEQUENCE FOR: " + i->first);
            }break;
              
              //We've already loaded the sequence and we just tried to load a GFF3?
            case GFF3:{
              cFastaFile gff3_in(file_name.c_str(), ios_base::in);
              ASSERT(!gff3_in.fail(), "Could not open GFF file: " + file_name);
              string line;
              while (!gff3_in.eof() && getline(gff3_in,line))
              {
                if(line.find("##FASTA") != string::npos)
                {
                  while (!gff3_in.eof() && getline(gff3_in,line))
                  {
                    //The line we're on should be a sequence name line
                    if(line[0] == '>')
                    {
                      //Remove the '>' from the line.
                      line.replace(0,1,"");
                      
                      //Is this line an EXACT match for the current SEQ_ID we're worried about?
                      if(line == i->first)
                      {
                        ERROR(file_name + "\nANOTHER FILE HAS ALREADY LOADED THE SEQUENCE FOR: " + i->first);
                      }
                    }
                  }
                }                
              }
            }break;
              
            default:
              break;
          }
        }
      }
      
     //ASSERT(m_seq_id_loaded[i->first] == i->second, file_name + "\nANOTHER FILE HAS ALREADY LOADED REFERENCE INFORMATION FOR: " + i->first);
    }
  }
  
  void cReferenceSequences::Verify()
  {
    bool Error = false;
    stringstream ss;    
    for (vector<cAnnotatedSequence>::iterator itr= this->begin(); itr != this->end(); itr++) {
      cAnnotatedSequence& as = *itr;
      if (!as.get_sequence_length()) {
        ss << "No sequence loaded for: " << as.m_seq_id << endl;
        Error = true;
      }
      else if ((uint32_t)as.m_length != as.get_sequence_length())  {
        ss << "FEATURES AND SEQUENCE DON'T MATCH: " << as.m_seq_id << endl;
        ss << "LENGTH:\t" << (uint32_t)as.m_length << "\tVS:\t" << as.get_sequence_length() << " (Sequence Length)" << endl;
        Error = true;
      }
      else if(!as.get_sequence_length())  {
        ss << "LIKELY CAUSE:" << endl;
        ss << "No nucleotide sequence loaded for " << as.m_seq_id <<  endl;
        ss << "\nMake sure that the reference you are loading contains a sequence." << endl;
        ss << "If loading a Genbank file, make sure it has the associated sequence" << endl;
        ss << "or load the sequence as a seperate FASTA file." << endl;
      }
      else  {        
        ss << "LIKELY CAUSE:" << endl;
        ss << "Multiple files loaded referring to the same sequence. Make" << endl;
        ss << "sure that any features loaded describe the sequence you" << endl;
        ss << "actually loaded." << endl;
      }
    }
    if (Error) ERROR(ss.str());
    if (this->empty()) ERROR("Reference files were not loaded");
  }

  void cReferenceSequences::ReadFASTA(const string &file_name) {
    
    cFastaFile ff(file_name, ios_base::in);
    cFastaSequence on_seq;

    uint32_t on_seq_id = this->size();

    while ( ff.read_sequence(on_seq) ) {
      
      // @JEB sorting the input files by type seems like a better way of doing this,
      //    but we can't do it just by looking at file name endings...
      //
      // If we have a dummy sequence id from loading a BULL feature file before a FASTA,
      // then fill it in and update the index
      if (m_seq_id_to_index.count(BULL_DUMMY_SEQ_ID)) {
        uint32_t seq_index = this->m_seq_id_to_index[BULL_DUMMY_SEQ_ID];
        this->m_seq_id_to_index.erase(BULL_DUMMY_SEQ_ID);
        (*this)[seq_index].m_seq_id = on_seq.m_name;
        m_seq_id_to_index[(*this)[seq_index].m_seq_id] = seq_index;       
      }
      else {
        this->add_new_seq(on_seq.m_name);
      }
      
      // copy the info over (could define an assignment operator...)
      cAnnotatedSequence& this_seq = (*this)[on_seq.m_name];
      this_seq.m_fasta_sequence = on_seq;
      this_seq.m_seq_id = on_seq.m_name;
      this_seq.m_description = on_seq.m_description;
      this_seq.m_length = on_seq.m_sequence.size();
    }
  }
    
  void cReferenceSequences::ReadFASTA(cFastaFile& ff) 
  {
    cFastaSequence on_seq;

    while ( ff.read_sequence(on_seq) ) {
      this->add_new_seq(on_seq.m_name);
      cAnnotatedSequence& new_seq = (*this)[on_seq.m_name];
      new_seq.m_fasta_sequence = on_seq;
      new_seq.m_seq_id = on_seq.m_name;
      new_seq.m_length = on_seq.m_sequence.size();
    }

  }

  void cReferenceSequences::WriteFASTA(const std::string &file_name, bool verbose) {
    
    if(verbose)cout << "Writing FASTA" << endl << "\t" << file_name << endl;
    
    cFastaFile ff(file_name, ios_base::out);
    for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
      if(it_as->m_length)ff.write_sequence(it_as->m_fasta_sequence);
    }
    
    if(verbose)cout << "\t**FASTA Complete**" << endl;
  }

  void cReferenceSequences::WriteFASTA(cFastaFile& ff, bool verbose) {

    if(verbose)cout << "Writing FASTA" << endl;
    
    for(vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
      if(it_as->m_length)ff.write_sequence(it_as->m_fasta_sequence);
    }
    
    if(verbose)cout << "\t**FASTA Complete**" << endl;
    
  }

  /*! ReadGFF abides by the following format:
   *  http://www.sequenceontology.org/gff3.shtml 
   */
  void cReferenceSequences::ReadGFF(const string& file_name)
  {
    cFastaFile in(file_name.c_str(), ios_base::in);
    ASSERT(!in.fail(), "Could not open GFF file: " + file_name);
    //! Step 1: Header //! cFastaFile removes first line, so no header...
    string line;
    map<string, vector<cSequenceFeature*> > id_features_table; 
    while (!in.eof() && getline(in,line)) {
      cSequenceFeaturePtr fp(new cSequenceFeature);
      cSequenceFeature& feature = *fp;

  //! Step 2: Check for GFF Tags, reads FASTA from here.
      // We are concerned about a couple of ## GFF tags
      if (line[0] == '#') {
        if (line.find("##species") != string::npos) {
          //TODO @GRC
          continue;
        }
        else if (line.find("##sequence-region") != string::npos) {
          // Line of form ##sequence-region seqid start end
          stringstream ls(line);
          string x, seq_id, start, end;
          ls >> x >> seq_id >> start >> end;
          this->add_new_seq(seq_id);
          (*this)[seq_id].m_length = from_string<uint32_t>(end);
          continue;
        }
        // Find embedded fasta file
        else if (line.find("##FASTA") != string::npos) {
          /* Things admittedly get a bit hairy right here, you have to
            take the next line (The one with ">XXXXXX") and set it as the
            current line for cFastaFile, after you've read the fasta you
            can exit the function since nothing should exist after.
            */
          getline(in,line);
          in.set_current_line(line);
          this->ReadFASTA(in);
          continue;
        } else {
          continue;
        }
      }
  /*! Step 3: Split line on tabs("\t") until last column, grab last column until endl("\n"),
      the default for getline(). !*/
      
      stringstream ss(line);
      string seq_id, start, end, strand;

      // Handle columns up to the last one, "attributes"
      // Column 1: "seqid"
      getline(ss, seq_id, '\t');
      // Column 2: "source"
      getline(ss, feature["source"], '\t');
      // Column 3: "type"
      getline(ss, feature["type"], '\t');
      // Column 4: "start"
      getline(ss, start, '\t');
      feature.m_location.set_start_1( from_string<uint32_t>(start));
      // Column 5: "end"
      getline(ss, end, '\t');
      feature.m_location.set_end_1( from_string<uint32_t>(end));
      // Column 6: "score"
      getline(ss, feature["score"], '\t');
      // Column 7: "strand"
      getline(ss, strand, '\t');
      feature.m_location.set_strand(0);
      if (strand == "+")
        feature.m_location.set_strand(1);
      else if (strand == "-")
        feature.m_location.set_strand(-1);        
      // Column 8: "phase"
      getline(ss, feature["phase"], '\t');
      // Column 9: "attributes"
      string raw_attributes;

      // Handle parsing the attributes
      getline(ss, raw_attributes);
      vector<string> attributes = split(raw_attributes, ";");

      //Split attribute's key and value by "="
      for (vector<string>::iterator itr = attributes.begin(); itr != attributes.end(); itr++) {
        string& attribute = *itr;
        vector<string> key_value = split(attribute,"=");
        string& key = key_value.front();
        string& value = key_value.back();
      //! Case 2: Multiple values for given key, split by ","
        feature.m_gff_attributes[key] = split(value, ",");
        // unescape special characters after splitting
        for (uint32_t i=0; i<feature.m_gff_attributes[key].size(); i++)  {
          feature.m_gff_attributes[key][i] = GFF3UnescapeString(feature.m_gff_attributes[key][i]);
        }
      }
            
      // Load certain information into the main hash, so breseq knows to use it
      if (feature.m_gff_attributes.count("Note"))
        feature["product"] = join(feature.m_gff_attributes["Note"], ",");
      
      if (feature.m_gff_attributes.count("Alias"))
        feature["accession"] = join(feature.m_gff_attributes["Alias"], ",");

      if (feature.m_gff_attributes.count("Name"))
        feature["name"] = join(feature.m_gff_attributes["Name"], ",");
      
      if (feature.m_gff_attributes.count("Pseudo"))
        feature.m_pseudo = from_string<bool>(feature.m_gff_attributes["Pseudo"][0]);
      
      // Load translation table over to the gene
      if (feature.m_gff_attributes.count("transl_table"))
        feature["transl_table"] = feature.m_gff_attributes["transl_table"][0];

    
      //! Step 4: Determine if sequence already exists (find or create if not found)
      this->add_new_seq(seq_id);
      
      // Handle features that cross the origin by adding them twice (only one will be written)
      // @GRC Should this be handled with cLocation rather then a 'bonus feature'
      if (feature.get_end_1() > (*this)[seq_id].m_length) {
        cSequenceFeaturePtr bonus_circular_feature_ptr(new cSequenceFeature);
        *bonus_circular_feature_ptr = feature;
        bonus_circular_feature_ptr->m_location.set_start_1(bonus_circular_feature_ptr->get_start_1() + 1 - (*this)[seq_id].m_length);
        bonus_circular_feature_ptr->m_location.set_end_1(bonus_circular_feature_ptr->get_end_1() - (*this)[seq_id].m_length);
        (*this)[seq_id].feature_push_front(bonus_circular_feature_ptr);
      }
    
      // If this is a landmark "region" corresponding to the entire chromosome grab extra information
      if ((feature["type"] == "region") && (feature.get_start_1() == 1) && (feature.get_end_1() == (*this)[seq_id].m_length)) {
        if (feature.m_gff_attributes.count("Is_circular"))
          (*this)[seq_id].m_is_circular = (feature.m_gff_attributes["Is_circular"][0] == "true");

        if (feature.m_gff_attributes.count("Note"))
          (*this)[seq_id].m_description = feature.m_gff_attributes["Note"][0];
      }

      // Add feature / determine if it needs to be added as a sub_location to another feature.
      if (feature.m_gff_attributes.count("ID")) {
        string id = feature.m_gff_attributes["ID"][0];
        if (id_features_table.count(id))  {
          vector<cSequenceFeature*> features = id_features_table[id];
          bool is_sub_location = false;
          
          //Search to see if 'type' matches also.
          for (uint32_t i = 0; i < features.size();  ++i) {
            if ((*features[i])["type"] == feature["type"]) {
              //Update cLocation to include this sublocation.
              features[i]->m_location.add_sub_location(feature.m_location);
              is_sub_location = true;
              break;
            }
          }

          if (is_sub_location) {
            continue; //Don't add to feature lists.
          } else { 
            id_features_table[id].push_back(&feature);
          }

        }
      }

      (*this)[seq_id].feature_push_back(fp);

    }
    
    // sort because they may now be out of order
    for (vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
      // sort the list
      it_as->m_features.sort();
      it_as->m_genes.sort();
      it_as->m_repeats.sort();
    }
  }
  
  
/*! WriteGFF abides by the following format:
  http://www.sequenceontology.org/gff3.shtml !*/
void cReferenceSequences::WriteGFF( const string &file_name, bool verbose ){

  cFastaFile out(file_name.c_str(), ios_base::out);

  ASSERT(!out.fail(), "Failed to open file " + file_name);
    
  //Notify the user of output
  if(verbose)
  {
    cout << "Writing GFF3" << endl;
    cout << "\t" << file_name << endl;
  }
    
  //! Step 1: Header
  out << "##gff-version 3" << endl;

  for (vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
    if(it_as->m_length)out << "##sequence-region" << "\t" << it_as->m_seq_id << "\t" << "1" << "\t" << it_as->m_length << endl;
  }
  
  //! Step 2: Features
  for (vector<cAnnotatedSequence>::iterator it_as = this->begin(); it_as < this->end(); it_as++) {
    
    // sort the list
    it_as->m_features.sort();
    
    for (cSequenceFeatureList::iterator it = it_as->m_features.begin(); it != it_as->m_features.end(); it++) {
      cSequenceFeature& feat = **it;

      // skip first example of doubly-loaded feature that overlaps circular genome
      if (feat.get_start_1() < 1)
        continue;

      enum {SEQID = 0, SOURCE, TYPE, START_1, END_1, SCORE, STRAND, PHASE, ATTRIBUTES};
      vector<string> columns(9, ".");

      columns[SEQID] = it_as->m_seq_id;
      if (feat.count("source")) columns[SOURCE] = feat["source"];
      if (feat.count("type"))   columns[TYPE] = feat["type"];
      //Note: Will add start_1 and end_1 below.
      if (feat.count("score"))  columns[SCORE] = feat["score"];
      columns[STRAND] = feat.is_top_strand() ? "+" : "-";
      if (feat.count("phase"))  columns[PHASE] = feat["phase"];

      //Attributes
      vector<string> attributes;
      map<string,vector<string> >::const_iterator itr;
      for (itr = feat.m_gff_attributes.begin(); itr != feat.m_gff_attributes.end(); itr++) {
        const string& key = itr->first;
        const vector<string>& values = itr->second;
        
        vector<string> s;
        for (vector<string>::const_iterator it = values.begin(); it != values.end(); it++) {
          s.push_back(GFF3EscapeString(*it));
        }
        
        if(s.size())attributes.push_back(key + "=" + join(s, ","));
      }
      columns[ATTRIBUTES] = join(attributes, ";");

      const vector<cLocation> locations = feat.m_location.get_all_sub_locations();
      for (uint32_t i = 0; i < locations.size(); ++i) {
        columns[START_1] = to_string(locations[i].get_start_1());
        columns[END_1]   = to_string(locations[i].get_end_1());

        out << join(columns, "\t") << endl;

      }


    }
  }
  //! Step 3: Fasta, one command writes out all sequences
  out << "##FASTA" << endl;
  this->WriteFASTA(out);
  
  out.close();
    
  //Notify the user of output
  if(verbose)
  {
    cout << "\t**GFF3 Complete**" << endl;
  }
}

void cReferenceSequences::ReadGenBank(const string& in_file_name) {

  ifstream in(in_file_name.c_str(), ios_base::in);
  ASSERT(!in.fail(), "Could not open GenBank file: " + in_file_name);

  while (ReadGenBankFileHeader(in)) {
    
    cAnnotatedSequence& s = this->back();
    
    // add a 'region' feature for GFF3 output
    cSequenceFeaturePtr f(new cSequenceFeature);
    (*f)["type"] = "region";
    f->m_location.set_strand( 1);
    f->m_location.set_start_1( 1);
    f->m_location.set_end_1( s.m_length);
    
    if (s.m_is_circular)
      f->m_gff_attributes["Is_circular"].push_back("true");
    else
      f->m_gff_attributes["Is_circular"].push_back("false");
    
    f->m_gff_attributes["Note"].push_back(s.m_description);
    s.feature_push_back(f);
    
    ReadGenBankFileSequenceFeatures(in, s);
    ReadGenBankFileSequence(in, s);    
  }
}

bool cReferenceSequences::ReadGenBankFileHeader(ifstream& in) {

  //std::cout << "header" << std::endl;
  string line;
  bool found_LOCUS_line = false;
  cAnnotatedSequence* s = NULL;
  while (!in.eof()) {
    std::getline(in, line);
    string first_word = GetWord(line);
    RemoveLeadingTrailingWhitespace(line);

    // This is the first line
    if (first_word == "LOCUS") {

      // Example line 
      // LOCUS       EB03                 4629813 bp    DNA     circular BCT 13-APR-2007
      
      string w;
      w = GetWord(line);
      string seq_id = w;
      
      this->add_new_seq(seq_id);
      s = &((*this)[seq_id]);
      
      w = GetWord(line);
      s->m_length = atoi(w.c_str());
      
      w = GetWord(line);
      w = GetWord(line);
      w = GetWord(line);
      if (to_lower(w) == "circular")
        s->m_is_circular = true;

      // Should only be one line like this per record!
      ASSERT(!found_LOCUS_line, "Multiple LOCUS lines found in single GenBank record."); 
      found_LOCUS_line = true;
    }

    if (first_word == "DEFINITION") {
      ASSERT(s, "Missing LOCUS line before DEFINITION line in GenBank record.");
      s->m_description = line;
    }

    if (first_word == "FEATURES") break;
  }

  return (found_LOCUS_line);
}

/*! Load start, end, strand for feature from a GenBank coordinate string.
 *
 *  The string may cover multiple lines. Currently we do not handle discontinuous features,
 *  and features that cross the origin of a circular chromosome are returned with end < start. 
 */
void cSequenceFeature::ReadGenBankCoords(string& s, ifstream& in) {

// Typical coordinate strings:
//   1485..1928
//   complement(2644..3159)
//  
// Example from a circular genome:
//   complement(join(7504..7835,1..163))
//
  
  // Read through all parentheses and across multiple lines if necessary
  int32_t parentheses_level = 0;
  size_t parenthesis_pos = s.find_first_of("()");
  
  while(parenthesis_pos != string::npos) {
  if (s.at(parenthesis_pos) == '(') {
      parentheses_level++;
    } else {
      parentheses_level--;
    }
    parenthesis_pos = s.find_first_of("()", parenthesis_pos+1);
  }
    
  // Add in remaining lines
  string s2;
  while ((parentheses_level > 0) && !in.eof()) {
    std::getline(in, s2);
    RemoveLeadingTrailingWhitespace(s2);
  
    size_t on_pos = s2.find_first_of("()");
    while(on_pos != string::npos) {
  
      if (s2.at(on_pos) == '(') {
        parentheses_level++;
      } else {
        parentheses_level--;
      }
      on_pos = s2.find_first_of("()", on_pos+1);
    }
    s += s2;
  }
  ASSERT(parentheses_level == 0, "Unmatched parenthesis in feature location.");
  
  // Clean <> and space characters.
  s = substitute(s, "<", "");
  s = substitute(s, ">", "");
  s = substitute(s, " ", "");
  s = substitute(s, "\t", "");
  
  m_location = ParseGenBankCoords(s);

  return;
}
  
// Parses a sub-part of the full location string
cLocation cSequenceFeature::ParseGenBankCoords(string& s, int8_t in_strand)
{
  cLocation loc;
  
  //Look for an operation.
  //complement().
  if (s.find("complement(") == 0) {
    uint32_t n = string("complement(").size();
    string value = s.substr(n, s.size() - n - 1);
    loc = ParseGenBankCoords(value, -1 * in_strand);
  }
  //join()
  else if (s.find("join(") == 0) {
    uint32_t n = string("join(").size();
    vector<string> tokens = split(s.substr(n, s.size() - n - 1), ",");
    
    for (uint32_t i = 0; i < tokens.size(); ++i) {
      cLocation sub_loc = ParseGenBankCoords(tokens[i], in_strand);
      loc.m_sub_locations.push_back( sub_loc);
    }
    
    //Get cLocation::start,end from outer sub_locations.
    loc.set_start_1(loc.m_sub_locations.front().get_start_1());
    loc.set_end_1(loc.m_sub_locations.back().get_end_1());
  }
  else {
    //Split on .. (usually 2)
    vector<string> tokens = split(s, "..");
    loc.set_strand(in_strand);
    loc.set_start_1(atoi(tokens.front().c_str()));
    loc.set_end_1(atoi(tokens.back().c_str()));
  }
  
  return loc;
}

void cSequenceFeature::ReadGenBankTag(std::string& tag, std::string& s, std::ifstream& in) {

  // delete leading slash
  tag.erase(0,1);

  // erase through the equals on this line
  int pos = s.find("=");
  s.erase(0,pos+1);

  // If there is a quote then we need to continue to the last quote
  size_t first_quote_pos = s.find("\"");

  if (first_quote_pos == string::npos) {
    (*this)[tag] = s;
    return;
  }

  s.erase(0,first_quote_pos+1);

  size_t second_quote_pos = s.find("\"");

  // One liner
  if (second_quote_pos != string::npos) {
    s.erase(second_quote_pos,s.length());
    (*this)[tag] = s;
    return;
  }

  // If the value is quoted, we have to read additional lines until end quote
  string value = s;

  bool found_last_quote = false;
  while (!found_last_quote && !in.eof()) {
    std::getline(in, s);
    RemoveLeadingTrailingWhitespace(s);

    second_quote_pos = s.find("\"");
    if (second_quote_pos != string::npos) {
      s.erase(second_quote_pos,s.length());
      found_last_quote = true;
    }

    if (tag != "translation") value += " ";
    value += s;
  }
  assert(found_last_quote);

  (*this)[tag] = value;
}

void cReferenceSequences::ReadGenBankFileSequenceFeatures(std::ifstream& in, cAnnotatedSequence& s) {
  //std::cout << "features" << std::endl;
  cSequenceFeature* current_feature(NULL);
  cSequenceFeatureList all_features; // make preliminary list then add once entries are complete
  string line;
  while (!in.eof()) {
    getline(in, line);

    //cout << line << endl;
    string first_word = GetWord(line);

    //std::cout << first_word << "::" << line << std::endl;

    // Done with this section...
    if (first_word == "ORIGIN") break;

    // Major tag = new feature or information block
    if (first_word[0] != '/') {

      if (first_word != "BASE") {
        cSequenceFeaturePtr new_feature(new cSequenceFeature);
        all_features.push_back( new_feature );
        current_feature = new_feature.get();
        (*current_feature)["type"] = first_word;
        // parse the rest of the line
        std::string coord_s = GetWord(line);
        current_feature->ReadGenBankCoords(coord_s, in);
      }
    }
    // Minor tag = information about current feature
    else {
      ASSERT(current_feature, "No current feature.");

      // Remove leading slash
      current_feature->ReadGenBankTag(first_word, line, in); // reads multi-line entries
    }
  }


  for (cSequenceFeatureList::iterator it = all_features.begin(); it != all_features.end(); it++) {
    cSequenceFeature& feature = **it;

    // common changes for any type
    // use /note as the product if there is not product
    if (feature.SafeGet("product") == "")
    {
      feature["product"] = feature.SafeGet("note");
    }

    if (feature["type"] == "repeat_region" || feature["type"] == "mobile_element") {
      // Don't add unnamed ones to the list...
      //if (it->SafeGet("mobile_element") == "") continue;

      feature["name"] = "repeat_region";
      
      // Give the repeat region a name if NOTHING else can be found
      if (feature.SafeGet("note") != "")
        feature["name"] = feature.SafeGet("note");

      // E. coli case:
      if (feature.count("mobile_element") || 
          feature.count("mobile_element_type"))
      {
        if (feature.count("mobile_element")) {  
          feature["name"] = feature["mobile_element"];
        }
        if (feature.count("mobile_element_type")) {
          feature["name"] = feature["mobile_element_type"];
        }

        string& name = feature["name"];

        // remove prefix
        int pos = name.find("insertion sequence:");
        if (pos != -1)
          name.erase(pos,pos+19);

        // remove suffix if "IS(\d)}
        pos = name.find_first_of("0123456789");
        if (pos != -1) {
          int found = name.find_first_not_of("0123456789", pos);
          if (found != -1) {
            name.erase(found, name.length());
          }
        }
      }      

      // S. cerevisiae case
      if (feature.SafeGet("rpt_family") != "")
        feature["name"] = feature["rpt_family"];

      //std::cerr << (*it).SafeGet("mobile_element") << " " << (*it).SafeGet("name") << std::endl;

      if (feature.SafeGet("product") == "")
        feature["product"] = "repeat region";
    }
    else
    {
      // Add information
      if (feature.SafeGet("gene") != "")
        feature["name"] = feature["gene"];

      if (feature.SafeGet("transl_table") != "")
        feature["transl_table"] = feature["transl_table"];
      
      if ( (feature.SafeGet("name") == "") && (feature.SafeGet("locus_tag") != "") )
        feature["name"] = feature["locus_tag"];

      //std::cerr << (*it).SafeGet("name") << " " << (*it).SafeGet("gene") << " " << (*it).SafeGet("locus_tag") << std::endl;

      feature["accession"] = feature.SafeGet("locus_tag");

      // Is this feature marked as pseudo?
      feature.m_pseudo = (feature.count("pseudo") != 0);
    }
        
    // transfer to GFF
    feature["phase"] = "0";
    if (feature.SafeGet("locus_tag") != "")
      feature.m_gff_attributes["ID"] = make_vector<string>(feature["locus_tag"]);
    if (feature.SafeGet("product") != "")
      feature.m_gff_attributes["Note"] = make_vector<string>(feature["product"]);
    if (feature.m_pseudo) 
      feature.m_gff_attributes["Pseudo"].push_back("true");
    if (feature.SafeGet("accession") != "")
      feature.m_gff_attributes["Alias"] = make_vector<string>(feature["accession"]);
    if (feature.SafeGet("name") != "")
      feature.m_gff_attributes["Name"] = make_vector<string>(feature["name"]);

    if (feature.SafeGet("transl_table") != "")
      feature.m_gff_attributes["transl_table"] = make_vector<string>(feature["transl_table"]);

    
    // add an extra copy of the feature if it crosses the origin of a circular chromosome
    if (feature.get_end_1() < feature.get_start_1()) {
      cSequenceFeaturePtr bonus_circular_feature(new cSequenceFeature);
      *bonus_circular_feature = feature;
      bonus_circular_feature->m_location = cLocation(
                                                     bonus_circular_feature->get_start_1() - s.m_length + 1, 
                                                     feature.m_location.get_end_1() + s.m_length,
                                                     bonus_circular_feature->m_location.get_strand()
                                                     );
      s.feature_push_front( bonus_circular_feature );
    }
    
    s.feature_push_back(*it);
  }
  
  // sort the list
  s.m_features.sort();
  s.m_genes.sort();
  s.m_repeats.sort();
}

void cReferenceSequences::ReadGenBankFileSequence(std::ifstream& in, cAnnotatedSequence& s) {

  s.m_fasta_sequence.m_name = s.m_seq_id;

  string line;
  while (!in.eof()) {
    getline(in, line);
    string first_word = GetWord(line);
    RemoveLeadingTrailingWhitespace(line);

    //cerr << first_word << "::" << line << std::endl;

    if (first_word == "//") break;
    for(string::iterator i=line.begin(); i!=line.end(); ++i) {
      if (*i == ' ') continue;
      *i = toupper(*i);

      // Scrub nonstandard characters
      if ((*i != 'A' ) && (*i != 'T' ) && (*i != 'C' ) && (*i != 'G' )) {
        *i = 'N';
      }

      s.m_fasta_sequence.m_sequence += *i;
    }
  }

  //cout << s.m_sequence << std::endl;
}

void cReferenceSequences::ReadBull(const string& file_name) {

  // This should only work if we have zero or one reference sequence, because it does not
  // give a seq_id, we assume the FASTA before or the FASTA after is paired with it.
  ASSERT(this->size() <= 1, "Bull format currently only works with one reference sequence");

  // if we encounter the BULL file first, create an empty sequence
  if (this->size() == 0)
  {
    // if one already exists, throw an error
    ASSERT(!this->m_seq_id_to_index.count(BULL_DUMMY_SEQ_ID), 
           "Two BULL formatted feature files loaded in a row. Ambiguous assignment to FASTA files.");
    this->add_new_seq(BULL_DUMMY_SEQ_ID);
  }
    
  ifstream in(file_name.c_str());
  if(!in.good()) WARN("Could not open file:" + file_name);

  char line[10];

  cSequenceFeatureList all_features;

  while ( !in.eof() ) {

    string line;
    getline(in, line);
    vector<string> s = split_on_whitespace(line);

    if (s.size() == 0) continue;
    ASSERT(s.size() == 3, "Expected 3 words in line, found " + to_string(s.size()) + " in line:\n" + line);
    
    cSequenceFeaturePtr current_feature(new cSequenceFeature);
    
    current_feature->m_location.set_start_1( from_string<uint32_t>(s[0]));
    current_feature->m_location.set_end_1( from_string<uint32_t>(s[1]));
    current_feature->m_location.set_strand( (current_feature->get_start_1() <= current_feature->get_end_1()) ? 1 : -1);
    
    // transfer to GenBank
    (*current_feature)["name"] = s[2];
    (*current_feature)["type"] = "CDS";
    (*current_feature)["product"] = "";
    (*current_feature)["accession"] = s[2];
    
    // transfer to GFF3
    (*current_feature)["ID"] = s[2];
    (*current_feature)["Note"] = "";
    (*current_feature)["Alias"] = s[2];
    (*current_feature)["Name"] = s[2];
    
    all_features.push_back(current_feature);
  }
  in.close();
  
  cAnnotatedSequence& s = this->back();
  for (cSequenceFeatureList::iterator it = all_features.begin(); it != all_features.end(); it++) {
    cSequenceFeature& feat = **it;
    s.m_features.push_back(*it);
  }
}

/*! Returns the nucleotide sequence of a typical copy of a repeat.
 *
 */  
string cReferenceSequences::repeat_family_sequence(
                                                   const string &repeat_name, 
                                                   int8_t strand, 
                                                   int32_t &region_pos, 
                                                   string* seq_id,
                                                   cSequenceFeature* repeat_feature_picked
                                                   )
{  
  counted_ptr<cSequenceFeature> picked_rep(NULL);
  cAnnotatedSequence* picked_seq(NULL);
  
  map<string, uint32_t> repeat_sequence_count;
  map<string, uint32_t> repeat_sequence_pos;  
  map<uint32_t, uint32_t> repeat_size_count;  
  map<uint32_t, uint32_t> repeat_size_pos;
  map<uint32_t, string> loaded_elements_positions;
  set<string> loaded_elements;
  uint32_t largest_name = 0;
  
  // Loop through all reference sequences
  for (vector<cAnnotatedSequence>::iterator itr_seq = this->begin(); itr_seq != this->end(); itr_seq++) {
    cAnnotatedSequence& this_seq = *itr_seq;
    cSequenceFeatureList& repeats = this_seq.m_repeats;
    
    for (cSequenceFeatureList::iterator itr_rep = repeats.begin(); itr_rep != repeats.end(); itr_rep++) {
      cSequenceFeature& rep = **itr_rep;
      
      // This section is for error information later
      loaded_elements.insert(rep.SafeGet("name"));
      loaded_elements_positions[rep.get_start_1()] = rep.SafeGet("name");
      if(largest_name < rep.SafeGet("name").size())largest_name = rep.SafeGet("name").size();
      // End error info
      
      if(rep.SafeGet("name") != repeat_name)
        continue;
      
      // Stores all the sequences that match so we can compare them.
      if (region_pos < 0)
      {
        string adjSeq = this_seq.get_sequence_1(rep.get_start_1(), rep.get_end_1());
        if (strand != rep.get_strand())
          adjSeq = reverse_complement(adjSeq);
        
        repeat_size_count[rep.get_end_1() - rep.get_start_1() + 1]++;        
        repeat_size_pos[rep.get_end_1() - rep.get_start_1() + 1] = rep.get_start_1();
        repeat_sequence_count[adjSeq]++;
        repeat_sequence_pos[adjSeq] = rep.get_start_1();
      }
    }
  }
  
  // This will set the region_pos to whichever sequence has the most number of copies.
  // If number of copies is not greater than 1, we will use the most common size instead.
  if((region_pos < 0) && (repeat_sequence_count.size()))
  {    
    if (max_element(repeat_sequence_count.begin(), repeat_sequence_count.end(), map_comp_second<string, uint32_t>)->second > 1) {
      region_pos = repeat_sequence_pos[max_element(repeat_sequence_count.begin(), repeat_sequence_count.end(), map_comp_second<string, uint32_t>)->first];  }
    else  {
      region_pos = repeat_size_pos[max_element(repeat_size_count.begin(), repeat_size_count.end(), map_comp_second<uint32_t, uint32_t>)->first];  }
  }
  
  // Loop through all reference sequences
  for (vector<cAnnotatedSequence>::iterator itr_seq = this->begin(); itr_seq != this->end(); itr_seq++) {
    cAnnotatedSequence& this_seq = *itr_seq;
    cSequenceFeatureList& repeats = this_seq.m_repeats;
    
    for (cSequenceFeatureList::iterator itr_rep = repeats.begin(); itr_rep != repeats.end(); itr_rep++) {
      cSequenceFeature& rep = **itr_rep;
      
      if(rep.SafeGet("name") != repeat_name)
        continue;
      
      if (region_pos == rep.get_start_1()) {
        picked_seq = &this_seq;
        picked_rep = *itr_rep;
        break;
      }
    }
  }
  
  stringstream ssRepList;
  ssRepList << "LOADED REPEAT ELEMENTS:";  
  stringstream ssRepPos;
  ssRepPos << "LOADED REPEAT ELEMENTS POSITIONS:";
  
  for(set<string>::iterator k = loaded_elements.begin(); k != loaded_elements.end(); k++)  {
    ssRepList << "\n\t" << *k;  }
  
  for(map<uint32_t, string>::iterator jk = loaded_elements_positions.begin(); jk != loaded_elements_positions.end(); jk++) 
  {
    ssRepPos << "\n\t" << (*jk).second;
    for(uint32_t ik = (*jk).second.size(); ik < largest_name; ik++)  {
      ssRepPos << " ";  }    
    ssRepPos << "\t" << to_string((*jk).first); 
  }
    
  if (region_pos != -1)
    ASSERT(picked_rep.get() != NULL, "No valid " + repeat_name + " repeat element found for position " + to_string(region_pos) + "\n" + ssRepPos.str());

  ASSERT(picked_rep.get() != NULL, "No valid " + repeat_name + " found.\n" + ssRepList.str());
    
  string repeat_seq = picked_seq->get_sequence_1(picked_rep->get_start_1(), picked_rep->get_end_1());
  if (strand != picked_rep->get_strand())
    repeat_seq = reverse_complement(repeat_seq);
  
  if (seq_id)
    *seq_id = picked_seq->m_seq_id;
  if (repeat_feature_picked)
    *repeat_feature_picked = *picked_rep;
  return repeat_seq;  
}

/*! Find the closest edge of a repeat in the specified direction within the specified distance
 */
cSequenceFeaturePtr cReferenceSequences::find_closest_repeat_region_boundary(int32_t position, cSequenceFeatureList& repeat_list, int32_t& max_distance, int32_t direction)
{
  cSequenceFeaturePtr repeat_ptr(NULL);

  for (cSequenceFeatureList::iterator it = repeat_list.begin(); it != repeat_list.end(); ++it) {
    cSequenceFeaturePtr test_repeat_ptr = *it;
    
    // Distance from the appropriate end of the repeat
    int32_t test_distance = abs(static_cast<int32_t>(((direction == -1) ? position - test_repeat_ptr->get_end_1() : test_repeat_ptr->get_start_1() - position)));
    
    // We want the closest one without going over that is within max_distance
    if ( (test_distance >= 0) && (test_distance <= max_distance) ) {
      repeat_ptr = test_repeat_ptr;
      max_distance = test_distance;
    }
  }
  return repeat_ptr;
}

/*! Returns the last feature encountered that overlaps a position
 */
cSequenceFeaturePtr cReferenceSequences::get_overlapping_feature(cSequenceFeatureList& feature_list, int32_t pos)
{
  cSequenceFeaturePtr feature_ptr(NULL);
  for (cSequenceFeatureList::iterator it = feature_list.begin(); it != feature_list.end(); ++it) {
    cSequenceFeaturePtr test_feature_ptr = *it;
    if (pos >= test_feature_ptr->get_start_1() && pos <= test_feature_ptr->get_end_1())
      feature_ptr = test_feature_ptr;
  }
  return feature_ptr;
}

void cReferenceSequences::find_nearby_genes(
                                            cSequenceFeatureList& gene_list,
                                            int32_t pos_1,
                                            int32_t pos_2,
                                            vector<Gene>& within_genes,
                                            vector<Gene>& between_genes,
                                            vector<Gene>& inside_left_genes,
                                            vector<Gene>& inside_right_genes,
                                            Gene& prev_gene,
                                            Gene& next_gene
                                            )
{
  for (cSequenceFeatureList::iterator it = gene_list.begin(); it != gene_list.end(); ++it)
  {
    cSequenceFeature& test_gene_feat = **it;
    Gene test_gene = Gene(test_gene_feat); // up-cast, should be a better way to do this @JEB

    if (test_gene.get_end_1() < pos_1)
      prev_gene = test_gene;

    if (  (test_gene.get_start_1() <= pos_1) && (test_gene.get_end_1() >= pos_1)
       && (test_gene.get_start_1() <= pos_2) && (test_gene.get_end_1() >= pos_2) )
    {
      within_genes.push_back(test_gene);
    }
    else if ( (test_gene.get_start_1() <= pos_1) && (test_gene.get_end_1() >= pos_1) )
    {
      inside_left_genes.push_back(test_gene);
    }
    else if ( (test_gene.get_start_1() <= pos_2) && (test_gene.get_end_1() >= pos_2) )
    {
      inside_right_genes.push_back(test_gene);
    }
    else if ( (test_gene.get_start_1() >= pos_1) && (test_gene.get_end_1() <= pos_2) )
    {
      between_genes.push_back(test_gene);
    }
    // We've passed the changes, so it is in the previous intergenic space
    if (test_gene.get_start_1() > pos_2)
    {
      next_gene = test_gene;
      break;
    }
  }
}

// Retrieved from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi on 04-21-2012
vector<string> cReferenceSequences::translation_tables = make_vector<string>
  ("") // 0 placeholder to make vector 1-indexed
  ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 1
  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG") // 2
  ("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 3
  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 4
  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG") // 5
  ("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 6
  ("") // 7 deleted
  ("") // 8 deleted
  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG") // 9
  ("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 10
  ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 11
  ("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 12
  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG") // 13
  ("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG") // 14
  ("FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 15
  ("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 16
  ("") // 17 deleted
  ("") // 18 deleted
  ("") // 19 deleted
  ("") // 20 deleted
  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG") // 21
  ("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 22
  ("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 23
  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG") // 24
;
  
  vector<string> cReferenceSequences::initiation_codon_translation_tables = make_vector<string>
  ("") // 0 placeholder to make vector 1-indexed
  ("FFLMSSSSYY**CC*WLLLMPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 1
  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRMMMMTTTTNNKKSS**VVVMAAAADDEEGGGG") // 2
  ("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 3
  ("FFMMSSSSYY**CCWWLLLMPPPPHHQQRRRRMMMMTTTTNNKKSSRRVVVMAAAADDEEGGGG") // 4
  ("FFLMSSSSYY**CCWWLLLLPPPPHHQQRRRRMMMMTTTTNNKKSSSSVVVMAAAADDEEGGGG") // 5
  ("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 6
  ("") // 7 deleted
  ("") // 8 deleted
  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVMAAAADDEEGGGG") // 9
  ("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 10
  ("FFLMSSSSYY**CC*WLLLMPPPPHHQQRRRRMMMMTTTTNNKKSSRRVVVMAAAADDEEGGGG") // 11
  ("FFLLSSSSYY**CC*WLLLMPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 12
  ("FFLMSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVMAAAADDEEGGGG") // 13
  ("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG") // 14
  ("FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 15
  ("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 16
  ("") // 17 deleted
  ("") // 18 deleted
  ("") // 19 deleted
  ("") // 20 deleted
  ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVMAAAADDEEGGGG") // 21
  ("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") // 22
  ("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRMIIMTTTTNNKKSSRRVVVMAAAADDEEGGGG") // 23
  ("FFLMSSSSYY**CCWWLLLMPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVMAAAADDEEGGGG") // 24
  ;
  
map<string,uint16_t> cReferenceSequences::codon_to_aa_index = make_map<string,uint16_t>
  ("TTT",  0)("TTC",  1)("TTA",  2)("TTG",  3)
  ("TCT",  4)("TCC",  5)("TCA",  6)("TCG",  7)
  ("TAT",  8)("TAC",  9)("TAA", 10)("TAG", 11)
  ("TGT", 12)("TGC", 13)("TGA", 14)("TGG", 15)
  ("CTT", 16)("CTC", 17)("CTA", 18)("CTG", 19)
  ("CCT", 20)("CCC", 21)("CCA", 22)("CCG", 23)
  ("CAT", 24)("CAC", 25)("CAA", 26)("CAG", 27)
  ("CGT", 28)("CGC", 29)("CGA", 30)("CGG", 31)
  ("ATT", 32)("ATC", 33)("ATA", 34)("ATG", 35)
  ("ACT", 36)("ACC", 37)("ACA", 38)("ACG", 39)
  ("AAT", 40)("AAC", 41)("AAA", 42)("AAG", 43)
  ("AGT", 44)("AGC", 45)("AGA", 46)("AGG", 47)
  ("GTT", 48)("GTC", 49)("GTA", 50)("GTG", 51)
  ("GCT", 52)("GCC", 53)("GCA", 54)("GCG", 55)
  ("GAT", 56)("GAC", 57)("GAA", 58)("GAG", 59)
  ("GGT", 60)("GGC", 61)("GGA", 62)("GGG", 63)
;

char cReferenceSequences::translate_codon(string seq, uint32_t translation_table, uint32_t codon_pos_1)
{  
  ASSERT(seq.size()==3, "Attempt to translate codon without three bases.");
  ASSERT(translation_table <= cReferenceSequences::translation_tables.size(), "Unknown translation table requested.");
  string& tt = (codon_pos_1 == 1) 
    ? cReferenceSequences::initiation_codon_translation_tables[translation_table]
    : cReferenceSequences::translation_tables[translation_table]
  ;
  ASSERT(tt.size() == 64, "Unknown translation table requested.");
  
  return (cReferenceSequences::codon_to_aa_index.count(seq) == 0) 
    ? '?' 
    : tt[cReferenceSequences::codon_to_aa_index[seq]];
}

void cReferenceSequences::annotate_1_mutation(cDiffEntry& mut, uint32_t start, uint32_t end, bool repeat_override, bool ignore_pseudogenes)
{
  // this could be moved to the object
  string intergenic_seperator = "/";

  // initialize everything, even though we don"t always use it
  mut["locus_tag"] = "";
  mut["aa_position"] = "";
  mut["aa_ref_seq"] = "";
  mut["aa_new_seq"] = "";
  mut["codon_position"] = "";
  mut["codon_ref_seq"] = "";
  mut["codon_new_seq"] = "";
  mut["gene_name"] = "";
  mut["gene_position"] = "";
  mut["gene_product"] = "";
  mut["_gene_product_hide"] = "";
  mut["gene_list"] = ""; //#affected genes

  string seq_id = mut["seq_id"];

  //or die "Unknown seq_id in reference sequence info: $seq_id\n";

  cSequenceFeatureList& gene_list_ref = (*this)[seq_id].m_genes;
  cSequenceFeatureList& repeat_list_ref = (*this)[seq_id].m_repeats;

  int32_t size = end - start + 1;

  Gene prev_gene, next_gene;
  vector<Gene> within_genes;
  vector<Gene> between_genes;
  vector<Gene> inside_left_genes;
  vector<Gene> inside_right_genes;

  cSequenceFeaturePtr repeat_ptr(NULL);
  if (repeat_override)
  {
    assert(start == end);
    repeat_ptr = get_overlapping_feature(repeat_list_ref, start);
    if (repeat_ptr.get() != NULL)
    {
      Gene within_gene(*repeat_ptr);
      within_genes.push_back(within_gene);
    }
  }

  if (repeat_ptr.get() == NULL)
    find_nearby_genes(gene_list_ref, start, end, within_genes, between_genes, inside_left_genes, inside_right_genes, prev_gene, next_gene);

  // Mutation is intergenic
  if (within_genes.size() + between_genes.size() + inside_left_genes.size() + inside_right_genes.size() == 0)
  {
    mut["snp_type"] = "intergenic";

    mut["gene_name"] += (prev_gene.name.size() > 0) ? prev_gene.name : "–"; //en-dash
    mut["gene_name"] += intergenic_seperator;
    mut["gene_name"] += (next_gene.name.size() > 0) ? next_gene.name : "–"; //en-dash

    mut["locus_tag"] += (prev_gene.get_locus_tag().size() > 0) ? prev_gene.get_locus_tag() : "–"; //en-dash
    mut["locus_tag"] += intergenic_seperator;
    mut["locus_tag"] += (next_gene.get_locus_tag().size() > 0) ? next_gene.get_locus_tag() : "–"; //en-dash

    if (prev_gene.name.size() > 0)
    {
      mut["gene_position"] += "intergenic (";
      mut["gene_position"] += prev_gene.is_top_strand() ? "+" : "-"; //hyphen
      mut["gene_position"] += to_string(start - prev_gene.get_end_1());
    }
    else
    {
      mut["gene_position"] += "intergenic (–";
    }
    mut["gene_position"] += intergenic_seperator;
    if (next_gene.name.size() > 0)
    {
      mut["gene_position"] += next_gene.is_top_strand() ? "-" : "+"; //hyphen
      mut["gene_position"] += to_string(next_gene.get_start_1() - end);
    }
    else
    {
      mut["gene_position"] += "–"; //en-dash
    }
    mut["gene_position"] += ")";

    mut["gene_product"] += (prev_gene.name.size() > 0) ? prev_gene.product : "–"; //en-dash
    mut["gene_product"] += intergenic_seperator;
    mut["gene_product"] += (next_gene.name.size() > 0) ? next_gene.product : "–"; //en-dash

    return;
  }
  // Mutation is completely within genes
  else if (within_genes.size() > 0)
  {
    /// TODO: It can be within multiple genes, in which case we need to annotate
    /// the change it causes in each reading frame UGH! YUCKY!
    /// FOR NOW: just take the first of the within genes...
    Gene gene = within_genes[0];
    mut["gene_name"] = gene.name;
    mut["gene_product"] = gene.product;
    if (gene.get_locus_tag().size()) {
      mut["locus_tag"] = gene.get_locus_tag();
    }

    //#added for gene table
    mut["gene_list"] = gene.name;

    int32_t within_gene_start = gene.is_top_strand() ? gene.get_start_1() : gene.get_end_1();

    if (start == end)
    {
      mut["gene_position"] = to_string(abs(static_cast<int32_t>(start) - within_gene_start) + 1);
    }
    else
    {
      uint32_t gene_start = abs(static_cast<int32_t>(start) - within_gene_start) + 1;
      uint32_t gene_end = abs(static_cast<int32_t>(end) - within_gene_start) + 1;
      mut["gene_position"] = (gene_start < gene_end)
        ? to_string(gene_start) + "-" + to_string(gene_end)  //hyphen
        : to_string(gene_end) + "-" + to_string(gene_start); //hyphen
    }

    string gene_nt_size = to_string(gene.get_end_1() - gene.get_start_1() + 1);
    string gene_strand_char = gene.is_top_strand() ? ">" : "<"; 
    mut["gene_strand"] = gene_strand_char;
    
    // ...but the gene is a pseudogene or not a protein coding gene
    if (!ignore_pseudogenes && gene.pseudogene)
    {
      mut["snp_type"] = "pseudogene";
      mut["gene_position"] = "pseudogene (" + mut["gene_position"] + "/" + gene_nt_size + " nt)";
      return;
    }
    else if (gene.type != "CDS")
    {
      mut["snp_type"] = "noncoding";
      mut["gene_position"] = "noncoding (" + mut["gene_position"] + "/" + gene_nt_size + " nt)";
      return;
    }

    //#only add gene information to SNPs and RA mutations that don"t include indels...
    if ((mut._type != SNP) && !((mut._type == RA) && (mut["ref_base"] != ".") && (mut["new_base"] != ".")))
    {
      mut["gene_position"] = "coding (" + mut["gene_position"] + "/" + gene_nt_size + " nt)";
      return;
    }

    // this is for RA...
    if (!mut.entry_exists("ref_seq")) mut["ref_seq"] = mut["ref_base"];
    if (!mut.entry_exists("new_seq")) mut["new_seq"] = mut["new_base"];

    // determine the old and new translation of this codon
    mut["aa_position"] = to_string((from_string<uint32_t>(mut["gene_position"]) - 1) / 3 + 1); // 1 indexed
    uint32_t codon_pos_1 = int(abs(static_cast<int32_t>(start) - within_gene_start)) % 3 + 1;
    mut["codon_position"] = to_string<int32_t>(codon_pos_1); // 1 indexed
    
    string& ref_string = (*this)[seq_id].m_fasta_sequence.m_sequence;
    string codon_seq = gene.is_top_strand()
      ? ref_string.substr(gene.get_start_1() + 3 * (from_string<uint32_t>(mut["aa_position"]) - 1) - 1, 3)
      : reverse_complement(ref_string.substr(gene.get_end_1() - 3 * from_string<uint32_t>(mut["aa_position"]), 3));

    mut["codon_ref_seq"] = codon_seq;
    mut["aa_ref_seq"] = translate_codon(mut["codon_ref_seq"], gene.translation_table, codon_pos_1);

    mut["codon_new_seq"] = codon_seq;
    //#remember to revcom the change if gene is on opposite strand
    mut["codon_new_seq"][from_string<uint32_t>(mut["codon_position"]) - 1] = gene.is_top_strand() 
      ? mut["new_seq"][0]
      : reverse_complement(mut["new_seq"])[0];
    mut["aa_new_seq"] =  translate_codon(mut["codon_new_seq"], gene.translation_table, codon_pos_1);

    mut["snp_type"] = (mut["aa_ref_seq"] != mut["aa_new_seq"]) ? "nonsynonymous" : "synonymous";
  }

  //The mutation actually contains several genes
  else if (between_genes.size() + inside_left_genes.size() + inside_right_genes.size() > 0)
  {
    vector<Gene> gene_list;
    gene_list.insert( gene_list.end(), inside_left_genes.begin(), inside_left_genes.end() );
    gene_list.insert( gene_list.end(), between_genes.begin(), between_genes.end() );
    gene_list.insert( gene_list.end(), inside_right_genes.begin(), inside_right_genes.end() );

    vector<string> gene_name_list;
    vector<string> locus_tag_list;
    for (vector<Gene>::iterator it=inside_left_genes.begin(); it != inside_left_genes.end(); it++)
    {
      gene_name_list.push_back("[" + it->name + "]");

      if (it->get_locus_tag().size()) {
        locus_tag_list.push_back("[" + it->get_locus_tag() + "]");
      }

    }
    for (vector<Gene>::iterator it=between_genes.begin(); it != between_genes.end(); it++)
    {
      gene_name_list.push_back(it->name);

      if (it->get_locus_tag().size()) {
        locus_tag_list.push_back("[" + it->get_locus_tag() + "]");
      }

    }
    for (vector<Gene>::iterator it=inside_right_genes.begin(); it != inside_right_genes.end(); it++)
    {
      gene_name_list.push_back("[" + it->name + "]");

      if (it->get_locus_tag().size()) {
        locus_tag_list.push_back("[" + it->get_locus_tag() + "]");
      }

    }
    
    // We ended up calling this function a lot.
    string sJoinedGeneList = join(gene_name_list, ", ");
    
    mut["gene_product"] += sJoinedGeneList;
    mut["_gene_product_hide"] = "<i>" + mut["gene_product"] + "</i>";
    
    // @MDS0004 - This number here will condense the number of genes listed.
    // Change to show how many appear in HTML. This is a javascript call.
    //
    // Non-Button Toggle - Replace <input>
    // <a href=\"javascript:hideTog('" + sID + "');\">Toggle Genes</a>
    if(gene_list.size() > 15)
    {
      string sID = "gene_hide_" + to_string(mut._type) + "_" + mut._id;
      string sButtonID = sID + "_button";
      mut["_gene_product_hide"] =
        "<i title=\"" + sJoinedGeneList + "\"><b>" +
        to_string(gene_list.size()) +
        " genes</b> <noscript>Javascript Disabled: All genes shown.<br></noscript></i>"
        + "<noscript>" + mut["gene_product"] + "</noscript>"
        + "<div id=\"" + sID + "\" class=\"hidden\">"
        + mut["gene_product"] + "</div>"
        + "<input id=\"" + sButtonID + "\" type=\"button\" onclick=\"hideTog('" + sID + "');showTog('" + sButtonID + "')\" value=\"Show\" />";
      mut["gene_product"] = "<i>";
      mut["gene_product"] += "<b>" + to_string(gene_list.size()) + " genes</b><BR>";
      mut["gene_product"] += sJoinedGeneList;
      mut["gene_product"] += "</i>";
    }

    if (gene_name_list.size() == 1) {
      mut["gene_name"] = gene_name_list[0];
    } else {
      mut["gene_name"] = gene_name_list.front() + "–" + gene_name_list.back();  //en-dash
    }

    if (locus_tag_list.size() == 1) {
      mut["locus_tag"] = locus_tag_list[0];
    }
    else if (locus_tag_list.size() > 1) {
      mut["locus_tag"] = locus_tag_list.front() + "–" + locus_tag_list.back();  //en-dash
    }

  }
}

void cReferenceSequences::annotate_mutations(cGenomeDiff& gd, bool only_muts, bool ignore_pseudogenes)
{
  //keep track of other mutations that affect SNPs
  //because we may double-hit a codon

  //TODO: the proper way to do this is to create list of SNPs that have been hit
  // hashed by gene protein accession ID and AA position within gene
  // and have the annotation point to them (and back at them)
  // so that the codon will be correctly updated with all changes and we can notify the
  // changes that their SNP_type is not really SNP, but multiple hit SNP.
  
  vector<cDiffEntry*> snp_muts;

  diff_entry_list_t muts = gd.show_list();
  for (diff_entry_list_t::iterator it=muts.begin(); it!=muts.end(); it++)
  {
    cDiffEntry& mut= **it;

    if (only_muts && !(mut.is_mutation())) continue;

    switch (mut._type)
    {
      case SNP:{
        mut["_ref_seq"] = get_sequence_1(mut["seq_id"], from_string<uint32_t>(mut["position"]), from_string<int32_t>(mut["position"]));
        annotate_1_mutation(mut, from_string<int32_t>(mut["position"]), from_string<int32_t>(mut["position"]), false, ignore_pseudogenes);
        snp_muts.push_back(&mut);
      } break;
        
      case SUB:{
        annotate_1_mutation(mut, from_string<int32_t>(mut["position"]), from_string<int32_t>(mut["position"]) + from_string<int32_t>(mut["size"]) - 1);
      } break;
        
      case DEL:{
        annotate_1_mutation(mut, from_string<int32_t>(mut["position"]), from_string<int32_t>(mut["position"]) + from_string<int32_t>(mut["size"]) - 1);
      } break;
        
      case INS:{
        annotate_1_mutation(mut, from_string<int32_t>(mut["position"]), from_string<int32_t>(mut["position"]));
      } break;
        
      case CON:{
        annotate_1_mutation(mut, from_string<int32_t>(mut["position"]), from_string<int32_t>(mut["position"]) + from_string<int32_t>(mut["size"]) - 1);
      } break;
        
      case MOB:{
        annotate_1_mutation(mut, from_string<int32_t>(mut["position"]), from_string<int32_t>(mut["position"]) + from_string<int32_t>(mut["duplication_size"]) - 1);
      } break;
        
      case INV:{
        annotate_1_mutation(mut, from_string<int32_t>(mut["position"]), from_string<int32_t>(mut["position"]));
        mut["gene_name_1"] = mut["gene_name"];
        mut["gene_product_1"] = mut["gene_product"];
        annotate_1_mutation(mut, from_string<int32_t>(mut["position"]) + from_string<int32_t>(mut["size"])-1, from_string<int32_t>(mut["position"]) + from_string<int32_t>(mut["size"])-1);
        mut["gene_name_2"] = mut["gene_name"];
        mut["gene_product_2"] = mut["gene_product"];
        mut.erase("gene_name");
        mut.erase("gene_product");
      } break;
        
      case AMP:{
        annotate_1_mutation(mut, from_string<int32_t>(mut["position"]), from_string<int32_t>(mut["position"]) + from_string<int32_t>(mut["size"]) - 1);
      } break;
        
      case JC:{
        cDiffEntry side_1;
        side_1[SEQ_ID] = mut["side_1_seq_id"];
        annotate_1_mutation(side_1, from_string<int32_t>(mut["side_1_position"]), from_string<int32_t>(mut["side_1_position"]), true);
        //copy over entries with prefix
        for(diff_entry_map_t::iterator it=side_1.begin(); it!=side_1.end(); it++)
        {
          mut["_side_1"+ it->first] = it->second;
        }
        
        cDiffEntry side_2;
        side_2[SEQ_ID] = mut["side_2_seq_id"];
        annotate_1_mutation(side_2, from_string<int32_t>(mut["side_2_position"]), from_string<int32_t>(mut["side_2_position"]), true);
        //copy over entries with prefix
        for(diff_entry_map_t::iterator it=side_2.begin(); it!=side_2.end(); it++)
        {
          mut["_side_2"+ it->first] = it->second;
        }
      } break;
        
      case RA:{
        annotate_1_mutation(mut, from_string<int32_t>(mut["position"]), from_string<int32_t>(mut["position"]));
      } break;
        
      case MC:{
        annotate_1_mutation(mut, from_string<int32_t>(mut["start"]), from_string<int32_t>(mut["end"]));
      } break;
        
      case CN:{
        annotate_1_mutation(mut, from_string<int32_t>(mut["start"]), from_string<int32_t>(mut["end"]));
      } break;
        
      default:{
      } break;
    }
  }
  
  //scan snps to see if they affect the same codon
    for (size_t i = 0; i < snp_muts.size(); i++){
      for (size_t j = i + 1; j < snp_muts.size(); j++){
        if ((*snp_muts[i])["codon_position"] == "" ||
            (*snp_muts[j])["codon_position"] == ""
            ){
          continue;
        }
        
        int32_t i_position = from_string<int32_t>((*snp_muts[i])["position"]);
        int32_t j_position = from_string<int32_t>((*snp_muts[j])["position"]);
        int32_t snp_distance = i_position - j_position;
        
        int32_t i_codon_position = from_string<int32_t>((*snp_muts[i])["codon_position"]);
        int32_t j_codon_position = from_string<int32_t>((*snp_muts[j])["codon_position"]);
        
        if (i_position > j_position){
          if (snp_distance <= 2 && j_codon_position < i_codon_position){
            string new_codon = (*snp_muts[i])["codon_new_seq"];
            const char new_char = (*snp_muts[j])["new_seq"][0];
            new_codon[j_codon_position - 1] = new_char;
            (*snp_muts[i])["codon_new_seq"] = new_codon;
            (*snp_muts[j])["codon_new_seq"] = new_codon;
          }
        }
        else{
          if (snp_distance >= -2 && j_codon_position > i_codon_position){
            string new_codon = (*snp_muts[i])["codon_new_seq"];
            const char new_char = (*snp_muts[j])["new_seq"][0];
            new_codon[j_codon_position - 1] = new_char;
            (*snp_muts[i])["codon_new_seq"] = new_codon;
            (*snp_muts[j])["codon_new_seq"] = new_codon;
          }
        }
      }
    }//for
}

void cReferenceSequences::polymorphism_statistics(Settings& settings, Summary& summary)
{
  string reference_fasta_file_name = settings.reference_fasta_file_name;
  vector<string> seq_ids = this->seq_ids();

  // some local variable lookups for convenience
  double log10_ref_length = log(this->total_length()) / log(10);

  //
  // Replacement for below
  //
  // ToDo: This should really make a different column for each input read set.
  //
  string coverage_fn = settings.file_name(settings.unique_only_coverage_distribution_file_name, "@", "");
  string outputdir = dirname(coverage_fn) + "/";
  string count_file_name = outputdir + "error_counts.tab";

  ifstream COUNT(count_file_name.c_str());
  assert(COUNT.is_open());
  string count_header_line, count_covariate_line;
  getline(COUNT, count_covariate_line); // ignore the first line

  // we parse the covariates to know how many different qualities we are dealing with...
  covariates_used_t         covariate_used;
  covariates_max_t          covariate_max;
  covariates_enforce_max_t  covariate_enforce_max;
  covariates_offset_t       covariate_offset;
  bool                      per_position;
  cErrorTable::read_covariates(count_covariate_line, covariate_used, covariate_max, covariate_enforce_max, covariate_offset, per_position);

  getline(COUNT, count_header_line); // second line is the header

  vector<string> count_header_list = split(count_header_line, "\t");

  uint32_t count_column = UNDEFINED_UINT32;
  uint32_t quality_column = UNDEFINED_UINT32;
  for (uint32_t i = 0; i < count_header_list.size(); i++)
  {
    if (count_header_list[i] == "quality")
      quality_column = i;
    else if (count_header_list[i] == "count")
      count_column = i;
  }

  ASSERT( (quality_column != UNDEFINED_UINT32) && (count_column != UNDEFINED_UINT32),
          "'quality' and 'count' columns not found in file: " + count_file_name);

  vector<uint32_t> quality_count_list(covariate_max[k_quality]);
  string line;
  while (COUNT.good())
  {
    getline(COUNT, line);
    vector<string> line_list = split(line,  "\t");
    if (line_list.size() < 2) break; // empty line
    uint32_t count = from_string<uint32_t>(line_list[count_column]);
    uint32_t quality = from_string<uint32_t>(line_list[quality_column]);
    quality_count_list[quality] += count;
  }
  COUNT.close();

  string genome_error_counts_file_name = settings.genome_error_counts_file_name;

  ofstream GEC(genome_error_counts_file_name.c_str());
  assert(GEC.is_open());
  for (uint32_t i = 1; i < quality_count_list.size(); i++)
  {
    uint32_t val = 0;
    val = quality_count_list[i];
    GEC << val << endl;
  }
  GEC.close();

  string polymorphism_statistics_input_file_name = settings.polymorphism_statistics_input_file_name;
  string polymorphism_statistics_output_file_name = settings.polymorphism_statistics_output_file_name;

  /// Load the older GenomeDiff and add new fields
  string ra_mc_genome_diff_file_name = settings.ra_mc_genome_diff_file_name;
  cGenomeDiff gd(ra_mc_genome_diff_file_name);

  string polymorphism_statistics_r_script_file_name = settings.polymorphism_statistics_r_script_file_name;
  string polymorphism_statistics_r_script_log_file_name = settings.polymorphism_statistics_r_script_log_file_name;
  uint32_t total_reference_length = summary.sequence_conversion.total_reference_sequence_length;

  string command = "R --vanilla total_length=" + to_string<uint32_t>(total_reference_length) +
    " in_file="   + cString(polymorphism_statistics_input_file_name).escape_shell_chars() +
    " out_file="  + cString(polymorphism_statistics_output_file_name).escape_shell_chars() +
    " qual_file=" + cString(genome_error_counts_file_name).escape_shell_chars() +
    " < "         + cString(polymorphism_statistics_r_script_file_name).escape_shell_chars() +
    " > "         + cString(polymorphism_statistics_r_script_log_file_name).escape_shell_chars();
  SYSTEM(command, false, false, false); //NOTE: Not escaping shell characters here.

  // Read R file and add new results corresponding to all columns
  ifstream ROUT(polymorphism_statistics_output_file_name.c_str());
  assert(ROUT.is_open()); // or die "Could not find file: $polymorphism_statistics_output_file_name";
  string header;
  getline(ROUT, header);
  vector<string> header_list = split(header, "\t");

  cGenomeDiff new_gd;

  diff_entry_list_t muts = gd.evidence_list();
  for (diff_entry_list_t::iterator it=muts.begin(); it!=muts.end(); it++)
  {
    cDiffEntry& mut= **it;

    // lines only exist for RA evidence
    if (mut._type != RA)
    {
      new_gd.add(mut);
      continue;
    }

    // lines only exist for polymorphisms
    if ((mut[FREQUENCY] == "1") || (mut[FREQUENCY] == "0"))
    {
      new_gd.add(mut);
      continue;
    }

    string line;
    getline(ROUT, line);
    vector<string> line_list = split(line, "\t");

    for (uint32_t j = 0; j < header_list.size(); j++)
    {
      assert(line_list.size() > j); // die "Incorrect number of items on line:\n$line"
      mut[header_list[j]] = line_list[j];
    }

    // Evalue cutoff again (in case we are only running this part)
    if (double_from_string(mut[POLYMORPHISM_QUALITY]) < settings.polymorphism_log10_e_value_cutoff)
      add_reject_reason(mut, "EVALUE");

    // Frequency cutoff
    if ( (from_string<double>(mut["frequency"]) < settings.polymorphism_frequency_cutoff)
      || (from_string<double>(mut["frequency"]) > 1-settings.polymorphism_frequency_cutoff) )
      add_reject_reason(mut, "POLYMORPHISM_FREQUENCY_CUTOFF");

    // Minimum coverage on both strands
    double polymorphism_coverage_limit_both_bases = settings.polymorphism_minimum_new_coverage_each_strand;
    bool passed = true;
    vector<string> top_bot = split(mut["ref_cov"], "/");
    double top = from_string<double>(top_bot[0]);
    double bot = from_string<double>(top_bot[1]);
    passed = passed && (top >= polymorphism_coverage_limit_both_bases);
    passed = passed && (bot >= polymorphism_coverage_limit_both_bases);

    // Test strand criterion
    top_bot = split(mut["new_cov"], "/");
    top = from_string<double>(top_bot[0]);
    bot = from_string<double>(top_bot[1]);
    passed = passed && (top >= polymorphism_coverage_limit_both_bases);
    passed = passed && (bot >= polymorphism_coverage_limit_both_bases);
    if (!passed)
      add_reject_reason(mut, "POLYMORPHISM_STRAND");
    
    // Test bias criteria
    if (settings.polymorphism_bias_p_value_cutoff && (from_string<double>(mut["ks_quality_p_value"]) < settings.polymorphism_bias_p_value_cutoff))
      add_reject_reason(mut, "KS_QUALITY_P_VALUE");
    if (settings.polymorphism_bias_p_value_cutoff && (from_string<double>(mut["fisher_strand_p_value"]) < settings.polymorphism_bias_p_value_cutoff))
      add_reject_reason(mut, "FISHER_STRAND_P_VALUE");

    ////// Optionally, ignore if in a homopolymer stretch longer than this
    if (settings.polymorphism_reject_homopolymer_length)
    {
      uint32_t test_length = 20;
      string seq_id = mut["seq_id"];
      uint32_t end_pos = from_string<uint32_t>(mut["position"]);
      uint32_t start_pos = end_pos - test_length + 1;
      if (start_pos < 1) start_pos = 1;
      string bases = this->get_sequence_1(seq_id, start_pos, end_pos);

      uint32_t same_base_length = 0;
      string first_base = bases.substr(end_pos - start_pos, 1);
      for (uint32_t j = end_pos; j >= start_pos; j--)
      {
        string this_base = bases.substr(j - start_pos, 1);
        if (first_base != this_base) break;
        same_base_length++;
      }

      if (same_base_length >= settings.polymorphism_reject_homopolymer_length)
      {
        add_reject_reason(mut, "HOMOPOLYMER_STRETCH");
      }
    }

    if (
        mut.number_reject_reasons() > 0
        && (double_from_string(mut[POLYMORPHISM_QUALITY]) > settings.mutation_log10_e_value_cutoff)
        && (from_string<double>(mut[FREQUENCY]) > 0.5)
        )
    {
      mut[FREQUENCY] = "1";
      mut.erase(REJECT);

      // FIX -- need to re-evaluate whether it would have been accepted as a normal mutation
      // This is NOT the right quality being used here. Need a separate quality for consensus call and polymorphism call!
      if (double_from_string(mut[POLYMORPHISM_QUALITY]) < settings.mutation_log10_e_value_cutoff)
        add_reject_reason(mut, "EVALUE");
    }

    new_gd.add(mut);

  }

  ROUT.close();

  /// Write out the file which now has much more data
  string polymorphism_statistics_ra_mc_genome_diff_file_name = settings.polymorphism_statistics_ra_mc_genome_diff_file_name;
  new_gd.write(polymorphism_statistics_ra_mc_genome_diff_file_name);

}

std::string GetWord(std::string &s) {
  RemoveLeadingWhitespace(s);
  int found = s.find_first_of(" =\n\r");
  std::string w = s.substr(0,found);
  s.erase(0, found);
  return w;
}


void RemoveLeadingWhitespace(std::string &s) {
  int found = s.find_first_not_of(" \t");
  s.erase(0,found);
}

void RemoveLeadingTrailingWhitespace(std::string &s) {
  int found = s.find_first_not_of(" \t\n\r");
  s.erase(0,found);
  found = s.find_last_not_of(" \t\n\r");
  s.erase(found+1,s.length());
}

  
/*
 Count the number of nonmatching bases in the alignment between the read and reference sequences
   Unmatched (padded) bases at the end of a read count as mismatches.
 */
uint32_t alignment_mismatches(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info)
{
  bool verbose = false;
  uint32_t mismatches = 0;

  uint32_t index = a.reference_target_id();
  const string& const_ref_string = ref_seq_info[index].m_fasta_sequence.m_sequence;
  string ref_string = const_ref_string.substr(a.reference_start_0(), a.reference_match_length());
  uint32_t ref_pos = 0;

  string read_string = a.read_char_sequence().substr(a.query_start_0(), a.query_match_length());
  uint32_t read_pos = 0;

  uint32_t* cigar_list = a.cigar_array(); // cigar array for this alignment

  if (verbose)
  {
    cout << a.read_name() << endl;
  }

  for (uint32_t i = 0; i < a.cigar_array_length(); i++)
  {
    char op = cigar_list[i] & BAM_CIGAR_MASK;
    uint32_t len = cigar_list[i] >> BAM_CIGAR_SHIFT;

    // soft padding counts as a mismatch
    if (op == BAM_CSOFT_CLIP)
    {
      mismatches += len;
    }
    else if (op == BAM_CDEL)
    {
      mismatches += len;
      ref_pos += len;
    }
    else if (op == BAM_CINS)
    {
      mismatches += len;
      read_pos += len;
    }
    else if (op == BAM_CMATCH)
    {
      for (uint32_t j = 0; j < len; j++)
      {
        if (verbose)
        {
          cout << "REF: " << ref_pos << "  " << ref_string[ref_pos] << endl;
          cout << "READ: " << read_pos << "  " << read_string[read_pos] << endl;
        }
        if (ref_string[ref_pos] != read_string[read_pos])
        {
          mismatches++;
        }
        read_pos++;
        ref_pos++;
      }
    }
  }

  return mismatches;
}

string shifted_cigar_string(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info)
{
  bool verbose = false;

  string ref_seq = ref_seq_info[a.reference_target_id()].m_fasta_sequence.m_sequence;
  uint32_t ref_seq_index = 0;
  string read_seq = a.read_char_sequence();
  uint32_t read_seq_index = 0;
  vector<pair<char,uint16_t> > cigar_pair_array = a.cigar_pair_array();

  if (verbose)
  {
    cout << a.read_name() << endl;

    if (a.read_name() == "GW1ULQG02EIUY7")
    {
      cout << "debug" << endl;
    }
    cout << a.cigar_string() << endl;
  }

  for (vector<pair<char,uint16_t> >::iterator it = cigar_pair_array.begin(); it != cigar_pair_array.end(); it++)
  {
    char op = it->first;
    uint16_t len = it->second;

    if (op == 'D')
    {
      char base = ref_seq[ref_seq_index];
      bool all_same_base = true;
      for (uint32_t j = 1; j < len; j++)
        all_same_base = all_same_base && (ref_seq[ref_seq_index + j] == base);

      if (all_same_base)
      {
        uint32_t shift_amount = 0;
        while (ref_seq[ref_seq_index + len + shift_amount] == base)
          shift_amount++;

        if (shift_amount > 0)
        {
          if (verbose)
          {
            cout << "Shifting D: " << a.read_name() << endl;
          }

          if (it != cigar_pair_array.begin())
            (it - 1)->second += shift_amount;
          else
            cigar_pair_array.insert(it, make_pair<char,uint16_t>('M', shift_amount));

          if (it != cigar_pair_array.end())
            (it + 1)->second -= shift_amount;
          else
            cigar_pair_array.push_back(make_pair<char,uint16_t>('M', shift_amount));
        }

        ref_seq_index += len + shift_amount;
      }
      else
      {
        ref_seq_index += len;
      }
    }
    else if (op == 'I')
    {
      char base = read_seq[read_seq_index];
      bool all_same_base = true;
      for (uint32_t j = 1; j < len; j++)
        all_same_base = all_same_base && (read_seq[read_seq_index + j] == base);

      if (all_same_base)
      {
        uint32_t shift_amount = 0;
        while (read_seq[read_seq_index + len + shift_amount] == base)
          shift_amount++;

        if (shift_amount > 0)
        {
          if (verbose)
          {
            cout << "Shifting I: " << a.read_name() << endl;
          }

          if (it != cigar_pair_array.begin())
            (it - 1)->second += shift_amount;
          else
            cigar_pair_array.insert(it, make_pair<char,uint16_t>('M', shift_amount));

          if (it != cigar_pair_array.end())
            (it + 1)->second -= shift_amount;
          else
            cigar_pair_array.push_back(make_pair<char,uint16_t>('M', shift_amount));
        }

        read_seq_index += len + shift_amount;
      }
      else
      {
        read_seq_index += len;
      }
    }
    else
    {
      ref_seq_index += len;
      read_seq_index += len;
    }
  }

  stringstream shifted_cigar_string_ss;
  for (vector<pair<char,uint16_t> >::iterator it = cigar_pair_array.begin(); it != cigar_pair_array.end(); it++)
    shifted_cigar_string_ss << it->second << it->first;

  string shifted_cigar_string = shifted_cigar_string_ss.str();
  if (verbose)
  {
    cout << shifted_cigar_string << endl;
  }

  return shifted_cigar_string;
}


} // breseq namespace

