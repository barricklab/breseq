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

#include "libbreseq/contingency_loci.h"

using namespace std;

namespace breseq {
  
/*! Convenience wrapper around the identify_mutations_pileup class.
 */
// The file "loci" is a white space separated file where each line documents a homopolymer repeat. The second column must be the coordinate of the first base of the repeat, and the fifth column must be the name of the gene that the repeat is located in.
// The output file is a space separated file where each line contains the counts of the reads of each length at each locus. The last value of the line is the coordinate of the repeat.
// Strict is 1 when the program ignores reads that do not match exactly 5 bases around the homopolymer repeat, and 0 when the program does not ignore them.  
void analyze_contingency_loci(
                              const string& bam,
                              const string& fasta,
                              const vector<string>& ref_seq_file_names,
                              const string& output,
                              const string& loci,
                              int strict
                              ) {
  
  cout << "Loading reference sequences..." << endl;
  
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(ref_seq_file_names);

  cout << "Finding repeats..." << endl;
  homopolymer_repeat_list hr;
  identify_homopolymer_repeats(hr, ref_seq_info);  
  contingency_loci_pileup clp(bam,fasta, loci, strict);
  
  cout << "Analyzing alignments..." << endl;
  // For each repeat
  for( size_t i=0; i<hr.size(); i++ ){
    
    // Makes the region string
    stringstream ss; 
    ss << hr[i].start << "-" << hr[i].start + hr[i].length - 1; 
    
    // string region= "NC_012660:456-466:T";
    string region = hr[i].seq_id + ":" + ss.str() + ":" + hr[i].base;
    
    clp.analyze_contingency_locus(region);
  }
      
  // Output the histograms
  clp.printStats( output, ref_seq_info );
  
}
    
        
    
void analyze_contingency_loci_significance(const string& output, const vector<string>& strainfiles){
    
    vector<string> locifiles;
    
    for( unsigned int i=0; i<strainfiles.size(); i++ ){
        
        ifstream in( strainfiles[i].c_str() );
        string line;
        getline(in,line);
        vector<string> tokens = split( line, "\t");
        
        int seqid_index = 0;
        while ( strcmp(tokens[seqid_index].c_str(),"seq_id") != 0 ) {
            seqid_index += 1;
        }
        //
        while( !in.eof() ){
            FILE* locus = fopen( (tokens[seqid_index] + tokens[seqid_index+1]).c_str(), "ap" );
            
            // If the locus file has not yet already been accounted for, this will add it to the locifiles list
            bool included = false;
            for(unsigned long j=0; j<locifiles.size(); j++ ){
                if( strcmp( locifiles[j].c_str(), (tokens[seqid_index] + tokens[seqid_index+1]).c_str() ) == 0 ){
                    included = true;
                }
            }
            if(included == false){
                locifiles.push_back( (tokens[seqid_index] + tokens[seqid_index+1]).c_str() );
                fprintf(locus,"sample\tlength\tcount");
            }
            
            // Add the information from strain strainfiles[i] to the locus file
            getline(in,line);
            tokens = split( line, "\t");
            for(int j=0; j<seqid_index; j++ ){
                fprintf( locus, "%u\t%i\t%s", i, j+1, tokens[j].c_str() );
            }
            fclose(locus);
        }    
    }
    
    // Run the R script to calculate the p-value for each locus, then print them
    FILE* out = fopen( output.c_str(), "w");
    for( unsigned int i=0; i<locifiles.size(); i++ ){
        
        SYSTEM( "cp " + cString(locifiles[i]) + " example_data.txt" );
        SYSTEM( "Rscript R_fitting_commands.txt" );
      
        //@JEB this is not the right script -- some more integration needs to be done
      
        ifstream ROUT( "results.txt");
        string line;
        getline(ROUT,line);
        float p_value = atof(line.c_str());
        if( p_value < 0.5/locifiles.size() ){
            fprintf( out, "%s %f", locifiles[i].c_str(), p_value );
        }
    }
    
}



void identify_homopolymer_repeats(homopolymer_repeat_list& hr, const cReferenceSequences& ref_seqs)
{
  //For each sequence
  for( size_t i=0; i<ref_seqs.size(); i++ ){
    
    char base = ' ';
    int rnumber = 0;
    
    // For each nucleotide in the sequence
    for( size_t j=0; j<ref_seqs[i].m_fasta_sequence.m_sequence.size(); j++ ){
      if( ref_seqs[i].m_fasta_sequence.m_sequence[j] == base ){
        rnumber++;
      }
      else{
        if( (base != 'N') && (rnumber >= 8) ){
          homopolymer_repeat r;
          r.seq_id = ref_seqs[i].m_fasta_sequence.m_name;
          r.start = j+1-rnumber;
          r.base = base;
          r.length = rnumber;
          hr.push_back(r);
          //printf( "%i:%c\n", r.start, r.base );
        }
        base = ref_seqs[i].m_fasta_sequence.m_sequence[j];
        rnumber = 1;
      }
    }
    // This checks if the last nucleotide was part of a repeat
    if( rnumber >= 8 ){
      homopolymer_repeat r;
      r.seq_id = ref_seqs[i].m_fasta_sequence.m_name;
      r.start = ref_seqs[i].m_fasta_sequence.m_sequence.size()+1-rnumber;
      r.base = base;
      r.length = rnumber;
      hr.push_back(r);
      //printf( "%i:%c/n", r.start, r.base );
      
    }
    
  }
}


/*! Constructor.
 */
contingency_loci_pileup::contingency_loci_pileup(
                                                 const string& bam,
                                                 const string& fasta,
                                                 const string& loci,
                                                 int s
                                                 )
: pileup_base(bam, fasta), strict(s)

{
  cout << "Reading indices...\n";
  readIndices(indices, names ,loci);
  
  //cout << "Filename: " << split( bam, "." )[0].append(".tam") << "\n";
  //tf.open_write( split( bam, "." )[0].append(".tam"), fasta );
  set_print_progress(true);
}

// The string "region" must have fields separated by ':', where the first field is the name of the genome, the second field is the range of coordinates spanned by the region, and the third field is the base that is repeated. EG: "NC_012660:456-466:T"
void contingency_loci_pileup::analyze_contingency_locus(const string& region) {
  
  //cout << "========================" << endl;
  //cout << region << "\n";
  //cout << "========================" << endl;

  vector<string> r = split( region, ":");
  vector<string> coords = split( r[1], "-" );
  
  //cout << "BASE: " << r[2] << "\n";
  current_region.base = r[2][0];
  current_region.start = atoi( coords[0].c_str() );
  current_region.length = atoi( coords[1].c_str() ) - current_region.start + 1;
  
  // Proceeds if it is a contingecy locus, or there the file for the loci "NA" was given
  bool contingency = 0;
  for( size_t i=0; i<indices.size(); i++ ){
    if( static_cast<uint32_t>(indices[i]) == current_region.start ){
      //cout << "locus: " << indices[i] << "\n";
      //getchar();
      contingency = true;
    }
  }
  if( contingency == false && indices.size() != 0 ){
    return;
  }
  
  // Creates a repeats_stats object and pushes it into the class's array of
  repeats.push_back( repeat_stats( region ) );
  do_fetch( r[0].append(":").append(r[1]) );
  
  // Constructs histogram from the values added from fetch_callback
  vector<double> histogram(1,0);
  int count = 0;
  //printf( "Freqs.size(): %i\n", repeats.back().freqs.size() );
  for( size_t i=0; i<repeats.back().freqs.size(); i++ ){
    //printf( "%f ", repeats.back().freqs[i] );
    if( histogram.size() < static_cast<size_t>(repeats.back().freqs[i])+1 ){
      histogram.resize( static_cast<size_t>(repeats.back().freqs[i])+1 );
    }
    histogram[ int(repeats.back().freqs[i]) ]++;
  }
  //Normalizes it
  /*
   for( int i=0; i<histogram.size(); i++ ){
   histogram[i] = histogram[i]/count;
   }*/
  repeats.back().freqs = histogram;
  
}

/*! Called for each alignment.
 Counts the lengths of the repeats in each of the appropriate alignments and stores it into the appropriate repeats_stats in the array
 */ 
void contingency_loci_pileup::fetch_callback(const alignment_wrapper& a) {
  
  bool verbose = false;
  vector<pair<char,uint16_t> > cigar_pair = a.cigar_pair_array();
  int dist_to_first_match = 0;
  
  // optional guards
  //if (!a.proper_pair()) return;
  
  if (a.is_redundant()) return;

  // unwrap the whole alignment into a couple of strings with a vector of the reference positions
  // this code is very much like alignment_output

  
  string read_sequence = a.read_char_sequence();
  int32_t read_length = a.read_length();
  
  const char* ref_sequence = get_refseq(a.reference_target_id());
  int32_t ref_length = target_length(a.reference_target_id());

  // decide which side is better anchored with longer match
  uint32_t u_ref_start_1, u_ref_end_1;
  a.reference_bounds_1(u_ref_start_1, u_ref_end_1);
  int32_t ref_start_1 = static_cast<int32_t>(u_ref_start_1);
  int32_t ref_end_1 = static_cast<int32_t>(u_ref_end_1);
  
  int32_t homopolymer_start_1 = current_region.start;
  int32_t homopolymer_end_1 = current_region.start + current_region.length -1;
  
  //if (a.read_name() == "1:1841572")
  //  cout << "test";
  
  // there is more length at the start of the read or at the end?
  if ( homopolymer_start_1 - ref_start_1 > ref_end_1 - homopolymer_end_1 )
  {
    
    int32_t on_ref_pos_1 = a.reference_start_1();
    int32_t on_read_pos_1 = 1;
    
    // First we parse through to the part of the read at the left boundary of the repeat
    for( size_t i=0; i<cigar_pair.size(); i++ ) {
      char op = cigar_pair[i].first;
      uint16_t len = cigar_pair[i].second;
      
      if( op == 'S' )
        on_read_pos_1 += len;
    
      else if( op == 'M' ) {
        on_read_pos_1 += len;
        on_ref_pos_1 += len;  
      }
      
      // insertion wrt ref
      else if( op == 'I' ) 
        on_read_pos_1 += len;
      
      // deletion wrt ref
      else if( op == 'D' )
        on_ref_pos_1 += len;        
      else
        ERROR("Unknown CIGAR operation."); 
      
      if (on_ref_pos_1 >= homopolymer_start_1) 
      {
        // back up to right position
        on_read_pos_1 -= on_ref_pos_1 - homopolymer_start_1;
        on_ref_pos_1 = homopolymer_start_1;        
        break;
      }
    }
    
    if (verbose) cout << a.read_name() << endl;
    
    // Now find the length of the homopolymer repeat
    int32_t on_read_check_i = on_read_pos_1;
    int32_t on_ref_check_i = on_ref_pos_1;
    int32_t homopolymer_length = 0;
    
    if (verbose) cout << on_read_check_i << " " << on_ref_check_i << endl;
    if (verbose) cout << read_sequence[on_read_check_i-1] << endl;
    
    while( (on_read_check_i < static_cast<int32_t>(read_sequence.size())) 
          && (read_sequence[on_read_check_i-1] == current_region.base) ) {
      homopolymer_length++;
      on_read_check_i++;
    }
 
    while( (on_ref_check_i <= ref_length)
          && (ref_sequence[on_ref_check_i-1] == current_region.base) ) {
      on_ref_check_i++;
    }

    int32_t matched_right = 0;
    while( (on_ref_check_i <= ref_length) && (on_read_check_i <= read_length) 
          && (read_sequence[on_read_check_i-1] == ref_sequence[on_ref_check_i-1]) ) {
      matched_right++;
      on_read_check_i++;
      on_ref_check_i++;
    }
    
    on_read_check_i = on_read_pos_1-1;
    on_ref_check_i = on_ref_pos_1-1;
    
    int32_t matched_left = 0;
    while( (on_ref_check_i > 0) && (on_read_check_i > 0) 
          && (read_sequence[on_read_check_i-1] == ref_sequence[on_ref_check_i-1]) ) {
      matched_left++;
      on_read_check_i--;
      on_ref_check_i--;
    }
    
    if (verbose) cout << matched_left << "|" << homopolymer_length << "|" << matched_right << endl;
    
    if ((matched_left >= 5) && (matched_right >= 5)) {
      if (verbose) cout << "passed" << endl;
      repeats.back().freqs.push_back(homopolymer_length);
    }
  } else {
    
    int32_t on_ref_pos_1 = a.reference_end_1();
    int32_t on_read_pos_1 = a.read_length();
    
    // First we parse through to the part of the read at the left boundary of the repeat
    for(int32_t i=cigar_pair.size()-1; i>=0; i-- ) {
      char op = cigar_pair[i].first;
      uint16_t len = cigar_pair[i].second;
      
      if( op == 'S' )
        on_read_pos_1 -= len;
      
      else if( op == 'M' ) {
        on_read_pos_1 -= len;
        on_ref_pos_1 -= len;  
      }
      
      // insertion wrt ref
      else if( op == 'I' ) 
        on_read_pos_1 -= len;
      
      // deletion wrt ref
      else if( op == 'D' )
        on_ref_pos_1 -= len;        
      else
        ERROR("Unknown CIGAR operation."); 
      
      if (on_ref_pos_1 <= homopolymer_end_1) 
      {
        // back up to right position
        on_read_pos_1 += homopolymer_end_1 - on_ref_pos_1;
        on_ref_pos_1 = homopolymer_end_1;        
        break;
      }
    }
    
    if (verbose) cout << a.read_name() << endl;
    
    // Now find the length of the homopolymer repeat
    int32_t on_read_check_i = on_read_pos_1;
    int32_t on_ref_check_i = on_ref_pos_1;
    int32_t homopolymer_length = 0;
    
    if (verbose) cout << on_read_check_i << " " << on_ref_check_i << endl;
    if (verbose) cout << read_sequence[on_read_check_i-1] << endl;
    
    while( (on_read_check_i > 0) 
          && (read_sequence[on_read_check_i-1] == current_region.base) ) {
      homopolymer_length++;
      on_read_check_i--;
    }
    
    while( (on_ref_check_i > 0)
          && (ref_sequence[on_ref_check_i-1] == current_region.base) ) {
      on_ref_check_i--;
    }
    
    uint32_t matched_left = 0;
    while( (on_ref_check_i > 0) && (on_read_check_i > 0) 
          && (read_sequence[on_read_check_i-1] == ref_sequence[on_ref_check_i-1]) ) {
      matched_left++;
      on_read_check_i--;
      on_ref_check_i--;
    }
    
    on_read_check_i = on_read_pos_1+1;
    on_ref_check_i = on_ref_pos_1+1;
    
    uint32_t matched_right = 0;
    while( (on_ref_check_i <= ref_length) && (on_read_check_i <= read_length) 
          && (read_sequence[on_read_check_i-1] == ref_sequence[on_ref_check_i-1]) ) {
      matched_right++;
      on_read_check_i++;
      on_ref_check_i++;
    }
    
    if (verbose) cout << "rev|" << matched_left << "|" << homopolymer_length << "|" << matched_right << endl;
    
    if ((matched_left >= 5) && (matched_right >= 5)) {
      if (verbose) cout << "passed" << endl;
      repeats.back().freqs.push_back(homopolymer_length);
    }
  }
    
}

void contingency_loci_pileup::printStats(const string& output, cReferenceSequences& ref_seq_info)
{
  ofstream out(output.c_str());
  ASSERT(!out.fail(), "Could not open output file: " + output); 

  //Finds the maximum encountered
  size_t maxsize = 0;
  for( size_t i=0; i<repeats.size(); i++ ){
    if( maxsize < repeats[i].freqs.size() ){
      maxsize = repeats[i].freqs.size();
    }
  }
  
  //
  // Create header
  //
  
  vector<string> header_list;
  for( size_t j=1; j<maxsize; j++ ) {
    header_list.push_back(to_string(j) + "-bp");
  }
  
  // if we are in locus mode
  if( indices.size() )
    header_list.push_back("locus");
    
  header_list.push_back("seq_id");
  header_list.push_back("start");
  header_list.push_back("end");
  header_list.push_back("repeat_base");
  header_list.push_back("gene_position");
  header_list.push_back("gene");
  header_list.push_back("gene_product");
  
  out << join(header_list, "\t") << endl;
     
  //
  // Line for each repeat
  //
  
  for( size_t i=0; i<repeats.size(); i++ ) {
    
    vector<string> line_list;
    
    vector<string> description_list = split_on_any(repeats[i].region, ":-");
    
    cout << i << " " << repeats[i].region << endl;
    
    // all of the base count columns
    for( size_t j=1; j<maxsize; j++ )
      line_list.push_back( (j<repeats[i].freqs.size()) ? to_string<int32_t>(static_cast<int32_t>(repeats[i].freqs[j])) : "0");
    
    // Checks if it is a contingency loci. If so, prints out the name of the locus
    for( size_t j=1; j<=indices.size(); j++ ) {
      if( atoi( description_list[1].c_str() ) == indices[j] ) {
        line_list.push_back(names[j]);
        break;
      }
    }
    
    line_list.push_back(description_list[0]);
    line_list.push_back(description_list[1]);
    line_list.push_back(description_list[2]);
    line_list.push_back(description_list[3]);

    
    cDiffEntry de(SUB);
    de[SEQ_ID] = description_list[0];
    de[POSITION] = description_list[1];
    de["size"] = to_string(from_string<int32_t>(description_list[2]) - from_string<int32_t>(description_list[1]) + 1);
    de["new_seq"] = "A"; // Dummy value - not used
    ref_seq_info.annotate_1_mutation(de, from_string<int32_t>(de[POSITION]), from_string<int32_t>(de[POSITION]) + from_string<int32_t>(de["size"]) - 1);
    
    line_list.push_back("\"" + de["gene_strand"] + " " + de["gene_position"] + "\"");
    line_list.push_back("\"" + de["gene_name"] + "\"");
    line_list.push_back("\"" + de["gene_product"] + "\"");

    out << join(line_list, "\t") << endl;
  }
  out.close();
}

void contingency_loci_pileup::readIndices( vector<int>& indices, vector<string>&  names, const string& loci)
{
  
  if (loci == "") return;
  
  if( loci.compare("NA") ){
    indices = vector<int>();
  }
  
  ifstream file( loci.c_str() );
  vector<int> ind;
  vector<string> n;
  int i=0;
  int x = 0;
  string s;
  while( !file.eof() ){
    if( i%5 == 1 ){
      file >> x;
      ind.push_back(x);
    }
    else if( i%5 == 4 ){
      file >> s;
      n.push_back(s);
    }
    else{
      file >> s; 
    }
    i++;
  }
  
  file.close();
  names = n;
  indices =  ind;
}

  
  
  
  
} // namespace breseq

