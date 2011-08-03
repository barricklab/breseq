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

#include "breseq/contingency_loci.h"

using namespace std;

namespace breseq {

    vector<int> readIndices(){
        ifstream file( "cloci.txt");
        vector<int> ind;
        int i=0;
        int x = 0;
        string s;
        
        while( !file.eof() ){
            if( i%4 == 1 ){
                file >> x;
                ind.push_back(x);
            }
            else{
                file >> s; 
            }
            i++;
        }
        
        file.close();
        cout << "DONE!";
        return ind;
    }
    
    
    
    
    
                            
/*! Convenience wrapper around the identify_mutations_pileup class.
 */
void analyze_contingency_loci(
								const string& bam,
								const string& fasta,
								const string& output
 ) {
       
	// do the mutation identification:
    cReferenceSequences ref_seqs;
    LoadFeatureIndexedFastaFile(ref_seqs, "", fasta);
  
    cout << "Finding repeats...\n";
    homopolymer_repeat_list hr;
    identify_homopolymer_repeats(hr, ref_seqs);
  
    contingency_loci_pileup clp(bam,fasta);
    

    
    cout << "Analyzing alignments...\n";
    // For each repeat
    for( int i=0; i<hr.size(); i++ ){
        //cout << ".";
        //cout.flush();
        // Makes the region string
        stringstream ss;
        ss << hr[i].start;  ss << "-"; ss << hr[i].start + hr[i].length; 
        string region = hr[i].seq_id.append(":").append( ss.str() ).append(":") + hr[i].base;
        // string region= "NC_012660:456-466:T";
        cout << i << " ";
        cout.flush();

        clp.analyze_contingency_locus(region);
    }
    
    
    // Output the histograms
    cout << "This should be printing. \n";
    clp.printStats( output );
      
}
    


void identify_homopolymer_repeats(homopolymer_repeat_list& hr, const cReferenceSequences& ref_seqs)
{
    //For each sequence
    for( int i=0; i<ref_seqs.size(); i++ ){
        
        char base = ' ';
        int rnumber = 0;
        
        // For each nucleotide in the sequence
        for( int j=0; j<ref_seqs[i].m_fasta_sequence.m_sequence.size(); j++ ){
            if( ref_seqs[i].m_fasta_sequence.m_sequence[j] == base ){
                rnumber++;
            }
            else{
                if( rnumber >= 8 ){
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
															const string& fasta
                                                )
: pileup_base(bam, fasta)

{    //fastaf = fasta;
    cout << "Filename: " << split( bam, "." )[0].append(".tam") << "\n";
    cout << "Reading indices...\n";
    indices = readIndices();    
    cout << "Done!\n";
    tf.open_write( split( bam, "." )[0].append(".tam"), fasta );
    set_print_progress(true);
}

void contingency_loci_pileup::analyze_contingency_locus(const string& region) {
    
    //cout << region << "\n";
    vector<string> r = split( region, ":");
    vector<string> coords = split( r[1], "-" );
    
    //cout << "BASE: " << r[2] << "\n";
    current_region.base = r[2][0];
    current_region.start = atoi( coords[0].c_str() );
    current_region.length = atoi( coords[1].c_str() ) - current_region.start;
    
    //////// If the region is a contigency locus
    
    /*for( int i=0; i<indices.size(); i++ ){
        if( indices[i] == current_region.start ){
            return;
        }
    }*/
    
    int boo = 0;
    for( int i=0; i<indices.size(); i++ ){
        if( indices[i] == current_region.start ){
            boo = 1;
        }
    }
    if( boo == 0 ){
        return;
    }
    ///// If the region is a contingency locus
    
    // Creates a repeats_stats object and pushes it into the class's array of
    repeats.push_back( repeat_stats( region ) );
    do_fetch( r[0].append(":").append(r[1]) );
    
    //printf( "Printing counts: \n" );
    
    // Constructs histogram from the values added from fetch_callback
    vector<double> histogram(1,0);
    int count = 0;
    //printf( "Freqs.size(): %i\n", repeats.back().freqs.size() );
    for( int i=0; i<repeats.back().freqs.size(); i++ ){
        //printf( "%f ", repeats.back().freqs[i] );
        if( histogram.size() < repeats.back().freqs[i]+1 ){
            histogram.resize( repeats.back().freqs[i]+1 );
            histogram[ repeats.back().freqs[i] ] = 1;
            count++;
        }
    
        else{
            //printf( "%f ", repeats.back().freqs[i] );
            histogram[ int(repeats.back().freqs[i]) ]++;
            count++;
        }
    }
    //Normalizes it
    /*
    for( int i=0; i<histogram.size(); i++ ){
        histogram[i] = histogram[i]/count;
    }*/
    /*
    cout << "Histogram: ";
    for( int i=0; i<repeats.back().freqs.size(); i++ ){
        cout << repeats.back().freqs[i] << " ";
    }
    cout << "\n";
    getchar();*/
    
    //getchar();
    repeats.back().freqs = histogram;
    /*
    cout << "Histogram:\n";
    for( int i=0; i<histogram.size(); i++ ){
        printf( "%f ", histogram[i] );
    }
    cout << "\n";*/
    //getchar();*/

}

/*! Called for each alignment.
<<<<<<< local
    Counts the lengths of the repeats in each of the appropriate alignments and stores it into the appropriate repeats_stats in the array
 */ 
void contingency_loci_pileup::fetch_callback(const alignment_wrapper& a) {
    
    vector<pair<char,uint16_t> > cigar_pair = a.cigar_pair_array();
    //cout << "Alignment: " << a.read_name() << "\n";
    int dist_to_first_match = 0;
    
    // Goes through cigar to see if the alignment is appropriate. Then counts the length of the repeat
    for( int i=0; i<cigar_pair.size(); i++ ){
        if( cigar_pair[i].first == 'S' ){
            dist_to_first_match += cigar_pair[i].second;
        }
        // If it begins the matching
        if( cigar_pair[i].first == 'M' ){
            
            // If the alignment extends 5 bases before the region
            if( current_region.start - a.reference_start_1() >= 5 ){
                
                // If the matched/deleted portion extends 5 bases after the end of the repeat, count the length of the homopolymer repeat in the alignment
                if( ( (int)a.reference_end_1() - (int)(current_region.start+current_region.length)) > 5 ){

                    string sequence = a.read_char_sequence();
                    char base = current_region.base;

                    
                    //Find where in the read to begin counting repeats from the first match
                    int ref_bypass_length = current_region.start - a.reference_start_1();
                    int read_bypass_length = 0;
                    int cigar_start_pos;

                    
                    for( int j=0; j+i<cigar_pair.size(); j++ ){
                        if( cigar_pair[i+j].first == 'M' ){
                            if( ref_bypass_length - cigar_pair[i+j].second >= 0 ){
                                read_bypass_length += cigar_pair[i+j].second;
                                ref_bypass_length -= cigar_pair[i+j].second;
                            }
                            else{
                                read_bypass_length += ref_bypass_length;
                                cigar_start_pos = j+i;
                                break;
                            }
                            
                        }
                        
                        if( cigar_pair[i+j].first == 'D' ){
                            if( ref_bypass_length - cigar_pair[i+j].second >= 0 ){
                                ref_bypass_length -= cigar_pair[i+j].second;
                            }
                            else{
                                cigar_start_pos = j+i;
                                break;
                            }

                        }
                        
                        if( cigar_pair[i+j].first == 'I' ){
                            read_bypass_length += cigar_pair[i+j].second;
                        }
                        
                        cigar_start_pos = j+i;
                        
                    }
                    int read_repeat_length = 0;
                    int ref_repeat_length = current_region.length;

                    for( int j=cigar_start_pos; j<cigar_pair.size(); j++ ){
                        
                        if( cigar_pair[j].first == 'M' && j == cigar_start_pos ){
                            if( ref_repeat_length - (cigar_pair[j].second - ref_bypass_length) >= 0 ){
                                ref_repeat_length -= (cigar_pair[j].second - ref_bypass_length);
                                read_repeat_length += (cigar_pair[j].second - ref_bypass_length);
                            }
                            else{
                                read_repeat_length += ref_repeat_length;

                                break;
                            }
                            
                        }
                        
                        else if( cigar_pair[j].first == 'M' ){
                            if( ref_repeat_length - cigar_pair[j].second >= 0 ){
                                ref_repeat_length -=  cigar_pair[j].second;
                                read_repeat_length += cigar_pair[j].second;
                            }
                            else{
                                read_repeat_length += ref_repeat_length;
                                break;
                            }
                            
                        }
                        
                        else if( cigar_pair[j].first == 'D' && j == cigar_start_pos ){
                            if( ref_repeat_length - (cigar_pair[j].second - ref_bypass_length) >= 0 ){
                                ref_repeat_length -= (cigar_pair[j].second - ref_bypass_length);
                                                            }
                            else{
                                break;
                            }
                            
                        }
                        
                        else if( cigar_pair[i+j].first == 'D' ){
                            if( ref_repeat_length - cigar_pair[i+j].second >= 0 ){
                                ref_repeat_length -= cigar_pair[i+j].second;
                            }
                            else{
                                break;
                            }
                            
                        }
                        
                        else if( cigar_pair[i+j].first == 'I' ){
                            read_repeat_length += cigar_pair[i+j].second;
                            
                            

                        }
                        
                    }
                    

                    //cout << "Read repeat length: " << read_repeat_length << "\n";
                    //getchar();
                    //Counts how long the repeat lasts
                    
                    int repeat_length = 0;
                    //for( int j=read_bypass_length+dist_to_first_match; j<read_bypass_length+dist_to_first_match+read_repeat_length; j++ ){
                    //cout << "begin match at: " << read_bypass_length+dist_to_first_match << "\n";
                    for( int j=0; j<read_repeat_length; j++ ){
                        if( sequence[j+read_bypass_length+dist_to_first_match] == current_region.base ){
                            repeat_length++;
                        }
                        else{
                            return;
                        }

                    }
                    
                    //getchar();

                    repeats.back().freqs.push_back(repeat_length);
                    alignment_list as;
                    cout << "Is something happening???\n";
                    
                    bam_alignment* b = new bam_alignment( a  );
                    counted_ptr<bam_alignment> bp(b);
                    as.push_back(bp);
                    
                    tf.write_alignments( 1, as, NULL );

                    return;
                }
            }
        }
    }
}

      




} // namespace breseq

