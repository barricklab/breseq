
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

#include "libbreseq/rna_seq.h"

using namespace std;

namespace breseq {

string good_gene_name(const string& gene, const string& product) {
  string ggn;
  ggn += gene;
  ggn += "|";
  ggn += product;
  ggn = substitute(ggn, " ", "_");
  return ggn;
}
  
void RNASeq::tam_to_gene_counts(cReferenceSequences& ref_seq_info, const string& fasta_file_name, const vector<string>& tam_files, const string& out_tab_file_name, bool verbose)
{
  int32_t num_datasets = tam_files.size();
  
  ofstream output_tab_file(out_tab_file_name.c_str());
  ASSERT(output_tab_file.good(), "Error opening output file: " + out_tab_file_name);

  //initialize list to zeros
  map<string,vector<int32_t> > gene_counts;
  
  for(vector<cAnnotatedSequence>::iterator it_as = ref_seq_info.begin(); it_as < ref_seq_info.end(); it_as++) {
    for (cSequenceFeatureList::iterator itr_feat=it_as->m_genes.begin(); itr_feat!=it_as->m_genes.end(); ++itr_feat) {
      Gene g(**itr_feat);
      string ggn = good_gene_name(g.name, g.product);
      gene_counts[ggn] = vector<int32_t>(num_datasets, 0);
    }
  }
  
  size_t i = 0;
  for(vector<string>::const_iterator tfi=tam_files.begin(); tfi != tam_files.end(); ++tfi) {
    cout << "SAM File: " << *tfi << endl;
    
    uint32_t num_unmapped = 0;
    uint32_t num_not_in_genes = 0;
    uint32_t num_in_genes = 0;
    uint32_t num_multiple_mapped = 0;
    
    tam_file tf(*tfi, fasta_file_name, ios_base::in);

    alignment_list al;
    
    cSequenceFeaturePtr best_gene(NULL);
    size_t pre_file_total = 0;

    while(tf.read_alignments(al)) {
      
      bam_alignment* a = al.front().get();

      
      pre_file_total++;
      if (pre_file_total % 100000 == 0) 
        cout << "  Aligned reads processed: " << pre_file_total << endl;
      
      if (verbose) {
        cout << " >" << a->read_name() << endl;
        cout << "  " << a->reference_target_id() << ":" << a->reference_start_1() << "-" << a->reference_end_1() 
          << "(" << a->strand() << ")" << endl;
      }
      
      
      // Ignore multiple alignments
      if (al.size() != 1) {
        if (verbose)
          cout << "  Multiple alignments (not counted)"<< endl;
        
        num_multiple_mapped++;
        continue;
      }
      
      if (a->unmapped()) {
        
        if (verbose)
          cout << "  Unmapped (not counted)"<< endl;

        num_unmapped++; 
        continue;
      }
      
      uint32_t tid = a->reference_target_id();
      
      cAnnotatedSequence& this_ref_seq = ref_seq_info[a->reference_target_id()];
      
      int32_t best_overlap = 0;
      for (cSequenceFeatureList::iterator itr_feat=this_ref_seq.m_genes.begin(); itr_feat!=this_ref_seq.m_genes.end(); ++itr_feat) {

        cSequenceFeature& feat = **itr_feat;
        
        if (feat.m_location.strand != a->strand())
          continue;
        
        int32_t ulength = overlap_length(feat.m_location.start, feat.m_location.end, a->reference_start_1(), a->reference_end_1());
        
        if (ulength > best_overlap) {
          best_overlap = ulength;
          best_gene = *itr_feat;
        }
      }
      
      if (best_gene.get()) {
        Gene g(*best_gene);
        
        if (verbose)
          cout << "  Best gene: " << g.name << " " << g.product << " " << best_gene->m_location.start << "-" << best_gene->m_location.end 
            << "(" << static_cast<int32_t>(best_gene->m_location.strand) << ")" << endl;
        

        // because we print tab-delimited, they can really screw up our columns
        string safe_gene_name = substitute(g.name, "\t", " ");
        string ggn = good_gene_name(g.name, g.product);
        gene_counts[ggn][i]++;
        
        //cerr << g.name << endl;
        num_in_genes++;
      }
      else {
        if (verbose)
          cout << "  Not in gene (not counted)"<< endl;
        
        num_not_in_genes++;
      }
    }
    
    cout << "  Unmapped: " << num_unmapped << endl;
    cout << "  Multiple mapped: " << num_multiple_mapped << endl;
    cout << "  Not in genes: " << num_not_in_genes << endl;
    cout << "  In genes: " << num_in_genes << endl;
    
    i++;
  }
  
  // Print out
  
  // One column for each data set
  
  output_tab_file << "gene";
  for(vector<string>::const_iterator it=tam_files.begin(); it!=tam_files.end(); it++) {
    output_tab_file << "\t" << *it;
  }
  output_tab_file << endl;
  
  // Each data row
  for(map<string,vector<int32_t> >::iterator mit=gene_counts.begin(); mit!=gene_counts.end(); mit++) {
    
    output_tab_file << mit->first;
    
    for(vector<int32_t>::iterator cit=mit->second.begin(); cit!=mit->second.end(); cit++) {
      output_tab_file << "\t" << *cit;
    }
    output_tab_file << endl;
  }
}

} // namespace breseq
