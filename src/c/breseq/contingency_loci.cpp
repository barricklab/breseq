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
  
  homopolymer_repeat_list hr;
  identify_homopolymer_repeats(hr, ref_seqs);
  
  contingency_loci_pileup clp(bam,fasta);
  
  //for each homopolymer repeat {
    string region= "NC_012660:456-466";
    clp.analyze_contingency_locus(region);
  //}
  
}

void identify_homopolymer_repeats(homopolymer_repeat_list& hr, const cReferenceSequences& ref_seqs)
{

}
  

/*! Constructor.
 */
contingency_loci_pileup::contingency_loci_pileup(
															const string& bam,
															const string& fasta
                                                )
: pileup_base(bam, fasta)
{
  set_print_progress(true);
}

void contingency_loci_pileup::analyze_contingency_locus(const string& region) {
  do_fetch(region);
}

/*! Called for each alignment.
 */
void contingency_loci_pileup::fetch_callback(const alignment_wrapper& a) {
}
  
  




} // namespace breseq

