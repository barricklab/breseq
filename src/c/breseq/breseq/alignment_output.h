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

#ifndef _BRESEQ_ALIGNMENT_OUTPUT_H_
#define _BRESEQ_ALIGNMENT_OUTPUT_H_

#include "common.h"
#include "pileup_base.h"
#include "pileup.h"
#include "alignment.h"


using namespace std;

namespace breseq {

/*! This class is a FACTORY for generating HTML alignments
 */

//! returns more information about aligned reads given a sequence id string.
typedef struct {
    string seq_id;
    uint32_t length;
    string read_sequence;
    string qual_sequence;
    string aligned_bases;
    string aligned_quals;
    uint32_t reference_start;
    uint32_t reference_end;
    uint32_t start;
    uint32_t end;
    int32_t strand;
    bool    updated; //whether the read was updated at this pileup iteration already
}struct_aligned_read;


//!Helper struct for aligned_refs
typedef struct {
    uint32_t start;
    uint32_t end;
    string aligned_bases;
    string aligned_quals;
    uint32_t reference_length;
    string reference_name;
    char base;
    
}struct_aligned_reference;

typedef struct {
    string aligned_bases;
}struct_aligned_annotation;
    
class alignment_output_pileup : public pileup_base {
public:
    //! Constructor.
    alignment_output_pileup(const string& bam, const string& fasta,
                            const uint32_t maximum_to_align);
    //! Destructor.
    virtual ~alignment_output_pileup();
    //! Called for each genome position.
    virtual void pileup_callback(const pileup& aligned_reference);
    //! Called for each aligned read.
    virtual void fetch_callback(const alignment& a);

    //!Helper struct for aligned_reads


    map<string, struct_aligned_read> aligned_reads;



    vector<struct_aligned_reference> aligned_references;


    struct_aligned_annotation aligned_annotation;

    uint32_t unique_start;
    uint32_t unique_end;
    uint32_t total_reads;
    uint32_t processed_reads;
    uint32_t maximum_to_align;


    uint32_t insert_start;
    uint32_t insert_end;

    uint32_t last_pos;
    uint32_t max_indel;
    
    
    char base;
};

class alignment_output {
private:
    alignment_output_pileup m_alignment_output_pileup_object;

public:
    //! Constructor.
    alignment_output(string bam, string fasta, uint32_t in_maximum_to_align);
    //! Output an HTML alignment.
    string html_alignment(const string region);
    void create_alignment(const string bam, const string fasta, const string region);


};


}//end namespace breseq
#endif
