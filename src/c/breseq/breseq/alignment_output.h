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


class alignment_output {
public:
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
        bool updated; //whether the read was updated at this pileup iteration already
    }Aligned_Read;

    //! returns more information about an aligned references
    typedef struct {
        uint32_t start;
        uint32_t end;
        string aligned_bases;
        string aligned_quals;
        uint32_t reference_length;
        string reference_name;
        char base;
    }Aligned_Reference;

    typedef struct {
        string aligned_bases;
    }Aligned_Annotation;

    typedef map<string, Aligned_Read> Aligned_Reads;
    typedef vector<Aligned_Reference> Aligned_References;

private:
    //! Builds Aligned_Reads, Aligned_References and Aligned_Annotation
    class Alignment_Output_Pileup : public pileup_base {
    public:
        //! Constructor.
        Alignment_Output_Pileup(const string& bam, const string& fasta,
                                const uint32_t maximum_to_align);
        //! Destructor.
        virtual ~Alignment_Output_Pileup();
        //! Called for each genome position.
        virtual void pileup_callback(const pileup& aligned_reference);
        //! Called for each aligned read.
        virtual void fetch_callback(const alignment& a);

        Aligned_Reads aligned_reads;
        Aligned_References aligned_references;
        Aligned_Annotation aligned_annotation;
        uint32_t unique_start; //used in create alignment and passed to fetch
        uint32_t unique_end; //used in create alignment and passed to fetch
        uint32_t total_reads;
        uint32_t processed_reads;
        uint32_t maximum_to_align;

        uint32_t insert_start;
        uint32_t insert_end;

        uint32_t last_pos;
        uint32_t max_indel;
        char base;
    };
    //!Helper struct for set_quality_range
    typedef struct {
        vector<uint32_t> qual_to_color_index;
        vector<uint32_t> qaul_cutoffs;
    }Quality_Range;

    Alignment_Output_Pileup m_alignment_output_pileup;
    Aligned_Reads m_aligned_reads;
    Aligned_References m_aligned_references;
    Aligned_Annotation m_aligned_annotation;
    Quality_Range m_quality_range;

public:
    //! Constructor.
    alignment_output(string bam, string fasta, uint32_t in_maximum_to_align);
    //! Output an HTML alignment.
    string html_alignment(const string region);
    void create_alignment(const string bam, const string fasta, const string region);
    void set_quality_range();
private:
    string create_header_string();
    static bool sort_by_aligned_bases(const pair<string, Aligned_Read> a, const pair<string,Aligned_Read> b)
    {
        return (a.second.aligned_bases > b.second.aligned_bases);
    }
};






}//end namespace breseq
#endif
