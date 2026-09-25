/*****************************************************************************

 AUTHORS

   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors

 LICENSE AND COPYRIGHT

   Copyright (c) 2008-2010 Michigan State University
   Copyright (c) 2011-2025 The University of Texas at Austin
   Copyright (c) 2025-     Michigan State University

   breseq is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2, or (at your option) any later version.

   SPDX-License-Identifier: GPL-2.0-or-later

*****************************************************************************/

#ifndef _BRESEQ_MUTATION_PREDICTOR_H_
#define _BRESEQ_MUTATION_PREDICTOR_H_

#include "common.h"

#include "alignment.h"
#include "reference_sequence.h"
#include "genome_diff.h"
#include "settings.h"
#include "fastq.h"
#include "candidate_junctions.h"

using namespace std;

namespace breseq {

  //! FALLBACK largest run of coverage that a CN copy-number-0 region is allowed to bridge when
  //! joining missing-coverage fragments (see merge_MC_fragments_spanned_by_CN).
  //!
  //! The bound is really one CNery tile -- the HMM decides copy number over a whole window, so an
  //! island it called deleted through is narrower than one window by construction -- and
  //! merge_MC_fragments_spanned_by_CN reads that width per region from the CN entry's tile_size.
  //! This literal is only what it falls back to for a CN entry that carries no tile_size, i.e. one
  //! read back from a hand-made genome diff. It is deliberately generous: real islands measure
  //! 21-62 bp, and being wrong here costs a merge that should have happened, not a false one.
  const uint32_t kMaxMCMergeIslandBases = 200;

  //! Join missing-coverage fragments that a single CN copy-number-0 region spans into one MC.
  //!
  //! The MC state machine has no gap tolerance -- one column above the propagation cutoff ends the
  //! entry (identify_mutations.cpp) -- so a few reads mismapped into the middle of a real deletion
  //! split it into fragments. Each fragment then seeds predictMCtoDELbyHomology's register search
  //! with its own (too short) extent, and every fragment misses. A CN region called at copy number
  //! zero across the whole span is independent evidence that the island is not real coverage.
  //!
  //! Extents come from the MC fragments, never from the CN region: CN boundaries land on CNery's
  //! tiling while MC boundaries are per-base, and the deletion callers want the precise ones.
  void merge_MC_fragments_spanned_by_CN(cGenomeDiff& gd, uint32_t fallback_max_island = kMaxMCMergeIslandBases);

  //! Options for MutationPredictor::predict_gene_conversions. Its own struct rather than fields
  //! of Settings so that gdtools NORMALIZE, which has no Settings of its own, and the pipeline
  //! can both drive it.
  struct cGeneConversionOptions {
    //! A CON is only made when its identical stretch explains at least this many mutations.
    int32_t minimum_mutations;
    //! The converted recipient must match the donor exactly, without gaps, over a stretch at
    //! least this long that contains every replaced mutation.
    int32_t minimum_identical_length;
    //! Length of the exact k-mer, taken across each mutation, used to locate candidate donors.
    int32_t seed_length;
    bool    verbose;

    cGeneConversionOptions()
      : minimum_mutations(1)
      , minimum_identical_length(50)
      , seed_length(20)
      , verbose(false)
    {}
  };

	class MutationPredictor
	{
	public:

		static cReferenceSequences ref_seq_info;

		MutationPredictor(cReferenceSequences& ref_seq_info);


    void prepare_junctions(Settings& settings, Summary& summary, cGenomeDiff& gd);

    // Functions that handle specific predictions
    void predictMCplusJCtoDEL(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc);
    //! Deletions between two near-identical copies of an UNANNOTATED repeat, where no split-read
    //! junction can exist. Reads data/reference.bam directly, including the redundantly placed reads
    //! that base calling discards, and falls back to RA evidence when there is no BAM. A discordant
    //! pair is used when one is present, but only to seed the search and as extra evidence.
    void predictMCtoDELbyHomology(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& mc, diff_entry_list_t& dp, diff_entry_list_t& ra);
    void predictJCplusJCtoMOB(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc);
    void predictJCtoINSorSUBorDEL(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc, bool use_redundant_sides = false);
    void predictRAtoSNPorDELorINSorSUB(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& ra, diff_entry_list_t& mc );

    //! Tandem amplification: a junction pointing back upstream on the same strand, whose span a CN
    //! region calls at two or more copies. The junction gives base-resolution endpoints and proves
    //! the copies are adjacent; CN supplies the copy number, which a junction alone cannot.
    //!
    //! Deliberately separate from predictJCtoINSorSUBorDEL, which caps a junction-derived event at
    //! kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation while normalize_INS_to_AMP requires a repeat
    //! unit LARGER than that same constant -- so no junction can reach AMP through that path. This
    //! one takes only junctions above the cap, leaving that path's behavior untouched.
    void predictJCplusCNtoAMP(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& cn);

    //! Amplification between an existing copy of a repeat family and a new copy of it. The existing
    //! copy is an annotated repeat_region at one end of the CN region; the new copy exists only in
    //! this sample, so it is never annotated and can only be anchored by the MOB or the single
    //! junction that reports it. The junction at the end lying inside the existing copy has a
    //! redundant side and is discarded upstream, which is why coverage has to close the gap.
    //!
    //! This is recombination between two copies that both exist by the time it happens, so the tag is
    //! always between=, never mediated= -- mediated= describes a different event, where a new element
    //! is laid down with each new copy, and it would make APPLY insert one element too many.
    //!
    //! The duplicated block is written to END on a repeat copy, so repeating it reproduces the
    //! element between every pair of copies, exactly as the real amplification does. When the copy at
    //! the high-coordinate end is the NEW one, it is not in the reference to point at, so the AMP is
    //! placed within= the MOB that reports it and begins inside that element. Coordinates are taken
    //! at the highest positions that give the same result, the convention gdtools NORMALIZE enforces
    //! so two runs of the same event compare equal.
    void predictCNplusRepeatToAMP(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& cn);
    
    void filter_JC_indel_evidence(Settings& settings, Summary& summary, cGenomeDiff& gd, diff_entry_list_t& jc, diff_entry_list_t& mc);


    // Be safe against having mutations on the same reference with a full deletion of that reference
    void remove_mutations_on_deleted_reference_sequences(Settings& settings, Summary& summary, cGenomeDiff& gd);
    
    // Set ignore flag on evidence before we predict mutations
    void ignore_evidence_near_contig_ends(Settings& settings, Summary& summary, cGenomeDiff& gd);
    void remove_mutations_near_contig_ends(Settings& settings, Summary& summary, cGenomeDiff& gd);

    // Remove SC evidence whose clip window is already explained by a predicted mutation
    void remove_soft_clipping_near_mutations(Settings& settings, cGenomeDiff& gd);

    // Attach a discordant-pair (DP) evidence item as evidence for a mutation when its sides
    // (seq_id/position/strand) exactly match a JC that already supports that mutation.
    void add_matching_DP_evidence(cGenomeDiff& gd);

    // Attach a pair-distance (PD) evidence item as evidence for a mutation when its sides exactly
    // match a JC that already supports that mutation.
    void add_matching_PD_evidence(cGenomeDiff& gd);

    //! Name the mobile element a pair-distance (PD) insertion is consistent with, where the data
    //! supports naming one. Annotation only -- no mutation is created and no gate changes.
    void annotate_PD_with_repeat_family(cGenomeDiff& gd);

    // Attach a copy number (CN) evidence item as evidence for a DEL mutation that already explains
    // it. Never used to predict a deletion -- only to stop a copy-number-zero region describing a
    // deletion breseq already called from being reported a second time as unassigned evidence.
    void add_matching_CN_evidence(cGenomeDiff& gd);

    // Remove any DP item that a PD item already describes. PD sees the whole distribution where DP
    // sees only its discordant tail, so where both fire on one breakpoint the DP is redundant and
    // its coordinates are the worse of the two.
    void remove_DP_superseded_by_PD(Settings& settings, Summary& summary, cGenomeDiff& gd);

    // Combine a DP item with a MOB whose supporting JC unique side matches the DP unique side (within
    // a small window + same strand), even when their IS (redundant) sides are on different copies of
    // the element. Attaches the DP as evidence and, when the DP localizes a specific IS copy whose
    // sequence differs from the family consensus, sets mob_region= to override the MOB's default
    // consensus IS.
    void combine_DP_with_MOB_by_unique_side(cGenomeDiff& gd);

    // Cleans up indel positions by shifting them and adds addition fiels for simple sequence repeats 
    void normalize_and_annotate_tandem_repeat_mutations(Settings& settings, Summary& summary, cGenomeDiff& gd);
    
    // Cleans up predictions of large INS to make them AMP
    void normalize_INS_to_AMP(Settings& settings, Summary& summary, cGenomeDiff& gd);

    //! Replace groups of small mutations with the gene conversion (CON) that explains them.
    //!
    //! Recombination between two near-identical copies of a sequence -- paralogous genes, rRNA
    //! operons, copies of a mobile element -- overwrites a tract of the recipient copy with the
    //! donor copy's sequence. Read mapping sees that as a run of SNP/INS/DEL/SUB calls in the
    //! recipient, one at every position where the copies differed across the tract. Those calls
    //! reproduce the right genome, but they describe one event as many, and gdtools COUNT and any
    //! phylogeny built on them treat each as an independent mutation.
    //!
    //! Every consensus SNP/INS/DEL/SUB mutation is applied to its reference sequence once, giving
    //! the mutated genome with a map from each base back to its reference coordinate. A k-mer of
    //! that mutated sequence across each mutation is looked for everywhere in the reference, on
    //! both strands, and every hit is extended without gaps in both directions until the mutated
    //! sequence and the candidate donor disagree. A mutation is explained by a hit when its whole
    //! footprint lies inside that identical stretch: because the comparison is against the
    //! MUTATED sequence, a mutation can only sit inside it if the donor really carries the new
    //! base(s) and the reference did not. There is no notion of a cluster: a stretch absorbs
    //! every explained mutation it reaches, however far apart they are, which is what an rRNA
    //! operon needs. The best hit (most mutations explained, then the longest identical stretch)
    //! is accepted when it explains at least minimum_mutations and the identical stretch is at
    //! least minimum_identical_length long.
    //!
    //! The CON written has the MAXIMAL extent: its tract is the whole identical stretch, every
    //! base the converted recipient shares with the donor. Every tract inside that stretch that
    //! still covers the mutations applies to the same genome, so the widest one is the canonical
    //! form, and it is the only form that says which strand the donor is on when the conversion
    //! changed a single base. A stretch never covers a base that belongs to any other mutation
    //! (one that is polymorphic, pinned by within/before/mediated/between, of another type, or a
    //! CON already made), so tracts never overlap and a de novo SNP inside a converted region
    //! yields two CON entries around a SNP rather than a CON with a SNP inside it (VALIDATE would
    //! reject that: 'within' is not allowed on a CON).
    //!
    //! A hit whose donor overlaps the recipient stretch it explains is discarded. That is what an
    //! insertion or deletion of one unit of a tandem repeat looks like when the mutated sequence
    //! is matched against the same locus shifted by a unit, and those are already normalized as
    //! INS/DEL with repeat_* annotations. A donor that is the neighbouring, non-overlapping copy
    //! of a tandem array is kept. Copies of an equally good donor (a mobile element family) are
    //! broken deterministically: forward strand first, then reference order, then the lowest
    //! coordinate. Circular sequences are not searched across the origin.
    //!
    //! Returns the number of CON entries created. gd is left sorted. The donor is read from the
    //! reference as APPLY reads it, so the caller's APPLY-based self-check is a genuine test.
    int32_t predict_gene_conversions(cGenomeDiff& gd, const cGeneConversionOptions& options);

    //! The pipeline's gene conversion step. Runs predict_gene_conversions on the mutations the RA
    //! evidence produced, then reads the alignments to do what those mutations cannot: inside a
    //! converted tract breseq calls nothing, because the reads there either tie between the two
    //! copies and are discarded as redundant, or match the donor and are placed there. The SNPs
    //! only ever mark the edges. So every tract is extended along the reads that align to both
    //! copies (see extend_gene_conversion_by_reads), missing-coverage evidence left over from the
    //! deletion predictors is tested for the same signature, and the MC, DP, CN and RA evidence a
    //! conversion produces is attached to the CON so it stops being reported as unassigned.
    void predictGeneConversions(Settings& settings, Summary& summary, cGenomeDiff& gd,
                                diff_entry_list_t& mc, diff_entry_list_t& dp, diff_entry_list_t& cn);

    //! Attach the evidence a gene conversion leaves to its CON: missing coverage inside the tract or
    //! inside the donor (reads moved from one copy to the other), discordant pairs joining the flank
    //! just outside the tract to the corresponding end of the donor, copy number zero over the tract
    //! or a gain over the donor, and read alignment evidence inside the tract that reads the donor
    //! base. Claimed MC items are removed from the working list.
    void attach_gene_conversion_evidence(Settings& settings, Summary& summary, cGenomeDiff& gd,
                                         cDiffEntry& con, diff_entry_list_t& mc,
                                         diff_entry_list_t& dp, diff_entry_list_t& cn);

    // Master function
		void predict(Settings& settings, Summary& summary, cGenomeDiff& gd);

    // Helper functions
		static bool sort_by_hybrid(const counted_ptr<cDiffEntry>& a, const counted_ptr<cDiffEntry>& b);
		static bool sort_by_reject_score(const counted_ptr<cDiffEntry>& a, const counted_ptr<cDiffEntry>& b);
		static bool sort_by_pos(const counted_ptr<cDiffEntry>& a, const counted_ptr<cDiffEntry>& b);

	private:

    // Helper functions for tandem repeats
		cFeatureLocation* within_repeat(string seq_id, int32_t position);
    cFeatureLocation* near_repeat_boundary(string seq_id, int32_t position, int32_t direction);
    void find_repeat_unit(string& mutation_sequence, uint32_t& repeat_size, string& repeat_sequence);
    
    void normalizeINSposition(cAnnotatedSequence& ref_seq, cDiffEntry& de, string& repeat_sequence);
    void normalizeDELposition(cAnnotatedSequence& ref_seq, cDiffEntry& de, string& repeat_sequence);

    
    
    uint32_t find_original_num_repeat_units(cAnnotatedSequence& ref_seq, int32_t position, string& repeat_sequence);
    
    void combine_newly_adjacent_mutations(cGenomeDiff& gd);
    
    void assign_before_within_to_mutations(cGenomeDiff& gd);

    
	}; // class MutationPredictor
  
  
  // effects of base substitutions
  // List is 4 times as long as the genome, 
  // with a slot or all possible nucleotide changes (including no change) at a position
  // Each subset of entries is the effect of a change to ('A')('T')('C')('G'), respectively.
  enum BaseSubstitutionEffect {
    intergenic_base_substitution,
    noncoding_base_substitution,
    pseudogene_base_substitution,
    synonymous_coding_base_substitution,
    nonsynonymous_coding_base_substitution,
    nonsense_coding_base_substitution,
    unknown_coding_base_substitution,  // for dealing with degenerate bases
    no_change_base_substitution
  };
  typedef vector<BaseSubstitutionEffect> SequenceBaseSubstitutionEffects;
  
  inline BaseSubstitutionEffect string_to_bse(const string& s) {
    if (s == "intergenic") return intergenic_base_substitution;
    if (s == "noncoding") return noncoding_base_substitution;
    if (s == "pseudogene") return pseudogene_base_substitution;
    if (s == "synonymous") return synonymous_coding_base_substitution;
    if (s == "nonsynonymous") return nonsynonymous_coding_base_substitution;
    if (s == "nonsense") return nonsense_coding_base_substitution;
    if (s == "unknown") return unknown_coding_base_substitution;
    if (s == "no_change") return no_change_base_substitution;

    ERROR("Attempt to convert unrecognized string to base substititon effect: " + s);
    return unknown_coding_base_substitution;
  }

  inline string bse_to_string(BaseSubstitutionEffect bse) {
    switch (bse){
        case intergenic_base_substitution: return "intergenic";
        case noncoding_base_substitution: return "noncoding";
        case pseudogene_base_substitution: return "pseudogene";
        case synonymous_coding_base_substitution: return "synonymous";
        case nonsynonymous_coding_base_substitution: return "nonsynonymous";
        case nonsense_coding_base_substitution: return "nonsense";
        case unknown_coding_base_substitution: return "unknown";
        case no_change_base_substitution: return "no_change";
    }
    return "";
  }
  
  
  enum BaseType {
    intergenic_base,
    noncoding_base,
    pseudogene_base,
    protein_base
  };
  
  // 
  // List is as long as the genome
  enum BaseCDSStrand {
    no_CDS,
    forward,
    reverse,
    conflict
  };
  typedef vector<BaseCDSStrand> SequenceBaseCDSStrands;
  
  // Organizes information about a single position for easier counting
  struct BaseSubstitutionEffectPositionInfo {
    base_char                       m_base_char;
    BaseType                        m_base_type;
    vector<BaseSubstitutionEffect>  m_base_substitution_effect; // For change to ('A')('T')('C')('G')
    BaseCDSStrand                   m_base_cds_strand;          // For change to ('A')('T')('C')('G')
  };
  
  //<! Class for annotating effects of base substitution mutations
  //   that properly accounts for overlapping features.
  class BaseSubstitutionEffects {    
  public:
        
    static string                   separator;
    static vector<base_char>        base_char_list;
    static map<base_char, uint8_t>  base_char_to_base_index;
    static map<string,uint8_t>      base_change_to_base_pair_change_index;
    static vector<string>           base_change_list;
    static vector<string>           base_pair_change_list;
    static map<string,string>       base_change_to_base_pair_change;
    static vector<string>           base_change_type_list;
    static vector<string>           base_type_list;
    static map<BaseSubstitutionEffect,BaseType>  snp_type_to_base_type;
    
    map<string,SequenceBaseSubstitutionEffects> m_bse;
    map<string,SequenceBaseCDSStrands> m_bcs;
    
    void initialize_from_sequence(cReferenceSequences& ref_seq_info);
    
    BaseSubstitutionEffectPositionInfo position_info_1(
                                                       cReferenceSequences& ref_seq_info, 
                                                       string seq_id, 
                                                       uint32_t pos_1
                                                       );
    
  }; // class BaseSubstitutionEffects

  
  class BaseSubstitutionEffectCounts  {
  public:
    void initialize_possible_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse);
    void change_position_1_possible_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse, string seq_id, uint32_t pos_1, int32_t inc);
    void change_position_1_observed_totals(cReferenceSequences& ref_seq_info, BaseSubstitutionEffects& bse, string seq_id, uint32_t pos_1, string new_base, int32_t inc);

    static vector<string> base_pair_change_count_list; 
    static vector<string> base_change_type_count_list;

    struct BaseTypeCounts : public map<string,int32_t> {
      
      BaseTypeCounts()
      { 
        (*this)["nt"] = 0; (*this)["gc"] = 0; (*this)["at"] = 0;
      }
    };
    
    // map by base pair change
    struct BasePairChangeCounts : public map<string,int32_t> {
      
      BasePairChangeCounts()
      { 
        for (vector<string>::iterator it = BaseSubstitutionEffects::base_pair_change_list.begin();
             it != BaseSubstitutionEffects::base_pair_change_list.end(); ++it) {
          (*this)[*it] = 0;
        }
        (*this)["TOTAL"] = 0;
      }
    }; //class BaseSubstitutionEffectCounts
    
    map<BaseType, BaseTypeCounts > m_base_counts;
    map<string, BasePairChangeCounts > m_possible_base_pair_change_counts;
    map<string, BasePairChangeCounts > m_observed_base_pair_change_counts;

  };
  
  void MutationCountFile(
                         cReferenceSequences& ref_seq_info, 
                         vector<cGenomeDiff>& genome_diffs, 
                         string& output_file_name,
                         string& detailed_output_file_name,
                         bool base_substitution_statistics,
                         bool count_polymorphisms,
                         bool calculate_genome_size,
                         bool verbose
                         );
  
  
} // namespace breseq

#endif
