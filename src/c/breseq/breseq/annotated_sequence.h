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

#ifndef _BRESEQ_ANNOTATED_SEQUENCE_H_
#define _BRESEQ_ANNOTATED_SEQUENCE_H_

#include "breseq/common.h"

#include "breseq/fasta.h"
#include "breseq/alignment.h"
#include "breseq/genome_diff.h"

using namespace std;

namespace breseq {
	
	/*! Interface for loading sequences and sequence features from GenBank files.
  */
  
  
	/*! Sequence Feature class.

	 */ 

  // Currently everything is stored as strings...
  typedef string sequence_feature_key_t;
  typedef string sequence_feature_value_t;
  typedef map<sequence_feature_key_t, sequence_feature_value_t> sequence_feature_map_t; //!< Diff entry key-value map.
  
  class cSequenceFeature : public sequence_feature_map_t {
    
    public:
    
      // Could add accessors that convert strings to numbers...
      uint32_t m_start, m_end;
      int8_t m_strand;
    
      cSequenceFeature() {}
      cSequenceFeature(cSequenceFeature* _in) : sequence_feature_map_t(*_in) {
        m_start = _in->m_start;
        m_end = _in->m_end;
        m_strand = _in->m_strand;
      }
	  cSequenceFeature operator=(cSequenceFeature* _in) {
        m_start = _in->m_start;
        m_end = _in->m_end;
        m_strand = _in->m_strand;
        sequence_feature_map_t::operator=(*_in);
        return this;
      }
      
      //<! Safe accessor that returns empty string if not defined. 
      std::string SafeGet(sequence_feature_key_t in_key) { 
        sequence_feature_map_t::const_iterator it = this->find(in_key);
        if (it == this->end()) return std::string("");
        return it->second;
      }
      
      void ReadCoords(string& s, ifstream& in);
      void ReadTag(string& tag, string& s, ifstream& in);
  };

  

	/*! Sequence class.
	 */   
   
  class cAnnotatedSequence {
    
    public:      
      uint32_t m_length;
      string m_definition, m_version, m_seq_id;
    
      cFastaSequence m_fasta_sequence;            //!< Nucleotide sequence
      vector<cSequenceFeature> m_features;  //!< List of sequence features

    public:
    
      //Constructor for empty object
      cAnnotatedSequence() : 
        m_length(0), 
        m_definition("na"), 
        m_version("na"), 
        m_seq_id("na"),
        m_features(0) {} ;
    
      // Utility to get yop strand sequence
      string get_sequence(uint32_t start_1, uint32_t end_1) 
      {
        return m_fasta_sequence.m_sequence.substr(start_1 - 1, end_1 - start_1 + 1);
      }
  };

  
  /*! Reference Sequences
   
   Holds sequences and features for ALL reference sequences.
	 */ 
  
  class cReferenceSequences : public vector<cAnnotatedSequence> {
  protected:
    map<string,int> m_seq_id_to_index; // for looking up sequences by seq_id
    
  public:
    
    cReferenceSequences() {};    
    
    //!< Write a tab delimited feature 
    void WriteFeatureTable(const string &file_name);
    
    //!< Read a tab delimited feature
    void ReadFeatureTable(const string &file_name);
    
    //!< Write FASTA file       
    void WriteFASTA(const string &file_name);
    
    //!< Read FASTA file       
    void ReadFASTA(const std::string &file_name);
    
    //!< Write a tab delimited GFF3 file
    void WriteGFF( const string &file_name );  
      
    //!< Convert 
    uint32_t seq_id_to_index(const string& seq_id) 
      { assert(m_seq_id_to_index.count(seq_id)); return m_seq_id_to_index[seq_id]; };

    //!< Utility to get sequences by seq_id
    string get_sequence(const string& seq_id, uint32_t start_1, uint32_t end_1) 
    {
      return (*this)[seq_id_to_index(seq_id)].get_sequence(start_1, end_1);
    }

	void annotate_1_mutation(diff_entry& mut, uint32_t start, uint32_t end, bool repeat_override)
	{/*
		// this could be moved to the object
		my $intergenic_seperator = "/";

		// initialize everything, even though we don't always use it
		$mut->{aa_position} = "";
		$mut->{aa_ref_seq} = "";
		$mut->{aa_new_seq} = "";
		$mut->{codon_position} = "";
		$mut->{codon_ref_seq} = "";
		$mut->{codon_new_seq} = "";
		$mut->{gene_name} = "";
		$mut->{gene_position} = "";
		$mut->{gene_product} = "";
		@{$mut->{gene_list}} = (); #affected genes

		my $seq_id = $mut->{seq_id};

		my $gene_list_ref = $ref_seq_info->{gene_lists}->{$seq_id};
		my $repeat_list_ref = $ref_seq_info->{repeat_lists}->{$seq_id};
		my $ref_string = $ref_seq_info->{ref_strings}->{$seq_id};

		die "Unknown seq_id in reference sequence info: $seq_id\n" if ((!defined $gene_list_ref) || (!defined $repeat_list_ref) || (!defined $ref_string));

		my $size = $end - $start + 1;

		my ($prev_gene, $next_gene) = (undef, undef);
		my @within_genes = ();
		my @between_genes = ();
		my @inside_left_genes = ();
		my @inside_right_genes = ();

		my $repeat_region;
		if ($repeat_override)
		{
			die if ($start != $end);
			$repeat_region = get_overlapping_feature($repeat_list_ref, $start);
			push @within_genes, $repeat_region if ($repeat_region);
		}

		if (!$repeat_region)
		{
			my ($within_genes_list_ref, $between_genes_list_ref, $inside_left_genes_list_ref, $inside_right_genes_list_ref);
			($prev_gene, $next_gene, $within_genes_list_ref, $between_genes_list_ref, $inside_left_genes_list_ref, $inside_right_genes_list_ref)
				= find_nearby_genes($gene_list_ref, $start, $end);

			@within_genes = @$within_genes_list_ref;
			@between_genes = @$between_genes_list_ref;
			@inside_left_genes = @$inside_left_genes_list_ref;
			@inside_right_genes = @$inside_right_genes_list_ref;
		}


		## Mutation is intergenic
		if (scalar(@within_genes) + scalar(@between_genes) + scalar(@inside_left_genes) + scalar (@inside_right_genes)== 0)
		{
			$mut->{snp_type} = "intergenic";

			$mut->{gene_name} .= (defined $prev_gene) ? $prev_gene->{name} : "–";
			$mut->{gene_name} .= $intergenic_seperator;
			$mut->{gene_name} .= (defined $next_gene) ? $next_gene->{name} : "–";

			if (defined $prev_gene)
			{
				$mut->{gene_position} .= "intergenic (";
				$mut->{gene_position} .= ($prev_gene->{strand} == +1) ? "+" : "-";
				$mut->{gene_position} .= $start - $prev_gene->{end};
			}
			else
			{
				$mut->{gene_position} .= "intergenic (–";
			}
			$mut->{gene_position} .= $intergenic_seperator;
			if (defined $next_gene)
			{
				$mut->{gene_position} .= ($next_gene->{strand} == +1) ? "-" : "+";
				$mut->{gene_position} .= $next_gene->{start} - $end;
			}
			else
			{
				$mut->{gene_position} .= "–";
			}
			$mut->{gene_position} .= ")";

			$mut->{gene_product} .= (defined $prev_gene) ? $prev_gene->{product} : "–";
			$mut->{gene_product} .= $intergenic_seperator;
			$mut->{gene_product} .= (defined $next_gene) ? $next_gene->{product} : "–";

			return $mut;
		}
		// Mutation is completely within genes
		elsif (scalar @within_genes > 0)
		{
			/// TODO: It can be within multiple genes, in which case we need to annotate
			/// the change it causes in each reading frame UGH! YUCKY!
			/// FOR NOW: just take the first of the within genes...
			my $gene = $within_genes[0];
			$mut->{gene_name} = $gene->{name};
			$mut->{gene_product} = $gene->{product};

			#added for gene table
			@{$mut->{gene_list}} = ($gene->{name});

			my $within_gene_start = ($gene->{strand} == +1) ? $gene->{start} : $gene->{end};

			if ($start == $end)
			{
				$mut->{gene_position} = abs($start-$within_gene_start) + 1;
			}
			else
			{
				my $gene_start = abs($start-$within_gene_start) + 1;
				my $gene_end = abs($end-$within_gene_start) + 1;
				$mut->{gene_position} = ($gene_start < $gene_end) ? "$gene_start–$gene_end" : "$gene_end–$gene_start";
			}

			my $gene_nt_size = $gene->{end} - $gene->{start} + 1;

			## ...but the gene is a pseudogene or not a protein coding gene
			if ($gene->{pseudogene})
			{
				$mut->{snp_type} = "pseudogene";
				$mut->{gene_position} = "pseudogene ($mut->{gene_position}/$gene_nt_size nt)";
				return $mut;
			}
			elsif ($gene->{type} ne "protein")
			{
				$mut->{snp_type} = "noncoding";
				$mut->{gene_position} = "noncoding ($mut->{gene_position}/$gene_nt_size nt)";
				return $mut;
			}

			#only add gene information to SNPs and RA mutations that don't include indels...
			if (($mut->{type} ne 'SNP') && !(($mut->{type} eq 'RA') && ($mut->{ref_base} ne '.') && ($mut->{new_base} ne '.')))
			{
				$mut->{gene_position} = "coding ($mut->{gene_position}/$gene_nt_size nt)";
				return $mut;
			}

			## this is for RA...
			$mut->{ref_seq} = $mut->{ref_base} if (!defined $mut->{ref_seq});
			$mut->{new_seq} = $mut->{new_base} if (!defined $mut->{new_seq});

			## determine the old and new translation of this codon
			$mut->{aa_position} = int(($mut->{gene_position}-1)/3) + 1; ## 1 indexed
			$mut->{codon_position} = abs($start-$within_gene_start) % 3 + 1; ## 1 indexed

			my $codon_seq = ($gene->{strand} == +1) ?

				substr($ref_string, $gene->{start} + 3 * ($mut->{aa_position}-1) - 1, 3) :
				Breseq::Shared::revcom(substr($ref_string, $gene->{end} - 3 * $mut->{aa_position}, 3));

				#$ref_seq->trunc($gene->{start} + 3 * ($mut->{aa_position}-1),$gene->{start} + 3 * $mut->{aa_position} - 1) :
				#$ref_seq->trunc($gene->{end} - 3 * $mut->{aa_position}+1,$gene->{end} - 3 * ($mut->{aa_position}-1))->revcom;

			##Debug
			##print "$mut->{aa_position} $mut->{codon_position} $gene->{start} $gene->{end} $codon_seq\n";

			$mut->{codon_ref_seq} = $codon_seq;
			$mut->{aa_ref_seq} = bridge_translate($mut->{codon_ref_seq});
			#$mut->{aa_ref_seq} = $codon_seq->translate( undef, undef, undef, 11 )->seq();

			$mut->{codon_new_seq} = $codon_seq;
			#remember to revcom the change if gene is on opposite strand
			substr($mut->{codon_new_seq}, $mut->{codon_position} - 1, 1) = ($gene->{strand} == +1) ? $mut->{new_seq} : Breseq::Shared::revcom($mut->{new_seq});
			$mut->{aa_new_seq} =  bridge_translate($mut->{codon_new_seq});
			#$codon_seq->seq($mut->{codon_new_seq});
			#$mut->{aa_new_seq} = $codon_seq->translate( undef, undef, undef, 11 )->seq();

			$mut->{snp_type} = ($mut->{aa_ref_seq} ne $mut->{aa_new_seq}) ? "nonsynonymous" : "synonymous";
		}

		##The mutation actually contains several genes
		elsif (scalar(@between_genes) + scalar(@inside_left_genes) + scalar (@inside_right_genes) > 0)
		{
			my @gene_list = ( map({ "<i>[" . $_->{name} . "]</i>" } @inside_left_genes),
							  map({ "<i>" . $_->{name} . "</i>" } @between_genes),
							  map({ "<i>[" . $_->{name} ."]</i>" } @inside_right_genes) );


			#added for gene table
			@{$mut->{gene_list}} = ( map({ $_->{name} } @inside_left_genes),
									 map({ $_->{name} } @between_genes),
									 map({ $_->{name} } @inside_right_genes) );

			$mut->{gene_product} = join (", ", @gene_list);

			if (scalar @gene_list == 1)
			{
				$mut->{gene_name} = $gene_list[0];
			}
			else
			{
				$mut->{gene_name} = $gene_list[0] . "–" . $gene_list[-1];
			}
		}

		return $mut;*/
	}

	void annotate_mutations(genome_diff& gd, bool only_muts)
	{
		//keep track of other mutations that affect SNPs
		//because we may double-hit a codon

		//TODO: the proper way to do this is to create list of SNPs that have been hit
		// hashed by gene protein accession ID and AA position within gene
		// and have the annotation point to them (and back at them)
		// so that the codon will be correctly updated with all changes and we can notify the
		// changes that their SNP_type is not really SNP, but multiple hit SNP.

		/*my $snp_hits_hash;

		MUT: foreach my $mut ($gd->list)
		{
			next MUT if ($only_muts && (length($mut->{type}) != 3));

			if ($mut->{type} eq 'SNP')
			{
				$mut->{_ref_seq} = get_sequence($ref_seq_info, $mut->{seq_id}, $mut->{position}, $mut->{position});
				annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position});
			}
			elsif ($mut->{type} eq 'SUB')
			{
				annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position} + $mut->{size} - 1);
			}
			elsif ($mut->{type} eq 'DEL')
			{
				annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position} + $mut->{size} - 1);
			}
			elsif ($mut->{type} eq 'INS')
			{
				annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position});
			}
			elsif ($mut->{type} eq 'CON')
			{
				annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position} + $mut->{size} - 1);
			}
			elsif ($mut->{type} eq 'MOB')
			{
				annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position} + $mut->{duplication_size}-1);
			}
			elsif ($mut->{type} eq 'INV')
			{
				annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position});
				($mut->{gene_name_1}, $mut->{gene_product_1})  = ($mut->{gene_name}, $mut->{gene_product});
				annotate_1_mutation($ref_seq_info, $mut, $mut->{position} + $mut->{size}-1, $mut->{position} + $mut->{size}-1);
				($mut->{gene_name_2}, $mut->{gene_product_2})  = ($mut->{gene_name}, $mut->{gene_product});
				delete $mut->{gene_name};
				delete $mut->{gene_product};
			}
			elsif ($mut->{type} eq 'AMP')
			{
				annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position} + $mut->{size} - 1);
			}
			elsif ($mut->{type} eq 'JC')
			{
				annotate_1_mutation($ref_seq_info, $mut->{_side_1}, $mut->{side_1_position}, $mut->{side_1_position}, 1);
				annotate_1_mutation($ref_seq_info, $mut->{_side_2}, $mut->{side_2_position}, $mut->{side_2_position}, 1);
			}
			elsif ($mut->{type} eq 'RA')
			{
				annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position});
			}
			elsif ($mut->{type} eq 'MC')
			{
				annotate_1_mutation($ref_seq_info, $mut, $mut->{start}, $mut->{end});
			}
		}*/
	}
    
    map<string,int32_t> seq_order;
    map<string,string> trims;
    map<string,string> ref_strings;

    map<string, vector<cSequenceFeature> > repeat_lists;
    static cSequenceFeature* find_closest_repeat_region(uint32_t position, vector<cSequenceFeature>& repeat_list_ref, uint32_t max_distance, int32_t direction);
    vector<string> seq_ids; //< @GRC need to implement, get_keys(m_seq_id_to_index)?
  };
  
  /*! Helper function for creating cReferenceSequences
   */
  
  void LoadGenBankFile(cReferenceSequences& s, const vector<string>& in_file_names);
  bool LoadGenBankFileHeader(ifstream& in, cReferenceSequences& s);
  void LoadGenBankFileSequenceFeatures(ifstream& in, cAnnotatedSequence& s);
  void LoadGenBankFileSequence(ifstream& in, cAnnotatedSequence& s);
  
  void LoadFeatureIndexedFastaFile(cReferenceSequences& s, const string &in_feature_file_name, const string &in_fasta_file_name);
  
  void LoadBullFile(cReferenceSequences& s, const vector<string>& in_file_names);
  void LoadBullFeatureFile(ifstream& in, cAnnotatedSequence& s);

  /*! Utility functions.
  */
    
  std::string GetWord(string &s);
  void RemoveLeadingWhitespace(string &s);
  void RemoveLeadingTrailingWhitespace(string &s);

  uint32_t alignment_mismatches(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info);
  string shifted_cigar_string(const alignment_wrapper& a, const cReferenceSequences& ref_seq_info);

} // breseq namespace

#endif
