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

#ifndef _BRESEQ_GENOME_DIFF_H_
#define _BRESEQ_GENOME_DIFF_H_

#include "common.h"
#include "file_parse_errors.h"
#include "genome_diff_entry.h"

using namespace std;

namespace breseq {

class cReferenceSequences;
class Settings;
class cFeatureLocation;
class MutationTableOptions;

extern const int32_t kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation;
extern const int32_t kBreseq_large_mutation_size_cutoff;
  
/*! Genome Diff class.
 
 //  		Genome Diff files are tab delimitted. The first column defines the type of entry on a given line.
 //  		The second and third columns are type-nonspecific (id, parents), followed by type-specific
 //  		columns, then an arbitrary number of columns of more detailed data in a key=value format.
 
 */
class cGenomeDiff
{
public:
  
  enum group { MUTATIONS = 0, EVIDENCE, VALIDATION }; 

  typedef string key_t; 
  typedef vector<string> list_t;

  static uint32_t s_input_order_counter;
  
  //!---- Variables ---- !//
protected:  	
  string _file_path;                   //!< File name associated with this diff.
  diff_entry_list_t _entry_list;      //!< All diff entries.
  uint32_t _unique_id_counter;        //!< Smallest available id.
  map<string,bool> unique_id_used;

public:
  // @JEB should make this protected and add an accessor
  //! Metadata kept in .gd files
  struct Metadata
  {
    Metadata() : version("1.0"), time(-1.0), input_order(s_input_order_counter++) {}

    string version;
    string title;
    string author;
    string created;
    string program;
    string command;
    vector<string> ref_seqs;
    vector<string> read_seqs;
    vector<string> adapter_seqs;
    map<string,string> adapters_for_reads; // Keeps track of what adaptors belong to what read files
    vector<vector<string> > reads_by_pair; // Keeps track of pairs R1 and R2.
    map<string,string> breseq_data; // Use this to write values from pipeline to gd
    string treatment;
    string population;
    double time;
    string clone;
    uint32_t input_order;
  } metadata;
  
  
  //! ---- Constructor / Destructor ---- !//

  //! Constructor.
  cGenomeDiff() : _unique_id_counter(0) { }

  //! Constructor from file
  cGenomeDiff(const string& filename);

  //! Constructor that replaces ::merge(1,2) function
  cGenomeDiff(cGenomeDiff& merge1, cGenomeDiff& merge2, bool new_id=true, bool verbose=false);

  //! Destructor.
  ~cGenomeDiff() { }
  
  
  //!---- Accessors ---- !//

  string get_file_path() const {return _file_path;}
  string get_file_name() const {return path_to_filename(_file_path);}
  
  string get_title() const {return metadata.title;}

  void set_title(const string& in_title) 
  { metadata.title = in_title; replace(metadata.title.begin(), metadata.title.end(), ' ', '_'); }
  
  void add_breseq_data(const key_t &key, const string& value)
    { this->metadata.breseq_data.insert(pair<string,string>(key, value)); }

  string get_breseq_data(const key_t &key)
  { if ( this->metadata.breseq_data.find(key) != this->metadata.breseq_data.end() ) return this->metadata.breseq_data[key]; return ""; }
  
  
  //!---- Input and Output ---- !//
  
  //! Read a genome diff from a file.
  cFileParseErrors read(const string& filename, bool suppress_errors = false);
  
  //! Check to see if genome diff is valid with reference sequences
  cFileParseErrors valid_with_reference_sequences(cReferenceSequences& ref_seq, bool suppress_errors = false);
    
  //! Write the genome diff to a file.
  void write(const string& filename, bool include_unprintable_fields=false);
  
  
  //!---- Adding and Removing Entries ---- !//
  
  //! Helper function to find next unused id
  uint32_t new_unique_id();

  //! Add a new unique id to this entry
  void assign_unique_id_to_entry(cDiffEntry &de);
  
  //! Add an item to this genome diff. Returns pointer to new copy of item.
  diff_entry_ptr_t add(const cDiffEntry& item, bool assign_unique_id=true);
  
  //! Removes an entry, with properly updating mutations that pointed to it
  diff_entry_list_t::iterator remove(diff_entry_list_t::iterator remove_it);
  
  //! Remove mutations, evidence, or validation.
  void remove_group(cGenomeDiff::group group);
  
  //! Remove all of a specific type of entry
  void remove_type(gd_entry_type _type);
  
  void remove_all_but_mutations_and_unknown();
  
  //! Remove all mutations except for deletions of the entire sequence
  void remove_mutations_on_deleted_reference_sequence(const string& seq_id, const int32_t sequence_length);

  
  //!---- Accessing Entries ---- !//
  
  diff_entry_ptr_t find_by_id(string _id);
  
  //! Retrieve cDiffEntrys that match given type(s) 
  const diff_entry_list_t get_const_list() const { return _entry_list; }
  diff_entry_list_t get_list(const vector<gd_entry_type>& types = vector<gd_entry_type>()) const;
  
  // Sometimes we need to remove entries outside of the class methods... use this for that.
  diff_entry_list_t* get_mutable_list_ptr() {return &_entry_list; };
  
  void set_list(diff_entry_list_t& in_list) {  _entry_list = in_list; }
  
  //! retrieve cDiffEntrys that match given type(s) and do not have 'no_show' set
  diff_entry_list_t show_list(const vector<gd_entry_type>& types = vector<gd_entry_type>()) const;
  
  //! Gets parent of entry (entry using it as evidence), if there is one
  // @JEB: should return a list, there can be multiple inheritance
  diff_entry_list_t using_evidence_list(const cDiffEntry& evidence) const;
  
  //! Returns _entry_list with matching item._evidence
  //  (That is, any entries using 'item' as evidence)
  diff_entry_list_t in_evidence_list(const cDiffEntry& item) const;
  
  diff_entry_list_t mutation_list();
  diff_entry_list_t evidence_list();
  diff_entry_list_t validation_list();
  
  //! Removes all GD entries that aren't used as evidence.
  void filter_not_used_as_evidence(bool verbose=false);
  
  //! Remove items used as evidence by any mutations out of input list
  diff_entry_list_t filter_used_as_evidence(const diff_entry_list_t& list) const;
  
  //! Helper function for returning subsets below
  bool mutation_in_entry_of_type(cDiffEntry mut, const gd_entry_type type);
  bool mutation_deleted(cDiffEntry mut) { return mutation_in_entry_of_type(mut, DEL); }
  bool mutation_unknown_or_missing_coverage(cDiffEntry mut) { return (mutation_in_entry_of_type(mut, UN) || mutation_in_entry_of_type(mut, MC)); }

  //! Removes all annotation and other information, leaving only specs
  void strip_to_spec(); // strips all items
  
  //!---- Set Operations ---- !//
  
  //! Subtract mutations using gd_ref as reference.
  void set_subtract(cGenomeDiff& gd_ref, bool phylogeny_id_aware, bool frequency_aware, bool verbose=false);

  void set_intersect(cGenomeDiff& gd_ref, bool verbose=false);
  
  void set_union(const cGenomeDiff& gd_ref, bool evidence_mode, bool phylogeny_aware, bool verbose=false);
  
  //! Merge GenomeDiff information using gd_new as potential new info.
  void merge(const cGenomeDiff& merge_gd, bool reassign_ids=false, bool phylogeny_id_aware = false, bool verbose=false);
  
  //! fast merge, doesn't compare entries, but does renumber
  void merge_preserving_duplicates(const cGenomeDiff& gd);
  
  //! Helper function for fixing IDs after a set operation
  void reassign_unique_ids();
  
  //!---- Sorting Items in Genome Diff ---- !//
  
  static bool diff_entry_ptr_compare_sort(const diff_entry_ptr_t& a, const diff_entry_ptr_t& b);
  static bool diff_entry_ptr_sort(const diff_entry_ptr_t& a, const diff_entry_ptr_t& b);
  static bool diff_entry_ptr_sort_apply_order(const diff_entry_ptr_t& a, const diff_entry_ptr_t& b);
  
  // Normal sort. Used for printing, merging, and compare.
  void sort() { _entry_list.sort(diff_entry_ptr_sort); }

  // Sort -- taking into account 'before' and 'within' tags
  void sort_apply_order();
  
  // Test if entry with first_id is applied before second_id using BEFORE entries
  bool applied_before_id(const string& first_id, const string& second_id);
  
  bool still_duplicates_considering_within(const cDiffEntry& a, const cDiffEntry& b);
  void sort_and_check_for_duplicates(cFileParseErrors* file_parse_errors = NULL);
  
  //! Reconciles mutations that were predicted two different ways due to user evidence, leading to invalid GD
  void reconcile_mutations_predicted_two_ways();
  
  //!---- Simulating and Applying Mutations ---- !//
  
  //! Call to generate random mutations.
  void random_mutations(string exclusion_file,
                        string type,
                        uint32_t n_muts,
                        uint32_t buffer,
                        cReferenceSequences& ref_seq_info,
                        bool verbose = false);

  void mutations_to_evidence(cReferenceSequences &ref_seq, bool remove_mutations = true);
  
  // Helper function for apply_to_sequences
  void shift_positions(cDiffEntry& item, cReferenceSequences& ref_seq_info, bool verbose=false);

  // For constructing the sequence a MOB replaces things with
  string mob_replace_sequence(cReferenceSequences& ref_seq_info, 
                              cDiffEntry& mut, 
                              string* picked_seq_id = NULL, 
                              cFeatureLocation* picked_sequence_feature = NULL
                              );
  
  //! Call to apply Genome Diff to sequences
  void apply_to_sequences(cReferenceSequences &ref_seq_info,
                          cReferenceSequences& new_ref_seq_info,
                          bool verbose=false,
                          int32_t slop_distance=10,
                          int32_t size_cutoff_AMP_becomes_INS_DEL_mutation = kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation
                          );
  
  //! Remove mutations that overlap MASK items in another GD
  void mask_mutations(cGenomeDiff& mask_gd, bool mask_only_small, bool verbose);

  void filter_to_within_region(cReferenceSequences& ref_seq_info, const string& region);
  
  //! Shift mutations to preferred descriptions
  void normalize_mutations(cReferenceSequences &ref_seq, Settings& settings, bool verbose = false);
  
  bool read_counts_for_entry(const cDiffEntry& de, double& new_read_count, double& total_read_count);
  
  //!---- Comparing known lists of mutations/evidence to test files ---- !//
  
  static cGenomeDiff check(cGenomeDiff& ctrl, cGenomeDiff& test, bool verbose = false);
  
  static cGenomeDiff check_evidence(cReferenceSequences& sequence,
                                       uint32_t buffer,
                                       uint32_t shorten_length,
                                       cGenomeDiff& ctrl,
                                       cGenomeDiff& test,
                                       bool jc_only_accepted,
                                       bool verbose = false);
  static void write_jc_score_table(cGenomeDiff& compare, string table_file_path, bool verbose = false); 

  static void tabulate_mutation_frequencies_from_multiple_gds(cGenomeDiff& master_gd,
                                                     vector<cGenomeDiff>& gd_list,
                                                     vector<string> &title_list,
                                                     bool phylogeny_aware = false,
                                                     bool verbose = false);

  
  //!---- Format Conversion Functions: Member ---- !//

  // ! VCF files
  void read_vcf(const string& filename);
  void write_vcf(const string& filename, cReferenceSequences& ref_seq_info);

  //! GVF files
  void write_gvf(const string& filename, cReferenceSequences& ref_seq_info, bool snv_only = false);
  
  //! JSON files
  void write_json(const string &jsonfile);

  
  //!---- Format Conversion Functions: Static Convenience ---- !//

  //! Convert GD to TSV input file
  //!
  static void write_separated_values_file(
                        string& output_file_name,
                        const char* separator,
                        vector<cGenomeDiff>& gd_list,
                        bool preserve_evidence,
                        bool verbose = false
                       );
  
  //! Write a text table
  static void write_table_file(
                      string& output_file_name,
                      const char* separator,
                      cGenomeDiff& gd,
                      const vector<string>& gd_titles,
                      const MutationTableOptions& mutation_table_options
                      );

  
  //! Convert genome diff to GVF
  static void GD2GVF( const string& gdfile, const string& gvffile, cReferenceSequences& ref_seq_info, bool snv_only = false )
    { cGenomeDiff gd(gdfile); gd.write_gvf(gvffile, ref_seq_info, snv_only); }
  
  //! Convert VCF to genome diff
  static void GD2VCF( const string &gdfile, const string & vcffile, cReferenceSequences& ref_seq_info)
    { cGenomeDiff gd(gdfile); gd.write_vcf(vcffile, ref_seq_info); }

  static void VCF2GD( const string& vcffile, const string& gdfile )
  { cGenomeDiff gd; gd.read_vcf(vcffile); gd.write(gdfile); }
  
  //! Convert GD to PHYLIP/FASTA input file
  //! format must be PHYLIP or FASTA
  //! 
  static void write_genotype_sequence_file(
                           const string& format,
                           const string& output_file_name,
                           cGenomeDiff& master_gd, 
                           vector<cGenomeDiff>& gd_list,
                           cReferenceSequences& ref_seq_info,
                           bool missing_as_ancestral = false,
                           bool verbose = false);
  
  
  //! Convert GD to Circos files
  static void GD2Circos(const vector<string> &gd_file_names,
                        const vector<string> &reference_file_names,
                        const string &circos_directory,
                        double distance_scale,
                        double feature_scale);

  //! Convert genome diff to OLI format
  static void GD2OLI( const vector<string> &gd_file_names, 
                      const vector<string> &reference_file_names, 
                      const string& output_file_name,
                      const uint32_t large_size_cutoff,
                      const bool phylogeny_aware);
  
  // For creating coverage graphs in R
  static void GD2COV( const vector<string> &gd_file_names,
                      const vector<string> &reference_file_names,
                      const string& output_file_name,
                      const uint32_t chunk_size = 100);

  //! Functions for dealing with lists of Genome Diffs
  static void sort_gd_list_by_treatment_population_time(vector<cGenomeDiff>& genome_diffs);
  
};
  
  

}
#endif
