#ifndef _BRESEQ_OUTPUT_H_
#define _BRESEQ_OUTPUT_H_
#include "breseq/common.h"
#include "breseq/settings.h"
#include "breseq/annotated_sequence.h"
#include "breseq/genome_diff.h"
#include "breseq/settings.h"

namespace breseq
{
///Temporary structs 
struct Options {
  bool repeat_header;
};
struct Mutation{};
struct Interval{};
struct Reference{};
class tSettings:Settings{
  string print_run_name;
};

//TEST
/// sub html_index
///         my (string file_name, Settings settings, Summary summary, $ref_seq_info, $gd) = @_;
void html_index(string file_name, Settings settings, Summary summary,
                cReferenceSequences ref_seq_info, genome_diff gd);
// sub html_marginal_predictions
//         my (string file_name, Settings settings, $summary, $ref_seq_info, $gd) = @_;
void html_marginal_predictions(string file_name, Settings settings,Summary summary,
                               cReferenceSequences ref_seq_info, genome_diff gd);
// sub html_header
//         my ($title) = @_;
void html_header(string title);
// sub html_footer
// sub html_compare
//         my (Settings settings, string file_name, $title, $gd, $one_ref_seq, $gd_name_list_ref, $options) = @_;
void html_compare(Settings settings, string file_name, string title, 
                  genome_diff gd, bool one_ref_seq, vector<string> gd_name_list_ref, Options options); 
//? string gd_name_list_ref
// if (defined $gd_name_list_ref)
//                 @freq_header_list = @$gd_name_list_ref;
//? Options options
//$output_str.= $header_str if (($row_num != 0) && (defined $options->{repeat_header}) && ($row_num % $options->{repeat_header} == 0));

// sub html_compare_polymorphisms
//         my (Settings settings, string file_name, $title, $list_ref) = @_;
void html_compare_polymorphisms(Settings settings, string file_name, string title,
                                vector <string> list_ref);
//? vector <string> list_ref
// sub html_statistics
//         my (string file_name, Settings settings, $summary, $ref_seq_info) = @_;
void html_statistics(string file_name, Settings settings, Summary summary, 
                     cReferenceSequences ref_seq_info);
/// sub breseq_header_string
///         my (Settings settings) = @_;
string breseq_header_string(Settings settings);

// sub html_genome_diff_item_table_string
//         my (Settings settings, $gd, $list_ref) = @_;
void html_genome_diff_item_table_string(Settings settings, genome_diff gd, 
                                        vector <string> list_ref);
//? vector <string> list_ref

// sub formatted_mutation_annotation
//         my ($mut) = @_;
void formatted_mutation_annotation(Mutation mut);
//                 sub to_underline_red_codon
//                         my ($mut, $codon_key) = @_;

/// sub html_mutation_table_string
///         our (Settings settings, $gd, $list_ref, $relative_link, $legend_row, $one_ref_seq, $gd_name_list_ref, $options) = @_;
//html_mutation_table_string($settings, $gd, \@muts, undef, undef, $one_ref_seq, $gd_name_list_ref, $options)
string html_mutation_table_string(Settings settings, genome_diff gd, genome_diff::entry_list_t muts,
                                  bool legend_row, bool one_ref_seq, vector<string> gd_name_list_ref,
                                  Options options) 
{return "html_mutation_table_string:needs implementaion";}
//html_mutation_table_string($settings, $gd, \@muts, $relative_path, undef, $one_ref_seq)
string html_mutation_table_string(Settings settings, genome_diff gd, genome_diff::entry_list_t muts, 
                                  string relative_path, bool legend_row, bool one_ref_seq)
{return "html_mutation_table_string:needs implementation";}
//html_mutation_table_string($settings, $gd, $list_ref)
string html_mutation_table_string(Settings settings, genome_diff gd, vector <string> list_ref)
{return "html_mutation_table_string:needs implementaion";}
                               
/// sub html_read_alignment_table_string
///         my ($list_ref, $relative_link, $title, $show_reject_reason) = @_;
void html_read_alignment_table_string(vector <string> list_ref, string relative_link,
                                      string title, bool show_reject_reason);
/// sub html_missing_coverage_table_string
///         my ($list_ref, $relative_link, $title, $show_reject_reason) = @_;
//html_missing_coverage_table_string($list_ref, undef, undef, 1)
string html_missing_coverage_table_string(vector <string> list_ref, string relative_link, string title,
                                        bool show_reject_reason);
//html_missing_coverage_table_string(\@mc, $relative_path, "Unassigned missing coverage evidence...")
string html_missing_coverage_table_string(genome_diff::entry_list_t muts, string relative_path, string title)
{return "needs implementaion";}

/// sub html_new_junction_table_string
///   our ($list_ref, $relative_link, $title, $show_reject_reason) = @_;

//html_new_junction_table_string($list_ref, undef, undef, 1)
string html_new_junction_table_string(vector<cReferenceSequences> list_ref, string relative_link, string title,
                                    bool show_reject_reason);
//html_new_junction_table_string(\@jc, $relative_path, "Marginal new junction evidence...")
//html_new_junction_table_string(\@jcu, $relative_path, "Unassigned new junction evidence...")
string html_new_junction_table_string(genome_diff::entry_list_t jc, string relative_path, string title);


// sub html_evidence_file_name
//         my ($interval) = @_;
void html_evidence_file_name(Interval interval);
// sub html_evidence_file
//         my (Settings settings, $gd, $interval) = @_;
void  html_evidence_file(Settings settings, genome_diff gd, Interval interval);
// sub decode_reject_reason
//         my ($reject) = @_;
void decode_reject_reason(string reject);
// sub create_evidence_files
//         my (Settings settings, $gd) = @_;
//         sub add_evidence
//                 my ($evidence_file_name_key, $evidence_item) = @_;
void create_evidence_files(Settings settings, genome_diff gd);
// sub save_text_deletion_file
//         my ($deletion_file_name, $deletions_ref) = @_;
void save_text_deletion_file(string deletion_file_name, vector<Reference> deletions_ref);
// sub draw_coverage
//         my (Settings settings, $ref_seq_info, $gd) = @_;
void draw_coverage(Settings settings, cReferenceSequences ref_seq_info, genome_diff gd);
// sub record_time
//         my ($name) = @_;
void record_time(string name);
// sub time2string
//     my ($seconds) = @_;
void time2string(uint32_t seconds);
// sub save_statistics
//         my (string file_name, $data) = @_;
void save_statistics(string file_name, string data);
// sub load_statistics
//         my (string file_name) = @_;
void load_statistics(string file_name);
// sub nonbreaking
//         my ($text) = @_;
string nonbreaking(string text);
// sub htmlize
//         my ($text) = @_;
string htmlize(string text);
  
}// end breseq namespace
#endif