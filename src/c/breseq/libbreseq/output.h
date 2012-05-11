#ifndef _BRESEQ_OUTPUT_H_
#define _BRESEQ_OUTPUT_H_
#include "libbreseq/common.h"
#include "libbreseq/settings.h"
#include "libbreseq/reference_sequence.h"
#include "libbreseq/genome_diff.h"
//#include "cgicc/Cgicc.h"
namespace breseq
{

/*-----------------------------------------------------------------------------
 *  Diff_Entry Keywords - Specific to breseq output! 
 *  Others are defined in genome_diff.*
 *-----------------------------------------------------------------------------*/
extern const char* ALIGNMENT_EMPTY_CHANGE_LINE;
extern const char* ALIGNMENT_OVERLAP;
extern const char* BAM_PATH;
extern const char* DELETED;
extern const char* FASTA_PATH;
extern const char* FILE_NAME;
extern const char* FISHER_STRAND_P_VALUE;
extern const char* FLANKING_LEFT;
extern const char* GENES;
extern const char* GENE_NAME;
extern const char* GENE_POSITION;
extern const char* GENE_PRODUCT;
extern const char* GHOST_END;
extern const char* GHOST_SEQ_ID_END;
extern const char* GHOST_SEQ_ID_START;
extern const char* GHOST_START;
extern const char* GHOST_STRAND_END;
extern const char* GHOST_STRAND_START;
extern const char* INSERT_END;
extern const char* INSERT_START;
extern const char* ITEM;
extern const char* KS_QUALITY_P_VALUE;
extern const char* MC_SIDE_1;
extern const char* MC_SIDE_2;
extern const char* NEW_SEQ;
extern const char* NO_SHOW;
extern const char* PLOT;
extern const char* PREFIX;
extern const char* SIZE;
extern const char* TRUNCATE_END;
extern const char* TRUNCATE_START;
extern const char* _COVERAGE_PLOT_FILE_NAME;
extern const char* _EVIDENCE_FILE_NAME;
extern const char* _NEW_JUNCTION_EVIDENCE_FILE_NAME;
extern const char* _SIDE_1_EVIDENCE_FILE_NAME; 
extern const char* _SIDE_2_EVIDENCE_FILE_NAME;
extern const char* SIDE_1_OVERLAP;
extern const char* SIDE_1_JC;
extern const char* SIDE_2_JC;

namespace output
{


/*-----------------------------------------------------------------------------
 *  HTML Attribute Keywords
 *-----------------------------------------------------------------------------*/
extern const char* ALIGN_CENTER;
extern const char* ALIGN_RIGHT;
extern const char* ALIGN_LEFT;

/*-----------------------------------------------------------------------------
 *  Utilities for Encoding HTML
 *-----------------------------------------------------------------------------*/

  //! Wraps input in <i></i> tags which renders as italic text
  inline string i(const string& input) {return "<i>"+input+"</i>";}
  //! Wraps input in <b></b> tags which renders as bold text
  inline string b(const string& input) {return "<b>"+input+"</b>";}
  //! Wraps input in <a></a> tags which defines an anchor,
  //used to create a link to another target document.
  inline string a(const string& target, const string& input) 
    {return "<a href=\"" + target +"\">"+input+"</a>";} 
  //! Wraps input in <th></th> tags	which defines a header cell		
  inline string th(const string& input = "") { return "<th>"+input+"</th>";}
  inline string th(const string& attributes, const string& input) 
    {return "<th " + attributes + ">" + input + "</th>";}
  //! Wraps input in <td></td> tags	which defines a standar cell
  inline string start_td(const string& attributes)
    {return "<td " + attributes + ">";}  
  inline string td(const string& attributes, const string& input)
    {return "<td " + attributes + ">" + input + "</td>";}
  inline string td(const string& input ="")
  {return "<td>" + input + "</td>";}
  //! Wraps input in <tr></tr> tags which define a row in an HTML table	
  inline string start_tr(const string& attributes = "") 
    {return "<tr " + attributes + ">";}
  inline string tr(const string& input) {return "<tr>" + input + "</tr>";}
  inline string tr(const string& attributes, const string& input)
    {return start_tr(attributes) + input + "</tr>";}
  inline string start_table(const string& attributes)
    {return "<table " + attributes + ">";}
  inline string end_table()
    {return "</table>";}
  //! Wraps input in <font></font> tags
  inline string font(const string& attributes, const string& input) 
    {return "<font " + attributes + ">" + input + "</font>";}
  inline string html_footer()
    {return "</body></html>";}
  inline string div(const string& attributes, const string& input)
    {return "<div " + attributes + ">" + input + "</div>";}
  inline string img(const string& target)
    {return "<img src=\"" + target + "\" />";}
  inline string h1(const string & input)
    {return "<h1>" + input + "</h1>";} 
  inline string h2(const string & input)
    {return "<h2>" + input + "</h2>";} 
  
  //! Encodes dash, en dash and spaces to HTML
  string nonbreaking(const string& input);
  //! Encodes en dash
  string htmlize(const string& input);
  //! Adds commas to large numbers (ex 1000 to 1,000)
  string commify(const string& input);

  //Specific to Breseq HTML Files
  string html_header(const string& title, const Settings& settings);
  string breseq_header_string(const Settings& settings);

/*-----------------------------------------------------------------------------
 * HTML FILES 
 *-----------------------------------------------------------------------------*/
// Convenience structure for passing many options
struct MutationTableOptions {
  
  MutationTableOptions() 
  : repeat_header(0)
  , legend_row(false)
  , one_ref_seq(false)
  , shade_frequencies(false)
  {}
  
  uint32_t repeat_header;
  bool legend_row; 
  bool one_ref_seq;
  bool shade_frequencies;
  vector<string> gd_name_list_ref;
  string relative_link;
  
};
  
void html_index(const string& file_name, const Settings& settings, Summary& summary,
                cReferenceSequences& ref_seq_info, cGenomeDiff& gd);
void mark_gd_entries_no_show(const Settings& settings, cGenomeDiff& gd);
  
void html_marginal_predictions(const string& file_name, const Settings& settings, Summary& summary,
                               cReferenceSequences& ref_seq_info, cGenomeDiff& gd);
void html_statistics(const string& file_name, const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info);
  
void html_compare(
                  const Settings& settings,
                  const string &file_name, 
                  const string &title, 
                  cGenomeDiff& gd,
                  MutationTableOptions& mt_options
                  ); 
  
/*-----------------------------------------------------------------------------
 * HTML TABLE STRINGS
 *-----------------------------------------------------------------------------*/
  

  
struct Html_Mutation_Table_String : public string
{
  public:
    //!Constructors
    Html_Mutation_Table_String(
                               const Settings& settings,
                               cGenomeDiff& gd,
                               diff_entry_list_t& list_ref,
                               MutationTableOptions& options
                               );
    
   
    //! Main Build Object
    //!Factory Methods
    void Header_Line(bool print_main_header = true);
    void Item_Lines();
    //!Helper Functions
    string freq_to_string(const string& freq);//!< Used in Item_Lines()
    string freq_cols(vector<string> freq_list);//!< Used in Item_Lines()
    size_t total_cols; //!< Shared between Factory Methods, set in Header_Line()

    //!Parameters
    Settings settings;
    cGenomeDiff gd;
    diff_entry_list_t list_ref;
    MutationTableOptions options;
};


string html_missing_coverage_table_string(
                                          diff_entry_list_t& list_ref,
                                          bool show_details,
                                          const string& title = "Missing coverage evidence...",
                                          const string& relative_link=""
                                          );

string html_read_alignment_table_string(
                                        diff_entry_list_t& list_ref,
                                        bool show_details,
                                        const string& title = "Read alignment evidence...",
                                        const string& relative_link = ""
                                        );

string html_new_junction_table_string(diff_entry_list_t& jc,
                                      bool show_details,
                                      const string& title = "New junction evidence",
                                      const string& relative_link = ""
                                      );
  
string html_copy_number_table_string(
                                     diff_entry_list_t& list_ref, 
                                     bool show_details, 
                                     const string& title = "Copy number evidence", 
                                     const string& relative_link = ""
                                     );


string html_genome_diff_item_table_string(const Settings& settings, cGenomeDiff& gd, 
                                        diff_entry_list_t& list_ref);
string html_deletion_coverage_values_table_string(const Settings& settings, cReferenceSequences& ref_seq_info, Summary& summary);
/*-----------------------------------------------------------------------------
 * Helper Functions For Tables 
 *-----------------------------------------------------------------------------*/
string formatted_mutation_annotation(const cDiffEntry& mut);
string to_underline_red_codon(const cDiffEntry& mut,const string& codon_key);
string decode_reject_reason(const string & reject);

/*-----------------------------------------------------------------------------
 *  Create_Evidence_Files
 *-----------------------------------------------------------------------------*/
struct Evidence_Files
{
  class Evidence_Item : public cDiffEntry
  {
  public:
    Evidence_Item(diff_entry_map_t& _fields, diff_entry_ptr_t _item, diff_entry_ptr_t _parent_item)
    : cDiffEntry(_fields), item(_item), parent_item(_parent_item) {};
    
    diff_entry_ptr_t item;
    diff_entry_ptr_t parent_item;
  };

  Evidence_Files(const Settings& settings, cGenomeDiff& gd);
  
  vector<Evidence_Item> evidence_list;
  
  private:
  
    string html_evidence_file_name(Evidence_Item& evidence_item);
    void add_evidence(const string& file_name, diff_entry_ptr_t item,
                      diff_entry_ptr_t parent_item, map<string,string>& fields);
    string file_name(Evidence_Item& evidence_item);
    void html_evidence_file(const Settings& settings, cGenomeDiff& gd, Evidence_Item& item);

};

/*-----------------------------------------------------------------------------
 *  Output Utilities
 *-----------------------------------------------------------------------------*/

// sub draw_coverage
void draw_coverage(Settings& settings, cReferenceSequences& ref_seq_info, cGenomeDiff& gd);


}// end output namespace 
}// end breseq namespace
#endif
