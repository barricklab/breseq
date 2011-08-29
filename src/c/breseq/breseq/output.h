#ifndef _BRESEQ_OUTPUT_H_
#define _BRESEQ_OUTPUT_H_
#include "breseq/common.h"
#include "breseq/settings.h"
#include "breseq/annotated_sequence.h"
#include "breseq/genome_diff.h"
namespace breseq
{

/*-----------------------------------------------------------------------------
 *  TEMPORARY STRUCTS, STILL NEED IMPLEMENTATION.
 *-----------------------------------------------------------------------------*/
struct Options {
  bool repeat_header;
};

/*-----------------------------------------------------------------------------
 *  Diff_Entry Keywords 
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
// For JC
extern const char* SIDE_1_OVERLAP;
extern const char* SIDE_1_POSITION;
extern const char* SIDE_1_SEQ_ID;
extern const char* SIDE_1_STRAND;
extern const char* SIDE_2_POSITION;
extern const char* SIDE_2_SEQ_ID;
extern const char* SIDE_2_STRAND;
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
  //! Wraps input in <font></font> tags
  inline string font(const string& attributes, const string& input) 
    {return "<font " + attributes + ">" + input + "</font>";}
  inline string html_footer()
    {return "</html>";}
  inline string div(const string& attributes, const string& input)
    {return "<div " + attributes + ">" + input + "</div>";}
  inline string img(const string& target)
    {return "<img src=\"" + target + "\" />";}
  inline string h1(const string & input)
    {return "<h1>" + input + "</h1>";} 
  
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
void html_index(const string& file_name, const Settings& settings, Summary& summary,
                cReferenceSequences& ref_seq_info, genome_diff& gd);
void html_marginal_predictions(const string& file_name, const Settings& settings, Summary& summary,
                               cReferenceSequences& ref_seq_info, genome_diff& gd);
void html_statistics(const string& file_name, const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info);
void html_compare(const Settings& settings,const string &file_name, const string &title, 
                  genome_diff& gd, bool one_ref_seq, vector<string>& gd_name_list_ref, Options& options); 
void html_compare_polymorphisms(const Settings& settings, const string& file_name, const string& title,
                                vector<string>& list_ref);
/*-----------------------------------------------------------------------------
 * HTML TABLE STRINGS
 *-----------------------------------------------------------------------------*/
struct Html_Mutation_Table_String : public string
{
  public:
    //!Constructors
    Html_Mutation_Table_String(
                               const Settings& settings,
                               genome_diff& gd,
                               diff_entry_list& list_ref,
  			                       vector<string>& gd_name_list_ref,
                               Options& options,
                               bool legend_row = false, 
                               bool one_ref_seq = false,
  			                       const string& relative_link = "" 
                               );
    
    Html_Mutation_Table_String(
                               const Settings& settings,
                               genome_diff& gd,
                               diff_entry_list& list_ref,
  			                       const string& relative_path = "", 
                               bool legend_row = false, 
                               bool one_ref_seq = false
                               );
    
   
    //! Main Build Object
    //!Factory Methods
    void Header_Line();
    void Item_Lines();
    //!Helper Functions
    string freq_to_string(const string& freq);//!< Used in Item_Lines()
    string freq_cols(vector<string> freq_list);//!< Used in Item_Lines()
    size_t total_cols; //!< Shared between Factory Methods, set in Header_Line()

    //!Parameters
    Settings settings;
    genome_diff gd;
    diff_entry_list list_ref;
    bool legend_row;
    bool one_ref_seq;
    vector<string> gd_name_list_ref; 
    Options options;
    string relative_link;
};


string html_missing_coverage_table_string
  (diff_entry_list& list_ref, 
   bool show_reject_reason,
   const string& title = "Missing coverage evidence...",
   const string& relative_link="");

string html_read_alignment_table_string  
  (diff_entry_list& list_ref, 
   bool show_reject_reason,
   const string& title = "Read alignment evidence...",
   const string& relative_link = "");

string html_new_junction_table_string
  (diff_entry_list& jc,
   bool show_reject_reason,
   const string& title= "New junction evidence",
   const string& relative_link = "");

string html_genome_diff_item_table_string(const Settings& settings, genome_diff& gd, 
                                        diff_entry_list& list_ref);
/*-----------------------------------------------------------------------------
 * Helper Functions For Tables 
 *-----------------------------------------------------------------------------*/
string formatted_mutation_annotation(const diff_entry& mut);
string to_underline_red_codon(const diff_entry& mut,const string& codon_key);
string decode_reject_reason(const string & reject);

/*-----------------------------------------------------------------------------
 *  Create_Evidence_Files
 *-----------------------------------------------------------------------------*/
struct Evidence_Files
{
  struct Evidence_Item:diff_entry
  {
  	diff_entry parent_item;
  	diff_entry item;
  };

  struct Type_Not_Equal
  {
    string type;
    Type_Not_Equal(string in_type)
      : type(in_type)
      {}
    bool operator() (diff_entry_ptr item) { return (*item)._type != type;} 
  };

  Evidence_Files(const Settings& settings, genome_diff& gd);
  
  vector<Evidence_Item> evidence_list;
  
  private:
  
    string html_evidence_file_name(Evidence_Item& evidence_item);
    void add_evidence(const string& file_name, diff_entry& item,
                      diff_entry& parent_item, map<string,string> fields);
    string file_name(Evidence_Item& evidence_item);
    void html_evidence_file(const Settings& settings, genome_diff& gd, Evidence_Item& item);

};

/*-----------------------------------------------------------------------------
 *  Output Utilities
 *-----------------------------------------------------------------------------*/

// sub draw_coverage
void draw_coverage(Settings& settings, cReferenceSequences& ref_seq_info, genome_diff& gd);


}// end output namespace 
}// end breseq namespace
#endif
