#include "breseq/output.h"
#include "breseq/anyoption.h"
#include "breseq/alignment_output.h"
#include "breseq/coverage_output.h"

using namespace std;
namespace breseq
{

/*-----------------------------------------------------------------------------
 *  Diff_Entry Keywords 
 *-----------------------------------------------------------------------------*/
const char* ALIGNMENT_EMPTY_CHANGE_LINE="alignment_empty_change_line";
const char* ALIGNMENT_OVERLAP="alignment_overlap";
const char* BAM_PATH="bam_path";
const char* DELETED="deleted";
const char* FASTA_PATH="fasta_path";
const char* FILE_NAME="file_name";
const char* FISHER_STRAND_P_VALUE="fisher_strand_p_value";
const char* FLANKING_LEFT="flanking_left";
const char* GENES="genes";
const char* GENE_NAME="gene_name";
const char* GENE_POSITION="gene_position";
const char* GENE_PRODUCT="gene_product";
const char* GHOST_END="ghost_end";
const char* GHOST_SEQ_ID_END="ghost_seq_id_end";
const char* GHOST_SEQ_ID_START="ghost_seq_id_start";
const char* GHOST_START="ghost_start";
const char* GHOST_STRAND_END="ghost_strand_end";
const char* GHOST_STRAND_START="ghost_strand_start";
const char* INSERT_END="insert_end";
const char* INSERT_START="insert_start";
const char* ITEM="item";
const char* KS_QUALITY_P_VALUE="ks_quality_p_value";
const char* MC_SIDE_1="mc_side_1";
const char* MC_SIDE_2="mc_side_2";
const char* NEW_SEQ="new_seq";
const char* NO_SHOW="no_show";
const char* PLOT="plot";
const char* PREFIX="prefix";
const char* SIZE="size";
const char* TRUNCATE_END="truncate_end";
const char* TRUNCATE_START="truncate_start";
const char* _COVERAGE_PLOT_FILE_NAME="_coverage_plot_file_name";
const char* _EVIDENCE_FILE_NAME="_evidence_file_name";
const char* _NEW_JUNCTION_EVIDENCE_FILE_NAME="_new_junction_evidence_file_name";
const char* _SIDE_1_EVIDENCE_FILE_NAME="_side_1_evidence_file_name";
const char* _SIDE_2_EVIDENCE_FILE_NAME="_side_2_evidence_file_name";
//For JC
const char* SIDE_1_OVERLAP="side_1_overlap";
const char* SIDE_1_POSITION="side_1_position";
const char* SIDE_1_SEQ_ID="side_1_seq_id";
const char* SIDE_1_STRAND="side_1_strand";
const char* SIDE_2_POSITION="side_2_position";
const char* SIDE_2_SEQ_ID="side_2_seq_id";
const char* SIDE_2_STRAND="side_2_strand";
const char* SIDE_1_JC="side_1_jc";
const char* SIDE_2_JC="side_2_jc";
const char* SIDE_2_OVERLAP="side_2_overlap";
namespace output
{

/*
 * =====================================================================================
 *        Class:  HTML
 *  Description:  
 * =====================================================================================
 */
/*-----------------------------------------------------------------------------
 *  HTML Attribute Keywords
 *-----------------------------------------------------------------------------*/
const char* ALIGN_CENTER="align=\"center\"";
const char* ALIGN_RIGHT="align=\"right\"";
const char* ALIGN_LEFT="align=\"left\"";

/*-----------------------------------------------------------------------------
 *  HTML Utility for printing numbers
 *-----------------------------------------------------------------------------*/
string commify(const string& input)
{
  string retval;
  string temp = input;
  reverse(temp.begin(), temp.end());
  for (size_t i = 0; i < input.size(); i++) {
    if ((i+1)%3 > 0) {
      retval.push_back(temp[i]);
    } else  {
      retval.push_back(temp[i]);
      if (i+1 != input.size()) {
        retval.push_back(',');
      }
    }
  }
  reverse(retval.begin(),retval.end());
   return retval;
}
/*-----------------------------------------------------------------------------
 *  HTML Utility for Encoding HTML
 *-----------------------------------------------------------------------------*/
string nonbreaking(const string& input)
{   
  string retval = input;
  
  /* substitute nonbreaking en dash */
  retval = substitute(retval, "–", "&#8211;");
  
  /* substitute nonbreaking hyphen */
  retval = substitute(retval, "-", "&#8209;");

  /* substitute nonbreaking space */
  retval = substitute(retval, " ", "&nbsp;");

  return retval;
}
/*-----------------------------------------------------------------------------
 *  HTML Utility for Encoding HTML
 *-----------------------------------------------------------------------------*/
string htmlize (const string& input) 
{
  string retval = input;
    
  /* substitute nonbreaking en dash */
  retval = substitute(retval, "–", "&#8211;");

  return retval;
}

/*-----------------------------------------------------------------------------
 *  These style definitions are included between the HTML <head> 
 *  tags of every genereated .html page.
 *-----------------------------------------------------------------------------*/
string header_style_string() 
{
  stringstream ss(ios_base::out | ios_base::app);
  ss << "body {font-family: sans-serif; font-size: 11pt;}"                 << endl;
  ss << "th {background-color: rgb(0,0,0); color: rgb(255,255,255);}"      << endl;
  ss << "table {background-color: rgb(1,0,0); color: rgb(0,0,0);}"         << endl;
  ss << "tr {background-color: rgb(255,255,255);}"                         << endl;
  ss << ".mutation_in_codon {color:red; text-decoration : underline;}}"    << endl;
  ss << ".mutation_header_row {background-color: rgb(0,130,0);}"           << endl;
  ss << ".read_alignment_header_row {background-color: rgb(255,0,0);}"     << endl;
  ss << ".missing_coverage_header_row {background-color: rgb(0,100,100);}" << endl;
  ss << ".new_junction_header_row {background-color: rgb(0,0,155);}"       << endl;
  ss << ".alternate_table_row_0 {background-color: rgb(255,255,255);}"     << endl;
  ss << ".alternate_table_row_1 {background-color: rgb(230,230,245);}"     << endl;
  ss << ".polymorphism_table_row {background-color: rgb(160,255,160);}"    << endl;
  ss << ".highlight_table_row {background-color: rgb(192,255,255);}"       << endl;
  ss << ".reject_table_row {background-color: rgb(255,200,165);}"          << endl;
  ss << ".information_table_row {background-color: rgb(200,255,255);}"     << endl;
  ss << ".junction_repeat {background-color: rgb(255,165,0)}"              << endl;
  ss << ".junction_gene {}" << endl;
  
return ss.str();
}



void html_index(const string& file_name, const Settings& settings, Summary& summary,
                cReferenceSequences& ref_seq_info, genome_diff& gd)
{
  (void)summary; //TODO: unused?
  
  // Create Stream and Confirm It's Open
  ofstream HTML(file_name.c_str());

  if (!HTML.good()) {
    cerr << "Could not open file: " << file_name << endl;
    assert(HTML.good());
  }
  
  // Build HTML Head
  HTML << html_header("BRESEQ :: Mutation Predictions", settings);
  HTML << breseq_header_string(settings) << endl;
  HTML << "<p>" << endl;
// #   ###
// #   ## Mutation predictions
// #   ###
  HTML << "<!--Mutation Predictions -->" << endl;
  diff_entry_list muts = gd.list(make_list<string>(SNP)(INS)(DEL)(SUB)(MOB)(AMP));
  
  string relative_path = settings.local_evidence_path;
  
  if(!relative_path.empty())
    relative_path += "/";
  
  //Determine if more than one reference sequence is used
  bool one_ref_seq;

  if (get_keys<string,string>(ref_seq_info.ref_strings).size() == 1)
    one_ref_seq = true;
  else
    one_ref_seq = false;

  //Build Mutation Predictions table
  HTML << "<p>" << endl;
  HTML << Html_Mutation_Table_String(settings, gd, muts, relative_path, false, one_ref_seq) << endl;

// #   ###
// #   ## Unassigned evidence
// #   ###
  HTML << "<!--Unassigned evidence-->" << endl;
  
  diff_entry_list mc = gd.filter_used_as_evidence(gd.list(make_list<string>("MC")));
  
  if (mc.size() > 0) {
    HTML << "<p>" << html_missing_coverage_table_string(mc, false, "Unassigned missing coverage evidence", relative_path);
  }

  diff_entry_list jc = gd.filter_used_as_evidence(gd.list(make_list<string>("JC")));

  jc.remove_if(diff_entry::field_exists("no_show"));
  
  //Don't show junctions for circular chromosomes
  if (!settings.hide_circular_genome_junctions) {
    jc.remove_if(diff_entry::field_exists("circular_chromosome")); 
  }
   
  diff_entry_list jcu = jc;
  jcu.remove_if(diff_entry::field_exists("reject"));

  if (jcu.size() > 0) {
    HTML << "<p>" << endl;
    HTML << html_new_junction_table_string(jcu, false, "Unassigned new junction evidence...", relative_path);
  }

  HTML << "</html>";
  HTML.close();
}



void html_marginal_predictions(const string& file_name, const Settings& settings,Summary& summary,
                               cReferenceSequences& ref_seq_info, genome_diff& gd)
{
  (void)summary; //TODO: unused?

  
  // Create Stream and Confirm It's Open
  ofstream HTML(file_name.c_str());
  
  if(!HTML.good()) {
    cerr << "Could not open file: " <<  file_name << endl;
    assert(HTML.good());
  }

  // Build HTML Head
  HTML << html_header("BRESEQ :: Marginal Predictions",settings); 
  HTML << breseq_header_string(settings) << endl;
  HTML << "<p>" << endl;
  
  string relative_path = settings.local_evidence_path;
  
  if (!relative_path.empty())
    relative_path += "/";
  
  //Determine if more than one reference sequence is used
  bool one_ref_seq;
  if (get_keys<string,string>(ref_seq_info.ref_strings).size() == 1)
    one_ref_seq = true;
  else
    one_ref_seq = false; 

// #   ###
// #   ## Marginal evidence
// #   ###
  diff_entry_list ra = gd.filter_used_as_evidence(gd.list(make_list<diff_entry::key_t>("RA")));

// #   ## don't print ones that overlap predicted deletions or were marked to not show
// #   @ra = grep { !$_->{deleted} && !$_->{no_show} } @ra;
  //Don't print ones that overlap predicted deletions or were marked to not show.
  ra.remove_if(diff_entry::fields_exist(make_list<diff_entry::key_t>("deleted")("no_show")));
  
  if (ra.size() > 0) {
    HTML << "<p>" << endl;
    HTML << html_read_alignment_table_string(ra, false, "Marginal read alignment evidence...",
      relative_path) << endl;
  }

//
//
//
// #   
// #   my @jc = $gd->filter_used_as_evidence($gd->list('JC'));
    diff_entry_list jc = gd.filter_used_as_evidence(gd.list(make_list<string>("JC")));
// #   @jc = grep { !$_->{no_show} } @jc;
    jc.remove_if(diff_entry::field_exists("no_show"));
// #   @jc = grep { $_->{reject} } @jc;
    jc.remove_if(not1(diff_entry::field_exists("reject")));
// #   if (scalar @jc > 0)
// #   { 
// #     ## sort by score, not by position (the default order)...
// #     @jc = sort { -($a->{pos_hash_score} <=> $b->{pos_hash_score}) || -($a->{min_overlap_score} <=> $b->{min_overlap_score})  || ($a->{total_reads} <=> $a->{total_reads}) } @jc;
// #     print HTML p . html_new_junction_table_string(\@jc, $relative_path, "Marginal new junction evidence..."); 
// #   }
     if (jc.size()) {
       //Sort by score, not by position (the default order)...
       jc.sort(diff_entry::by_scores(
         make_list<diff_entry::key_t>("pos_hash_score")("min_overlap_score")("total_reads"))); 
      
       HTML << "<p>" << endl;
       HTML << html_new_junction_table_string(jc, false, "Marginal new junction evidence...", relative_path);
     }

    HTML <<  "</HTML>";
    HTML.close();
}

string html_header (const string& title, const Settings& settings)
{
  stringstream ss(ios_base::out | ios_base::app);  
  
  ss << "<html>" << endl;   
  ss << "<title>" << title;
  if (!settings.print_run_name.empty()) {
    ss << " :: " << settings.print_run_name;
  }
  ss << "</title>" << endl;
  
  ss << "<head>" << endl;
  ss << "<style type = \"text/css\">" << endl;
  ss << header_style_string() << endl;
  ss << "</style>" << endl;
  ss << "</head>" << endl;
 
  return ss.str();
}



void html_compare(Settings& settings,const string &file_name, const string &title, genome_diff& gd,
                  bool one_ref_seq, vector<string>& gd_name_list_ref, Options& options)
{
  // Create stream and confirm it's open
  ofstream HTML(file_name.c_str());
  
  if(!HTML.good()) {
    cerr << "Could not open file: " <<  file_name << endl;
    assert(HTML.good());
  }

  //Build html head
  HTML << "<html>" << endl;   
  HTML << "<title>" << title << "</title>" << endl;
  HTML << "<head>" << endl;
  HTML << "<style type=\"text/css\">" << endl;
  HTML << header_style_string() << endl;
  HTML << "</style>" << endl;
  HTML << "</head>" << endl;
  
  diff_entry_list muts = gd.mutation_list();

  HTML << Html_Mutation_Table_String(settings, gd, muts, gd_name_list_ref, options, false, one_ref_seq, "");
  HTML << "</html>";
  HTML.close();
}

void html_compare_polymorphisms(Settings& settings, const string& file_name, const string& title, diff_entry_list& list_ref)
{
  (void)settings; //TODO: unused?

  // Create stream and confirm it's open
  ofstream HTML(file_name.c_str());
  
  if(!HTML.good()) {
    cerr << "Could not open file: " <<  file_name << endl;
    assert(HTML.good());
  }

  //Build html head
  HTML << "<html>" << endl;
  HTML << "<title>" << title << "</title>" << endl;
  HTML << "<head>" << endl;
  HTML << "<style type =\"text/css\">" << endl;
  HTML << header_style_string() << endl;
  HTML << "</style>" << endl;
  HTML << "</head>" << endl;
  HTML << html_read_alignment_table_string(list_ref, true); 
  HTML << "</html>" << endl;
  HTML.close();
}
  

void html_statistics(const string &file_name, const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info)
{  
  // Create stream and confirm it's open
  ofstream HTML(file_name.c_str());
  
  if(!HTML.good()) {
    cerr << "Could not open file: " <<  file_name << endl;
    assert(HTML.good());
  }

  //Build html head
  HTML << html_header("BRESEQ :: Summary Statistics", settings);
  HTML << breseq_header_string(settings) << endl;
  HTML << "<p>" << endl;

  //Write read file information
  //HTML << "<!-- Write fastq read file informations -->" << endl;
  HTML << "<table border=\"0\" cellspace=\"1\" cellpadding=\"5\">" << endl;
  HTML << "<tr>" << th() << th("fastq read file") << th("reads") << 
                    th("bases") << th("longest") << "</tr>" << endl;
  for(cReadFiles::const_iterator it=settings.read_files.begin(); it!=settings.read_files.end(); it++)
  {
    const AnalyzeFastq& s = summary.sequence_conversion.reads[it->m_base_name];
    
    HTML << "<tr>";
    HTML << td( a(Settings::relative_path( 
                      settings.file_name(settings.error_rates_plot_file_name, "#", it->m_base_name), settings.output_path
                                          ), 
                "errors" 
                )
              );
    HTML << td(it->m_base_name);
    HTML << td(ALIGN_RIGHT, commify(to_string(s.num_reads)));
    HTML << td(ALIGN_RIGHT, commify(to_string(s.num_bases)));
    HTML << td(ALIGN_RIGHT, to_string(s.max_read_length) + "&nbsp;bases");

    HTML << "</tr>";
  }
  
  HTML << "<tr class=\"highlight_table_row\">";
  HTML << td();
  HTML << td(b("total"));
  HTML << td(ALIGN_RIGHT , b(commify(to_string(summary.sequence_conversion.num_reads))) );
  HTML << td(ALIGN_RIGHT , b(commify(to_string(summary.sequence_conversion.num_bases))) );
  HTML << td(b(commify(to_string(summary.sequence_conversion.max_read_length))) + "&nbsp;bases");
  HTML << "</tr></table>";
  
  //Write reference sequence information
  //HTML << "<!-- Write reference sequence information -->" << endl;
  HTML << "<p>" << endl;
  HTML << "<table border=\"0\" cellspacing=\"1\" cellpadding=\"5\" >" << endl;
  HTML << "<tr>" << th() << 
                    th() << 
                    th("reference sequence") << 
                    th("length") << 
                    th(ALIGN_LEFT, "description") << 
          "</tr>" << endl;
             
  size_t total_length = 0;
  for(cReferenceSequences::iterator it=ref_seq_info.begin(); it!=ref_seq_info.end(); it++)
  {
    total_length += it->m_length;
    
    HTML << td( a(Settings::relative_path( 
                                          settings.file_name(settings.overview_coverage_plot_file_name, "@", it->m_seq_id), settings.output_path
                                          ), 
                  "coverage" 
                  )
               );
    HTML << td( a(Settings::relative_path( 
                                          settings.file_name(settings.unique_only_coverage_plot_file_name, "@", it->m_seq_id), settings.output_path
                                          ), 
                  "distribution" 
                  )
               ); 
    HTML << td(it->m_seq_id);
    HTML << td(ALIGN_RIGHT, commify(to_string(it->m_length)));
    HTML << td(it->m_definition);
  }  
  
  HTML << "<tr class=\"highlight_table_row\">";
  HTML << td();
  HTML << td();
  HTML << td(b("total"));
  HTML << td(ALIGN_RIGHT, b(commify(to_string(total_length))) );
  HTML << td();
  HTML << "</tr>" << endl;

// # //TODO @JEB Summary
// #   ## junction only reference sequences
// #   foreach my $seq_id (@{$ref_seq_info->{junction_only_seq_ids}})
// #   {
// #     my $c = $summary->{sequence_conversion}->{reference_sequences}->{$seq_id};
// #     print HTML Tr(
// #       td({-colspan=>"2", -align=>"center"}, "junction&nbsp;only"), 
// #       td($seq_id), 
// #       td({-align=>"right"},commify($c->{length})), 
// #       td($c->{definition})
// #     );
// #     $total_length+= $c->{length};
// #   }
// #   
// #   print HTML end_table();
  HTML << "</table>" << endl;  
  
  // Write Execution Times
  const vector<ExecutionTime>& times = settings.execution_times;
  // HTML << "<!-- Write Times -->" << endl;
  HTML << "<p>"  << endl;
  HTML << h1("Execution Times") << endl;
  HTML << start_table("width=\"100%\" border=\"1\" cellspacing=\"0\" cellpadding=\"3\"") << endl;
  HTML << "<tr>" << th("Step") << th("Start") << th("End") << th("Elapsed") << "</tr>" << endl; 
  double total_time_elapsed = 0; 

  for (vector<ExecutionTime>::const_iterator itr = times.begin(); itr != times.end(); itr ++) {  
    const ExecutionTime& t = (*itr);
    
    if (t._message.empty()) continue;

    HTML << "<tr>";
    HTML << td(t._message);
    HTML << td(nonbreaking(t._formatted_time_start));
    HTML << td(nonbreaking(t._formatted_time_end));
    HTML << td(nonbreaking(t._formatted_time_elapsed));
    HTML << "</tr>" << endl; 

    total_time_elapsed += t._time_elapsed;    
  }

  HTML << "<tr class=\"highlight_table_row\">"<< endl;
  HTML << "<td colspan=\"3\" >" << b("Total") << "</td>" << endl;
  HTML << "<td>" << (b(nonbreaking(Settings::elapsedtime2string(total_time_elapsed)))) << "</td>" << endl;
  HTML << "</tr>" << endl;

  HTML << "</table>";
  HTML.close();
}

string breseq_header_string(const Settings& settings)
{
  stringstream ss(ios_base::out | ios_base::app);
  
  //copy over the breseq_graphic which we need if it doesn't exist
  if (!file_exists(settings.breseq_small_graphic_to_file_name.c_str())) {
    _system("cp " + settings.breseq_small_graphic_from_file_name + " " + settings.breseq_small_graphic_to_file_name);
  }
  
  ss << "<table width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"3\">" << endl;
  ss << "<tr>" << endl;
  ss << td(a(settings.website, img(Settings::relative_path(settings.breseq_small_graphic_to_file_name, settings.output_path))));
  ss << endl;
  ss << start_td("width=\"100%\"") << endl;
  ss << settings.byline << endl;
  ss << "<br>";
  ss << a(Settings::relative_path(settings.index_html_file_name, settings.output_path), "mutation predictions"); 
  ss << " | " << endl;
  ss << a(Settings::relative_path(settings.marginal_html_file_name, settings.output_path), "marginal predictions");
  ss << " | " << endl;
  ss << a(Settings::relative_path(settings.summary_html_file_name, settings.output_path), "summary statistics");
  ss << " | " << endl;
  ss << a(Settings::relative_path(settings.final_genome_diff_file_name, settings.output_path), "genome diff");
  ss << " | " << endl;
  ss << a(Settings::relative_path(settings.log_file_name, settings.output_path), "command line log");
  ss << endl;
  ss << "</td></tr></table>" << endl;
  return ss.str();
}


string html_genome_diff_item_table_string(const Settings& settings, genome_diff& gd, diff_entry_list& list_ref)
{
  if(list_ref.empty()) return "";

  diff_entry& first_item = *list_ref.front();
  //mutation
  if(first_item._type.length() == 3)
  {
    diff_entry_list new_list(list_ref.begin(), list_ref.end());
    return Html_Mutation_Table_String(settings, gd, new_list); 
  }
  //evidence
  else
  {
    if(first_item._type == "MC")
    {
      return html_missing_coverage_table_string(list_ref, true, "", "" );
    }
    else if(first_item._type == "RA")
    {
      return html_read_alignment_table_string(list_ref, true);
    }
    else if(first_item._type == "JC")
    {
      return html_new_junction_table_string(list_ref,false, "", "" );
    }
  }  
  return "";
}

/*-----------------------------------------------------------------------------
 *  FORMATTED_MUTATION_ANNOTATION
 *-----------------------------------------------------------------------------*/
string formatted_mutation_annotation(const diff_entry& mut)
{
  stringstream ss;

  // additional formatting for some variables
  if((mut.entry_exists("snp_type")) && (mut["snp_type"] != "intergenic") &&
     (mut["snp_type"] != "noncoding") && (mut["snp_type"] != "pseudogene"))
  {    
    ss << mut["aa_ref_seq"] << mut["aa_position"] << mut["aa_new_seq"];

    string codon_ref_seq = to_underline_red_codon(mut, "codon_ref_seq");
    string codon_new_seq = to_underline_red_codon(mut, "codon_new_seq");
    
    ss << "&nbsp;" << codon_ref_seq << "&rarr;" << codon_new_seq << "&nbsp;";  
  }
  else // mut[SNP_TYPE] == "NC"
  {
    ss << nonbreaking(mut[GENE_POSITION]); 
  }
  return ss.str(); 
}

/*-----------------------------------------------------------------------------
 *  Helper function for formatted_mutation_annotation
 *-----------------------------------------------------------------------------*/
string to_underline_red_codon(const diff_entry& mut, const string& codon_key)
{
  if (!mut.entry_exists(codon_key) || 
      !mut.entry_exists("codon_position") ||
      mut["codon_position"] == "") {
    return "";
  }

  stringstream ss; //!< codon_string
  
  string codon_ref_seq = mut[codon_key];
  for (size_t i = 0; i < codon_ref_seq.size(); i++) {

    if (i == from_string(mut["codon_position"])) {
      ss << font("class=\"mutation_in_codon\"", codon_ref_seq.substr(i,1));
    }
    else 
    {
      ss << codon_ref_seq[i];
    }
  }
  return ss.str();
}

string html_read_alignment_table_string(diff_entry_list& list_ref, bool show_reject_reason, const string& title, const string& relative_link)
{
  stringstream ss; //!< Main build object for function
  stringstream ssf; //!< Used for formatting strings
  
  ss << start_table("border=\"0\" cellspacing=\"1\", cellpadding=\"3\"") << endl;
  
  // Evidence hyperlinks will be the first column if they exist
  bool link;
  if (list_ref.front().get() != 0 && (*list_ref.front()).entry_exists(_EVIDENCE_FILE_NAME)) {
    link = true;
  } else {
    link = false;
  }

  //Determine Number of Columns in Table
  size_t total_cols = link ? 11 : 10;
  
  //Create Column Titles
  //seq_id/position/change/freq/score/cov/annotation/genes/product
  if (title != "") {
    ss << tr(th("colspan=\"" + to_string(total_cols) + 
                "\" align=\"left\" class=\"read_alignment_header_row\"", title)) << endl;
    ss << "<tr>" << endl;
  }
  
  if (link) {
    ss << th("&nbsp;") << endl;
  }
  ss << th("seq&nbsp;id") << endl;
  ss << th("colspan=\"2\"", "position") << endl;
  ss << th("change")     << endl <<
        th("freq")       << endl <<
        th("score")      << endl <<
        th("cov")        << endl <<
        th("annotation") << endl <<
        th("genes")       << endl;
  
  ss << th("width=\"100%\"", "product") << endl;
  ss << "</tr>" << endl;
  
  //Loop through list_ref to build table rows
  for (diff_entry_list::iterator itr = list_ref.begin();
       itr != list_ref.end(); itr ++) {  
    diff_entry& c = **itr;
    bool is_polymorphism = false;
    if (c.entry_exists(FREQUENCY) && (from_string<double>(c[FREQUENCY]) != 1.0)) {
      is_polymorphism = true;
    }
    
    string row_class = "normal_table_row";
    
    if (is_polymorphism) {
      row_class = "polymorphism_table_row";
    }
    
    // Start building a single table row
    ss << start_tr("class=\"" + row_class + "\"");
    if (link) {
      ss << td(a(relative_link + c[_EVIDENCE_FILE_NAME], "*"));
    }

    string fisher_p_value;
    if (c.entry_exists("fisher_strand_p_value") &&
       (c["fisher_strand_p_value"] != "ND")) 
    {
      ssf.precision(1);
      ssf << scientific << from_string<double>(c["fisher_strand_p_value"]); //TODO Confirm
      fisher_p_value = nonbreaking("&nbsp;(" + ssf.str() + ")");
      //Clear Formated String Stream
      ssf.str("");
      ssf.clear();
     }

    ss << td(ALIGN_CENTER, nonbreaking(c[SEQ_ID]));
    ss << td(ALIGN_RIGHT, commify(c["position"]));
    ss << td(ALIGN_RIGHT, c["insert_position"]);
    ss << td(ALIGN_CENTER, c["ref_base"] + "&rarr;" + c["new_base"]); // "Change" Column
    ssf.width(4);
    ssf.precision(1);
    ssf << fixed << from_string<double>(c["frequency"]) * 100 << "%" << endl;
    ss << td(ALIGN_RIGHT, ssf.str());
    //Clear formated string stream
    ssf.str("");
    ssf.clear();

    if (is_polymorphism) {
      // display extra score data for polymorphisms...
      string log_fisher = "999";
      string log_ks = "999";
         
      if (c.entry_exists(FISHER_STRAND_P_VALUE) &&
          from_string<double>(c[FISHER_STRAND_P_VALUE]) > 0 )
        log_fisher = to_string(log(from_string<double>(c[FISHER_STRAND_P_VALUE])));

      if (c.entry_exists(KS_QUALITY_P_VALUE) &&
          from_string<double>(c[KS_QUALITY_P_VALUE]) > 0 ) 
        log_ks = to_string(log(from_string<double>(c[KS_QUALITY_P_VALUE]))/log(10));

      ssf.precision(1);
      ssf << fixed << c["polymorphism_quality"] << " " <<
                      log_fisher << " " <<
                      log_ks;
      ss << td(ALIGN_RIGHT, nonbreaking(ssf.str()));
      //Clear formated string stream
      ssf.str("");
      ssf.clear();

    }

    else {
      ssf.precision(1);
      ssf << fixed << c["quality"] << endl;
      ss << td(ALIGN_RIGHT, nonbreaking(ssf.str()));
      ssf.str("");
      ssf.clear();
    }  
    // Build "cov" column value
    vector<string> temp_cov = split(c[TOT_COV], "/");
    string top_cov = temp_cov[0];
    string bot_cov = temp_cov[1];
    string total_cov = to_string(from_string<uint32_t>(top_cov) + 
                                 from_string<uint32_t>(bot_cov));
    ss << td(ALIGN_CENTER, total_cov);// "Cov" Column
    ss << td(ALIGN_CENTER, nonbreaking(formatted_mutation_annotation(c))); //"Annotation" Column
    ss << td(ALIGN_CENTER, i(nonbreaking(c[GENE_NAME])));
    ss << td(ALIGN_LEFT, htmlize(c[GENE_PRODUCT]));
    ss << "</tr>" << endl;

    if (show_reject_reason) 
    {
      vector<string> reject_reasons = c.get_reject_reasons();
      for (vector<string>::iterator itr = reject_reasons.begin(); itr != reject_reasons.end(); itr ++) 
      {  
        string& reject = (*itr);
       ss << tr("class=\"reject_table_row\"", 
                td("colspan=\"" + to_string(total_cols) + "\"",
                   "Rejected: " + decode_reject_reason(reject)));
      }
    
      /* Fisher Strand Test */
      if (c.entry_exists("fisher_strand_p_value")) 
      {
        ssf.precision(2);
        ssf << scientific << from_string<float>(c["fisher_strand_p_value"]);
        string fisher_strand_p_value = ssf.str();
        
        //Clear formated string stream
        ssf.str("");
        ssf.clear();

        ss << tr("class=\"information_table_row\"", 
                 td("colspan=\"total_cols\"",
                    "Strands of reads supporting (+/-):&nbsp;&nbsp;" +
                    b("new") + " base " + "(" + c[NEW_COV] + ")" + ":&nbsp;&nbsp;" +
                    b("ref") + " base " + "(" + c[REF_COV] + ")" + ":&nbsp;&nbsp;" +
                    b("total") + " (" + c[TOT_COV] + ")")); 
        ss << tr("class=\"information_table_row\"", 
                 td("colspan=\"" + to_string(total_cols) + "\"",
                    "Fisher's exact test strand distribution" +
                    i("p") + "-value = " +fisher_strand_p_value));
      } //end fisher_strand_p_value

      /* Kolmogorov-Smirnov Test */
      if (c.entry_exists("ks_quality_p_value")) {
      ssf.precision(2);
      ssf << scientific << (c.entry_exists(KS_QUALITY_P_VALUE) ? 
        from_string<float>(c["ks_quality_p_value"]) :
        0);
      string ks_quality_p_value = ssf.str();
      
      //Clear formated string stream
      ssf.str("");
      ssf.clear();

      ss << tr("class=\"information_table_row\"", 
               td("colspan=\"" + to_string(total_cols) + "\"",
                  " Kolmogorov-Smirnov test that lower quality scores support polymorphism than reference" + //TODO @ JEB grammar?
                  i("p") + "-value = " +ks_quality_p_value));
      } //end ks_quality_p_value
    } // end show_reject_reason
  } // end list_ref loop

  ss << "</table>" << endl;
  return ss.str();
}

string html_missing_coverage_table_string(diff_entry_list& list_ref, bool show_reject_reason, const string& title, const string& relative_link)
{
  stringstream ss; //!< Main Build Object in Function
  
  ss << endl;
  ss << start_table("border=\"0\" cellspacing=\"1\" cellpadding=\"3\" width=\"100%\"") << endl;
  
  bool coverage_plots = ((list_ref.front()).get() != NULL && (*list_ref.front()).entry_exists("_EVIDENCE_FILE_NAME"));
  
  bool link = ((*list_ref.front()).entry_exists("_side_1_evidence_file_name")) && 
              ((*list_ref.front()).entry_exists("_side_2_evidence_file_name"));
  
  size_t total_cols = link ? 11 : 8;

  if (title != "") {
    ss << "<tr>" << th("colspan=\"" + to_string(total_cols) + "\" align=\"left\" class=\"missing_coverage_header_row\"", title) << "</tr>" << endl;   
  }
  
  ss << "<tr>";
    
  if (link) {
    ss << th("&nbsp;") <<  th("&nbsp;");
    if (coverage_plots) {
      ss << th("&nbsp;");
    }
  }

  ss << th("seq&nbsp;id") << endl <<
        th("start")       << endl <<
        th("end")         << endl <<
        th("size")        << endl <<
        th("&larr;cov")   << endl <<
        th("cov&rarr;")   << endl <<
        th("gene")        << endl;
  
  ss << th("width=\"100%\"", "description") << endl;
  ss << "</tr>" << endl;

      for (diff_entry_list::iterator itr = list_ref.begin(); itr != list_ref.end(); itr ++) {  
        diff_entry& c =  **itr;

        ss << "<tr>" << endl;
      if (link) {
        ss << td(a(relative_link + c[_SIDE_1_EVIDENCE_FILE_NAME], "*")) << endl;
        ss << td(a(relative_link + c[_SIDE_2_EVIDENCE_FILE_NAME], "*")) << endl;
        
        if (coverage_plots  ) {
          ss << td(a(relative_link + c[_EVIDENCE_FILE_NAME], "*")) << endl;
        }

      }

      string start = c[START];
      if (from_string<uint32_t>(c[START_RANGE]) > 0) 
      {
        start += "–" + 
          to_string(from_string<uint32_t>(c[START]) + 
                    from_string<uint32_t>(c[START_RANGE]));
      }
      string end = c[END];
      if (from_string<uint32_t>(c[END_RANGE]) > 0) 
      {
         end += "" + 
           to_string(from_string<uint32_t>(c[END]) -
                     from_string<uint32_t>(c[END_RANGE]));
      }

      string size = to_string(from_string<uint32_t>(c[END]) - from_string<uint32_t>(c[START]) + 1);
        
      if ((from_string<uint32_t>(c[END_RANGE]) > 0) ||
          (from_string<uint32_t>(c[START_RANGE]) > 0)) {
       
        uint32_t size_value = 
        from_string<uint32_t>(c[END]) - 
        from_string<uint32_t>(c[START]) + 1 -
        from_string<uint32_t>(c[END_RANGE]) -
        from_string<uint32_t>(c[START_RANGE]);
       
        size = to_string(size_value) + "–" + size; 
      }
       
      ss << td(nonbreaking(c[SEQ_ID])) << endl;
      ss << td(ALIGN_RIGHT, nonbreaking(start)) << endl;
      ss << td(ALIGN_RIGHT, nonbreaking(end)) << endl;
      ss << td(ALIGN_RIGHT, nonbreaking(size)) << endl;
      ss << td(ALIGN_CENTER, nonbreaking(c[LEFT_OUTSIDE_COV] + "[" + c[LEFT_INSIDE_COV] + "]")) <<endl;
      ss << td(ALIGN_CENTER, nonbreaking(c[RIGHT_INSIDE_COV] + "[" + c[RIGHT_OUTSIDE_COV] + "]")) << endl;
      ss << td(ALIGN_CENTER, i(nonbreaking(c[GENE_NAME]))) << endl;
      ss << td(ALIGN_LEFT, htmlize(c[GENE_PRODUCT])) << endl;
      ss << "</tr>" << endl;
// #     
// #     if ($show_reject_reason)
// #     {
// #       foreach my $reject (GenomeDiff::get_reject_reasons($c))
// #       {
// #         $output_str.= Tr({-class=>'reject_table_row'}, td({-colspan => $total_cols}, "Rejected: " . decode_reject_reason($reject)));
// #       }
// #     }
    if (show_reject_reason && c.entry_exists(REJECT)) {
      vector<string> reject_reasons = c.get_reject_reasons();
      for (vector<string>::iterator itr = reject_reasons.begin();
           itr != reject_reasons.end(); itr ++) {  
        string& reject = (*itr);
        ss << tr("class=\"reject_table_row\"", 
                 td("colspan=\"" + to_string(total_cols) + "\"",
                    "Rejected: " + decode_reject_reason(reject)));  
      }
    }
  }
  ss << "</table>" << endl;
  return ss.str();
}

string html_new_junction_table_string(diff_entry_list& list_ref, bool show_reject_reason, const string& title, const string& relative_link)
{
  stringstream ss; //!<< Main Build Object for Function
  diff_entry& test_item = *list_ref.front();

  bool link = (test_item.entry_exists(_SIDE_1_EVIDENCE_FILE_NAME) &&
               test_item.entry_exists(_SIDE_2_EVIDENCE_FILE_NAME) &&
               test_item.entry_exists(_NEW_JUNCTION_EVIDENCE_FILE_NAME));

  ss << start_table("border=\"0\" cellspacing=\"1\" cellpadding=\"3\"") << endl;
  size_t total_cols = link ? 10 : 8;
  
  if (title != "") {
    ss << tr(th("colspan=\"" + to_string(total_cols) + "\" align=\"left\" class=\"new_junction_header_row\"", title)) << endl;
  }
// #     
// #   #####################
// #   #### HEADER LINE ####
// #   #####################
  ss << "<!-- Header Lines for New Junction -->" << endl; 
  ss << "<tr>" << endl;
    
  if (link) {
    ss << th("colspan=\"2\"", "&nbsp;") << endl;
  }
  ss << th("seq&nbsp;id") << endl <<
        th("position")    << endl <<
        th("overlap")     << endl <<
        th("reads")       << endl <<
        th("score")       << endl <<
        th("annotation")  << endl <<
        th("gene")        << endl;
  
  ss << th("width=\"100%\"","product") << endl;
  ss << "</tr>" << endl;
  ss << endl;
// #   
// #   ####################
// #   #### ITEM LINES ####
// #   ####################
  ss << "<!-- Item Lines for New Junction -->" << endl;
// #   
// #   ## the rows in this table are linked (same background color for every two)
  uint32_t row_bg_color_index = 0; 
   
  for (diff_entry_list::iterator itr = list_ref.begin(); itr != list_ref.end(); itr ++) 
  {  
    diff_entry& c = **itr;
// #     ##############
// #     ### Side 1 ###
// #     ##############
    ss << "<!-- Side 1 Item Lines for New Junction -->" << endl;
// #     my $key = 'side_1';     
// #     my $annotate_key = "junction_" . $c->{"$key\_annotate_key"};
    string key = "side_1";
    string annotate_key = "junction_" + c[key + "_annotation_key"];
// #     $output_str.= start_Tr({-class=> "mutation_table_row_$row_bg_color_index"});
    ss << start_tr("class=\"mutation_table_row_" +
                  to_string(row_bg_color_index) +"\"") << endl;
// #     $output_str.= td({-rowspan=>2}, a({-href=>"$relative_link$c->{_new_junction_evidence_file_name}"}, "*")) if ($link); 
// #     { 
     if (link) {
      ss << td("rowspan=\"2\"", 
              a(relative_link + c[_NEW_JUNCTION_EVIDENCE_FILE_NAME], "*" )) << endl;
     }

     { // Begin hiding data for side 1

// #       $output_str.= td({-rowspan=>1}, a({-href=>"$relative_link$c->{_side_1_evidence_file_name}"}, "?")) if ($link); 
      if (link) {   
        ss << td("rowspan=\"1\"", 
                a(relative_link + c[_SIDE_1_EVIDENCE_FILE_NAME], "?")) << endl;
      }
// #       $output_str.= td({-rowspan=>1, -class=>"$annotate_key"}, nonbreaking($c->{"$key\_seq_id"}));      
      ss << td("rowspan=\"1\" class=\"" + annotate_key + "\"",
            nonbreaking(c[key + "_seq_id"])) << endl;
// #       $output_str.= td( {-align=>"center", -class=>"$annotate_key"}, ($c->{"$key\_strand"} == +1) ? $c->{"$key\_position"} . "&nbsp;=": "=&nbsp;" . $c->{"$key\_position"} );
      //? if (c[key + "_strand"] == "+1") 
      if (from_string<int32_t>(c[key + "_strand"]) == 1) { 
        ss << td("align=\"center\" class-\"" + annotate_key +"\"",
                c[key + "_position"] + "&nbsp;=");
      } else {
        ss << td("align=\"center\" class-\"" + annotate_key +"\"",
                "=&nbsp;" + c[key + "_position"]);
      }
// #       $output_str.= td( {-rowspan=>2, -align=>"center"}, $c->{overlap} );
      ss << td("rowspan=\"2\" align=\"center\"", c["overlap"]) << endl;
// #       $output_str.= td( {-rowspan=>2, -align=>"center"}, $c->{total_reads} );
      ss << td("rowspan=\"2\" align=\"center\"", c["total_reads"]) << endl;
// #       $output_str.= td( {-rowspan=>2, -align=>"center"}, b("&lt;" . $c->{pos_hash_score} . "&gt;") . br . $c->{min_overlap_score} );
      ss << td("rowspan=\"2\" align=\"center\"", 
               b("&lt;" + c["pos_hash_score"] + "&gt;") + 
                 "<br>" + c["min_overlap_score"]) << endl;
// #       $output_str.= td( {-align=>"center", -class=>"$annotate_key"}, nonbreaking($c->{"_$key"}->{gene_position}) );
      ss << td("align=\"center\" class=\"" + annotate_key + "\"", 
              nonbreaking(c["_" + key] + c[GENE_POSITION])) << endl;
// #       $output_str.= td( {-align=>"center", -class=>"$annotate_key"}, i(nonbreaking($c->{"_$key"}->{gene_name})) );
      ss << td("align=\"center\" class=\"" + annotate_key + "\"", 
              nonbreaking(c["_" + key] + c[GENE_NAME])) << endl;
// #       $output_str.= td( {-class=>"$annotate_key"}, htmlize($c->{"_$key"}->{gene_product}) );
      ss << td("class=\"" + annotate_key + "\"",
              htmlize(c["_" + key] + c[GENE_PRODUCT])) << endl;  
// #     }
    } // End hiding data for side 1
// #     $output_str.= end_Tr;
    ss << "</tr>" << endl;


// #     ##############
// #     ### Side 2 ###
// #     ##############
    ss << "<!-- Side 2 Item Lines for New Junction -->" << endl;
// #     $key = 'side_2';
    key = "side_2";
// #     $annotate_key = "junction_" . $c->{"$key\_annotate_key"};
    annotate_key = "junction_" + c[key + "_annotate_key"];
// #     $output_str.= start_Tr({-class=> "mutation_table_row_$row_bg_color_index"});    
    ss << start_tr("class=\"mutation_table_row_" + 
                  to_string(row_bg_color_index) + "\"") << endl;
// #     {
    { //Begin hiding data for side 2
// #       $output_str.= td({-rowspan=>1}, a({-href=>"$relative_link$c->{_side_2_evidence_file_name}"}, "?")) if ($link); 
      if (link) {
        ss << td("rowspan=\"1\"", 
                a("href=\"" + relative_link + 
                  c[_SIDE_2_EVIDENCE_FILE_NAME] + "\"", "?"));
      }
// #       $output_str.= td({-rowspan=>1, -class=>"$annotate_key"}, nonbreaking($c->{"$key\_seq_id"}));    
      ss << td("rowspan=\"1\" class=\"" + annotate_key + "\"",
              nonbreaking(c[key + "_seq_id"])) << endl;
// #       $output_str.= td( {-align=>"center", -class=>"$annotate_key"}, ($c->{"$key\_strand"} == +1) ? $c->{"$key\_position"} . "&nbsp;=": "=&nbsp;" . $c->{"$key\_position"} );
// #
      //? if (c[key + "_strand"] == "+1") 
      if (from_string<int32_t>(c[key + "_strand"]) == 1) { 
        ss << td("align=\"center\" class-\"" + annotate_key +"\"",
                c[key + "_position"] + "&nbsp;=") << endl;
      } else {
        ss << td("align=\"center\" class-\"" + annotate_key +"\"",
                "=&nbsp;" + c[key + "_position"]) << endl;
      } 
// #       $output_str.= td( {-align=>"center", -class=>"$annotate_key"}, nonbreaking($c->{"_$key"}->{gene_position}) );
      ss << td("align=\"center\" class=\"" + annotate_key + "\"",
              nonbreaking(c["_" + key] + c[GENE_POSITION])) << endl;
// #       $output_str.= td( {-align=>"center", -class=>"$annotate_key"}, i(nonbreaking($c->{"_$key"}->{gene_name})) );
      ss << td("align=\"center\" class=\"" + annotate_key + "\"",
              i(nonbreaking(c["_" + key] + c[GENE_NAME]))) << endl;
// #       $output_str.= td( {-class=>"$annotate_key"}, htmlize($c->{"_$key"}->{gene_product}) );
      ss << td("class=\"" + annotate_key + "\"",
              htmlize(c["_" + key] + c[GENE_PRODUCT])) << endl;
// #     }     
    } //end hiding data for side 2
// #     $output_str.= end_Tr;
  ss << "</tr>" << endl;
// #     
// #     if ($show_reject_reason)
// #     {
// #       foreach my $reject (GenomeDiff::get_reject_reasons($c))
// #       {
// #         $output_str.= Tr({-class=>'reject_table_row'}, td({-colspan => $total_cols}, "Rejected: " . decode_reject_reason($reject)));
// #       }
// #     }
  if (show_reject_reason && c.entry_exists("reject")) {
    genome_diff gd;

    vector<string> reject_reasons = c.get_reject_reasons();
    
    for (vector<string>::iterator itr = reject_reasons.begin();
         itr != reject_reasons.end(); itr++) {
      string& reject(*itr);
    
      ss << tr("class=\"reject_table_row\"",
              td("colspan=\"" + to_string(total_cols) + "\"",
                "Rejected: " + decode_reject_reason(reject))) << endl;
    }

  }
// #     
// #     $row_bg_color_index = ($row_bg_color_index+1)%2;
  row_bg_color_index = (row_bg_color_index+1) % 2;//(row_bg_color_index) % 2; 
// #   }
  }// End list_ref Loop
  ss << "</table>" << endl;
  return ss.str();
}

string decode_reject_reason(const string& reject)
{
  
  if (reject == "NJ")
  {
    return "Position hash score below cutoff.";
  }
  else if (reject == "EVALUE")
  {
    return "E-value exceeds prediction threshold.";
  }
  else if (reject == "STRAND")
  {
    return "Prediction not supported by reads on both strands.";
  }
  else if (reject == "FREQ")
  {
    return "Prediction has fr==uency below cutoff threshold.";
  }  
  else if (reject == "COV")
  {
    return "Prediction has coverage below cutoff threshold.";
  }
  else if (reject == "BIAS_P_VALUE")
  {
    return "Prediction has biased strand and/or quality scores supporting polymorphism.";
  }
  else if (reject == "KS_QUALITY_P_VALUE")
  {
    return "Prediction has significantly lower quality scores supporting polymorphism compared to reference.";
  }
  else if (reject == "KS_QUALITY_P_VALUE_UNUSUAL_POLY")
  {
    return "Prediction has biased quality score distribution for polymorphism bases.";
  }
  else if (reject == "KS_QUALITY_P_VALUE_UNUSUAL_REF")
  {
    return "Prediction has biased quality score distribution for new bases.";
  }
  else if (reject == "KS_QUALITY_P_VALUE_UNUSUAL_NEW")
  {
    return "Prediction has biased quality score distribution for ref bases.";
  }
  else if (reject == "KS_QUALITY_P_VALUE_UNUSUAL_ALL")
  {
    return "Prediction has biased quality score distribution for all bases.";
  }
  else if (reject == "FISHER_STRAND_P_VALUE")
  {
    return "Prediction has biased read strand distribution supporting polymorphism.";
  }  
  else if (reject == "POLYMORPHISM_STRAND")
  {
    return "Polymorphism prediction not supported by minimum number of reads on both strands.";
  }
  else if (reject == "POLYMORPHISM_FREQUENCY_CUTOFF")
  {
    return "Polymorphism does not pass frequency cutoff.";
  }
  else if (reject == "HOMOPOLYMER_STRETCH")
  {
    return "Polymorphism is in a homopolymer stretch.";
  }
  
  return "";
}
// # 

/*
 * =====================================================================================
 *        Class:  Evidence_Files
 *  Description:  
 * =====================================================================================
 */
void EvidenceFiles::EvidenceItem::MUT::init_MUT(diff_entry_ptr& item, const Settings& settings, genome_diff& gd)
{ 
  bam_path = settings.reference_bam_file_name;
  bam_path = settings.reference_fasta_file_name;
  prefix = type = item->_type;
  seq_id = item->_id;
  start = from_string<uint32_t>((*item)[START]);
  end = start;

  //Determine What kind of mutation it is.
  if (item->_type == INS) {  
    insert_start = 1;
    insert_end = (*item)[NEW_SEQ].length();
  } else if (item->_type == DEL) {
    genome_diff gd;
    diff_entry_list mutation_list = gd.mutation_evidence_list(*item);
    //Only do deletions if they have within_read evidence
    if (count_if(mutation_list.begin(), mutation_list.end(), diff_entry::is_type("RA")) > 0)
      return;
    else
      end = start + from_string<uint32_t>((*item)[SIZE]) - 1;
  } else if (item->_type == SUB) {
    end = start + (*item)[NEW_SEQ].length() - 1;
  }
  item = item;
  parent_item = item;
  quality_score_cutoff = settings.base_quality_cutoff;
  evidence_list = gd.mutation_evidence_list(*parent_item);

  file_name = this->html_evidence_file_name();
  //Set item's file name
  item->_fields[_EVIDENCE_FILE_NAME] = file_name;
  output_path = settings.evidence_path + "/" + file_name;
}

void EvidenceFiles::EvidenceItem::RA::init_RA(diff_entry_ptr& item, const Settings& settings, genome_diff& gd)
{
  bam_path = settings.reference_bam_file_name;
  bam_path = settings.reference_fasta_file_name;
  seq_id = item->_id;
  start = from_string<uint32_t>((*item)[POSITION]);
  end = from_string<uint32_t>((*item)[POSITION]);
  insert_start = from_string<uint32_t>((*item)[INSERT_POSITION]);
  insert_end = from_string<uint32_t>((*item)[INSERT_POSITION]);
  parent_item = item;
  item = item;
  prefix = item->_type;
  quality_score_cutoff = settings.base_quality_cutoff;
  evidence_list = gd.mutation_evidence_list(*parent_item);
  
  file_name = this->html_evidence_file_name();
  //Set item's file name
  item->_fields[_EVIDENCE_FILE_NAME] = file_name;
  output_path = settings.evidence_path + "/" + file_name;
}


void EvidenceFiles::EvidenceItem::MC::init_MC(diff_entry_ptr& item, const Settings& settings, genome_diff& gd)
{
  //Determine if a parent exists
  diff_entry_ptr parent_item = gd.parent(*item);
  
  if (parent_item.get() == NULL) parent_item = item;

  //Side 1
  side_1_evidence.bam_path = settings.reference_bam_file_name;
  side_1_evidence.fasta_path  = settings.reference_fasta_file_name;
  side_1_evidence.seq_id = item->_id;
  side_1_evidence.start = (from_string<uint32_t>((*item)[START]) - 1);
  side_1_evidence.end = (from_string<uint32_t>((*item)[START]) - 1);
  side_1_evidence.parent_item = parent_item;
  side_1_evidence.item = item;
  side_1_evidence.prefix = "MC_SIDE_1";
  side_1_evidence.file_name = side_1_evidence.html_evidence_file_name();
  side_1_evidence.item->_fields[_SIDE_1_EVIDENCE_FILE_NAME] = side_1_evidence.file_name;
  side_1_evidence.quality_score_cutoff = settings.base_quality_cutoff;
  side_1_evidence.evidence_list = gd.mutation_evidence_list(*parent_item);
  //Side 2
  side_2_evidence.bam_path = settings.reference_bam_file_name;
  side_2_evidence.fasta_path  = settings.reference_fasta_file_name;
  side_2_evidence.seq_id = item->_id;
  side_2_evidence.start = (from_string<uint32_t>((*item)[END]) + 1);
  side_2_evidence.end = (from_string<uint32_t>((*item)[END]) + 1);
  side_2_evidence.parent_item = parent_item;
  side_2_evidence.item = item;
  side_2_evidence.prefix = "MC_SIDE_2";
  side_2_evidence.file_name = side_2_evidence.html_evidence_file_name();
  side_2_evidence.item->_fields[_SIDE_2_EVIDENCE_FILE_NAME] = side_2_evidence.file_name;
  side_2_evidence.quality_score_cutoff = settings.base_quality_cutoff;
  side_2_evidence.evidence_list = gd.mutation_evidence_list(*parent_item);
  //Evidence_Files
  evidence.seq_id = item->_id;
  evidence.start = from_string<uint32_t>((*item)[START]);
  evidence.end = from_string<uint32_t>((*item)[END]);
  evidence.parent_item = parent_item;
  evidence.item = item;
  evidence.prefix = "MC_PLOT";
  evidence.plot = (*item)[_COVERAGE_PLOT_FILE_NAME];
  evidence.file_name = evidence.html_evidence_file_name(); 
  evidence.item->_fields[_COVERAGE_PLOT_FILE_NAME];
  evidence.evidence_list = gd.mutation_evidence_list(*parent_item);
}



void EvidenceFiles::EvidenceItem::JC::init_JC(diff_entry_ptr& item, const Settings& settings, genome_diff& gd)
{
  diff_entry_ptr parent_item = gd.parent(*item);
  if (parent_item.get() == NULL) parent_item = item;

  //Define here for cleaner declarations below.
  uint32_t start = UINT_MAX;
  uint32_t end = UINT_MAX;
  //Regenereate the alignment overlap from the junction key
  uint32_t alignment_overlap = from_string<uint32_t>((*item)[ALIGNMENT_OVERLAP]);
  uint32_t flanking_left = from_string<uint32_t>((*item)[FLANKING_LEFT]);
   
  if ( alignment_overlap == 0) {
    start = flanking_left;
    end = (flanking_left + 1);
  } else if (alignment_overlap > 0) {
    start = (flanking_left + 1);
    end = (flanking_left + alignment_overlap);
  } else { //item->{overlap} < 0);
    start = (flanking_left + 1);
    end = (flanking_left - alignment_overlap);
  }
  assert(start != UINT_MAX || end != UINT_MAX);//TODO @JEB possible negative value?
  
  //Define here for cleaner declarations below.
  uint32_t side_1_overlap = from_string<uint32_t>((*item)[SIDE_1_OVERLAP]);
  uint32_t side_1_position = from_string<uint32_t>((*item)[SIDE_1_POSITION]);
  string side_1_strand = (*item)[SIDE_1_STRAND];
  string side_1_seq_id = (*item)[SIDE_1_SEQ_ID];
  uint32_t side_2_overlap = from_string<uint32_t>((*item)[SIDE_2_OVERLAP]);
  uint32_t side_2_position = from_string<uint32_t>((*item)[SIDE_2_POSITION]);
  string side_2_strand = (*item)[SIDE_2_STRAND];
  string side_2_seq_id = (*item)[SIDE_2_SEQ_ID];
  uint32_t truncate_start = (flanking_left + 1 + abs(alignment_overlap - side_2_overlap));

  //Evidence 
  evidence.bam_path = settings.junction_bam_file_name;
  evidence.fasta_path = settings.candidate_junction_fasta_file_name;
  evidence.seq_id = item->_id;
  evidence.start = start;
  evidence.end = end;
  evidence.parent_item = parent_item;
  evidence.item = item;
  evidence.prefix = "JC";
  evidence.file_name = evidence.html_evidence_file_name();
  evidence.item->_fields[_NEW_JUNCTION_EVIDENCE_FILE_NAME] = evidence.file_name;
  evidence.quality_score_cutoff = settings.base_quality_cutoff;
  evidence.evidence_list = gd.mutation_evidence_list(*parent_item);
  
    //Extra information
    evidence.alignment_empty_change_line = true;
    //Reference info
    evidence.alignment_reference_info_side_1.truncate_end = flanking_left + side_1_overlap;
    evidence.alignment_reference_info_side_1.ghost_end = side_1_position;
    evidence.alignment_reference_info_side_1.ghost_strand = side_1_seq_id;
    evidence.alignment_reference_info_side_1.ghost_seq_id = side_1_seq_id;
  
    evidence.alignment_reference_info_side_2.truncate_start = flanking_left + side_1_overlap;
    evidence.alignment_reference_info_side_2.ghost_start = side_2_position;
    evidence.alignment_reference_info_side_2.ghost_strand = side_2_seq_id;
    evidence.alignment_reference_info_side_2.ghost_seq_id = side_2_seq_id;
  
  //This is the mothership file taht we show first when clicking on evidence from a mutation...
  item->_fields[_EVIDENCE_FILE_NAME] = item->_fields[_NEW_JUNCTION_EVIDENCE_FILE_NAME];


  //Side 1
  side_1_evidence.bam_path = settings.junction_bam_file_name;
  side_1_evidence.fasta_path = settings.candidate_junction_fasta_file_name;
  side_1_evidence.seq_id = (*item)[SIDE_1_SEQ_ID];
  side_1_evidence.start = from_string<uint32_t>((*item)[SIDE_1_POSITION]);
  side_1_evidence.end = from_string<uint32_t>((*item)[SIDE_1_POSITION]);
  side_1_evidence.parent_item = parent_item;
  side_1_evidence.item = item;
  side_1_evidence.prefix = "JC_SIDE_1_" + (*item)[SIDE_2_SEQ_ID] + "_" + (*item)[SIDE_2_POSITION]; // Need to be unique
  side_1_evidence.file_name = side_1_evidence.html_evidence_file_name();
  side_1_evidence.item->_fields[_SIDE_1_EVIDENCE_FILE_NAME] = side_1_evidence.file_name;
  side_1_evidence.quality_score_cutoff = settings.base_quality_cutoff;
  side_1_evidence.evidence_list = gd.mutation_evidence_list(*parent_item);
  //Side 2
  side_2_evidence.bam_path = settings.junction_bam_file_name;
  side_2_evidence.fasta_path = settings.candidate_junction_fasta_file_name;
  side_2_evidence.seq_id = (*item)[SIDE_2_SEQ_ID];
  side_2_evidence.start = from_string<uint32_t>((*item)[SIDE_2_POSITION]);
  side_2_evidence.end = from_string<uint32_t>((*item)[SIDE_2_POSITION]);
  side_2_evidence.parent_item = parent_item;
  side_2_evidence.item = item;
  side_2_evidence.prefix = ("JC_SIDE_2_" + (*item)[SIDE_1_POSITION] + "_" + (*item)[SIDE_1_POSITION]); // Need to be unique
  side_2_evidence.file_name = side_2_evidence.html_evidence_file_name();
  side_2_evidence.item->_fields[_SIDE_2_EVIDENCE_FILE_NAME] = side_2_evidence.file_name;
  side_2_evidence.quality_score_cutoff = settings.base_quality_cutoff;
  side_2_evidence.evidence_list = gd.mutation_evidence_list(*parent_item);

}

void EvidenceFiles::EvidenceItem::MC::Evidence::create_html_file(const Settings& settings, genome_diff& gd)
{
 if (true) { //TODO @JEB settings.verbose
   cerr << "Creating evidence file: " << this->file_name << endl;   
 }
  ofstream HTML(this->output_path.c_str());

 if (!HTML.good()) {
    cerr << "Could not open file: " << this->output_path << endl;
    assert(HTML.good());
  }
  
  // Build HTML Head
  HTML << html_header("BRESEQ :: Results", settings);
  
  vector<string> types = make_list<string>("RA")("MC")("JC");
 
  for (vector<string>::iterator itr = types.begin(); itr != types.end(); itr ++) 
  {  
    string& type = (*itr);

    this->evidence_list.remove_if(diff_entry::is_not_type(type));   

    if(this->evidence_list.empty()) continue;

    HTML << html_genome_diff_item_table_string(settings, gd, this->evidence_list);
    HTML << "<p>";
  }
  //Create MC Evidence Plot
  HTML << div(ALIGN_CENTER, img(this->plot));

  HTML << "</html>";
  HTML.close();
}


void EvidenceFiles::EvidenceItem::BaseEvidence::create_html_file(const Settings&settings, genome_diff& gd)
{
 if (true) { //TODO @JEB settings.verbose
   cerr << "Creating evidence file: " << this->file_name << endl;   
 }

 ofstream HTML(this->output_path.c_str());

 if (!HTML.good()) {
    cerr << "Could not open file: " << this->output_path << endl;
    assert(HTML.good());
  }
  
  // Build HTML Head
  HTML << html_header("BRESEQ :: Results", settings);
  
  vector<string> types = make_list<string>("RA")("MC")("JC");
 
  for (vector<string>::iterator itr = types.begin(); itr != types.end(); itr ++) 
  {  
    string& type = (*itr);

    this->evidence_list.remove_if(diff_entry::is_not_type(type));   

    if(this->evidence_list.empty()) continue;

    HTML << html_genome_diff_item_table_string(settings, gd, this->evidence_list);
    HTML << "<p>";
  }
  stringstream ss;   
  ss << this->seq_id << ":" << this->start;
  ss << (this->insert_start == UINT_MAX ? "" :to_string(this->insert_start));
  ss << "-" << this->end;
  ss << (this->insert_end == UINT_MAX ? "" : to_string(this->insert_end));
  cerr << "Creating read alignment for region " << ss.str() << endl;
  
  //TODO @GRC settings.maximum_reads_to_align 
  alignment_output ao(this->bam_path, this->fasta_path, 200 , 
    settings.base_quality_cutoff);
   
  HTML << ao.html_alignment(ss.str());

  

  HTML << endl << "</html>";
  HTML.close();

}


string EvidenceFiles::EvidenceItem::BaseEvidence::html_evidence_file_name()
{
  stringstream ss;

  ss << prefix;
  ss << "_" << seq_id;
  ss << "_" << start;
  if (this->insert_start != UINT_MAX) 
    ss << "." << this->insert_start;
  ss << "_" << end;
  if (this->insert_end != UINT_MAX)
    ss << "." << this->insert_end;
  ss << "_alignment.html";
  
  return ss.str();
}
  

void EvidenceFiles::initEvidenceItems(const Settings& settings, genome_diff& gd)
{  
  // Fasta and BAM files for making alignments.
  string reference_bam_file_name = settings.reference_bam_file_name;
  string reference_fasta_file_name = settings.reference_fasta_file_name;

  // hybrids use different BAM files for making the alignments!!!
  string junction_bam_file_name = settings.junction_bam_file_name;
  string junction_fasta_file_name = settings.candidate_junction_fasta_file_name;

  create_path(settings.evidence_path);

  // Handle MC 
  // We make alignments of two regions for deletions: upstream and downstream edges.
  diff_entry_list MC_items = gd.list(make_list<string>(MC));
  MC_items.remove_if(diff_entry::field_exists(NO_SHOW));

  diff_entry_list::iterator itr_item; 
  for (itr_item = MC_items.begin(); itr_item != MC_items.end(); itr_item++) {
    counted_ptr<EvidenceItem::MC> mc;

    mc->init_MC(*itr_item, settings, gd);
    
    this->m_MC_items.push_back(mc);
  }
  
  //Handle Mutations 
  diff_entry_list MUT_items = gd.list(make_list<string>(SNP)(INS)(DEL)(SUB));
  MUT_items.remove_if(diff_entry::field_exists(NO_SHOW));

  for (itr_item = MUT_items.begin(); itr_item != MUT_items.end(); itr_item++) {
    
    //Build mutation and add to 
    counted_ptr<EvidenceItem::MUT> mut;
    
    mut->init_MUT(*itr_item, settings, gd);
    
    this->m_MUT_items.push_back(mut);
    
    //Add evidence to RA items as well
    diff_entry_list MUT_evidence = gd.mutation_evidence_list(**itr_item);
    for (diff_entry_list::iterator evidence = MUT_evidence.begin();
         evidence != MUT_evidence.end(); evidence++) {
      if ((*evidence)->_type != RA)
        continue;
      else
        (**evidence)[_EVIDENCE_FILE_NAME] = (**itr_item)[_EVIDENCE_FILE_NAME];

    }
  }

  // Handle RA 
  
  diff_entry_list RA_items = gd.filter_used_as_evidence(gd.list(make_list<string>(RA)));
  RA_items.remove_if(diff_entry::field_exists(NO_SHOW));
  for (itr_item = RA_items.begin(); itr_item != RA_items.end(); itr_item++) {
    counted_ptr<EvidenceItem::RA> ra;

    ra->init_RA(*itr_item, settings, gd);

    this->m_RA_items.push_back(ra);
  }
  
  // This additional information is used for the complex reference line.
  // Note that it is completely determined by the original candidate junction sequence 
  // positions and overlap: alignment_pos and alignment_overlap.
  
  diff_entry_list JC_items = gd.list(make_list<string>(JC));
  JC_items.remove_if(diff_entry::field_exists(NO_SHOW));

  for (itr_item = JC_items.begin(); itr_item != JC_items.end(); itr_item ++) 
  {  
    counted_ptr<EvidenceItem::JC> jc;
    jc->init_JC(*itr_item, settings, gd);

    this->m_JC_items.push_back(jc);
  }

}

/*-----------------------------------------------------------------------------
 *  Create the HTML Evidence File
 *-----------------------------------------------------------------------------*/


void EvidenceFiles::htmlOutput(const Settings& settings, genome_diff& gd)
{
  this->initEvidenceItems(settings, gd);
  

  //MC Evidence
  for (list<counted_ptr<EvidenceItem::MC> >::iterator itr = m_MC_items.begin();
       itr != m_MC_items.end(); itr++) {
    counted_ptr<EvidenceItem::MC> mc = *itr;
    
    mc->evidence.create_html_file(settings, gd);
    mc->side_1_evidence.create_html_file(settings, gd);
    mc->side_2_evidence.create_html_file(settings, gd);

  }

  //RA Evidence
  for (list<counted_ptr<EvidenceItem::RA> >::iterator itr = m_RA_items.begin();
       itr != m_RA_items.end(); itr++) {
    counted_ptr<EvidenceItem::RA> ra = *itr;

    ra->create_html_file(settings, gd);

  }

  //JC Evidence 
  for (list<counted_ptr<EvidenceItem::JC> >::iterator itr = m_JC_items.begin();
       itr != m_JC_items.end(); itr++) {
    counted_ptr<EvidenceItem::JC> jc = *itr;

    jc->evidence.create_html_file(settings, gd);
    jc->side_1_evidence.create_html_file(settings, gd);
    jc->side_2_evidence.create_html_file(settings, gd);

  }

  //MUT Evidence
  for (list<counted_ptr<EvidenceItem::MUT> >::iterator itr = m_MUT_items.begin();
       itr != m_MUT_items.end(); itr ++) {
    counted_ptr<EvidenceItem::MUT> mut = *itr;

    mut->create_html_file(settings, gd);
  }

}

/*-----------------------------------------------------------------------------
 *  //End Create_Evidence_Files
 *-----------------------------------------------------------------------------*/

void draw_coverage(Settings& settings, cReferenceSequences& ref_seq_info, genome_diff& gd)
{  
  coverage_output co(
                     settings.reference_bam_file_name, 
                     settings.reference_fasta_file_name, 
                     settings.coverage_plot_r_script_file_name, 
                     settings.coverage_plot_path
                     );
  co.output_format("png");
  
  create_path(settings.coverage_plot_path);
  string coverage_plot_path = settings.coverage_plot_path;
  
  // Coverage overview plots of entire reference sequences
  for (cReferenceSequences::iterator it = ref_seq_info.begin(); it != ref_seq_info.end(); ++it)
  {
    cAnnotatedSequence& seq = *it;
    string region = seq.m_seq_id + ":" + "1" + "-" + to_string(seq.m_length);
    string this_complete_coverage_text_file_name = settings.file_name(settings.overview_coverage_plot_file_name, "@", seq.m_seq_id);
    
    cerr << "Creating coverage plot for region: " << region << endl;
    co.plot(region, this_complete_coverage_text_file_name);
   }
  
  // Zoom-in plots of individual deletions
  vector<string> mc_types = make_list<string>("MC");
	diff_entry_list mc = gd.list(mc_types);
  for (diff_entry_list::iterator it=mc.begin(); it!=mc.end(); it++)
  {
    diff_entry_ptr& item = *it;
    uint32_t start = from_string<uint32_t>((*item)[START]);
    uint32_t end = from_string<uint32_t>((*item)[END]);
    uint32_t size = end - start + 1;
    
    uint32_t _shaded_flanking = floor(static_cast<double>(size) / 10.0);
    if (_shaded_flanking < 100) _shaded_flanking = 100;
    co.shaded_flanking(_shaded_flanking);
    
    string region = (*item)[SEQ_ID] + ":" + (*item)[START] + "-" + (*item)[END];
    string coverage_plot_file_name = settings.evidence_path + "/" + (*item)[SEQ_ID] + "_" + (*item)[START] + "-" + (*item)[END];

    cerr << "Creating coverage plot for region: " << region << endl;
    co.plot(region, coverage_plot_file_name);
  }
}

/*
 * =====================================================================================
 *        Class:  Html_Mutation_Table_String
 *  Description:  
 * =====================================================================================
 */
Html_Mutation_Table_String::Html_Mutation_Table_String(
                                                       const Settings& settings,
                                                       genome_diff& gd,
                                                       diff_entry_list& list_ref,
                                                       vector<string>& gd_name_list_ref,
                                                       Options& options,
                                                       bool legend_row, 
                                                       bool one_ref_seq,
                                                       const string& relative_link
                                                       )
  : string()
  , total_cols(0)
  , settings(settings)
  , gd(gd)
  , list_ref(list_ref)
  , legend_row(legend_row)
  , one_ref_seq(one_ref_seq)
  , gd_name_list_ref(gd_name_list_ref)
  , options(options)
  , relative_link(relative_link)
{
  (*this) += "<!--Output Html_Mutation_Table_String-->\n";
  (*this) += "<table border=\"0\" cellspacing=\"1\" cellpadding=\"3\">\n";
  
  this->Header_Line();
  this->Item_Lines();
}

Html_Mutation_Table_String::Html_Mutation_Table_String(
                                                       const Settings& settings,
                                                       genome_diff& gd,
                                                       diff_entry_list& list_ref,
  			                                               const string& relative_path, 
                                                       bool legend_row, 
                                                       bool one_ref_seq
                                                       )
  : string()
  , total_cols(0)
  , settings(settings)
  , gd(gd)
  , list_ref(list_ref)
  , legend_row(legend_row)
  , one_ref_seq(one_ref_seq)
  , relative_link(relative_path)

{


  
  (*this) += "<!--Output Html_Mutation_Table_String-->\n";
  (*this) += "<table border=\"0\" cellspacing=\"1\" cellpadding=\"3\">\n";
  vector<string> gd_name_list_ref;
  this->gd_name_list_ref = gd_name_list_ref;
  
  Options options;
  options.repeat_header = false;
  this->options = options;
  
  this->Header_Line();
  this->Item_Lines();
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Html_Mutation_Table_String
 *      Method:  Html_Mutation_Table_String :: Header_Line
 * Description:  
 *--------------------------------------------------------------------------------------
 */
void Html_Mutation_Table_String::Header_Line()
{
// #   #####################
// #   #### HEADER LINE ####
// #   #####################
  stringstream ss(ios_base::out | ios_base::app); //<! Main Build Object for Function
  string header_text = ((list_ref.size() > 1) ? "Predicted mutations" : "Predicted mutation");
// #
// #   # There are three possibilities for the frequency column(s)
// #   # (1) We don't want it at all. (Single genome no poly prediction)   
// #   my @freq_header_list = ();
  vector<string> freq_header_list = make_list<string>("");
// #   # (2) We want multiple columns because we are comparing genomes.
// #   if (defined $gd_name_list_ref)
// #   {
// #     @freq_header_list = @$gd_name_list_ref;
// #   }

// #   # (3) We want a single column (polymorphism prediction)
// #   elsif ($settings->{polymorphism_prediction})
// #   {
// #     @freq_header_list = ("freq");
// #   }
  if (gd_name_list_ref.size() > 0) {
    freq_header_list = gd_name_list_ref;
  } 
  else if(settings.polymorphism_prediction) {
    freq_header_list = make_list<string>("freq");
  }

  if (settings.lenski_format) {
    vector<string> header_list = split(freq_header_list.front(), "|");
    size_t header_rows = header_list.size() - 1; //!< -1 is necessary for C++

    total_cols = 7 + freq_header_list.size() ;
    if(!one_ref_seq) total_cols += 1; 
    if(!settings.no_evidence) total_cols += 1;

    for (size_t i = 1; i <= header_rows; i++) {
     ss << "<!-- Header Line -->" << endl;
     ss << "<tr>" << endl;
      if(!settings.no_evidence)
        ss << th("evidence") << endl;
      if(!one_ref_seq)
       ss << th(nonbreaking("seq id")) << endl;

      ss << th( (header_rows == i) ? "position" : "") << endl;
      ss << th( (header_rows == i) ? "mutation" : "") << endl;
      ss << th( (header_rows == i) ? "annotation" : "") << endl;
      ss << th( (header_rows == i) ? "gene" : "") << endl;
      ss << th("width=\"100%\"", (header_rows == i) ? "description" : "") << endl;
      for (vector<string>::iterator itr = freq_header_list.begin(); itr != freq_header_list.end(); itr++) {
        string& freq_header_item(*itr);        

        vector<string> header_list = split(freq_header_item, "|");        
        string this_header_string = header_list[i-1];
        while(this_header_string.find("_") != string::npos)
        {
          size_t pos = this_header_string.find("_");
          this_header_string.replace(pos, 1, "&nbsp;");//TODO confim "_" is 1 char
        }
        string this_header_string_1 = header_list.front();
        string this_header_string_2 = header_list[1];

        string color = "black";  
  
        if( this_header_string_1 == "UC" )
        color = "gray";
        else if( this_header_string_1 == "clade_1" )
        color = "green";  
        else if( this_header_string_1 == "clade_2" )
          color = "blue";  
        else if( this_header_string_1 == "clade_3" &&
                 this_header_string_2 == "ZDB483" ||
                 this_header_string_2 == "ZDB30" )
          color = "orange";
        else if( this_header_string_1 == "clade_3" )
          color = "red";

        ss << th("style=\"background-color:" + color + "\"", this_header_string) << endl;  
        ss << th(this_header_string) << endl;
      }    
      ss << th(header_rows == i ? "position" : "") << endl;
      ss << "</tr>" << endl;
    }
  } else {
    total_cols = 5 + freq_header_list.size() - 1;
    if (!one_ref_seq) {
    total_cols += 1;
    }
    if (!settings.no_evidence) {
      total_cols += 1;
    }

    ss <<  "<tr>" << endl;
   if (!settings.no_evidence) {
     ss << th("evidence") << endl;
   } 
   if(!one_ref_seq) {
     ss << th(nonbreaking("seq id")) << endl;
   }

   ss << th("position") << endl;
   ss << th("mutation") << endl;

   if(!freq_header_list[0].empty()) {
     for (vector<string>::iterator itr = freq_header_list.begin() ;
          itr != freq_header_list.end() ; itr++) {
       string& freq_header_item = *itr;
       ss << th(freq_header_item) << endl;  
     }
   }
 
   ss << th("annotation") << endl;
   ss << th("gene") << endl;
   ss << th("width=\"100%\"","description") << endl;
   ss << "</tr>" << endl; 
  }

  if(!settings.no_header) {
          (*this) += tr(th("colspan=\"" + to_string(total_cols) +
                "\" align=\"left\" class=\"mutation_header_row\"", header_text));
  }

 ss << endl;
 (*this) += ss.str(); 
}
//===============================================================================
//       CLASS: Html_Mutation_Table_String
//      METHOD: Item_Lines
// DESCRIPTION: 
//===============================================================================
void Html_Mutation_Table_String::Item_Lines()
{
// #   ####################
// #   #### ITEM LINES ####
// #   ####################
// # 
  size_t row_num = 0;

  stringstream ss(ios_base::out | ios_base::app); 
  ss << "<!-- Item Lines -->" << endl;
  for (diff_entry_list::iterator itr = list_ref.begin();
       itr != list_ref.end(); itr ++) { 
    diff_entry& mut = (**itr);
    if ((row_num != 0) && (options.repeat_header != 0) && (row_num % options.repeat_header == 0))
    {
      Header_Line();
    }
    row_num++;
        
    // Build Evidence Column
    string evidence_string;
    if (!settings.no_evidence) {
      bool already_added_RA = false;
       
      diff_entry_list mutation_evidence_list = gd.mutation_evidence_list(mut);
      
      for (diff_entry_list::iterator itr = mutation_evidence_list.begin();
           itr != mutation_evidence_list.end(); itr ++) {  
        diff_entry& evidence_item = **itr;

        if (evidence_item._type == RA) {
          if (already_added_RA) 
            continue;
          else 
            already_added_RA = true;
        }
        
        if (!evidence_string.empty()) evidence_string += "&nbsp;";
        evidence_string += a(relative_link + evidence_item[_EVIDENCE_FILE_NAME], evidence_item._type);
      }
    }
  
    string row_class = "normal_table_row";
// #     
// #     # There are three possibilities for the frequency column(s)
// #     # (1) We don't want it at all. (Single genome no poly prediction)   
// #     my @freq_list = ();
    vector<string> freq_list = make_list<string>("");
  //TODO @JEB 
// #     # (2) We want multiple columns because we are comparing genomes.
// #     if (defined $gd_name_list_ref) {
// #       #"_freq_[name]" keys were made within GenomeDiff structure        
// #       @freq_list = map { $mut->{"frequency_$_"} } @$gd_name_list_ref; 
// #       $row_class = "alternate_table_row_" . ($row_num % 2);
// #     }
   
// #     # (3) We want a single column (polymorphism prediction)
// #     elsif ($settings->{polymorphism_prediction}) {      
// #       if ((defined $mut->{frequency}) && ($mut->{frequency} != 1)) {
// #         $row_class = "polymorphism_table_row";  
// #       }     
// #       push @freq_list, $mut->{frequency};
// #     }
// #     
    if (settings.polymorphism_prediction) {
      if(mut.entry_exists("frequency") && from_string<size_t>(mut["frequency"]) != 1) {
        string row_class = "polymorphism_table_row";
      }
      freq_list.push_back(mut["frequency"]);
    }
      
    // ### marshal cells defined depending on mutation type
    string cell_seq_id = nonbreaking(mut[SEQ_ID]);
    string cell_position = commify(mut[POSITION]);
    string cell_mutation;
    string cell_mutation_annotation = nonbreaking(formatted_mutation_annotation(mut));
    string cell_gene_name = i(nonbreaking(mut[GENE_NAME]));
    string cell_gene_product = htmlize(mut["gene_product"]);

    if (mut._type == SNP) {
      cell_mutation = mut["_ref_seq"] + "&rarr;" + mut[NEW_SEQ];
    } else if (mut._type == INS) {
      cell_mutation = "+";
      cell_mutation += mut[NEW_SEQ];
    } else if (mut._type == DEL) {
      cell_mutation = nonbreaking("&Delta;") + commify(mut["size"]) + " bp";
      string annotation_str;
      annotation_str = mut.entry_exists("mediated") ? mut["mediated"] + "-mediated"  : ""; 
      if(annotation_str.empty()) {
        annotation_str = nonbreaking(mut["gene_position"]);
      } 
      cell_mutation_annotation =  nonbreaking(annotation_str);
    } else if (mut._type == SUB) {
      cell_mutation = nonbreaking(mut["size"] + "bp&rarr;" + mut["new_seq"]);
    } else if (mut._type == CON) {
      cell_mutation = nonbreaking(mut["size"] + "bp&rarr;" + mut["region"]);
    } else if (mut._type == MOB) {
      stringstream s(ios_base::out | ios_base::app);
      
      stringstream s_start(ios_base::out | ios_base::app);
      if (mut.entry_exists("ins_start")) {
        s_start << "+" << mut["ins_start"];
      }
      if (mut.entry_exists("del_start")) {
        s_start << "&Delta;" << mut["del_start"];
      }
      if (!(s_start.str()).empty()) {
        s << s_start << " :: ";
      }

      s << mut["repeat_name"] << " (";
      s << (mut["strand"] == "+1" ? "+" : (mut["strand"] == "-1" ? "&minus;" : "?"));
      s << ")";

      stringstream s_end(ios_base::out | ios_base::app);
      if (mut.entry_exists("del_end")) {
        s_end << "&Delta;" + mut["del_end"];
      }
      if (mut.entry_exists("ins_end")) {
        s_end << "+" << mut["ins_end"];
      }
      if (!(s_end.str()).empty()) {
        s << " :: " << s_end;
      }
      
      //dup_str not necessary
      stringstream s_dup(ios_base::out | ios_base::app);
      if (from_string<int>(mut["duplication_size"]) >= 0) {
        s_dup << "+" << mut["duplication_size"];
      } else {
        s_dup << "&Delta;" << abs(from_string<int>(mut["duplication_size"]));
      }
      s << s_dup.str() << "bp";
      cell_mutation = nonbreaking(s.str());
    } else if (mut._type == INV) {
      cell_mutation = nonbreaking(commify(mut["size"]) + " bp inversion");
      cell_gene_name = i(nonbreaking(mut["gene_name_1"])) + "&darr;" +
                       i(nonbreaking(mut["gene_name_2"]));
      cell_gene_product = htmlize(mut["gene_product_1"]) + "&darr;" + 
                          htmlize(mut["gene_product_2"]);
    } else if (mut._type == AMP) {
      cell_mutation = nonbreaking(commify(mut["size"]) + "bp x " + mut["new_copy_number"]);
      cell_mutation_annotation = 
        from_string<uint8_t>(mut["new_copy_number"]) == 2 ?
          "duplication" : "amplification";
    }
    // ###### PRINT THE TABLE ROW ####
    ss << endl << "<!-- Print The Table Row -->" << endl; 
    ss << start_tr("class=\"" + row_class + "\"") << endl;

    if (!settings.no_evidence) {
      ss << td(ALIGN_CENTER, evidence_string) << "<!-- Evidence -->" << endl;
    }
    if (!one_ref_seq) {
      ss << td(ALIGN_CENTER, cell_seq_id) << "<!-- Seq_Id -->" << endl;
    }
    ss << td(ALIGN_CENTER, cell_position) << "<!-- Position -->" << endl;

    ss << td(ALIGN_CENTER, cell_mutation) << "<!-- Cell Mutation -->" << endl;
    if (settings.lenski_format) {
      ss << "<!-- Lenski_Format -->" << endl;
      ss << td(ALIGN_CENTER, cell_mutation_annotation) << endl;
      ss << td(ALIGN_CENTER, cell_gene_name) << endl;
      ss << td(ALIGN_CENTER, cell_gene_product) << endl;
    }
    
    //Need if statement for C++
    if (freq_list.size() >= 1 && !freq_list[0].empty()) {
      ss << freq_cols(freq_list) << endl;
    } 
    if (settings.lenski_format) {
      ss << "<!-- Lenski Format -->" << endl;
      ss << td(ALIGN_CENTER, cell_position) << endl;
    } else {
      ss << td(ALIGN_CENTER, cell_mutation_annotation) << endl;
      ss << td(ALIGN_CENTER, cell_gene_name) << endl;
      ss << td(ALIGN_CENTER, cell_gene_product) << endl;
    }
    ss << "</tr>" << endl;
    
    ss << "<!-- End Table Row -->" << endl;
  } //##### END TABLE ROW ####
  
  if (legend_row) {
    ss << "<tr>" << endl;
    ss << td("colspan=\"" + to_string(total_cols) + "\"", 
                    b("Evidence codes: RA = read alignment, MC = missing coverage, JC = new junction"));
    ss << "</tr>" << endl;
  }
  ss << "</table>";
  
  (*this) += ss.str();
}
// # 
//===============================================================================
//       CLASS: Html_Mutation_Table_String
//      METHOD: freq_to_string  
// DESCRIPTION: Helper function used in Item_Lines
//===============================================================================
string Html_Mutation_Table_String::freq_to_string(const string& freq)
{
  if (freq == "H") {
    return "H";
  }
  if (from_string<double>(freq) == 0.0) {
    return "";
  }
  stringstream ss;
  if (from_string<double>(freq) == 1.0 || freq.empty()) {
    ss << "100%";
  } 
  else {
    float conv_freq = from_string<float>(freq) * 100;
    ss.width(4);
    ss.precision(1);
    ss << conv_freq;
  }
  return ss.str(); 
}

//===============================================================================
//       CLASS: Html_Mutation_Table_String
//      METHOD: freq_cols  
// DESCRIPTION: Helper function used in Item_Lines
//===============================================================================
string Html_Mutation_Table_String::freq_cols(vector<string> freq_list)
{
  stringstream ss(ios_base::out | ios_base::app);
  for (vector<string>::iterator itr = freq_list.begin();
       itr != freq_list.end(); itr ++) {  
    string& freq = (*itr);
    if (settings.shade_frequencies) {
      string bgcolor;
      if (freq == "1") {
       bgcolor = "Blue";
      }
      if (!bgcolor.empty()) {
        ss << td("align=\"right\" bgcolor=\"" + bgcolor +"\"", "&nbsp;"); //TODO Check
      } 
      else {
        ss << td(ALIGN_RIGHT,"&nbsp;");
      }
    }
    else {
      ss << td (ALIGN_RIGHT, freq_to_string(freq));
    }
  }
  return ss.str();
}


}//end namespace output
}//end namespace breseq

