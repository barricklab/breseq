#include "breseq/output.h"
#include "breseq/anyoption.h"
#include "breseq/alignment_output.h"


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
  /* substitute nonbreaking hyphen */
  substitute(retval, "-", "&#8209;");

  /* substitute nonbreaking en dash */
  substitute(retval, "–", "&#8211;");

  /* substitute nonbreaking space */
  substitute(retval, " ", "&nbsp;");

  return retval;
}
/*-----------------------------------------------------------------------------
 *  HTML Utility for Encoding HTML
 *-----------------------------------------------------------------------------*/
string htmlize (const string& input) 
{
  string retval = input;
  //substitute nonbreaking en dash
  while (retval.find("–") != string::npos) {
    size_t pos = retval.find("–");
    retval.replace(pos, 3, "&#8211;");//!< En Dash is 3 chars, not 1
  }
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
  
  ss << breseq_header_string(settings) << endl;
  ss << "<p>" << endl;
 
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

  HTML << "<!-- Write fastq read file informations -->" << endl;
  HTML << "<table border=\"0\" cellspace=\"1\" cellpadding=\"5\">" << endl;
  HTML << "<tr>" << th("fastq read file") << th("reads") << 
                    th("bases") << th("longest") << "</tr>" << endl;
// # //TODO @JEB Summary
// #   foreach my $read_file ($settings->read_files)
// #   {
// #     my $c = $summary->{sequence_conversion}->{reads}->{$read_file};
// #     print HTML Tr(
// #       td(a({-href=>$settings->html_path('error_rates_plot_file_name', {'#'=>$read_file})}, "errors")),
// #       td($read_file), 
// #       td({-align=>"right"}, commify($c->{num_reads})), 
// #       td({-align=>"right"},commify($c->{num_bases})), 
// #       td($c->{max_read_length} . "&nbsp;bases"), 
// #     );
// #   }
// #   print HTML Tr({-class=>"highlight_table_row"}, 
// #     td(),
// #     td(b("total")), 
// #     td({-align=>"right"},b(commify($summary->{sequence_conversion}->{num_reads}))), 
// #     td({-align=>"right"},b(commify($summary->{sequence_conversion}->{num_bases}))), 
// #     td(b($summary->{sequence_conversion}->{max_read_length} . "&nbsp;bases")), 
// #   );
// #   print HTML end_table();
// # 
// #   ## Write reference sequence information
// #   print HTML p; 
  
  //Write reference sequence information
  HTML << "<!-- Write reference sequence information -->" << endl;
  HTML << "<p>" << endl;
  HTML << "<table border=\"0\" cellspacing=\"1\" cellpadding=\"5\" >" << endl;
  HTML << "<tr>" << th() << 
                    th() << 
                    th("reference sequences") << 
                    th("length") << 
                    th("description") << 
          "</tr>" << endl;
             
  size_t total_length = 0;
  
// # //TODO @JEB Summary
// #   foreach my $seq_id (@{$ref_seq_info->{seq_ids}})
// #   {   
// #     my $c = $summary->{sequence_conversion}->{reference_sequences}->{$seq_id};
// #     
// #     print HTML Tr(
// #       td(a({-href=>$settings->html_path('coverage_plot_file_name', {'@'=>$seq_id})}, "coverage")), 
// #       td(a({-href=>$settings->html_path('unique_only_coverage_plot_file_name', {'@'=>$seq_id})}, "distribution")), 
// #       td($seq_id), 
// #       td({-align=>"right"},commify($c->{length})), 
// #       td($c->{definition})
// #     );
// #     $total_length+= $c->{length};
// #   } 
// #   
// #   print HTML Tr({-class=>"highlight_table_row"},
// #     td(),
// #     td(),
// #     td(b("total")), 
// #     td(b({-align=>"right"},commify($total_length))), 
// #     td()
// #   );
  //TODO move this into above loop when summary is implemented
  
  HTML << "<tr class=\"highlight_table_row\">";
  HTML << td();
  HTML << td();
  HTML << td(b("total"));
  HTML << td("<b align=\"right\">" + commify(to_string(total_length)) + "</b>");
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

// #     
// #     
  vector<ExecutionTime> times = settings.execution_times;
// #   ## Write Times
// HTML << "<!-- Write Times -->" << endl;
  HTML << "<p>"  << endl;
  HTML << h1("Execution Times") << endl;
  HTML << start_table("width=\"100%\" border=\"1\" cellspacing=\"0\" cellpadding=\"3\"") << endl;
  HTML << "<tr>" << th("Step") << th("Start") << th("End") << th("Elapsed") << "</tr>" << endl; 
  time_t total_time_elapsed = 0; 

  for (vector<ExecutionTime>::iterator itr = times.begin();
       itr != times.end(); itr ++) {  
    ExecutionTime& t = (*itr);
    
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
  HTML << "<td>" << (b(nonbreaking(to_string(&total_time_elapsed)))) << "</td>" << endl;
  HTML << "</tr>" << endl;

  HTML << "</table>";
  HTML.close();
}

string breseq_header_string(const Settings& settings)
{
  stringstream ss(ios_base::out | ios_base::app);
// #   
// #   #copy over the breseq_graphic
  string breseq_graphic_from_file_name = settings.breseq_small_graphic_from_file_name;
  string breseq_graphic_to_file_name= settings.breseq_small_graphic_to_file_name;
// #   if (!-e $breseq_graphic_to_file_name)
// #   {
// #     copy($breseq_graphic_from_file_name, $breseq_graphic_to_file_name);
// #   }
  if (breseq_graphic_to_file_name.empty()) {
    breseq_graphic_to_file_name = breseq_graphic_to_file_name;
  }
  ss << "<table width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"3\">" << endl;
  ss << "<tr>" << endl;
  ss << td(a(settings.website, img(settings.html_path(settings.breseq_small_graphic_to_file_name))));
  ss << endl;
  ss << start_td("width=\"100%\"") << endl;
// #   $output_string .= $settings->{byline} /TODO @JEB set byline to breseq version
  ss << settings.byline << endl;
  ss << "<br>";
  ss << a(settings.html_path("index_html_file_name"), "mutation predictions"); 
  ss << " | " << endl;
  ss << a(settings.html_path("marginal_html_file_name"), "marginal predictions");
  ss << " | " << endl;
  ss << a(settings.html_path("summary_html_file_name"), "summary statistics");
  ss << " | " << endl;
  ss << a(settings.html_path("final_genome_diff_file_name"), "genome diff");
  ss << " | " << endl;
  ss << a(settings.html_path("log_file_name"), "command line log");
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
  ss << tr(th("colspan=\"" + to_string(total_cols) + 
              "\" align=\"left\" class=\"read_alignment_header_row\"", title)) << endl;
  ss << "<tr>" << endl;
  
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
  ss << "<tr>" << th("colspan=\"" + to_string(total_cols) + "\" align=\"left\" class=\"missing_coverage_header_row\"", title) << "</tr>" << endl;   
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
         end += "–" + 
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

string html_new_junction_table_string(diff_entry_list& list_ref, bool show_reject_reason, const string& title, const string& relative_link
                               )
{
  stringstream ss; //!<< Main Build Object for Function
  diff_entry& test_item = *list_ref.front();

  bool link = (test_item.entry_exists(_SIDE_1_EVIDENCE_FILE_NAME) &&
               test_item.entry_exists(_SIDE_2_EVIDENCE_FILE_NAME) &&
               test_item.entry_exists(_NEW_JUNCTION_EVIDENCE_FILE_NAME));

  ss << start_table("border=\"0\" cellspacing=\"1\" cellpadding=\"3\"") << endl;
  size_t total_cols = link ? 10 : 8;
  ss << tr(th("colspan=\"" + to_string(total_cols) + "\" align=\"left\" class=\"new_junction_header_row\"", title)) << endl;
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
Evidence_Files::Evidence_Files(const Settings& settings, genome_diff& gd)
{  
  // Fasta and BAM files for making alignments.
  string reference_bam_file_name = settings.reference_bam_file_name;
  string reference_fasta_file_name = settings.reference_fasta_file_name;

  // hybrids use different BAM files for making the alignments!!!
  string junction_bam_file_name = settings.junction_bam_file_name;
  string junction_fasta_file_name = settings.candidate_junction_fasta_file_name;

  // We make alignments of two regions for deletions: upstream and downstream edges.
  diff_entry_list items_MC = gd.list(make_list<string>(MC));
  for (diff_entry_list::iterator itr = items_MC.begin(); itr != items_MC.end(); itr ++) 
  {  
    diff_entry& item = **itr;
    if (item.entry_exists(NO_SHOW)) continue;
     
    counted_ptr<diff_entry> parent_item = gd.parent(item);
    if (parent_item.get() == NULL)
      parent_item = *itr;

   add_evidence(_SIDE_1_EVIDENCE_FILE_NAME,
                item,
                 *parent_item,
                make_map<string,string>
                (BAM_PATH, reference_bam_file_name)
                (FASTA_PATH, reference_fasta_file_name)
                (PREFIX, "MC_SIDE_1")
                (SEQ_ID, item[SEQ_ID])            
                (START, to_string(from_string<uint8_t>(item[START]) - 1))
                (END,  to_string(from_string<uint8_t>(item[START]) - 1)));

    add_evidence(_SIDE_2_EVIDENCE_FILE_NAME,
                item,
                 *parent_item,
                make_map<string,string>
                (BAM_PATH, reference_bam_file_name)
                (FASTA_PATH, reference_fasta_file_name)
                (PREFIX, "MC_SIDE_2")
                (SEQ_ID, item[SEQ_ID])            
                (START, to_string(from_string<uint8_t>(item[END]) + 1))
                (END,  to_string(from_string<uint8_t>(item[END]) + 1)));

    add_evidence(_EVIDENCE_FILE_NAME,
                item,
                 *parent_item,
                make_map<string,string>
                (BAM_PATH, junction_bam_file_name)
                (FASTA_PATH, junction_fasta_file_name)
                (PREFIX, "MC_PLOT")
                (SEQ_ID, item[SEQ_ID])            
                (START, to_string(from_string<uint8_t>(item[END]) + 1))
                (END,  to_string(from_string<uint8_t>(item[END]) + 1))
                (PLOT, item[_COVERAGE_PLOT_FILE_NAME]));
    
  } // mc_item list
  
  diff_entry_list items_SNP_INS_DEL_SUB = gd.list(make_list<string>(SNP)(INS)(DEL)(SUB));

  for (diff_entry_list::iterator itr = items_SNP_INS_DEL_SUB.begin(); itr != items_SNP_INS_DEL_SUB.end(); itr ++) 
  {  
    diff_entry& item = **itr;
    diff_entry_list mutation_evidence_list = gd.mutation_evidence_list(item);
    if (item.entry_exists(NO_SHOW)) continue;

    // #this reconstructs the proper columns to draw
    uint32_t start = from_string<uint32_t>(item[POSITION]);
    uint32_t end = start;
    uint32_t insert_start = 0;
    uint32_t insert_end = 0;

    if (item._type == INS) 
    {
      insert_start = 1;
      insert_end = item[NEW_SEQ].size();
    }
    else if (item._type == DEL) 
    {
      bool has_ra_evidence;
      for (diff_entry_list::iterator itr = mutation_evidence_list.begin(); itr != mutation_evidence_list.end(); itr ++) 
      {  
        diff_entry& evidence_item = **itr;
        if (evidence_item._type == RA) has_ra_evidence = true;
      }
      if(!has_ra_evidence) continue;  

      end = start + from_string<uint32_t>(item[SIZE]) - 1;
    }

    else if (item._type == SUB ) 
    {
      end = start + item[NEW_SEQ].size() - 1;
    }

    add_evidence(_EVIDENCE_FILE_NAME,
                 item,
                 item,
                 make_map<string,string>
                 (BAM_PATH, reference_bam_file_name)
                 (FASTA_PATH, reference_fasta_file_name)
                 (SEQ_ID, item[SEQ_ID])            
                 (START, to_string(start))
                 (END,  to_string(end))
                 (INSERT_START, to_string(insert_start))
                 (INSERT_END, to_string(insert_end))
                 (PREFIX, item._type));


    /* Add evidence to RA items as well */
    for (diff_entry_list::iterator itr = mutation_evidence_list.begin(); itr != mutation_evidence_list.end(); itr ++) 
    {  
      diff_entry& evidence_item = **itr;
      if (evidence_item._type != RA) continue;
      evidence_item[_EVIDENCE_FILE_NAME] = item[_EVIDENCE_FILE_NAME];  
    }
  }
  

  // Still create files for RA evidence that was not good enough to predict a mutation from
  diff_entry_list ra_list = gd.filter_used_as_evidence(gd.list(make_list<string>(RA)));
  
  for (diff_entry_list::iterator itr = ra_list.begin(); itr != ra_list.end(); itr ++) 
  {  
    diff_entry& item = **itr;
    if (item.entry_exists(NO_SHOW)) continue;

    add_evidence(_EVIDENCE_FILE_NAME,
                 item,
                 item,
                 make_map<string,string>
                 (BAM_PATH, reference_bam_file_name)
                 (FASTA_PATH, reference_fasta_file_name)
                 (SEQ_ID, item[SEQ_ID])            
                 (START, item[POSITION])
                 (END, item[POSITION])
                 (INSERT_START, item[INSERT_POSITION])
                 (INSERT_END, item[INSERT_POSITION])
                 (PREFIX, item._type));
  }
  // This additional information is used for the complex reference line.
  // Note that it is completely determined by the original candidate junction sequence 
  // positions and overlap: alignment_pos and alignment_overlap.
  
  diff_entry_list items_JC = gd.list(make_list<string>(JC));

  for (diff_entry_list::iterator itr = items_JC.begin(); itr != items_JC.end(); itr ++) 
  {  
    diff_entry& item = **itr;
    if (item.entry_exists(NO_SHOW)) continue;

    diff_entry_ptr parent_item = gd.parent(item);

    if(parent_item.get() == NULL) {
      parent_item = *itr;
    }

    uint32_t start;
    uint32_t end;

    if (from_string<uint32_t>(item[ALIGNMENT_OVERLAP]) == 0) 
    {
      start = from_string<uint32_t>(item[FLANKING_LEFT]);
      end = from_string<uint32_t>(item[FLANKING_LEFT]) + 1;
    }
    else if (from_string <uint32_t>(item[ALIGNMENT_OVERLAP]) > 0) 
    {
      start = from_string<uint32_t>(item[FLANKING_LEFT]) + 1;
      end = from_string<uint32_t>(item[FLANKING_LEFT]) + 
            from_string<uint32_t>(item[ALIGNMENT_OVERLAP]);
    }
    else //if (from_string <uint32_t>(item[ALIGNMENT_OVERLAP]) > 0) 
    {
      start = from_string<uint32_t>(item[FLANKING_LEFT]) + 1;
      end = from_string<uint32_t>(item[FLANKING_LEFT]) - from_string<uint32_t>(item[ALIGNMENT_OVERLAP]);
    }
    
    add_evidence(_NEW_JUNCTION_EVIDENCE_FILE_NAME,
                  item,
                  *parent_item,
                  make_map<string,string>
                 (BAM_PATH, junction_bam_file_name)
                 (FASTA_PATH, junction_fasta_file_name)
                 (SEQ_ID, item["key"])              // @JEB TODO: Need to update this to SEQ_ID when making evidence file
                 (START, to_string(start))
                 (END, to_string(end))
                 (PREFIX, JC)
                 (ALIGNMENT_EMPTY_CHANGE_LINE, "1")
                 (TRUNCATE_END, to_string(
                                          from_string<uint32_t>(item[FLANKING_LEFT]) +
                                          from_string<uint32_t>(item[SIDE_1_OVERLAP])
                                         ))
                 (GHOST_END, item[SIDE_1_POSITION])
                 (GHOST_STRAND_END, item[SIDE_1_STRAND])
                 (GHOST_SEQ_ID_END, item[SIDE_1_SEQ_ID])
                 (TRUNCATE_START, to_string(
                                            from_string<uint32_t>(item[FLANKING_LEFT]) + 
                                            1 +
                                            abs(from_string<int32_t>(item[ALIGNMENT_OVERLAP]))
                                           ))
                 (GHOST_START, item[SIDE_2_POSITION])
                 (GHOST_STRAND_START, item[SIDE_2_STRAND])
                 (GHOST_SEQ_ID_START, item[SIDE_2_SEQ_ID]));



    // this is the flagship file that we show first when clicking on evidence from a mutation...
    item[_EVIDENCE_FILE_NAME] = item[_NEW_JUNCTION_EVIDENCE_FILE_NAME];
    string side_1_key_str;
    add_evidence(_SIDE_1_EVIDENCE_FILE_NAME,
                 item,
                 *parent_item,
                 make_map<string,string>
                 (BAM_PATH, reference_bam_file_name)
                 (FASTA_PATH, reference_fasta_file_name)
                 (SEQ_ID, item[SIDE_1_SEQ_ID])            
                 (START, item[SIDE_1_POSITION])
                 (END, item[SIDE_1_POSITION])
                 (PREFIX, "JC_SIDE_1" + '_' +
                          item[SIDE_2_SEQ_ID] + '_' +
                          item[SIDE_2_POSITION] + '_' +
                          item[SIDE_2_POSITION]
                 )); 

    add_evidence(_SIDE_2_EVIDENCE_FILE_NAME,
                 item,
                 *parent_item,
                 make_map<string,string>
                 (BAM_PATH, reference_bam_file_name)
                 (FASTA_PATH, reference_fasta_file_name)
                 (SEQ_ID, item[SIDE_2_SEQ_ID])            
                 (START, item[SIDE_2_POSITION])
                 (END, item[SIDE_2_POSITION])
                 (PREFIX, "JC_SIDE_2" + '_' +
                          item[SIDE_1_SEQ_ID] + '_' +
                          item[SIDE_1_POSITION] + '_' +
                          item[SIDE_1_POSITION]
                 ));
  }

  // now create evidence files
  create_path(settings.evidence_path);
  for (vector<Evidence_Item>::iterator itr = evidence_list.begin(); itr != evidence_list.end(); itr ++) 
  {  
    Evidence_Item& e = (*itr);

    if (settings.verbose) {
      cerr << "Creating evidence file: " + e[FILE_NAME] << endl;   
    }
    html_evidence_file(settings, gd, e);
  }
}

/*-----------------------------------------------------------------------------
 *  Helper Function For Create_Evidence_Files()
 *-----------------------------------------------------------------------------*/
void Evidence_Files::add_evidence(const string& evidence_file_name_key, diff_entry& item,
                                  diff_entry& parent_item, map<string,string> fields)
{
  Evidence_Item evidence_item;
  evidence_item._fields = fields;
  evidence_item.item = item;
  evidence_item.parent_item = parent_item;

  evidence_item[FILE_NAME] = file_name(evidence_item);
  evidence_item[evidence_file_name_key] = evidence_item[FILE_NAME];
  
  evidence_list.push_back(evidence_item);
}
/*-----------------------------------------------------------------------------
 *  Helper Function For Create_Evidence_Files()
 *-----------------------------------------------------------------------------*/
string Evidence_Files::file_name(Evidence_Item& evidence_item)
{
  stringstream ss(ios_base::out | ios_base::app);

  ss << evidence_item[PREFIX];
  ss << "_" << evidence_item[SEQ_ID];
  ss << "_" << evidence_item[START];
  ss << evidence_item.entry_exists(INSERT_START) ? "." + evidence_item[INSERT_START] : "";
  ss << "_" << evidence_item[END];
  ss << evidence_item.entry_exists(INSERT_END) ? "." + evidence_item[INSERT_END] : "";
  ss << "_alignment.html";
  
  return ss.str();
}


/*-----------------------------------------------------------------------------
 *  Create the HTML Evidence File
 *-----------------------------------------------------------------------------*/
// # 
// # 
void 
Evidence_Files::html_evidence_file (
                    const Settings& settings, 
                    genome_diff& gd, 
                    Evidence_Item& item
                   )
{
  item["output_path"] = settings.evidence_path + "/" + item[FILE_NAME];

  
  // Create Stream and Confirm It's Open
  ofstream HTML(item["output_path"].c_str());
  
  if (!HTML.good()) {
    cerr << "Could not open file: " << item["output_path"] << endl;
    assert(HTML.good());
  }
  
  // Build HTML Head
  HTML << html_header("BRESEQ :: Results", settings);
  
  // print a table for the main item
  // followed by auxiliary tables for each piece of evidence

  diff_entry parent_item = item.parent_item;
  diff_entry_list evidence_list = gd.mutation_evidence_list(parent_item);

  vector<string> types = make_list<string>("RA")("MC")("JC");
  
  for (vector<string>::iterator itr = types.begin(); itr != types.end(); itr ++) 
  {  
    string& type = (*itr);

    diff_entry_list this_evidence_list = evidence_list;
    this_evidence_list.remove_if(diff_entry::is_not_type(type));   
    
    if(this_evidence_list.empty()) continue;

    HTML << html_genome_diff_item_table_string(settings, gd, this_evidence_list);
    HTML << "<p>"; 
  }

  
  if (item.entry_exists(PLOT) && !item[PLOT].empty()) {
    HTML << div(ALIGN_CENTER, img(item[PLOT]));
  } else {
    stringstream ss;   
    ss << item[SEQ_ID] << ":" << item[START];
    if (item[INSERT_START].size() > 0) {
      ss << "." << item[INSERT_START];
    }
    ss << "-" << item[END];
    if (item[INSERT_END].size() > 0) {
      ss << "." << item[INSERT_END];
    }
    cerr << "Creating read alignment for region " << ss.str() << endl;

    if (settings.base_quality_cutoff != 0) {
      item["BASE_QUALITY_CUTOFF"] = to_string(settings.base_quality_cutoff);
    }
    
    alignment_output ao(item[BAM_PATH], item[FASTA_PATH], from_string<uint32_t>(ss.str()), settings.base_quality_cutoff);
     
    HTML << ao.html_alignment(ss.str());

  }
  HTML << endl << "</html>";
  HTML.close();
}


/*-----------------------------------------------------------------------------
 *  //End Create_Evidence_Files
 *-----------------------------------------------------------------------------*/
void save_text_deletion_file(string deletion_file_name, diff_entry_list& deletions_ref)
{
	ofstream DEL(deletion_file_name.c_str()); //or die "Could not open: $deletion_file_name";
	DEL << "seq_id\tstart\tend\n";
	for (diff_entry_list::iterator d = deletions_ref.begin(); d != deletions_ref.end(); d++)
		DEL << (**d)["seq_id"] << "\t" << (**d)["start"] << "\t" << (**d)["end"] << "\n";
	DEL.close();
}

void draw_coverage(Settings& settings, cReferenceSequences* ref_seq_info, genome_diff& gd)
{
	vector<string> mc_types = make_list<string>("MC");
	diff_entry_list mc = gd.list(mc_types);
	string drawing_format = "png";

	{
		create_path("coverage_plot_path");
		string coverage_plot_path = settings.coverage_plot_path;
		string deletions_text_file_name = settings.deletions_text_file_name;
		save_text_deletion_file(deletions_text_file_name, mc);

		for (uint32_t i = 0; i < ref_seq_info->size(); i++)
		{
			string seq_id = (*ref_seq_info)[i].m_seq_id;
			string this_complete_coverage_text_file_name = settings.file_name("complete_coverage_text_file_name", "@", seq_id);
			string command = settings.lib_path + "plot_coverage --drawing-format " + drawing_format + " -t " + coverage_plot_path + " -p " + settings.coverage_plot_path + " -i " + deletions_text_file_name + " -c " + this_complete_coverage_text_file_name + " --seq_id=" + seq_id;
			assert(system(command.c_str()) >= 0);

			// need to assign link names that correspond to what the R script is doing
			diff_entry_list this_deletions;
			diff_entry_list::iterator it;
			
			if (seq_id.size() > 0)
				for (it = mc.begin(); it != mc.end(); it++)
					if ((**it)["seq_id"] == seq_id)
						this_deletions.push_back(*it);

			uint32_t j = 1;
			for (it = this_deletions.begin(); it != this_deletions.end(); it++)
				(**it)["_coverage_plot_file_name"] = seq_id + to_string(j++) + "." + drawing_format;
		}
    remove_file(deletions_text_file_name);
	}
}

vector<ExecutionTime> execution_times;
string record_time(string name)
{
	time_t this_time = time(NULL);
// #   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($this_time); 
	string formatted_time = to_string(this_time);
	ExecutionTime new_time;
	new_time._name = name;
	new_time._formatted_time = formatted_time;
	new_time._time = this_time;
	new_time._time_elapsed = 0;
	new_time._formatted_time_elapsed = "";

	//if we had a previous time
	time_t time_since_last;
	//string time_since_last_formatted;
	if (execution_times.size() > 0)
	{
		time_since_last = this_time - execution_times[execution_times.size() - 1]._time;
		new_time._time_elapsed = time_since_last;
		new_time._formatted_time_elapsed = to_string(time_since_last);
	}
	execution_times.push_back(new_time);
	return formatted_time;
}
// # 
// # sub time2string
// # {
// #     my ($seconds) = @_;
// #     # Convert seconds to days, hours, minutes, seconds
// #     my @parts = gmtime($seconds);
// #     my $ret = '';
// #     if(sprintf("%4d",$parts[7])>0)
// #     {
// #         $ret .= sprintf("%4d",$parts[7]);
// #         $ret .= sprintf(" %s",($parts[7]>1)?'days':'day');
// #     }
// #     if(sprintf("%4d",$parts[2])>0)
// #     {
// #         $ret .= sprintf("%4d",$parts[2]);
// #         $ret .= sprintf(" %s",($parts[2]>1)?'hours':'hour');
// #     }
// #     if(sprintf("%4d",$parts[1])>0)
// #     {
// #         $ret .= sprintf("%4d",$parts[1]);
// #         $ret .= sprintf(" %s",($parts[1]>1)?'minutes':'minute');
// #     }
// #     if(sprintf("%4d",$parts[0])>0)
// #     {
// #         $ret .= sprintf("%4d",$parts[0]);
// #         $ret .= sprintf(" %s",($parts[0]>1)?'seconds':'second');
// #     }
// #     return $ret;
// # }
// # 
// # 
// # ## Rather than dealing with our own formats for saving 
// # ## many different kinds of summary statistics,
// # ## just use freeze and thaw to directly save and load
// # ## the perl data structures (generally hashes).
// # use Storable;
// # sub save_statistics
// # {
// #   my ($file_name, $data) = @_;
// #   store $data, "$file_name";
// # }
// # 
// # sub load_statistics
// # {
// #   my ($file_name) = @_;
// #   return retrieve($file_name);
// # }
// # 
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
// #   my $header_text = "Predicted mutation";
// #   $header_text .= "s" if (scalar @$list_ref > 1);
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
  
// #   
// #   my $total_cols; // Class parameter
  
// #   my $header_str = '';

// #   
// #   if ($settings->{lenski_format}) {
  if (settings.lenski_format) {
// #     my @header_list = split /\|/, $freq_header_list[0];
    vector<string> header_list = split(freq_header_list.front(), "|");
// #     my $header_rows = scalar @header_list;
    size_t header_rows = header_list.size() - 1; //!< -1 is necessary for C++
// #     
// #     $total_cols = 7 + scalar @freq_header_list;
    total_cols = 7 + freq_header_list.size() ;
// #     $total_cols += 1 if (!$one_ref_seq);
    if(!one_ref_seq) total_cols += 1; 
// #     $total_cols += 1 if (!$settings->{no_evidence});  ## evidence column 
    if(!settings.no_evidence) total_cols += 1;
// #       
// #     for (my $i=1; $i<=$header_rows; $i++)
// #     { 
    for (size_t i = 1; i <= header_rows; i++) {
     ss << "<!-- Header Line -->" << endl;
     ss << "<tr>" << endl;
// #       $header_str.= th("evidence") if (!$settings->{no_evidence}); 
      if(!settings.no_evidence)
        ss << th("evidence") << endl;
// #       $header_str.= th(nonbreaking("seq id")) if (!$one_ref_seq); 
      if(!one_ref_seq)
       ss << th(nonbreaking("seq id")) << endl;
// #     
// #       $header_str .= th(
// #         [
// #           ($header_rows == $i) ? "position" : "",
// #           ($header_rows == $i) ? "mutation" : "",
// #           ($header_rows == $i) ? "annotation" : "", 
// #           ($header_rows == $i) ? "gene" : ""
// #         ]
// #       );
      ss << th( (header_rows == i) ? "position" : "") << endl;
      ss << th( (header_rows == i) ? "mutation" : "") << endl;
      ss << th( (header_rows == i) ? "annotation" : "") << endl;
      ss << th( (header_rows == i) ? "gene" : "") << endl;
// #       $header_str.= th({-width => "100%"}, ($header_rows == $i) ? "description" : ""); 
      ss << th("width=\"100%\"", (header_rows == i) ? "description" : "") << endl;
// #       foreach my $freq_header_item (@freq_header_list) {
      for (vector<string>::iterator itr = freq_header_list.begin();
           itr != freq_header_list.end(); itr++) {
        string& freq_header_item(*itr);        
// #         
// #         my @header_list = split /\|/, $freq_header_item; 
        vector<string> header_list = split(freq_header_item, "|");        
// #         my $this_header_string = $header_list[$i-1];
        string this_header_string = header_list[i-1];
// #         $this_header_string =~ s/_/&nbsp;/g;
        while(this_header_string.find("_") != string::npos)
        {
          size_t pos = this_header_string.find("_");
          this_header_string.replace(pos, 1, "&nbsp;");//TODO confim "_" is 1 char
        }
// #         my $this_header_string_1 = $header_list[0];
        string this_header_string_1 = header_list.front();
// #         my $this_header_string_2 = $header_list[1];
        string this_header_string_2 = header_list[1];
// # 
// #         my $color = "black";
        string color = "black";  
// #         
// #         if ($this_header_string_1 eq "UC") {
// #           $color = "gray";
// #         } elsif ($this_header_string_1 eq "clade_1") {
// #           $color = "green";
// #         } elsif ($this_header_string_1 eq "clade_2") {
// #           $color = "blue";
// #         } elsif ( ($this_header_string_1 eq "clade_3") 
// #           && ( ($this_header_string_2 eq "ZDB483") 
// #           || ($this_header_string_2 eq "ZDB30") ) )
// #         {
// #           $color = "orange";
// #         } elsif ($this_header_string_1 eq "clade_3") {
// #           $color = "red";
// #         }   
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
//#     
// #         $header_str .= th({style=>"background-color:$color"}, $this_header_string);
        ss << th("style=\"background-color:" + 
                         color + "\"", this_header_string) << endl;  
// #         #$header_str .= th($this_header_string);
        ss << th(this_header_string) << endl;
// #       }
      }
// #     

// #       $header_str .= th(
// #         [
// #           ($header_rows == $i) ? "position" : "",
// #         ]
// #       );    
      ss << th(header_rows == i ? "position" : "") << endl;
// #       $header_str.= end_Tr; 
      ss << "</tr>" << endl;
// #     }
    }
// #     
// #   } else {
  } else {
// #     $total_cols = 5 + scalar @freq_header_list;
    total_cols = 5 + freq_header_list.size() - 1;
// #     $total_cols += 1 if (!$one_ref_seq);
    if (!one_ref_seq) {
    total_cols += 1;
    }
// #     $total_cols += 1 if (!$settings->{no_evidence});  ## evidence column   
    if (!settings.no_evidence) {
      total_cols += 1;
    }
// #       
// #     $header_str.= start_Tr();
    ss <<  "<tr>" << endl;
// #     $header_str.= th("evidence") if (!$settings->{no_evidence}); 
   if (!settings.no_evidence) {
     ss << th("evidence") << endl;
   } 
// #     $header_str.= th(nonbreaking("seq id")) if (!$one_ref_seq); 
   if(!one_ref_seq) {
     ss << th(nonbreaking("seq id")) << endl;
   }
// #     
    
// #     $header_str .= th(
// #       [
// #         "position",
// #         "mutation",
// #       ]
// #     );
   ss << th("position") << endl;
   ss << th("mutation") << endl;
// #     
// #     foreach my $freq_header_item (@freq_header_list) {
// #       $header_str .= th( [$freq_header_item] );
// #     } 
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
// #   
// #   if (!$settings->{no_header})
// #   {
// #     $output_str.= Tr(th({-colspan => $total_cols, -align => "left", -class=>"mutation_header_row"}, $header_text));
// #   }
  if(!settings.no_header) {
          (*this) += tr(th("colspan=\"" + to_string(total_cols) +
                "\" align=\"left\" class=\"mutation_header_row\"", header_text));
  }
// #   $output_str.= $header_str; //  
// #
 ss << endl;
 (*this) += ss.str(); 
  //(*this) += ss.str();
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
// #   my $row_num = 0;
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
      if (!evidence_string.empty()) evidence_string += "&nbsp;"; //TODO Confirm "if statement" needed?
      evidence_string += a(relative_link + evidence_item[_EVIDENCE_FILE_NAME],
        evidence_item._type);
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
  //
// #     
// #     ## marshal cells defined depending on mutation type
// #     # $evidence_string
   //TODO double check with perl
    string cell_seq_id = nonbreaking(mut[SEQ_ID]);
    string cell_position = commify(mut[POSITION]);
    string cell_mutation;
    string cell_mutation_annotation = nonbreaking(formatted_mutation_annotation(mut));
    string cell_gene_name = i(nonbreaking(mut[GENE_NAME]));
    string cell_gene_product = htmlize(mut["gene_product"]);
// #         
// #     if ($mut->{type} eq 'SNP') {
// #       $cell_mutation = "$mut->{_ref_seq}&rarr;$mut->{new_seq}";
// #     } elsif ($mut->{type} eq 'INS') {
// #       $cell_mutation = "+$mut->{new_seq}";
// #     } elsif ($mut->{type} eq 'DEL') {
// #       $cell_mutation = nonbreaking("&Delta;" . commify($mut->{size}) . " bp");
// #       my $annotation_str = '';
// #       $annotation_str = "between $mut->{between}" if ($mut->{between});
// #       $annotation_str = "$mut->{mediated}-mediated" if ($mut->{mediated});
// #       $annotation_str = nonbreaking($mut->{gene_position}) if (!$annotation_str); 
// #       $cell_mutation_annotation = nonbreaking($annotation_str);
// #     } elsif ($mut->{type} eq 'SUB') {
// #       $cell_mutation = nonbreaking("$mut->{size} bp&rarr;$mut->{new_seq}");
// #     } elsif ($mut->{type} eq 'CON') {
// #       $cell_mutation = nonbreaking("$mut->{size} bp&rarr;$mut->{region}");      

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
// #     } elsif ($mut->{type} eq 'MOB') {
// #       my $s;
// #       my $s_start = '';
// #       $s_start .= "+" . $mut->{ins_start} if ($mut->{ins_start});
// #       $s_start .= "&Delta;" . $mut->{del_start} if ($mut->{del_start});
// #       $s.= $s_start . " :: " if ($s_start);
// # 
// #       $s .= "$mut->{repeat_name} (";
// #       $s .= (($mut->{strand}==+1) ? '+' : (($mut->{strand}==-1) ? '&minus;' : '?'));
// #       $s .= ")";
// # 
// #       my $s_end = '';
// #       $s_end .= "&Delta;" . $mut->{del_end} if ($mut->{del_end});
// #       $s_end .= "+" . $mut->{ins_end} if ($mut->{ins_end});
// #       $s.= " :: " . $s_end if ($s_end);
// # 
// #       my $dup_str = ($mut->{duplication_size} >= 0) ? "+$mut->{duplication_size}" : "&Delta;" . abs($mut->{duplication_size});
// #       $s .= " ($dup_str) bp";     
// #       $cell_mutation = nonbreaking($s);
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
// #     } elsif ($mut->{type} eq 'INV') {
// #       $cell_mutation =  nonbreaking(commify($mut->{size}) . " bp inversion");
// #       $cell_gene_name = i(nonbreaking($mut->{gene_name_1})) . "&darr;" . i(nonbreaking($mut->{gene_name_2}));
// #       $cell_gene_product = htmlize($mut->{gene_product_1}) . "&darr;" . htmlize($mut->{gene_product_2});
    } else if (mut._type == INV) {
      cell_mutation = nonbreaking(commify(mut["size"]) + " bp inversion");
      cell_gene_name = i(nonbreaking(mut["gene_name_1"])) + "&darr;" +
                       i(nonbreaking(mut["gene_name_2"]));
      cell_gene_product = htmlize(mut["gene_product_1"]) + "&darr;" + 
                          htmlize(mut["gene_product_2"]);
// #     } elsif ($mut->{type} eq 'AMP') {       
// #       $cell_mutation = nonbreaking(commify($mut->{size}) . " bp x $mut->{new_copy_number}");
// #       $cell_mutation_annotation = ($mut->{new_copy_number} == 2) ? "duplication" : "amplification";       
// #     }
    } else if (mut._type == AMP) {
      cell_mutation = nonbreaking(commify(mut["size"]) + "bp x " + mut["new_copy_number"]);
      cell_mutation_annotation = 
        from_string<uint8_t>(mut["new_copy_number"]) == 2 ?
          "duplication" : "amplification";
    }
//  ###### PRINT THE TABLE ROW ####
    ss << endl << "<!-- Print The Table Row -->" << endl; 
    ss << start_tr("class=\"" + row_class + "\"") << endl;

    if (!settings.no_evidence) {
      ss << td(ALIGN_CENTER, evidence_string) << "<!-- Evidence -->" << endl;
    }
    if (!one_ref_seq) {
      ss << td(ALIGN_CENTER, cell_seq_id) << "<!-- Seq_Id -->" << endl;
    }
    ss << td(ALIGN_CENTER, cell_position) << "<!-- Position -->" << endl;
// #     $output_str.= td({align=>"center"}, $cell_mutation);
// #     if ($settings->{lenski_format}) {
// #       $output_str.= td({align=>"center"}, $cell_mutation_annotation);
// #       $output_str.= td({align=>"center"}, $cell_gene_name);
// #       $output_str.= td({align=>"left"}, $cell_gene_product);
// #     }
    ss << td(ALIGN_CENTER, cell_mutation) << "<!-- Cell Mutation -->" << endl;
    if (settings.lenski_format) {
      ss << "<!-- Lenski_Format -->" << endl;
      ss << td(ALIGN_CENTER, cell_mutation_annotation) << endl;
      ss << td(ALIGN_CENTER, cell_gene_name) << endl;
      ss << td(ALIGN_CENTER, cell_gene_product) << endl;
    }
// #     $output_str.= freq_cols(@freq_list);
    //Need if statement for C++
    if (freq_list.size() >= 1 && !freq_list[0].empty()) {
      ss << freq_cols(freq_list) << endl;
    }
// #     
// #     if ($settings->{lenski_format}) {
// #       $output_str.= td({align=>"right"}, $cell_position);
// #     } else {
// #       $output_str.= td({align=>"center"}, $cell_mutation_annotation);     
// #       $output_str.= td({align=>"center"}, $cell_gene_name);
// #       $output_str.= td({align=>"left"}, $cell_gene_product);
// #     }
// #     $output_str.= end_Tr;   
    if (settings.lenski_format) {
      ss << "<!-- Lenski Format -->" << endl;
      ss << td(ALIGN_CENTER, cell_position) << endl;
    } else {
      ss << td(ALIGN_CENTER, cell_mutation_annotation) << endl;
      ss << td(ALIGN_CENTER, cell_gene_name) << endl;
      ss << td(ALIGN_CENTER, cell_gene_product) << endl;
    }
    ss << "</tr>" << endl;
    
// #     ##### END TABLE ROW ####
    ss << "<!-- End Table Row -->" << endl;
// #   }
  }
// #   
// #   if ($legend_row) {
// #     $output_str.= start_Tr(); 
// #     $output_str.= td({-colspan=>$total_cols}, b("Evidence codes: RA = read alignment, MC = missing coverage, JC = new junction"));
// #     $output_str.= end_Tr; 
// #   }
  if (legend_row) {
    ss << "<tr>" << endl;
    ss << td("colspan=\"" + to_string(total_cols) + "\"", 
                    b("Evidence codes: RA = read alignment, MC = missing coverage, JC = new junction"));
    ss << "</tr>" << endl;
  }
    
// #   
// #   $output_str.= end_table;  
// # }
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
// #     sub freq_to_string
// #     {
// #       my ($freq) = @_;
// #       
// #       $freq = 1 if (!defined $freq);
// #       return 'H' if ($freq eq 'H');
// #       return '?' if ($freq eq '?');
// #       return '' if ($freq == 0);
// #       my $frequency_string;
// #       if ($freq == 1) {
// #         $frequency_string = "100%";
// #       }
// #       else { 
// #         $frequency_string = sprintf("%4.1f%%", $freq*100);
// #       }
// #       return $frequency_string;
// #     }
  ///TODO check that my change of freq=1 is okay
  if (freq == "H") {
    return "H";
  }
  if (from_string<uint8_t>(freq) == 0) {
    return "";
  }
  stringstream ss;
  if (from_string<uint8_t>(freq) == 1 || freq.empty()) {
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
// #     sub freq_cols:
// #     {
// #       my @freq_list = @_;
// #       my $output_str = '';
// #       
// #       foreach my $freq (@freq_list) {
// #         if ($settings->{shade_frequencies}) {
// #           my $bgcolor;
// #           $bgcolor = 'Blue' if ($freq == 1);
// #           if (defined $bgcolor) {
// #             $output_str .= td({align=>"right", bgcolor=>$bgcolor}, '&nbsp;');
// #           }
// #           else {
// #             $output_str .= td({align=>"right"}, '&nbsp;');
// #           }
// #         }
// #         else {
// #           $output_str .= td({align=>"right"}, freq_to_string($freq));
// #         }
// #       }
// #       return $output_str;
// #     }
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

