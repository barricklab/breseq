#include "libbreseq/output.h"
#include "libbreseq/anyoption.h"
#include "libbreseq/alignment_output.h"
#include "libbreseq/coverage_output.h"

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
const char* GENE_PRODUCT_HIDE="_gene_product_hide";
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
const char* SIDE_1_OVERLAP="side_1_overlap";
const char* CORRECTED_KEY="corrected_key";
  
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
string htmlize(const string& input) 
{
  string retval = input;
    
  /* substitute nonbreaking en dash */
  retval = substitute(retval, "–", "&#8211;");
  
  /* substitute nonbreaking hyphen */
  retval = substitute(retval, "-", "&#8209;");

  return retval;
}

/*-----------------------------------------------------------------------------
 *  These style definitions are included between the HTML <head> 
 *  tags of every genereated .html page.
 *-----------------------------------------------------------------------------*/
string header_style_string() 
{
  stringstream ss;
  ss << "body {font-family: sans-serif; font-size: 11pt;}"                 << endl;
  ss << "th {background-color: rgb(0,0,0); color: rgb(255,255,255);}"      << endl;
  ss << "table {background-color: rgb(1,0,0); color: rgb(0,0,0);}"         << endl;
  ss << "tr {background-color: rgb(255,255,255);}"                         << endl;
  ss << ".mutation_in_codon {color:red; text-decoration : underline;}"     << endl;
  ss << ".mutation_header_row {background-color: rgb(0,130,0);}"           << endl;
  ss << ".read_alignment_header_row {background-color: rgb(255,0,0);}"     << endl;
  ss << ".missing_coverage_header_row {background-color: rgb(0,100,100);}" << endl;
  ss << ".new_junction_header_row {background-color: rgb(0,0,155);}"       << endl;
  ss << ".copy_number_header_row {background-color: rgb(153,102,0);}"      << endl;
  ss << ".alternate_table_row_0 {background-color: rgb(255,255,255);}"     << endl;
  ss << ".alternate_table_row_1 {background-color: rgb(235,235,235);}"     << endl;
  ss << ".gray_table_row {background-color: rgb(230,230,245);}"            << endl;
  ss << ".polymorphism_table_row {background-color: rgb(160,255,160);}"    << endl;
  ss << ".highlight_table_row {background-color: rgb(192,255,255);}"       << endl;
  ss << ".reject_table_row {background-color: rgb(255,200,165);}"          << endl;
  ss << ".information_table_row {background-color: rgb(200,255,255);}"     << endl;
  ss << ".junction_repeat {background-color: rgb(255,165,0)}"              << endl;
  ss << ".junction_gene {}"                                                << endl;
  ss << ".hidden { display: none; }"                                       << endl;
  ss << ".unhidden { display: block; }"                                    << endl;
  
return ss.str();
}
  
/*-----------------------------------------------------------------------------
 *  Javascript to include inside the head.
 *-----------------------------------------------------------------------------*/
string javascript_string() 
{
  stringstream ss;
  ss << "<script type=\"text/javascript\">"                                     << endl;
  ss << "  function hideTog(divID) {"                                           << endl;
  ss << "    var item = document.getElementById(divID);"                        << endl;
  ss << "    if (item) {"                                                       << endl;
  ss << "      item.className=(item.className=='hidden')?'unhidden':'hidden';"  << endl;
  ss << "    }"                                                                 << endl;
  ss << "  }"                                                                   << endl;
  ss << "  function showTog(butID) {"                                           << endl;
  ss << "    var button = document.getElementById(butID);"                      << endl;
  ss << "    if (button) {"                                                     << endl;
  ss << "      button.value=(button.value=='Show')?'Hide':'Show';"              << endl;
  ss << "    }"                                                                 << endl;
  ss << "  }"                                                                   << endl;
  ss << "</script>"                                                             << endl;
  
  return ss.str();
}



void html_index(const string& file_name, const Settings& settings, Summary& summary,
                cReferenceSequences& ref_seq_info, cGenomeDiff& gd)
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
  diff_entry_list_t muts = gd.show_list(make_vector<gd_entry_type>(SNP)(INS)(DEL)(SUB)(MOB)(AMP));
  
  string relative_path = settings.local_evidence_path;
  
  if(!relative_path.empty())
    relative_path += "/";
  
  //Determine if more than one reference sequence is used
  bool one_ref_seq;

  if (ref_seq_info.size() == 1)
    one_ref_seq = true;
  else
    one_ref_seq = false;

  //Build Mutation Predictions table
  HTML << "<p>" << endl;
  MutationTableOptions mt_options;
  mt_options.relative_link = relative_path;
  mt_options.one_ref_seq = one_ref_seq;
  HTML << Html_Mutation_Table_String(settings, gd, muts, mt_options) << endl;

  //////
  // Unassigned evidence
  //////

  diff_entry_list_t cn = gd.filter_used_as_evidence(gd.show_list(make_vector<gd_entry_type>(CN)));
  cn.remove_if(cDiffEntry::rejected()); 
  
  if (cn.size() > 0) {
    HTML << "<p>" << html_copy_number_table_string(cn, false, "Unassigned copy number evidence", relative_path);
  }
  
  diff_entry_list_t ra = gd.filter_used_as_evidence(gd.show_list(make_vector<gd_entry_type>(RA)));
  ra.remove_if(cDiffEntry::rejected()); 
  ra.remove_if(cDiffEntry::field_exists("mixed"));
  if (ra.size() > 0) {
    HTML << "<p>" << html_read_alignment_table_string(ra, false, "Unassigned read alignment evidence", relative_path);
  }
  
  diff_entry_list_t mc = gd.filter_used_as_evidence(gd.show_list(make_vector<gd_entry_type>(MC)));
  mc.remove_if(cDiffEntry::rejected()); 

  if (mc.size() > 0) {
    HTML << "<p>" << html_missing_coverage_table_string(mc, false, "Unassigned missing coverage evidence", relative_path);
  }

  diff_entry_list_t jc = gd.filter_used_as_evidence(gd.show_list(make_vector<gd_entry_type>(JC)));
  jc.remove_if(cDiffEntry::rejected()); 

  //Don't show junctions for circular chromosomes
  if (settings.hide_circular_genome_junctions) {
    jc.remove_if(cDiffEntry::field_exists("circular_chromosome")); 
  }
   
  diff_entry_list_t jcu = jc;
  if (jcu.size() > 0) {
    HTML << "<p>" << endl;
    HTML << html_new_junction_table_string(jcu, false, "Unassigned new junction evidence...", relative_path);
  }
  
  HTML << html_footer();
  HTML.close();
}


void mark_gd_entries_no_show(const Settings& settings, cGenomeDiff& gd)
{
  /////
  // RA evidence
  //////
  
  vector<gd_entry_type> ra_types = make_vector<gd_entry_type>(RA);
  list<counted_ptr<cDiffEntry> > ra = gd.filter_used_as_evidence(gd.list(ra_types));
  ra.remove_if(not1(cDiffEntry::field_exists("reject")));
  ra.remove_if(cDiffEntry::field_exists("deleted"));
  ra.sort(cDiffEntry::by_scores(make_vector<string>("genotype_quality")));
  
  // Filtering for polymorphism predictions

  /* @JEB experimenting with removal
  // Require a certain amount of coverage for the variant
  diff_entry_list_t new_ra = gd.filter_used_as_evidence(gd.list(ra_types));
  for (diff_entry_list_t::iterator item = new_ra.begin(); item != new_ra.end(); item++)
  {
    vector<string> top_bot = split((**item)["tot_cov"], "/");
    uint32_t top = from_string<int32_t>(top_bot[0]);
    uint32_t bot = from_string<int32_t>(top_bot[1]);
    if (top + bot <= 2) (**item)["no_show"] = "1";
  }
  ra.remove_if(cDiffEntry::no_show());
  */
  
  // We mark polymorphism entries above some limit as no_show 
  // to prevent creating too many alignments.
  uint32_t i=0;
  uint32_t j=0;
  for(diff_entry_list_t::iterator item = ra.begin(); item != ra.end(); item++)
  {
    cDiffEntry& ra_item = **item;
    // This is a consensus prediction
    if ( (!ra_item.entry_exists(FREQUENCY)) || (from_string<double>(ra_item[FREQUENCY]) == 1.0 ) ) {
      if (j++ > settings.max_rejected_read_alignment_evidence_to_show)
        ra_item[NO_SHOW] = "1";
    }
    // This is a polymorphism prediction
    else if (i++ > settings.max_rejected_polymorphisms_to_show)
      ra_item[NO_SHOW] = "1";
  }
  
  /////
  // JC evidence
  //////
  
  vector<gd_entry_type> jc_types = make_vector<gd_entry_type>(JC);
  
  diff_entry_list_t jc = gd.filter_used_as_evidence(gd.list(jc_types));
  
  // This is the correct way to erase from a list! @MDS
  for (diff_entry_list_t::iterator it = jc.begin(); it != jc.end(); )
  {
    if (((*it)->number_reject_reasons() == 0))
      it = jc.erase(it);
    else
      it++;
  }
  
  jc.sort(cDiffEntry::by_scores(make_vector<diff_entry_key_t>("pos_hash_score")("total_non_overlap_reads")));
  
  i = 0;
  for (diff_entry_list_t::iterator it = jc.begin(); it != jc.end(); it++ )
  {
    if (i++ >= settings.max_rejected_junction_evidence_to_show)
      (**it)["no_show"] = "1";
  }


}
  
void html_marginal_predictions(const string& file_name, const Settings& settings,Summary& summary,
                               cReferenceSequences& ref_seq_info, cGenomeDiff& gd)
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
  if (ref_seq_info.size() == 1)
    one_ref_seq = true;
  else
    one_ref_seq = false; 

  // ###
  // ## Marginal evidence
  // ###
  
  // RA evidence
  vector<gd_entry_type> ra_types = make_vector<gd_entry_type>(RA);
  list<counted_ptr<cDiffEntry> > mixed_ra = gd.filter_used_as_evidence(gd.show_list(ra_types));
  mixed_ra.remove_if(not1(cDiffEntry::field_exists("mixed")));
  if (mixed_ra.size() > 0) {
    HTML << "<p>" << endl;
    HTML << html_read_alignment_table_string(mixed_ra, false, "Possible mixed read alignment sites...", relative_path) << endl;
  }
  
  list<counted_ptr<cDiffEntry> > ra = gd.filter_used_as_evidence(gd.show_list(ra_types));
  ra.remove_if(not1(cDiffEntry::field_exists("reject")));
  ra.sort(cDiffEntry::by_scores(make_vector<string>("quality")));
  if (ra.size() > 0) {
    HTML << "<p>" << endl;
    HTML << html_read_alignment_table_string(ra, false, "Marginal read alignment evidence...", relative_path) << endl;
  }
  
  vector<gd_entry_type> jc_types = make_vector<gd_entry_type>(JC);
  diff_entry_list_t jc = gd.filter_used_as_evidence(gd.show_list(jc_types));
  jc.remove_if(not1(cDiffEntry::field_exists("reject")));
  if (jc.size()) {
    //Sort by score, not by position (the default order)...
    jc.sort(cDiffEntry::by_scores(
     make_vector<diff_entry_key_t>("pos_hash_score")("total_non_overlap_reads")));

    HTML << "<p>" << endl;
    HTML << html_new_junction_table_string(jc, false, "Marginal new junction evidence...", relative_path);
  }

  HTML << html_footer();
  HTML.close();
}

string html_header (const string& title, const Settings& settings)
{
  stringstream ss(ios_base::out | ios_base::app);  
  
  ss << "<!DOCTYPE html" << endl;
	ss << "PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"" << endl;
  ss << "\"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">" << endl;
  ss << "<html xmlns=\"http://www.w3.org/1999/xhtml\" lang=\"en-US\" xml:lang=\"en-US\">" << endl;
  
  ss << "<html>" << endl;  
  ss << "<head>" << endl;
  ss << "<title>";
  if (!settings.print_run_name.empty()) {
    ss << settings.print_run_name << " :: ";
  }
  ss << title;
  ss << "</title>" << endl;
  
  ss << "<style type = \"text/css\">" << endl;
  ss << header_style_string() << endl;
  ss << "</style>" << endl;
  ss << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\" />" << endl;
  ss << javascript_string() << endl;
  ss << "</head>" << endl;
  ss << "<body>" << endl;
  return ss.str();
}



void html_compare(
                  const Settings& settings,
                  const string &file_name, 
                  const string &title, 
                  cGenomeDiff& gd,
                  MutationTableOptions& mt_options
                  )
{
  // Create stream and confirm it's open
  ofstream HTML(file_name.c_str());
  
  if(!HTML.good()) {
    cerr << "Could not open file: " <<  file_name << endl;
    assert(HTML.good());
  }

  diff_entry_list_t list_ref = gd.mutation_list();
  
  //Build html head
  HTML << html_header(title, settings);
  HTML << Html_Mutation_Table_String(settings, gd, list_ref, mt_options);
  HTML << html_footer();
  HTML.close();
}
  

void html_statistics(const string &file_name, const Settings& settings, Summary& summary, cReferenceSequences& ref_seq_info)
{  
  // Create stream and confirm it's open
  ofstream HTML(file_name.c_str());
  ASSERT(HTML.good(), "Could not open file: " + file_name);

  //Build html head
  HTML << html_header("BRESEQ :: Summary Statistics", settings);
  HTML << breseq_header_string(settings) << endl;
  HTML << "<p>" << endl;

  ////
  // Write read file information
  ////
  
  HTML << h2("Read File Information") << endl;
  HTML << start_table("border=\"0\" cellspace=\"1\" cellpadding=\"5\"") << endl;
  HTML << start_tr() << th() << th("read file") << th("reads") << 
                    th("bases") << th("longest") << th("mapped") << "</tr>" << endl;
  for(cReadFiles::const_iterator it=settings.read_files.begin(); it!=settings.read_files.end(); it++)
  {
    const Summary::AnalyzeFastq& s = summary.sequence_conversion.reads[it->m_base_name];
    const Summary::AlignmentResolution::ReadFile& rf = summary.alignment_resolution.read_file[it->m_base_name];
    HTML << start_tr();
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
    double percent_mapped = 100 * (1.0 - static_cast<double>(rf.num_unmatched_reads) / rf.num_total_reads);
    HTML << td(ALIGN_RIGHT, to_string(percent_mapped, 1) + "%");
    HTML << end_tr();
  }
  
  HTML << start_tr("class=\"highlight_table_row\"");
  HTML << td();
  HTML << td(b("total"));
  HTML << td(ALIGN_RIGHT , b(commify(to_string(summary.sequence_conversion.num_reads))) );
  HTML << td(ALIGN_RIGHT , b(commify(to_string(summary.sequence_conversion.num_bases))) );
  HTML << td(b(commify(to_string(summary.sequence_conversion.max_read_length))) + "&nbsp;bases");
  double total_percent_mapped = 100 * (1.0 - static_cast<double>(summary.alignment_resolution.total_unmatched_reads) / summary.alignment_resolution.total_reads);
  HTML << td(ALIGN_RIGHT, to_string(total_percent_mapped, 1) + "%");
  HTML << end_tr();
  HTML << end_table();
  
  ////
  // Write reference sequence information
  ////
  
  HTML << h2("Reference Sequence Information") << endl;
  HTML << "<p>" << endl;
  HTML << "<table border=\"0\" cellspacing=\"1\" cellpadding=\"5\" >" << endl;
  HTML << "<tr>" << th() << 
                    th() << 
                    th("seq id") << 
                    th("length") <<
                    th(ALIGN_CENTER, "fit mean") <<
                    th(ALIGN_CENTER, "fit dispersion") <<
 //                   th(ALIGN_CENTER, "nbinom mean") <<
 //                   th(ALIGN_CENTER, "nbinom size") <<
                    th(ALIGN_LEFT, "description") <<
          "</tr>" << endl;
             
  size_t total_length = 0;
  for(cReferenceSequences::iterator it=ref_seq_info.begin(); it!=ref_seq_info.end(); it++)
  {
    total_length += it->m_length;
    
    bool fragment_with_fit_coverage = (summary.unique_coverage[it->m_seq_id].nbinom_mean_parameter != 0);
    bool fragment_with_no_coverage = (summary.unique_coverage[it->m_seq_id].average == 0);

    HTML << (fragment_with_fit_coverage ? "<tr>" : "<tr class=\"gray_table_row\">");
    HTML << td( a(Settings::relative_path( 
                                          settings.file_name(settings.overview_coverage_plot_file_name, "@", it->m_seq_id), settings.output_path
                                          ), 
                  "coverage" 
                  )
               );
    
    // There may be absolutely no coverage and no graph will exist...
    if (!fragment_with_no_coverage)
    {
      HTML << td( a(Settings::relative_path( 
                                          settings.file_name(settings.unique_only_coverage_plot_file_name, "@", it->m_seq_id), settings.output_path
                                          ), 
                  "distribution" 
                  )
               ); 
    }
    else
    {
      HTML << td(nonbreaking("none aligned"));
    }
    
    HTML << td(it->m_seq_id);
    HTML << td(ALIGN_RIGHT, commify(to_string(it->m_length)));
    
    HTML << td(ALIGN_CENTER, to_string(summary.unique_coverage[it->m_seq_id].nbinom_mean_parameter, 1));
    HTML << td(ALIGN_CENTER, to_string(summary.unique_coverage[it->m_seq_id].nbinom_dispersion));
    //HTML << td(ALIGN_CENTER, to_string(summary.unique_coverage[it->m_seq_id].nbinom_mean_parameter, 1));
    //HTML << td(ALIGN_CENTER, to_string(summary.unique_coverage[it->m_seq_id].nbinom_size_parameter, 1));

    HTML << td(it->m_description);
    HTML << "</tr>";
  }  
  
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
  
  HTML << "<tr class=\"highlight_table_row\">";
  HTML << td();
  HTML << td();
  HTML << td(b("total"));
  HTML << td(ALIGN_RIGHT, b(commify(to_string(total_length))) );
  HTML << td();
  HTML << td();
  HTML << td();
  HTML << "</tr>" << endl;
  HTML << "</table>" << endl;

  HTML << "<p>" << "<b>fit dispersion</b> is the ratio of the variance to the mean for the negative binomial fit. It is =1 for Poisson and >1 for over-dispersed data." << endl;
  
  ////
  // Junction evidence
  ////
  
  if (!settings.no_junction_prediction)
  {
    HTML << h2("Junction Evidence") << endl;
    HTML << "<p>" << endl;
    
    //stringstream num;
    //num << scientific << setprecision(2) << settings.junction_accept_pr;
    //HTML << "E-value cutoff: " <<  num.str() << endl;
    //HTML << "<p>" << endl;
    HTML << start_table("border=\"0\" cellspacing=\"1\" cellpadding=\"5\"") << endl;
    HTML << start_tr() <<
    th("reference sequence") << 
    th("pr(no read start)") <<
//    th("pos hash score cutoff") <<
//    th("distance score cutoff") <<
    end_tr() << endl;
    
    size_t total_length = 0;
    for(cReferenceSequences::iterator it=ref_seq_info.begin(); it!=ref_seq_info.end(); it++)
    {
      HTML << start_tr() << endl;
      HTML << td(it->m_seq_id);
      // this score is always an integer for now, even though it is typed as a double
      bool fragment_with_fit_coverage = (summary.unique_coverage[it->m_seq_id].nbinom_mean_parameter != 0);

      if (!fragment_with_fit_coverage) {
        HTML << td(ALIGN_CENTER, "NA"); 
      } else {
        HTML << td(ALIGN_CENTER, to_string(summary.preprocess_error_count[it->m_seq_id].no_pos_hash_per_position_pr, 5));
//        HTML << td(ALIGN_CENTER, to_string(summary.alignment_resolution.pos_hash_cutoffs[it->m_seq_id].back()
//        HTML << td(ALIGN_CENTER, to_string(summary.alignment_resolution.distance_cutoffs[it->m_seq_id]));
      }
      HTML << end_tr();
    }  
    HTML << end_table();
    
    HTML << "<p>" << "pr(no read start) is the probability that there will not be an aligned read whose first base matches a given position on a given strand." << endl;
  }
  
  //
  // Read alignment evidence options
  //
  double polymorphism_log10_e_value_cutoff;               // Default = mutation_log10_e_value_cutoff = 10
  double polymorphism_bias_p_value_cutoff;                // Default = 0.05
  double polymorphism_frequency_cutoff;                   // Default = 0.1 for mixed base | 0.0 for polymorphism
  int32_t polymorphism_minimum_new_coverage_each_strand; // Default = 1
  uint32_t polymorphism_reject_homopolymer_length;        // Default = 0 (OFF)
  bool no_indel_polymorphisms;                            // Default = false
  
  
  {
    HTML << h2("Read Alignment Evidence") << endl;
    HTML << "<p>" << endl;
    HTML << start_table("border=\"0\" cellspacing=\"1\" cellpadding=\"5\"") << endl;
    HTML << tr(th("option") + th("value"));

    string mode = "Consensus";
    if (settings.polymorphism_prediction) {
      mode = "Full Polymorphism";
    } else {
      mode = "Consensus/Mixed Base";
    }
    HTML << tr(td("Mode") + td(mode));
    HTML << tr(td("Ploidy") + td("1 (haploid)"));
    
    HTML << tr(td("Consensus mutation E-value cutoff") + td(s(settings.mutation_log10_e_value_cutoff)));


    HTML << tr(td("Polymorphism E-value cutoff") + td(s(settings.polymorphism_log10_e_value_cutoff)));
    HTML << tr(td("Polymorphism frequency cutoff")
               + td((settings.polymorphism_frequency_cutoff == 0) ? "OFF" : s(settings.polymorphism_frequency_cutoff))
               );
    HTML << tr(td("Polymorphism minimum coverage each strand") 
               + td((settings.polymorphism_minimum_new_coverage_each_strand == 0) ? "OFF" : s(settings.polymorphism_minimum_new_coverage_each_strand))
               );
    HTML << tr(td("Polymorphism homopolymer length cutoff") 
               + td((settings.polymorphism_reject_homopolymer_length == 0) ? "OFF" : s(settings.polymorphism_reject_homopolymer_length))
               );
    HTML << tr(td("Polymorphism bias cutoff") 
               + td((settings.polymorphism_bias_p_value_cutoff == 0) ? "OFF" : s(settings.polymorphism_bias_p_value_cutoff))
               );
    HTML << tr(td("Polymorphism predict indel polymorphisms") 
               + td((settings.polymorphism_reject_homopolymer_length == 0) ? "OFF" : s(settings.polymorphism_reject_homopolymer_length))
               );

    HTML << end_table();    
  }
  
  ////
  // Custom deletion coverage cutoffs
  ////
  
  if (settings.deletion_coverage_propagation_cutoff || settings.deletion_coverage_seed_cutoff)
    HTML << html_deletion_coverage_values_table_string(settings, ref_seq_info, summary);
  
  
  ////
  // Write Execution Times
  ////
  
  const vector<ExecutionTime>& times = settings.execution_times;
  // HTML << "<!-- Write Times -->" << endl;
  HTML << "<p>"  << endl;
  HTML << h2("Execution Times") << endl;
  HTML << start_table("border=\"0\" cellspacing=\"1\" cellpadding=\"5\"") << endl;
  HTML << "<tr>" << th("step") << th("start") << th("end") << th("elapsed") << "</tr>" << endl; 
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
  HTML << html_footer();
  HTML.close();
}

string html_deletion_coverage_values_table_string(const Settings& settings, cReferenceSequences& ref_seq_info, Summary& summary)
{
  stringstream ss; //!< Main Build Object in Function

  ss << "<p>";
  ss << h2("Coverage Values");
  ss << start_table("border=\"1\" cellspacing=\"0\" cellpadding=\"3\"");
  // Table Header
  ss << "<tr>";
  ss << th("Parameter");

  cReferenceSequences::iterator it;
  for(it = ref_seq_info.begin(); it != ref_seq_info.end(); it++)
      ss << th(it->m_seq_id);

  ss << "</tr>";

  // Table Rows  s
  ss << "<tr>";
  ss << td("Calculated average");
  for(it = ref_seq_info.begin(); it != ref_seq_info.end(); it++)
    ss << td(to_string(summary.unique_coverage[it->m_seq_id].average));
  ss << "</tr>";

  ss << "<tr>";
  ss << td("Calculated propagation cutoff");
  for(it = ref_seq_info.begin(); it != ref_seq_info.end(); it++)
    ss << td(to_string(summary.unique_coverage[it->m_seq_id].deletion_coverage_propagation_cutoff));
  ss << "</tr>";

  ss << "<tr>";
  ss << td("Calculated seed cutoff");
  for(it = ref_seq_info.begin(); it != ref_seq_info.end(); it++)
    ss << td(to_string(summary.unique_coverage[it->m_seq_id].deletion_coverage_seed_cutoff));
  ss << "</tr>";


  ss << "<tr>";
  ss << td("User defined propagation cutoff");
  if (settings.deletion_coverage_propagation_cutoff < 1)
    ss << td(to_string(settings.deletion_coverage_propagation_cutoff) + "%");
  else
    ss << td(to_string(settings.deletion_coverage_propagation_cutoff));

  for(size_t i = 0; i < ref_seq_info.size() - 1; i++)
    ss << td();
  ss << "</tr>";

  ss << "<tr>";
  ss << td("User defined seed cutoff");
  if (settings.deletion_coverage_seed_cutoff < 1)
    ss << td(to_string(settings.deletion_coverage_seed_cutoff) + "%");
  else
    ss << td(to_string(settings.deletion_coverage_seed_cutoff));

  for(size_t i = 0; i < ref_seq_info.size() - 1; i++)
    ss << td();
  ss << "</tr>";



  ss << "</table>";

  return ss.str();
}

string breseq_header_string(const Settings& settings)
{
  stringstream ss(ios_base::out | ios_base::app);
  
  //copy over the breseq_graphic which we need if it doesn't exist - don't show command
  if (!file_exists(settings.breseq_small_graphic_to_file_name.c_str())) {
    copy_file(settings.breseq_small_graphic_from_file_name, settings.breseq_small_graphic_to_file_name);
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


string html_genome_diff_item_table_string(const Settings& settings, cGenomeDiff& gd, diff_entry_list_t& list_ref)
{
  if(list_ref.empty()) return "";

  cDiffEntry& first_item = *list_ref.front();
  //mutation
  if(first_item.is_mutation())
  {
    MutationTableOptions mt_options;
    mt_options.detailed=true;
    return Html_Mutation_Table_String(settings, gd, list_ref, mt_options); 
  }
  //evidence
  else
  {
    if(first_item._type == MC)
    {
      return html_missing_coverage_table_string(list_ref, true);
    }
    else if(first_item._type == RA)
    {
      return html_read_alignment_table_string(list_ref, true);
    }
    else if(first_item._type == JC)
    {
      return html_new_junction_table_string(list_ref,true);
    }
    else if(first_item._type == CN)
    {
      return html_copy_number_table_string(list_ref,true);
    }
  }  
  return "";
}

/*-----------------------------------------------------------------------------
 *  FORMATTED_MUTATION_ANNOTATION
 *-----------------------------------------------------------------------------*/
string formatted_mutation_annotation(const cDiffEntry& mut)
{
  stringstream ss;
  
  // additional formatting for some variables
  if((mut.entry_exists("snp_type")) && (mut.get("snp_type") != "intergenic") &&
     (mut.get("snp_type") != "noncoding") && (mut.get("snp_type") != "pseudogene"))
  {    
    ss << mut.get("aa_ref_seq") << mut.get("aa_position") << mut.get("aa_new_seq");
    string codon_ref_seq = to_underline_red_codon(mut, "codon_ref_seq");
    string codon_new_seq = to_underline_red_codon(mut, "codon_new_seq");
    ss << "&nbsp;(" << codon_ref_seq << "&rarr;" << codon_new_seq << ")&nbsp;";  
  }
  else // mut[SNP_TYPE] == "NC"
  {
    if(mut.entry_exists(GENE_POSITION)){
      ss << nonbreaking(mut.get(GENE_POSITION));
    }
  }
  return ss.str(); 
}

/*-----------------------------------------------------------------------------
 *  Helper function for formatted_mutation_annotation
 *-----------------------------------------------------------------------------*/
string to_underline_red_codon(const cDiffEntry& mut, const string& codon_key)
{
  if (!mut.entry_exists(codon_key) || 
      !mut.entry_exists("codon_position") ||
      mut.get("codon_position") == "") {
    return "";
  }

  stringstream ss; //!< codon_string
  
  string codon_ref_seq = mut.get(codon_key);
  // codon_position is 1-indexed
  uint32_t codon_position = from_string<int32_t>(mut.get("codon_position"));
  for (size_t i = 0; i < codon_ref_seq.size(); i++) {
    if (i+1 == codon_position) {
      ss << font("class=\"mutation_in_codon\"", codon_ref_seq.substr(i,1));
    }
    else 
    {
      ss << codon_ref_seq[i];
    }
  }
  return ss.str();
}

string html_read_alignment_table_string(diff_entry_list_t& list_ref, bool show_details, const string& title, const string& relative_link)
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
        th("reads")      << endl <<
        th("annotation") << endl <<
        th("genes")      << endl;
  
  ss << th("width=\"100%\"", "product") << endl;
  ss << "</tr>" << endl;
  
  //Loop through list_ref to build table rows
  for (diff_entry_list_t::iterator itr = list_ref.begin();
       itr != list_ref.end(); itr ++) {  
    cDiffEntry& c = **itr;
    
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
      string log_fisher = "";
      string log_ks = "";
         
      if (c.entry_exists(FISHER_STRAND_P_VALUE) &&
          from_string<double>(c[FISHER_STRAND_P_VALUE]) > 0 )
        log_fisher = to_string(log(from_string<double>(c[FISHER_STRAND_P_VALUE])));

      if (c.entry_exists(KS_QUALITY_P_VALUE) &&
          from_string<double>(c[KS_QUALITY_P_VALUE]) > 0 ) 
        log_ks = to_string(log(from_string<double>(c[KS_QUALITY_P_VALUE]))/log(10));

      ssf << fixed << setprecision(1);
      ssf << from_string<double>(c["quality"]);
      if (log_fisher != "") ssf << " " << log_fisher;
      if (log_ks != "") ssf << " " << log_ks;
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
    ss << td(ALIGN_CENTER, formatted_mutation_annotation(c)); //"Annotation" Column DON'T call nonbreaking on the whole thing
    ss << td(ALIGN_CENTER, i(nonbreaking(c[GENE_NAME])));
    ss << td(ALIGN_LEFT, htmlize(c[GENE_PRODUCT]));
    ss << "</tr>" << endl;

    if (show_details) 
    {
      vector<string> reject_reasons = c.get_reject_reasons();
      for (vector<string>::iterator itr = reject_reasons.begin(); itr != reject_reasons.end(); itr ++) 
      {  
        string& reject = (*itr);
       ss << tr("class=\"reject_table_row\"", 
                td("colspan=\"" + to_string(total_cols) + "\"",
                   "Rejected: " + decode_reject_reason(reject)));
      }
      
      // Show information about the strands supporting the change
      ss << tr("class=\"information_table_row\"", 
               td("colspan=\"" + to_string(total_cols) + "\"",
                  "Reads supporting (aligned to +/- strand):&nbsp;&nbsp;" +
                  b("new") + " base " + "(" + c[NEW_COV] + ")" + ":&nbsp;&nbsp;" +
                  b("ref") + " base " + "(" + c[REF_COV] + ")" + ":&nbsp;&nbsp;" +
                  b("total") + " (" + c[TOT_COV] + ")")); 
    
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
                 td("colspan=\"" + to_string(total_cols) + "\"",
                    "Fisher's exact test strand distribution " +
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
                    " Kolmogorov-Smirnov test that lower quality scores support polymorphism than reference " +
                    i("p") + "-value = " +ks_quality_p_value));
      } //end ks_quality_p_value
    } // end show_details
  } // end list_ref loop

  ss << "</table>" << endl;
  return ss.str();
}

string html_missing_coverage_table_string(diff_entry_list_t& list_ref, bool show_details, const string& title, const string& relative_link)
{
  ASSERT(list_ref.front().get(), "No items in table");
  
  stringstream ss; //!< Main Build Object in Function
  
  ss << endl;
  ss << start_table("border=\"0\" cellspacing=\"1\" cellpadding=\"3\" width=\"100%\"") << endl;
  
  bool coverage_plots = list_ref.front()->entry_exists(_EVIDENCE_FILE_NAME);
    
  bool link = ( list_ref.front()->entry_exists(_SIDE_1_EVIDENCE_FILE_NAME) && 
                list_ref.front()->entry_exists(_SIDE_2_EVIDENCE_FILE_NAME) );
  
  size_t total_cols = link ? 11 : 8;

  if (title != "") {
    ss << "<tr>" << th("colspan=\"" + to_string(total_cols) + "\" align=\"left\" class=\"missing_coverage_header_row\"", title) << "</tr>" << endl;   
  }
  
  ss << "<tr>";
    
  if (link) 
  {
    ss << th("&nbsp;") <<  th("&nbsp;");
    if (coverage_plots) 
    {
      ss << th("&nbsp;");
    }
  }

  ss << th("seq&nbsp;id") << endl <<
        th("start")       << endl <<
        th("end")         << endl <<
        th("size")        << endl <<
        th("&larr;reads")   << endl <<
        th("reads&rarr;")   << endl <<
        th("gene")        << endl;
  
  ss << th("width=\"100%\"", "description") << endl;
  ss << "</tr>" << endl;

  for (diff_entry_list_t::iterator itr = list_ref.begin(); itr != list_ref.end(); itr ++) 
  {  
    cDiffEntry& c =  **itr;

    ss << "<tr>" << endl;
    if (link) 
    {
      ss << td(a(relative_link + c[_SIDE_1_EVIDENCE_FILE_NAME], "*")) << endl;
      ss << td(a(relative_link + c[_SIDE_2_EVIDENCE_FILE_NAME], "*")) << endl; 
      
      if (coverage_plots) 
        ss << td(a(relative_link + c[_EVIDENCE_FILE_NAME], "&divide;")) << endl;
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
      
    if (
        (from_string<uint32_t>(c[END_RANGE]) > 0) 
        || (from_string<uint32_t>(c[START_RANGE]) > 0)
        ) 
    {
     
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
    ss << td(ALIGN_CENTER, nonbreaking(c[LEFT_OUTSIDE_COV] + " [" + c[LEFT_INSIDE_COV] + "]")) <<endl;
    ss << td(ALIGN_CENTER, nonbreaking( "[" + c[RIGHT_INSIDE_COV] + "] " + c[RIGHT_OUTSIDE_COV])) << endl;
    ss << td(ALIGN_CENTER, i(nonbreaking(c[GENE_NAME]))) << endl;
    ss << td(ALIGN_LEFT, htmlize(c[GENE_PRODUCT])) << endl;
    ss << "</tr>" << endl;

    if (show_details && c.entry_exists(REJECT)) 
    {
      vector<string> reject_reasons = c.get_reject_reasons();
      for (vector<string>::iterator itr = reject_reasons.begin(); itr != reject_reasons.end(); itr ++) 
      {  
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
  
// helper function

string string_to_fixed_digit_string(string s, uint32_t precision = 2)
{
  if (s == "NA")
    return "NA";
  double value = from_string<double>(s);
  stringstream ss;
  ss << fixed << setprecision(precision) << value;
  return ss.str();
}
  
string html_new_junction_table_string(diff_entry_list_t& list_ref, bool show_details, const string& title, const string& relative_link)
{
  stringstream ss; //!<< Main Build Object for Function
  cDiffEntry& test_item = *list_ref.front();

  bool link = (test_item.entry_exists(_SIDE_1_EVIDENCE_FILE_NAME) &&
               test_item.entry_exists(_SIDE_2_EVIDENCE_FILE_NAME) &&
               test_item.entry_exists(_NEW_JUNCTION_EVIDENCE_FILE_NAME));

  ss << start_table("border=\"0\" cellspacing=\"1\" cellpadding=\"3\"") << endl;
  size_t total_cols = link ? 11 : 9;
  
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
        th("reads&nbsp;(cov)") << endl <<
        th("reads&nbsp;(cov)") << endl <<
        th("score")       << endl <<
        th("skew")        << endl <<
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
  
  for (diff_entry_list_t::iterator itr = list_ref.begin(); itr != list_ref.end(); itr ++) 
  {  
    cDiffEntry& c = **itr;
// #     ##############
// #     ### Side 1 ###
// #     ##############
    ss << "<!-- Side 1 Item Lines for New Junction -->" << endl;
    
    string key = "side_1";
    string annotate_key = "junction_" + c[key + "_annotate_key"];
    ss << start_tr("class=\"mutation_table_row_" + to_string(row_bg_color_index) +"\"") << endl;

     if (link) {
      ss << td("rowspan=\"2\"", 
              a(relative_link + c[_NEW_JUNCTION_EVIDENCE_FILE_NAME], "*" )) << endl;
     }

     { // Begin hiding data for side 1

      if (link) {   
        ss << td("rowspan=\"1\"", 
                a(relative_link + c[_SIDE_1_EVIDENCE_FILE_NAME], "?")) << endl;
      }
      ss << td("rowspan=\"1\" class=\"" + annotate_key + "\"",
            nonbreaking(c[key + "_seq_id"])) << endl;
       
      if (from_string<int32_t>(c[key + "_strand"]) == 1) { 
        ss << td("align=\"center\" class=\"" + annotate_key +"\"",
                c[key + "_position"] + "&nbsp;=");
      } else {
        ss << td("align=\"center\" class=\"" + annotate_key +"\"",
                "=&nbsp;" + c[key + "_position"]);
      }
      
      ss << td("align=\"center\" class=\"" + annotate_key +"\"",
                c[key + "_read_count"] + " (" + string_to_fixed_digit_string(c[key + "_coverage"], 3) + ")" );
      
      
      //no longer print overlap
      //ss << td("rowspan=\"2\" align=\"center\"", c["overlap"]) << endl;
      ss << td("rowspan=\"2\" align=\"center\"", 
               c["new_junction_read_count"] + " (" + string_to_fixed_digit_string(c["new_junction_coverage"], 3) + ")" ) << endl;
      ss << td("rowspan=\"2\" align=\"center\"", 
               c["pos_hash_score"] + "/" +  c["max_pos_hash_score"]) << endl;
      ss << td("rowspan=\"2\" align=\"center\"", 
               c["neg_log10_pos_hash_p_value"]) << endl;
                
                
               //" (" + c["max_left"] + "/" + c["max_right"] + ")") << endl;
      ss << td("align=\"center\" class=\"" + annotate_key + "\"", 
              nonbreaking(c["_" + key + GENE_POSITION])) << endl;
      ss << td("align=\"center\" class=\"" + annotate_key + "\"", 
              i(nonbreaking(c["_" + key + GENE_NAME]))) << endl;
      ss << td("class=\"" + annotate_key + "\"",
              htmlize(c["_" + key + GENE_PRODUCT])) << endl;  
    } // End hiding data for side 1
    ss << end_tr() << endl;


// #     ##############
// #     ### Side 2 ###
// #     ##############
    ss << "<!-- Side 2 Item Lines for New Junction -->" << endl;
    key = "side_2";
    annotate_key = "junction_" + c[key + "_annotate_key"];
    ss << start_tr("class=\"mutation_table_row_" + 
                  to_string(row_bg_color_index) + "\"") << endl;

    { //Begin hiding data for side 2
      if (link) {
        ss << td("rowspan=\"1\"", 
                a(relative_link +  c[_SIDE_2_EVIDENCE_FILE_NAME], "?"));
      }
      ss << td("rowspan=\"1\" class=\"" + annotate_key + "\"",
              nonbreaking(c[key + "_seq_id"])) << endl;

      if (from_string<int32_t>(c[key + "_strand"]) == 1) { 
        ss << td("align=\"center\" class=\"" + annotate_key +"\"",
                c[key + "_position"] + "&nbsp;=") << endl;
      } else {
        ss << td("align=\"center\" class=\"" + annotate_key +"\"",
                "=&nbsp;" + c[key + "_position"]) << endl;
      } 
      
      ss << td("align=\"center\" class=\"" + annotate_key +"\"",
               c[key + "_read_count"] + " (" + string_to_fixed_digit_string(c[key + "_coverage"], 3) + ")" );
      
      ss << td("align=\"center\" class=\"" + annotate_key + "\"",
              nonbreaking(c["_" + key + GENE_POSITION])) << endl;
      ss << td("align=\"center\" class=\"" + annotate_key + "\"",
              i(nonbreaking(c["_" + key + GENE_NAME]))) << endl;
      ss << td("class=\"" + annotate_key + "\"",
              htmlize(c["_" + key + GENE_PRODUCT])) << endl;
    } //end hiding data for side 2
    
  ss << end_tr() << endl;

  /* Extra debug output
  if (show_details) {
    ss << tr(   td("colspan=\"" + to_string(total_cols) + "\"",
                "Continuation Left: " + c["continuation_left"] + "&nbsp;&nbsp;Continuation Right: " + c["continuation_right"] )) << endl;
  }
  */
    
  if (show_details && c.entry_exists("reject")) {
    cGenomeDiff gd;

    vector<string> reject_reasons = c.get_reject_reasons();
    
    for (vector<string>::iterator itr = reject_reasons.begin();
         itr != reject_reasons.end(); itr++) {
      string& reject(*itr);
    
      ss << tr("class=\"reject_table_row\"",
              td("colspan=\"" + to_string(total_cols) + "\"",
                "Rejected: " + decode_reject_reason(reject))) << endl;
    }

  }
  row_bg_color_index = (row_bg_color_index+1) % 2;//(row_bg_color_index) % 2; 

  }// End list_ref Loop
  ss << "</table>" << endl;
  return ss.str();
}
  
string html_copy_number_table_string(diff_entry_list_t& list_ref, bool show_details, const string& title, const string& relative_link)
{
  stringstream ss; //!<< Main Build Object for Function
  cDiffEntry& test_item = *list_ref.front();
  
  bool link = test_item.entry_exists(_EVIDENCE_FILE_NAME);
  
  ss << start_table("border=\"0\" cellspacing=\"1\" cellpadding=\"3\"") << endl;
  size_t total_cols = link ? 9 : 8;
  
  if (title != "") {
    ss << tr(th("colspan=\"" + to_string(total_cols) + "\" align=\"left\" class=\"copy_number_header_row\"", title)) << endl;
  }
  // #   #####################
  // #   #### HEADER LINE ####
  // #   #####################
  ss << "<tr>" << endl;
  
  if (link) {
    ss << th("&nbsp;") << endl;
  }
  ss << th("seq&nbsp;id") << endl;
  ss << th("start") << endl;
  ss << th("end") << endl;
  ss << th("tile&nbsp;size") << endl;
  ss << th("copy&nbsp;number") << endl;
  ss << th("rel&nbsp;coverage") << endl;
  ss << th("gene") << endl;
  ss << th("width=\"100%\"","product") << endl;
  ss << "</tr>" << endl;
  ss << endl;
  // #   ####################
  // #   #### ITEM LINES ####
  // #   ####################
  
  for (diff_entry_list_t::iterator itr = list_ref.begin(); itr != list_ref.end(); itr ++) 
  {  
    cDiffEntry& c = **itr;
    
    ss << start_tr("class=\"normal_table_row\"") << endl;
    
    if (link)
      ss << td(a(relative_link + c[_EVIDENCE_FILE_NAME], "*" )) << endl;
    
    ss << td(ALIGN_LEFT, nonbreaking(c[SEQ_ID])) << endl;
    ss << td(ALIGN_LEFT, nonbreaking(c[START])) << endl;
    ss << td(ALIGN_LEFT, nonbreaking(c[END])) << endl;

    ss << td(ALIGN_CENTER, nonbreaking(c["tile_size"])) << endl;
    ss << td(ALIGN_CENTER, nonbreaking(c["copy_number"])) << endl;
    
    stringstream num; 
    num << fixed << setprecision(2) << from_string<double>(c["relative_coverage"]);
    ss << td(ALIGN_CENTER, num.str()) << endl;
    
    ss << td(ALIGN_CENTER, i(nonbreaking(c[GENE_NAME])));
    ss << td(ALIGN_LEFT, htmlize(c[GENE_PRODUCT]));
    
    
    if (show_details && c.entry_exists("reject")) {
      cGenomeDiff gd;
      
      vector<string> reject_reasons = c.get_reject_reasons();
      
      for (vector<string>::iterator itr = reject_reasons.begin();
           itr != reject_reasons.end(); itr++) {
        string& reject(*itr);
        
        ss << tr("class=\"reject_table_row\"",
                 td("colspan=\"" + to_string(total_cols) + "\"",
                    "Rejected: " + decode_reject_reason(reject))) << endl;
      }
    }
    ss << end_tr();
  }// End list_ref Loop
  
  ss << end_table() << endl;
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
Evidence_Files::Evidence_Files(const Settings& settings, cGenomeDiff& gd)
{  
  // Fasta and BAM files for making alignments.
  string reference_bam_file_name = settings.reference_bam_file_name;
  string reference_fasta_file_name = settings.reference_fasta_file_name;

  // hybrids use different BAM files for making the alignments!!!
  string junction_bam_file_name = settings.junction_bam_file_name;
  string junction_fasta_file_name = settings.candidate_junction_fasta_file_name;

  
  // We make alignments of two regions for deletions: upstream and downstream edges.
  diff_entry_list_t items_MC = gd.show_list(make_vector<gd_entry_type>(MC));
  //cerr << "Number of MC evidence items: " << items_MC.size() << endl;

  
  for (diff_entry_list_t::iterator itr = items_MC.begin(); itr != items_MC.end(); itr ++) 
  {  
    diff_entry_ptr_t item(*itr);
     
    diff_entry_ptr_t parent_item(gd.parent(*item));
    if (parent_item.get() == NULL)
      parent_item = *itr;

   add_evidence(_SIDE_1_EVIDENCE_FILE_NAME,
                item,
                parent_item,
                make_map<string,string>
                (BAM_PATH, reference_bam_file_name)
                (FASTA_PATH, reference_fasta_file_name)
                (PREFIX, "MC_SIDE_1")
                (SEQ_ID, (*item)[SEQ_ID])            
                (START, to_string(from_string<uint32_t>((*item)[START]) - 1))
                (END,  to_string(from_string<uint32_t>((*item)[START]) - 1)));

    add_evidence(_SIDE_2_EVIDENCE_FILE_NAME,
                item,
                parent_item,
                make_map<string,string>
                (BAM_PATH, reference_bam_file_name)
                (FASTA_PATH, reference_fasta_file_name)
                (PREFIX, "MC_SIDE_2")
                (SEQ_ID, (*item)[SEQ_ID])            
                (START, to_string(from_string<uint32_t>((*item)[END]) + 1))
                (END,  to_string(from_string<uint32_t>((*item)[END]) + 1)));
    
    add_evidence(_EVIDENCE_FILE_NAME,
                item,
                parent_item,
                make_map<string,string>
                (BAM_PATH, reference_bam_file_name)
                (FASTA_PATH, reference_fasta_file_name)
                (PREFIX, "MC_PLOT")
                (SEQ_ID, (*item)[SEQ_ID])            
                (START, (*item)[START])
                (END,  (*item)[END])
                (PLOT, (*item)[_COVERAGE_PLOT_FILE_NAME])); // filled by draw_coverage
    
  } // mc_item list
  
  
  
  diff_entry_list_t items_SNP_INS_DEL_SUB = gd.show_list(make_vector<gd_entry_type>(SNP)(INS)(DEL)(SUB));
  //cerr << "Number of SNP_INS_DEL_SUB evidence items: " << items_SNP_INS_DEL_SUB.size() << endl;

  for (diff_entry_list_t::iterator itr = items_SNP_INS_DEL_SUB.begin(); itr != items_SNP_INS_DEL_SUB.end(); itr ++) 
  {  
    diff_entry_ptr_t item = *itr;
    diff_entry_list_t mutation_evidence_list = gd.mutation_evidence_list(*item);

    // #this reconstructs the proper columns to draw
    uint32_t start = from_string<uint32_t>((*item)[POSITION]);
    uint32_t end = start;
    uint32_t insert_start = 0;
    uint32_t insert_end = 0;

    if (item->_type == INS) 
    {
      insert_start = 1;
      insert_end = (*item)[NEW_SEQ].size();
    }
    else if (item->_type == DEL) 
    {
      bool has_ra_evidence = false;
      for (diff_entry_list_t::iterator itr = mutation_evidence_list.begin(); itr != mutation_evidence_list.end(); itr ++) 
      {  
        cDiffEntry& evidence_item = **itr;
        if (evidence_item._type == RA) has_ra_evidence = true;
      }
      if(!has_ra_evidence) continue;  

      end = start + from_string<uint32_t>((*item)[SIZE]) - 1;
    }

    else if (item->_type == SUB ) 
    {
      end = start + (*item)[NEW_SEQ].size() - 1;
    }

    add_evidence(_EVIDENCE_FILE_NAME,
                 item,
                 item,
                 make_map<string,string>
                 (BAM_PATH, reference_bam_file_name)
                 (FASTA_PATH, reference_fasta_file_name)
                 (SEQ_ID, (*item)[SEQ_ID])            
                 (START, to_string(start))
                 (END,  to_string(end))
                 (INSERT_START, to_string(insert_start))
                 (INSERT_END, to_string(insert_end))
                 (PREFIX, to_string((*item)._type)));


    // Add evidence to RA items as well
    for (diff_entry_list_t::iterator itr = mutation_evidence_list.begin(); itr != mutation_evidence_list.end(); itr ++) 
    {  
      cDiffEntry& evidence_item = **itr;
      if (evidence_item._type != RA) continue;
      evidence_item[_EVIDENCE_FILE_NAME] = (*item)[_EVIDENCE_FILE_NAME];  
    }
  }
  
  
  
  
  // Still create files for RA evidence that was not good enough to predict a mutation from
  diff_entry_list_t items_RA = gd.filter_used_as_evidence(gd.show_list(make_vector<gd_entry_type>(RA)));
  //cerr << "Number of RA evidence items: " << items_RA.size() << endl;

  for (diff_entry_list_t::iterator itr = items_RA.begin(); itr != items_RA.end(); itr ++) 
  {  
    diff_entry_ptr_t item = *itr;

    add_evidence(_EVIDENCE_FILE_NAME,
                 item,
                 item,
                 make_map<string,string>
                 (BAM_PATH, reference_bam_file_name)
                 (FASTA_PATH, reference_fasta_file_name)
                 (SEQ_ID, (*item)[SEQ_ID])            
                 (START, (*item)[POSITION])
                 (END, (*item)[POSITION])
                 (INSERT_START, (*item)[INSERT_POSITION])
                 (INSERT_END, (*item)[INSERT_POSITION])
                 (PREFIX, to_string(item->_type)));
  }
  // This additional information is used for the complex reference line.
  // Note that it is completely determined by the original candidate junction sequence 
  // positions and overlap: alignment_pos and alignment_overlap.
  
  
  diff_entry_list_t items_JC = gd.show_list(make_vector<gd_entry_type>(JC));
  //cerr << "Number of JC evidence items: " << items_JC.size() << endl;

  for (diff_entry_list_t::iterator itr = items_JC.begin(); itr != items_JC.end(); itr ++) 
  {  
    diff_entry_ptr_t item = *itr;

    diff_entry_ptr_t parent_item(gd.parent(*item));
    if (parent_item.get() == NULL)
      parent_item = *itr;

    uint32_t start = 0;
    uint32_t end = 0;

    if (from_string<int32_t>((*item)[ALIGNMENT_OVERLAP]) == 0) 
    {
      start = from_string<uint32_t>((*item)[FLANKING_LEFT]);
      end = from_string<uint32_t>((*item)[FLANKING_LEFT]) + 1;
    }
    else if (from_string <int32_t>((*item)[ALIGNMENT_OVERLAP]) > 0) 
    {
      start = from_string<uint32_t>((*item)[FLANKING_LEFT]) + 1;
      end = from_string<uint32_t>((*item)[FLANKING_LEFT]) + 
            from_string<int32_t>((*item)[ALIGNMENT_OVERLAP]);
    }
    else //if (from_string <uint32_t>((*item)[ALIGNMENT_OVERLAP]) > 0) 
    {
      start = from_string<uint32_t>((*item)[FLANKING_LEFT]) + 1;
      end = from_string<uint32_t>((*item)[FLANKING_LEFT]) - from_string<int32_t>((*item)[ALIGNMENT_OVERLAP]);
    }
    
    // The "key"/ID is set early in breseq.  It must remain unique and unchanging
    // through the run so we know what we're referencing.  Because we derive values
    // from the name (like here), sometimes those values won't match the new resolved
    // values for each side.  This is a dirty fix so that the evidence files
    // will be operating with the correct positions for each side.  Search @MDS0001
    // to find out where we finally access this modified information.
    JunctionInfo juncInfo((*item)["key"]);
    juncInfo.sides[0].redundant = from_string<int32_t>((*item)[SIDE_1_OVERLAP]);
    juncInfo.sides[0].position = from_string<int32_t>((*item)[SIDE_1_POSITION]);
    juncInfo.sides[1].position = from_string<int32_t>((*item)[SIDE_2_POSITION]); 
    
    add_evidence(_NEW_JUNCTION_EVIDENCE_FILE_NAME,
                  item,
                  parent_item,
                  make_map<diff_entry_key_t,diff_entry_value_t>
                 (BAM_PATH, junction_bam_file_name)
                 (FASTA_PATH, junction_fasta_file_name)
                 (SEQ_ID, (*item)["key"])
                 (START, to_string(start))
                 (END, to_string(end))
                 (PREFIX, "JC")
                 (ALIGNMENT_EMPTY_CHANGE_LINE, "1")
                 (CORRECTED_KEY, juncInfo.junction_key())
                 );
    
    // this is the flagship file that we show first when clicking on evidence from a mutation...
    (*item)[_EVIDENCE_FILE_NAME] = (*item)[_NEW_JUNCTION_EVIDENCE_FILE_NAME];
    string side_1_key_str = "JC_SIDE_1_" + (*item)[SIDE_2_SEQ_ID] + "_" + (*item)[SIDE_2_POSITION] + "_" + (*item)[SIDE_2_POSITION];
    add_evidence(_SIDE_1_EVIDENCE_FILE_NAME,
                 item,
                 parent_item,
                 make_map<string,string>
                 (BAM_PATH, reference_bam_file_name)
                 (FASTA_PATH, reference_fasta_file_name)
                 (SEQ_ID, (*item)[SIDE_1_SEQ_ID])            
                 (START, (*item)[SIDE_1_POSITION])
                 (END, (*item)[SIDE_1_POSITION])
                 (PREFIX, side_1_key_str)
                 ); 

    string side_2_key_str = "JC_SIDE_2_" + (*item)[SIDE_1_SEQ_ID] + "_" + (*item)[SIDE_1_POSITION] + "_" + (*item)[SIDE_1_POSITION];
    add_evidence(_SIDE_2_EVIDENCE_FILE_NAME,
                 item,
                 parent_item,
                 make_map<string,string>
                 (BAM_PATH, reference_bam_file_name)
                 (FASTA_PATH, reference_fasta_file_name)
                 (SEQ_ID, (*item)[SIDE_2_SEQ_ID])            
                 (START, (*item)[SIDE_2_POSITION])
                 (END, (*item)[SIDE_2_POSITION])
                 (PREFIX, side_2_key_str)
                 );
  }
  
  // Copy number evidence
  diff_entry_list_t items_CN = gd.filter_used_as_evidence(gd.show_list(make_vector<gd_entry_type>(CN)));
  
  for (diff_entry_list_t::iterator itr = items_CN.begin(); itr != items_CN.end(); itr ++) 
  {  
    diff_entry_ptr_t item = *itr;
    
    add_evidence(_EVIDENCE_FILE_NAME,
                 item,
                 item,
                 make_map<string,string>
                 (BAM_PATH, reference_bam_file_name)
                 (FASTA_PATH, reference_fasta_file_name)
                 (PREFIX, "CN_PLOT")
                 (SEQ_ID, (*item)[SEQ_ID])            
                 (START, (*item)[START])
                 (END,  (*item)[END])
                 (PLOT, (*item)[_COVERAGE_PLOT_FILE_NAME])); // filled by draw_coverage

  }

  
  // now create evidence files
  create_path(settings.evidence_path);
  //cerr << "Total number of evidence items: " << evidence_list.size() << endl;
  
  for (vector<Evidence_Item>::iterator itr = evidence_list.begin(); itr != evidence_list.end(); itr ++) 
  {  
    Evidence_Item& e = (*itr);
    //cerr << "Creating evidence file: " + e[FILE_NAME] << endl;   
    html_evidence_file(settings, gd, e);
  }
}

/*-----------------------------------------------------------------------------
 *  Helper Function For Create_Evidence_Files()
 *-----------------------------------------------------------------------------*/

string Evidence_Files::html_evidence_file_name(Evidence_Item& evidence_item)
{
  
  //set up the file name
  string s = evidence_item[PREFIX];  
  s += "_";
  s += evidence_item[SEQ_ID];  
  s += "_";
  s += evidence_item[START];
  
  if (evidence_item.entry_exists(INSERT_START))
  {
    s += ".";
    s += evidence_item[INSERT_START];
  }
  
  s += "_";
  s += evidence_item[END];
  
  if (evidence_item.entry_exists(INSERT_END))
  {
    s += ".";
    s += evidence_item[INSERT_END];
  }
  
  s += "_alignment.html";  

  return s;
}
  
  
void Evidence_Files::add_evidence(const string& evidence_file_name_key, diff_entry_ptr_t item,
                                  diff_entry_ptr_t parent_item, diff_entry_map_t& fields)
{
  Evidence_Item evidence_item(fields, item, parent_item);
  evidence_item[FILE_NAME] = html_evidence_file_name(evidence_item);
  
  // this is added to the actual genome diff entry so that we know where to link
  (*item)[evidence_file_name_key] = evidence_item[FILE_NAME];
  
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
  if(evidence_item.entry_exists(INSERT_START))
    ss << "." + evidence_item[INSERT_START];
  ss << "_" << evidence_item[END];
  if(evidence_item.entry_exists(INSERT_END))
    ss << "." + evidence_item[INSERT_END];
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
                    cGenomeDiff& gd, 
                    Evidence_Item& item
                   )
{  
  string output_path = settings.evidence_path + "/" + item[FILE_NAME];

  // Create Stream and Confirm It's Open
  ofstream HTML(output_path.c_str());
  
  if (!HTML.good()) {
    cerr << "Could not open file: " << item["output_path"] << endl;
    assert(HTML.good());
  }
  
  // Build HTML Head
  HTML << html_header("BRESEQ :: Evidence", settings);
  
  // print a table for the main item
  // followed by auxiliary tables for each piece of evidence

  diff_entry_ptr_t parent_item = item.parent_item;
  diff_entry_list_t parent_list;
  parent_list.push_back(parent_item);
  HTML << html_genome_diff_item_table_string(settings, gd, parent_list);
  HTML << "<p>";
  
  diff_entry_list_t evidence_list = gd.mutation_evidence_list(*parent_item);

  vector<gd_entry_type> types = make_vector<gd_entry_type>(RA)(MC)(JC);
  
  for (vector<gd_entry_type>::iterator itr = types.begin(); itr != types.end(); itr ++)
  {  
    const gd_entry_type type = *itr;
    diff_entry_list_t this_evidence_list = evidence_list;
    this_evidence_list.remove_if(cDiffEntry::is_not_type(type));   
    
    if(this_evidence_list.empty()) continue;

    HTML << html_genome_diff_item_table_string(settings, gd, this_evidence_list);
    HTML << "<p>"; 
  }
  
  if (item.entry_exists(PLOT) && !item[PLOT].empty())
    HTML << div(ALIGN_CENTER, img(item[PLOT]));
  else
  {
    bool bCorrectExist = false;
    if(item.entry_exists(CORRECTED_KEY))
      bCorrectExist = true;
        
    stringstream ss;
    stringstream sc;
    
    ss << item[SEQ_ID] << ":" << item[START];
    if(bCorrectExist)sc << item[CORRECTED_KEY] << ":" << item[START];
    
    if (item[INSERT_START].size() > 0)
    {
      ss << "." << item[INSERT_START];
      if(bCorrectExist)sc << "." << item[INSERT_START];
    }
    
    ss << "-" << item[END];
    if(bCorrectExist)sc << "-" << item[END];
    
    if (item[INSERT_END].size())
    {
      ss << "." << item[INSERT_END];
      if(bCorrectExist)sc << "." << item[INSERT_END];
    }
    cerr << "Creating read alignment for region: " << ss.str() << endl;

    if (settings.base_quality_cutoff != 0)
      item["base_quality_cutoff"] = to_string(settings.base_quality_cutoff);
    
    alignment_output ao(item[BAM_PATH], item[FASTA_PATH], settings.maximum_reads_to_align, settings.base_quality_cutoff);

    HTML << ao.html_alignment(ss.str(), sc.str());

  }
  HTML << html_footer();
  HTML.close();
}


/*-----------------------------------------------------------------------------
 *  //End Create_Evidence_Files
 *-----------------------------------------------------------------------------*/

void draw_coverage(Settings& settings, cReferenceSequences& ref_seq_info, cGenomeDiff& gd)
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
  
  // Zoom-in plots of individual deletions and copy number variation
  vector<gd_entry_type> mc_types = make_vector<gd_entry_type>(MC)(CN);
	diff_entry_list_t mc = gd.show_list(mc_types);
  for (diff_entry_list_t::iterator it=mc.begin(); it!=mc.end(); it++)
  {
    diff_entry_ptr_t& item = *it;
    uint32_t start = from_string<uint32_t>((*item)[START]);
    uint32_t end = from_string<uint32_t>((*item)[END]);
    uint32_t size = end - start + 1;
    
    uint32_t _shaded_flanking = static_cast<uint32_t>(floor(static_cast<double>(size) / 10.0));
    if (_shaded_flanking < 100) _shaded_flanking = 100;
    co.shaded_flanking(_shaded_flanking);
    
    string region = (*item)[SEQ_ID] + ":" + (*item)[START] + "-" + (*item)[END];
    string coverage_plot_file_name = settings.evidence_path + "/" + (*item)[SEQ_ID] + "_" + (*item)[START] + "-" + (*item)[END] + "." + co.output_format();

    string link_coverage_plot_file_name = Settings::relative_path(coverage_plot_file_name, settings.evidence_path);    
    (*item)[_COVERAGE_PLOT_FILE_NAME] = link_coverage_plot_file_name;
    
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
                                                       cGenomeDiff& gd,
                                                       diff_entry_list_t& list_ref,
                                                       MutationTableOptions& options
                                                       )
  : string()
  , total_cols(0)
  , settings(settings)
  , gd(gd)
  , list_ref(list_ref)
  , options(options)
{
  (*this) += "<!--Output Html_Mutation_Table_String-->\n";
  (*this) += "<table border=\"0\" cellspacing=\"1\" cellpadding=\"3\">\n";
  
  
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
void Html_Mutation_Table_String::Header_Line(bool print_main_header)
{
  // #####################
  // #### HEADER LINE ####
  // #####################
  
  string header_text = ((list_ref.size() > 1) ? "Predicted mutations" : "Predicted mutation");

  stringstream ss; //<! Main Build Object for Function

  // There are three possibilities for the frequency column(s)
  // (1) We don't want it at all. (Single genome no poly prediction)
  vector<string> freq_header_list;

  if (options.gd_name_list_ref.size() > 1) {
    freq_header_list = options.gd_name_list_ref;
  } 
  else if(settings.polymorphism_prediction) {
    freq_header_list = make_vector<string>("freq");
  }

  if (settings.lenski_format) {
    vector<string> header_list = split(freq_header_list.front(), "|");
    size_t header_rows = header_list.size() - 1; //!< -1 is necessary for C++

    total_cols = 7 + freq_header_list.size() ;
    if(!options.one_ref_seq) total_cols += 1; 
    if(!settings.no_evidence) total_cols += 1;

    for (size_t i = 1; i <= header_rows; i++) {
     ss << "<!-- Header Line -->" << endl;
     ss << "<tr>" << endl;
      if(!settings.no_evidence)
        ss << th("evidence") << endl;
      if(!options.one_ref_seq)
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
                 (this_header_string_2 == "ZDB483" ||
                 this_header_string_2 == "ZDB30") )
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
    total_cols = 5 + freq_header_list.size();
    if (!options.one_ref_seq)
      total_cols += 1;

    if (!settings.no_evidence)
      total_cols += 1;


    ss <<  "<tr>" << endl;
    if (!settings.no_evidence)
      ss << th("evidence") << endl;

    if(!options.one_ref_seq)
      ss << th(nonbreaking("seq id")) << endl;

    ss << th("position") << endl;
    ss << th("mutation") << endl;

    if(freq_header_list.size() > 0) {
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

  if(print_main_header && (header_text != ""))
    (*this) += tr(th("colspan=\"" + to_string(total_cols) + "\" align=\"left\" class=\"mutation_header_row\"", header_text));
  
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

  stringstream ss; 
  ss << "<!-- Item Lines -->" << endl;
  (*this) += ss.str();
  ss.str("");
  for (diff_entry_list_t::iterator itr = list_ref.begin(); itr != list_ref.end(); itr ++) { 
    cDiffEntry& mut = (**itr);

    if ((row_num != 0) && (options.repeat_header != 0) && (row_num % options.repeat_header == 0))
    {
      Header_Line(false); // don't print main header again
    }
    row_num++;
        
    // Build Evidence Column
    string evidence_string;
    if (!settings.no_evidence) {
      bool already_added_RA = false;
       
      diff_entry_list_t mutation_evidence_list = gd.mutation_evidence_list(mut);
      
      for (diff_entry_list_t::iterator evitr = mutation_evidence_list.begin(); evitr != mutation_evidence_list.end(); evitr ++) {  
        cDiffEntry& evidence_item = **evitr;

        if (evidence_item._type == RA) {
          if (already_added_RA) 
            continue;
          else 
            already_added_RA = true;
        }
        
        if (!evidence_string.empty()) evidence_string += "&nbsp;";
        evidence_string += a(options.relative_link + evidence_item[_EVIDENCE_FILE_NAME], to_string(evidence_item._type));
      }
    }
  
    string row_class = "normal_table_row";
    if (options.gd_name_list_ref.size() > 1) {
      row_class = "alternate_table_row_" + to_string(row_num % 2);
    }
    
    // There are three possibilities for the frequency column(s)
    // (1) We don't want it at all. (Single genome no poly prediction)   
    vector<string> freq_list;
   
    // (2) We are in compare mode and need a column for each file
    if (options.gd_name_list_ref.size() > 1) {
      //for each gd
      for (uint32_t j = 0; j < options.gd_name_list_ref.size(); j++) {
        string base_name = options.gd_name_list_ref[j];
        string key = "frequency_" + base_name;
        freq_list.push_back(mut.count(key) ? mut[key] : "0");
      }
    }
    
    // (3) We want a single column (polymorphism prediction)
    if (settings.polymorphism_prediction) {
      // polymorphisms get highlighted
      if(mut.entry_exists(FREQUENCY) && (from_string<double>(FREQUENCY) != 1)) {
        row_class = "polymorphism_table_row";
        freq_list.push_back(mut[FREQUENCY]);
      }
      else // frequencies of other entries assumed to be 1.00
      {
        freq_list.push_back("1");
      }
    }
      
    // ### marshal cells defined depending on mutation type
    string cell_seq_id = nonbreaking(mut[SEQ_ID]);
    string cell_position = commify(mut[POSITION]);
    string cell_mutation;
    string cell_mutation_annotation = formatted_mutation_annotation(mut); // Do NOT make nonbreaking
    string cell_gene_name = i(nonbreaking(mut[GENE_NAME]));
    string cell_gene_product = htmlize(mut[GENE_PRODUCT]);
    
    // @MDS0004 - If the product contains more than a set number of genes
    // replace the name with the one that hides it with javascript.
    if(mut.count(GENE_PRODUCT_HIDE))
      if(mut[GENE_PRODUCT_HIDE].size() > 0)
        cell_gene_product = htmlize(mut[GENE_PRODUCT_HIDE]);

    // build 'mutation' column = description of the genetic change    
    switch (mut._type)
    {
      case SNP:{
        cell_mutation = mut["_ref_seq"] + "&rarr;" + mut[NEW_SEQ];
      } break;
        
      case INS:{
        cell_mutation = "+";
        cell_mutation += mut[NEW_SEQ];
      } break;
        
      case DEL:{
        cell_mutation = nonbreaking("&Delta;" + commify(mut["size"]) + " bp");
        string annotation_str;
        
        // special annotation for mediated- and between repeat elements
        if (mut.entry_exists("mediated")) 
          annotation_str = mut["mediated"] + "-mediated"; 
        if (mut.entry_exists("between")) 
          annotation_str = "between " + mut["between"];
        // default
        if(annotation_str.empty()) {
          annotation_str = nonbreaking(mut["gene_position"]);
        } 
        cell_mutation_annotation =  nonbreaking(annotation_str);
      } break;
        
      case SUB:{
        if (options.detailed || (mut["new_seq"].size() <= 4)) {
          cell_mutation = nonbreaking(mut["size"] + " bp&rarr;" + mut["new_seq"]);
        } else {
          cell_mutation = nonbreaking(mut["size"] + " bp&rarr;" + s(mut["new_seq"].size()) + " bp");
        }
      } break;
        
      case CON:{
        cell_mutation = nonbreaking(mut["size"] + " bp&rarr;" + mut["region"]);
      } break;
        
      case MOB:{
        stringstream s;
        
        stringstream s_start;
        if (mut.entry_exists("ins_start")) {
          s_start << "+" << mut["ins_start"];
        }
        if (mut.entry_exists("del_start")) {
          s_start << "&Delta;" << mut["del_start"];
        }
        if (!(s_start.str()).empty()) {
          s << s_start.str() << " :: ";
        }
        
        s << mut["repeat_name"] << " (";
        switch (from_string<int32_t>(mut["strand"]))
        {
          case -1:
            s << "–";
            break;
          case 0:
            s << "?";
            break;
          case +1:
            s << "+";
            break;
        }
        s << ")";
        
        if (from_string<int32_t>(mut["duplication_size"]) > 0) {
          s << " +" << mut["duplication_size"] << " bp";
        } else if (from_string<int32_t>(mut["duplication_size"]) < 0) {
          s << " &Delta;" << abs(from_string(mut["duplication_size"])) << " bp";
        }
        
        stringstream s_end;
        if (mut.entry_exists("del_end")) {
          s_end << " &Delta;" << mut["del_end"];
        }
        if (mut.entry_exists("ins_end")) {
          s_end << " +" << mut["ins_end"];
        }
        if (!(s_end.str()).empty()) {
          s << " ::" << s_end.str();
        }
        
        cell_mutation = nonbreaking(s.str());
      } break;
        
      case INV:{
        cell_mutation = nonbreaking(commify(mut["size"]) + " bp inversion");
        cell_gene_name = i(nonbreaking(mut["gene_name_1"])) + "&darr;" +
        i(nonbreaking(mut["gene_name_2"]));
        cell_gene_product = htmlize(mut["gene_product_1"]) + "&darr;" + 
        htmlize(mut["gene_product_2"]);
      } break;
        
      case AMP:{
        cell_mutation = nonbreaking(commify(mut["size"]) + " bp x " + mut["new_copy_number"]);
        cell_mutation_annotation = 
        from_string<uint32_t>(mut["new_copy_number"]) == 2 ?
        "duplication" : "amplification";
      } break;
        
      default:
        break;
    }
    // ###### PRINT THE TABLE ROW ####
    ss << endl << "<!-- Print The Table Row -->" << endl; 
    ss << start_tr("class=\"" + row_class + "\"") << endl;

    if (!settings.no_evidence) {
      ss << td(ALIGN_CENTER, evidence_string) << "<!-- Evidence -->" << endl;
    }
    if (!options.one_ref_seq) {
      ss << td(ALIGN_CENTER, cell_seq_id) << "<!-- Seq_Id -->" << endl;
    }
    ss << td(ALIGN_RIGHT, cell_position) << "<!-- Position -->" << endl;
    
    ss << td(ALIGN_CENTER, cell_mutation) << "<!-- Cell Mutation -->" << endl;
          
    if (settings.lenski_format) {
      ss << "<!-- Lenski_Format -->" << endl;
      ss << td(ALIGN_CENTER, cell_mutation_annotation) << endl;
      ss << td(ALIGN_CENTER, cell_gene_name) << endl;
      ss << td(ALIGN_LEFT, cell_gene_product) << endl;
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
      ss << td(ALIGN_LEFT, cell_gene_product) << endl;
    }
    ss << end_tr() << endl;

    ss << "<!-- End Table Row -->" << endl;
    
    (*this) += ss.str();
    ss.str("");
    
  } //##### END TABLE ROW ####
  
  if (options.legend_row) {
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
  if (freq == "?")
    return "?";
  
  if (freq == "D")
    return "&Delta;";
  
  if (freq == "H")
    return "H";

  if (from_string<double>(freq) == 0.0)
    return "";

  stringstream ss;
  if (from_string<double>(freq) == 1.0 || freq.empty())
    ss << "100%";

  else {
    double conv_freq = from_string<double>(freq) * 100;
    ss.width(4);
    ss.setf(ios_base::fixed);
    ss.precision(1);
    ss << conv_freq << "%";
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
  stringstream ss;
  for (vector<string>::iterator itr = freq_list.begin();
       itr != freq_list.end(); itr ++) {  
    string& freq = (*itr);
    if (options.shade_frequencies) {
      string bgcolor;
      if (freq == "1") {
       bgcolor = "Blue";
      }
      if (!bgcolor.empty()) {
        ss << td("align=\"right\" bgcolor=\"" + bgcolor +"\"", "&nbsp;");
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

