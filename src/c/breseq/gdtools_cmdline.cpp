/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2012 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#include "libbreseq/common.h"
#include "libbreseq/settings.h"
#include "libbreseq/anyoption.h"
#include "libbreseq/mutation_predictor.h"
#include "libbreseq/genome_diff.h"
#include "libbreseq/reference_sequence.h"
#include "libbreseq/output.h"

using namespace breseq;
using namespace output;

int gdtools_usage()
{
  UserOutput uout(0);
  uout("Usage: gdtools [COMMAND] [OPTIONS]");
  uout("Manipulate Genome Diff (*.gd) files using the following commands.");

  uout("General:");
  uout << "apply                  apply mutations to a sequence" << endl;
  uout << "compare                compare control versus test mutations" << endl;
  uout << "not-evidence           remove evidence not used by any mutations" << endl;
  uout << "annotate               get additional mutation information from sequence" << endl;
  uout << "normalize              normalize mutations to a sequence" << endl;
  uout << "filter                 remove mutations given filtering expressions" << endl;
  uout << "merge                  combine multiple GD files" << endl;
  uout << "header                 create or add header entries" << endl;

  uout("Format Conversions:");
  uout << "gd2gvf                 GD to Genome Variant Format(GVF)" << endl;
  uout << "vcf2gd                 Variant Call Format(VCF) to GD" << endl;
  uout << "gd2circos              GD to Circos Data" << endl;

  uout("Set Operations:");
  uout << "subtract               remove mutations" << endl;
  uout << "intersect              locate equal mutations" << endl;
  uout << "union                  combine mutations, removing duplicates" << endl;

  uout("TACC:");
  uout << "download               download reference and read files given appropriate GD header info" << endl;
  uout << "runfile                create a commands file and launcher script for use on TACC" << endl;
  
  return 0;
}

int do_intersection(int argc, char *argv[])
{
  AnyOption options("gdtools INTERSECT [-o output.gd] input1.gd input2.gd ...");
  options("output,o",  "Output Genome Diff file name", "output.gd");
  options("verbose,v", "Verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("Creates a new Genome Diff file with mutations that are present in ALL input Genome Diff files.");

  if (options.getArgc() < 2) {
    options.addUsage("");
    options.addUsage("Must provide at least two input Genome Diff files.");
    options.printUsage();
    return -1;
  }

  UserOutput uout("INTERSECT");

  uout("Reading input GD files") << options.getArgv(0) << endl;
  cGenomeDiff gd1(options.getArgv(0));

  for(int32_t i = 1; i < options.getArgc(); ++i) {
    uout << options.getArgv(i) << endl;
    cGenomeDiff gd2(options.getArgv(i));
    gd1.set_intersect(gd2, options.count("verbose"));
  }

  uout("Writing output GD file", options["output"]);
  gd1.write(options["output"]);

  return 0;
}

int do_union(int argc, char *argv[])
{
  AnyOption options("gdtools UNION [-o output.gd] input1.gd input2.gd ...");
  options("output,o",  "output GD file name", "output.gd");
  options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("\tCreates a GD file of mutations and associated evidence that exist in each input GD file.");
  options.addUsage("Duplicate mutations are merged into a single mutation.");

  if (options.getArgc() < 2) {
    options.addUsage("");
    options.addUsage("Must provide at least two input Genome Diff files.");
    options.printUsage();
    return -1;
  }
  UserOutput uout("UNION");

  uout("Reading input GD files") << options.getArgv(0) << endl;
  cGenomeDiff gd1(options.getArgv(0));

  for(int32_t i = 1; i < options.getArgc(); ++i) {
    uout << options.getArgv(0) << endl;
    cGenomeDiff gd2(options.getArgv(i));
    gd1.set_union(gd2, options.count("verbose"));
  }

  uout("Assigning unique IDs");
  gd1.assign_unique_ids();

  uout("Writing output GD file", options["output"]); 
  gd1.write(options["output"]);

  return 0;
}

int do_apply(int argc, char *argv[])
{
  AnyOption options("gdtools APPLY [ -o output.gff3 -f GFF3 ] -r reference.gbk [-r reference2.gbk ] input.gd");
  options("output,o",    "output file Default (output.*)");
  options("format,f",    "output file's format (Options: FASTA, GFF3)", "FASTA");
  options("reference,r", "Genbank, GFF3, or FASTA reference sequence file");
  options("verbose,v",   "Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Input a single Genome Diff, and as many reference files");
  options.addUsage("as you like.  Using the Genome Diff, this will apply all");
  options.addUsage("the mutations to the reference sequences, output is to");
  options.addUsage("a single file that includes all the references in the.");
  options.addUsage("requested format.");

  string format = to_upper(options["format"]);

  if ((format != "FASTA") && (format != "GFF3")) {
    options.addUsage("");
    options.addUsage("Did not recognize format: " + options["format"]);
    options.printUsage();
    return -1;
  }
  
  string output;
  if (options.count("output")) {
    output = options["output"];
  }
  else {
    output = "output." + to_lower(format);
  }  

  if (!options.getArgc()) {
    options.addUsage("");
    options.addUsage("No additional files provided.");
    options.printUsage();
    return -1;
  }

  UserOutput uout("APPLY");

  cReferenceSequences ref_seq_info;

  const vector<string>& ref_paths = from_string<vector<string> >(options["reference"]);
  uout("Reading input reference files", ref_paths);
  ref_seq_info.LoadFiles(ref_paths);
  cReferenceSequences new_ref_seq_info = cReferenceSequences::deep_copy(ref_seq_info);

  uout("Reading input GD file", options.getArgv(0));
  cGenomeDiff gd(options.getArgv(0));


  //Check to see if every item in the loaded .gd is applicable to the reference file.
  ASSERT(gd.is_valid(ref_seq_info, options.count("verbose")), "Reference file and GenomeDiff file don't match.");
  gd.apply_to_sequences(ref_seq_info, new_ref_seq_info, options.count("verbose"));

  uout("Writing output file in " + format + " format");

  if (format == "FASTA") {
    uout << options["fasta"] << endl;
    new_ref_seq_info.WriteFASTA(output, options.count("verbose"));
  }
  else if (format == "GFF3") {
    uout << options["gff3"] << endl;
    new_ref_seq_info.WriteGFF(output, options.count("verbose"));
  }

  return 0;
}

int do_subtract(int argc, char *argv[])
{
  AnyOption options("gdtools SUBTRACT [-o output.gd] input.gd subtract1.gd [subtract2.gd ...]");
  options("output,o",  "output GD file", "output.gd");
  options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Creates a new Genome Diff file of mutations from the input file that are present after removing mutations present in any of the subtracted Genome Diff files.");

  if (options.getArgc() < 2) {
    options.addUsage("");
    options.addUsage("Input Genome Diff and at least one Genome Diff to subtract not provided.");
    options.printUsage();
    return -1;
  }

  const bool verbose = options.count("verbose");

  UserOutput uout("SUBTRACT");

  uout("Starting with input GD file", options.getArgv(0));
  cGenomeDiff gd1(options.getArgv(0));

  uout("Subtracting mutations from files");
  for (int32_t i = 1; i < options.getArgc(); ++i) {
    uout << options.getArgv(i) << endl;
    cGenomeDiff gd2(options.getArgv(i));
    gd1.set_subtract(gd2, verbose);
  }
  
  uout("Writing output GD file", options["output"]);
  gd1.write(options["output"]);
  
  return 0;
}

int do_merge(int argc, char *argv[])
{
  AnyOption options("gdtools MERGE [-o output.gd] input1.gd input2.gd ...");
  options("output,o",     "output GD file", "output.gd");
  options("unique,u",     "unique entries only", TAKES_NO_ARGUMENT);  
  options("id,i",         "reorder IDs", TAKES_NO_ARGUMENT);
  options("verbose,v",    "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Input as many GenomeDiff files as you want, and have them");
  options.addUsage("merged together into a single GenomeDiff file specified.");
  options.addUsage("Unique IDs will remain unique across files, any IDs that");
  options.addUsage("aren't unique will get new ones.  Mutations that use those");
  options.addUsage("IDs will be properly updated to acknowledge the change.");
  
  if (options.getArgc() < 2) {
    options.addUsage("");
    options.addUsage("At least two input Genome Diff files must be provided.");
    options.printUsage();
    return -1;
  }

  const bool verbose = options.count("verbose");

  UserOutput uout("MERGE");
  
  uout("Merging input GD files") << options.getArgv(0) << endl;
  cGenomeDiff gd1(options.getArgv(0));

  //Load all the GD files that were input.
  for(int32_t i = 1; i < options.getArgc(); i++)
  {
    uout << options.getArgv(i) << endl;;
    cGenomeDiff gd2(options.getArgv(i));
    gd1.merge(gd2, options.count("unique"), options.count("id"), options.count("verbose"));
  }
  
  uout("Writing output GD file", options["output"]);

  gd1.write(options["output"]);
  
  return 0;
}

int do_weights(int argc, char* argv[])
{
  AnyOption options("gdtools WEIGHTS [-o output.gd input1.gd input2.gd ...");
  options("output,o",    "output GD file", "output.gd");
  options("verbose,v",   "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("\tCreates a GD file of mutations and associated evidence that exist in each input GD file.");
  options.addUsage("Duplicate mutations are merged into a single mutation with a 'weight' field added to the");
  options.addUsage("entry. The 'weight' field is the inverse of the frequency of that mutation occuring accross");
  options.addUsage("all input GD files. Unique mutations will therefor have a 'weight' of 1");

  if (options.getArgc() < 2) {
    options.addUsage("");
    options.addUsage("At least two input Genome Diff files must be provided.");
    options.printUsage();
    return -1;
  }


  const bool verbose = options.count("verbose");

  UserOutput uout("WEIGHTS");

  uout("Readings input GD files"); 
  uout << options.getArgv(0) << endl;
  cGenomeDiff gd1(options.getArgv(0));

  //Load all the GD files that were input.
  for(int32_t i = 1; i < options.getArgc(); i++)
  {
    uout << options.getArgv(i) << endl;
    cGenomeDiff gd2(options.getArgv(i));
    gd1.fast_merge(gd2);
  }


  uout("Calculating weights");
  //Store entries as keys and pointers to matching entries.
  diff_entry_list_t muts = gd1.mutation_list();

  bool (*comp_fn) (const diff_entry_ptr_t&, const diff_entry_ptr_t&) = diff_entry_ptr_sort;
  typedef map<diff_entry_ptr_t, vector<diff_entry_ptr_t>, bool(*)(const diff_entry_ptr_t&, const diff_entry_ptr_t&)> diff_entry_map_t;

  diff_entry_map_t weights_table(comp_fn);

  for (diff_entry_list_t::iterator it = muts.begin(); 
       it != muts.end(); ++it) {
    weights_table[*it].push_back(*it);
  }

  //Given the number of pointers, determine frequency and add 'weight' field to entries.
  for (diff_entry_map_t::iterator it = weights_table.begin();
       it != weights_table.end(); ++it) {
    float weight = 1.0f / static_cast<float>(it->second.size());
    uout << weight << '\t' << it->first->to_spec() << endl;
    
    for (uint32_t i  = 0; i < it->second.size(); ++i) {
      (*it->second[i])["weight"] = to_string(weight);
    }

    if (verbose) {
      cout << "\tEntry: " << *it->first << endl;
      cout << "\tMatches with:" << endl;
      for (uint32_t i  = 0; i < it->second.size(); ++i) {
        cout << "\t\t[entry]: " << (*it->second[i]) << endl;
      }
      cout << endl;
      cout << endl;
    }

  }
  uout("Assigning unique IDs");
  gd1.assign_unique_ids();
  
  uout("Writing output GD file", options["output"]);
  gd1.write(options["output"]);

  return 0;
}


int do_compare(int argc, char *argv[])
{
  AnyOption options("gdtools COMPARE [-o output.gd] control.gd test.gd");
  options("output,o",     "output GD file", "comp.gd");
  options("reference,r",  "reference sequence file");
  options("evidence",     "compare evidence", TAKES_NO_ARGUMENT);
  options("jc-buffer",    "length of sequence segment to compare for JC evidence", 50);
  options("plot-jc",      "plot JC Precision versus Score, argument is a prefix for the file paths");
  options("verbose,v",    "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("   Given a control input GD file with valid mutations, compare each additional");
  options.addUsage("test GD file's mutations and determine if the mutations are true-postivie,");
  options.addUsage("false-negative or false-positive, where:");
  options.addUsage("");
  options.addUsage("       TP, true-positive  : mutation exists in both control and test GD files.");
  options.addUsage("       FN, false-negative : mutation exists in the control GD file but not in the test GD files.");
  options.addUsage("       FP, false-positive : mutation exists in one of the test GD files but not in the control GD file.");
  options.addUsage("");
  options.addUsage("The control and test GD files are merged into the given output GD file, for each");
  options.addUsage("mutation a 'compare' field will be added with a TP, FN, FP value. ");

  if (options.getArgc() != 2) {
    options.addUsage("");
    options.addUsage("Exactly two input Genome Diff files must be provided.");
    options.printUsage();
    return -1;
  }

  if ((options.count("evidence") || options.count("plot-jc")) && !options.count("reference")) {
    options.addUsage("");
    options.addUsage("Reference file is needed to compare evidence.");
    options.printUsage();
    return -1;
  }


  UserOutput uout("COMPARE");

  uout("Reading input GD files");
  uout << "Control: " << options.getArgv(0) << endl;
  cGenomeDiff ctrl(options.getArgv(0));
  uout << "Test:    " << options.getArgv(1) << endl;
  cGenomeDiff test(options.getArgv(1));


  uout("Comparing control vs test GD file");
  cGenomeDiff comp;
  if (options.count("evidence") || options.count("plot-jc")) {

    cReferenceSequences ref;
    ref.LoadFiles(from_string<vector<string> >(options["reference"]));
    cGenomeDiff* gds[2] = {&ctrl, &test};
    for (uint32_t i = 0; i < 2; ++i) {
      diff_entry_list_t muts = gds[i]->mutation_list();
      diff_entry_list_t evid = gds[i]->evidence_list();

      if (muts.size() && evid.empty()) {
        uout("Converting mutations to evidence for file: " + gds[i]->file_name());
        gds[i]->mutations_to_evidence(ref);
      }

    }

    uout("Comparing evidence");
    comp = cGenomeDiff::compare_evidence(ref, un(options["jc-buffer"]), ctrl, test, options.count("verbose"));

    if (options.count("plot-jc")) {
      uout("Evaluating results to plot data.");
      string prefix = options["plot-jc"];
      if (prefix.rfind('/') != string::npos) {
        size_t pos = prefix.rfind('/');
        string dir_path = prefix.substr(0, pos);
        uout << "Creating directory: " + dir_path << endl;
        create_path(dir_path);
      }

      string table_path = prefix + ".table.txt";
      uout << "Creating table: " + table_path << endl;
      cGenomeDiff::write_jc_score_table(comp, table_path, options.count("verbose"));
      string cv_exe = comp.metadata.author == "tophat" ? "tophat" : "breseq";
      uint32_t cutoff = cv_exe == "tophat" ? 0 : 3;

      string plot_jc_score_script_name = "/plot_jc_scores.r";
      string plot_jc_score_script_path = DATADIR + plot_jc_score_script_name;

      uout << "Creating plots: " + prefix + ".preciscion.png and " << prefix + ".sensitivity.png" << endl;;
      string cmd = plot_jc_score_script_path + " " + table_path + " " + prefix + " " + s(cutoff) + " " + cv_exe;
      SYSTEM(cmd, true, false, true);
    }

  } else {
    uout("Comparing mutations");
    comp = cGenomeDiff::compare(ctrl, test, options.count("verbose"));
  }

  uout("Assigning unique IDs to Genome Diff entries");
  comp.assign_unique_ids();

  uout("Writing output GD file", options["output"]);
  comp.write(options["output"]);

  return 0;
}

int do_gd2gvf( int argc, char* argv[]){
  AnyOption options("gdtools GD2GVF [-o output.gvf] input.gd"); 

  options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("output,o","name of output file", "output.gvf")
    ("input,i","gd file to convert") 
    ("snv-only", "only output SNV entries", TAKES_NO_ARGUMENT)
    ;
  options.processCommandArgs( argc,argv);

  options.addUsage("");
  options.addUsage("Creates a Genome Variant Format (GVF) file of mutations present in a input GD file.");
  options.addUsage("GVF files are a type of GFF3 file for describing genomic changes.");
  
  if( options.getArgc() != 1 ){
    options.addUsage("");
    options.addUsage("Provide a single input Genome Diff file.");
    options.printUsage();
    return -1;
  }
  
  try{
      GDtoGVF( options.getArgv(0), options["output"], options.count("snv-only") );
  } 
  catch(...){ 
      return -1; // failed 
  }
  
  return 0;
}

int do_vcf2gd( int argc, char* argv[])
{
  AnyOption options("VCF2GD [-o output.gd] intput.vcf");
  options("output,o","name of output Genome Diff file", "output.gd");
  options.processCommandArgs( argc,argv);

  options.addUsage("");
  options.addUsage("Creates a GD file of mutations present in a input Variant Call Format (VCF) file.");
  
  if( options.getArgc() != 1 ){
    options.addUsage("");
    options.addUsage("Provide a single input VCF file.");
    options.printUsage();
    return -1;
  }
  string input = options.getArgv(0);
  
  cGenomeDiff gd = cGenomeDiff::from_vcf(input);
  const string &file_name = options.count("output") ?
        options["output"] :
        cString(input).remove_ending("vcf") + "gd";

  gd.write(file_name);

  return 0;
}

int do_gd2circos(int argc, char *argv[]){
  
  AnyOption options("gdtools GD2CIRCOS -r <reference> [-r <reference2> ...] -o <output_dir> input1.gd [input2.gd ...]");
  
  options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("reference,r", "reference file(s) specified in gd file")
    ("output,o", "name of directory to save Circos files")
    ("distance,d", "the distance from the center the first axis will be in proportion to the default size", "1.0")
    ("feature,f", "the scale of the features in proportion to the default size", "1.0")
    ;
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Creates text files which can be read by Circos to create a visual representation of a GD file. Executes Circos and saves images in the output directory.");
  
  if (!options.count("output")){
    options.addUsage("");
    options.addUsage("No output provided.");
    options.printUsage();
    return -1;
  }
  if (!options.count("reference")){
    options.addUsage("");
    options.addUsage("No reference provided.");
    options.printUsage();
    return -1;
  }
  
  double distance_scale;
  double feature_scale;
  
  distance_scale = from_string<double>(options["distance"]);
  feature_scale = from_string<double>(options["feature"]);
  
  vector<string> gd_names;
  for (int32_t i = 0; i < options.getArgc(); i++){
    gd_names.push_back(options.getArgv(i));
  }
  
  if (gd_names.size() == 0){
    options.addUsage("");
    options.addUsage("No input provided.");
    options.printUsage();
    return -1;
  }
  
  GDtoCircos(gd_names, 
             from_string<vector<string> >(options["reference"]),
             options["output"],
             distance_scale,
             feature_scale);
  
  return 0;
}

int do_mira2gd(int argc, char* argv[]){
  //unsupported before it ever saw the light of day.
  AnyOption options("gdtools MIRA2GD [-o output.gd] input.mira");
  
  options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("output,o", "name of gd file to save", "output.gd")
    ("input,i", "mira feature analysis file to convert") 
    ;
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Creates a GD file from a MIRA feature analysis file. Be sure to normalize the GD created afterward.");
  
  if( options.getArgc() != 1 ){
    options.addUsage("");
    options.addUsage("Provide a single input MIRA file.");
    options.printUsage();
    return -1;
  }
  string input = options.getArgv(0);
  
  try{
      MIRAtoGD( input, options["output"]);
  } 
  catch(...){ 
      return -1; // failed 
  }
  
  return 0;
}

int do_not_evidence(int argc, char *argv[])
{
  AnyOption options("NOT_EVIDENCE [-o output.gd] input.gd");
  options("output,o","output GD file", "output.gd");
  options("id,i","Reorder IDs (Flag)", TAKES_NO_ARGUMENT);
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Takes a GenomeDiff file and removes all of the entries that");
  options.addUsage("are NOT used as evidence by a mutation.  Outputs to a new");
  options.addUsage("GenomeDiff file if specified.  If no output is specified,");
  options.addUsage("verbose will still inform what evidence isn't being used.");
  
  if( options.getArgc() != 1 ){
    options.addUsage("");
    options.addUsage("Provide a single input Genome Diff file.");
    options.printUsage();
    return -1;
  }
  string input = options.getArgv(0);
  
  cGenomeDiff gd1(input);
  
  gd1.filter_not_used_as_evidence(options.count("verbose"));
  
  if(options.count("id"))
  {
    cGenomeDiff gd2;
    gd2.merge(gd1, true, true);
    gd2.write(options["output"]);
  }
  else  {
    gd1.write(options["output"]);  
  }
  
  return 0;
}

int do_annotate(int argc, char* argv[])
{
  AnyOption options("gdtools ANNOTATE [-o annotated.* --html] -r reference.gbk input.1.gd [input.2.gd ... ]");
  
  options
  ("help,h", "produce advanced help message", TAKES_NO_ARGUMENT)
  ("output,o", "path to output file with added mutation data. (DEFAULT: annotated.gd OR annotated.html")
  ("reference,r", "reference sequence in GenBank flatfile format (REQUIRED)")
  ("ignore-pseudogenes", "treats pseudogenes as normal genes for calling AA changes", TAKES_NO_ARGUMENT)
  ("html", "generate HTML output instead of GenomeDiff (REQUIRED)", TAKES_NO_ARGUMENT)
  ;
  options.addUsage("");
  options.addUsage("If multiple GenomeDiff input files are provided, then they are");
  options.addUsage("merged and the frequencies from each file shown for each mutation.");

  options.processCommandArgs(argc, argv);
  
  UserOutput uout("ANNOTATE");
  
  bool html_output_mode = options.count("html");
  string output_file_name;
  if (options.count("output")) 
    output_file_name = options["output"];
  else
    output_file_name = (html_output_mode) ? "annotated.html" : "annotated.gd";
  
  vector<string> gd_path_names;
  for (int32_t i = 0; i < options.getArgc(); ++i) {
    gd_path_names.push_back(options.getArgv(i));
  }
  
  if ( (gd_path_names.size() == 0) 
    || !options.count("reference")){
    options.printUsage();
    return -1;
  }
  
  // more than one file was provided as input
  bool compare_mode = (gd_path_names.size() > 1);
  
  vector<string> gd_base_names;
  for (uint32_t i = 0; i < gd_path_names.size(); i++){    
    cString s(gd_path_names[i]);
    s = s.get_base_name_no_extension();
    gd_base_names.push_back(s);
  }
  
  // First use merge to produce a file with a line for each mutation
  cGenomeDiff gd;
  vector<diff_entry_list_t> mut_lists;
  vector<cGenomeDiff> gd_list;

  for (uint32_t i = 0; i < gd_path_names.size(); i++){
    uout("Reading input GD file",gd_path_names[i]);
    cGenomeDiff single_gd(gd_path_names[i]);
    //remove the evidence to speed up the merge
    //single_gd.remove(cGenomeDiff::EVIDENCE);
    //No, can't remove UN evidence!
    gd_list.push_back(single_gd);
    mut_lists.push_back(single_gd.mutation_list());
    gd.merge(single_gd);
  }
  gd.sort();
  
  // Then add frequency columns for all genome diffs
  if (compare_mode) {
    
    diff_entry_list_t de_list = gd.mutation_list();
    bool found = false;
    
    for (diff_entry_list_t::iterator it = de_list.begin(); it != de_list.end(); it++) { 

      diff_entry_ptr_t& this_mut = *it;
    
      // for each genome diff compared
      for (uint32_t i=0; i<mut_lists.size(); i++) { 
		
        string freq_key = "frequency_" + gd_base_names[i];
        (*this_mut)[freq_key] = "0";

        diff_entry_list_t& mut_list = mut_lists[i];
        if (mut_list.size() == 0) 
          continue; 
        
        bool found = false;
        
        // for top mutation in this genomedff (they are sorted by position)
        diff_entry_ptr_t check_mut;
        check_mut = mut_list.front();        
        
        // we found the exact same mutation
        if ( (check_mut.get() != NULL) && (*check_mut == *this_mut) ) {
          
          if (check_mut->count(FREQUENCY))
            (*this_mut)[freq_key] = (*check_mut)[FREQUENCY];
          else
            (*this_mut)[freq_key] = "1";
          
          // remove the item
          mut_list.pop_front();
          continue;
        }
          
        if (gd_list[i].mutation_deleted(*this_mut)) {
          (*this_mut)[freq_key] = "D";
          continue;
        }
        
        if (gd_list[i].mutation_unknown(*this_mut)) {
          (*this_mut)[freq_key] = "?";
          continue;
        }
      }
    }
  }
  

  vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
  uout("Reading input reference sequence files") << reference_file_names << endl;
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(reference_file_names);
    
  uout("Annotating mutations");
  ref_seq_info.annotate_mutations(gd, true, options.count("ignore-pseudogenes"));
    
  if (html_output_mode) {
    
    uout("Writing output HTML file", output_file_name);
    
    Settings settings;
    // No evidence needs to be transferred to options and initialized correctly within breseq
    settings.no_evidence = true;
    
    MutationTableOptions mt_options;
    if (compare_mode)
      mt_options.repeat_header = true;
    mt_options.one_ref_seq = ref_seq_info.size() == 1;
    mt_options.gd_name_list_ref = gd_base_names;
    mt_options.repeat_header = 10;
    
    html_compare(settings, output_file_name, "Mutation Comparison", gd, mt_options);
        
  } else {
    uout("Writing output Genome Diff file", options["output"]);
    gd.write(output_file_name);
  }
  
  return 0;
}

int do_count(int argc, char* argv[])
{
  AnyOption options("gdtools COUNT [-o count.csv] -r reference.gbk input.1.gd [input.2.gd ... ]");
  options.addUsage("");
  options.addUsage("Counts the numbers of mutations and other statistics for each input GenomeDiff file.");
  options
  ("help,h", "produce advanced help message", TAKES_NO_ARGUMENT)
  ("output,o", "path to output file with added mutation data.", "count.csv")
  ("reference,r", "reference sequence in GenBank flatfile format (REQUIRED)")
  ("ignore-pseudogenes", "treats pseudogenes as normal genes for calling AA changes", TAKES_NO_ARGUMENT)
  ("base-substitution-statistics,b", "calculate detailed base substitution statistics", TAKES_NO_ARGUMENT)
  
  ;
  options.processCommandArgs(argc, argv);
  
  UserOutput uout("COUNT");
  
  string output_file_name = options["output"];
  
  vector<string> gd_path_names;
  for (int32_t i = 0; i < options.getArgc(); ++i) {
    gd_path_names.push_back(options.getArgv(i));
  }
  
  if ( (gd_path_names.size() == 0) 
      || !options.count("reference")){
    options.printUsage();
    return -1;
  }
    
  vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
  uout("Reading input reference sequence files") << reference_file_names << endl;
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(reference_file_names);
  
  // Load and annotate genome diffs
  vector<cGenomeDiff> genome_diffs; 

  for (uint32_t i=0; i<gd_path_names.size(); i++) {
    cGenomeDiff gd(gd_path_names[i]);
    uout("Annotating mutations " + gd_path_names[i]);
    ref_seq_info.annotate_mutations(gd, true, options.count("ignore-pseudogenes"));
    genome_diffs.push_back(gd);
  }
  
  MutationCountFile(
                    ref_seq_info, 
                    genome_diffs, 
                    output_file_name, 
                    options.count("base-substitution-statistics") 
                    );
  
  return 0;
}

int do_normalize_gd(int argc, char* argv[])
{
  AnyOption options("gdtools NORMALIZE [-o output.gd] -r reference.gbk input.gd");
  options
  ("output,o"     , "output Genome Diff file.", "output.gd")
  ("reference,r"  , "input reference file.")
  ("verbose,v"    , "Verbose Mode (Flag)", TAKES_NO_ARGUMENT);

  options.processCommandArgs(argc, argv);
  options.addUsage("");
  options.addUsage("Creates a GD file of mutations that have been normalized to the input reference files.");

  if (!options.count("reference")) {
    options.addUsage("");
    options.addUsage("No reference provided.");
    options.printUsage();
    return -1;
  }

  if( options.getArgc() != 1 ){
    options.addUsage("");
    options.addUsage("Provide a single input Genome Diff file.");
    options.printUsage();
    return -1;
  }
  string input = options.getArgv(0);

  bool verbose = options.count("verbose");

  UserOutput uout("NORMALIZE");

  vector<string> rfns = from_string<vector<string> >(options["reference"]);
  uout("Reading in reference sequence files") << rfns << endl;
  cReferenceSequences rs;
  rs.LoadFiles(rfns);

  cGenomeDiff gd(input);
  
  uout("Normalizing mutations");
  gd.normalize_to_sequence(rs);
  diff_entry_list_t muts = gd.mutation_list();

  for (diff_entry_list_t::iterator it = muts.begin();
       it != muts.end(); it++) {
    diff_entry_ptr_t mut = *it;
    if (mut->entry_exists("norm") && (*mut)["norm"] == "is_not_valid"){
      uout << "INVALID_MUTATION\t" << **it << endl; 
      (*mut)["comment_out"] = "True";
    } 
  }

  uout("Writing output GD file", options["output"]);
  gd.write(options["output"]);

  return 0;
}

int do_filter_gd(int argc, char* argv[]) {
  AnyOption options("gdtools FILTER  [-o output.gd] -f filter1 [-f filter2] [-m SNP] input.gd");
  options("output,o", "Output Genome Diff file.", "output.gd");
  options("mut_type,m", "Only consider this mutation type for filtering.");
  options("filter,f", "Filters to apply when selecting Genome Diff entries.");
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("Creates a GD file of mutations that evaluate as false to the input filtering expressions.");

  UserOutput uout("FILTER");
  
  if (options.getArgc() != 1) {
    options.addUsage("");
    options.addUsage("Provide a single Genome Diff input file.");
    options.printUsage();
    return -1;
  }

  vector<string> filters = from_string<vector<string> >(options["reference"]);
  if (!options.count("filter")) {
    options.addUsage("");
    options.addUsage("You must supply filter arguments.");
    options.printUsage();
    return -1;
  }
  
  assert(filters.size());

  cGenomeDiff gd(options["input"]);

  diff_entry_list_t muts;
  if (!options.count("mut_type")) {
    muts = gd.mutation_list();
  }
  else {
    vector<string> mut_strs= from_string<vector<string> >(options["mut_type"]);
    assert (mut_strs.size());
    vector<gd_entry_type> mut_types;
    for (size_t i = 0; i < mut_strs.size(); ++i) {
      for (size_t j = 0; j < gd_entry_type_lookup_table.size(); ++j) {
        if (gd_entry_type_lookup_table[j] == mut_strs[i]) {
          mut_types.push_back((gd_entry_type)j);
          break;
        }
      }
    }
    muts = gd.list(mut_types);
  }
  if (!muts.size()) {
    printf("No mutations found.\n");
    return 0;
  }

  cGenomeDiff output_gd;
  const vector<string> evals = make_vector<string>("==")("!=")("<=")(">=")("<")(">");
  for (diff_entry_list_t::iterator it = muts.begin(); it != muts.end(); ++it) {
    cDiffEntry mut = **it;

    vector<string> reasons;
    for (vector<string>:: const_iterator jt = filters.begin(); jt != filters.end(); ++jt) {
      bool is_filtered = false;
      // Parse filter string.
      string filter = *jt;
      for (size_t i = 0; i < filter.size(); ++i) {
        if (filter[i] == ' ') {
          filter.erase(i,1);
        }
      }
      string key   = "";
      size_t eval  = string::npos;
      string value = "";
      for (size_t i = 0; i < evals.size(); ++i) {
        if (filter.find(evals[i]) != string::npos) {
          size_t pos = filter.find(evals[i]);
          eval = i;
          string segment = filter;
          segment.erase(pos, evals[i].size());
          key = segment.substr(0, pos);
          value = segment.substr(pos);
          break;
        }
      }
      assert(eval != string::npos);
      assert(key.size());
      assert(value.size());

      if (mut.count(key)) {
        if (value.find_first_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ") == string::npos) {
          switch(eval)
          {
            case 0: if (from_string<float>(mut[key]) == from_string<float>(value)) is_filtered = true; break;
            case 1: if (from_string<float>(mut[key]) != from_string<float>(value)) is_filtered = true; break;
            case 2: if (from_string<float>(mut[key]) <= from_string<float>(value)) is_filtered = true; break;
            case 3: if (from_string<float>(mut[key]) >= from_string<float>(value)) is_filtered = true; break;
            case 4: if (from_string<float>(mut[key]) <  from_string<float>(value)) is_filtered = true; break;
            case 5: if (from_string<float>(mut[key]) >  from_string<float>(value)) is_filtered = true; break;
          };

          if (is_filtered) {
            reasons.resize(reasons.size() + 1);
            sprintf(reasons.back(), "%s %s %s", key.c_str(), evals[eval].c_str(), value.c_str());
          }
        }
      }
    }
    if (reasons.size()) {
      printf("Filtered[%s]: %s\n", join(reasons, ", ").c_str(), mut.to_string().c_str());
      mut["filtered"] = join(reasons, ", ");
      mut["comment_out"] = "true";
    }
    output_gd.add(mut);
  }
  uout("Writing output GD file", options["output"]);
  output_gd.write(options["output"]);
  return 0;
}

int do_rand_muts(int argc, char *argv[])
{
  AnyOption options("Usage: gdtools RANDOM-MUTATIONS -r <reference> -o <output.gd> -t <type>");  
  options("reference,r","Reference file");  
  options("output,o","Output file");
  options("type,t","Type of mutation to generate");
  options("exclude,e","Exclusion file");
  options("number,n","Number of mutations to generate", static_cast<uint32_t>(1000));
  options("buffer,b","Buffer distance between mutations", static_cast<uint32_t>(50));
  options("seq,s","Sequence to use from reference");  
  options("seed","Seed for the random number generator");
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Using -reference, this command will generate a --number of");
  options.addUsage("mutations, and space them based on the supplied --length.");  
  options.addUsage("Not supplying --seq will use the first sequence in the reference.");
  options.addUsage("");
  options.addUsage("Required fields are -r, -o, and -t.");
  options.addUsage("Valid types: SNP, INS, DEL, MOB, AMP, RMD");
  options.addUsage("INS:1-10 will generate insertions of size 1 to 10.");
  options.addUsage("DEL:1-10 will generate deletions of size 1 to 10.");
  
  if(argc <= 1)  {
    options.printUsage();
    return -1;  }
  
  if (!options.count("reference") || !file_exists(options["reference"].c_str())) {
    options.addUsage("");
    options.addUsage("You must supply the --reference option for input.");
    options.addUsage("If you feel you've received this message in error, please");
    options.addUsage("check to see that the file exists.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("You must supply the --output option for output.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("type")) {
    options.addUsage("");
    options.addUsage("You must supply the --type option so that we can choose");
    options.addUsage("which mutation to generate.");
    options.printUsage();
    return -1;
  }
  
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(from_string<vector<string> >(options["reference"]));
  
  if(options.count("exclude"))  {
    options["exclude"];  }
  
  int ref_seq_id = 0;
  if(options.count("seq"))  {
    ref_seq_id = ref_seq_info.seq_id_to_index(options["seq"]);  }
  
  if (options.count("seed")) {
    cSimFastqSequence::SEED_VALUE = un(options["seed"]);
  }

  cGenomeDiff gd;
  gd.add_breseq_data("COMMAND", join(vector<string>(argv, argv + argc), " "));
  gd.random_mutations(
      options["exclude"],
      options["type"],
      un(options["number"]),
      un(options["buffer"]),
      ref_seq_info[ref_seq_id],
      options.count("verbose")
      );

  gd.mutations_to_evidence(ref_seq_info, false);

  gd.assign_unique_ids();
  
  gd.write(options["output"]);
  
  return 0;
}

int do_mutations_to_evidence(int argc, char *argv[])
{
  AnyOption options("Usage: breseq MUTATIONS-TO-EVIDENCE -r <reference> -o <output.gd> input.gd");  
  options("reference,r","Reference file");  
  options("output,o","Output file");
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  if(argc == 1)  {
    options.printUsage();
    return -1;  }
  
  if (!options.count("reference")) {
    options.addUsage("");
    options.addUsage("You must supply the --reference option.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("You must supply the --output option.");
    options.printUsage();
    return -1;
  }
  
  vector<string> gd_path_names;
  for (int32_t i = 0; i < options.getArgc(); ++i) {
    gd_path_names.push_back(options.getArgv(i));
  }
  if (gd_path_names.size() != 1) {
    options.addUsage("");
    options.addUsage("You must supply exactly one input genome diff.");
    options.printUsage();
  }
  
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(from_string<vector<string> >(options["reference"]));
  
  cGenomeDiff gd(gd_path_names[0]);
  gd.mutations_to_evidence(ref_seq_info);
  gd.write(options["output"]);
  
  return 0;
}



int do_header(int argc, char* argv[]) {
  AnyOption options("gdtools HEADER [-o output.gd] [-r reference] file1.fastq file2.fastq ...");
  options("output,o",      "output GD file");
  options("reference,r",   "reference file");
  options("tag,t",         "header tag to add to Genome Diff file, input as <key>=<value> will produce #=<key> <value>");
  options("input,i",       "modify input Genome Diff file inplace");
  options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Create or add '#=<TAG>' entries to the header of a GenomeDiff file,");
  options.addUsage("the -r argument will be added as #=REFSEQ and the *.fastq arguments");
  options.addUsage("will be added as #=READSEQ");

  if (!options.count("reference") && !options.count("tag") && !options.getArgc()) {
    options.printUsage();
    return -1;
  }

  if (!options.count("input") && !options.count("output")) {
    options.addUsage("");
    options.addUsage("No input or output Genome Diff file provided.");
    options.printUsage();
    return -1;
  }

  UserOutput uout("HEADER");

  cGenomeDiff gd;
  if (options.count("input")) {
    uout("Reading in " + options["input"]);
    gd.read(options["input"]);
  }

  if (options.count("reference")) {
    gd.metadata.ref_seqs.clear();
    uout("Adding #=REFSEQ header info for");
    vector<string> refs = from_string<vector<string> >(options["reference"]);
    for (uint32_t i = 0; i < refs.size(); ++i) {
      uout << refs[i] << endl;
      gd.metadata.ref_seqs.push_back(refs[i]);
    }
  }

  if (options.getArgc()) {
    gd.metadata.read_seqs.clear();
    uout("Adding #=READSEQ header info for");
    for (int32_t i = 0; i < options.getArgc(); ++i) {
      uout << options.getArgv(i) << endl;
      gd.metadata.read_seqs.push_back(options.getArgv(i));
    }
  }

  if (options.count("tag")) {
    uout("Adding user defined tags");
    vector<string> tags = from_string<vector<string> >(options["tag"]);
    for (uint32_t i = 0; i < tags.size(); ++i) {
      cKeyValuePair kvp(tags[i], '=');
      uout << "#=" << to_upper(kvp.get_key()) << '\t' << kvp.get_value() << endl;
      gd.add_breseq_data(to_upper(kvp.get_key()), kvp.get_value());
    }
  }

  string output_path = "";
  if (options.count("input") && options.count("output")) {
    output_path = options["output"];
  } 
  else 
  if (options.count("output")) {
    output_path = options["output"];
  }
  else
  if (options.count("input")) {
    output_path = options["input"];
  }

  uout("Writing output GD file", output_path );
  gd.write(output_path);


  return 0;
}

int do_download(int argc, char *argv[])
{
  stringstream ss;
  ss << "Usage: breseq DOWNLOAD -l <user:password> -d <download_dir> -g <genome_diff_dir>\n";
  ss << "Usage: breseq DOWNLOAD -l <user:password> -d <download_dir> <file1.gd file2.gd file3.gd ...>\n";

  AnyOption options(ss.str());
  options("login,l",           "Login user:password information for private server access.");
  options("download-dir,d",    "Output directory to download file to.", "02_Downloads");
  options("genome-diff-dir,g", "Directory to searched for genome diff files.", "01_Data");
  options("test"           ,   "Test urls in genome diff files, doesn't download the file.", TAKES_NO_ARGUMENT);
  options("reference-only",    "Only downloads the reference sequence files for this file.", TAKES_NO_ARGUMENT);

  options.processCommandArgs(argc, argv);

  options.addUsage("\nExamples:");
  options.addUsage("  breseq DOWNLOAD -l john:1234 -d downloads -g data");
  options.addUsage("  breseq DOWNLOAD -l john:1234 -d downloads 1B4.gd GRC2000.gd");

  //! Step: Confirm genome diff files have been input.
  list<string> file_names;
  if (options.getArgc()) {
    const size_t n = options.getArgc();
    for (size_t i = 0; i < n; ++i)
    file_names.push_back(options.getArgv(i));
  } else {
    const string &data_dir = cString(options["genome-diff-dir"]).trim_ends_of('/');
    if (ifstream(data_dir.c_str()).good()) {
      const string cmd = cString("ls %s/*.gd", data_dir.c_str());
      SYSTEM_CAPTURE(back_inserter(file_names), cmd, true);
    }
  }

  if (file_names.empty()) {
    options.addUsage("\nERROR: You must input genome diff files or a directory to search for genome diff files.");
    options.printUsage();
    return -1;
  }

  printf("\n++Starting download.\n");

  string download_dir = cString(options["download-dir"]).trim_ends_of('/');
  if (!options.count("test")) {
    create_path(download_dir);
  }

  //! Step: Create lookup table to store c_style format strings for urls and file paths.
  map<string, map<string, string> > lookup_table;
  //Url formats.
    lookup_table["GENBANK"]
        ["url_format"] = "http://www.ncbi.nlm.nih.gov/sviewer/?db=nuccore&val=%s&report=gbwithparts&retmode=text";
    lookup_table["SRA"]
        ["url_format"] = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s.fastq.gz";
    lookup_table["BARRICKLAB-PUBLIC"]
        ["url_format"] = "http://barricklab.org/%s";
    lookup_table["BARRICKLAB-PRIVATE"]
        ["url_format"] = options.count("login") ? "ftp://" + options["login"] + "@backup.barricklab.org/%s" : "";

  //File path formats.
    lookup_table["GENBANK"]
        ["file_path_format"] = download_dir + "/%s.gbk";
    lookup_table["SRA"]
        ["file_path_format"] = download_dir + "/%s.fastq.gz";
    lookup_table["BARRICKLAB-PUBLIC"]
        ["file_path_format"] = download_dir + "/%s";
    lookup_table["BARRICKLAB-PRIVATE"]
        ["file_path_format"] = download_dir + "/%s";

  /*! Wrather than crash; gather [gd_file_name][reason] = error_value.
      and output at the end of the function. */
  map<string, map<string, string> > error_report;

  for (;file_names.size(); file_names.pop_front()) {
    const string &file_name = file_names.front();
    cout << endl << "Parsing file: " << file_name << endl;
    cGenomeDiff gd(file_name);
    const vector<string> &refs  = gd.metadata.ref_seqs;
    const vector<string> &reads = gd.metadata.read_seqs;

    list<string> seqs_kv_pairs;
    copy(refs.begin(), refs.end(), back_inserter(seqs_kv_pairs));
    if (!options.count("reference-only")) {
      copy(reads.begin(), reads.end(), back_inserter(seqs_kv_pairs));
    }

    for (;seqs_kv_pairs.size(); seqs_kv_pairs.pop_front()) {
      const cKeyValuePair seq_kvp(seqs_kv_pairs.front(), ':');
      if (!seq_kvp.valid()) {
        error_report[file_name]["NOT_KEY_VALUE_PAIR"] = seq_kvp;
        continue;
      }

      const string &key = to_upper(seq_kvp.get_key());
      const string &value = cString(seq_kvp.get_value()).trim_ends_of('/');

      if (!lookup_table.count(key)) {
        error_report[file_name]["INVALID_KEY"] = key;
        continue;
      }

      //! Step: Get file path and check if it has already been downloaded or is empty.
      const string &base_name = cString(value).get_base_name();
      string file_path = "";
      sprintf(file_path, lookup_table[key]["file_path_format"].c_str(), base_name.c_str());
      assert(file_path.size());

      bool is_downloaded =
          ifstream(file_path.c_str()).good() && !file_empty(file_path.c_str());

      bool is_gzip =
          cString(file_path).ends_with(".gz");

      const string &gunzip_path = is_gzip ?
          cString(file_path).remove_ending(".gz") : "";

      if (is_gzip) {
        if (ifstream(gunzip_path.c_str()).good() && !file_empty(gunzip_path.c_str())) {
         is_downloaded = true;
         file_path = gunzip_path;
        }
      }

      if (is_downloaded) {
        printf("File:%s already downloaded in directory %s.\n",
               file_path.c_str(), download_dir.c_str());
        if (!options.count("test")) continue;
      }

      //! Step: Get url and download, gunzip if necessary.
      string url = "";
      const char *url_format = lookup_table[key]["url_format"].c_str();
      if (key == "GENBANK") {
        sprintf(url, url_format, value.c_str());
      }
      else if (key == "SRA") {
        const string &first  = value.substr(0,6), second = value.substr(0,9), third  = value;
        sprintf(url, url_format, first.c_str(), second.c_str(), third.c_str());
      }
      else if (key == "BARRICKLAB-PUBLIC") {
        sprintf(url, url_format, value.c_str());
      }
      else if (key == "BARRICKLAB-PRIVATE") {
        ASSERT(options.count("login"), "Provide the login option (-l user:password) to access private files.");
        sprintf(url, url_format, value.c_str());
      }

      string wget_cmd = options.count("test") ?
            cString("wget --spider \"%s\"", url.c_str()) :
            cString("wget -O %s \"%s\"",file_path.c_str(), url.c_str());

      const bool is_url_error = SYSTEM(wget_cmd, false, false, false) != 0;
      if (is_url_error) {
        error_report[file_name]["INVALID_URL"]  = url;
        continue;
      }

      if (is_gzip && !is_url_error) {
        string gunzip_cmd = "";
        sprintf(gunzip_cmd, "gunzip %s", file_path.c_str());
        if (options.count("test")) {
          cout << gunzip_cmd << endl;
        } else {
          const bool is_gunzip_error = SYSTEM(gunzip_cmd, false, false, false) != 0;
          if (is_gunzip_error) {
            error_report[file_name]["CORRUPT_GZIP_FILE"] = file_path;
          } else {
            file_path = gunzip_path;
          }
        }
      }

      /*! Step: Confirm files are in proper format (mainly want to check that we haven't
      just downloaded an html error page.) */
      if (!options.count("test")) {
        ifstream in(file_path.c_str());
        string first_line = "";
        std::getline(in, first_line);
        if (cString(file_path).ends_with(".gbk")) {
          if (first_line.find("LOCUS") != 0) {
            error_report[file_name]["NOT_GBK_FILE"] = file_path;
          }
        }
        else if(cString(file_path).ends_with(".fastq")) {
          if (first_line[0] != '@') {
            error_report[file_name]["NOT_FASTQ_FILE"] = file_path;
          }
        }
      }
    }//End sequences loop.
  }//End genome diff file_names loop.

  //! Step: Output error_report.
  if (error_report.size()) {
    printf("\n\nERROR: The following problems were encountered:\n");
    map<string, map<string, string> >::const_iterator i = error_report.begin();
    for (;i != error_report.end(); ++i) {
      cout << endl << "Genome diff file: " << i->first << endl;
      map<string, string>::const_iterator j = i->second.begin();
      for (;j != i->second.end(); ++j) {
        printf("\t[%s]: %s\n", j->first.c_str(), j->second.c_str());
      }
    }
  }


  return 0;
}

int do_runfile(int argc, char *argv[])
{
  stringstream ss;
  ss << "Usage: breseq RUNFILE -e <executable> -d <downloads dir> -o <output dir> -l <error log dir> -r <runfile name> -g <genome diff data dir>\n";
  ss << "Usage: breseq RUNFILE -e <executable> -d <downloads dir> -o <output dir> -l <error log dir> -r <runfile name> <file1.gd file2.gd file3.gd ...>";
  AnyOption options(ss.str());
  options("executable,e",     "Executable program to run, add extra options here.", "breseq");
  options("options",          "Options to be passed to the executable.");

  options("data_dir,g",       "Directory to searched for genome diff files.", "01_Data");
  options("downloads_dir,d",  "Downloads directory where read and reference files are located.", "02_Downloads");
  options("output_dir,o",     "Output directory for commands within the runfile.", "03_Output");
  options("log_dir,l",        "Directory for error log file that captures the executable's stdout and sterr.", "04_Logs");
  options("runfile,r",        "Name of the run file to be output.", "commands");
  options("launcher",         "Name of the launcher script run on TACC", "launcher");

  options("job",              "Job name. Alters add suffix to runfile and launcher script.");
  options("email",            "Email address to be added to launcher script. Informs you when a job starts and finishes.");

  options("aln",              "Pass alignment files from breseq's pipeline as an argument.", TAKES_NO_ARGUMENT);

  options.addUsage("\n");
  options.addUsage("***Reminder: Create the error log directory before running TACC job.");
  options.addUsage("\n");
  options.addUsage("Examples:");
  options.addUsage("\tCommand: breseq runfile -o 1B4_Mutated -l 1B4_Mutated_Errors 1B4.gd");
  options.addUsage("\t  Output: breseq -o 1B4_Mutated -r NC_012660.1.gbk SRR172993.fastq >& 1B4_Mutated_Errors/1B4.errors.txt");
  options.addUsage("\n");
  options.addUsage("\tCommand: breseq runfile -d 02_Downloads -l 04_Errors -g 01_Data");
  options.addUsage("\t  Output: breseq -o 1B4 -r 02_Downloads/NC_012660.1.gbk 02_Downloads/SRR172993.fastq >& 04_Errors/1B4.errors.txt");
  options.addUsage("\t  Output: breseq -o ZDB111 -r 02_Downloads/REL606.5.gbk 02_Downloads/SRR098039.fastq >& 04_Errors/ZDB111.errors.txt");
  options.processCommandArgs(argc, argv);

  //! Step: Confirm genome diff files have been input.
  list<string> file_names;
  if (options.getArgc()) {
    const size_t n = options.getArgc();
    for (size_t i = 0; i < n; ++i)
    file_names.push_back(options.getArgv(i));
  } else {
    const string &data_dir = cString(options["data_dir"]).trim_ends_of('/');
    if (ifstream(data_dir.c_str()).good()) {
      const string &cmd = cString("ls %s/*.gd", data_dir.c_str());
      SYSTEM_CAPTURE(back_inserter(file_names), cmd, true);
    }
  }

  if (file_names.empty()) {
    options.addUsage("\nERROR: You must input genome diff files or a directory to search for genome diff files.");
    options.printUsage();
    return -1;
  }

  const string &downloads_dir =
      cString(options["downloads_dir"]).trim_ends_of('/');

  map<string, map<string, string> > lookup_table;
  //Paths where .gbk and .fastq files are located.
  lookup_table["GENBANK"]
      ["download_path_format"] = downloads_dir + "/%s.gbk";
  lookup_table["SRA"]
      ["download_path_format"] = downloads_dir + "/%s.fastq";
  lookup_table["BARRICKLAB-PUBLIC"]
      ["download_path_format"] = downloads_dir + "/%s";
  lookup_table["BARRICKLAB-PRIVATE"]
      ["download_path_format"] = downloads_dir + "/%s";
  lookup_table["LOCAL"]
      ["download_path_format"] = downloads_dir + "/%s";

  const string &exe = options["executable"];

  string job           = "";
  string log_dir       = cString(options["log_dir"]).trim_ends_of('/');
  string output_dir    = cString(options["output_dir"]).trim_ends_of('/');
  string runfile_path  = cString(options["runfile"]).trim_ends_of('/');
  string launcher_path = cString(options["launcher"]).trim_ends_of('/');

  if (options.count("job")) {
    job = options["job"];
    runfile_path  += "_" + job;
    launcher_path += "_" + job;
  }
  create_path(log_dir.c_str());

  const string &log_path_format = log_dir + "/%s.log.txt";

  ofstream runfile(runfile_path.c_str());
  size_t n_cmds = 0;
  for (;file_names.size(); file_names.pop_front()) {
    const string &file_name = file_names.front();
    cout << endl << "Parsing file: " << file_name << endl;
    cGenomeDiff gd(file_name);
    const vector<string> &refs  = gd.metadata.ref_seqs;
    const vector<string> &reads = gd.metadata.read_seqs;

    list<string> seq_kv_pairs;
    copy(refs.begin(), refs.end(), back_inserter(seq_kv_pairs));
    copy(reads.begin(), reads.end(), back_inserter(seq_kv_pairs));
    assert(seq_kv_pairs.size());

    //! Step: Begin building command line.
    stringstream ss;

    //! Part 1: Executable and options to pass to it if given by user.
    ss << exe;

    if (options.count("options")) {
      ss << " " << options["options"];
    }

    //! Part 2: Pipeline's output path.
    ss << " -o " << output_dir + "/" + gd.metadata.run_name;

    size_t n_refs = refs.size();
    for (;seq_kv_pairs.size(); seq_kv_pairs.pop_front()) {
      const cKeyValuePair seq_kvp(seq_kv_pairs.front(), ':');

      cString download_path = "";
        if (seq_kvp.valid()) {

        const string &key = to_upper(seq_kvp.get_key());
        const string &value = cString(seq_kvp.get_value()).trim_ends_of('/');
        const string &base_name = cString(value).get_base_name();

        if (!lookup_table.count(key)) {
          WARN("File: " + file_name + " has invalid key: " + key);
          continue;
        }
        download_path = cString(lookup_table[key]["download_path_format"].c_str(), base_name.c_str());
      } else {
        download_path = seq_kvp;
      }

      //! Part 3: Reference argument path(s).
      if (n_refs) {
        ss << " -r " << download_path;
        n_refs--;
      } else {
      //! Part 4: Read arguement path(s).
        if (!options.count("aln")) {
          if (download_path.ends_with(".gz")) download_path.remove_ending(".gz");
          ss << " " << download_path;
        }
      }
    }
    if (options.count("aln")) {
      ss << " -a " << options["output_dir"] << "/breseq/" << gd.metadata.run_name << "/03_candidate_junctions/best.sam";
    }
    //! Part 5: Error log path.
    ss << " >& " << cString(log_path_format.c_str(), gd.metadata.run_name.c_str());

    //! Step: Output to file.
    cout << ss.str() << endl;
    runfile << ss.str() << endl;
    ++n_cmds;
  }

  //! Step: Create launcher script.
  /*Note: For lonestar we are under the current assumption that a 4way 12 will
    run 3 breseq jobs. On Ranger a 16way 16 will run 16 breseq jobs.*/
  //Ranger:   '/share/home/$NUM/$USER'
  //Lonestar: '/home1/$NUM/$USER'
  size_t tasks = 0, nodes = 0;
  const cString &home_path = SYSTEM_CAPTURE("echo $HOME", true);
  // RANGER
  if (home_path.starts_with("/share")) {
    tasks = 16;
    nodes = static_cast<size_t>(ceilf(static_cast<float>(n_cmds) / 16.f) * 16);
  }
  // Default to LONESTAR
  else {
    if (!home_path.starts_with("/home1/")) {
      WARN("TACC system not determined, defaulting to Lonestar.");
    }
    tasks = 4;
    nodes = static_cast<size_t>(ceilf(static_cast<float>(n_cmds) / 3.f) * 12);
  }
  assert(tasks || nodes);

  const string &pwd = SYSTEM_CAPTURE("pwd", true);

  job = job.size() ? job : "breseq";

  ofstream launcher(launcher_path.c_str());
  // #$ Parameters.
  fprintf(launcher, "#!/bin/csh\n");
  fprintf(launcher, "#$ -N %s\n", job.c_str());
  fprintf(launcher, "#$ -pe %uway %u\n", tasks, nodes);
  fprintf(launcher, "#$ -q normal\n");
  fprintf(launcher, "#$ -o %s.o$JOB_ID\n", job.c_str());
  fprintf(launcher, "#$ -l h_rt=14:00:00\n");
  fprintf(launcher, "#$ -V\n");
  fprintf(launcher, "#$ -cwd\n");

  if (options.count("email")) {
    fprintf(launcher, "#$ -M %s\n", options["email"].c_str());
  }

  fprintf(launcher, "#$ -m be\n");
  fprintf(launcher, "#$ -A breseq\n");
  fprintf(launcher, "\n");

  // Set environmental variables.
  fprintf(launcher, "module load launcher\n");
  fprintf(launcher, "setenv EXECUTABLE     $TACC_LAUNCHER_DIR/init_launcher\n");
  fprintf(launcher, "setenv CONTROL_FILE   %s\n", runfile_path.c_str());
  fprintf(launcher, "setenv WORKDIR        %s\n", pwd.c_str());
  fprintf(launcher, "\n");

  // Job submission.
  fprintf(launcher, "cd $WORKDIR/\n");
  fprintf(launcher, "$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE\n");

  launcher.close();
  SYSTEM("chmod +x " + launcher_path, true);

  return 0;
}

int main(int argc, char* argv[]) {
	//Extract the sub-command argument.
	string command;
	char* argv_new[argc];
	int argc_new = argc - 1;

  //Print out our generic header.
  Settings::command_line_run_header();
                                      
	if (argc > 1) {

		command = argv[1];
		argv_new[0] = argv[0];
		for (int32_t i = 1; i < argc; i++)
			argv_new[i] = argv[i + 1];

	} else {
    gdtools_usage();
    return -1; 
	}

	//Pass the command to the proper handler.
	command = to_upper(command);

  // Genome Diff Commands:
  if (command == "APPLY") {
    return do_apply(argc_new, argv_new);    
  } else if (command == "COMPARE") {
    return do_compare(argc_new, argv_new);
  } else if (command == "NOT-EVIDENCE") {        //TODO merge with FILTER
    return do_not_evidence(argc_new, argv_new);
  } else if (command == "ANNOTATE") {
    return do_annotate(argc_new, argv_new);
  } else if (command == "COUNT") {
    return do_count(argc_new, argv_new);
  } else if (command == "NORMALIZE") {
    return do_normalize_gd(argc_new, argv_new);
  } else if (command == "FILTER") {
    return do_filter_gd(argc_new, argv_new);

  } else if (command == "INTERSECT") {
    return do_intersection(argc_new, argv_new);
  } else if (command == "UNION") {
    return do_union(argc_new, argv_new);
  } else if (command == "SUBTRACT") {
    return do_subtract(argc_new, argv_new);    
  } else if (command == "MERGE") {
    return do_merge(argc_new, argv_new);    
  } else if (command == "WEIGHTS") {
    return do_weights(argc_new, argv_new);
  } else if (command == "GD2GVF") {             
    return do_gd2gvf(argc_new, argv_new);
  } else if (command == "VCF2GD") {             
    return do_vcf2gd( argc_new, argv_new);
  } else if (command == "GD2CIRCOS"){
    return do_gd2circos(argc_new, argv_new);
  } else if(command == "MIRA2GD"){
    return do_mira2gd(argc_new, argv_new);
  } else if(command == "HEADER"){
    return do_header(argc_new, argv_new);
  } else if ((command == "RANDOM-MUTATIONS") || (command == "RAND-MUTS")) {
    return do_rand_muts(argc_new, argv_new);
  } else if (command == "MUTATIONS-TO-EVIDENCE") {
    return do_mutations_to_evidence(argc_new, argv_new);
  } else if (command == "DOWNLOAD") {
    return do_download(argc_new, argv_new);
  } else if (command == "RUNFILE") {
    return do_runfile(argc_new, argv_new);
  } else {
    cout << "Unrecognized command: " << command << endl;
    gdtools_usage();
  }
  return 0;

}

