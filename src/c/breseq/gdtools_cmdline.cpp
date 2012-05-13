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

#include "libbreseq/common.h"
#include "libbreseq/settings.h"
#include "libbreseq/anyoption.h"
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

  uout("Format Conversions:");
  uout << "gd2gvf                 GD to Genome Variant Format(GVF)" << endl;
  uout << "vcf2gd                 Variant Call Format(VCF) to GD" << endl;
  uout << "gd2circos              GD to Circos Data" << endl;

  uout("Set Operations:");
  uout << "subtract               remove mutations" << endl;
  uout << "intersect              locate equal mutations" << endl;
  uout << "union                  combine mutations, removing duplicates" << endl;
  
  return 0;
}

int do_intersection(int argc, char *argv[])
{
  AnyOption options("INTERSECT -o <output.gd> <input1.gd input2.gd ...>");
  options("output,o",  "output GD file name");
  options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("\tCreates a GD file of mutations and associated evidence that exists across all input GD files.");

  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("No output provided.");
    options.printUsage();
    return -1;
  }  

  if (!options.getArgc()) {
    options.addUsage("");
    options.addUsage("No additional files provided.");
    options.printUsage();
    return -1;
  }

  if (options.getArgc() < 2) {
    options.addUsage("");
    options.addUsage("Not enough additional files provided.");
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
  AnyOption options("UNION -o <output.gd> <input1.gd input2.gd ...>");
  options("output,o",  "output GD file name");
  options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("\tCreates a GD file of mutations and associated evidence that exist in each input GD file.");
  options.addUsage("Duplicate mutations are merged into a single mutation.");

  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("No output provided.");
    options.printUsage();
    return -1;
  }  

  if (!options.getArgc()) {
    options.addUsage("");
    options.addUsage("No additional files provided.");
    options.printUsage();
    return -1;
  }

  if (options.getArgc() < 2) {
    options.addUsage("");
    options.addUsage("Not enough additional files provided.");
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
  AnyOption options("APPLY -o <output> -f <fasta/gff3> -r <reference> <input.gd>");
  options("output,o",    "output file");
  options("format,f",    "output file's format (Options: fasta, gff3)");
  options("reference,r", ".gbk/.gff3/.fasta/.bull reference sequence file");
  options("verbose,v",   "Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Input a single GenomeDiff, and as many reference files");
  options.addUsage("as you like.  Using the GenomeDiff, this will apply all");
  options.addUsage("the mutations to the reference sequences, output is to");
  options.addUsage("a single file that includes all the references.");
  
  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("No output provided.");
    options.printUsage();
    return -1;
  }  

  string format;
  if (!options.count("format")) {
    options.addUsage("");
    options.addUsage("No format provided.");
    options.printUsage();
    return -1;
  } else {
    format = to_upper(options["format"]);
  }

  if ((format != "FASTA") && (format != "GFF3")) {
    options.addUsage("");
    options.addUsage("Did not recognize format: " + options["format"]);
    options.printUsage();
    return -1;
  }

  if (!options.getArgc()) {
    options.addUsage("");
    options.addUsage("No additional files provided.");
    options.printUsage();
    return -1;
  }

  UserOutput uout("APPLY");

  cReferenceSequences ref_seq_info;
  cReferenceSequences new_ref_seq_info;

  const vector<string>& ref_paths = from_string<vector<string> >(options["reference"]);
  uout("Reading input reference files", ref_paths);
  ref_seq_info.LoadFiles(ref_paths);
  new_ref_seq_info.LoadFiles(ref_paths);


  uout("Reading input GD file", options.getArgv(0));
  cGenomeDiff gd(options.getArgv(0));


  //Check to see if every item in the loaded .gd is applicable to the reference file.
  ASSERT(gd.is_valid(ref_seq_info, options.count("verbose")), "Reference file and GenomeDiff file don't match.");
  gd.apply_to_sequences(ref_seq_info, new_ref_seq_info, options.count("verbose"));

  uout("Writing output file in " + format + " format");

  if (format == "FASTA") {
    uout << options["fasta"] << endl;
    new_ref_seq_info.WriteFASTA(options["output"], options.count("verbose"));
  }
  else if (format == "GFF3") {
    uout << options["gff3"] << endl;
    new_ref_seq_info.WriteGFF(options["output"], options.count("verbose"));
  }

  return 0;
}

int do_subtract(int argc, char *argv[])
{
  AnyOption options("SUBTRACT [-o output.gd] input.gd minus1.gd [minus2.gd ...]");
  options("output,o",  "output GD file");
  options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Creates a GD file of mutations from the input.gd after removing mutations present in additional minus#.gd GD files.");
  
  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("No output provided.");
    options.printUsage();
    return -1;
  }

  if (options.getArgc() < 2) {
    options.addUsage("");
    options.addUsage("Not enough input files provided.");
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
  AnyOption options("MERGE -o <output.gd> <input1.gd input2.gd ...>");
  options("output,o",     "output GD file");
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
  
  if(argc == 1)  {
    options.printUsage();
    return -1;  }
  
  if (!options.getArgc()) {
    options.addUsage("");
    options.addUsage("No input provided.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("no output provided.");
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
  AnyOption options("WEIGHTS -o <output.gd> <input1.gd input2.gd ...>");
  options("output,o",    "output GD file");
  options("verbose,v",   "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("\tCreates a GD file of mutations and associated evidence that exist in each input GD file.");
  options.addUsage("Duplicate mutations are merged into a single mutation with a 'weight' field added to the");
  options.addUsage("entry. The 'weight' field is the inverse of the frequency of that mutation occuring accross");
  options.addUsage("all input GD files. Unique mutations will therefor have a 'weight' of 1");

  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("No output provided.");
    options.printUsage();
    return -1;
  }

  if (!options.getArgc()) {
    options.addUsage("");
    options.addUsage("No additional input provided.");
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
  AnyOption options("COMPARE -o <output.gd> <control.gd> <test.gd>");
  options("output,o",  "output GD file");
  options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
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

  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("No output provided.");
    options.printUsage();
    return -1;
  }

  if (options.getArgc() < 2) {
    options.addUsage("");
    options.addUsage("Not enough input file provided.");
    options.printUsage();
    return -1;
  }

  UserOutput uout("COMPARE");

  uout("Reading input GD files");
  uout << "Control: " << options.getArgv(0) << endl;
  cGenomeDiff gd1(options.getArgv(0));
  uout << "Test:    " << options.getArgv(1) << endl;
  cGenomeDiff gd2(options.getArgv(1));


  uout("Comparing control vs test GD file");
  gd1.compare(gd2, options.count("verbose"));

  uout("Assigning unique IDs");
  gd1.assign_unique_ids();

  uout("Writing output GD file", options["output"]);
  gd1.write(options["output"]);

  return 0;
}

int do_convert_gvf( int argc, char* argv[]){
  AnyOption options("GD2GVF --input <input.gd> --output <output.gvf>"); 

  options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("input,i","gd file to convert") 
    ("output,o","name of output file")
    ("snv-only", "only output SNV entries", TAKES_NO_ARGUMENT)
    ;
  options.processCommandArgs( argc,argv);

  options.addUsage("");
  options.addUsage("Creates a Genome Variant Format(GVF) file of mutations present in a input GD file.");
  
  if( !options.count("input") || !options.count("output") ){
      options.printUsage(); return -1;
  }
  
  try{
      GDtoGVF( options["input"], options["output"], options.count("snv-only") );
  } 
  catch(...){ 
      return -1; // failed 
  }
  
  return 0;
}

int do_convert_gd( int argc, char* argv[])
{
  AnyOption options("VCF2GD --input <intput.vcf> --output <output.gd>");
  options("input,i","gd file to convert");
  options("output,o","name of output file");
  options.processCommandArgs( argc,argv);

  options.addUsage("");
  options.addUsage("Creates a GD file of mutations present in a input Variant Call Format(VCF) file.");
  
  if(!options.count("input") && !options.count("output")){
    options.printUsage();
    return -1;
  }
  if(!options.count("input")) {
    options.printUsage();
    return -1;
  }
  
  cGenomeDiff gd = cGenomeDiff::from_vcf(options["input"]);
  const string &file_name = options.count("output") ?
        options["output"] :
        cString(options["input"]).remove_ending("vcf") + "gd";

  gd.write(file_name);

  return 0;
}

int do_convert_circos(int argc, char *argv[]){
  
  AnyOption options("GD2CIRCOS -r <reference> [-r <reference2> ...] -o <output_dir> input1.gd [input2.gd ...]");
  
  options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("reference,r", "reference file(s) specified in gd file")
    ("output,o","name of directory to save Circos files")
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
  
  try{
    GDtoCircos(gd_names, 
               from_string<vector<string> >(options["reference"]),
               options["output"]);
  } 
  catch(...){ 
      return -1; // failed 
  }
  
  return 0;
}

int do_convert_mira(int argc, char* argv[]){
  //unsupported before it ever saw the light of day.
  AnyOption options("MIRA2GD --input <input.gd> --output <output_dir>");
  
  options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("input,i","mira feature analysis file to convert") 
    ("output,o","name of gd file to save")
    ;
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Creates a GD file from a MIRA feature analysis file. Be sure to normalize the GD created afterward.");
  
  if (!options.count("input") && !options.count("output")){
    options.printUsage();
    return -1;
  }
  if (!options.count("input")){
    options.printUsage();
    return -1;
  }
  
  try{
      MIRAtoGD( options["input"], options["output"]);
  } 
  catch(...){ 
      return -1; // failed 
  }
  
  return 0;
}

int do_not_evidence(int argc, char *argv[])
{
  AnyOption options("NOT_EVIDENCE -g <genomediff.gd> -o <output.gd>");
  options("genomediff,g","input GD files");
  options("output,o","output GD file");
  options("id,i","Reorder IDs (Flag)", TAKES_NO_ARGUMENT);
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Takes a GenomeDiff file and removes all of the entries that");
  options.addUsage("are NOT used as evidence by a mutation.  Outputs to a new");
  options.addUsage("GenomeDiff file if specified.  If no output is specified,");
  options.addUsage("verbose will still inform what evidence isn't being used.");
  
  if(argc == 1)  {
    options.printUsage();
    return -1;  }
  
  if (!options.count("genomediff")) {
    options.addUsage("");
    options.addUsage("You must supply the --genomediff option for input.");
    options.printUsage();
    return -1;
  }
  
  cGenomeDiff gd1(options["genomediff"]);
  
  gd1.filter_not_used_as_evidence(options.count("verbose"));
  
  if(options.count("output"))
  {
    if(options.count("id"))
    {
      cGenomeDiff gd2;
      gd2.merge(gd1, true, true);
      gd2.write(options["output"]);
    }
    else  {
      gd1.write(options["output"]);  }
  }
  
  return 0;
}

int do_annotate(int argc, char* argv[])
{
  AnyOption options("ANNOTATE [-o annotated.* --html] -r reference.gbk input.1.gd [input.2.gd ... ]");
  
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
    uout("Writing output HTML file", options["output"]);
    gd.write(output_file_name);
  }
  
  return 0;
}

int do_count(int argc, char* argv[])
{
  AnyOption options("COUNT [-o count.csv] -r reference.gbk input.1.gd [input.2.gd ... ]");
  options.addUsage("");
  options.addUsage("Counts the numbers of mutations and other statistics for each input GenomeDiff file.");
  options
  ("help,h", "produce advanced help message", TAKES_NO_ARGUMENT)
  ("output,o", "path to output file with added mutation data.", "count.csv")
  ("reference,r", "reference sequence in GenBank flatfile format (REQUIRED)")
  ("ignore-pseudogenes", "treats pseudogenes as normal genes for calling AA changes", TAKES_NO_ARGUMENT)
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
  
  // Could be a parameter > this is a "large" mutation, <= this is an "small" mutation
  int32_t large_size_cutoff = 50;
    
  vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
  uout("Reading input reference sequence files") << reference_file_names << endl;
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(reference_file_names);
  //ref_seq_info.total_length();
  
  // Need to create base substitution file...
  //my $bsf;
  //my $bs_totals;
  //if ($base_substitution_file)
  //{
  //  $bsf = GenomeDiff::BaseSubstitutionFile->read(input_file => $base_substitution_file, nt_sequence => $ref_seq_info->{ref_strings}->{$main_seq_id});
  //  $bs_totals = $bsf->totals;
  //  print Dumper($bs_totals);
  //}
  
  // Load and annotate genome diffs
  vector<cGenomeDiff> genome_diffs; 

  for (uint32_t i=0; i<gd_path_names.size(); i++) {
    cGenomeDiff gd(gd_path_names[i]);
    uout("Annotating mutations " + gd_path_names[i]);
    ref_seq_info.annotate_mutations(gd, true, options.count("ignore-pseudogenes"));
    genome_diffs.push_back(gd);
  }
  
  // Figure out the names of all "repeat" columns
  map<string,bool> mob_name_hash;
  map<string,bool> con_name_hash;
  
  for (vector<cGenomeDiff>::iterator it=genome_diffs.begin(); it != genome_diffs.end(); ++it) {
    cGenomeDiff &gd = *it;
        
    diff_entry_list_t muts = gd.mutation_list();
    for (diff_entry_list_t::iterator it=muts.begin(); it != muts.end(); ++it) {
      cDiffEntry& mut = **it;
      if (mut._type == MOB) {
        mob_name_hash[mut["repeat_name"]] = true;
      }
      if (mut._type == CON) {
        if (mut.entry_exists("mediated"))
            mob_name_hash[mut["mediated"]] = true;
      }
      if (mut._type == DEL) {
        if (mut.entry_exists("mediated"))
          mob_name_hash[mut["mediated"]] = true;
      }
    }
  }
  vector<string> mob_name_list = map_keys_to_list<string,bool>(mob_name_hash);
  sort(mob_name_list.begin(), mob_name_list.end());
  vector<string> con_name_list = map_keys_to_list<string,bool>(con_name_hash);
  sort(con_name_list.begin(), con_name_list.end());
  
  vector<string> column_headers;
  column_headers.push_back("sample");
  column_headers.push_back("total");
  column_headers.push_back("base_substitution");
  column_headers.push_back("small_indel");
  column_headers.push_back("large_deletion");
  column_headers.push_back("large_insertion");
  column_headers.push_back("large_amplification");
  column_headers.push_back("large_substitution");
  column_headers.push_back("mobile_element_insertion");
  column_headers.push_back("gene_conversion");
  column_headers.push_back("deleted_bp");
  column_headers.push_back("inserted_bp");
  column_headers.push_back("repeat_inserted_bp");
  column_headers.push_back("called_bp");
  column_headers.push_back("total_bp");
  
  vector<string> header_snp_types = prefix_each_in_vector(snp_types, "base_substitution.");
  column_headers.insert(column_headers.end(),header_snp_types.begin(), header_snp_types.end());
  
  vector<string> header_mob_name_list = prefix_each_in_vector(mob_name_list, "mobile_element.");
  column_headers.insert(column_headers.end(),header_mob_name_list.begin(), header_mob_name_list.end());
  
  vector<string> header_con_name_list = prefix_each_in_vector(con_name_list, "gene_conversion.");
  column_headers.insert(column_headers.end(),header_con_name_list.begin(), header_con_name_list.end());
  
  ofstream output_file(output_file_name.c_str());
  ASSERT(output_file.good(), "Error writing to file: " + output_file_name);
  
  uout("Calculating base substitution effects in reference sequences") << endl;
  BaseSubstitutionEffects my_bsf(ref_seq_info);
  
  /*    
   if ($base_substitution_file) {
   
   foreach my $snp_type (@GenomeDiff::BaseSubstitutionFile::bsf_snp_types) {
   foreach my $bp_change (@GenomeDiff::BaseSubstitutionFile::bp_change_label_list, 'TOTAL') {
   push @column_headers, "POSSIBLE.$snp_type\.$bp_change";
   }
   }
   foreach my $snp_type (@GenomeDiff::BaseSubstitutionFile::bsf_snp_types) {
   foreach my $bp_change (@GenomeDiff::BaseSubstitutionFile::bp_change_label_list, 'TOTAL') {
   push @column_headers, "OBSERVED.$snp_type\.$bp_change";
   }
   }
   }
   */
  
  output_file << join(column_headers, ",") << endl;

  /*
   #deep copy sums
   my $this_bs_totals;
   if ($base_substitution_file) {
   $this_bs_totals = dclone($bs_totals);
   }
   */

  for (vector<cGenomeDiff>::iterator it=genome_diffs.begin(); it != genome_diffs.end(); ++it) {
    cGenomeDiff &gd = *it;
    uout("Counting mutations " + gd.metadata.run_name);
    
    // Zero out counts
    int32_t total_deleted = 0;
    int32_t total_inserted = 0;
    int32_t total_repeat_inserted = 0;
    int32_t total_bp = ref_seq_info.total_length();
    
    // Complicated map storing a bunch of counts
    map<string,map<string,int32_t> > count;
    
    for(vector<string>::const_iterator snp_type = snp_types.begin(); snp_type != snp_types.end(); ++snp_type) {
      count["type"][*snp_type] = 0;
    }
    for(vector<string>::const_iterator mob_name = mob_name_list.begin(); mob_name != mob_name_list.end(); ++mob_name) {
      count["mob"][*mob_name] = 0;
    }
    
    count["base_substitution"][""] = 0;
    count["large_deletion"][""] = 0;
    count["small_indel"][""] = 0;
    count["large_insertion"][""] = 0;
    count["large_amplification"][""] = 0;
    count["large_substitution"][""] = 0;
    count["gene_conversion"][""] = 0;
    count["mobile_element_insertion"][""] = 0;
    
    //my $this_bs_counts = [];

    diff_entry_list_t mut_list = gd.mutation_list();
    for (diff_entry_list_t::iterator it=mut_list.begin(); it != mut_list.end(); ++it) {		
      
      cDiffEntry& mut = **it;
      
      if (mut.is_marked_deleted()) continue;
          
      //count SNPs and examples of mobile element insertions
      if (mut._type == SNP) {
        count["base_substitution"][""]++;
        string base_change = mut[REF_BASE] + " " + mut[NEW_BASE];
        //my $bp_change = $GenomeDiff::BaseSubstitutionFile::base_to_bp_change->{$base_change};	
        //if ($base_substitution_file) {
        //  $this_bs_counts = $bsf->add_bp_change_to_totals($this_bs_counts, $mut->{position}, $bp_change);					
        //}
        count["type"][mut["snp_type"]]++;
      }
 
      if (mut._type == DEL) {
        total_deleted += from_string<int32_t>(mut[SIZE]);
        
        if (mut.entry_exists("mediated"))
          count["mob"][mut["mediated"]]++;
        
        if (from_string<int32_t>(mut[SIZE]) > large_size_cutoff)
          count["large_deletion"][""]++;
        else
          count["small_indel"][""]++;
      }
    
      if (mut._type == INS) {
        int32_t ins_size = mut[NEW_SEQ].size();

        total_inserted += ins_size;

        if (mut.entry_exists("mediated"))
          count["mob"][mut["mediated"]]++;
        
        if (ins_size > large_size_cutoff)
          count["large_insertion"][""]++;
        else
          count["small_indel"][""]++;
      }
      
      if (mut._type == SUB) {
        int32_t old_size = from_string<int32_t>(mut[SIZE]);
        int32_t new_size = mut[NEW_SEQ].size();
        
        if (new_size - old_size > 0)
          total_inserted += new_size - old_size;
        else
          total_deleted += old_size - new_size;
          
        if (mut.entry_exists("mediated"))
          count["mob"][mut["mediated"]]++;
        
        if (abs(new_size - old_size) > large_size_cutoff)
          count["large_substitution"][""]++;
        else
          count["small_indel"][""]++;
      }

      if (mut._type == CON) {
        // TODO: Need size change?
        
        if (mut.entry_exists("mediated"))
          count["con"][mut["mediated"]]++;
        count["gene_conversion"][""]++;
      }
   
      if (mut._type == MOB) {
        int32_t rpos = -1;
        string repeat_seq = ref_seq_info.repeat_family_sequence(mut["repeat_name"], 1, rpos);
        
        int32_t this_length = repeat_seq.size();
        total_inserted += this_length;
        total_repeat_inserted += this_length;
        //print "Repeat $mut->{repeat_name} $this_length bp\n";
        
        count["mob"][mut["repeat_name"]]++;
        count["mobile_element_insertion"][""]++;

      } 
      
      if (mut._type == AMP) {
        int32_t this_size = mut.mutation_size_change(ref_seq_info);
        total_inserted += this_size;
        
        if (this_size > large_size_cutoff)
          count["large_amplification"][""]++;
        else
          count["small_indel"][""]++;
      }
    }
      
    // statistics for UN
    int32_t un_bp = 0;
    
    diff_entry_list_t un_list = gd.list(make_vector<gd_entry_type>(UN));
    for (diff_entry_list_t::iterator it=un_list.begin(); it!= un_list.end(); ++it) {
      cDiffEntry& un = **it;
      un_bp += from_string<int32_t>(un[END]) - from_string<int32_t>(un[START]) + 1;

      /*
      if ($base_substitution_file)
      {
        my $ref_string = $ref_seq_info->{ref_strings}->{$main_seq_id};
        
        for (my $p = $un->{start}; $p <= $un->{end}; $p++)
        {
          ### subroutine to subtract certain position info
          $bsf->subtract_position_1_from_totals($this_bs_totals, $p);					
        }
      }
       */
      
    }
    
    int32_t called_bp = total_bp - un_bp;
 
    vector<string> this_columns;
 
    this_columns.push_back(gd.metadata.run_name);
    this_columns.push_back(to_string(mut_list.size()));
    this_columns.push_back(to_string(count["base_substitution"][""]));
    this_columns.push_back(to_string(count["small_indel"][""]));
    this_columns.push_back(to_string(count["large_deletion"][""]));
    this_columns.push_back(to_string(count["large_insertion"][""]));
    this_columns.push_back(to_string(count["large_amplification"][""]));
    this_columns.push_back(to_string(count["large_substitution"][""]));
    this_columns.push_back(to_string(count["mobile_element_insertion"][""]));
    this_columns.push_back(to_string(count["gene_conversion"][""]));
    this_columns.push_back(to_string(total_deleted));
    this_columns.push_back(to_string(total_inserted));
    this_columns.push_back(to_string(total_repeat_inserted));
    this_columns.push_back(to_string(called_bp));
    this_columns.push_back(to_string(total_bp));
    
    vector<string> snp_type_counts = map_key_list_to_values_as_strings(count["type"], snp_types);
    this_columns.insert(this_columns.end(),snp_type_counts.begin(), snp_type_counts.end());
    
    vector<string> mob_type_counts = map_key_list_to_values_as_strings(count["mob"], mob_name_list);
    this_columns.insert(this_columns.end(),mob_type_counts.begin(), mob_type_counts.end());
    
    vector<string> con_type_counts = map_key_list_to_values_as_strings(count["con"], con_name_list);
    this_columns.insert(this_columns.end(),con_type_counts.begin(), con_type_counts.end());
 
   /*
   if ($base_substitution_file) {			
   for (my $i = 0; $i < scalar @$this_bs_totals; $i++) {
   foreach my $bp_change (@GenomeDiff::BaseSubstitutionFile::bp_change_label_list, 'TOTAL') {
   push @this_columns, (defined $this_bs_totals->[$i]->{$bp_change}) ? $this_bs_totals->[$i]->{$bp_change} : '0';
   }
   }
   
   
   for (my $i = 0; $i < scalar @$this_bs_counts; $i++) {
   foreach my $bp_change (@GenomeDiff::BaseSubstitutionFile::bp_change_label_list, 'TOTAL') {
   push @this_columns, (defined $this_bs_counts->[$i]->{$bp_change}) ? $this_bs_counts->[$i]->{$bp_change} : '0';
   }
   }
   }
    */
    output_file << join(this_columns, ",") << endl;
  }	
  
  return 0;
}

int do_normalize_gd(int argc, char* argv[])
{
  AnyOption options("NORMALIZE -o <output.gd> -r <reference> <input1.gd input2.gd input3.gd ...>");
  options
  ("output,o"     , "output GD file.")
  ("reference,r"  , "input reference file.")
  ("verbose,v"    , "Verbose Mode (Flag)", TAKES_NO_ARGUMENT);

  options.processCommandArgs(argc, argv);
  options.addUsage("");
  options.addUsage("Creates a GD file of mutations that have been normalized to the input reference files.");

  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("No output provided.");
    options.printUsage();
    return -1;
  }

  if (!options.count("reference")) {
    options.addUsage("");
    options.addUsage("No reference provided.");
    options.printUsage();
    return -1;
  }

  if (!options.getArgc()) {
    options.addUsage("");
    options.addUsage("No additional input provided.");
    options.printUsage();
    return -1;
  }

  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("You must supply the --output option for output.");
    options.printUsage();
    return -1;
  }

  bool verbose = options.count("verbose");

  UserOutput uout("NORMALIZE");

  vector<string> rfns = from_string<vector<string> >(options["reference"]);
  uout("Reading in reference sequence files") << rfns << endl;
  cReferenceSequences rs;
  rs.LoadFiles(rfns);

  uout("Merging input GD files") <<  options.getArgv(0) << endl;
  cGenomeDiff gd(options.getArgv(0));

  for (int32_t i = 1; i < options.getArgc(); ++i) {
    uout << options.getArgv(i) << endl;
    cGenomeDiff temp_gd(options.getArgv(i));
    gd.merge(temp_gd, false, false, options.count("verbose"));
  }

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
  AnyOption options("FILTER -i <input.gd> -o <output.gd> -m <mutation type> <filters>");
  options("input,i", "Genome diff file to be filtered.");
  options("output,o", "Filter genome diff file.");
  options("mut_type,m", "Only consider this mutation type for filtering.");
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("Creates a GD file of mutations that evaluate as false to the input filtering expressions.");

  if (!options.count("input")) {
    options.addUsage("");
    options.addUsage("You must supply an input genome diff file to filter with option -g");
    options.printUsage();
    return -1;
  }

  if (!options.getArgc()) {
    options.addUsage("");
    options.addUsage("You must supply filter arguments.");
    options.printUsage();
    return -1;
  }

  // Use original file as output if argument is not passed.
  if (!options.count("output")) {
    options["output"] = options["input"];
  }

  UserOutput uout("FILTER");

  list<string> filters;
  for(int32_t i = 0; i < options.getArgc(); ++i) {
    filters.push_back(options.getArgv(i));
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
    for (list<string>:: const_iterator jt = filters.begin(); jt != filters.end(); ++jt) {
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
    return do_convert_gvf(argc_new, argv_new);
  } else if (command == "VCF2GD") {             
    return do_convert_gd( argc_new, argv_new);
  } else if (command == "GD2CIRCOS"){
    return do_convert_circos(argc_new, argv_new);
  } else if(command == "MIRA2GD"){
    return do_convert_mira(argc_new, argv_new);
  }

  return 0;

}

