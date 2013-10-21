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
  uout << "VALIDATE               check formating of input files" << endl;
  uout << "APPLY                  apply mutations to a sequence" << endl;
  uout << "ANNOTATE (or COMPARE)  annotate the effects of mutations and compare multiple samples" << endl;
  uout << "CHECK                  compare control versus test mutations" << endl;
  //uout << "normalize              normalize mutations to a sequence" << endl;
  //uout << "header                 create or add header entries" << endl;
  //uout << "mRNA-Stability         determine mRNA free energy difference of mutations" << endl;

  uout("Set and Filtering Operations:");
  uout << "SUBTRACT               remove mutations in one file from another" << endl;
  uout << "INTERSECT              keep shared mutations in two files" << endl;
  uout << "UNION                  combine mutations, removing duplicates" << endl;
  uout << "MERGE                  combine mutations, preserving duplicates" << endl;
  uout << "FILTER                 remove mutations given filtering expressions" << endl;
  uout << "NOT-EVIDENCE           remove evidence not used by any mutations" << endl;
    
  uout("Creating test data:");
  uout << "SIMULATE-MUTATIONS       create a file containing random mutations" << endl;
    
  uout("Format Conversions:");
  uout << "GD2VCF                 GD to Variant Call Format (VCF)" << endl;
  uout << "GD2GVF                 GD to Genome Variation Format (GVF)" << endl;
  //uout << "vcf2gd                 Variant Call Format(VCF) to GD" << endl;
  uout << "GD2CIRCOS              GD to Circos Data" << endl;

  uout("TACC Utilities:");
  uout << "DOWNLOAD               download reference and read files given appropriate GD header info" << endl;
  uout << "RUNFILE                create a commands file and launcher script for use on TACC" << endl;
  
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
  gd1.reassign_unique_ids();

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

  bool (*comp_fn) (const diff_entry_ptr_t&, const diff_entry_ptr_t&) = cGenomeDiff::diff_entry_ptr_sort;
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
  gd1.reassign_unique_ids();
  
  uout("Writing output GD file", options["output"]);
  gd1.write(options["output"]);

  return 0;
}

int do_check_plot(int argc, char *argv[])
{
    AnyOption options("gdtools CHECK-PLOT [-o output] test1.gd [test2.gd ...]");
    options("output,o",     "Prefix for output files prefix_ ", "output");
    options("verbose,v",    "verbose mode", TAKES_NO_ARGUMENT);
    options.processCommandArgs(argc, argv);
    
    options.addUsage("");
    options.addUsage("Creates sensitivity and precision plots for Genome Diff files output by CHECK, ");
    options.addUsage("which contain TP|FP|FN information and a score value for mutations or evidence.");

    if (options.getArgc() < 1) {
        options.addUsage("");
        options.addUsage("At least one input Genome Diff files must be provided.");
        options.printUsage();
        return -1;
    }
    
    UserOutput uout("CHECK-PLOT");
    
    string gd_file_name = options.getArgv(0);
    uout("Loading Genome Diff file:" + gd_file_name);
    cGenomeDiff merged_gd(gd_file_name);
    merged_gd.metadata.breseq_data.erase("TP|FN|FP");
    
    // Load and quick merge genome diff files (not removing duplicates)
    for (int32_t i=1; i<options.getArgc(); i++) {
        string gd_file_name = options.getArgv(i);
        uout("Loading Genome Diff file:" + gd_file_name);
        cGenomeDiff gd(gd_file_name);
        
        // Merge in a way that preserves duplicates
        merged_gd.fast_merge(gd);
    }
    
    
    string prefix = options["output"];

    // Write out the merged file
    string merged_gd_path = prefix + ".gd";
    uout << "Creating merged Genome Diff: " + merged_gd_path << endl;
    merged_gd.write(merged_gd_path);
    
    string table_path = prefix + ".table.txt";
    uout << "Creating table: " + table_path << endl;
    cGenomeDiff::write_jc_score_table(merged_gd, table_path, options.count("verbose"));
    string cv_exe = merged_gd.metadata.author == "tophat" ? "tophat" : "breseq";
    uint32_t cutoff = cv_exe == "tophat" ? 0 : 3;
    
    string plot_jc_score_script_name = "/plot_jc_scores.r";
    string plot_jc_score_script_path = DATADIR + plot_jc_score_script_name;
    
    uout << "Creating plots: " + prefix + ".precision.png and " << prefix + ".sensitivity.png" << endl;
    string cmd = plot_jc_score_script_path + " " + table_path + " " + prefix + " " + s(cutoff) + " " + cv_exe;
    SYSTEM(cmd, true, false, true);
    
    return 0;
}

// Check format of Genome Diff and report all format errors
int do_validate(int argc, char *argv[])
{
    AnyOption options("gdtools VALIDATE input1.gd [input2.gd]");
    //options("reference,r",  "reference sequence file");
    options.addUsage("");
    options.addUsage("Validates whether the format of the input Genome Diff files is correct.");

    options.processCommandArgs(argc, argv);

    // Simply read files. Really need some way to change deadly errors to warnings in read function.
    for (int32_t i=0; i<options.getArgc(); i++) {
        string gd_file_name = options.getArgv(i);
        cerr <<  "Validating Genome Diff format for file: " << gd_file_name << endl;
        cGenomeDiff gd;
        cFileParseErrors pe = gd.read(gd_file_name, true);
        if (pe._errors.size() == 0) {
            cerr << "  FORMAT OK" << endl;
        } else {
            pe.print_errors(false);
        }
    }    
    return 0;
}

int do_check(int argc, char *argv[])
{
  AnyOption options("gdtools CHECK [-o output.gd] control.gd test.gd");
  options("output,o",     "output GD file", "comp.gd");
  options("reference,r",  "reference sequence file");
  options("evidence",     "compare evidence", TAKES_NO_ARGUMENT);
  options("jc-buffer",    "length of sequence segment to compare for JC evidence", 50);
  options("jc-shorten",   "length to shorten control segments by when comparing JC evidence for overlap", 5);
  options("plot-jc",      "plot JC Precision versus Score, argument is a prefix for the file paths");
  options("verbose,v",    "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("Compare a control input GD file with known valid mutations and a test GD file");
  options.addUsage("");
  options.addUsage("The control and test GD files are merged into the given output GD file. For each");
  options.addUsage("mutation a 'compare' field will be added with a TP, FN, FP value where: ");
  options.addUsage("");
  options.addUsage("       TP, true-positive  : mutation exists in both control and test GD files.");
  options.addUsage("       FN, false-negative : mutation exists in the control GD file but not in the test GD files.");
  options.addUsage("       FP, false-positive : mutation exists in one of the test GD files but not in the control GD file.");
  options.addUsage("");
    
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


  UserOutput uout("CHECK");

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

    uout("Comparing evidence");
    comp = cGenomeDiff::validate_evidence(ref, un(options["jc-buffer"]), un(options["jc-shorten"]), ctrl, test, options.count("verbose"));

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

      uout << "Creating plots: " + prefix + ".precision.png and " << prefix + ".sensitivity.png" << endl;
      string cmd = plot_jc_score_script_path + " " + table_path + " " + prefix + " " + s(cutoff) + " " + cv_exe;
      SYSTEM(cmd, true, false, true);
    }

  } else {
    uout("Comparing mutations");
    comp = cGenomeDiff::validate(ctrl, test, options.count("verbose"));
  }

  uout("Assigning unique IDs to Genome Diff entries");
  comp.reassign_unique_ids();

  uout("Writing output GD file", options["output"]);
  comp.write(options["output"]);

  return 0;
}

int do_gd2vcf( int argc, char* argv[])
{
    AnyOption options("gdtools GD2VCF [-o output.vcf] input.gd"); 
    
    options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("reference,r",  "reference sequence file (REQUIRED)")
    ("output,o","name of output file", "output.vcf")
    ;
    options.processCommandArgs( argc,argv);
    
    options.addUsage("");
    options.addUsage("Creates a Variant Call Format (VCF) file of mutations present in an input Genome Diff file.");
    options.addUsage("VCF is a community format that can be loaded into viewers and used as input to other programs.");
    
    if( options.getArgc() != 1 ){
        options.addUsage("");
        options.addUsage("You must provide exactly one input Genome Diff file.");
        options.printUsage();
        return -1;
    }
    
    if (!options.count("reference")) {
        options.addUsage("");
        options.addUsage("You must provide a reference sequence file (-r).");
        options.printUsage();
        return -1;
    }
        
    cReferenceSequences ref_seq_info;
    ref_seq_info.LoadFiles(from_string<vector<string> >(options["reference"]));
    
    cGenomeDiff::GD2VCF( options.getArgv(0), options["output"], ref_seq_info );
    
    return 0;
}


int do_gd2gvf( int argc, char* argv[])
{
  AnyOption options("gdtools GD2GVF [-o output.gvf] input.gd"); 

  options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("reference,r","Reference file")
    ("output,o","name of output file", "output.gvf")
    ("snv-only","only include SNP/SNV entries in output", TAKES_NO_ARGUMENT)
    ;
  options.processCommandArgs( argc,argv);

  options.addUsage("");
  options.addUsage("Creates a Genome Variation Format (GVF) file of mutations present in an input Genome Diff file.");
  
  if( options.getArgc() != 1 ){
    options.addUsage("");
    options.addUsage("Provide a single input Genome Diff file.");
    options.printUsage();
    return -1;
  }
    
  vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(reference_file_names);
    
  cGenomeDiff::GD2GVF( options.getArgv(0), options["output"], ref_seq_info, options.count("snv-only") );
  
  return 0;
}

int do_vcf2gd( int argc, char* argv[])
{
  AnyOption options("VCF2GD [-o output.gd] input.vcf");
  options("output,o","name of output Genome Diff file", "output.gd");
  options.processCommandArgs( argc,argv);

  options.addUsage("");
  options.addUsage("Creates a GD file of mutations present in an input Variant Call Format (VCF) file.");
  
  if( options.getArgc() != 1 ){
    options.addUsage("");
    options.addUsage("Provide a single input VCF file.");
    options.printUsage();
    return -1;
  }
  
  cGenomeDiff::VCF2GD(options.getArgv(0), options["output"]);
  
  return 0;
}

int do_gd2circos(int argc, char *argv[])
{
  
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
  
  cGenomeDiff::GD2Circos(gd_names, 
             from_string<vector<string> >(options["reference"]),
             options["output"],
             distance_scale,
             feature_scale);
  
  return 0;
}

int do_mira2gd(int argc, char* argv[])
{
  //unsupported before it ever saw the light of day.
  AnyOption options("gdtools MIRA2GD [-o output.gd] input.mira");
  
  options
    ("help,h", "produce this help message", TAKES_NO_ARGUMENT)
    ("output,o", "name of gd file to save", "output.gd")
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
  
  ERROR("Not implemented.");
  //cGenomeDiff::MIRA2GD( input, options["output"]);
  
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
  AnyOption options("gdtools ANNOTATE/COMPARE [-o annotated.html] -r reference.gbk input.1.gd [input.2.gd ... ]");
  
  options
  ("help,h", "produce advanced help message", TAKES_NO_ARGUMENT)
  ("output,o", "path to output file with added mutation data. (DEFAULT: annotated.html OR annotated.gd")
  ("reference,r", "reference sequence in GenBank flatfile format (REQUIRED)")
  ("ignore-pseudogenes", "treats pseudogenes as normal genes for calling AA changes", TAKES_NO_ARGUMENT)
  ("gd", "generate GenomeDiff output instead of HTML", TAKES_NO_ARGUMENT)
  ;
  options.addUsage("");
  options.addUsage("If multiple GenomeDiff input files are provided, then they are");
  options.addUsage("merged and the frequencies from each file are shown for each mutation.");

  options.processCommandArgs(argc, argv);
  
  UserOutput uout("ANNOTATE");
  
  bool html_output_mode = !options.count("gd");
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
      uint32_t this_mut_position = from_string<uint32_t>((*this_mut)[POSITION]);
    
      // for each genome diff compared
      for (uint32_t i=0; i<mut_lists.size(); i++) { 
		
        string freq_key = "frequency_" + gd_base_names[i];
        (*this_mut)[freq_key] = "0";

        diff_entry_list_t& mut_list = mut_lists[i];
        if (mut_list.size() == 0) 
          continue; 
        
        bool found = false;
        

          
        // We have some problems when there are multiple INS after the same position in a genomediff...
        // they get merged to one, e.g.
        // INS	177	1750	T7_WT_Genome	24198	T	frequency=0.1310
        // INS	178	1751	T7_WT_Genome	24198	T	frequency=0.0250
        while (mut_list.size() && (from_string<uint32_t>((*mut_list.front())[POSITION]) < this_mut_position)) {
          cout << (*mut_list.front())[POSITION] << " < " << this_mut_position << endl;
          mut_list.pop_front();
        }
        if (mut_list.size() == 0) 
          continue; 
        // End code for multiple INS problem.  
          
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
  ("verbose,v", "produce output for each mutation counted.", TAKES_NO_ARGUMENT)
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
    ref_seq_info.annotate_mutations(gd, true, options.count("ignore-pseudogenes"), options.count("verbose"));
    genome_diffs.push_back(gd);
  }
  
  MutationCountFile(
                    ref_seq_info, 
                    genome_diffs, 
                    output_file_name, 
                    options.count("base-substitution-statistics"),
                    options.count("verbose")
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

int do_filter_gd(int argc, char* argv[]) 
{
  AnyOption options("gdtools FILTER  [-o output.gd] -f filter1 [-f filter2] [-m SNP] input.gd");
  options("output,o", "Output Genome Diff file.", "output.gd");
  options("mut_type,m", "Only consider this mutation type for filtering.");
  options("filter,f", "Filters to apply when selecting Genome Diff entries. Mutations that match this condition will be removed from the resulting Genome Diff file. Enclose the value of this parameter in quotes, e.g. -f \"frequency<=0.05\".");
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

  if (!options.count("filter")) {
    options.addUsage("");
    options.addUsage("You must supply filter arguments.");
    options.printUsage();
    return -1;
  }
  vector<string> filters = from_string<vector<string> >(options["filter"]);

  assert(filters.size());

  uout("Reading input GD file") << options.getArgv(0) << endl;
  cGenomeDiff gd(options.getArgv(0));

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
    
      ASSERT(key.size() && value.size(), "Error in format of filter: " + filter);
        
      // special case of frequency missing if 100%
      if ( (key=="frequency") && (mut.count(key)==0) ) {
          mut[key] = "1";
      }  
        
        
      if (mut.count(key)) {
          
        // Numeric
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
        }
        // String
        else {
          switch(eval)
          {
            case 0: if (mut[key] == value) is_filtered = true; break;
            case 1: if (mut[key] != value) is_filtered = true; break;
            case 2: if (mut[key] <= value) is_filtered = true; break;
            case 3: if (mut[key] >= value) is_filtered = true; break;
            case 4: if (mut[key] <  value) is_filtered = true; break;
            case 5: if (mut[key] >  value) is_filtered = true; break;
          };
        }
        if (is_filtered) {
          reasons.resize(reasons.size() + 1);
          sprintf(reasons.back(), "%s %s %s", key.c_str(), evals[eval].c_str(), value.c_str());
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

int do_simulate_mutations(int argc, char *argv[])
{
  AnyOption options("Usage: gdtools SIMULATE-MUTATIONS [-n 100] -r <reference> -o <output.gd> -t <type>");  
  options("reference,r","Reference file");  
  options("output,o","Output file");
  options("type,t","Type of mutation to generate");
  options("number,n","Number of mutations to generate", static_cast<uint32_t>(1000));
  options("exclude,e","Exclusion file containing regions no mutations will overlap, usually repeats");
  options("buffer,b","Buffer distance between mutations and excluded intervals", static_cast<uint32_t>(50));
  options("seq,s","Reference sequence id to use from reference file");  
  options("seed","Seed for the random number generator");
  options("verbose,v","Verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Using -reference, this command will generate a --number of");
  options.addUsage("mutations, and space them based on the supplied --buffer.");  
  options.addUsage("Not supplying --seq will use the first sequence in the reference.");
  options.addUsage("");
  options.addUsage("Required fields are -r, -o, and -t.");
  options.addUsage("Valid types: SNP, INS, DEL, MOB, AMP, RMD");
  options.addUsage(" INS:1-10 will generate insertions of size 1 to 10.");
  options.addUsage(" DEL:1-10 will generate deletions of size 1 to 10.");
  options.addUsage(" MOB:1-10 will generate insertions of random repeat regions with");
  options.addUsage("          target site duplications of 1 to 10 bases.");
  options.addUsage(" AMP:100-1000,2-5 will generate tandem amplifications with sizes of");
  options.addUsage("          100-1000 bases and final copy numbers of 2-5.");

  options.addUsage("");

  options.addUsage("An exclusion file of repeats should be generated using a MUMmer command like this:");
  options.addUsage("  mummer -maxmatch -b -c -l 36 REL606.fna REL606.fna > reference.exclude.coords");
  options.addUsage("The resulting reference.exclude.coords file can be used with the --exclude option.");

  if(argc <= 1)  {
    options.printUsage();
    return -1;  }
  
  if (!options.count("reference") || !file_exists(options["reference"].c_str())) {
    options.addUsage("");
    options.addUsage("You must supply the --reference option for input.");
    options.addUsage("If you feel you've received this message in error, please");
    options.addUsage("check to see that the reference file you supplied exists.");
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
  vector<string> cmd_args(argv, argv + argc);
  cmd_args.insert(cmd_args.begin() + 1, "random-mutations");
  gd.add_breseq_data("COMMAND", join(cmd_args, " "));
  gd.random_mutations(
      options["exclude"],
      options["type"],
      un(options["number"]),
      un(options["buffer"]),
      ref_seq_info[ref_seq_id],
      options.count("verbose")
      );

  gd.mutations_to_evidence(ref_seq_info, false);

  gd.reassign_unique_ids();
  
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

int do_header(int argc, char* argv[]) 
{
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

  options("name",             "Name for job.", "breseq");
  options("email",            "Email address to be added to launcher script. Informs you when a job starts and finishes.");
  options("time",             "Alloted time for job.", "14:00:00");

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

  string name          = options["name"];
  string log_dir       = cString(options["log_dir"]).trim_ends_of('/');
  string output_dir    = cString(options["output_dir"]).trim_ends_of('/');
  string runfile_path  = cString(options["runfile"]).trim_ends_of('/');
  string launcher_path = cString(options["launcher"]).trim_ends_of('/');

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

    if (refs.size() == 0) {
      cerr << ">> Skipping file because no #REFSEQ= header lines found." << endl ;
      continue;  
    }
    
    if (reads.size() == 0) {
      cerr << ">> Skipping file because no #READSEQ= header lines found." << endl ;
      continue;  
    }
        
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

  ofstream launcher(launcher_path.c_str());
  // #$ Parameters.
  fprintf(launcher, "#!/bin/csh\n");
  fprintf(launcher, "#$ -N %s\n", name.c_str());
  fprintf(launcher, "#$ -pe %uway %u\n", tasks, nodes);
  fprintf(launcher, "#$ -q normal\n");
  fprintf(launcher, "#$ -o %s.o$JOB_ID\n", name.c_str());
  fprintf(launcher, "#$ -l h_rt=%s\n", options["time"].c_str());
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
  fprintf(launcher, "setenv WORKDIR        .\n");
  fprintf(launcher, "\n");

  // Job submission.
  fprintf(launcher, "cd $WORKDIR/\n");
  fprintf(launcher, "$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE\n");

  launcher.close();
  SYSTEM("chmod +x " + launcher_path, true);

  return 0;
}

int do_mrna_stability(int argc, char *argv[])
{
  AnyOption options("gdtools mrna_stability [-o output.gd] [-r reference.gbk] input1.gd");
  options("output,o",  "Output fasta file for vienna RNA fold", "output.fa");
  options("reference,r", "Genbank, GFF3, or FASTA reference sequence file");
  options("flanking_sequence,f", "Number of bases on either side of synonymous SNP to consider for mRNA folding", "15");
  options("verbose,v", "Verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("Creates a fasta file with posistion of Synonymous SNP,"); 
  options.addUsage("base change, and mRNA sequence of location plus/minus given flanking sequence");
  options.addUsage("To be used in vienna RNA fold");
                                            
  if (options.getArgc() != 1) {//@ded any interest in allowing multiple inserts and merge with frequencies before calcuation?
    options.addUsage("");
    options.addUsage("Must provide exactly 1 input.gd.");
    options.printUsage();
    return -1;
  }

    UserOutput uout("ANNOTATE");
//@ded html not an option, and annotated GD file is not a desired output.
/*    
    bool html_output_mode = options.count("html");
    string output_file_name;
    if (options.count("output")) 
        output_file_name = options["output"];
    else
        output_file_name = /#(html_output_mode) ? "annotated.html" :@ded no html#/ "annotated.gd";
*/    
    vector<string> gd_path_names;
    for (int32_t i = 0; i < options.getArgc(); ++i) {
        gd_path_names.push_back(options.getArgv(i));
    }

//@ded uneeded for stability as only allowing 1 input file
// checks that more than 1 input and at least 1 reference is given
/*    if ( (gd_path_names.size() == 0) 
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
*/     
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
/*    
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
@ded not allowing for multiple inputs*/     
    vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
    uout("Reading input reference sequence files") << reference_file_names << endl;
    cReferenceSequences ref_seq_info;
    ref_seq_info.LoadFiles(reference_file_names);
    
    uout("Annotating mutations");
    ref_seq_info.annotate_mutations(gd, true, options.count("ignore-pseudogenes"));
//@ded end of used annotation commands
    uout("Annotation complete. Writing synonymous SNP mRNA sequences.");
    
/*    if (html_output_mode) {
        
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
        uout("Writing output Genome Diff file", options["output"]);//not needed either
//        gd.write(output_file_name);
    //}*/    
  diff_entry_list_t muts = gd.show_list();
  stringstream ss; //fasta build file 
  uint32_t flanking_sequence = from_string<uint32_t>(options["flanking_sequence"]);
    
  for (diff_entry_list_t::iterator it=muts.begin(); it!=muts.end(); it++)
  { 
    cDiffEntry& mut= **it;  
    if (mut._type == SNP && mut["snp_type"] == "synonymous") {
    // determine region to get sequence of
      const uint32_t mutation_stability_start = from_string<uint32_t>(mut["position"]) - flanking_sequence;
      const uint32_t mutation_stability_end = from_string<uint32_t>(mut["position"]) + flanking_sequence;
      
    // determine reference and mutation sequence
      string syn_ref_seq, mut_ref_seq;
      syn_ref_seq = mut["gene_strand"] == ">" 
          ? ref_seq_info.get_sequence_1(mut["seq_id"], mutation_stability_start, mutation_stability_end)
          : reverse_complement(ref_seq_info.get_sequence_1(mut["seq_id"], mutation_stability_start, mutation_stability_end));
      mut_ref_seq = syn_ref_seq;
      mut_ref_seq.replace (flanking_sequence,1,mut["gene_strand"] == ">" 
          ? mut["new_seq"]
          : reverse_complement(mut["new_seq"]));         

    //construct fasta file
      ss << ">" << mut["position"] << "_" << mut["codon_ref_seq"] << "-R" << endl;
      ss << syn_ref_seq << endl;
      ss << ">" << mut["position"] << "_" << mut["codon_new_seq"] << "-M" << endl;
      ss << mut_ref_seq << endl;
    }
  }  
    
//write fasta file
  ofstream myfile(options["output"].c_str());
  myfile << ss.str();
  return 0;                                            
}

int do_translate_proteome(int argc, char *argv[])
{
    AnyOption options("gdtools PROTEOME [-o output.fna] -r reference.gbk");
    options("output,o",  "Base name for FASTA file output", "output");
    options("reference,r", "Genbank, GFF3, or FASTA reference sequence file");
    options("verbose,v", "Verbose mode", TAKES_NO_ARGUMENT);
    options.processCommandArgs(argc, argv);
    
    options.addUsage("");
    options.addUsage("Creates two FASTA files, one normal and one with protein translating amber stops as X,"); 
    
    cString output_base_name = options["output"];
    
    if (!options.count("reference")) {
        options.printUsage();
        return -1;
    }
    
    cReferenceSequences ref_seqs;
    ref_seqs.LoadFiles(from_string<vector<string> >(options["reference"]));
    
    cString nfile_name(output_base_name + ".normal.fna");
    ofstream nfile(nfile_name.c_str());
    cString afile_name(output_base_name + ".amber.fna");
    ofstream afile(afile_name.c_str());
    cString ufile_name(output_base_name + ".unique.fna");
    ofstream ufile(ufile_name.c_str());
    cString bfile_name(output_base_name + ".amber.tab");
    ofstream bfile(bfile_name.c_str());
    
    string n_translation_table = cReferenceSequences::translation_tables[11];    
    string n_translation_table_1 = cReferenceSequences::initiation_codon_translation_tables[11];
    
    string a_translation_table = n_translation_table;
    string a_translation_table_1 = n_translation_table_1;
    a_translation_table[cReferenceSequences::codon_to_aa_index["TAG"]] = 'X';
    a_translation_table_1[cReferenceSequences::codon_to_aa_index["TAG"]] = 'X';

    
    cSequenceFeaturePtr feature_ptr(NULL);
    
    uint32_t total_protein_count = 0;
    uint32_t skipped_protein_count = 0;
    uint32_t amber_terminated_protein_count = 0;
    //list< pair<string,uint32_t> > protein_backup_codon_length_list;    
    
    for(vector<cAnnotatedSequence>::iterator its = ref_seqs.begin(); its !=ref_seqs.end(); its++) {
    
        cAnnotatedSequence& on_seq = *its;
        cSequenceFeatureList& feature_list = on_seq.m_features;

        for (cSequenceFeatureList::iterator it = feature_list.begin(); it != feature_list.end(); ++it) {
            cSequenceFeature& on_feature = **it;
            
            if (on_feature["type"] == "CDS") {
                
                total_protein_count++;
                
                if (on_feature.get_strand() != 0) {                
                    string n_protein = cReferenceSequences::translate_protein(on_seq, on_feature.m_location, n_translation_table, n_translation_table_1);
                    string a_protein = cReferenceSequences::translate_protein(on_seq, on_feature.m_location, a_translation_table, a_translation_table_1);
                    
                    int32_t length_diff = a_protein.size() - n_protein.size();
                    
                    nfile << ">" << on_feature["name"] << endl << n_protein << endl;
                    afile << ">" << on_feature["name"] << "_amber_" << length_diff << "_bp_to_new_stop" << endl << a_protein << endl;
                    
                    ufile << ">" << on_feature["name"] << endl << n_protein << endl;
                    if (length_diff) {
                        amber_terminated_protein_count++;
                        ufile << ">" << on_feature["name"] << "_amber_" << length_diff << "_bp_to_new_stop" << endl << a_protein << endl;
 
                        bfile << on_feature["name"] << "\t" << length_diff << "\t" << on_feature["locus_tag"] <<  "\t" << on_feature["product"] << endl;
                        //protein_backup_codon_length_list.push_back(make_pair(on_feature["name"],length_diff));
                    }
                }
                else
                {
                    skipped_protein_count++;
                }
                
            }
        }
    }
    
    cout << "Amber terminated proteins: " << amber_terminated_protein_count << endl;
    cout << "Skipped proteins: " << skipped_protein_count << endl;
    cout << "Total proteins: " << total_protein_count << endl;
    
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
  if (command == "VALIDATE") {
    return do_validate(argc_new, argv_new);
  } else if (command == "APPLY") {
    return do_apply(argc_new, argv_new);    
  } else if (command == "COMPARE") {
    return do_annotate(argc_new, argv_new);
  } else if (command == "CHECK") {
      return do_check(argc_new, argv_new);  
  } else if (command == "CHECK-PLOT") {
      return do_check_plot(argc_new, argv_new);
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
  } else if (command == "GD2VCF") {             
    return do_gd2vcf(argc_new, argv_new);
  } else if (command == "VCF2GD") {             
    return do_vcf2gd( argc_new, argv_new);
  } else if (command == "GD2GVF") {             
    return do_gd2gvf( argc_new, argv_new);
  } else if (command == "GD2CIRCOS"){
    return do_gd2circos(argc_new, argv_new);
  } else if(command == "MIRA2GD"){
    return do_mira2gd(argc_new, argv_new);
  } else if(command == "HEADER"){
    return do_header(argc_new, argv_new);
  } else if ((command == "RANDOM-MUTATIONS") || (command == "SIMULATE-MUTATIONS")) {
    return do_simulate_mutations(argc_new, argv_new);
  } else if (command == "MUTATIONS-TO-EVIDENCE") {
    return do_mutations_to_evidence(argc_new, argv_new);
  } else if (command == "DOWNLOAD") {
    return do_download(argc_new, argv_new);
  } else if (command == "RUNFILE") {
    return do_runfile(argc_new, argv_new);
  } else if (command == "MRNA-STABILITY") {
    return do_mrna_stability(argc_new, argv_new);
  } else if (command == "PROTEOME") {
      return do_translate_proteome(argc_new, argv_new);
  } else {
    cout << "Unrecognized command: " << command << endl;
    gdtools_usage();
  }
  return 0;

}

