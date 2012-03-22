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
#include "libbreseq/annotated_sequence.h"

using namespace breseq;
int do_intersection (int argc, char* argv[]);
int do_mutate       (int argc, char* argv[]);
int do_subtract     (int argc, char* argv[]);
int do_merge        (int argc, char* argv[]);
int do_compare      (int argc, char* argv[]);
int do_convert_gvf  (int argc, char* argv[]);
int do_convert_gd   (int argc, char* argv[]);
int do_not_evidence (int argc, char* argv[]);
int do_annotate     (int argc, char* argv[]);
int do_normalize_gd (int argc, char* argv[]);
int do_filter_gd    (int argc, char* argv[]);
int do_weights      (int argc, char* argv[]);




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
    //TODO for genome diff // Gives default usage in this case.
    return -1; 
	}

	//Pass the command to the proper handler.
	command = to_upper(command);

  // Genome Diff Commands:
  if (command == "APPLY") {
    return do_mutate(argc_new, argv_new);    
  } else if (command == "COMPARE") {
    return do_compare(argc_new, argv_new);
  } else if (command == "NOT-EVIDENCE") {
    return do_not_evidence(argc_new, argv_new);
  } else if (command == "ANNOTATE") {
    return do_annotate(argc_new, argv_new);
  } else if (command == "NORMALIZE") {
    return do_normalize_gd(argc_new, argv_new);
  } else if (command == "FILTER") {
    return do_filter_gd(argc_new, argv_new);

  } else if (command == "INTERSECT") {
    return do_intersection(argc_new, argv_new);
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
  }

  return 0;

}




int do_intersection(int argc, char *argv[])
{
  AnyOption options("INTERSECT -o <output.gd> <file1.gd file2.gd file3.gd ...>");
  options("output,o", "output intersection mutations to this file");
  options.processCommandArgs(argc, argv);

  typedef vector<string> string_vector_t;
  string_vector_t gd_file_names;
  //Get GD Files
  for (int32_t i = 0; i < options.getArgc() ; i++) {
    const string &file_name = options.getArgv(i);
    gd_file_names.push_back(file_name);
  }

  if (!options.count("output") ||
      gd_file_names.size() < 2) {
    printf("\n");
    printf("ERROR!\n");
    printf("Ouput:%s\n", options["output"].c_str());
    printf("GenomeDiff file count: %i\n", static_cast<int>(gd_file_names.size()));
    printf("\n");
    options.printUsage();
    return -1;
  }

  cGenomeDiff first_gd(*gd_file_names.begin());
  const diff_entry_list_t &first_mutations = first_gd.mutation_list();

  //Strip counted_ptr
  typedef set<cDiffEntry> diff_entry_set_t;
  diff_entry_set_t first_mutations_set;
  for (diff_entry_list_t::const_iterator it = first_mutations.begin();
       it != first_mutations.end(); it++) {
    first_mutations_set.insert(**it);
  }

  /*! Step: Compare to other genome diffs and reassign first_mutations_set when applicable.
    Stop when there is no intersection or we run out of file_names
  */
  string_vector_t::const_iterator it_file_name = gd_file_names.begin();
  advance(it_file_name,1);
  while(it_file_name != gd_file_names.end()) {

    cGenomeDiff second_gd(*it_file_name);
    const diff_entry_list_t &second_mutations = second_gd.mutation_list();

    //Strip counted_ptr
    diff_entry_set_t second_mutations_set;
    for (diff_entry_list_t::const_iterator it = second_mutations.begin();
         it != second_mutations.end(); it++) {
      second_mutations_set.insert(**it);
    }

    //Find intersection
    diff_entry_set_t intersecting_mutations_set;
    set_intersection(first_mutations_set.begin(), first_mutations_set.end(),
                     second_mutations_set.begin(), second_mutations_set.end(),
                     inserter(intersecting_mutations_set, intersecting_mutations_set.begin()));

    if (intersecting_mutations_set.empty()) {
      printf("No intersecting mutations were found across the current file: %s\n",
             it_file_name->c_str());
      return -1;
    }


    first_mutations_set = intersecting_mutations_set;

    it_file_name++;
  }

  printf("Found %i intersecting mutations across files: %s.\n",
         static_cast<int>(first_mutations_set.size()), join(gd_file_names, ", ").c_str());

  cGenomeDiff intersecting_gd;

  for (diff_entry_set_t::iterator it = first_mutations_set.begin();
       it != first_mutations_set.end(); it++) {
    //cDiffEntry de = *it;
    intersecting_gd.add(*it);
  }

  intersecting_gd.write(options["output"]);

	return 0;
}

int do_mutate(int argc, char *argv[])
{
  AnyOption options("APPLY -g <file.gd> -r <reference>");
  options("genomediff,g", "genome diff file");
  options("reference,r",".gbk/.gff3/.fasta/.bull reference sequence file");
  options("fasta,f","output FASTA file");
  options("gff3,3","output GFF3 file");
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Input a single GenomeDiff, and as many reference files");
  options.addUsage("as you like.  Using the GenomeDiff, this will apply all");
  options.addUsage("the mutations to the reference sequences, output is to");
  options.addUsage("a single file that includes all the references.");
  
  if(argc == 1)  {
    options.printUsage();
    return -1;  }
  
  if (!options.count("genomediff") ||
      !options.count("reference")) {
    options.addUsage("");
    options.addUsage("You must supply BOTH the --genomediff and --reference");
    options.addUsage("options for input.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("gff3") && !options.count("fasta")) {
    options.addUsage("");
    options.addUsage("You must supply at least one of the --fasta or --gff3 options for output.");
    options.printUsage();
    return -1;
  }  
  
  cGenomeDiff gd(options["genomediff"]);
  cReferenceSequences ref_seq_info;
  cReferenceSequences new_ref_seq_info;
  ref_seq_info.LoadFiles(from_string<vector<string> >(options["reference"]));
  new_ref_seq_info.LoadFiles(from_string<vector<string> >(options["reference"]));

  //Check to see if every item in the loaded .gd is
  // applicable to the reference file.
  ASSERT(gd.is_valid(ref_seq_info, options.count("verbose")), "Reference file and GenomeDiff file don't match.");
  
  gd.apply_to_sequences(ref_seq_info, new_ref_seq_info, options.count("verbose"));

  if (options.count("fasta"))
    new_ref_seq_info.WriteFASTA(options["fasta"], options.count("verbose"));
  if (options.count("gff3"))
    new_ref_seq_info.WriteGFF(options["gff3"], options.count("verbose"));

  return 0;
}

int do_subtract(int argc, char *argv[])
{
  AnyOption options("SUBTRACT -1 <file.gd> -2 <file.gd> -o <output.gd>");
  options("input1,1","input GD file 1");
  options("input2,2","input GD file 2");
  options("output,o","output GD file");
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Input two GenomeDiff files and compare mutations.");
  options.addUsage("Mutations that appear in BOTH will be subtracted from");
  options.addUsage("--input1 and results will output to a new file specified.");
  
  if(argc == 1)  {
    options.printUsage();
    return -1;  }
  
  if (!options.count("input1") ||
      !options.count("input2")) {
    options.addUsage("");
    options.addUsage("You must supply BOTH the --input1 and --input2 options");
    options.addUsage("for input.");
    options.printUsage();
    return -1;
  }
  
  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("You must supply the --output option for output.");
    options.printUsage();
    return -1;
  }
  
  cGenomeDiff gd1(options["input1"]);
  cGenomeDiff gd2(options["input2"]);
  
  gd1.subtract(gd2, options.count("verbose"));
  
  gd1.write(options["output"]);
  
  return 0;
}

int do_merge(int argc, char *argv[])
{
  AnyOption options("MERGE -o <output.gd> <file1.gd file2.gd file3.gd ...>");
  options("output,o",     "Output GD file");
  options("unique,u",     "Unique Entries Only (Flag)", TAKES_NO_ARGUMENT);  
  options("id,i",         "Reorder IDs (Flag)", TAKES_NO_ARGUMENT);
  options("verbose,v",    "Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
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
    options.addUsage("You must supply the --output option for output.");
    options.printUsage();
    return -1;
  }

  const bool verbose = options.count("verbose");
  
  cGenomeDiff gd1(options.getArgv(0));
  
  //Load all the GD files that were input.
  for(uint32_t i = 1; i < options.getArgc(); i++)
  {
    cGenomeDiff gd2(options.getArgv(i));
    gd1.merge(gd2, options.count("unique"), options.count("id"), options.count("verbose"));
  }
  
  gd1.write(options["output"]);
  
  return 0;
}

int do_weights(int argc, char* argv[])
{
  AnyOption options("WEIGHTS -o <output.gd> <file1.gd file2.gd file3.gd ...>");
  options("output,o",    "Output GD file");
  options("verbose,v",   "Verbose mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);


  if (!options.getArgc()) {
    options.addUsage("");
    options.addUsage("No input provided.");
    options.printUsage();
    return -1;
  }

  if (!options.count("output")) {
    options.addUsage("");
    options.addUsage("You must supply the --output option for output.");
    options.printUsage();
    return -1;
  }

  const bool verbose = options.count("verbose");

  cGenomeDiff gd1(options.getArgv(0));

  //Load all the GD files that were input.
  for(uint32_t i = 1; i < options.getArgc(); i++)
  {
    cGenomeDiff gd2(options.getArgv(i));
    gd1.fast_merge(gd2);
  }

  //Store entries as keys and pointers to matching entries.
  diff_entry_list_t muts = gd1.mutation_list();
  map<cDiffEntry, vector<diff_entry_ptr_t> > weights_table;

  for (diff_entry_list_t::iterator it = muts.begin(); 
       it != muts.end(); ++it) {
    cDiffEntry mut = **it;
    weights_table[mut].push_back(*it);
  }

  //Given the number of pointers, determine frequency and add 'weight' field to entries.
  for (map<cDiffEntry, vector<diff_entry_ptr_t> >::iterator it = weights_table.begin();
       it != weights_table.end(); ++it) {
    float weight = 1.0f / static_cast<float>(it->second.size());
    
    for (uint32_t i  = 0; i < it->second.size(); ++i) {
      (*it->second[i])["weight"] = to_string(weight);
    }

    if (verbose) {
      cout << "\tEntry: " << it->first << endl;
      cout << "\tMatches with:" << endl;
      for (uint32_t i  = 0; i < it->second.size(); ++i) {
        cout << "\t\t[entry]: " << (*it->second[i]) << endl;
      }
      cout << endl;
      cout << endl;
    }

  }
  
  gd1.write(options["output"]);

  return 0;
}


int do_compare(int argc, char *argv[])
{
  AnyOption options("COMPARE-GD -c <control.gd> -t <test.gd> -o <output.gd>");
  options("control,c", "control genome diff file, mutations within are assumed to be accurate");
  options("test,t",    "test genome diff file, mutations not been checked for accuracy");
  options("output,o",  "output compared genome diff file name");
  options.processCommandArgs(argc, argv);


  if (!options.count("control") ||
      !options.count("test")    ||
      !options.count("output")) {
    options.printUsage();
    return -1;
  }

  cGenomeDiff control_gd(options["control"]);
  cGenomeDiff test_gd(options["test"]);

  cGenomeDiff compare_gd = cGenomeDiff::compare_genome_diff_files(control_gd, test_gd);

  compare_gd.write(options["output"]);


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
  AnyOption options("ANNOTATE -r reference.gbk  -i input.gd -o annotated.gd");
  
  options
  ("help,h", "produce advanced help message", TAKES_NO_ARGUMENT)
  // convert to basing everything off the main output path, so we don't have to set so many options
  ("input,i", "path to input genome diff (REQUIRED)")
  ("output,o", "path to output genome diff with added mutation data (REQUIRED)")
  ("reference,r", "reference sequence in GenBank flatfile format (REQUIRED)")
  ("ignore-pseudogenes", "treats pseudogenes as normal genes for calling AA changes", TAKES_NO_ARGUMENT)
  ;
  options.processCommandArgs(argc, argv);
  
  if ( !options.count("input")
    || !options.count("output")
    || !options.count("reference")
      ) {
    options.printUsage();
    return -1;
  }
  
  cGenomeDiff gd(options["input"]);
  
  vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(reference_file_names);
  ref_seq_info.annotate_mutations(gd, false, options.count("ignore-pseudogenes"));
  gd.write(options["output"]);
  
  return 0;
}

int do_normalize_gd(int argc, char* argv[])
{
  AnyOption options("NORMALIZE-GD -g <input.gd> -r <reference> -o <output.gd>");
  options
  ("genome_diff,g", "Input genome diff file.")
  ("reference,r"  , "Input reference file.")
  ("output,o"     , "Output normalized genome diff file.")
  ("verbose,v"    , "Verbose Mode (Flag)", TAKES_NO_ARGUMENT);

  options.processCommandArgs(argc, argv);
  options.addUsage("");
  options.addUsage("Takes a genome diff file to be normalized.");
  options.addUsage("Then normalizes the diff entries to the reference.");
  options.addUsage("Outputs a normalized genome diff file.");
  
  if(argc == 1)  {
    options.printUsage();
    return -1;  }

  if (!options.count("genome_diff")) {
    options.addUsage("");
    options.addUsage("You must supply the --genome_diff option for input.");
    options.printUsage();
    return -1;
  }

  if (!options.count("reference")) {
    options.addUsage("");
    options.addUsage("You must supply the --reference option for input.");
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

  vector<string> rfns = from_string<vector<string> >(options["reference"]);
  cReferenceSequences rs;
  rs.LoadFiles(rfns);

  cGenomeDiff gd(options["genome_diff"]);

  printf("\n++Normalizing genome diff file: %s to reference file: %s \n\n",
          options["genome_diff"].c_str(), join(rfns,",").c_str());

  gd.normalize_to_sequence(rs);
  const diff_entry_list_t &muts = gd.mutation_list();

  cGenomeDiff new_gd;
  for (diff_entry_list_t::const_iterator it = muts.begin();
       it != muts.end(); it++) {
    cDiffEntry de = **it;

    if (de.entry_exists("norm") && de["norm"] == "is_not_valid"){
      printf("\tINVALID_MUTATION:%s\n", de.to_string().c_str());
      de["comment_out"] = "True";
    } 
    
    new_gd.add(de);
  }
  new_gd.write(options["output"]);

  printf("\n++Normilization completed.\n\n");

  return 0;
}

int do_filter_gd(int argc, char* argv[]) {
  AnyOption options("FILTER -g <input.gd> -o <output.gd> -m <mutation type> <filters>");
  options("genome_diff,g", "Genome diff file to be filtered.");
  options("output_gd,o", "Filter genome diff file.");
  options("mut_type,m", "Only consider this mutation type for filtering.");
  options.processCommandArgs(argc, argv);

  if (!options.count("genome_diff")) {
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
  if (!options.count("output_gd")) {
    options["output_gd"] = options["genome_diff"];
  }

  list<string> filters;
  for(int32_t i = 0; i < options.getArgc(); ++i) {
    filters.push_back(options.getArgv(i));
  }
  assert(filters.size());

  cGenomeDiff gd(options["genome_diff"]);

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
      mut["zfiltered"] = join(reasons, ", ");
      mut["comment_out"] = "true";
    }
    output_gd.add(mut);
  }
  output_gd.write(options["output_gd"]);
  return 0;
}


