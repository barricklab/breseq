/*****************************************************************************

 AUTHORS

   Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com> and other contributors

 LICENSE AND COPYRIGHT

   Copyright (c) 2008-2010 Michigan State University
   Copyright (c) 2011-2025 The University of Texas at Austin
   Copyright (c) 2025-     Michigan State University

   breseq is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2, or (at your option) any later version.

   SPDX-License-Identifier: GPL-2.0-or-later

*****************************************************************************/

#include "libbreseq/common.h"
#include "libbreseq/settings.h"
#include "libbreseq/anyoption.h"
#include "libbreseq/mutation_predictor.h"
#include "libbreseq/flagged_regions.h"
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
  uout << "VALIDATE               check formatting of input files" << endl;
  uout << "APPLY                  apply mutations to a sequence" << endl;
  uout << "ANNOTATE (or COMPARE)  annotate the effects of mutations and compare multiple samples" << endl;
	uout << "MUTATIONS              (re)predict mutations from evidence" << endl;
  uout << "CHECK                  compare control versus test mutations" << endl;
  uout << "NORMALIZE              normalize mutation positions and annotations" << endl;
  //uout << "header                 create or add header entries" << endl;

  uout("Set and Filtering Operations:");
  uout << "SUBTRACT               remove mutations in one file from another" << endl;
  uout << "INTERSECT              keep shared mutations in two files" << endl;
  uout << "UNION/MERGE            combine mutations, removing duplicates" << endl;
	uout << "FILTER/REMOVE          remove mutations matching specified conditions" << endl;
	uout << "MASK                   remove mutation predictions in masked regions" << endl;
  uout << "NOT-EVIDENCE           remove evidence not used by any mutations" << endl;
	
  uout("Format Conversions:");
  uout << "GD2VCF                 GD to Variant Call Format (VCF)" << endl;
	uout << "VCF2GD                 Variant Call Format(VCF) to GD" << endl;
  uout << "GD2GVF                 GD to Genome Variation Format (GVF)" << endl;
	uout << "MUMMER2MASK            Create a mask GD file from MUMmer output" << endl;

	uout("Analysis:");
	uout << "COUNT                  count statistics for different types of mutations" << endl;
	uout << "PHYLOGENY              create maximum parsimony tree from mutations (requires PHYLIP)" << endl;
	
  uout("TACC Utilities:");
  uout << "DOWNLOAD               download reference and read files from GD header info" << endl;

  return 0;
}

int do_intersection(int argc, char *argv[])
{
  AnyOption options("gdtools INTERSECT [-o output.gd] input1.gd input2.gd ...");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("output,o",  "Output Genome Diff file name", "output.gd");
  options("verbose,v", "Verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("Creates a new Genome Diff file with mutations that are present in ALL input Genome Diff files.");

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
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

int do_merge(int argc, char *argv[])
{
  AnyOption options("gdtools MERGE/UNION [-o output.gd] input1.gd input2.gd ...");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("output,o",  "output GD file name", "output.gd");
	options("preserve-evidence,e", "By default evidence items with two-letter codes are removed (RA, JC, MC, ...). Supply this option to retain them. ", TAKES_NO_ARGUMENT);
	options("phylogeny-aware,p", "Do not consider mutations that differ in their 'phylogeny_id' fields equivalent.", TAKES_NO_ARGUMENT);
	options("population-aware,q", "Do not consider mutations in GD files with different POPULATION metadata values equivalent.", TAKES_NO_ARGUMENT);
  options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("Creates a GD file with entries that exist in any input GD file.");
	options.addUsage("Duplicate entries are merged into a single entry.");
	options.addUsage("");
	options.addUsage("Header information will be inherited from the first input file,");
	options.addUsage("so this function can also be used to transfer metadata to a new file.");

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
	/*
  if (options.getArgc() < 1) {
    options.addUsage("");
    options.addUsage("Must provide at least two input Genome Diff files.");
    options.printUsage();
    return -1;
  }
	 
	 */
  UserOutput uout("MERGE");
	cout << endl << "    Preserving: " << (!options.count("evidence") ? "Mutations (3-letter codes)" : "Evidence (2-letter codes)") << endl;

  uout("Reading input GD files") << endl;
	
	cGenomeDiff gd1(options.getArgv(0));
	uout << "    " << options.getArgv(0) << endl;
  for(int32_t i = 1; i < options.getArgc(); ++i) {
    uout << "    " << options.getArgv(i) << endl;
    cGenomeDiff gd2(options.getArgv(i));
    gd1.set_union(gd2, options.count("preserve-evidence"), options.count("phylogeny-aware"), options.count("population-aware"), options.count("verbose"));
  }

  uout("Assigning unique IDs");
	gd1.reassign_unique_ids();

  uout("Writing output GD file", options["output"]); 
  gd1.write(options["output"]);

  return 0;
}

int do_apply(int argc, char *argv[])
{
  AnyOption options("gdtools APPLY [ -o output.gff3 -f GFF3 -s seq_id ] -r reference.gbk input.gd");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("output,o",    "Output file name (DEFAULT=output.*)");
  options("format,f",    "Output file format (Options: FASTA, GENBANK, GFF3)", "FASTA");
  options("reference,r", "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
  options("seq-id,s",    "Sequence ID to keep in output. If this argument is provided, other sequences are deleted after the APPLY. May be provided multiple times.");
	options("polymorphism-mode,p",  "Apply all mutations in GD file regardless of their frequency. By default (without this option) only mutations with 100% frequency are applied.", TAKES_NO_ARGUMENT);
	options("applied-gd,a",  "Output file name for GD with mutations updated to coordinates in the output sequences.");
  options("verbose,v",   "Verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
	options.addUsage("Input a single GenomeDiff, and as many reference files as you like.  Using the GenomeDiff, this will apply all the mutations to the reference sequences, output is to a single file that includes all the references in the requested format.");
	options.addUsage("");
	options.addUsage("The input GenomeDiff file is expected to only have consensus mutations. Polymorphic mutations are ignored unless the --polymorphism-mode flag is supplied.");
	options.addUsage("");
	options.addUsage("The --apply-gd option causes a GenomeDiff file to be output that is the input GenomeDiff file with positions of mutations shifted to where they occur in the output sequence. It also has applied_seq_id, applied_start, and applied_end fields defining the changed bases.");
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
  if (!options.count("reference")) {
    options.addUsage("");
    options.addUsage("No reference files provided.");
    options.printUsage();
  }
	

	string format = to_upper(options["format"]);
	if (format == "GFF") format = "GFF3"; // allow synonym
  if ((format != "FASTA") && (format != "GENBANK") && (format != "GFF3")) {
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
		if (format=="FASTA") output = "output.fasta";
		if (format=="GFF3") output = "output.gff";
		if (format=="GENBANK") output = "output.gbk";
  }  

  if (!options.getArgc()) {
    options.addUsage("");
    options.addUsage("No input files provided.");
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
  gd.valid_with_reference_sequences(ref_seq_info);
	
	// Set all frequencies to one
	if (options.count("polymorphism-mode")) {
		diff_entry_list_t mutation_list = gd.mutation_list();
		for (diff_entry_list_t::iterator itr_mut = mutation_list.begin(); itr_mut != mutation_list.end(); itr_mut++)
		{
			cDiffEntry& mut(**itr_mut);
			
			if ((mut.count(FREQUENCY)!=0) && (from_string<double>(mut[FREQUENCY]) != 1.0))
				mut[FREQUENCY] = "1";
		}
	}
	
	gd.apply_to_sequences(ref_seq_info, new_ref_seq_info, options.count("verbose"));

	// If seq-id is present keep only certain sequence ids
	// Keep the order of the --seq-id|-s to allow user to specify this
	if (options.count("seq-id")) {
		const vector<string>& seq_ids = from_string<vector<string> >(options["seq-id"]);
		cReferenceSequences kept_ref_seq_info;
		set<string> used_ids_se;
		
		for (vector<string>::const_iterator i=seq_ids.begin(); i != seq_ids.end(); i++) {
			if (used_ids_se.find(*i) == used_ids_se.end()) {
				
				if ( !new_ref_seq_info.seq_id_exists(*i) ) {
					ERROR("Specified--seq-id|-s does not exist in loaded files: " + *i);
				} else {
					kept_ref_seq_info.push_back(new_ref_seq_info[*i]);
				}
			} else {
				WARN("Duplicate --seq-id|-s encountered and ignored: " + *i);
			}
		}
		new_ref_seq_info = kept_ref_seq_info;
	}
	
  uout("Writing output file in " + format + " format");

  if (format == "FASTA") {
    new_ref_seq_info.WriteFASTA(output);
  }
  else if (format == "GENBANK") {
    new_ref_seq_info.WriteGenBank(output);
  }
  else if (format == "GFF3") {
    new_ref_seq_info.WriteGFF(output);
  }

	if (options.count("applied-gd")) {
		uout("Writing applied GenomeDiff file:" + options["applied-gd"]);
		gd.write(options["applied-gd"]);
	}
	
  return 0;
}

int do_subtract(int argc, char *argv[])
{
  AnyOption options("gdtools SUBTRACT [-o output.gd] input.gd subtract1.gd [subtract2.gd ...]");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("output,o",  "output GD file", "output.gd");
	options("phylogeny-aware,p", "Do not consider mutations that differ in their 'phylogeny_id' fields equivalent.", TAKES_NO_ARGUMENT);
	options("population-aware,q", "Do not consider mutations in GD files with different POPULATION metadata values equivalent.", TAKES_NO_ARGUMENT);
	options("frequency-aware,f", "Use the frequencies of mutations when performing the subtraction. Normally an input mutation is removed if it appears at any frequency in a subtracted file. In this mode its frequency is reduced by the frequency in each subtracted file. If the resulting frequency is zero or below, then the mutation is removed.", TAKES_NO_ARGUMENT);

  options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
	options.addUsage("Creates a new Genome Diff file that contains all mutation entries that are still");
	options.addUsage("present in the input file after removing mutations that are in any of the");
	options.addUsage("subtracted Genome Diff files. All evidence and validation entries are");
	options.addUsage("retained in the output file.");

  if (options.getArgc() < 2) {
    options.addUsage("");
    options.addUsage("Input Genome Diff and at least one Genome Diff to subtract not provided.");
    options.printUsage();
    return -1;
  }
	
	if (options.count("help")) {
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
    gd1.set_subtract(gd2, options.count("phylogeny-aware"), options.count("population-aware"), options.count("frequency-aware"), verbose);
  }
  
  uout("Writing output GD file", options["output"]);
  gd1.write(options["output"]);
  
  return 0;
}

int do_weights(int argc, char* argv[])
{
  AnyOption options("gdtools WEIGHTS [-o output.gd input1.gd input2.gd ...");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("output,o",    "output GD file", "output.gd");
  options("verbose,v",   "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("\tCreates a GD file of mutations and associated evidence that exist in each input GD file.");
  options.addUsage("Duplicate mutations are merged into a single mutation with a 'weight' field added to the");
  options.addUsage("entry. The 'weight' field is the inverse of the frequency of that mutation occurring across");
  options.addUsage("all input GD files. Unique mutations will therefor have a 'weight' of 1");

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
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
    gd1.merge_preserving_duplicates(gd2);
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

// Check format of Genome Diff and report all format errors
int do_validate(int argc, char *argv[])
{
	AnyOption options("gdtools VALIDATE -r reference.gbk input1.gd [input2.gd]");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("reference,r",  "File containing reference sequences in GenBank, GFF3, or FASTA format. If provided, will validate seq_ids and positions in the GD file using these.  Option may be provided multiple times for multiple files. (OPTIONAL)");
	options("verbose,v",  "Verbose mode. Outputs additional information about progress. (OPTIONAL)");
	options.addUsage("");
	options.addUsage("Validates whether the format of the input Genome Diff files is correct.");

	options.processCommandArgs(argc, argv);

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
	bool verbose = options.count("verbose");
	
	if (options.getArgc() == 0) {
			options.printUsage();
			return -1;
	}
	
	cerr << "Running gdtools VALIDATE on " << options.getArgc() << " file(s)..." << endl;
	
	// Further read in reference sequences and check IDs of mutations if user provided
	// Give a warning if they weren't provided to remind the user this is possible.
	cReferenceSequences ref;
	if (options.count("reference")) {
			ref.LoadFiles(from_string<vector<string> >(options["reference"]));
	} else {
			cerr << "  WARNING: No reference sequence(s) provided (-r option). Genomic coordinates will not be checked." << endl;
	}
	
	
	uint32_t num_files_processed = 0;
	uint32_t num_files_with_errors = 0;
	
	// Simply read files.
	for (int32_t i=0; i<options.getArgc(); i++) {
			string gd_file_name = options.getArgv(i);
			num_files_processed++;
		
			cerr << endl << "File: " << gd_file_name << endl;

			cGenomeDiff gd;
			cFileParseErrors pe = gd.read(gd_file_name, true);

			// we can't test the next part if there are fatal errors
			if (pe.fatal()) {
				pe.print_errors(false);
				num_files_with_errors++;
				continue;
			}
		
			cFileParseErrors pe2;
			if (options.count("reference")) {
					pe2 = gd.valid_with_reference_sequences(ref, true);
			}
			
			if (pe._errors.size() + pe2._errors.size() == 0) {
				cerr << "  FORMAT OK" << endl;
			} else {
				cerr << endl << "File: " << gd_file_name << endl;
				pe.print_errors(false);
				pe2.print_errors(false);
				num_files_with_errors++;
			}
	}
	
	cerr << endl;
	cerr <<     "::: Summary :::" << endl;
	cerr <<     "  Total files processed:         " << setw(6) << num_files_processed << endl;
	if (num_files_with_errors == 0) {
			cerr << "* No formatting errors found!" << endl;
	} else {
			cerr << "  Files with formatting errors:  " << setw(6) << num_files_with_errors << endl;
	}
	
	return 0;
}

int do_check(int argc, char *argv[])
{
  AnyOption options("gdtools CHECK [-o output.gd] control.gd test.gd");
	options("help,h", "display detailed help message", TAKES_NO_ARGUMENT);
  options("output,o",         "output GD file", "comp.gd");
  options("reference,r",      "file containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
  options("evidence",         "compare evidence (only JC at present)", TAKES_NO_ARGUMENT);
  options("jc-buffer",        "when comparing JC evidence, length of sequence segment to compare for JC evidence", 50);
  options("jc-shorten",       "when comparing JC evidence, length to shorten control segments by when comparing JC evidence for overlap", 5);
  options("jc-only-accepted", "when comparing JC evidence, do not score/count rejected items", TAKES_NO_ARGUMENT);
  options("verbose,v",        "verbose mode", TAKES_NO_ARGUMENT);
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
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
  if (options.getArgc() != 2) {
    options.addUsage("");
    options.addUsage("Exactly two input Genome Diff files must be provided.");
    options.printUsage();
    return -1;
  }

  if (options.count("evidence") && !options.count("reference")) {
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
  if (options.count("evidence")) {

    cReferenceSequences ref;
    ref.LoadFiles(from_string<vector<string> >(options["reference"]));

    uout("Comparing evidence");
    comp = cGenomeDiff::check_evidence(ref, un(options["jc-buffer"]), un(options["jc-shorten"]), ctrl, test, options.count("jc-only-accepted"), options.count("verbose"));

  } else {
    uout("Comparing mutations");
    comp = cGenomeDiff::check(ctrl, test, options.count("verbose"));
  }

  uout("Assigning unique IDs to Genome Diff entries");
  comp.reassign_unique_ids();

  uout("Writing output GD file", options["output"]);
  comp.write(options["output"]);

  return 0;
}

int do_convert(int argc, char* argv[], string forced_format = "")
{
	string usage_command = forced_format.empty() ? "CONVERT -f <format>" : forced_format;
	AnyOption options("gdtools " + usage_command + " [-o <output_file> -r <refseq>] <input_file>");

	options("help,h", "display detailed help message", TAKES_NO_ARGUMENT);
	if (forced_format.empty())
		options("format,f","file format to output: GD, VCF, GVF, or JSON (REQUIRED)");
	options("output,o","name of output file (DEFAULT = <input_file>.*)", "");
		options("reference,r",  "file containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED for VCF input or with annotate option)");
	options("annotate,a","annotate mutations and evidence in input file", TAKES_NO_ARGUMENT);
	options("snv-only","only include SNP/SNV entries in GVF output", TAKES_NO_ARGUMENT);

	options.processCommandArgs( argc,argv);

	options.addUsage("");
	options.addUsage("Convert a single file from/to GenomeDiff format.");
	options.addUsage("Allowed input formats: GD, VCF");
	options.addUsage("");

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}

	UserOutput uout(usage_command);

	string output_format;
	if (!forced_format.empty()) {
		output_format = to_upper(forced_format);
	} else if (options.count("format")) {
		output_format = to_upper(options["format"]);
	} else {
		options.addUsage("");
		options.addUsage("You must supply the output file format option (--format,-f).");
		options.printUsage();
		return -1;
	}

	if( options.getArgc() != 1 ){
		options.addUsage("");
		options.addUsage("You must provide exactly one input Genome Diff file.");
		options.printUsage();
		return -1;
	}
	string input_file_name = options.getArgv(0);
	
	// Set default output file name
	string output_file_name = options["output"];
	if (output_file_name.size() == 0) {
		string default_ending = to_lower(output_format);
		output_file_name = cString(input_file_name).get_base_name_no_extension();
		output_file_name += "." + default_ending;
	}
	
	cReferenceSequences ref_seq_info;
	if ( (output_format == "VCF") || (output_format == "GVF") || options.count("annotate") ) {
		if (!options.count("reference")) {
			options.addUsage("");
			if (output_format == "VCF") {
				options.addUsage("You must provide a reference sequence file (-r) for VCF output.");
			} else if (output_format == "GVF") {
				options.addUsage("You must provide a reference sequence file (-r) for GVF output.");
			} else if (options.count("annotate")) {
				options.addUsage("You must provide a reference sequence file (-r) with the annotate (-a) option.");
			}
			options.printUsage();
			return -1;
		}
		ref_seq_info.LoadFiles(from_string<vector<string> >(options["reference"]));
	}

	// Load the input file
	
	// Sniff the first line of the file to determine format
	ifstream sniff(input_file_name);
	ASSERT(sniff.is_open(), "Could not open input file: " + input_file_name);
	string first_line;
	getline(sniff, first_line);
	first_line = to_upper(first_line);
	sniff.close();
	
	// Load as GenomeDiff
	uout("Reading file", input_file_name);
	cGenomeDiff gd;

	if (first_line.find("##FILEFORMAT=VCF") != string::npos) {
		uout("Input format", "VCF");
		gd.read_vcf(input_file_name);
	} else if (first_line.find("#=GENOME_DIFF") != string::npos) {
		uout("Input format", "GD");
		gd.read(input_file_name);
	} else {
		ERROR("Could not determine file format of input file. Expected first line to contain:\n   #=GENOME_DIFF for GenomeDiff\n   ##FILEFORMAT=VCF for VCF");
	}
	
	// Annotate if requested - not just mutations
	if (options.count("annotate")) {
		uout("Annotating mutations");
		ref_seq_info.annotate_mutations(gd, false, false);
	}
					
	// Write in format requested
	uout("Writing file", output_file_name);
	uout("Output format", output_format);
	if ( (output_format == "GENOMEDIFF") || (output_format == "GD") ) {
		gd.write(output_file_name);
	} else if (output_format == "VCF") {
		gd.write_vcf(output_file_name, ref_seq_info);
	} else if (output_format == "GVF") {
		gd.write_gvf(output_file_name, ref_seq_info, options.count("snv-only"));
	} else if (output_format == "JSON") {
		gd.write_json(output_file_name);
	}
	
	return 0;
}

int do_not_evidence(int argc, char *argv[])
{
  AnyOption options("NOT_EVIDENCE [-o output.gd] input.gd");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("output,o","output GD file", "output.gd");
  options("id,i","Reorder IDs (Flag)", TAKES_NO_ARGUMENT);
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
	
  options.addUsage("");
  options.addUsage("Takes a GenomeDiff file and removes all of the entries that");
  options.addUsage("are NOT used as evidence by a mutation.  Outputs to a new");
  options.addUsage("GenomeDiff file if specified.  If no output is specified,");
  options.addUsage("verbose will still inform what evidence isn't being used.");
	
	options.processCommandArgs(argc, argv);

	if (options.count("help")) {
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
  
  cGenomeDiff gd1(input);
  
  gd1.filter_not_used_as_evidence(options.count("verbose"));
  
  if(options.count("id"))
  {
    cGenomeDiff gd2;
    gd2.merge(gd1, true);
    gd2.write(options["output"]);
  }
  else  {
    gd1.write(options["output"]);  
  }
  
  return 0;
}


// How does polymorphism_search_found work?
//   If it is NULL then nothing happens and compare_mode is honored
//   If it is not null then it returns TRUE/FALSE whether polymorphisms were found and swtiches the output to compare mode
//      if there is only one file provided but it has polymorphisms. This is only needed for HTML output!
void load_merge_multiple_gd_files(cGenomeDiff& gd, vector<cGenomeDiff>& gd_list, vector<string>& gd_path_names, vector<string>& gd_titles, cReferenceSequences& ref_seq_info, bool strip_to_mutations_and_unknown, bool* polymorphism_search_found, bool compare_mode, AnyOption& options, UserOutput& uout)
{
	
	if (polymorphism_search_found != NULL) {
		*polymorphism_search_found = false; // handles putting in the polymorphism column if only one file provided
	}
	
	for (uint32_t i = 0; i < gd_path_names.size(); i++) {
		uout("Reading input GD file",gd_path_names[i]);
		cGenomeDiff single_gd(gd_path_names[i]);
		
		// This allows a switch to compare mode so the correct columns are present if polymorphisms are found so they can be displayed properly in HTML
		if ((polymorphism_search_found != NULL) && !compare_mode && (i==0)) {
			diff_entry_list_t muts = single_gd.mutation_list();
			for(diff_entry_list_t::iterator it=muts.begin(); it != muts.end(); it++) {
				cDiffEntry de = **it;
				if (de.count(FREQUENCY) && from_string<double>(de[FREQUENCY]) != 1.0) {
					*polymorphism_search_found = true;
					break;
				}
			}
		}
		
		if (strip_to_mutations_and_unknown) {
			single_gd.remove_all_but_mutations_and_unknown();
		}
		
		// Safety for downstream operations - must strip_to_mutations_and_unknown if in compare mode
		ASSERT(!compare_mode || strip_to_mutations_and_unknown, "Unable to merge multiple GD files and preserve evidence entries.");
		
		// Clean to the desired region here to avoid
		if (options.count("region")) {
			single_gd.filter_to_within_region(ref_seq_info, options["region"]);
		}
		
		// Decide whether to merge in a new column
		if ( (!options.count("collapse")) || (single_gd.mutation_list().size() > 0)) {
			gd.merge(single_gd, false, options.count("phylogeny-aware"), options.count("population-aware"));
			//cout << gd.get_list().size();
			gd_list.push_back(single_gd); // it's important to add a copy that has UN items intact
		}
		
		// Copy over the metadata if there is only one file
		if (gd_path_names.size() == 1) {
			gd.metadata = single_gd.metadata;
		}
	}
	
	// Sort the full merged list
	gd.sort();
	
	// Also sort the order of the genome diff files to take into account
	// population, time, clone, etc. information
	
	cGenomeDiff::sort_gd_list_by_treatment_population_time(gd_list);

	// Then add frequency columns for all genome diffs
	if (compare_mode || ( (polymorphism_search_found != NULL) && (*polymorphism_search_found) ) ) {
		uout("Tabulating frequencies of mutations across all files");
		cGenomeDiff::tabulate_mutation_frequencies_from_multiple_gds(gd, gd_list, gd_titles, options.count("phylogeny-aware"), options.count("population-aware"));
	}
	
}


int do_annotate(int argc, char* argv[])
{
  AnyOption options("gdtools ANNOTATE/COMPARE [-o annotated.html] -r reference.gbk input.1.gd [input.2.gd ... ]");
	
  options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("output,o", "Path to output file with added mutation data. (DEFAULT: output.*)");
	options("reference,r", "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
	options("format,f", "Type of output file to generate. See options below", "HTML");
	options("add-html-fields,a", "Add formatted fields that are used for generating HTML output. Only applicable to GD and JSON output formats", TAKES_NO_ARGUMENT);
	options("add-text-fields,b", "Add formatted fields in UTF-8 encoded text that are similar to those used in HTML output. Only applicable to GD, TSV, CSV, and JSON output formats", TAKES_NO_ARGUMENT);
	options("ignore-pseudogenes", "Treat pseudogenes as valid ORFs for calling amino acid substitutions", TAKES_NO_ARGUMENT);
	options("repeat-header", "In HTML mode, repeat the header line every this many rows (0=OFF)", "0");
	options("phylogeny-aware,p", "Do not consider mutations that differ in their 'phylogeny_id' fields equivalent.", TAKES_NO_ARGUMENT);
	options("population-aware,q", "Do not consider mutations in GD files with different POPULATION metadata values equivalent.", TAKES_NO_ARGUMENT);
	options("region,g", "Only show mutations that overlap this reference sequence region (e.g., REL606:64722-65312)");
	options("preserve-evidence,e", "By default evidence items with two-letter codes are removed (RA, JC, MC, ...). Supply this option to retain them. Only affects output in GD and JSON formats. This option can only be used with a single input GD file (i.e., not in COMPARE mode). ", TAKES_NO_ARGUMENT);
	options("collapse,c", "Do not show samples (columns) unless they have at least one mutation", TAKES_NO_ARGUMENT);
	options("inactivating-overlap-fraction", "Mutations within this fraction of the length of a gene from its beginning are assigned to the 'genes_inactivating' versus 'the genes_overlapping' list if they fulfill other criteria.", cReferenceSequences::k_inactivating_overlap_fraction);
	options("inactivating-size-cutoff", "INS, DEL, and SUB mutations in genes that are longer than this length cutoff are always assigned to the 'genes_inactivating' list, even if they are in-frame or in noncoding genes.", cReferenceSequences::k_inactivating_size_cutoff);
	options("promoter-distance", "Mutations upstream and within this distance of the beginning of a gene have it added to their 'genes_promoter' list.", cReferenceSequences::k_promoter_distance);
	options("genbank-field-for-seq-id", "Which GenBank header field will be used to assign sequence IDs. Valid choices are LOCUS, ACCESSION, and VERSION. The default is to check those fields, in that order, for the first one that exists. This must match the field used when the input GD file(s) were created.", "AUTOMATIC", ADVANCED_OPTION);
	
	options.addUsage("");
	options.addUsage("ANNOTATE mutations in one or more GenomeDiff files. If multiple input files are provided, then also COMPARE the frequencies for identical mutations across samples.");
  options.addUsage("");
  options.addUsage("Valid output formats:");
  options.addUsage("  HTML     Descriptive table viewable in a web browser");
  options.addUsage("  GD       GenomeDiff with added annotation of mutations");
	options.addUsage("  TABLE    Comma-separated values UTF-8 text file suitable for input into Excel or R");
	options.addUsage("             Simplified output with the same content/format as HTML output.");
	options.addUsage("  TSV/CSV  Tab- or comma-separated values file suitable for input into R");
	options.addUsage("             Detailed output with all GD fields. One line per mutation per sample.");
	options.addUsage("  PHYLIP   Alignment of genotypes in PHYLIP format for phylogenetic analysis");
	options.addUsage("  FASTA    Multi-FASTA alignment of genotypes for phylogenetic analysis");
	options.addUsage("  JSON     JavaScript object notation file suitable for parsing");
	options.addUsage("");
	options.addUsage("When multiple GD files are provided, the #=TITLE metadata line in each file is used to name each sample/column. If that information is not present, then the GD file name is used (removing the *.gd suffix).");
	options.addUsage("");
	
	options.addUsage("In output, frequencies of 'D' mean that this mutation occurs within");
	options.addUsageSameLine("a region that is deleted by a different mutation in the genome in question. Frequencies");
	options.addUsageSameLine("of '?' indicate that there were not enough aligned reads to call a mutation at this position");
	options.addUsageSameLine("in the genome in question (either for or against the mutation).");
	options.addUsage("");
	options.addUsage("PHYLIP/FASTA FORMAT");
	options.addUsage("");
	options.addUsage("In this output format, each column in the genotype 'sequence' that is created corresponds to");
	options.addUsageSameLine("a unique mutational event. For SNPs the base present at that position in each genome is shown.");
	options.addUsageSameLine("For other types of mutations, 'A' is used for the ancestral allele (e.g., no transposon insertion),");
	options.addUsageSameLine("and 'T' is used for the derived allele (e.g., new transposon copy inserted).");
	options.addUsageSameLine("'N' is used when it is ambiguous as to whether the mutation occurred in the lineage");
	options.addUsageSameLine("leading to this genome, either because the position is deleted or there are not sufficient");
	options.addUsageSameLine("reads aligned to a position to call a mutation (i.e., is inside an UN region).");
	options.addUsage("");
	options.addUsage("PHYLIP/FASTA output is designed to be input into other programs to create a phylogenetic tree by parsimony methods. In these analyses, be sure that you use basic parsimony that weights all genotype sequence changes equally (e.g., same transition/transversion rates).");
	options.addUsage("");
	options.addUsage("INPUT INTO EXCEL");
	options.addUsage("");
	options.addUsage("You can load files output in the TABLE or HTML formats into Excel. For loading TABLE output, open an existing workbook and 'Import' as a Text file. You MUST choose Unicode (UTF-8) as the file origin and a comma as the delimiter. For loading HTML output, choose 'Import' as an HTML file." );
	options.addUsage("");
	options.addUsage("MUTATION EFFECTS CLASSIFICATION");
	options.addUsage("");
	options.addUsage("Each mutation has a 'genes_overlapping' list assigned based on the genes it overlaps.");
	options.addUsage("");
	options.addUsage("If the mutation affects a position within --inactivating-overlap-fraction of the length of an overlapping gene from its start, the gene is moved to the 'genes_inactivated' list if it is also a MOB, INV, INT, SNP causing a nonsense mutation, or an INS, DEL, or SUB that results in a size change that is <= the --inactivating-size-cutoff and results in a frameshift or stop codon in a protein-coding gene or an INS, DEL, SUB with a size change > the --inactivating-size-cutoff for any type of gene, even if it is in-frame in a protein-coding gene.");
	options.addUsage("");
	options.addUsage("If there are no 'genes_overlapping' or 'genes_inactivated', a mutation has a 'genes_promoter' list assigned to all genes that are <= the --promoter_cutoff bp upstream of a gene. There can be multiple qualifying genes assigned to each list, but each gene will only be in one of the lists.");

  options.processCommandArgs(argc, argv);
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}

  string genbank_field_for_seq_id = to_upper(options["genbank-field-for-seq-id"]);
  if (   (genbank_field_for_seq_id != "AUTOMATIC")
      && (genbank_field_for_seq_id != "LOCUS")
      && (genbank_field_for_seq_id != "VERSION")
      && (genbank_field_for_seq_id != "ACCESSION")
      ) {
    options.addUsage("");
    options.addUsage("Value of --genbank-field-for-seq-id must be one of the following: AUTOMATIC, LOCUS, VERSION, ACCESSION.");
    options.addUsage("");
    options.addUsage("Value provided was: " + genbank_field_for_seq_id);
    options.printUsage();
    return -1;
  }

  UserOutput uout("ANNOTATE/COMPARE");
  
  string output_format = to_upper(options["format"]);
  string output_file_name;
  if (output_format == "HTML") {
    output_file_name = "output.html";
  } else if (output_format == "GD") {
    output_file_name = "output.gd";
  } else if (output_format == "PHYLIP") {
    output_file_name = "output.phylip";
	} else if (output_format == "FASTA") {
			output_file_name = "output.fa";
	} else if (output_format == "TSV") {
		output_file_name = "output.tsv";
	} else if ( (output_format == "CSV") || (output_format == "TABLE") ) {
		output_file_name = "output.csv";
  } else if (output_format == "JSON") {
		output_file_name = "output.json";
	} else {
		options.addUsage("");
		options.addUsage("OPTION ERROR:\nUnknown output format (--format|-f) of " + output_format + " requested. Valid choices are HTML, GD, TABLE, TSV, CSV, PHYLIP, FASTA, JSON");
		options.printUsage();
		return -1;
	}
  // User setting overrules default names
  if (options.count("output")) 
    output_file_name = options["output"];
    
						 
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
	
	// Perform some pre-checks to keep things in scope
	if ( ((output_format == "PHYLIP") || (output_format == "FASTA")) && !compare_mode) {
		options.addUsage("");
		options.addUsage("OPTION ERROR:\nYou must provide more than one input GD file in PHYLIP ot FASTA mode.");
		options.printUsage();
		return -1;
	}
	if (compare_mode && options.count("preserve-evidence")) {
		options.addUsage("");
		options.addUsage("OPTION ERROR:\nThe --preserve-evidence option can only be used with one input GD file.");
		options.printUsage();
		return -1;
	}
	
	// Give a warning if add-html-fields used when it isn't necessary
	if (options.count("add-html-fields") && !((output_format == "GD") || (output_format == "JSON")) ) {
		WARN("--add-html-fields option is ignored for the selected output format " + output_format);
	}
	
	// Give a warning if add-html-fields used when it isn't necessary
	if (options.count("add-text-fields") && !((output_format == "GD") || (output_format == "JSON") || (output_format == "TSV") || (output_format == "CSV")) ) {
		WARN("--add-text-fields option is ignored for the selected output format " + output_format);
	}
	
	// Load reference files
	vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
	uout("Reading input reference sequence files") << reference_file_names << endl;
	cReferenceSequences ref_seq_info;
	ref_seq_info.LoadFiles(reference_file_names, genbank_field_for_seq_id);

	cGenomeDiff gd;
	vector<cGenomeDiff> gd_list;
	vector<string> gd_titles;
	bool polymorphism_search_found(false);
	
  if (output_format == "HTML") {
		
		load_merge_multiple_gd_files(gd, gd_list, gd_path_names, gd_titles, ref_seq_info, true, &polymorphism_search_found, compare_mode, options, uout);

		uout("Annotating mutations");
		
		ref_seq_info.annotate_mutations(gd, false, options.count("ignore-pseudogenes"), compare_mode, kBreseq_large_mutation_size_cutoff, false, from_string<double>(options["inactivating-overlap-fraction"]), from_string<uint32_t>(options["inactivating-size-cutoff"]), from_string<uint32_t>(options["promoter-distance"]) );
		
    uout("Writing output HTML file", output_file_name);
		
    Settings settings;
    // No evidence needs to be transferred to options and initialized correctly within breseq
    settings.no_evidence = true;
    
    MutationTableOptions mt_options(settings);
		mt_options.repeat_header = from_string<int32_t>(options["repeat-header"]);

    if (polymorphism_search_found)
      mt_options.force_frequencies_for_one_reference = true;
    mt_options.one_ref_seq = ref_seq_info.size() == 1;
    mt_options.gd_name_list_ref = gd_titles;
		mt_options.force_show_sample_headers = gd_path_names.size() > 1;
    html_compare(settings, output_file_name, "Mutation Comparison", gd, mt_options);
        
  } else if (output_format == "GD") {
		
		load_merge_multiple_gd_files(gd, gd_list, gd_path_names, gd_titles, ref_seq_info, !options.count("preserve-evidence"), NULL, compare_mode, options, uout);
				
		uout("Annotating mutations");
		ref_seq_info.annotate_mutations(gd, false, options.count("ignore-pseudogenes"), compare_mode);
		
		uout("Writing output Genome Diff file", options["output"]);
		
		// Add extra HTML annotations
		if (options.count("add-html-fields")) {
			// Only defaults accessible - which include javascript output...
			Settings settings;
			MutationTableOptions mutation_table_options(settings);
			diff_entry_list_t muts = gd.mutation_list();
			for (diff_entry_list_t::iterator itr = muts.begin(); itr != muts.end(); itr ++) {
				cDiffEntry& mut = (**itr);
				add_html_fields_to_mutation(mut, mutation_table_options);
			}
		}
		
		// Add extra TEXT annotations
		if (options.count("add-text-fields")) {
			// Only defaults accessible - which include javascript output...
			Settings settings;
			MutationTableOptions mutation_table_options(settings);
			diff_entry_list_t muts = gd.mutation_list();
			for (diff_entry_list_t::iterator itr = muts.begin(); itr != muts.end(); itr ++) {
				cDiffEntry& mut = (**itr);
				add_text_fields_to_mutation(mut, mutation_table_options);
			}
		}
		
    gd.write(output_file_name);
  } else if (output_format == "PHYLIP") {
		uout("Writing output PHYLIP alignment file", options["output"]);
		
		load_merge_multiple_gd_files(gd, gd_list, gd_path_names, gd_titles, ref_seq_info, true, NULL, compare_mode, options, uout);
		gd.write_genotype_sequence_file("PHYLIP", output_file_name, gd, gd_list, ref_seq_info);
		
	} else if (output_format == "FASTA") {
		uout("Writing output PHYLIP alignment file", options["output"]);
	
		load_merge_multiple_gd_files(gd, gd_list, gd_path_names, gd_titles, ref_seq_info, true, NULL, compare_mode, options, uout);
		gd.write_genotype_sequence_file("FASTA", output_file_name, gd, gd_list, ref_seq_info);
		
	} else if (output_format == "TABLE") {

		// Force compare mode
		compare_mode = true;
		load_merge_multiple_gd_files(gd, gd_list, gd_path_names, gd_titles, ref_seq_info, true, NULL, compare_mode, options, uout);

		// Remove evidence
		gd.remove_group(cGenomeDiff::EVIDENCE);
		
		uout("Annotating mutations");
		ref_seq_info.annotate_mutations(gd, false, options.count("ignore-pseudogenes"), compare_mode);
		
		Settings settings;
		MutationTableOptions mutation_table_options(settings);

		uout("Writing output TABLE file", options["output"]);
		cGenomeDiff::write_table_file(output_file_name, ",", gd, gd_titles, mutation_table_options);
		
	} else if ( (output_format == "TSV") || (output_format == "CSV") ) {
		// Load the entire list and pass it to the writer
		
		vector<cGenomeDiff> ev_list;

		for (vector<string>::iterator it=gd_path_names.begin(); it!=gd_path_names.end(); it++) {
			uout("Reading/annotating input GD file",*it);

			cGenomeDiff this_gd(*it);
			cGenomeDiff this_ev(*it);
			
			ref_seq_info.annotate_mutations(this_gd, false, options.count("ignore-pseudogenes"), compare_mode);
			
			// Add extra TEXT annotations
			if (options.count("add-text-fields")) {
				// Only defaults accessible - which include javascript output...
				Settings settings;
				MutationTableOptions mutation_table_options(settings);
				diff_entry_list_t muts = this_gd.mutation_list();
				for (diff_entry_list_t::iterator itr = muts.begin(); itr != muts.end(); itr ++) {
					cDiffEntry& mut = (**itr);
					add_text_fields_to_mutation(mut, mutation_table_options);
				}
			}
			
			gd_list.push_back(this_gd);
			ev_list.push_back(this_ev);
		}
		

		
		if (output_format == "TSV") {
			uout("Writing output TSV file", options["output"]);
			cGenomeDiff::write_separated_values_file(output_file_name, "\t", gd_list, options.count("preserve-evidence"));
		} else if (output_format == "CSV") {
			uout("Writing output CSV file", options["output"]);
			cGenomeDiff::write_separated_values_file(output_file_name, ",", gd_list, options.count("preserve-evidence"));
		}
		
	} else if (output_format == "JSON") {
		uout("Writing output JSON file", options["output"]);
		
		load_merge_multiple_gd_files(gd, gd_list, gd_path_names, gd_titles, ref_seq_info, !options.count("preserve-evidence"), NULL, compare_mode, options, uout);
		
		uout("Annotating mutations");
		ref_seq_info.annotate_mutations(gd, false, options.count("ignore-pseudogenes"), compare_mode);
		
		uout("Writing output JSON file", output_file_name);
		
		// Add extra HTML annotations
		if (options.count("add-html-fields")) {
			// Only defaults accessible - which include javascript output...
			Settings settings;
			MutationTableOptions options(settings);
			
			diff_entry_list_t muts = gd.mutation_list();
			for (diff_entry_list_t::iterator itr = muts.begin(); itr != muts.end(); itr ++) {
				cDiffEntry& mut = (**itr);
				add_html_fields_to_mutation(mut, options);
			}
		}
		
		// Add extra TEXT annotations
		if (options.count("add-text-fields")) {
			// Only defaults accessible - which include javascript output...
			Settings settings;
			MutationTableOptions mutation_table_options(settings);
			diff_entry_list_t muts = gd.mutation_list();
			for (diff_entry_list_t::iterator itr = muts.begin(); itr != muts.end(); itr ++) {
				cDiffEntry& mut = (**itr);
				add_text_fields_to_mutation(mut, mutation_table_options);
			}
		}
		
		gd.write_json(output_file_name);
		
	}
	
  return 0;
}

int do_mutations(int argc, char* argv[])
{
	AnyOption options("gdtools MUTATIONS [-o output.gd] -r reference.gbk input.gd");
	
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("output,o", "Path to output GD file with mutations predicted from evidence.", "output.gd");
	options("reference,r", "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");

	options.addUsage("");
	options.addUsage("Predicts mutations from the evidence entries in the input GD file. Any mutation entries (three-letter code lines) already present in the input GD file will be removed.");
	
	options.processCommandArgs(argc, argv);
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
	UserOutput uout("MUTATIONS");
	
	// Check the input
	if (options.getArgc() != 1) {
		options.addUsage("");
		options.addUsage("Exactly one input Genome Diff file must be provided.");
		options.printUsage();
		return -1;
	}
	
	if (!options.count("reference")){
		options.printUsage();
		return -1;
	}
	
	cGenomeDiff gd(options.getArgv(0));
	gd.remove_group(cGenomeDiff::MUTATIONS);
	
	vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
	uout("Reading input reference sequence files") << reference_file_names << endl;
	cReferenceSequences ref_seq_info;
	ref_seq_info.LoadFiles(reference_file_names);
	
	uout("Predicting mutations");
	MutationPredictor mp(ref_seq_info);
	Settings settings;
	Summary summary;
	mp.predict(settings, summary, gd);
	
	uout("Writing output Genome Diff file", options["output"]);
	gd.write(options["output"]);
	
	return 0;
}


int do_phylogeny(int argc, char* argv[])
{
	AnyOption options("gdtools PHYLOGENY [-o phylogeny] -r reference.gbk input.1.gd [input.2.gd ... ]");
	
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("verbose,v", "produce output for each mutation counted.", TAKES_NO_ARGUMENT);
	options("output,o", "base name for output files.", "phylogeny");
	options("reference,r", "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
	//options("ignore-pseudogenes", "treats pseudogenes as normal genes for calling AA changes", TAKES_NO_ARGUMENT);
	options("missing-as-ancestral,a", "Count missing data (mutations in UN regions) as the ancestral allele rather than as an unknown allele (N).", TAKES_NO_ARGUMENT);
	options("phylogeny-aware,p", "Do not consider mutations that differ in their 'phylogeny_id' fields equivalent.", TAKES_NO_ARGUMENT);
	options("population-aware,q", "Do not consider mutations in GD files with different POPULATION metadata values equivalent.", TAKES_NO_ARGUMENT);
	options("dnapars-command,c", "Command to run for phylip 'dnapars' executable. It can include a path or just the command if the executable is located in your $PATH (Default: try 'phylip dnapars' and then 'dnapars')");

	options.addUsage("");
	options.addUsage("Uses PHYLIP to construct a phylogentic tree. If you are including an ancestor");
	options.addUsageSameLine("to root the tree, you should include it as the very first Genome Diff file.");
	options.addUsage("");
	options.addUsage("You MUST have PHYLIP installed in your $PATH such that executing 'phylip dnapars' or 'dnapars' will run this command. If you have it installed in a different way, try using the --dnapars-command to specify the executable (and if necessary, also the dnapars subcommand).");

	options.processCommandArgs(argc, argv);
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
	UserOutput uout("PHYLOGENY");
	
	// User setting overrules default names
	string output_base_name = "phylogeny";
	if (options.count("output")) {
		output_base_name = options["output"];
	}
	
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
	
	// First use merge to produce a file with a line for each mutation
	cGenomeDiff gd;
	vector<cGenomeDiff> gd_list;
	
	cGenomeDiff::sort_gd_list_by_treatment_population_time(gd_list);

	bool polymorphisms_found = false; // handled putting in the polymorphism column if only one file provided
	for (uint32_t i = 0; i < gd_path_names.size(); i++){
		uout("Reading input GD file",gd_path_names[i]);
		cGenomeDiff single_gd(gd_path_names[i]);
		gd_list.push_back(single_gd);
	}    

	uint32_t file_num = 1;
	vector<string> title_list;
	for (vector<cGenomeDiff>::iterator it=gd_list.begin(); it!= gd_list.end(); it++) {
		cGenomeDiff& single_gd = *it;
		gd.merge(single_gd, false, options.count("phylogeny-aware"), options.count("population-aware"));
		title_list.push_back(single_gd.get_title());
		single_gd.set_title("_" + to_string<uint32_t>(file_num++) + "_");
	}

	gd.sort();
	
	uout("Tabulating mutation frequencies across samples");

	// Then add frequency columns for all genome diffs
	vector<string> dummy_title_list;
	cGenomeDiff::tabulate_mutation_frequencies_from_multiple_gds(gd, gd_list, dummy_title_list, options.count("phylogeny-aware"), options.count("population-aware"));

	// Save merged GD file
	string merged_gd_file_name =  output_base_name + ".merged.gd";
	gd.write(merged_gd_file_name);
	
	vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
	uout("Reading input reference sequence files") << reference_file_names << endl;
	cReferenceSequences ref_seq_info;
	ref_seq_info.LoadFiles(reference_file_names);
	
	//uout("Annotating mutations");
	//ref_seq_info.annotate_mutations(gd, true);
	
	string phylip_input_file_name = output_base_name + ".genotypes.txt";
	uout("Writing output PHYLIP alignment file", phylip_input_file_name);
	gd.write_genotype_sequence_file("PHYLIP", phylip_input_file_name, gd, gd_list, ref_seq_info, options.count("missing-as-ancestral"));
	string phylip_script_file_name = output_base_name + ".phylip.commands.txt";
	ofstream phylip_script(phylip_script_file_name.c_str());
	phylip_script << phylip_input_file_name << endl;
	phylip_script << "V" << endl; // print only one tree!
	phylip_script << "1" << endl;
	phylip_script << "4" << endl; // print steps
	phylip_script << "Y" << endl;
	
	if (file_exists("outtree")) remove_file("outtree");
	if (file_exists("outfile")) remove_file("outfile");
	
	// Figure out what executable works or is specified for dnapars
	string dnapars_command;
	if (options.count("dnapars-command")) {
		dnapars_command = options["dnapars-command"];
		uout("DNAPARS command specified as: ", dnapars_command);
	} else {
		// See if `phylip` is installed
		dnapars_command = SYSTEM_CAPTURE("which phylip", true);
		if (dnapars_command.size() > 0) {
			// Need to add the subcommand here
			dnapars_command += " dnapars";
		}
		else {
			// See if `dnapars` is installed
			dnapars_command = SYSTEM_CAPTURE("which dnapars", true);
		}
		// Didn't find one? Bail
		if (dnapars_command.size() == 0) {
			ERROR("Could not find 'phylip' or 'dnapars' command in $PATH");
		}
		
		uout("DNAPARS command found as: ", dnapars_command);
	}
	
	uout("Running DNAPARS on", phylip_input_file_name);
	SYSTEM(dnapars_command + " < " + phylip_script_file_name, false, false, false);
	
	string phylip_original_tree_file_name = "outtree";
	string phylip_renamed_tree_file_name = output_base_name + ".tre";
	ofstream renamed_tree(phylip_renamed_tree_file_name.c_str());
	string slurped_file;
	ifstream original_tree(phylip_original_tree_file_name.c_str());
	
	// Read the entire file
	string line;
	while(original_tree){
			getline(original_tree, line);
			slurped_file += line;
	}
	
	// Replace all file names 
	file_num = 1;
	for (vector<string>::iterator it = title_list.begin(); it != title_list.end(); it++) {
			slurped_file = substitute(slurped_file, "_" + to_string<uint32_t>(file_num++) + "_", *it);
	}
	renamed_tree << slurped_file;
	
	string phylip_output_file_name = output_base_name + ".phylip.output.txt";
	
	// Changed to deleting the original tree with placeholder names
	SYSTEM("rm outtree");
	//string phylip_original_tree_save_file_name = output_base_name + ".original.phylip.tre";
	//SYSTEM("mv outtree " + phylip_original_tree_save_file_name);
	
	// Keep the phylip output though
	SYSTEM("mv outfile " + phylip_output_file_name);
	
	// Create mutation key file
	uout("Creating mutation key file");

	string mutation_key_file_name = output_base_name + ".mutation.key.txt";
	ofstream mutation_key(mutation_key_file_name.c_str());
	diff_entry_list_t mut_list = gd.mutation_list();
	uint32_t i=0;
	for(diff_entry_list_t::iterator it = mut_list.begin(); it != mut_list.end(); it++) {
			cDiffEntry& mut = **it;
			mutation_key << to_string(mut._type) + "-" + mut[POSITION] + "-" + mut[GENE_NAME] << endl;
	}
	
	// Create sample key file
	string sample_key_file_name = output_base_name + ".sample.key.txt";
	ofstream sample_key(sample_key_file_name.c_str());
	i=1;
	
	set<string> used_names;
	
	for(vector<string>::iterator it = title_list.begin(); it != title_list.end(); it++) {
		if (used_names.count(*it)) {
			WARN("Tree will have duplicate taxon names due to repeated GD title: " + *it);
		}
		used_names.insert(*it);
		sample_key << "_" + to_string<uint32_t>(i++) + "_" << "\t" << *it << endl;
	}
	
	return 0;
}


int do_count(int argc, char* argv[])
{
  AnyOption options("gdtools COUNT [-o count.csv] -r reference.gbk input.1.gd [input.2.gd ... ]");
  options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("verbose,v", "produce output for each mutation counted.", TAKES_NO_ARGUMENT);
	options("reference,r", "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
	options("output,o", "path to output CSV file with count data.", "count.csv");
	options("detailed-output,d", "path to optional output tab-delimited file with detailed information about all mutations (Default = OFF)");
	options("calculate-genome-size,s", "use APPLY to calculate final genome sizes", TAKES_NO_ARGUMENT);
	options("base-substitution-statistics,b", "calculate detailed base substitution statistics", TAKES_NO_ARGUMENT);
	options("count-polymorphisms,p", "count polymorphic mutations (those with frequencies < 1). (Default = FALSE)", TAKES_NO_ARGUMENT);

	options.addUsage("");
	options.addUsage("Counts the numbers of mutations and other statistics for each input GenomeDiff file.");
  options.addUsage("");
  options.addUsage("In the output \"small\" mutations are ≤ " + to_string<int32_t>(kBreseq_large_mutation_size_cutoff) + " bp. \"large\" mutations are >" + to_string<int32_t>(kBreseq_large_mutation_size_cutoff) + " bp");
	
	options.addUsage("For base substitutions overlapping multiple genes, the mutation is counted as the most significant effect it has on any of the genes it overlaps with the precedence: nonsense > nonsynonymous > synonymous.");
	
  options.processCommandArgs(argc, argv);
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
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
    ref_seq_info.annotate_mutations(gd, true, options.count("ignore-pseudogenes"), false, kBreseq_large_mutation_size_cutoff, options.count("verbose"));
    genome_diffs.push_back(gd);
  }
	
	

	
	string detailed_output_file_name = options.count("detailed-output") ? options["detailed-output"] : "";
  MutationCountFile(
                    ref_seq_info, 
                    genome_diffs, 
                    output_file_name,
										detailed_output_file_name,
                    options.count("base-substitution-statistics"),
										options.count("count-polymorphisms"),
										options.count("calculate-genome-size"),
                    options.count("verbose")
                    );
  
  return 0;
}

int do_normalize_gd(int argc, char* argv[])
{
  AnyOption options("gdtools NORMALIZE [-o output.gd] -r reference.gbk input.gd");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("verbose,v"         , "verbose mode (flag)", TAKES_NO_ARGUMENT);
	options("output,o"          , "output Genome Diff file.", "output.gd");
	options("reference,r"       , "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
	options("reassign-ids,s"    , "reassign ids to lowest numbers possible.", TAKES_NO_ARGUMENT);
	options("repeat-adjacent,a" , "mark repeat-region adjacent, mediated, and between mutations.", TAKES_NO_ARGUMENT);
	options("dont-check-apply,x" , "skip step that checks consistency of normalize using APPLY.", TAKES_NO_ARGUMENT);

	const int32_t kDistanceToRepeat = 20;
	
  options.addUsage("");
  options.addUsage("Creates a GD file of mutations that have been normalized to the input reference files. ");
	options.addUsage("");
  options.addUsage("This process involves (1) converting AMP mutations of ≤50 bp to indels, ");
  options.addUsageSameLine("(2) shifting INS and DEL mutations to the highest coordinates possible, (3)");
  options.addUsageSameLine("adding repeat_seq, repeat_length, repeat_ref_copies, and repeat_new_copies fields");
  options.addUsageSameLine("for INS and DEL mutations that are in tandem sequence repeats of ≥5 bases in the reference sequence,");
  options.addUsageSameLine("and (4) flagging SNP, INS, or DEL mutations with sizes ≤50 bp that are within 20 bp of the ends of");
  options.addUsageSameLine("annotated mobile_element copies in the reference genome with the field mobile_element_adjacent=1");
  options.addUsage("");
	options.addUsage("Optionally, assigns 'adjacent', 'mediated', or 'between' tags to mutations within");
	options.addUsageSameLine(to_string(kDistanceToRepeat) + " bp of annotated repeat regions");
	options.addUsageSameLine("to indicate these may be hotspots that experience elevated mutation rates. (They will be counted separately");
	options.addUsageSameLine("from other mutations in gdtools COUNT). This process preserves any previous versions of these tags.");
	options.addUsageSameLine("DEL mutations with a size < " + to_string(kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation) + " bp near ");
	options.addUsageSameLine("repeat_regions are treated as 'adjacent' rather than 'mediated'." );
	options.addUsage("");
	options.addUsage("Important: the 'adjacent', 'mediated', or 'between' tags will be added based on the current genome sequence");
	options.addUsageSameLine("when this mutation is applied. So they will consider new mobile element insertions that are applied before a");
	options.addUsageSameLine("mutation, but not those that are applied after, for example. The 'before' tag can be used to enforce a specific ordering");
	options.addUsageSameLine("of these mutations so that the mobile element insertion occurs first.");
	options.addUsage("");
	options.addUsage("The coordinates of DEL mutations that have 'mediated' or 'between' tags are NOT shifted, so that");
	options.addUsageSameLine("they will remain adjacent to any relevant mobile elements or repeats.");
	options.addUsage("");
	options.addUsage("Any mutations including 'no_normalize=1' in their definition will not be normalized.");

	options.processCommandArgs(argc, argv);

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
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

    
  vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
  uout("Reading input reference sequence files.") << reference_file_names << endl;
  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(reference_file_names);

  cGenomeDiff gd(input);
	
  Settings settings;
	cReferenceSequences new_ref_seq_info;

	// To annotate repeat adjacent mutations... we have to apply
	if (options.count("repeat-adjacent")) {
		uout("Annotating repeat_region mediated, between, and adjacent mutations.");
		
		// load new copy since apply will change things
		cGenomeDiff apply_gd(input);
		
		new_ref_seq_info = cReferenceSequences::deep_copy(ref_seq_info);
		apply_gd.apply_to_sequences(ref_seq_info, new_ref_seq_info, options.count("verbose"), kDistanceToRepeat, settings.size_cutoff_AMP_becomes_INS_DEL_mutation);
		
		
		// Now transfer between and mediated tags to ones with the same IDs
		diff_entry_list_t gd_mutation_list = gd.mutation_list();
		vector<string> transfer_list = make_vector<string>("between")("mediated")("adjacent");
		
		for(diff_entry_list_t::iterator gd_it = gd_mutation_list.begin(); gd_it != gd_mutation_list.end(); gd_it++) {
			
			diff_entry_ptr_t gd_mut = *gd_it;
			diff_entry_ptr_t apply_mut_p = apply_gd.find_by_id(gd_mut->_id);
			
			// Leftover debug code
			ASSERT(apply_mut_p->_id == gd_mut->_id, "Mutation lists not identical.");
			
			for(vector<string>::iterator key_it=transfer_list.begin(); key_it != transfer_list.end(); key_it++) {
				gd_mut->erase(*key_it);
				if (apply_mut_p->entry_exists(*key_it))
					(*gd_mut)[*key_it] = (*apply_mut_p)[*key_it];
			}
		}
	} else if (!options.count("dont-check-apply")) {
		
		new_ref_seq_info = cReferenceSequences::deep_copy(ref_seq_info);
		cGenomeDiff apply_gd(input);
		apply_gd.apply_to_sequences(ref_seq_info, new_ref_seq_info, false, kDistanceToRepeat, settings.size_cutoff_AMP_becomes_INS_DEL_mutation);
	}
	
	uout("Normalizing mutations.");
	gd.normalize_mutations(ref_seq_info, settings);
	
	if (options.count("reassign-ids")) {
		uout("Reassigning mutation and evidence ids.");
		gd.reassign_unique_ids();
	}
	
	if (!options.count("dont-check-apply")) {
		uout("Using APPLY to check that normalization didn't change the mutated sequence.");
		cReferenceSequences	verify_ref_seq_info = cReferenceSequences::deep_copy(ref_seq_info);
		cGenomeDiff verify_gd(input); // must load new copy or positions will be shifted by apply_to_sequences
		verify_gd.apply_to_sequences(ref_seq_info, verify_ref_seq_info, false, kDistanceToRepeat, settings.size_cutoff_AMP_becomes_INS_DEL_mutation);
		
		vector<string> seq_ids = verify_ref_seq_info.seq_ids();
		vector<string> new_seq_ids = new_ref_seq_info.seq_ids();
		
		for (vector<string>::const_iterator it = seq_ids.begin(); it != seq_ids.end(); it++)
		{
			if (new_ref_seq_info[*it].m_fasta_sequence.get_sequence() != verify_ref_seq_info[*it].m_fasta_sequence.get_sequence()) {
				WARN("Failed APPLY test. Discrepancies between sequences produced before and after NORMALIZE. Check ordering of mutations.");
			}
		}
	}
	
	uout("Writing output Genome Diff file", options["output"]);
	gd.write(options["output"]);

  return 0;
}

int do_remove_gd(int argc, char* argv[])
{
  AnyOption options("gdtools FILTER/REMOVE [-o output.gd] -c condition1 [-c condition2] [-m SNP] input.gd");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("output,o", "Output Genome Diff file.", "output.gd");
	options("condition,c", "Condition for removing entries from the input Genome Diff file. Enclose the value of this parameter in quotes if it includes spaces, e.g. -c \"frequency <= 0.05\". Both field names and field values are case-sensitive. You may include multiple conditions on the same command line. Only entries that satisfy ALL conditions will be removed.");
	options("preserve-evidence,e", "By default evidence items with two-letter codes are removed (RA, JC, MC, ...). Supply this option to retain them. ", TAKES_NO_ARGUMENT);
	options("renumber,n","Renumber IDs of entries that remain after performing filtering.", TAKES_NO_ARGUMENT);
	options("comment,m","Instead of removing filtered entries, comment them out in the output GD file.", TAKES_NO_ARGUMENT);
	options("verbose,v","Print information about why entries were removed to the console.", TAKES_NO_ARGUMENT);

  options.processCommandArgs(argc, argv);

  options.addUsage("");
  options.addUsage("Removes entries from a GD file for which ALL of the provided conditions evaluate to true.");
	options.addUsage("");
	options.addUsage("Operators recognized in conditions are: ==, >, >=, >, <=, !=.");
	options.addUsage("");
	options.addUsage("Use the special value UNDEFINED to match that a field does not exist for the given entry.");
	options.addUsage("");
  options.addUsage("Use the special field ASSIGNED to filter evidence entries based on whether any mutation references them in the parent-id field. This field has a value of 0 (for false) or 1 (for true)");
	options.addUsage("");
	options.addUsage("USAGE EXAMPLES");
	options.addUsage("");
	options.addUsage("Remove all but insertion and deletion mutations");
	options.addUsage("  gdtools REMOVE -c \"type != INS\" -c \"type != DEL\" ...");
	options.addUsage("");
	options.addUsage("Remove all but synonymous mutations (must run gd ANNOTATE first to add field)");
	options.addUsage("  gdtools ANNOTATE -f GD -r reference.gbk -o annotated.gd input.gd");
	options.addUsage("  gdtools REMOVE -c \"snp_type != synonymous\" -o final.gd annotated.gd");
	options.addUsage("");
	options.addUsage("Remove all but unassigned JC evidence (takes two commands). Needs -e option.");
	options.addUsage("  gdtools REMOVE -e -c \"ASSIGNED=1\" -o intermediate.gd input.gd");
	options.addUsage("  gdtools REMOVE -e -c \"type != JC\" -o final.gd intermediate.gd");
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
  UserOutput uout("REMOVE");
  
  if (options.getArgc() != 1) {
    options.addUsage("");
    options.addUsage("Provide a single Genome Diff input file.");
    options.printUsage();
    return -1;
  }

  if (!options.count("condition")) {
    options.addUsage("");
    options.addUsage("You must supply at least one condition argument.");
    options.printUsage();
    return -1;
  }
  vector<string> filters = from_string<vector<string> >(options["condition"]);

  assert(filters.size());

  uout("Reading input GD file", options.getArgv(0));
  cGenomeDiff gd(options.getArgv(0));

	// Remove evidence unless requesting to keep it
	if (!options.count("preserve-evidence")) {
		gd.remove_group(cGenomeDiff::EVIDENCE);
		//this_gd.remove_group(cGenomeDiff::VALIDATION);
	}
	
	// Load the normal conditions
	vector<string> filter_keys;
	vector<size_t> filter_evals;
	vector<string> filter_values;
	
	const vector<string> evals = make_vector<string>("==")("!=")("<=")(">=")("<")(">");
	for (vector<string>:: const_iterator jt = filters.begin(); jt != filters.end(); ++jt) {
			
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
				
				filter_keys.push_back(key);
				filter_evals.push_back(eval);
				filter_values.push_back(value);
				
				break;
			}
		}
	
		ASSERT(key.size() && value.size(), "Error in format of filter: " + filter);
	}
	
	// Add the special "assigned" field
	diff_entry_list_t entries = gd.get_list();

	// Create a map of all used evidence
	map<string,bool> used_as_evidence;
	for (diff_entry_list_t::iterator it = entries.begin(); it != entries.end(); it++) {
		cDiffEntry& de = **it;
		for (vector<string>::iterator ev_it = de._evidence.begin(); ev_it != de._evidence.end(); ev_it++) {
			used_as_evidence[*ev_it] = true;
		}
	}
		
	// Assign based on presence in map
	for (diff_entry_list_t::iterator it = entries.begin(); it != entries.end(); it++) {
		cDiffEntry& de = **it;
		if (used_as_evidence.count(de._id)) {
			de["ASSIGNED"] = "1";
		} else {
			de["ASSIGNED"] = "0";
		}
	}
	
	uout("Filtering GD entries");
	
	// Perform filtering
	cGenomeDiff new_gd;
	diff_entry_list_t* mutable_entry_list = gd.get_mutable_list_ptr();
	diff_entry_list_t::iterator it = mutable_entry_list->begin();
	bool advance_iterator;
  while (it != mutable_entry_list->end()) {
		advance_iterator = true;
    cDiffEntry& de = **it;

    vector<string> reasons;

		for (size_t i=0; i<filter_evals.size(); i++) {
		
			bool is_filtered = false;
			string& key = filter_keys[i];
			size_t eval = filter_evals[i];
			string& value = filter_values[i];
			
			//Note special case for matching UNDEFINED values
			string test_string = de.count(key) ? de[key] : "UNDEFINED";
			
			// special case of frequency missing if 100%
			if ( (key==FREQUENCY) && (test_string=="UNDEFINED") ) {
					test_string = "1";
			}
			
			// Special cases for non-map entries
			if (key=="type") {
				test_string = to_string(de._type);
			}
			
			// UNDEFINED only works with == and !=
			if (test_string == "UNDEFINED") {
				switch(eval)
				{
					case 0: if (test_string == value) is_filtered = true; break;
					case 1: if (test_string != value) is_filtered = true; break;
					
					default:
						ERROR("UNDEFINED can only be used with the == and !- comparison operators");
				};
			}
			// Numeric
			else if (value.find_first_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ") == string::npos) {
				switch(eval)
				{
					case 0: if (from_string<float>(test_string) == from_string<float>(value)) is_filtered = true; break;
					case 1: if (from_string<float>(test_string) != from_string<float>(value)) is_filtered = true; break;
					case 2: if (from_string<float>(test_string) <= from_string<float>(value)) is_filtered = true; break;
					case 3: if (from_string<float>(test_string) >= from_string<float>(value)) is_filtered = true; break;
					case 4: if (from_string<float>(test_string) <  from_string<float>(value)) is_filtered = true; break;
					case 5: if (from_string<float>(test_string) >  from_string<float>(value)) is_filtered = true; break;
				};
			}
			// String
			else {
				switch(eval)
				{
					case 0: if (test_string == value) is_filtered = true; break;
					case 1: if (test_string != value) is_filtered = true; break;
					case 2: if (test_string <= value) is_filtered = true; break;
					case 3: if (test_string >= value) is_filtered = true; break;
					case 4: if (test_string <  value) is_filtered = true; break;
					case 5: if (test_string >  value) is_filtered = true; break;
				};
			}
			if (is_filtered) {
				reasons.resize(reasons.size() + 1);
				sprintf(reasons.back(), "%s %s %s", key.c_str(), evals[eval].c_str(), value.c_str());
			}
				
    }
		
    if (reasons.size() == filters.size()) {
			if (options.count("verbose")) {
				cout << "Removed [" << join(reasons, ", ") << "]: " << de.as_string() << endl;
			}
      de["filtered"] = join(reasons, ", ");
			
			// Comment out or delete
			if (options.count("comment")) {
				de["comment_out"] = "true";
			} else {
				it = gd.remove(it);
				advance_iterator = false;
			}
    }
		
		if (advance_iterator) it++;
  }

	// Remove the special fields
	for (diff_entry_list_t::iterator it = entries.begin(); it != entries.end(); it++) {
		cDiffEntry& de = **it;
		de.erase("ASSIGNED");
	}
	
	// Fix the ids - if requested
	if(options.count("renumber"))
	{
		uout("Renumbering GD entries");

		cGenomeDiff gd2;
		gd2.metadata = gd.metadata;
		gd2.merge(gd, true);
		gd = gd2;
	}
	
  uout("Writing output GD file", options["output"]);
  gd.write(options["output"]);
  return 0;
}

int do_mummer2mask(int argc, char* argv[])
{
	AnyOption options("gdtools MUMMER2MASK [-o output.gd -p 36] -r reference.fna input.coords");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("reference,r","File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
	options("output,o", "Output Genome Diff file.", "output.gd");
	options("padding,p", "Additional padding to add to each end of every MASK region.", 0);
	options("merge,g", "Merge regions if they are within this distance (before adding padding).", 0);
	options("minimum,m", "Minimum size of a region after padding and merging to remain in the MASK list.", 0);
	options.addUsage("");
	options.addUsage("The input file of repeats should be generated using a MUMmer command like this:");
	options.addUsage("  mummer -maxmatch -b -c -l 36 reference.fna reference.fna > input.coords");
	options.addUsage("");
	options.addUsage("Creates a GD file with MASK entries from the output of MUMmer on the reference genome.");
	options.addUsage("This file can be used to remove mutation predictions that overlap the masked regions");
	options.addUsage("using the gdtools MASK command.");

	options.processCommandArgs(argc, argv);

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
	if (options.getArgc() != 1) {
		options.addUsage("");
		options.addUsage("Provide exactly one input.coords file generated by MUMmer.");
		options.printUsage();
		return -1;
	}
	
	if (!options.count("reference")) {
		options.addUsage("");
		options.addUsage("You must supply the --reference|-r option.");
		options.printUsage();
		return -1;
	}
	
	UserOutput uout("MUMMER2MASK");
	
	cReferenceSequences ref_seq_info;
	ref_seq_info.LoadFiles(from_string<vector<string> >(options["reference"]));
	
	cGenomeDiff gd;
	cFlaggedRegions flagged_regions;
	
	uout << "Reading MUMmer file: " << options.getArgv(0) << endl;
	flagged_regions.read_mummer(options.getArgv(0), ref_seq_info);
	
	uout << "Merging regions" << endl;
	flagged_regions.merge_within_distance(from_string<uint32_t>(options["merge"]));
	
	uout << "Adding padding to ends of regions" << endl;
	flagged_regions.add_padding_to_ends(from_string<int32_t>(options["padding"]), from_string<int32_t>(options["minimum"]));
	
	// Flagged regions has to keep track of the seq_id
	
	for (cReferenceSequences::iterator its=ref_seq_info.begin(); its != ref_seq_info.end(); its++) {
	
		cFlaggedRegions::regions_t regions = flagged_regions.all_regions(its->m_seq_id);
		for(cFlaggedRegions::regions_t::iterator it=regions.begin(); it!=regions.end(); it++ ) {
			cDiffEntry mask_entry(MASK);
			mask_entry[SEQ_ID] = its->m_seq_id;
			mask_entry[POSITION] = to_string<uint32_t>(it->first);
			mask_entry[SIZE] = to_string<uint32_t>(it->second - it->first + 1);
			gd.add(mask_entry);
		}
	}
	uout << "Writing GD file: " << options["output"] << endl;
	gd.write(options["output"]);
	
	return 0;
}


int do_mask_gd(int argc, char* argv[])
{
	AnyOption options("gdtools MASK  [-o output.gd] input.gd mask.gd");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("small,s", "Mask only 'small' mutations defined as: all SNP mutations; INS, DEL, and SUB mutations with sizes ≤ 20 bp; and all INS and DEL mutations causing expansion or contraction of simple sequence repeats (SSRs) with at least two repeats of the same unit of one to several bp and a total length of 5 bp in the reference genome. If these mutations are marked as 'mediated' or 'between' repeats, then they are NOT removed.", TAKES_NO_ARGUMENT);
	options("output,o", "Output Genome Diff file.", "output.gd");
	options("verbose,v","Verbose mode", TAKES_NO_ARGUMENT);

	options.addUsage("");
	options.addUsage("Creates a GD file where mutations in the input GD that are located within certain regions of the reference genome are removed. These regions are defined as MASK entries in the mask GD file. Mutations that overlap only masked reference bases (and therefore do not overlap any unmasked bases) are removed.");
	
	options.processCommandArgs(argc, argv);
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
	UserOutput uout("MASK");
	
	if (options.getArgc() != 2) {
		options.addUsage("");
		options.addUsage("Provide two Genome Diff input files as input.");
		options.printUsage();
		return -1;
	}
	
	uout("Reading input GD file") << options.getArgv(0) << endl;
	cGenomeDiff gd(options.getArgv(0));
	
	uout("Reading mask GD file") << options.getArgv(1) << endl;
	cGenomeDiff mask_gd(options.getArgv(1));
	
	cGenomeDiff new_gd(gd);
	new_gd.metadata = gd.metadata; // copy all of the important information
	
	// Mask mutations
	new_gd.mask_mutations(mask_gd, options.count("small"), options.count("verbose"));
	
	uout("Writing output GD file", options["output"]);
	new_gd.write(options["output"]);
	return 0;
}

int do_read_count(int argc, char* argv[])
{
	
	AnyOption options("gdtools READ-COUNT [-i input.gd -o output.gd] input1.fastq input2.fastq ... ");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("input,i", "Input GenomeDiff file", "input.gd");
	options("output,o", "Output GenomeDiff file containing read counts", "output.gd");
	options.addUsage("");
	options.addUsage("Counts the number of reads and bases in the input FASTQ files and adds them as ORIGINAL-READS and ORIGINAL-BASES header information in the output GenomeDiff file.");
	
	options.processCommandArgs( argc,argv);
	
	if (options.count("help")) {
		options.printAdvancedUsage();
		return -1;
	}
	
	if (options.getArgc() == 0) {
		options.addUsage("");
		options.addUsage("Must input at least one FASTQ file as an unnamed argument.");
		options.printUsage();
		return -1;
	}
	
	UserOutput uout("READ-COUNT");

	
	cGenomeDiff gd(options["input"]);
	uout << "Reading input GenomeDiff file: " << options["input"] << endl;
	uint64_t total_reads(0);
	uint64_t total_bases(0);
	
	vector<string> fastq_file_names;
	for (int32_t i = 0; i < options.getArgc(); i++)
	{
		uout << "Reading input FASTQ file: " << options.getArgv(i) << endl;;

		uint64_t num_original_reads, num_original_bases;
		uint32_t read_length_min, read_length_max;
		uint8_t min_quality_score, max_quality_score;
		string quality_format = cFastqQualityConverter::predict_fastq_file_format(options.getArgv(i), num_original_reads, num_original_bases, read_length_min, read_length_max, min_quality_score, max_quality_score);
		
		uout << "  Reads: " << num_original_reads << endl;
		uout << "  Bases: " << num_original_bases << endl;
		
		total_reads += num_original_reads;
		total_bases += num_original_bases;
	}
	
	gd.add_breseq_data("ORIGINAL-READS", to_string<uint64_t>(total_reads));
	gd.add_breseq_data("ORIGINAL-BASES", to_string<uint64_t>(total_bases));
	
	uout << "Writing output GenomeDiff file: " << options["output"] << endl;
	gd.write(options["output"]);
	
	uout << "SUCCESSFULLY COMPLETED" << endl;;
	return 0;
}



int do_simulate_mutations(int argc, char *argv[])
{
  AnyOption options("Usage: gdtools SIMULATE-MUTATIONS [-n 100] -r <reference> -o <output.gd> -t <type>");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("reference,r","File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
  options("output,o","Output file");
  options("type,t","Type of mutation to generate");
  options("number,n","Number of mutations to generate", static_cast<uint32_t>(1000));
  options("exclude,e","Exclusion file containing regions no mutations will overlap, usually repeats");
  options("buffer,b","Buffer distance between mutations and excluded intervals", static_cast<uint32_t>(50));
  options("seq,s","Reference sequence id to use from reference file");  
  options("seed","Seed for the random number generator");
  options("verbose,v","Verbose mode", TAKES_NO_ARGUMENT);
  
  options.addUsage("");
  options.addUsage("Using --reference, this command will generate a --number of");
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

	options.processCommandArgs(argc, argv);

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
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
      ref_seq_info,
      options.count("verbose")
      );

  gd.mutations_to_evidence(ref_seq_info, false);

  gd.reassign_unique_ids();
  
  gd.write(options["output"]);
  
  return 0;
}

int do_mutations_to_evidence(int argc, char *argv[])
{
  AnyOption options("Usage: gdtools MUTATIONS-TO-EVIDENCE -r <reference> -o <output.gd> input.gd");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("reference,r","File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");  
  options("output,o","Output file");
  options("verbose,v","Verbose Mode (Flag)", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
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
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
  options("output,o",      "output GD file");
  options("reference,r",   "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
  options("tag,t",         "header tag to add to GenomeDiff file, input as <key>=<value> will produce #=<key> <value>");
  options("input,i",       "modify input GenomeDiff file inplace");
  options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
  options.processCommandArgs(argc, argv);
  
  options.addUsage("");
  options.addUsage("Create or add '#=<TAG>' entries to the header of a GenomeDiff file,");
  options.addUsage("the -r argument will be added as #=REFSEQ and the *.fastq arguments");
  options.addUsage("will be added as #=READSEQ");

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
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

int do_reheader(int argc, char* argv[])
{
	AnyOption options("gdtools REHEADER [-o output.gd] header.gd input.gd");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("output,o", "output GD file", "output.gd");
	options("verbose,v", "verbose mode", TAKES_NO_ARGUMENT);
	options.processCommandArgs(argc, argv);
	
	options.addUsage("");
	options.addUsage("Outputs a new GD file with content from input.gd and header information from header.gd");
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
	if (options.getArgc() != 2) {
		cout << "Must supply exactly two GenomeDiff files as input." << endl;
		return -1;
	}
	
	UserOutput uout("REHEADER");
	
	cGenomeDiff header_gd(options.getArgv(0));
	cGenomeDiff input_gd(options.getArgv(1));

	input_gd.metadata = header_gd.metadata;
	
	uout("Writing output GD file", options["output"] );
	input_gd.write(options["output"]);
	
	return 0;
}

int do_download(int argc, char *argv[])
{
  stringstream ss;
  ss << "Usage: gdtools DOWNLOAD -l <user:password> -d <download_dir> -g <genome_diff_dir>\n";
  ss << "Usage: gdtools DOWNLOAD -l <user:password> -d <download_dir> <file1.gd file2.gd file3.gd ...>\n";

  AnyOption options(ss.str());
	options("help,h",            "Display detailed help message", TAKES_NO_ARGUMENT);
  options("login,l",           "Login user:password information for private server access.");
  options("download-dir,d",    "Output directory to download file to.", "02_Downloads");
  options("genome-diff-dir,g", "Directory to search for genome diff files.", "01_Data");
  options("test"           ,   "Test urls in genome diff files, doesn't download the file", TAKES_NO_ARGUMENT);
  options("reference-only",    "Only downloads the reference sequence files for this file", TAKES_NO_ARGUMENT);
	options("ungzip,z",					 "Decompress gzipped read files", TAKES_NO_ARGUMENT);

	options.addUsage("\nExamples:");
	options.addUsage("  gdtools DOWNLOAD -l john:1234 -d downloads -g data");
	options.addUsage("  gdtools DOWNLOAD -l john:1234 -d downloads 1B4.gd GRC2000.gd");

  options.processCommandArgs(argc, argv);

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
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
	
	string hostname = SYSTEM_CAPTURE("hostname -f", true);
  cout << endl << "Hostname:" << hostname << endl;
	bool onTACC = (hostname.find("tacc.utexas.edu") != string::npos);
	if (onTACC) cout << "Detected TACC system..." << endl;

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
    lookup_table["NCBI-SRA"]
        ["url_format"] = "%s";
    lookup_table["BARRICKLAB-PUBLIC"]
        ["url_format"] = "http://barricklab.org/%s";
		lookup_table["NO-DOWNLOAD"]
				["url_format"] = "no-download:%s";
		lookup_table["DEFAULT"]
				["url_format"] = "%s";
	
	if (!onTACC) {
		lookup_table["BARRICKLAB-PRIVATE"]
			["url_format"] = "scp stampede.tacc.utexas.edu:/corral-repl/utexas/breseq/%s";
	} else {
		lookup_table["BARRICKLAB-PRIVATE"]
			["url_format"] = "cp /corral-repl/utexas/breseq/%s";
	}
			
	// Old usage
	//    lookup_table["BARRICKLAB-PRIVATE"]
  //      ["url_format"] = options.count("login") ? "ftp://" + options["login"] + "@backup.barricklab.org/%s" : "";

  //File path formats.
    lookup_table["GENBANK"]
        ["file_path_format"] = download_dir + "/%s.gbk";
    lookup_table["NCBI-SRA"]
        ["file_path_format"] = download_dir + "/%s";;
    lookup_table["BARRICKLAB-PUBLIC"]
        ["file_path_format"] = download_dir + "/%s";
    lookup_table["BARRICKLAB-PRIVATE"]
        ["file_path_format"] = download_dir + "/%s";
		lookup_table["NO-DOWNLOAD"]
				["file_path_format"] = download_dir + "/%s";
		lookup_table["DEFAULT"]
				["file_path_format"] = download_dir + "/%s";
	
  /*! Rather than crash; gather [gd_file_name][reason] = error_value.
      and output at the end of the function. */
  map<string, map<string, string> > error_report;

  for (;file_names.size(); file_names.pop_front()) {
    const string &file_name = file_names.front();
    cout << endl << "Parsing file: " << file_name << endl;
    cGenomeDiff gd(file_name);
    const vector<string> &refs  = gd.metadata.ref_seqs;
    const vector<string> &reads = gd.metadata.read_seqs;
    const vector<string> &adapters = gd.metadata.adapter_seqs;
      
    list<string> seqs_kv_pairs;
    copy(refs.begin(), refs.end(), back_inserter(seqs_kv_pairs));
    if (!options.count("reference-only")) {
      copy(reads.begin(), reads.end(), back_inserter(seqs_kv_pairs));
      copy(adapters.begin(), adapters.end(), back_inserter(seqs_kv_pairs));
    }

    for (;seqs_kv_pairs.size(); seqs_kv_pairs.pop_front()) {
      const cKeyValuePair seq_kvp(seqs_kv_pairs.front(), ':');
      if (!seq_kvp.valid()) {
        error_report[file_name]["NOT_KEY_VALUE_PAIR"] = seq_kvp;
        continue;
      }

      string key = to_upper(seq_kvp.get_key());
			string value = cString(seq_kvp.get_value()).trim_ends_of('/');

			// In this case we have a URL http:// etc
      if (!lookup_table.count(key)) {
				key = "DEFAULT";
				value = seqs_kv_pairs.front();
			}
			
			if (key == "NO-DOWNLOAD") continue;

      //! Step: Get file path and check if it has already been downloaded or is empty.
      const string &base_name = cString(value).get_base_name();
      string file_path = "";
      sprintf(file_path, lookup_table[key]["file_path_format"].c_str(), base_name.c_str());
      assert(file_path.size());

      bool is_downloaded =
          ifstream(file_path.c_str()).good() && !file_empty(file_path.c_str());
			
			// Need to glob to check for matching filenames
			if (key == "NCBI-SRANCBI-SRA") {
				glob_t glob_result;
				glob((file_path + "*").c_str(),GLOB_TILDE,NULL,&glob_result);
				is_downloaded = (glob_result.gl_pathc > 0);
			}
			
      bool is_gzip =
          cString(file_path).ends_with(".gz");
			if (options.count("ungzip")) is_gzip = true;
			
      const string &gunzip_path = is_gzip ?
          cString(file_path).remove_ending(".gz") : "";

      if (options.count("ungzip")) {
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
      else if (key == "NCBI-SRA") {
        sprintf(url, url_format, value.c_str());
      }
      else if (key == "BARRICKLAB-PUBLIC") {
        sprintf(url, url_format, value.c_str());
      }
      else if (key == "BARRICKLAB-PRIVATE") {
        //ASSERT(options.count("login"), "Provide the login option (-l user:password) to access private files.");
        sprintf(url, url_format, value.c_str());
			} else if (key == "DEFAULT") { //case for normal URLs
				sprintf(url, url_format, value.c_str());
			}

			string wget_cmd;
			if (key == "BARRICKLAB-PRIVATE") {
				wget_cmd = url + " " + file_path;
				cout << wget_cmd << endl;
			} else if (key == "NCBI-SRA") {
				wget_cmd = "fastq-dump --split-files --outdir " + download_dir + " " + url;
		  } else {
				wget_cmd = options.count("test") ?
							cString("wget --spider \"%s\"", url.c_str()) :
							cString("wget -O %s \"%s\"",file_path.c_str(), url.c_str());
			}
			
				
      const bool is_url_error = SYSTEM(wget_cmd, false, false, false) != 0;
      if (is_url_error) {
        error_report[file_name]["INVALID_URL"]  = url;
        continue;
      }

      if (options.count("ungzip")) {
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
        breseq::getline(in, first_line);
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

int do_gd2coverage(int argc, char* argv[])
{
	
	AnyOption options("gdtools GD2COV [-o output -r reference.gbk] input1.gd input2.gd ... ");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("reference,r",  "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
	options("output,o","Base name of output files", "output");
	options("tile-size,t","Size of tiling for coverage graph", 1000);


	options.addUsage("");
	options.addUsage("Creates files for graphing duplicated and deleted regions across many evolved genomes.");

	options.processCommandArgs( argc,argv);

	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
	if (!options.count("reference")) {
		options.addUsage("");
		options.addUsage("You must provide a reference sequence file (-r).");
		options.printUsage();
		return -1;
	}
	
	vector<string> gd_file_names;
	for (int32_t i = 0; i < options.getArgc(); i++)
	{
		string file_name = options.getArgv(i);
		gd_file_names.push_back(file_name);
	}
	
	cGenomeDiff::GD2COV(gd_file_names, from_string<vector<string> >(options["reference"]), options["output"], from_string<int32_t>(options["tile-size"]));
	
	return 0;
}


int do_deleted_genes(int argc, char* argv[])
{
	
	AnyOption options("gdtools deleted-genes [-o output.csv -r reference.gbk] input1.gd input2.gd ... ");
	options("help,h", "Display detailed help message", TAKES_NO_ARGUMENT);
	options("reference,r",  "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)");
	options("output,o","Name of output table file", "output.csv");
	
	
	options.addUsage("");
	options.addUsage("Writes out a table of 0/1 for whether a big deletion (>50 nt) overlaps each gene.");
	
	options.processCommandArgs( argc,argv);
	
	if (options.count("help")) {
		options.printUsage();
		return -1;
	}
	
	if (!options.count("reference")) {
		options.addUsage("");
		options.addUsage("You must provide a reference sequence file (-r).");
		options.printUsage();
		return -1;
	}
		
	if (options.getArgc() == 0) {
		options.addUsage("");
		options.addUsage("You must provide at least one Genome Diff file.");
		options.printUsage();
		return -1;
	}

	// Load references
	vector<string> reference_file_names = from_string<vector<string> >(options["reference"]);
	cReferenceSequences ref_seq_info;
	ref_seq_info.LoadFiles(reference_file_names);
	
	// Handle each GD
	vector<string> gd_file_names;
	for (int32_t i = 0; i < options.getArgc(); i++)
	{
		string file_name = options.getArgv(i);
		
		cGenomeDiff gd(file_name);
		
		// Create flagged regions
		cFlaggedRegions fr;
		
		diff_entry_list_t muts = gd.mutation_list();
		for (diff_entry_list_t::iterator it = muts.begin(); it != muts.end(); it++) {
			cDiffEntry& mut = **it;
			
			if (mut._type == DEL) {
				cReferenceCoordinate start_1 = mut.get_reference_coordinate_start();
				cReferenceCoordinate end_1 = mut.get_reference_coordinate_end();
				if (end_1 - start_1 + 1 >= 50) {
					fr.flag_region(mut[SEQ_ID], start_1.get_position(), end_1.get_position());
				}
			}
		}
		
		ofstream out(options["output"].c_str());
		
		out << join(make_vector<string>("file")("title")("gene")("locus_tag")("deleted"), ",") << endl;
		for(cReferenceSequences::iterator its=ref_seq_info.begin(); its!=ref_seq_info.end(); its++) {
			
			cAnnotatedSequence& ref_seq = *its;
				
			for (cSequenceFeatureList::iterator it = ref_seq.m_genes.begin(); it != ref_seq.m_genes.end(); ++it) {
				cSequenceFeature& feat = **it;
				cGeneFeature gene = cGeneFeature(feat); // up-cast
			
				out << gd.get_file_name();
				out << "," << gd.get_title();
				out << "," << gene.name;
				out << "," << gene.m_gff_attributes["ID"];
				
				bool is_deleted = false;
				for (cFeatureLocationList::iterator itr = feat.m_locations.begin(); itr != feat.m_locations.end(); itr++) {
					is_deleted = is_deleted || fr.is_flagged(ref_seq.m_seq_id, itr->get_start_1(), itr->get_end_1());
				}
				
				if (is_deleted) {
					out << ",1";
				} else {
					out << ",0";
				}
				out << endl;
			}
			
		}
	
	}
	
	// Output loop
	
	
	return 0;
}


                                            
int main(int argc, char* argv[]) {
	
	signal(SIGSEGV, seg_fault_handler);
	Settings::set_global_paths();
	
	//Extract the sub-command argument.
	string command;
	vector<char*> argv_new(argc);
	int argc_new = argc - 1;

	if (argc > 1) {

		command = argv[1];
		argv_new[0] = argv[0];
		for (int32_t i = 1; i < argc; i++)
			argv_new[i] = argv[i + 1];

	} else {
		Settings::command_line_run_header();
    gdtools_usage();
    return -1; 
	}

	//Pass the command to the proper handler.
	command = to_upper(command);

	// gnu standard return version string
	if ( (argc_new == 1) && (command == "--VERSION") ) {
		string github_revision_string(GITHUB_REVISION_STRING);
		if (github_revision_string.length() > 0) github_revision_string = " " + github_revision_string;
		cout << "gdtools " << VERSION << github_revision_string << endl;
		return 0;
	}
	
	//Print out our generic header.
	Settings::command_line_run_header();
	
  // Genome Diff Commands:
  if (command == "VALIDATE") {
    return do_validate(argc_new, argv_new.data());
  } else if (command == "APPLY") {
    return do_apply(argc_new, argv_new.data());    
  } else if ( (command == "COMPARE") || (command == "ANNOTATE") ) {
    return do_annotate(argc_new, argv_new.data());
  } else if (command == "CHECK") {
    return do_check(argc_new, argv_new.data());
  } else if (command == "NOT-EVIDENCE") {        //TODO merge with FILTER
    return do_not_evidence(argc_new, argv_new.data());
	} else if (command == "MUTATIONS") {
		return do_mutations(argc_new, argv_new.data());
  } else if (command == "PHYLOGENY"){
    return do_phylogeny(argc_new, argv_new.data());
  } else if (command == "COUNT") {
    return do_count(argc_new, argv_new.data());
  } else if (command == "NORMALIZE") {
    return do_normalize_gd(argc_new, argv_new.data());
  } else if (command == "FILTER") {
		return do_remove_gd(argc_new, argv_new.data());
	} else if (command == "REMOVE") {
		return do_remove_gd(argc_new, argv_new.data());
	} else if (command == "INTERSECT") {
    return do_intersection(argc_new, argv_new.data());
	} else if ( (command == "MERGE") || (command == "UNION") ) {
		return do_merge(argc_new, argv_new.data());
  } else if (command == "SUBTRACT") {
    return do_subtract(argc_new, argv_new.data());
  } else if (command == "WEIGHTS") {
    return do_weights(argc_new, argv_new.data());
	} else if (command == "CONVERT") {
		return do_convert(argc_new, argv_new.data());
  } else if (command == "GD2VCF") {
    return do_convert(argc_new, argv_new.data(), "VCF");
  } else if (command == "VCF2GD") {
    return do_convert(argc_new, argv_new.data(), "GD");
  } else if (command == "GD2GVF") {
    return do_convert(argc_new, argv_new.data(), "GVF");
  } else if(command == "HEADER"){
    return do_header(argc_new, argv_new.data());
	} else if(command == "REHEADER"){
		return do_reheader(argc_new, argv_new.data());
  } else if ((command == "RANDOM-MUTATIONS") || (command == "SIMULATE-MUTATIONS")) {
    return do_simulate_mutations(argc_new, argv_new.data());
  } else if (command == "MUTATIONS-TO-EVIDENCE") {
    return do_mutations_to_evidence(argc_new, argv_new.data());
  } else if (command == "DOWNLOAD") {
    return do_download(argc_new, argv_new.data());
	} else if (command == "GD2COV") {
		return do_gd2coverage(argc_new, argv_new.data());
	} else if (command == "DELETED-GENES") {
		return do_deleted_genes(argc_new, argv_new.data());
	} else if (command =="MUMMER2MASK") {
		return do_mummer2mask(argc_new, argv_new.data());
	} else if (command =="MASK") {
		return do_mask_gd(argc_new, argv_new.data());
	} else if (command =="READ-COUNT") {
		return do_read_count(argc_new, argv_new.data());
  } else {
    cout << "Unrecognized command: " << command << endl;
    gdtools_usage();
  }
  return 0;

}

