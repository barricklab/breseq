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

#include "libbreseq/samtools_commands.h"

using namespace std;

// Pull in the C code for these functions included
// as the original samtools source files here

extern "C" {

// samtools index
int bam_index(int argc, char *argv[]);

// samtools sort
int bam_sort(int argc, char *argv[]);

// samtools import
int main_import(int argc, char *argv[]);

// samtools faidx
int faidx_main(int argc, char *argv[]);

}


namespace breseq {


void samtools_index(const string& bam_file_name) {
  
  cli_arguments args;
  args.push_back("index");
  args.push_back(bam_file_name);
  args.push_back(bam_file_name + ".bai");
  
  int argc;
  char** argv;
  args.get_args(&argc,&argv);
  
  string command = "[samtools] " + args.as_string();
  cout << command << endl;
  int return_code = bam_index(argc, argv);
  ASSERT(return_code == 0, "Error code returned by call to command: " + to_string(return_code));
  
}

void samtools_sort(const string& unsorted_bam_file_name, const string& sorted_bam_file_name) {
  
  cli_arguments args;
  args.push_back("sort");
  args.push_back("-o");
  args.push_back(sorted_bam_file_name);
  args.push_back("-T");
  args.push_back(sorted_bam_file_name);
  args.push_back(unsorted_bam_file_name);


  int argc;
  char** argv;
  args.get_args(&argc,&argv);
  
  string command = "[samtools] " + args.as_string();
  cout << command << endl;
  int return_code = bam_sort(argc, argv);
  ASSERT(return_code == 0, "Error code returned by call to command: " + to_string(return_code));
  
}


void samtools_import(const string& faidx_file_name, const string& sam_file_name, const string& output_bam_file_name) {
  
  cli_arguments args;
  args.push_back("import");
  args.push_back(faidx_file_name);
  args.push_back(sam_file_name);
  args.push_back(output_bam_file_name);
  
  int argc;
  char** argv;
  args.get_args(&argc,&argv);

  string command = "[samtools] " + args.as_string();
  cout << command << endl;
  int return_code = main_import(argc, argv);
  ASSERT(return_code == 0, "Error code returned by call to command: " + to_string(return_code));
  
}

void samtools_faidx(const string& fasta_file_name) {
  
  cli_arguments args;
  args.push_back("faidx");
  args.push_back(fasta_file_name);
  
  int argc;
  char** argv;
  args.get_args(&argc,&argv);
  
  string command = "[samtools] " + args.as_string();
  cout << command << endl;
  int return_code = faidx_main(argc, argv);
  ASSERT(return_code == 0, "Error code returned by call to command: " + to_string(return_code));
  
}
  
}
