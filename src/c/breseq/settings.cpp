/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2008-2010 Michigan State University
 Copyright (c) 2011-2022 The University of Texas at Austin
 
 breseq is free software; you can redistribute it and/or modify it under the
 terms the GNU General Public License as published by the Free Software
 Foundation; either version 1, or (at your option) any later version.
 
 *****************************************************************************/


#include "libbreseq/settings.h"

#include "libbreseq/anyoption.h"


using namespace std;

namespace breseq
{
  const int32_t kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation = 50;
  const double kBreseq_ignore_within_this_multiple_of_average_read_length_of_contig_end = 2.0;

  // small mutation <= kBreseq_large_mutation_size_cutoff < large_mutation
  const int32_t kBreseq_large_mutation_size_cutoff = kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation;
  
  const char* kBreseqAlignmentScoreBAMTag = "X5";
  const char* kBreseqBestAlignmentScoreBAMTag = "X6";

  string Settings::global_bin_path;
  string Settings::global_program_data_path;

  //! Multithreading
  ctpl::thread_pool Settings::pool(1);
  std::mutex Settings::lock;

  string Settings::output_divider("================================================================================");
  
	void cReadFiles::Init(const vector<string>& read_file_names, bool sam_files)
	{
		//clear any existing info
		this->clear();

		uint32_t on_id = 0;
		uint32_t on_error_group = 0;
		uint32_t on_paired_end_group = 0;

    //keep track of duplicates
    set<string> used_base_names;  // just filename base
    set<string> used_original_names;  // complete path

		for (vector<string>::const_iterator it = read_file_names.begin(); it < read_file_names.end(); it++)
		{
			cReadFile rf;
			rf.m_original_file_name = *it;
      
      // Check for exact duplicate paths, all but the first will be skipped with errors.
      if (used_original_names.find(rf.m_original_file_name) != used_original_names.end()) {
        WARN("Duplicate read file path provided: " + rf.m_original_file_name + "\nThis file will only be used as input once.");
        continue;
      }
        
      used_original_names.insert(rf.m_original_file_name);
      
			rf.m_paired_end_group = on_paired_end_group++;
			rf.m_error_group = on_error_group++;
			rf.m_id = on_id++;

			// create base name
			rf.m_base_name = rf.m_original_file_name;
			// - beginning path
			size_t pos = rf.m_base_name.rfind("/");
			if (pos != string::npos) rf.m_base_name.erase(0, pos + 1);
			// - trailing .fastq, must be at the end of the sequence
      if (!sam_files) {
        pos = rf.m_base_name.rfind(".gz");
        if ((pos != string::npos) && (pos == rf.m_base_name.size() - 3))
          rf.m_base_name.erase(pos, 3);
        pos = rf.m_base_name.rfind(".fastq");
        if ((pos != string::npos) && (pos == rf.m_base_name.size() - 6))
          rf.m_base_name.erase(pos, 6);
      } else {
        pos = rf.m_base_name.rfind(".sam");
        if ((pos != string::npos) && (pos == rf.m_base_name.size() - 4))
          rf.m_base_name.erase(pos, 4);
      }
      
      // This code gives a new base name when there are two files that have the exact same
      // filename but are in different directories. Ex: path1/fastq1.gz and path2/fastq1.gz
      string original_base_name = rf.m_base_name;
      uint32_t duplicate_index = 0;
      while ( used_base_names.find(rf.m_base_name) != used_base_names.end()) {
        duplicate_index++;
        rf.m_base_name = original_base_name + "_" + to_string(duplicate_index);
      }
      
      // Warn if there is a duplicate so they know what the new name means.
      if (rf.m_base_name != original_base_name) {
        WARN("Duplicate read base file name will be renamed in breseq output as follows:\n" + original_base_name + " => " + rf.m_base_name);
      }
      
      used_base_names.insert(rf.m_base_name);
			
      // set up the map for converting base names to fastq file names to be used
      read_file_to_fastq_file_name_map[rf.m_base_name] = rf.m_original_file_name;
      
			this->push_back(rf);
		}
	}
  
  string cReadFiles::base_name_to_read_file_name(const string& base_name)
  {    
    if (read_file_to_converted_fastq_file_name_map.count(base_name)) 
    {
      return read_file_to_converted_fastq_file_name_map[base_name];
    }
    
    assert(read_file_to_fastq_file_name_map.count(base_name));
    
    return read_file_to_fastq_file_name_map[base_name];
  }
  
  cReferenceSequenceSettings::cReferenceSequenceSettings(
                            cReferenceSequences& ref_seq_info,
                            vector<string>& reference_file_names,
                            vector<string>& contig_reference_file_names,
                            vector<string>& junction_only_reference_file_names  
                            )
  {    
    m_reference_file_names = reference_file_names;
    m_contig_reference_file_names = contig_reference_file_names;
    m_junction_only_reference_file_names = junction_only_reference_file_names;
    
    // Initialize the junction_only and other sets
    set<string> junction_only_file_name_set(m_junction_only_reference_file_names.begin(), m_junction_only_reference_file_names.end());
    set<string> contig_file_name_set(m_contig_reference_file_names.begin(), m_contig_reference_file_names.end());

    for(cReferenceSequences::iterator it=ref_seq_info.begin(); it!=ref_seq_info.end(); it++) {
      if ( junction_only_file_name_set.count(it->get_file_name()) ) {
        m_junction_only_seq_id_set.insert(it->m_seq_id);
      } else {
        m_call_mutations_seq_id_set.insert(it->m_seq_id);
        if ( contig_file_name_set.count(it->get_file_name()) ) {
          m_contig_seq_id_set.insert(it->m_seq_id);
          it->m_is_contig = true;
        }
      }
    }
    
    // Map by filename (only for contig_reference_file_names to id)
    map<string,uint32_t> contig_to_id_map;
    uint32_t on_coverage_group_index = 0; 

    for(cReferenceSequences::iterator it=ref_seq_info.begin(); it!=ref_seq_info.end(); it++) {
    
      // contig_reference sequence file -> lump with previous or start new
      if ( contig_file_name_set.count(it->get_file_name()) ) {
        uint32_t this_coverage_group_index;
        if (contig_to_id_map.count(it->get_file_name())) {
          this_coverage_group_index = contig_to_id_map[it->get_file_name()];
        } else {
          this_coverage_group_index = m_seq_ids_by_coverage_group.size();
          contig_to_id_map[it->get_file_name()] = this_coverage_group_index;
          m_seq_ids_by_coverage_group.push_back(vector<string>());
        }
        m_seq_ids_by_coverage_group[this_coverage_group_index].push_back(it->m_seq_id);
      } 
      // otherwise add new entry to list
      else {
        vector<string> v;
        v.push_back(it->m_seq_id);
        m_seq_ids_by_coverage_group.push_back(v);
      }
      
      // Add an entry to lookup the coverage group id frome the seq_id
      m_seq_id_to_coverage_group_map[it->m_seq_id] = m_seq_ids_by_coverage_group.size()-1;
    }
  }

  
  // Set up defaults and build paths without command-line arguments
  Settings::Settings(const string& _base_output_path)
  {
    // things work fine if this is empty ""
    this->base_output_path = _base_output_path;
    
    this->pre_option_initialize();
    
    // no command line arguments here
    
		this->post_option_initialize();
	}
  
  Settings::Settings(int argc, char* argv[])
  {
    this->pre_option_initialize(argc, argv);
    
    // setup and parse configuration options:
    AnyOption options("Usage: breseq -r reference.gbk [-r reference2.gbk ...] reads1.fastq [reads2.fastq.gz ...]");
    options.addUsage("");
    options.addUsage("Run the breseq pipeline for predicting mutations from haploid microbial re-sequencing data.");
    options.addUsage("");
    options.addUsage("FASTQ read files (which may be gzipped) are input as the last unnamed argument(s).");
    options.addUsage("");
    options.addUsage("Allowed Options");
    
    options
		("help,h", "Produce help message showing advanced options", TAKES_NO_ARGUMENT)
    ("reference,r", "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)")
    ("name,n", "Human-readable name of the analysis run for output (DEFAULT=<none>)", "")
    ("num-processors,j", "Number of processors to use in multithreaded steps", 1)
    //("verbose,v","Produce verbose output",TAKES_NO_ARGUMENT, ADVANCED_OPTION) @JEB - not consistently implemented
		("output,o", "Path to breseq output", ".")
    ("polymorphism-prediction,p", "The sample is not clonal. Predict polymorphic (mixed) mutations. Setting this flag changes from CONSENSUS MODE (the default) to POLYMORPHISM MODE", TAKES_NO_ARGUMENT);

    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Read File Options", ADVANCED_OPTION);
    options
    ("limit-fold-coverage,l", "Analyze a subset of the input FASTQ sequencing reads with enough bases to provide this theoretical coverage of the reference sequences. A value between 60 and 120 will usually speed up the analysis with no loss in sensitivity for clonal samples. The actual coverage achieved will be somewhat less because not all reads will map (DEFAULT=OFF)", "", ADVANCED_OPTION)
    ("nanopore-fast-basecalling,x", "Set recommended options for nanopore data processed in fast basecalling mode to rule out false-positive mutations due to high homopolymer indel error rates. Equivalent to --consensus-reject-indel-homopolymer-length 4 --polymorphism-reject-indel-homopolymer-length 4 consensus/polymorphism --polymorphism-no-indel --bowtie2-stage1 \"" + this->bowtie2_stage2 + "\" --bowtie2-stage2 \"\". If you provide any of these options on their own, then they will override these preset options. NOTE: For data processed in high-accuracy basecalling mode, this option is not necessary.", TAKES_NO_ARGUMENT)
    ("aligned-sam", "Input files are aligned SAM files, rather than FASTQ files. Junction prediction steps will be skipped. Be aware that breseq assumes: (1) Your SAM file is sorted such that all alignments for a given read are on consecutive lines. You can use 'samtools sort -n' if you are not sure that this is true for the output of your alignment program. (2) You EITHER have alignment scores as additional SAM fields with the form 'AS:i:n', where n is a positive integer and higher values indicate a better alignment OR it defaults to calculating an alignment score that is equal to the number of bases in the read minus the number of inserted bases, deleted bases, and soft clipped bases in the alignment to the reference. The default highly penalizes split-read matches (with CIGAR strings such as M35D303M65).", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("read-min-length", "Reads in the input FASTQ file that are shorter than this length will be ignored. (0 = OFF)", 18, ADVANCED_OPTION)
    ("read-max-same-base-fraction", "Reads in the input FASTQ file in which this fraction or more of the bases are the same will be ignored. (0 = OFF)", 0.9, ADVANCED_OPTION)
    ("read-max-N-fraction", "Reads in the input FASTQ file in which this fraction or more of the bases are uncalled as N will be ignored. (0 = OFF)", 0.5, ADVANCED_OPTION)
    ("long-read-trigger-length", "Mark a file as containing long reads and enable read splitting if the longest read has a length that is greater than or equal to this value. (0 = OFF)", 1000, ADVANCED_OPTION)
    ("long-read-split-length", "Split input reads in a file marked as having long reads into pieces that are at most this many bases long. Using values much larger than the default for this parameter will likely degrade the speed and accuracy of breseq because of how it performs mapping and analyzes split-read alignments. Filters such as --read-min-length are applied to split reads. (0 = OFF)", 200, ADVANCED_OPTION)
    ("long-read-distribute-remainder", "When splitting long reads, divide them into equal pieces that are less than the split length. If this option is not chosen (the default), reads will be split into chunks with exactly the split length and any remaining bases after the last chunk will be ignored.", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("genbank-field-for-seq-id", "Which GenBank header field will be used to assign sequence IDs. Valid choices are LOCUS, ACCESSION, and VERSION. The default is to check those fields, in that order, for the first one that exists. If you override the default, you will need to use the converted reference file (data/reference.gff) for further breseq and gdtools operations on breseq output!", "AUTOMATIC", ADVANCED_OPTION)
    ;
    
    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Reference File Options", ADVANCED_OPTION);
    options
    ("contig-reference,c", "File containing reference sequences in GenBank, GFF3, or FASTA format. The same coverage distribution will be fit to all of the reference sequences in this file simultaneously. This is appropriate when they are all contigs from a genome that should be present with the same copy number. Use of this option will improve performance when there are many contigs and especially when some are very short (≤1,000 bases).", NULL, ADVANCED_OPTION)
    ("junction-only-reference,s", "File containing reference sequences in GenBank, GFF3, or FASTA format. These references are only used for calling junctions with other reference sequences. An example of appropriate usage is including a transposon sequence not present in a reference genome. Option may be provided multiple times for multiple files.", NULL, ADVANCED_OPTION)
    ("targeted-sequencing,t", "Reference sequences were targeted for ultra-deep sequencing (using pull-downs or amplicons). Do not fit coverage distribution.", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("user-evidence-gd","User supplied Genome Diff file of JC and/or RA evidence items. The breseq output will report the support for these sequence changes even if they do not pass the normal filters for calling mutations in this sample.", "", ADVANCED_OPTION)
    ;
    
    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Read Alignment Options", ADVANCED_OPTION);
    options
    ("minimum-mapping-quality,m", "Ignore alignments with less than this mapping quality (MQ) when calling mutations. MQ scores are equal to -10log10(P), where P is the probability that the best alignment is not to the correct location in the reference genome. The range of MQ scores returned by bowtie2 is 0 to 255.", 0, ADVANCED_OPTION)
    ("base-quality-cutoff,b", "Ignore bases with quality scores lower than this value", 3, ADVANCED_OPTION)
    ("quality-score-trim", "Trim the ends of reads past any base with a quality score below --base-quality-score-cutoff.", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("require-match-length", "Only consider alignments that cover this many bases of a read", 0, ADVANCED_OPTION)
    ("require-match-fraction", "Only consider alignments that cover this fraction of a read", 0.9, ADVANCED_OPTION)
    ("maximum-read-mismatches", "Don't consider reads with this many or more bases or indels that are different from the reference sequence. Unaligned bases at the end of a read also count as mismatches. Unaligned bases at the beginning of the read do NOT count as mismatches. (DEFAULT=OFF)", "", ADVANCED_OPTION)
    ;
    
    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Bowtie2 Mapping/Alignment Options", ADVANCED_OPTION);
    options
    ("bowtie2-scoring", "All calls to bowtie2 must use the same commands for assigning scores to read alignments. Larger scores are assumed to be better by breseq. Each call to bowtie2 has this option added to its command line. (DEFAULT=\"" + this->bowtie2_scoring + "\")", "", ADVANCED_OPTION)
    ("bowtie2-stage1", "Settings for mapping criteria used for the stage 1 alignment. This step is normally meant for quickly aligning near-perfect matches. (DEFAULT=\"" + this->bowtie2_stage1 + "\")", "", ADVANCED_OPTION)
    ("bowtie2-stage2", "Settings for mapping criteria used for the stage 2 alignment. If set to the empty string \"\", then stage 2 alignment is skipped. This step is normally meant for exhaustively mapping reads that were unmapped by stage 1. (DEFAULT=\"" + this->bowtie2_stage2 + "\")", "", ADVANCED_OPTION)
    ("bowtie2-junction", "Settings for mapping criteria used in aligning reads to candidate junctions. (DEFAULT=\"" + this->bowtie2_junction + "\")", "", ADVANCED_OPTION)
    ;
    options.addUsage("In addition to these values, breseq automatically sets the seed size for bowtie2 read mapping (-L option) to a value that is scaled to the read length (r). This value is 0.5 * r for stage 1, 5 + 0.1 * r for stage 2, and 0.3 * r for junction mapping. In each case, it is bounded to the range [4,31] as required by bowtie2. Be warned that breseq internally rescores alignments with a scoring scheme setting +1 for match, -3 for mismatch, -2 for gap open, and -3 for gap extend for consistency when comparing alternative alignments present in the bowtie2 output.", ADVANCED_OPTION);
    
    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Junction (JC) Evidence Options", ADVANCED_OPTION);
    options
    ("no-junction-prediction", "Do not predict new sequence junctions", TAKES_NO_ARGUMENT)
    ("junction-indel-split-length", "Split read alignments on indels of this many or more bases. Indel mutations of this length or longer will be predicted by JC evidence and those that are shorter will be predicted from RA evience", 3, ADVANCED_OPTION)
    ("junction-alignment-pair-limit", "Only consider this many passed alignment pairs when creating candidate junction sequences (0 = DO NOT LIMIT)", 100000, ADVANCED_OPTION)
    ("junction-minimum-candidates", "Test at least this many of the top-scoring junction candidates, regardless of their length", 100, ADVANCED_OPTION)
    ("junction-maximum-candidates", "Test no more than this many of the top-scoring junction candidates (0 = DO NOT LIMIT)", 5000, ADVANCED_OPTION)
    ("junction-candidate-length-factor", "Accept top-scoring junction candidates to test until their cumulative length is this factor times the total reference sequence length (0 = DO NOT LIMIT)", 0.1, ADVANCED_OPTION)
    ("junction-minimum-candidate-pos-hash-score", "Minimum number of distinct spanning read start positions required to create a junction candidate for further testing", 2, ADVANCED_OPTION)
    ("junction-score-cutoff", "Maximum negative log10 probability of uneven coverage across a junction breakpoint to accept (0 = OFF)", 3.0, ADVANCED_OPTION)
    ("junction-minimum-pos-hash-score", "Minimum number of distinct spanning read start positions required to accept a junction (DEFAULT = consensus mode, 3; polymorphism mode, 3)", "", ADVANCED_OPTION)
    ("junction-minimum-side-match", "Minimum number of bases a read must extend past any overlap or read-only sequence at the breakpoint of a junction on each side to count as support for the junction (DEFAULT = consensus mode, 1; polymorphism mode, 6)", "", ADVANCED_OPTION)
    ("junction-minimum-pr-no-read-start-per-position", "Minimum probablilty assigned that no mapped read will start at a given position and strand for junction prediction", 0.1, ADVANCED_OPTION)
    ("junction-allow-suboptimal-matches", "Assign a read to the junction candidate with the most overall support as long as its match to this junction is better than to any location in the reference sequence, even if it matches a different junction candidate better. This behavior was the default before v0.35.0. It will align more reads to junctions but risks misassigning some reads to the wrong junction candidates. It is only recommended that you use this option in CONSENSUS mode", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ;
    
    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Missing Coverage (MC) Evidence Options", ADVANCED_OPTION);
    options
    ("deletion-coverage-seed-cutoff","Value for coverage below which MC are seeded", 0, ADVANCED_OPTION)
    ("deletion-coverage-propagation-cutoff","Value for coverage above which MC ends stop. 0 = calculated from coverage distribution", 0, ADVANCED_OPTION)
    ("call-mutations-overlapping-MC", "If provided, don't ignore mutations predicted from RA evidence that overlap MC evidence", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ;
    
    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Consensus Read Alignment (RA) Evidence Options", ADVANCED_OPTION);
    options
    ("consensus-score-cutoff", "Log10 E-value cutoff for consensus base substitutions and small indels (DEFAULT = 10)", "", ADVANCED_OPTION)
    ("consensus-frequency-cutoff", "Only predict consensus mutations when the variant allele frequency is above this value. (DEFAULT = consensus mode, 0.8; polymorphism mode, 0.8)", "", ADVANCED_OPTION)
    ("consensus-minimum-variant-coverage", "Only predict consensus mutations when at least this many reads support the mutation. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", ADVANCED_OPTION)
    ("consensus-minimum-total-coverage", "Only predict consensus mutations when at least this many reads total are aligned to a genome position. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", ADVANCED_OPTION)
    ("consensus-minimum-variant-coverage-each-strand", "Only predict consensus mutations when at least this many reads on each strand support the mutation. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", ADVANCED_OPTION)
    ("consensus-minimum-total-coverage-each-strand", "Only predict consensus mutations when at least this many reads on each strand are aligned to a genome position. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", ADVANCED_OPTION)
    ("consensus-reject-indel-homopolymer-length", "Reject insertion/deletion polymorphisms which could result from expansion/contraction of homopolymer repeats with this length or greater in the reference genome (0 = OFF) (DEFAULT = consensus mode, OFF; polymorphism mode, OFF) ", "", ADVANCED_OPTION)
    ("consensus-reject-surrounding-homopolymer-length", "Reject polymorphic base substitutions that create a homopolymer with this many or more of one base in a row. The homopolymer must begin and end after the changed base. For example, TATTT->TTTTT would be rejected with a setting of 5, but ATTTT->TTTTT would not. (0 = OFF) (DEFAULT = consensus mode, OFF; polymorphism mode, OFF)", "", ADVANCED_OPTION)
    ;

    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Polymorphism Read Alignment (RA) Evidence Options", ADVANCED_OPTION);
    options

    ("polymorphism-score-cutoff", "Log10 E-value cutoff for test of polymorphism vs no polymorphism (DEFAULT = consensus mode, 10; polymorphism mode, 2)", "", ADVANCED_OPTION)
    ("polymorphism-frequency-cutoff", "Only predict polymorphisms when the minor variant allele frequency is greater than this value. For example, a setting of 0.05 will reject all polymorphisms with a non-reference frequency of <0.05, and any variants with a non-reference frequency of ≥ 0.95 (which is 1 - 0.05) will be rejected as polymorphisms and instead predicted to be consensus mutations (DEFAULT = consensus mode, 0.2; polymorphism mode, 0.05)", "", ADVANCED_OPTION)
    ("polymorphism-minimum-variant-coverage", "Only predict polymorphisms when at least this many reads support each alternative allele. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", ADVANCED_OPTION)
    ("polymorphism-minimum-total-coverage", "Only predict polymorphisms when at least this many reads total are aligned to a genome position. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", ADVANCED_OPTION)
    ("polymorphism-minimum-variant-coverage-each-strand", "Only predict polymorphisms when at least this many reads on each strand support each alternative allele. (DEFAULT = consensus mode, 0; polymorphism mode, 2)", "", ADVANCED_OPTION)
    ("polymorphism-minimum-total-coverage-each-strand", "Only predict polymorphisms when at least this many reads on each strand are aligned to a genome position. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", ADVANCED_OPTION)
    ("polymorphism-bias-cutoff", "P-value criterion for Fisher's exact test for strand bias AND K-S test for quality score bias. (0 = OFF) (DEFAULT = consensus mode, OFF; polymorphism mode, OFF)", "", ADVANCED_OPTION)
    ("polymorphism-no-indels", "Do not predict insertion/deletion polymorphisms ≤" + to_string(kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation) + " bp from read alignment or new junction evidence", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("polymorphism-reject-indel-homopolymer-length", "Reject insertion/deletion polymorphisms which could result from expansion/contraction of homopolymer repeats with this length or greater in the reference genome (0 = OFF) (DEFAULT = consensus mode, OFF; polymorphism mode, 3) ", "", ADVANCED_OPTION)
    ("polymorphism-reject-surrounding-homopolymer-length", "Reject polymorphic base substitutions that create a homopolymer with this many or more of one base in a row. The homopolymer must begin and end after the changed base. For example, TATTT->TTTTT would be rejected with a setting of 5, but ATTTT->TTTTT would not. (0 = OFF) (DEFAULT = consensus mode, OFF; polymorphism mode, 5)", "", ADVANCED_OPTION)
    ;
    
    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Output Options", ADVANCED_OPTION);
    options
    ("max-flanking-columns", "Maximum number of columns in aligned reads to show flanking the region of interest in the HTML output for an evidence item (0=ALL)", 100, ADVANCED_OPTION)
    ("max-displayed-reads", "Maximum number of reads to display in the HTML output for an evidence item (0=ALL)", 100, ADVANCED_OPTION)
    ("brief-html-output", "Don't create detailed output files for evidence (no read alignments or coverage plots)", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("header-genome-diff,g", "Include header information from this GenomeDiff file in output.gd", "", ADVANCED_OPTION)
    ("no-javascript", "Don't include javascript in the HTML output", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ;
    
    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Pipeline Control Options", ADVANCED_OPTION);
    options
    ("skip-RA-MC-prediction", "Skip generating read alignment and missing coverage evidence.", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("skip-JC-prediction", "Skip generating new junction evidence.", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("skip-MC-prediction", "Skip generating missing coverage evidence.", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ;
    
    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Debugging Options", ADVANCED_OPTION);
    options
    ("keep-intermediates,k","Do not delete intermediate files.", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("per-position-file", "Create additional file of per-position aligned bases", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("junction-debug", "Output additional junction debugging files", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ;
    
    options.addUsage("", ADVANCED_OPTION);
    options.addUsage("Experimental Options (Use at your own risk)", ADVANCED_OPTION);
    options
    ("cnv","Do experimental copy number variation prediction",TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ("cnv-tile-size", "Tile size for copy number variation prediction", 500, ADVANCED_OPTION)
    ("cnv-ignore-redundant", "Only consider non-redundant coverage when using cnv", TAKES_NO_ARGUMENT, ADVANCED_OPTION)
    ;
    
    options.processCommandArgs(argc, argv);
    
    options.addUsage("");
		options.addUsage("Utility Command Usage: breseq [command] options ...");
    options.addUsage("  Sequence Utility Commands: CONVERT-FASTQ, CONVERT-REFERENCE, GET-SEQUENCE");
    options.addUsage("  Breseq Post-Run Commands: BAM2ALN, BAM2COV, CL-TABULATE");
    options.addUsage("  Other Commands: SOFT-CLIPPING");
    options.addUsage("");
    options.addUsage("For help using a utility command, type: breseq [command] "); 
    options.addUsage("");
    options.addUsage(output_divider);
    
    // make sure that the other config options are good:
    if (options.count("help"))
    {
      options.printAdvancedUsage();
      exit(-1);
    }
    
    //! Settings: Global Workflow and Output
    
    this->base_output_path = options["output"];
    
    // Do this so the default printed can be pretty
    if (this->base_output_path == "current directory")
      this->base_output_path = ".";
    
    vector<string> split_output_path = split(this->base_output_path, "\n");
    if (split_output_path.size() > 1) {
      options.addUsage("");
      options.addUsage("Only one output path can be specified (-o).");
      options.printUsage();
      exit(-1);
    }
    
    // Remaining command line items are read files
    // Read sequence file provided?
		if (options.getArgc() == 0)
		{
      options.addUsage("");
      options.addUsage("No read sequence files provided.");
      options.printUsage();
      exit(-1);		
    }
    for (int32_t i = 0; i < options.getArgc(); i++)
    {
      string read_file_name = options.getArgv(i);
      this->read_file_names.push_back(read_file_name);
    }
    this->aligned_sam_mode = options.count("aligned-sam");
    if (aligned_sam_mode) {
      cerr << "Input files are aligned SAM instead of FASTQ (--aligned-sam option)." << endl;
      cerr << "No junction prediction will take place." << endl;
      cerr << output_divider << endl;
      this->skip_new_junction_prediction = true;
    }
    this->read_files.Init(read_file_names, this->aligned_sam_mode);
    
    if (options.count("limit-fold-coverage")) {
      this->read_file_coverage_fold_limit = from_string<double>(options["limit-fold-coverage"]);
    }
    
    this->read_file_max_same_base_fraction = from_string<double>(options["read-max-same-base-fraction"]);
    this->read_file_max_N_fraction = from_string<double>(options["read-max-N-fraction"]);
    this->read_file_long_read_trigger_length = from_string<uint32_t>(options["long-read-trigger-length"]);
    this->read_file_long_read_split_length = from_string<uint32_t>(options["long-read-split-length"]);
    this->read_file_long_read_distribute_remainder = options.count("long-read-distribute-remainder");

    
    bool long_read_include_remainder;           // Default = false COMMAND-LINE OPTION
    uint32_t long_read_split_length;            // Default = 500 COMMAND-LINE OPTION
    
    this->genbank_field_for_seq_id = options["genbank-field-for-seq-id"];
    this->genbank_field_for_seq_id = to_upper(this->genbank_field_for_seq_id);
    if (   (this->genbank_field_for_seq_id != "AUTOMATIC")
        && (this->genbank_field_for_seq_id != "LOCUS")
        && (this->genbank_field_for_seq_id != "VERSION")
        && (this->genbank_field_for_seq_id != "ACCESSION")
        ) {
      options.addUsage("");
      options.addUsage("Value of --genbank-field-for-seq-id must be one of the following: AUTOMATIC, LOCUS, VERSION, ACCESSION.");
      options.addUsage("");
      options.addUsage("Value provided was: " + this->genbank_field_for_seq_id);
      options.printUsage();
      exit(-1);
    }
    
    
    // Reference sequence provided?
		if (options.count("reference") + options.count("contig-reference") + options.count("junction-only-reference") == 0)
		{
      options.addUsage("");
      options.addUsage("No reference sequences provided (-r|-c|-s).");
      options.printUsage();
      exit(-1);
		}
    
    if (options.count("reference")) {
      this->normal_reference_file_names = from_string<vector<string> >(options["reference"]);
    }
    this->all_reference_file_names.insert(  this->all_reference_file_names.end(), 
                                            this->normal_reference_file_names.begin(), 
                                            this->normal_reference_file_names.end() 
                                         ); 
                                        
    // Important to check for NULL before converting
    if (options.count("contig-reference")) {
      this->contig_reference_file_names = from_string<vector<string> >(options["contig-reference"]);
    }
    
    this->all_reference_file_names.insert(  this->all_reference_file_names.end(), 
                                            this->contig_reference_file_names.begin(), 
                                            this->contig_reference_file_names.end() 
                                         ); 
    
    // Important to check for NULL before converting
    if (options.count("junction-only-reference")) {
      this->junction_only_file_names = from_string<vector<string> >(options["junction-only-reference"]);
    }
    
    this->all_reference_file_names.insert(  this->all_reference_file_names.end(), 
                                            this->junction_only_file_names.begin(), 
                                            this->junction_only_file_names.end() 
                                         );    
    
    // Note that for full setup of the cReferenceSequenceSettings object,
    // we require a cReferenceSequences object to be instantiated outside of the cSettings object
    // and a subsequent call to init_reference_sequences()
    
    this->user_evidence_genome_diff_file_name = options["user-evidence-gd"];
    
    this->custom_run_name = options["name"];
    
    this->num_processors = from_string<int32_t>(options["num-processors"]);
    
    this->do_copy_number_variation = options.count("cnv");
    this->copy_number_variation_tile_size = from_string<uint32_t>(options["cnv-tile-size"]);
    this->ignore_redundant_coverage = options.count("cnv-ignore-redundant");
    
    this->verbose = options.count("verbose");
    
    //! Settings: Read Alignment and Candidate Junction Read Alignment

    this->require_match_length = from_string<uint32_t>(options["require-match-length"]);
    this->require_match_fraction = from_string<double>(options["require-match-fraction"]);
    if (options.count("maximum-read-mismatches")) {
        this->maximum_read_mismatches = from_string<int32_t>(options["maximum-read-mismatches"]);
    }
    
    //! Settings: Mutation Identification
    this->base_quality_cutoff = from_string<uint32_t>(options["base-quality-cutoff"]);
    if (options.count("quality-score-trim")) this->quality_score_trim = this->base_quality_cutoff;

    this->deletion_coverage_propagation_cutoff = from_string<double>(options["deletion-coverage-propagation-cutoff"]);
    ASSERT(this->deletion_coverage_propagation_cutoff >= 0, "Argument --deletion-coverage-propagation-cutoff must be >= 0")
    this->deletion_coverage_seed_cutoff = from_string<double>(options["deletion-coverage-seed-cutoff"]);
    ASSERT(this->deletion_coverage_propagation_cutoff >= 0, "Argument --deletion-coverage-seed-cutoff must be >= 0")
    
    this->call_mutations_overlapping_missing_coverage = options.count("call-mutations-overlapping-MC");
    
    //! Settings: bowtie2
    //  all have default that we only overwrite if present on command line
    if (options.count("bowtie2-scoring")) this->bowtie2_scoring = options["bowtie2-scoring"];
    if (options.count("bowtie2-stage1")) this->bowtie2_stage1 = options["bowtie2-stage1"];
    if (options.count("bowtie2-stage2")) this->bowtie2_stage2 = options["bowtie2-stage2"];
    if (options.count("bowtie2-junction")) this->bowtie2_junction = options["bowtie2-junction"];

    //! Settings: Junction Prediction
    this->skip_new_junction_prediction = this->skip_new_junction_prediction || options.count("no-junction-prediction");
    this->minimum_candidate_junctions = from_string<int32_t>(options["junction-minimum-candidates"]);
    this->maximum_candidate_junctions = from_string<int32_t>(options["junction-maximum-candidates"]);
    this->maximum_candidate_junction_length_factor = from_string<double>(options["junction-candidate-length-factor"]);
    this->minimum_candidate_junction_pos_hash_score = from_string<double>(options["junction-minimum-candidate-pos-hash-score"]);
    this->maximum_junction_sequence_passed_alignment_pairs_to_consider = from_string<uint64_t>(options["junction-alignment-pair-limit"]);
    this->junction_pos_hash_neg_log10_p_value_cutoff = from_string<double>(options["junction-score-cutoff"]);
    this->junction_allow_suboptimal_matches = options.count("junction-allow-suboptimal-matches");
    
    //! Settings: Pipeline Control
    this->skip_read_alignment_and_missing_coverage_prediction = options.count("skip-RA-MC-prediction");
    this->skip_new_junction_prediction = this->skip_new_junction_prediction || options.count("skip-JC-prediction");
    this->skip_missing_coverage_prediction = options.count("skip-MC-prediction");
    
    //! Settings: Debugging
    this->keep_all_intermediates = options.count("keep-intermediates");

    
    //
    // Set the read alignment evidence model
    // Different defaults for the three modes:
    //   Mixed base - consensus and strong polymorphism calls as marginal (new default mode)
    //   Polymorphism prediction - consensus and any mixed calls (--p option)
    //
    
    this->polymorphism_prediction = options.count("polymorphism-prediction");
    if (this->polymorphism_prediction) {
      
      this->minimum_mapping_quality = 0;
      
      this->mutation_log10_e_value_cutoff = 10;
      this->consensus_frequency_cutoff = 0.8; // zero is OFF - ensures any rejected poly with high freq move to consensus!
      this->consensus_minimum_variant_coverage = 0;
      this->consensus_minimum_total_coverage = 0;
      this->consensus_minimum_variant_coverage_each_strand = 0;
      this->consensus_minimum_total_coverage_each_strand = 0;
      this->consensus_reject_indel_homopolymer_length = 0; // zero is OFF!
      this->consensus_reject_surrounding_homopolymer_length = 0; // zero is OFF!
      
      this->polymorphism_log10_e_value_cutoff = 2;
      this->polymorphism_frequency_cutoff = 0.05;
      this->polymorphism_minimum_variant_coverage = 0;
      this->polymorphism_minimum_total_coverage = 0;
      this->polymorphism_minimum_variant_coverage_each_strand = 2;
      this->polymorphism_minimum_total_coverage_each_strand = 0;
      
      this->mixed_base_prediction = false;
      this->polymorphism_reject_indel_homopolymer_length = 3; // zero is OFF!
      this->polymorphism_reject_surrounding_homopolymer_length = 5; // zero is OFF!
      this->polymorphism_bias_p_value_cutoff = 0;
      this->polymorphism_no_indels = false;
      this->polymorphism_precision_decimal = 0.000001;
      this->polymorphism_precision_places = 8;
      
      this->minimum_alignment_resolution_pos_hash_score = 3;
      this->junction_minimum_side_match = 6;
      this->junction_pos_hash_neg_log10_p_value_cutoff = 0; // OFF
      
    }
    // This is strictly true if we are not in polymorphism mode...
    else /* if (this->mixed_base_prediction)*/ {

      this->minimum_mapping_quality = 0;
      
      this->mutation_log10_e_value_cutoff = 10;
      this->consensus_frequency_cutoff = 0.8;
      this->consensus_minimum_variant_coverage = 0;
      this->consensus_minimum_total_coverage = 0;
      this->consensus_minimum_variant_coverage_each_strand = 0;
      this->consensus_minimum_total_coverage_each_strand = 0;
      this->consensus_reject_indel_homopolymer_length = 0; // zero is OFF!
      this->consensus_reject_surrounding_homopolymer_length = 0; // zero is OFF!
      
      this->polymorphism_log10_e_value_cutoff = 10;
      this->polymorphism_frequency_cutoff = 0.2;
      this->polymorphism_minimum_variant_coverage = 0;
      this->polymorphism_minimum_total_coverage = 0;
      this->polymorphism_minimum_variant_coverage_each_strand = 0;
      this->polymorphism_minimum_total_coverage_each_strand = 0;

      this->polymorphism_reject_indel_homopolymer_length = 0; // zero is OFF!
      this->polymorphism_reject_surrounding_homopolymer_length = 0; // zero is OFF!
      this->polymorphism_bias_p_value_cutoff = 0;
      this->polymorphism_no_indels = false;
      this->polymorphism_precision_decimal = 0.000001;
      this->polymorphism_precision_places = 3;
      
      this->minimum_alignment_resolution_pos_hash_score = 3;
      this->junction_minimum_side_match = 1;
    }
    
    // en masse overrides come before specific overrides
    
    if (options.count("nanopore-fast-basecalling")) {
      this->consensus_reject_indel_homopolymer_length = 4;
      this->polymorphism_reject_indel_homopolymer_length = 4;
      this->polymorphism_no_indels = true;
      
      // Just do what is normally stage 2 alignment
      this->bowtie2_stage1 = this->bowtie2_stage2;
      this->bowtie2_stage2 = "";
    }
      
    // override the default settings
    
    if (options.count("minimum-mapping-quality")) {
      this->minimum_mapping_quality = from_string<int32_t>(options["minimum-mapping-quality"]);
    }
    
    if (options.count("junction-indel-split-length"))
      this->preprocess_junction_min_indel_split_length = from_string<int32_t>(options["junction-indel-split-length"]);
    
    if (options.count("consensus-score-cutoff"))
      this->mutation_log10_e_value_cutoff = from_string<double>(options["consensus-score-cutoff"]);
    if (this->mutation_log10_e_value_cutoff < 0) {
      options.addUsage("");
      options.addUsage("--consensus-score-cutoff must be ≥0");
      options.printUsage();
      exit(-1);
    }
    if (options.count("consensus-frequency-cutoff"))
      this->consensus_frequency_cutoff = from_string<double>(options["consensus-frequency-cutoff"]);
    
    if (options.count("consensus-minimum-variant-coverage"))
      this->consensus_minimum_variant_coverage = from_string<uint32_t>(options["consensus-minimum-variant-coverage"]);
    if (options.count("consensus-minimum-total-coverage"))
      this->consensus_minimum_total_coverage = from_string<uint32_t>(options["consensus-minimum-total-coverage"]);
    if (options.count("consensus-minimum-variant-coverage-each-strand"))
      this->consensus_minimum_variant_coverage_each_strand = from_string<uint32_t>(options["consensus-minimum-variant-coverage-each-strand"]);
    if (options.count("consensus-minimum-total-coverage-each-strand"))
      this->consensus_minimum_total_coverage_each_strand = from_string<uint32_t>(options["consensus-minimum-total-coverage-each-strand"]);
    if (options.count("consensus-reject-indel-homopolymer-length"))
      this->consensus_reject_indel_homopolymer_length = from_string<int32_t>(options["consensus-reject-indel-homopolymer-length"]);
    if (options.count("consensus-reject-surrounding-homopolymer-length"))
      this->consensus_reject_surrounding_homopolymer_length = from_string<int32_t>(options["consensus-reject-surrounding-homopolymer-length"]);
    
    if (options.count("polymorphism-frequency-cutoff"))
      this->polymorphism_frequency_cutoff = from_string<double>(options["polymorphism-frequency-cutoff"]);
    
    if (options.count("polymorphism-minimum-variant-coverage"))
      this->polymorphism_minimum_variant_coverage = from_string<uint32_t>(options["polymorphism-minimum-variant-coverage"]);
    if (options.count("polymorphism-minimum-total-coverage"))
      this->polymorphism_minimum_total_coverage = from_string<uint32_t>(options["polymorphism-minimum-total-coverage"]);
    if (options.count("polymorphism-minimum-variant-coverage-each-strand"))
      this->polymorphism_minimum_variant_coverage_each_strand = from_string<uint32_t>(options["polymorphism-minimum-variant-coverage-each-strand"]);
    if (options.count("polymorphism-minimum-total-coverage-each-strand"))
      this->polymorphism_minimum_total_coverage_each_strand = from_string<uint32_t>(options["polymorphism-minimum-total-coverage-each-strand"]);
    
    if (options.count("polymorphism-no-indels"))
      this->polymorphism_no_indels = options.count("polymorphism-no-indels");
    if (options.count("polymorphism-reject-indel-homopolymer-length"))
      this->polymorphism_reject_indel_homopolymer_length = from_string<int32_t>(options["polymorphism-reject-indel-homopolymer-length"]);
    if (options.count("polymorphism-reject-surrounding-homopolymer-length"))
      this->polymorphism_reject_surrounding_homopolymer_length = from_string<int32_t>(options["polymorphism-reject-surrounding-homopolymer-length"]);
    if (options.count("polymorphism-score-cutoff"))
      this->polymorphism_log10_e_value_cutoff = from_string<double>(options["polymorphism-score-cutoff"]);
    if (this->polymorphism_log10_e_value_cutoff < 0) {
      options.addUsage("");
      options.addUsage("--polymorphism-score-cutoff must be ≥0");
      options.printUsage();
      exit(-1);
    }
    if (options.count("polymorphism-bias-cutoff"))
      this->polymorphism_bias_p_value_cutoff = from_string<double>(options["polymorphism-bias-cutoff"]);

    // Warn of possibly confusing settings
    if (this->consensus_frequency_cutoff > 1 - this->polymorphism_frequency_cutoff) {
      cerr << endl << "=== WARNING ===" << endl;
      cerr << "Mutations with a frequency between " << (1 - this->polymorphism_frequency_cutoff) << " and " << this->consensus_frequency_cutoff << " will not be predicted with these settings. ";
      cerr << "Unless this is your intent, either set \'--consensus_frequency_cutoff " << (1 - this->polymorphism_frequency_cutoff) << "\' or set \'--polymorphism_frequency_cutoff " << (1 - this->consensus_frequency_cutoff) << "\' on the command line." << endl << endl;
    }
    
    // Junction options
    if (options.count("junction-minimum-pos-hash-score"))
      this->minimum_alignment_resolution_pos_hash_score = from_string<uint32_t>(options["junction-minimum-pos-hash-score"]);
    
    if (options.count("junction-minimum-side-match"))
      this->junction_minimum_side_match = from_string<int32_t>(options["junction-minimum-side-match"]);

    // Not a command line option... redundant with new setting but also separate.
    this->required_both_unique_length_per_side = this->junction_minimum_side_match;
    
    
    this->targeted_sequencing = options.count("targeted-sequencing");
    if (this->targeted_sequencing)
      this->skip_missing_coverage_prediction = true;
    
    this->print_mutation_identification_per_position_file = options.count("per-position-file");
    
    this->junction_debug = options.count("junction-debug");
    
    this->max_displayed_reads = from_string<int32_t>(options["max-displayed-reads"]);
    this->max_flanking_columns = from_string<int32_t>(options["max-flanking-columns"]);
    this->skip_alignment_or_plot_generation = options.count("brief-html-output");
    this->no_javascript = options.count("no-javascript");
    if (options.count("header-genome-diff"))
      this->header_genome_diff_file_name = options["header-genome-diff"];
  
		this->post_option_initialize();
    
    //////// Check here for any conflicting options
    
    /*
    if (this->user_evidence_genome_diff_file_name.size() && !this->polymorphism_prediction) {
      ERROR("You must run breseq in polymorphism mode (-p) when supplying --user-evidence-gd.");
    }
    */
    
    ////////
    
    
    
    // Log the command line
    time_t stamp_time = time(NULL);
    this->log(ctime(&stamp_time));	
    this->log(this->full_command_line + "\n");
	}
  
  void Settings::command_line_run_header()
  {
    cerr << output_divider << endl;
    fprintf(stderr, "%s  %s   %s\n", PACKAGE_STRING, GITHUB_REVISION_STRING, PACKAGE_URL);
    fprintf(stderr, "\n");
    fprintf(stderr, "Active Developers: Barrick JE, Deatherage DE\n");
    fprintf(stderr, "Contact:           <%s>\n", PACKAGE_BUGREPORT);
    fprintf(stderr, "\n");
    fprintf(stderr, "License GPLv2+: GNU GPL version 2 or later <https://gnu.org/licenses/gpl.html>");
    fprintf(stderr, "This is free software: you are free to change and redistribute it. \n");
    fprintf(stderr, "There is NO WARRANTY, to the extent permitted by law.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Copyright (c) 2008-2010 Michigan State University\n");
    fprintf(stderr, "Copyright (c) 2011-2022 The University of Texas at Austin\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "If you use breseq in your research, please cite:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Deatherage, D.E., Barrick, J.E. (2014) Identification of mutations\n"); 
    fprintf(stderr, "  in laboratory-evolved microbes from next-generation sequencing\n");
    fprintf(stderr, "  data using breseq. Methods Mol. Biol. 1151: 165–188.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "If you use structural variation (junction) predictions, please cite:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Barrick, J.E., Colburn, G., Deatherage D.E., Traverse, C.C.,\n");
    fprintf(stderr, "  Strand, M.D., Borges, J.J., Knoester, D.B., Reba, A., Meyer, A.G. \n");
    fprintf(stderr, "  (2014) Identifying structural variation in haploid microbial genomes \n");
    fprintf(stderr, "  from short-read resequencing data using breseq. BMC Genomics 15:1039.\n");
    cerr << output_divider << endl;
  }

  // Settings in here should be static (not set at the command line).
  // Options that are set are in Settings::Settings(int argc, char* argv[])
	void Settings::pre_option_initialize(int argc, char* argv[])
	{
    
    ////////////////////
    //! Data
    ////////////////////
    
    this->byline = "<b><i>breseq</i></b>&nbsp;&nbsp;version ";
    this->byline += PACKAGE_VERSION;
    this->byline += "&nbsp;&nbsp;";
    this->byline += GITHUB_REVISION_STRING;

    this->website = PACKAGE_URL;
    
    // Settings
    // Initialize all variables that do not have default initializers (non-strings mostly) 

    this->arguments = "";
    if (argc > 0) this->full_command_line = argv[0];
    for(int i=1; i<argc; i++)
    {
      if (!this->arguments.empty()) this->arguments += " ";
      this->arguments += argv[i];
    }
    this->full_command_line += " " + this->arguments;
   
    ////////////////////
    //! Settings
    ////////////////////
    
    //! Settings: Global Workflow and Output
    
    //! Read file options
    this->aligned_sam_mode  = false;
    this->read_file_coverage_fold_limit = 0.0;
    this->read_file_max_same_base_fraction = 0.9;
    this->read_file_read_length_min = 18;
    this->read_file_long_read_trigger_length = 1000;
    this->read_file_long_read_split_length = 200;
    this->read_file_long_read_distribute_remainder = false;
    
    //! Options that control which parts of the pipeline to execute
    this->skip_read_filtering = false;
    this->skip_new_junction_prediction = false;
		this->skip_read_alignment_and_missing_coverage_prediction = false;
		this->skip_missing_coverage_prediction = false;
    this->skip_alignment_or_plot_generation = false;
		this->do_copy_number_variation = false;
		this->do_periodicity = false;
    
    //! DEBUG options
    this->verbose = false;
		this->alignment_read_limit = 0;         
		this->resolve_alignment_read_limit = 0; 
    this->candidate_junction_read_limit = 0;
    this->no_unmatched_reads = false;
    this->keep_all_intermediates = false;

    //! Settings: Read Alignment and Candidate Junction Read Alignment
    this->require_match_length = 0;         
    this->require_match_fraction = 0.9; // CL value overrides
    this->maximum_read_mismatches = -1;

    this->bowtie2_junction_maximum_alignments_to_consider_per_read = 2000;
    this->bowtie2_genome_maximum_alignments_to_consider_per_read = 2000;

    this->bowtie2_scoring = "--ma 1 --mp 3 --np 0 --rdg 2,3 --rfg 2,3 --ignore-quals";
    //this->bowtie2_score_parameters = "--ma 1 --mp 3 --np 0 --rdg 2,3 --rfg 2,3";
    
    // Different values for these:
    // We report ALL alignments for junctions
    string bowtie2_junction_alignment_reporting_parameters = (bowtie2_junction_maximum_alignments_to_consider_per_read > 0)
    ? " -k " + to_string(this->bowtie2_junction_maximum_alignments_to_consider_per_read) : " -a";
    
    // Need to report all matches for proper handling of
    string bowtie2_genome_alignment_reporting_parameters = (bowtie2_genome_maximum_alignments_to_consider_per_read > 0)
    ? " -k " + to_string(this->bowtie2_genome_maximum_alignments_to_consider_per_read) : " -a";
        
    // Note: this leaves off -L, since it is set based on read length
    
    // Set these defaults
    this->bowtie2_stage1 = "--local -i S,1,0.25 --score-min L,1,0.9" + bowtie2_genome_alignment_reporting_parameters; // "-L 22 -i S,1,0.25 --score-min L,0,0.9";
    
    this->bowtie2_stage2 = "--local -i S,1,0.25 --score-min L,6,0.2" + bowtie2_genome_alignment_reporting_parameters; // "-L 9 -i C,1,0 --score-min L,6,0.2";
    
    //JEB: Changed 2015-07-13 to not give too many alignments for longer (150 bp) reads, but to keep thresholds similar for short (36 bp) reads
    this->bowtie2_junction = "--local -i S,1,0.25 --score-min L,1,0.70" + bowtie2_junction_alignment_reporting_parameters; // "-L 9 -i C,1,0 --score-min L,6,0.2";

    
    this->num_processors = 1;
    
    //! Settings: Candidate Junction Prediction
		this->preprocess_junction_min_indel_split_length = 3;
    this->required_both_unique_length_per_side_fraction = 0.2; 
    this->unmatched_end_length_factor =  1 - this->require_match_fraction;
    this->unmatched_end_minimum_read_length = 50;
    this->maximum_junction_sequence_negative_overlap_length_fraction = 0.4;
    this->maximum_junction_sequence_negative_overlap_length_minimum = 12;
    this->maximum_junction_sequence_positive_overlap_length_fraction = 0.4;
    this->maximum_junction_sequence_positive_overlap_length_minimum = 12;
    this->highly_redundant_junction_ignore_passed_pair_limit = 200;

    // Extra options that are mostly phased out because they are hard nucleotide cutoffs
    this->maximum_junction_sequence_insertion_length = 0;
    this->maximum_junction_sequence_overlap_length = 0;
    this->required_both_unique_length_per_side = 0;
    this->required_one_unique_length_per_side = 0;
    
    this->minimum_candidate_junction_pos_hash_score = 2;
    this->penalize_negative_junction_overlap = true;

    //! Settings: Alignment Resolution
    this->add_split_junction_sides = true;
    this->junction_pos_hash_neg_log10_p_value_cutoff = 3;
    this->minimum_alignment_resolution_pos_hash_score = 3;
    this->minimum_pr_no_read_start_per_position = 0.1;
    this->junction_allow_suboptimal_matches = false;

    //! Settings: Mutation Identification
    this->base_quality_cutoff = 3;
    this->quality_score_trim = 0;
    this->deletion_coverage_propagation_cutoff = 0;
    this->deletion_coverage_seed_cutoff = 0;
    this->call_mutations_overlapping_missing_coverage = false;
      
    this->polymorphism_prediction = false;
    this->mixed_base_prediction = true;
    
    this->consensus_minimum_variant_coverage = 0;
    this->consensus_minimum_total_coverage = 0;
    this->consensus_minimum_variant_coverage_each_strand = 0;
    this->consensus_minimum_total_coverage_each_strand = 0;
    
    this->polymorphism_log10_e_value_cutoff = this->mutation_log10_e_value_cutoff;
		this->polymorphism_bias_p_value_cutoff = 0;
		this->polymorphism_frequency_cutoff = 0.1;
		this->polymorphism_minimum_variant_coverage = 0;
    this->polymorphism_minimum_total_coverage = 0;
    this->polymorphism_minimum_variant_coverage_each_strand = 0;
    this->polymorphism_minimum_total_coverage_each_strand = 0;
    this->polymorphism_reject_indel_homopolymer_length = 0;
    this->polymorphism_reject_surrounding_homopolymer_length = 0;
		this->polymorphism_no_indels = false;
    
    //! Settings: Mutation Prediction
    this->size_cutoff_AMP_becomes_INS_DEL_mutation = kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation;
    this->ignore_within_this_multiple_of_average_read_length_of_contig_end = kBreseq_ignore_within_this_multiple_of_average_read_length_of_contig_end;
    
    //! Settings: Output
    this->max_flanking_columns = 100;
    this->max_displayed_reads = 100;
    this->alignment_mask_ref_matches = false;
    this->max_nucleotides_to_show_in_tables = 20;
		this->max_rejected_read_alignment_evidence_to_show = 20;
		this->max_rejected_junction_evidence_to_show = 10;
		this->hide_circular_genome_junctions = true;
    
		this->lenski_format = false;
		this->no_evidence = false;
    this->no_javascript = false;
    this->no_list_js = false;
    this->add_metadata_to_gd = true;
    
    ////////////////////
    //! File Paths
    ////////////////////
    
    this->bin_path = get_bin_path();
    this->program_data_path = get_program_data_path();
    //ASSERT(this->bin_path.size() > 0, "Uninitialized bin_path");
    ASSERT(this->program_data_path.size() > 0, "Uninitialized program_data_path");
  }

	void Settings::post_option_initialize()
	{     
    this->init_installed();
    
    ////////////////////
    //! Settings
    ////////////////////
    
		// problems if there are spaces b/c shell removes quotes before we know about them
		// thus require run names to only use underscores (but when printing output, remove).
    this->print_custom_run_name = substitute(this->custom_run_name, "_", " ");
    
    ////////////////////
    //! File Paths
    ////////////////////
    
		//// '#' replaced with read file name
		//// '@' replaced by seq_id of reference sequence

    // Do not allow the output path to have spaces. This will cause many invocations we can't trust?
    // ASSERT(this->base_output_path.find(" ") == string::npos, "Output path cannot contain spaces: \"" + this->base_output_path + "\".");
    
    //! Paths: Sequence conversion
		this->sequence_conversion_path = "01_sequence_conversion";
		if (this->base_output_path.size() > 0) this->sequence_conversion_path = this->base_output_path + "/" + this->sequence_conversion_path;
		this->sequence_conversion_done_file_name = this->sequence_conversion_path + "/sequence_conversion.done";

		this->converted_fastq_file_name = this->sequence_conversion_path + "/#.converted.fastq";
		this->unwanted_fasta_file_name = this->sequence_conversion_path + "/unwanted.fasta";
		this->reference_trim_file_name = this->sequence_conversion_path + "/@.trims";
		this->sequence_conversion_summary_file_name = this->sequence_conversion_path + "/summary.json";

    //! Paths: Read alignment
		this->reference_alignment_path = "02_reference_alignment";
		if (this->base_output_path.size() > 0) this->reference_alignment_path = this->base_output_path + "/" + this->reference_alignment_path;
		this->reference_alignment_done_file_name = this->reference_alignment_path + "/alignment.done";

		this->reference_hash_file_name = this->reference_alignment_path + "/reference";
		this->reference_sam_file_name = this->reference_alignment_path + "/#.reference.sam";

    this->stage1_reference_sam_file_name = this->reference_alignment_path + "/#.stage1.sam";
    this->stage1_unmatched_fastq_file_name = this->reference_alignment_path + "/#.stage1.unmatched.fastq";
    this->stage2_reference_sam_file_name = this->reference_alignment_path + "/#.stage2.matched.sam";

    //! Paths: Junction Prediction
		this->candidate_junction_path = "03_candidate_junctions";
		if (this->base_output_path.size() > 0) this->candidate_junction_path = this->base_output_path + "/" + this->candidate_junction_path;

		this->preprocess_junction_done_file_name = this->candidate_junction_path + "/preprocess_junction_alignment.done";
		this->preprocess_junction_best_sam_file_name = this->candidate_junction_path + "/best.sam";
		this->preprocess_junction_split_sam_file_name = this->candidate_junction_path + "/#.split.sam";
    
    this->coverage_junction_done_file_name = this->candidate_junction_path + "/coverage_junction_alignment.done";
		this->coverage_junction_best_bam_unsorted_file_name = this->candidate_junction_path + "/best.unsorted.bam";
		this->coverage_junction_best_bam_file_name = this->candidate_junction_path + "/best.bam";
		this->coverage_junction_distribution_file_name = this->candidate_junction_path + "/@.unique_only_coverage_distribution.tab";
		this->coverage_junction_plot_file_name = this->candidate_junction_path + "/@.coverage.pdf";
		this->coverage_junction_summary_file_name = this->candidate_junction_path + "/coverage.summary.json";
    this->coverage_junction_error_count_summary_file_name = this->candidate_junction_path + "/error_count.summary.json";

    this->candidate_junction_done_file_name = this->candidate_junction_path + "/candidate_junction.done";
		this->candidate_junction_summary_file_name = this->candidate_junction_path + "/candidate_junction_summary.json";
		this->candidate_junction_fasta_file_name = this->candidate_junction_path + "/candidate_junction.fasta";
    this->candidate_junction_detailed_file_name = this->candidate_junction_path + "/candidate_junction.detailed.txt";
		this->candidate_junction_faidx_file_name = this->candidate_junction_path + "/candidate_junction.fasta.fai";

    //! Paths: Junction Alignment
		this->candidate_junction_alignment_path = "04_candidate_junction_alignment";
		if (this->base_output_path.size() > 0) this->candidate_junction_alignment_path = this->base_output_path + "/" + this->candidate_junction_alignment_path;
		this->candidate_junction_alignment_done_file_name = this->candidate_junction_alignment_path + "/candidate_junction_alignment.done";

		this->candidate_junction_hash_file_name = this->candidate_junction_alignment_path + "/candidate_junction";
		this->candidate_junction_sam_file_name = this->candidate_junction_alignment_path + "/#.candidate_junction.sam";

    //! Paths: Alignment Resolution
		this->alignment_resolution_path = "05_alignment_correction";
		if (this->base_output_path.size() > 0) this->alignment_resolution_path = this->base_output_path + "/" + this->alignment_resolution_path;
		this->alignment_correction_done_file_name = this->alignment_resolution_path + "/alignment_resolution.done";

		this->resolved_reference_sam_file_name = this->alignment_resolution_path + "/reference.sam";
		this->resolved_junction_sam_file_name = this->alignment_resolution_path + "/junction.sam";
		this->alignment_resolution_summary_file_name = this->alignment_resolution_path + "/summary.json";
		this->jc_genome_diff_file_name = this->alignment_resolution_path + "/jc_evidence.gd";

    this->junction_debug_file_name = this->alignment_resolution_path + "/junction_debug.txt";

    //! Paths: BAM conversion
		this->bam_path = "06_bam";
		if (this->base_output_path.size() > 0) this->bam_path = this->base_output_path + "/" + this->bam_path;
		this->bam_done_file_name = this->bam_path + "/bam.done";

		this->reference_bam_unsorted_file_name = this->bam_path + "/reference.unsorted.bam";
		this->junction_bam_unsorted_file_name = this->bam_path + "/junction.unsorted.bam";
		this->junction_bam_file_name = this->bam_path + "/junction.bam";

		//! Paths: Error Calibration
		this->error_calibration_path = "07_error_calibration";
		if (this->base_output_path.size() > 0) this->error_calibration_path = this->base_output_path + "/" + this->error_calibration_path;
		this->error_counts_file_name = this->error_calibration_path + "/error_counts.tab";
		this->error_rates_done_file_name = this->error_calibration_path + "/error_rates.done";

		this->error_rates_file_name = this->error_calibration_path + "/error_rates.tab";
		this->error_counts_done_file_name = this->error_calibration_path + "/error_counts.done";
		this->coverage_file_name = this->error_calibration_path + "/@.coverage.tab";
		this->unique_only_coverage_distribution_file_name = this->error_calibration_path + "/@.unique_only_coverage_distribution.tab";
		this->error_rates_summary_file_name = this->error_calibration_path + "/summary.json";
		this->error_rates_base_qual_error_prob_file_name = this->error_calibration_path + "/base_qual_error_prob.#.tab";
		this->plot_error_rates_r_script_file_name = this->program_data_path + "/plot_error_rate.r";
		this->plot_error_rates_fit_r_script_file_name = this->error_calibration_path + "/fit.#.r_script";
		this->plot_error_rates_r_script_log_file_name = this->error_calibration_path + "/#.plot_error_rate.log";

		//! Paths: Mutation Identification
		this->mutation_identification_path = "08_mutation_identification";
		if (this->base_output_path.size() > 0) this->mutation_identification_path = this->base_output_path + "/" + this->mutation_identification_path;

    this->mutation_identification_done_file_name = this->mutation_identification_path + "/mutation_identification.done";
		this->mutation_identification_per_position_file_name = this->mutation_identification_path + "/per_position_file.tab";
		this->complete_coverage_text_file_name = this->mutation_identification_path + "/@.coverage.tab";
		this->genome_error_counts_file_name = this->mutation_identification_path + "/error_counts.tab";
		this->ra_mc_genome_diff_file_name = this->mutation_identification_path + "/ra_mc_evidence.gd";

    this->polymorphism_statistics_done_file_name = this->mutation_identification_path + "/polymorphism_statistics.done";
		this->polymorphism_statistics_input_file_name = this->mutation_identification_path + "/polymorphism_statistics_input.tab";
		this->polymorphism_statistics_output_file_name = this->mutation_identification_path + "/polymorphism_statistics_output.tab";
		this->polymorphism_statistics_r_script_file_name = this->program_data_path + "/polymorphism_statistics.r";
		this->polymorphism_statistics_r_script_log_file_name = this->mutation_identification_path + "/polymorphism_statistics_output.log";
		this->polymorphism_statistics_ra_mc_genome_diff_file_name = this->mutation_identification_path + "/ra_mc_evidence_polymorphism_statistics.gd";

    //! Paths: Copy Number Variation
		this->copy_number_variation_path = "09_copy_number_variation";
    if (this->base_output_path.size() > 0) this->copy_number_variation_path = this->base_output_path + "/" + this->copy_number_variation_path;
    this->copy_number_variation_done_file_name = this->copy_number_variation_path + "/copy_number_variation.done";
    
    this->tiled_complete_coverage_text_file_name = this->copy_number_variation_path + "/@.tiled.tab";
    this->tiled_for_edging_text_file_name = this->copy_number_variation_path + "/@.tiled_for_edging.tab";
    this->ranges_text_file_name = this->copy_number_variation_path + "/@.ranges.tab";
    this->cnv_history_text_file_name = this->copy_number_variation_path + "/@.history.tab";
    this->smoothed_ranges_text_file_name = this->copy_number_variation_path + "/@.smoothed_ranges.tab";
    this->final_cnv_text_file_name = this->copy_number_variation_path + "/@.cnv_final.tab";
    this->copy_number_variation_cn_genome_diff_file_name = this->copy_number_variation_path + "/@.cn_evidence.gd";
    
    //! Paths: Output
		this->output_path = "output";
		if (this->base_output_path.size() > 0) this->output_path = this->base_output_path + "/" + this->output_path;
		this->output_done_file_name = this->output_path + "/output.done";
    
    this->output_json_summary_file_name = this->output_path + "/summary.json";
    this->output_vcf_file_name = this->output_path + "/output.vcf";

		this->log_file_name = this->output_path + "/log.txt";
		this->index_html_file_name = this->output_path + "/index.html";
		this->summary_html_file_name = this->output_path + "/summary.html";
		this->marginal_html_file_name = this->output_path + "/marginal.html";

		this->local_evidence_path = "evidence";
		this->evidence_path = this->output_path + "/" + this->local_evidence_path;
		this->evidence_genome_diff_file_name = this->evidence_path + "/evidence.gd";
		this->final_genome_diff_file_name = this->output_path + "/output.gd";
    this->preannotated_genome_diff_file_name = this->evidence_path + "/preannotated.gd";
    this->annotated_genome_diff_file_name = this->evidence_path + "/annotated.gd";
    this->output_mutation_prediction_done_file_name = this->evidence_path + "/mutation_prediction.done";
    this->output_mutation_annotation_done_file_name = this->evidence_path + "/mutation_annotation.done";

		this->local_coverage_plot_path = "evidence";
		this->coverage_plot_path = this->output_path + "/" + this->local_coverage_plot_path;
    this->coverage_plot_r_script_file_name = this->program_data_path + "/plot_coverage.r";    
		this->overview_coverage_plot_file_name = this->coverage_plot_path + "/@.overview.png";

		this->output_calibration_path = this->output_path + "/calibration";
		this->unique_only_coverage_plot_file_name = this->output_calibration_path + "/@.unique_coverage.pdf";
		this->error_rates_plot_file_name = this->output_calibration_path + "/#.error_rates.pdf";

		this->breseq_small_graphic_from_file_name = this->program_data_path + "/breseq_small.png";
		this->breseq_small_graphic_to_file_name = this->output_path + "/" + this->local_evidence_path + "/breseq_small.png";
    
    //! Paths: Data
		this->data_path = "data";
		if (this->base_output_path.size() > 0) this->data_path = this->base_output_path + "/" + this->data_path;

		this->reference_bam_file_name = this->data_path + "/reference.bam";
		this->reference_fasta_file_name = this->data_path + "/reference.fasta";
		this->reference_faidx_file_name = this->data_path + "/reference.fasta.fai";
		this->reference_gff3_file_name = this->data_path + "/reference.gff3";
		this->unmatched_read_file_name = this->data_path + "/#.unmatched.fastq";
    this->data_vcf_file_name = this->data_path + "/output.vcf";
    this->data_genome_diff_file_name = this->data_path + "/output.gd";
    this->data_annotated_genome_diff_file_name = this->data_path + "/annotated.gd";
    this->data_json_summary_file_name = this->data_path + "/summary.json";

    //


    //! Paths: Experimental
    this->long_pairs_file_name = this->output_path + "/long_pairs.tab";

    
    //! Multithreading
    Settings::pool.resize(this->num_processors);
  }


	void Settings::init_installed()
	{

    // Save the path for reference
    char * pPath;
    pPath = getenv("PATH");
    this->installed["path"] = "";
    if (pPath!=NULL)
    {
      this->installed["path"] = pPath;
    }
    
    // detect bowtie2 and bowtie2-build system-wide install
    this->installed["bowtie2"] = SYSTEM_CAPTURE("which bowtie2", true);
    this->installed["bowtie2-build"] = SYSTEM_CAPTURE("which bowtie2-build", true);
    this->installed["bowtie2_version"] = "";
    this->installed["bowtie2_version_string"] = "";
    
    if (this->installed["bowtie2"].size() != 0) {
      string version_string = SYSTEM_CAPTURE(this->installed["bowtie2"] + " --version", true);
      size_t start_version_pos = version_string.find("version");
      if (start_version_pos != string::npos) {
        start_version_pos = version_string.find_first_not_of(" \t\r\n", start_version_pos+7);
        size_t end_version_pos = version_string.find_first_not_of("0123456789", start_version_pos+1);
        string new_version_string = version_string.substr(start_version_pos, end_version_pos - start_version_pos);
        
        start_version_pos = version_string.find_first_not_of(".", end_version_pos+1);
        end_version_pos = version_string.find_first_of(" \t\r\n", start_version_pos+1);
        
        new_version_string += "." + version_string.substr(start_version_pos, end_version_pos - start_version_pos);
        this->installed["bowtie2_version_string"] = new_version_string;
        
        // - instead of . appears with beta
        new_version_string = substitute(new_version_string, "-", ".");

        // beta counts as another sub version, just delete that string
        if (new_version_string.find("beta") != string::npos) {
          new_version_string = substitute(new_version_string, "beta", "");
        }
        
        // Split apart string
        vector<string> split_version = split(new_version_string, ".");
        
        // Pad things out to four items
        for(size_t i=split_version.size(); i<4; i++) {
          split_version.push_back("0");
        }
        
        // Only pay attention to four items!
        uint64_t numerical_version = 0;
        for(size_t i=0; i<4; i++) {
          numerical_version *= 1000;
          numerical_version += n(split_version[i]);
        }
        this->installed["bowtie2_version"] = s(numerical_version);
        //cout << this->installed["bowtie2_version"] << endl;
      }
      else {
        this->installed["bowtie2_version_error_message"] = version_string;
      }
    }
    
    //cerr << "Bowtie2: " << this->installed["bowtie2"] << " Bowtie2 version: " << this->installed["bowtie2_version_string"] << endl;

    
		this->installed["R"] = SYSTEM_CAPTURE("which R", true);
    this->installed["R_version_string"] = "";
    this->installed["R_version"] = "";
    
		if (this->installed["R"].size() > 0)
		{
			string R_version = SYSTEM_CAPTURE("R --version", true);
      
      // default if output does not match our pattern
      this->installed["R_version"] = "0";
      
      size_t start_version_pos = R_version.find("R version ");
      if (start_version_pos != string::npos)
      {
        start_version_pos+=10;
        size_t end_version_pos = R_version.find(" ", start_version_pos);
        
        if (end_version_pos == string::npos)
          end_version_pos = R_version.size();
        end_version_pos--;
        
        string version_string = R_version.substr(start_version_pos, end_version_pos - start_version_pos + 1);
        this->installed["R_version_string"] = version_string;
        vector<string> split_version_string = split(version_string, ".");
        if (split_version_string.size() == 3)
        {
          uint32_t R_numerical_version = from_string<uint32_t>(split_version_string[0]) * 1000000 + from_string<uint32_t>(split_version_string[1]) * 1000 + from_string<uint32_t>(split_version_string[2]);
          this->installed["R_version"] = to_string(R_numerical_version);
        }
      }
      else {
        this->installed["R_version_error_message"] = R_version;
        
      }
		}

    //cerr << "R: " << this->installed["R"] << " R version: " << this->installed["R_version_string"] << endl;
    
	}

	void Settings::check_installed()
	{
    // Developer's Note
    //
    // If you are running things through a debugger (like in XCode), your $PATH may need to be
    // set to include the paths where you have Bowtie2 and R installed within your IDE.
    
		bool good_to_go = true;

    
    if (this->installed["bowtie2"].size() == 0)
    {
      good_to_go = false;
      cerr << "---> ERROR Required executable \"bowtie2\" or \"bowtie2-build\" not found." << endl;
      cerr << "---> See http://bowtie-bio.sourceforge.net/bowtie2" << endl;
    }
    else if (this->installed.count("bowtie2_version_error_message")) {
      good_to_go = false;
      cerr << "---> ERROR Could not determine version of installed executable \"bowtie2\" or \"bowtie2-build\"." << endl;
      cerr << "---> For found executable installed at [" << this->installed["bowtie2"] << "]" << endl;
      cerr << "---> Commands \"bowtie2 --version\" or \"bowtie2-build --version\" returned:" << endl;
      cerr << this->installed["bowtie2_version_error_message"] << endl;
    }
    // version encoded in triplets of numbers
    else if (from_string<uint64_t>(this->installed["bowtie2_version"]) < 2000000007) {
      good_to_go = false;
      cerr << "---> ERROR Required executable \"bowtie2\" version 2.0.0-beta7 or later not found." << endl;
      cerr << "---> Your version is " << this->installed["bowtie2_version_string"] << "." << endl;
      cerr << "---> See http://bowtie-bio.sourceforge.net/bowtie2" << endl;
    }
    else if ((from_string<uint64_t>(this->installed["bowtie2_version"]) == 2000003000)
             ||  (from_string<uint64_t>(this->installed["bowtie2_version"]) == 2000004000)) {
      good_to_go = false;
      cerr << "---> ERROR bowtie2 versions 2.0.3 and 2.0.4 are known to have bugs in" << endl;
      cerr << "---> ERROR SAM output that can cause breseq to crash. Please upgrade." << endl;
      cerr << "---> Your version is " << this->installed["bowtie2_version_string"] << "." << endl;
      cerr << "---> See http://bowtie-bio.sourceforge.net/bowtie2" << endl;
    }
    else if (from_string<uint64_t>(this->installed["bowtie2_version"]) < 2001000000) {
      good_to_go = true;
      cerr << "---> WARNING bowtie2 versions before 2.1.0 may produce output that varies slightly" << endl;
      cerr << "---> WARNING from later versions. This may cause consistency tests to fail." << endl;
      cerr << "---> WARNING We strongly suggest that you update to bowtie2 2.1.0+." << endl;
      cerr << "---> Your version is " << this->installed["bowtie2_version_string"] << "." << endl;
      cerr << "---> See http://bowtie-bio.sourceforge.net/bowtie2" << endl;
    }
    else if (from_string<uint64_t>(this->installed["bowtie2_version"]) == 2003001000) {
      good_to_go = false;
      cerr << "---> ERROR bowtie2 version 2.3.1 is known to have a major bug in SAM" << endl;
      cerr << "---> ERROR output that will cause breseq to crash. Please upgrade/downgrade." << endl;
      cerr << "---> Your version is " << this->installed["bowtie2_version_string"] << "." << endl;
      cerr << "---> See http://bowtie-bio.sourceforge.net/bowtie2" << endl;
    }
    else {
      cerr << "---> bowtie2  :: version " << this->installed["bowtie2_version_string"] << " [" << this->installed["bowtie2"] << "]" << endl;
      // WARN if preferred version is not used
      if (from_string<uint64_t>(this->installed["bowtie2_version"]) != 2004005000) {
        cerr << "---> bowtie2  :: NOTE :: breseq output may vary slightly depending on your bowtie2 version," << endl;
        cerr << "---> bowtie2  :: NOTE :: and occasionally bowtie2 versions may have bugs that cause crashes." << endl;
        cerr << "---> bowtie2  :: NOTE :: bowtie2 version 2.4.5 is recommended with this breseq version." << endl;
      }
    }
		// R version 2.1 required
		if (this->installed["R"].size() == 0)
		{
			good_to_go = false;
			cerr << "---> ERROR Required executable \"R\" not found." << endl;
			cerr << "---> See http://www.r-project.org" << endl;
		}
    else if (this->installed.count("R_version_error_message")) {
      good_to_go = false;
      cerr << "---> ERROR Could not determine version of installed executable \"R\"." << endl;
      cerr << "---> For found executable installed at [" << this->installed["R"] << "]" << endl;
      cerr << "---> Command \"R --version\" returned:" << endl;
      cerr << this->installed["R_version_error_message"] << endl;
    }
    else if (from_string<uint32_t>(this->installed["R_version"]) < 2001000)
		{
      good_to_go = false;
      cerr << "---> ERROR Required executable \"R version 2.1.0 or later\" not found." << endl;
      cerr << "---> Your version is " << this->installed["R_version_string"] << "." << endl;
      cerr << "---> See http://www.r-project.org" << endl;
    }
    else {
      cerr << "---> R        :: version " << this->installed["R_version_string"] << " [" << this->installed["R"] << "]" << endl;
    }
    
    /*

    if (this->installed["samtools"].size() == 0)
    {
      good_to_go = false;
      cerr << "---> ERROR Required executable \"samtools\" not found." << endl;
      cerr << "---> This should have been installed by the breseq installer." << endl;
    }
    else if (this->installed.count("samtools_version_error_message")) {
      good_to_go = false;
      cerr << "---> ERROR Could not determine version of installed executable \"samtools\"." << endl;
      cerr << "---> For found executable installed at [" << this->installed["samtools"] << "]" << endl;
      cerr << "---> Command \"samtools\" returned:" << endl;
      cerr << this->installed["samtools_version_error_message"] << endl;
    }
    else {
      cout << "---> samtools :: version " << this->installed["samtools_version_string"] << " [" << this->installed["samtools"] << "]" << endl;
    }
    */
    
		if (!good_to_go) exit(0);
	}

	bool Settings::do_step(const string& done_key, const string& message)
	{
		string done_file_name = done_key;
		this->done_key_messages[done_key] = message;
    this->set_current_step_done_key(done_key);
		if (!file_exists(done_file_name.c_str()))
		{
      cerr << "+++   NOW PROCESSING " << message << endl;
      this->record_start_time(message);
      return true;
		}
    
		cerr << "--- ALREADY COMPLETE " << message << endl;

		ExecutionTime et;
    et.retrieve(done_file_name);
    done_key_intermediate_files = et._done_key_intermediate_files;
    // @JEB Should check for errors, such as incomplete reads...

    this->execution_times.push_back(et);

		return false;
	}

	void Settings::done_step(const string& done_key)
	{
    // Delete intermediate files
    if (!this->keep_all_intermediates) {
      for (vector<string>::iterator it= done_key_intermediate_files[done_key].begin(); it != done_key_intermediate_files[done_key].end(); it++) {
        remove_file(*it, false, true);
      }
    }
    
		// Create the done file with timing information
		string done_file_name = done_key;
		string message = this->done_key_messages[done_key];
    execution_times.back()._done_key_intermediate_files = done_key_intermediate_files;
		this->record_end_time(message);

		this->execution_times.back().store(done_file_name);
	}
  
  
  void Settings::log(const string& message)
  {    
    create_path(this->output_path);
    ofstream LOG(this->log_file_name.c_str(), ios::app);
    ASSERT(!LOG.fail(), "Failed to open log file: " + this->log_file_name);
    LOG << message << endl;
    LOG.close();
  }
  
  // Records an intermediate file that will be cleaned up
  // when a step is done if the user has not specified to 
  // keep intermediate files. Added to current step.
  void Settings::track_intermediate_file(const string& done_key, const string& file_path)
  {
    done_key_intermediate_files[done_key].push_back(file_path);
  }
  
} // breseq namespace

