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


#include "settings.h"

#include "anyoption.h"
#include "genome_diff.h"


using namespace std;

namespace breseq
{
  const int32_t kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation = 50;
  const double kBreseq_ignore_within_this_multiple_of_average_read_length_of_contig_end = 2.0;

  // small mutation <= kBreseq_large_mutation_size_cutoff < large_mutation
  const int32_t kBreseq_large_mutation_size_cutoff = kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation;
  
  const char* kBreseqAlignmentScoreBAMTag = "X5";
  const char* kBreseqBestAlignmentScoreBAMTag = "X6";
  // The two log-likelihood terms of the junction-vs-reference ratio for this read, carried into
  // the BAMs so that the junction read counting in the Output stage can recover the comparison the
  // resolution stage already made. Both are 100 * log10 P (see alignment_log10_likelihood(),
  // reference_sequence.cpp).
  //
  // X8 is this record's OWN likelihood, under the reference it is aligned to. X7 is the best
  // likelihood the same read achieved under the COMPETING hypothesis: on junction.bam that is the
  // reference genome, on reference.bam it is the candidate junctions. X7 is constant across a
  // read's records; X8 is not.
  //
  // These are NOT the alignment scores. Routing still uses AS/X5, which is quality-blind; only the
  // read weighting reads these. Absent tags mean "no competing hypothesis on record" -> full weight
  // toward whichever BAM the read sits in, which is the correct reading when a read aligned to only
  // one of the two references.
  const char* kBreseqOtherHypothesisLogLikelihoodBAMTag = "X7";
  const char* kBreseqOwnHypothesisLogLikelihoodBAMTag = "X8";

  string Settings::global_bin_path;
  string Settings::global_program_data_path;

  //! Multithreading
  ctpl::thread_pool Settings::pool(1);
  std::mutex Settings::lock;

  string Settings::output_divider("================================================================================");
  
  void cReadFileSets::Init(const vector<string>& read_file_names, bool sam_files, bool allow_paired)
  {
    this->clear();
    read_file_to_fastq_file_name_map.clear();
    read_file_to_converted_fastq_file_name_map.clear();

    uint32_t on_id = 0;
    uint32_t on_error_group = 0;
    uint32_t on_paired_end_group = 0;

    set<string> used_base_names;
    set<string> used_original_names;

    // Step 1: build a flat list of cReadFile objects with stripped base names.
    vector<cReadFile> all_files;

    for (vector<string>::const_iterator it = read_file_names.begin(); it < read_file_names.end(); it++)
    {
      cReadFile rf;
      rf.m_original_file_name = *it;

      if (used_original_names.find(rf.m_original_file_name) != used_original_names.end()) {
        WARN("Duplicate read file path provided: " + rf.m_original_file_name + "\nThis file will only be used as input once.");
        continue;
      }
      used_original_names.insert(rf.m_original_file_name);

      rf.m_paired_end_group = on_paired_end_group++;
      rf.m_error_group = on_error_group++;
      rf.m_id = on_id++;

      rf.m_base_name = rf.m_original_file_name;
      size_t pos = rf.m_base_name.rfind("/");
      if (pos != string::npos) rf.m_base_name.erase(0, pos + 1);
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

      string original_base_name = rf.m_base_name;
      uint32_t duplicate_index = 0;
      while (used_base_names.find(rf.m_base_name) != used_base_names.end()) {
        duplicate_index++;
        rf.m_base_name = original_base_name + "_" + to_string(duplicate_index);
      }
      if (rf.m_base_name != original_base_name) {
        WARN("Duplicate read base file name will be renamed in breseq output as follows:\n" + original_base_name + " => " + rf.m_base_name);
      }
      used_base_names.insert(rf.m_base_name);

      all_files.push_back(rf);
    }

    // Step 2: detect paired files (only when --paired-mapping is active).
    if (allow_paired) {
      // Two files form a pair when their base names have equal length and differ in exactly
      // one position where one has '1' and the other has '2'.
      map<string, size_t> base_name_to_index;
      for (size_t i = 0; i < all_files.size(); i++)
        base_name_to_index[all_files[i].m_base_name] = i;

      vector<bool> used(all_files.size(), false);

      for (size_t i = 0; i < all_files.size(); i++) {
        if (used[i]) continue;

        const string& name_a = all_files[i].m_base_name;
        vector<pair<size_t,size_t>> pair_candidates;  // (position_in_name, index_of_partner)

        for (size_t k = 0; k < name_a.size(); k++) {
          if (name_a[k] != '1' && name_a[k] != '2') continue;
          char other_digit = (name_a[k] == '1') ? '2' : '1';
          string candidate = name_a;
          candidate[k] = other_digit;
          auto it = base_name_to_index.find(candidate);
          if (it != base_name_to_index.end() && !used[it->second])
            pair_candidates.push_back({k, it->second});
        }

        if (pair_candidates.size() > 1) {
          WARN("Read file \"" + name_a + "\" could be paired with multiple other read files.\nIt will be treated as unpaired.");
        }

        cReadFileSet rfs;
        if (pair_candidates.size() == 1) {
          size_t k = pair_candidates[0].first;
          size_t j = pair_candidates[0].second;
          string canonical = name_a;
          canonical[k] = 'X';
          rfs.m_base_name = canonical;
          if (name_a[k] == '1') {
            rfs.m_files.push_back(all_files[i]);  // R1
            rfs.m_files.push_back(all_files[j]);  // R2
          } else {
            rfs.m_files.push_back(all_files[j]);  // R1
            rfs.m_files.push_back(all_files[i]);  // R2
          }
          used[i] = true;
          used[j] = true;
        } else {
          rfs.m_base_name = name_a;
          rfs.m_files.push_back(all_files[i]);
          used[i] = true;
        }

        // Populate the base-name → original file map for each individual file in this set.
        for (auto& rf : rfs.m_files)
          read_file_to_fastq_file_name_map[rf.m_base_name] = rf.m_original_file_name;

        this->push_back(rfs);
      }
    } else {
      // Without --paired-mapping, every file is its own single-file set.
      for (auto& rf : all_files) {
        cReadFileSet rfs;
        rfs.m_base_name = rf.m_base_name;
        rfs.m_files.push_back(rf);
        read_file_to_fastq_file_name_map[rf.m_base_name] = rf.m_original_file_name;
        this->push_back(rfs);
      }
    }
  }

  string cReadFileSets::base_name_to_read_file_name(const string& base_name)
  {
    if (read_file_to_converted_fastq_file_name_map.count(base_name))
      return read_file_to_converted_fastq_file_name_map[base_name];
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
		("help,h", "Produce help message showing advanced options (use --expert-help to also show expert options)", TAKES_NO_ARGUMENT)
    ("reference,r", "File containing reference sequences in GenBank, GFF3, or FASTA format. Option may be provided multiple times for multiple files (REQUIRED)")
    ("name,n", "Human-readable name of the analysis run for output (DEFAULT=<none>)", "")
    ("num-processors,j", "Number of processors to use in multithreaded steps", 1)
    //("verbose,v","Produce verbose output",TAKES_NO_ARGUMENT, NORMAL_OPTION) @JEB - not consistently implemented
		("output,o", "Path to breseq output", ".")
    ("polymorphism-prediction,p", "The sample is not clonal. Predict polymorphic (mixed) mutations. Setting this flag changes from CONSENSUS MODE (the default) to POLYMORPHISM MODE", TAKES_NO_ARGUMENT);

    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Read File Options", NORMAL_OPTION);
    options
    ("limit-fold-coverage,l", "Analyze a subset of the input FASTQ sequencing reads with enough bases to provide this theoretical coverage of the reference sequences. A value between 60 and 120 will usually speed up the analysis with no loss in sensitivity for clonal samples. The actual coverage achieved will be somewhat less because not all reads will map (DEFAULT=OFF)", "", NORMAL_OPTION)
    ("nanopore-fast,x", "Set recommended options for nanopore data processed in fast basecalling mode to rule out false-positive mutations caused by high homopolymer indel error rates. Equivalent to --consensus-reject-indel-homopolymer-length 4 --polymorphism-reject-indel-homopolymer-length 4 --polymorphism-no-indels. If you provide any of these options on their own, they override these preset values. Does not change read mapping settings. NOTE: For data processed in high-accuracy basecalling mode, this option is not necessary. (DEFAULT=OFF)", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("aligned-sam", "Input files are aligned SAM files, rather than FASTQ files. Junction prediction steps will be skipped. Be aware that breseq assumes: (1) Your SAM file is sorted such that all alignments for a given read are on consecutive lines. You can use 'samtools sort -n' if you are not sure that this is true for the output of your alignment program. (2) Any 'AS:i:n' alignment scores already present in your SAM file are IGNORED and overwritten. breseq rescores every alignment itself from the CIGAR and the reference, scoring +1 per matching base, -3 per mismatching base, -2 per gap opened and -3 per gap extended, with soft- and hard-clipped bases costing nothing. It does this so that alternative placements of the same read are always compared on one consistent scale. Note this scoring highly penalizes split-read matches (with CIGAR strings such as M35D303M65).", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("read-min-length", "Reads in the input FASTQ file that are shorter than this length will be ignored. (0 = OFF)", 18, NORMAL_OPTION)
    ("read-max-same-base-fraction", "Reads in the input FASTQ file in which this fraction or more of the bases are the same will be ignored. (0 = OFF)", 0.9, NORMAL_OPTION)
    ("read-max-N-fraction", "Reads in the input FASTQ file in which this fraction or more of the bases are uncalled as N will be ignored. (0 = OFF)", 0.5, NORMAL_OPTION)
    ("long-read-trigger-length", "Mark a file as containing long reads and enable read splitting if the longest read has a length that is greater than or equal to this value. (0 = OFF)", 1000, NORMAL_OPTION)
    ("long-read-split-length", "Split input reads in a file marked as having long reads into pieces that are at most this many bases long. Using values much larger than the default for this parameter will likely degrade the speed and accuracy of breseq because of how it performs mapping and analyzes split-read alignments. Filters such as --read-min-length are applied to split reads. (0 = OFF)", 200, NORMAL_OPTION)
    ("long-read-distribute-remainder", "When splitting long reads, divide them into equal pieces that are less than the split length. If this option is not chosen (the default), reads will be split into chunks with exactly the split length and any remaining bases after the last chunk will be ignored.", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("genbank-field-for-seq-id", "Which GenBank header field will be used to assign sequence IDs. Valid choices are LOCUS, ACCESSION, and VERSION. The default is to check those fields, in that order, for the first one that exists. If you override the default, you will need to use the converted reference file (data/reference.gff) for further breseq and gdtools operations on breseq output!", "AUTOMATIC", NORMAL_OPTION)
    ;
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Reference File Options", NORMAL_OPTION);
    options
    ("contig-reference,c", "File containing reference sequences in GenBank, GFF3, or FASTA format. The same coverage distribution will be fit to all of the reference sequences in this file simultaneously. This is appropriate when they are all contigs from a genome that should be present with the same copy number. Use of this option will improve performance when there are many contigs and especially when some are very short (≤1,000 bases).", NULL, NORMAL_OPTION)
    ("junction-only-reference,s", "File containing reference sequences in GenBank, GFF3, or FASTA format. These references are only used for calling junctions with other reference sequences. An example of appropriate usage is including a transposon sequence not present in a reference genome. Option may be provided multiple times for multiple files.", NULL, NORMAL_OPTION)
    ("targeted-sequencing,t", "Reference sequences were targeted for ultra-deep sequencing (using pull-downs or amplicons). Do not fit coverage distribution.", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("user-evidence-gd","User supplied Genome Diff file of JC and/or RA evidence items. The breseq output will report the support for these sequence changes even if they do not pass the normal filters for calling mutations in this sample.", "", NORMAL_OPTION)
    ;
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Read Alignment Options", NORMAL_OPTION);
    options
    ("minimum-mapping-quality,m", "Ignore alignments with less than this mapping quality (MQ) when calling mutations. MQ scores are equal to -10log10(P), where P is the probability that the best alignment is not to the correct location in the reference genome. The range of MQ scores returned by bowtie2 is 0 to 255.", 0, NORMAL_OPTION)
    ("base-quality-cutoff,b", "Ignore bases with quality scores lower than or equal to this value", 3, NORMAL_OPTION)
    ("quality-score-trim", "Trim the ends of reads past any base with a quality score lower than or equal to --base-quality-cutoff.", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("require-match-length", "Only consider alignments that cover this many bases of a read", 0, NORMAL_OPTION)
    ("require-match-fraction", "Only consider alignments that cover this fraction of a read (automatically lowered to 0.5 when --predict-soft-clipping is used, unless set explicitly)", 0.9, NORMAL_OPTION)
    ("maximum-read-mismatches", "Don't consider reads with this many or more bases or indels that are different from the reference sequence. Unaligned bases at the end of a read also count as mismatches. Unaligned bases at the beginning of the read do NOT count as mismatches. (DEFAULT=OFF)", "", NORMAL_OPTION)
    ("junction-indel-split-length", "Threshold (in bases) for splitting an indel found during candidate-junction preprocessing: an indel already present in a single alignment's own CIGAR that is this length or longer is split into separate junction-candidate-supporting alignment records (JC evidence), unless it is entirely a length-change to a reference homopolymer that was already this long or longer (in which case it is left as RA evidence). (DEFAULT = 5)", "", NORMAL_OPTION)
    ;
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Bowtie2 Mapping/Alignment Options", NORMAL_OPTION);
    options
    ("paired-mapping", "DEPRECATED: paired-mapping is now the default. Use --no-paired-mapping to opt out.", TAKES_NO_ARGUMENT, DEPRECATED_OPTION)
    ("no-paired-mapping", "Disable paired-end detection of read files whose names differ only in a '1'/'2' character. By default (paired-mapping ON) such files are grouped as a paired-end set for reporting and paired-aware FASTQ normalization; reads are still aligned independently in single-end mode. Pass this flag to treat every read file as an independent single-end set. (DEFAULT=OFF, i.e. paired-mapping is ON)", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("bowtie2-stage1", "Complete settings (scoring scheme, alignment mode, and mapping criteria) used for the stage 1 alignment. This step is normally meant for quickly aligning near-perfect matches. (DEFAULT=\"" + this->bowtie2_stage1 + "\")", "", EXPERT_OPTION)
    ("bowtie2-stage2", "Complete settings (scoring scheme, alignment mode, and mapping criteria) used for the stage 2 alignment. If set to the empty string \"\", then stage 2 alignment is skipped. This step is normally meant for exhaustively mapping reads that were unmapped by stage 1. (DEFAULT=\"" + this->bowtie2_stage2 + "\")", "", EXPERT_OPTION)
    ("bowtie2-junction", "Complete settings (scoring scheme, alignment mode, and mapping criteria) used in aligning reads to candidate junctions. (DEFAULT=\"" + this->bowtie2_junction + "\")", "", EXPERT_OPTION)
    ("local-mapping", "Use --local (instead of the default --end-to-end) alignment mode for the stage 1 bowtie2 alignment and the candidate-junction realignment pass, restoring the complete --local defaults for --bowtie2-stage1/--bowtie2-junction (--ma 1, with --score-min bounded at 0). --bowtie2-stage2 always uses --local regardless of this option. Has no effect on either setting for which the corresponding --bowtie2-* option is also given explicitly, which always takes precedence. (DEFAULT=OFF, i.e. --end-to-end mapping)", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ;
    options.addUsage("In addition to these values, breseq automatically sets the seed size for bowtie2 read mapping (-L option) to a value that is scaled to the read length (r). This value is 0.5 * r for stage 1 and 0.25 * r for both stage 2 and junction mapping. In each case, it is bounded to the range [9,32]. Be warned that breseq internally rescores alignments with a scoring scheme setting +1 for match, -3 for mismatch, -2 for gap open, and -3 for gap extend for consistency when comparing alternative alignments present in the bowtie2 output.", EXPERT_OPTION);
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Junction (JC) Evidence Options", NORMAL_OPTION);
    options
    ("no-junction-prediction", "Do not predict new sequence junctions", TAKES_NO_ARGUMENT)
    ("junction-alignment-pair-limit", "Only consider this many passed alignment pairs when creating candidate junction sequences (0 = DO NOT LIMIT)", 100000, NORMAL_OPTION)
    ("junction-minimum-candidates", "Test at least this many of the top-scoring junction candidates, regardless of their length", 100, NORMAL_OPTION)
    ("junction-maximum-candidates", "Test no more than this many of the top-scoring junction candidates (0 = DO NOT LIMIT)", 5000, NORMAL_OPTION)
    ("junction-candidate-length-factor", "Accept top-scoring junction candidates to test until their cumulative length is this factor times the total reference sequence length (0 = DO NOT LIMIT)", 0.1, NORMAL_OPTION)
    ("junction-minimum-candidate-pos-hash-score", "Minimum number of distinct spanning read start positions required to create a junction candidate for further testing", 2, NORMAL_OPTION)
    ("junction-score-cutoff", "Maximum negative log10 probability of uneven coverage across a junction breakpoint to accept (0 = OFF)", 3.0, NORMAL_OPTION)
    ("junction-minimum-pos-hash-score", "Minimum number of distinct spanning read start positions required to accept a junction (DEFAULT = consensus mode, 3; polymorphism mode, 3)", "", NORMAL_OPTION)
    ("junction-minimum-side-match", "Minimum number of bases a read must extend past any overlap or read-only sequence at the breakpoint of a junction on each side to count as support for the junction (DEFAULT = 1)", "", NORMAL_OPTION)
    ("junction-minimum-pr-no-read-start-per-position", "Minimum probablilty assigned that no mapped read will start at a given position and strand for junction prediction", 0.1, NORMAL_OPTION)
    ("junction-allow-suboptimal-matches", "Assign a read to the junction candidate with the most overall support as long as its match to this junction is better than to any location in the reference sequence, even if it matches a different junction candidate better. This behavior was the default before v0.35.0. It will align more reads to junctions but risks misassigning some reads to the wrong junction candidates. It is only recommended that you use this option in CONSENSUS mode", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("junction-no-read-weighting", "Count every junction-spanning read as full support, instead of splitting it between the junction and the reference according to how well its alignment scores tell them apart. By default breseq fits that split together with the junction frequency, so a read that barely distinguishes the two counts for less than one that clearly does.", TAKES_NO_ARGUMENT, EXPERT_OPTION)
    ("concordant-pairs-to-make-unique", "When paired-mapping is enabled (the default; disable with --no-paired-mapping), minimum number of concordant read pairs that must agree before a redundant junction side is reassigned to the specific repeat copy their mates map next to (and marked unique). Higher values make false disambiguations rarer. (DEFAULT = 3)", 3, EXPERT_OPTION)
    ;
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Missing Coverage (MC) Evidence Options", NORMAL_OPTION);
    options
    ("deletion-coverage-seed-cutoff","Value for coverage below which MC are seeded", 0, NORMAL_OPTION)
    ("deletion-coverage-propagation-cutoff","Value for coverage above which MC ends stop. 0 = calculated from coverage distribution", 0, NORMAL_OPTION)
    ("call-mutations-overlapping-MC", "If provided, don't ignore mutations predicted from RA evidence that overlap MC evidence", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ;

    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Discordant Pair (DP) Evidence Options", NORMAL_OPTION);
    options
    ("predict-discordant-pairs", "Predict discordant read-pair (DP) evidence for structural variants. This functionality is experimental and OFF by default; requires paired-mapping (the default). (DEFAULT = OFF)", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("discordant-pair-seed", "Minimum number of discordant read pairs within a paired-mapping-distance window required to seed a DP candidate region. (DEFAULT = 3)", 3, NORMAL_OPTION)
    ("discordant-pair-skew-cutoff", "Reference cutoff for the discordant-pair (DP) skew score. The DP skew marginalizes over the empirical concordant-pair crossing distribution when its reference sequence has at least 10^(cutoff+1) non-deletion positions; smaller references fall back to a negative-binomial fit (whose parametric tail is not capped near log10(N), but is conservative). (DEFAULT = 3.0)", 3.0, NORMAL_OPTION)
    ("discordant-pair-minimum-pairs", "Minimum number of read pairs shared by two DP candidate regions for the junction between them to be examined at all. (DEFAULT = 2)", 2, NORMAL_OPTION)
    ("discordant-pair-frequency-cutoff", "Only accept DP evidence when the lower confidence bound on its local variant frequency -- discordant pairs divided by discordant plus concordant pairs spanning the breakpoint -- is above this value. 0 = OFF. Defaults to --polymorphism-frequency-cutoff and follows it if you change it. (DEFAULT = consensus mode, 0.10; polymorphism mode, 0.05)", "", NORMAL_OPTION)
    ("discordant-pair-minimum-crossing", "Only apply the DP skew test when at least this many concordant pairs are expected to span a normal position at the breakpoint's local coverage. Below this the skew test has no power (short paired-mapping distance distributions leave almost no concordant pair spanning any position) and is not used to reject. 0 = always apply. (DEFAULT = 10.0)", 10.0, NORMAL_OPTION)
    // NOTE: pass sub-0.1 defaults as STRINGS. AnyOption stringifies a numeric default with
    // to_string(double), whose default precision is 1 decimal place, so a 0.05 passed as a double
    // silently becomes "0.1".
    ("discordant-pair-background-e-value-cutoff", "Reject DP evidence whose number of supporting read pairs is expected to arise this many times or more across all candidate junctions by chance, given the genome-wide background of spurious discordant pairs. 0 = OFF. (DEFAULT = 0.05)", "0.05", NORMAL_OPTION)
    ;

    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Missing Pair (MP) Evidence Options", NORMAL_OPTION);
    options
    ("predict-missing-pairs", "Predict missing read-pair (MP) evidence: places where reads pile up whose mates did not map anywhere, the signature of a novel sequence inserted into the genome. This functionality is experimental and OFF by default; requires paired-mapping (the default). (DEFAULT = OFF)", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("missing-pair-seed", "Minimum number of reads with unmapped mates within a paired-mapping-distance window required to seed an MP candidate region. (DEFAULT = 3)", 3, NORMAL_OPTION)
    ("missing-pair-seed-fraction", "Minimum fraction of the reads on one strand within a paired-mapping-distance window whose mates did not map, required to seed an MP candidate region. This is what distinguishes a real insertion from the roughly constant genome-wide background of reads whose mates failed to align -- without it a fixed count seeds continuously at any decent coverage and the whole covered region becomes one candidate. (DEFAULT = 0.25)", 0.25, NORMAL_OPTION)
    ("missing-pair-minimum-reads", "Only accept MP evidence supported by at least this many reads with unmapped mates. (DEFAULT = 3)", 3, NORMAL_OPTION)
    ("missing-pair-minimum-distinct", "Only accept MP evidence whose supporting reads start at at least this many distinct positions, so that PCR duplicates of one molecule cannot carry a prediction. (DEFAULT = 2)", 2, NORMAL_OPTION)
    // NOTE: pass sub-0.1 defaults as STRINGS -- see the note in the DP block above. Here the default is
    // mode-dependent, so it is registered empty and only read back when the user actually supplied it.
    ("missing-pair-frequency-cutoff", "Only accept MP evidence when the lower 95% confidence bound on its local variant frequency -- reads with unmapped mates divided by those plus the read pairs spanning the same position -- is at or above this value. 0 = OFF. Defaults to --polymorphism-frequency-cutoff, which tests the same kind of bound on the same kind of frequency, and follows it if you change it. (DEFAULT = consensus mode, 0.10; polymorphism mode, 0.05)", "", NORMAL_OPTION)
    ;

    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Pair Distance (PD) Evidence Options", NORMAL_OPTION);
    options
    ("predict-pair-distance", "Predict pair distance (PD) evidence: places where the read pairs whose unsequenced middle gap spans one point are collectively shifted to longer (a deletion) or shorter (an insertion) mapping distances. DP only sees pairs that are outliers INDIVIDUALLY -- wrong orientation, or a distance above the discordant cutoff -- and has no lower cutoff at all, so events of a few hundred bases are invisible to it in both directions. PD tests the collective shift instead. This functionality is experimental and OFF by default; requires paired-mapping (the default). (DEFAULT = OFF)", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("pair-distance-seed", "Minimum number of read pairs in the matching distribution tail whose gap covers a position, required to seed a PD candidate region. (DEFAULT = 3)", 3, NORMAL_OPTION)
    ("pair-distance-seed-z", "Minimum |z| of the rank-sum statistic over all read pairs whose gap covers a position, required to seed a PD candidate region. This is the 'average quantile per pair' test: under the null each covering pair's distance quantile is uniform, so the mean quantile is a sensitive detector of a shift too small to make any single pair an outlier. 0 = derive from the reference length, so that the genome-wide false-seed rate stays near 1%. (DEFAULT = 0, derived)", 0.0, NORMAL_OPTION)
    ("pair-distance-maximum-span", "Ignore read pairs mapping farther apart than this when seeding PD. Pairs this long are DP's business, and an unbounded distance would make the seeding window unbounded. 0 = derive (twice the paired-mapping distance cutoff). (DEFAULT = 0, derived)", 0, NORMAL_OPTION)
    ("pair-distance-minimum-pairs", "Only accept PD evidence supported by at least this many shifted read pairs. (DEFAULT = 3)", 3, NORMAL_OPTION)
    ("pair-distance-minimum-distinct", "Only accept PD evidence whose supporting pairs start at at least this many distinct positions, so that PCR duplicates of one molecule cannot carry a prediction. (DEFAULT = 2)", 2, NORMAL_OPTION)
    ("pair-distance-minimum-shift", "Only accept PD evidence whose estimated size shift is at least this many bases. 0 = derive from the width of the paired-mapping distance distribution. (DEFAULT = 0, derived)", 0, NORMAL_OPTION)
    // NOTE: pass sub-0.1 defaults as STRINGS -- see the note in the DP block above. Here the default is
    // mode-dependent, so it is registered empty and only read back when the user actually supplied it.
    ("pair-distance-tail-quantile", "Quantile of the paired-mapping distance distribution defining its tails, for the PD seed's counting test. (DEFAULT = 0.05)", "0.05", NORMAL_OPTION)
    ("pair-distance-frequency-cutoff", "Only accept PD evidence when the lower 95% confidence bound on its local variant frequency -- shifted read pairs divided by those plus the normal-distance pairs covering the same point -- is at or above this value. 0 = OFF. Defaults to --polymorphism-frequency-cutoff, which tests the same kind of bound on the same kind of frequency, and follows it if you change it. (DEFAULT = consensus mode, 0.10; polymorphism mode, 0.05)", "", NORMAL_OPTION)
    ;

    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Consensus Read Alignment (RA) Evidence Options", NORMAL_OPTION);
    options
    ("consensus-score-cutoff", "Log10 E-value cutoff for consensus base substitutions and small indels (DEFAULT = 10)", "", NORMAL_OPTION)
    ("consensus-frequency-cutoff", "Only predict consensus mutations when a 95% confidence bound on the variant allele frequency is at or above this value. Which bound is tested depends on the mode, because the two modes assume opposite defaults. In consensus mode the LOWER bound is tested, asking whether the variant is confidently the majority allele. In polymorphism mode the UPPER bound is tested, asking whether a fixed interpretation can still not be ruled out. (DEFAULT = consensus mode, 0.50; polymorphism mode, 0.95)", "", NORMAL_OPTION)
    ("consensus-minimum-variant-coverage", "Only predict consensus mutations when at least this many reads support the mutation. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", NORMAL_OPTION)
    ("consensus-minimum-total-coverage", "Only predict consensus mutations when at least this many reads total are aligned to a genome position. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", NORMAL_OPTION)
    ("consensus-minimum-variant-coverage-each-strand", "Only predict consensus mutations when at least this many reads on each strand support the mutation. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", NORMAL_OPTION)
    ("consensus-minimum-total-coverage-each-strand", "Only predict consensus mutations when at least this many reads on each strand are aligned to a genome position. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", NORMAL_OPTION)
    ("consensus-reject-indel-homopolymer-length", "Reject insertion/deletion polymorphisms which could result from expansion/contraction of homopolymer repeats with this length or greater in the reference genome (0 = OFF) (DEFAULT = consensus mode, OFF; polymorphism mode, OFF) ", "", NORMAL_OPTION)
    ("consensus-reject-surrounding-homopolymer-length", "Reject polymorphic base substitutions that create a homopolymer with this many or more of one base in a row. The homopolymer must begin and end after the changed base. For example, TATTT->TTTTT would be rejected with a setting of 5, but ATTTT->TTTTT would not. (0 = OFF) (DEFAULT = consensus mode, OFF; polymorphism mode, OFF)", "", NORMAL_OPTION)
    ;
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Polymorphism Read Alignment (RA) Evidence Options", NORMAL_OPTION);
    options

    ("polymorphism-score-cutoff", "Log10 E-value cutoff for test of polymorphism vs no polymorphism (DEFAULT = consensus mode, 10; polymorphism mode, 2)", "", NORMAL_OPTION)
    ("polymorphism-frequency-cutoff", "Only predict polymorphisms when the LOWER 95% confidence bound on the minor variant allele frequency is at or above this value, so a variant is never rejected merely for having shallow coverage -- only for being confidently below the cutoff. Variants that cannot be shown to lie below --consensus-frequency-cutoff are predicted as consensus mutations instead. (DEFAULT = consensus mode, 0.10; polymorphism mode, 0.05)", "", NORMAL_OPTION)
    ("polymorphism-minimum-variant-coverage", "Only predict polymorphisms when at least this many reads support each alternative allele. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", NORMAL_OPTION)
    ("polymorphism-minimum-total-coverage", "Only predict polymorphisms when at least this many reads total are aligned to a genome position. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", NORMAL_OPTION)
    ("polymorphism-minimum-variant-coverage-each-strand", "Only predict polymorphisms when at least this many reads on each strand support the non-reference allele. (DEFAULT = 2)", "", NORMAL_OPTION)
    ("polymorphism-minimum-total-coverage-each-strand", "Only predict polymorphisms when at least this many reads on each strand are aligned to a genome position. (DEFAULT = consensus mode, 0; polymorphism mode, 0)", "", NORMAL_OPTION)
    ("polymorphism-strand-bias-cutoff", "P-value criterion for Fisher's exact test for strand bias: reject a polymorphism whose variant reads are distributed across the two strands differently from the total read depth. (0 = OFF) (DEFAULT = 0.05)", "", NORMAL_OPTION)
    ("polymorphism-quality-bias-cutoff", "P-value criterion for the K-S test for base-quality bias: reject a polymorphism whose variant reads have a different base-quality distribution from the reference reads. This also rejects strand-balanced, well-supported calls, so it is OFF by default. (0 = OFF) (DEFAULT = OFF)", "", NORMAL_OPTION)
    ("polymorphism-bias-cutoff", "Shorthand that sets BOTH --polymorphism-strand-bias-cutoff and --polymorphism-quality-bias-cutoff. Either specific option overrides it. (0 = OFF)", "", EXPERT_OPTION)
    ("polymorphism-no-indels", "Do not predict insertion/deletion polymorphisms ≤" + to_string(kBreseq_size_cutoff_AMP_becomes_INS_DEL_mutation) + " bp from read alignment or new junction evidence", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("polymorphism-reject-indel-homopolymer-length", "Reject insertion/deletion polymorphisms which could result from expansion/contraction of homopolymer repeats with this length or greater in the reference genome (0 = OFF) (DEFAULT = OFF) ", "", NORMAL_OPTION)
    ("polymorphism-reject-surrounding-homopolymer-length", "Reject polymorphic base substitutions that create a homopolymer with this many or more of one base in a row. The homopolymer must begin and end after the changed base. For example, TATTT->TTTTT would be rejected with a setting of 5, but ATTTT->TTTTT would not. (0 = OFF) (DEFAULT = OFF)", "", NORMAL_OPTION)
    ;
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Soft Clipping (SC) Evidence Options (HIGHLY EXPERIMENTAL)", NORMAL_OPTION);
    options
    ("predict-soft-clipping", "Predict soft clipping (SC) evidence: positions where reads are unexpectedly soft clipped at their ends, which may indicate unannotated structural variation. This evidence type is experimental and disabled by default.", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("soft-clipping-minimum-bases", "Minimum number of soft-clipped bases at a read end to count as a soft-clipping event. This also sets how much aligned sequence a read must have on BOTH sides of a position to count against soft clipping there, so raising it well beyond the read length leaves no reads to compare against.", 12, NORMAL_OPTION)
    ("soft-clipping-score-cutoff", "Log10 E-value cutoff for soft-clipping evidence (DEFAULT = 3). 0 = OFF.", 3.0, NORMAL_OPTION)
    // NOTE: these are registered with an empty default and only read below when the
    // user actually passed them, so the exact values set in pre_option_initialize()
    // survive. Passing a numeric default here would route it through
    // to_string(double), which formats with only ONE decimal place -- 0.75 would
    // become "0.8" and 0.005 and 1e-5 would both become "0.0" (i.e. silently OFF).
    // This is the same pattern --polymorphism-frequency-cutoff uses.
    ("soft-clipping-consensus-base-fraction", "Fraction of clipped bases that must match for a clipped read to count as sharing the same clipped sequence as the others at its position. Multiplied by the number of bases compared and rounded DOWN, so the default of 0.95 allows one mismatch in 12 bases. Reads clipped for uninteresting reasons (adapter read-through, low-quality tails) differ from each other, while reads spanning a real breakpoint all carry the same donor sequence. (DEFAULT = 0.95). 0 = OFF, count every clipped read.", "", NORMAL_OPTION)
    ("soft-clipping-dispersion-trim-frequency", "Positions where at least this fraction of reads are clipped are excluded from the estimate of the background clipping-rate dispersion, so that genuine breakpoints do not inflate the null they are tested against (DEFAULT = 0.25). 0 = OFF.", "", NORMAL_OPTION)
    ("soft-clipping-maximum-dispersion", "Upper bound on the estimated dispersion (intra-class correlation) of the clipping rate. Caps how conservative the null can become on references too small to estimate it reliably (DEFAULT = 0.005). 0 = OFF, use a plain binomial null.", "", NORMAL_OPTION)
    ("soft-clipping-minimum-rate", "Lower bound on the estimated background clipping rate, preventing an unusually clean library from making every clipped position significant (DEFAULT = 1e-5). 0 = OFF.", "", NORMAL_OPTION)
    ("soft-clipping-minimum-read-count", "Minimum number of consensus-supporting clipped reads required for a position to be reported at all (DEFAULT = 3). 0 = OFF.", 3, NORMAL_OPTION)
    ("soft-clipping-frequency-cutoff", "Minimum fraction of reads that must be clipped at a position for the evidence to be accepted rather than rejected. Defaults to --polymorphism-frequency-cutoff, the same frequency used for predicting mutations, and follows it if you change it (DEFAULT = consensus mode, 0.10; polymorphism mode, 0.05). 0 = OFF.", "", NORMAL_OPTION)
    ("soft-clipping-consensus-fraction-cutoff", "Minimum fraction of the clipped reads at a position that must agree on the consensus clipped sequence for the evidence to be accepted rather than rejected (DEFAULT = 0.5). 0 = OFF.", "", NORMAL_OPTION)
    ;
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Copy numnber (CN) Evidence Options (HIGHLY EXPERIMENTAL)", NORMAL_OPTION);
    options
    ("predict-copy-number","Predict copy number variation evidence using CNery",TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ;
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Output Options", NORMAL_OPTION);
    options
    ("header-genome-diff,g", "Include header information from this GenomeDiff file in output.gd", "", NORMAL_OPTION)
    ("output-unmapped-reads", "Output unmapped reads to file: " + this->unmapped_reads_fastq_file_name, TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("no-evidence-html", "Don't create output files for evidence (e.g., read alignments and coverage plots)", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("zip-html", "DEPRECATED: HTML evidence is now bundled into a ZIP archive (evidence.zip) by default. Use --unzipped-html to opt out.", TAKES_NO_ARGUMENT, DEPRECATED_OPTION)
    ("unzipped-html", "Do not bundle evidence files into a ZIP archive. By default, evidence pages are bundled into evidence.zip and loaded from the archive via JavaScript in the main output files. Pass this flag to instead write the individual evidence files to disk uncompressed. (DEFAULT=OFF, i.e. HTML evidence is zipped)", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("no-javascript", "Don't include javascript in the HTML output", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("max-flanking-columns", "Maximum number of columns in aligned reads to show flanking the region of interest in the HTML output for an evidence item (0=ALL)", 100, NORMAL_OPTION)
    ("max-displayed-reads", "Maximum number of reads to display in the HTML output for an evidence item (0=ALL)", 100, NORMAL_OPTION)
    ;
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Masking Options", NORMAL_OPTION);
    options
    ("mask-gd", "Mask predicted mutations that overlap MASK entries in this GenomeDiff file (mutations are marked ignore=masked, not deleted)", "", NORMAL_OPTION)
    ("mask-mode", "Mode for masking mutations and evidence: 'ALL' masks all mutation and evidence types in masked regions; 'SMALL' masks only small mutations (SNP,DEL,INS,SUB <=" + to_string(kBreseq_large_mutation_size_cutoff) + " bp) and small evidence (RA,MC,CN <=" + to_string(kBreseq_large_mutation_size_cutoff) + " bp); JC and SC evidence are always treated as large and are never masked in SMALL mode", "ALL", NORMAL_OPTION)
    ;

    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Pipeline Control Options", NORMAL_OPTION);
    options
    ("skip-RA-MC-prediction", "Skip generating read alignment and missing coverage evidence.", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("skip-JC-prediction", "Skip generating new junction evidence.", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("skip-MC-prediction", "Skip generating missing coverage evidence.", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("skip-homologous-DEL-prediction", "Skip predicting large deletions between two near-identical copies of an unannotated repeat (paralogous genes, an rRNA operon) from missing coverage evidence. These deletions leave no new junction, so they are located by reading which copy the surviving reads came from.", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ;
    
    options.addUsage("", NORMAL_OPTION);
    options.addUsage("Debugging Options", NORMAL_OPTION);
    options
    ("keep-intermediates,k","Do not delete intermediate files.", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("per-position-file", "Create additional file of per-position aligned bases", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ("junction-debug", "Output additional junction debugging files", TAKES_NO_ARGUMENT, NORMAL_OPTION)
    ;

    // Basic option, deliberately listed last: shows the expert-level help tier.
    options
    ("expert-help,$", "Show the full help message, including expert-level options such as the raw bowtie2 command settings (--bowtie2-stage1/-stage2/-junction)", TAKES_NO_ARGUMENT)
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
    if (options.count("expert-help"))
    {
      options.printExpertUsage();
      exit(-1);
    }
    if (options.count("help"))
    {
      options.printNormalUsage();
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
    this->paired_mapping = !options.count("no-paired-mapping");

    // Discordant-pair (DP) evidence prediction is experimental and opt-in. It needs the paired-mapping
    // insert distribution + candidate regions, so it is ignored (with a warning) under --no-paired-mapping.
    this->predict_discordant_pairs = options.count("predict-discordant-pairs");
    if (this->predict_discordant_pairs && !this->paired_mapping) {
      cerr << "WARNING: --predict-discordant-pairs requires paired-mapping, but --no-paired-mapping was" << endl;
      cerr << "         given. No discordant-pair (DP) evidence will be predicted." << endl;
      cerr << output_divider << endl;
      this->predict_discordant_pairs = false;
    }

    // Missing-pair (MP) evidence prediction is experimental and opt-in. Like DP it needs the
    // paired-mapping insert distribution, and its qualifying reads only exist in a paired run.
    this->predict_missing_pairs = options.count("predict-missing-pairs");
    if (this->predict_missing_pairs && !this->paired_mapping) {
      cerr << "WARNING: --predict-missing-pairs requires paired-mapping, but --no-paired-mapping was" << endl;
      cerr << "         given. No missing-pair (MP) evidence will be predicted." << endl;
      cerr << output_divider << endl;
      this->predict_missing_pairs = false;
    }

    // Pair-distance (PD) evidence prediction is experimental and opt-in. Like DP and MP it needs the
    // paired-mapping insert distribution -- which IS its null model -- so it cannot run without it.
    this->predict_pair_distance = options.count("predict-pair-distance");
    if (this->predict_pair_distance && !this->paired_mapping) {
      cerr << "WARNING: --predict-pair-distance requires paired-mapping, but --no-paired-mapping was" << endl;
      cerr << "         given. No pair-distance (PD) evidence will be predicted." << endl;
      cerr << output_divider << endl;
      this->predict_pair_distance = false;
    }

    // Backward compatibility: --paired-mapping is now the default. Accept it, but warn.
    if (options.count("paired-mapping")) {
      cerr << "WARNING: The --paired-mapping option is DEPRECATED. Paired-mapping detection is now" << endl;
      cerr << "         enabled by default, so this flag has no effect. To turn paired-mapping OFF," << endl;
      cerr << "         use --no-paired-mapping instead." << endl;
      cerr << output_divider << endl;
    }
    this->read_file_sets.Init(read_file_names, this->aligned_sam_mode, this->paired_mapping);
    
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
    
    this->predict_copy_number = options.count("predict-copy-number");

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

    this->discordant_pair_seed = from_string<int32_t>(options["discordant-pair-seed"]);
    ASSERT(this->discordant_pair_seed >= 0, "Argument --discordant-pair-seed must be >= 0")
    this->discordant_pair_skew_cutoff = from_string<double>(options["discordant-pair-skew-cutoff"]);
    ASSERT(this->discordant_pair_skew_cutoff >= 0, "Argument --discordant-pair-skew-cutoff must be >= 0")
    this->discordant_pair_minimum_pairs = from_string<int32_t>(options["discordant-pair-minimum-pairs"]);
    ASSERT(this->discordant_pair_minimum_pairs >= 1, "Argument --discordant-pair-minimum-pairs must be >= 1")
    this->discordant_pair_minimum_crossing = from_string<double>(options["discordant-pair-minimum-crossing"]);
    ASSERT(this->discordant_pair_minimum_crossing >= 0, "Argument --discordant-pair-minimum-crossing must be >= 0")
    this->discordant_pair_background_e_value_cutoff = from_string<double>(options["discordant-pair-background-e-value-cutoff"]);
    ASSERT(this->discordant_pair_background_e_value_cutoff >= 0, "Argument --discordant-pair-background-e-value-cutoff must be >= 0")

    this->missing_pair_seed = from_string<int32_t>(options["missing-pair-seed"]);
    ASSERT(this->missing_pair_seed >= 0, "Argument --missing-pair-seed must be >= 0")
    this->missing_pair_seed_fraction = from_string<double>(options["missing-pair-seed-fraction"]);
    ASSERT((this->missing_pair_seed_fraction >= 0) && (this->missing_pair_seed_fraction <= 1),
           "Argument --missing-pair-seed-fraction must be in the range [0,1]")
    this->missing_pair_minimum_reads = from_string<int32_t>(options["missing-pair-minimum-reads"]);
    ASSERT(this->missing_pair_minimum_reads >= 0, "Argument --missing-pair-minimum-reads must be >= 0")
    this->missing_pair_minimum_distinct = from_string<int32_t>(options["missing-pair-minimum-distinct"]);
    ASSERT(this->missing_pair_minimum_distinct >= 0, "Argument --missing-pair-minimum-distinct must be >= 0")

    this->pair_distance_seed = from_string<int32_t>(options["pair-distance-seed"]);
    ASSERT(this->pair_distance_seed >= 0, "Argument --pair-distance-seed must be >= 0")
    this->pair_distance_tail_quantile = from_string<double>(options["pair-distance-tail-quantile"]);
    ASSERT((this->pair_distance_tail_quantile > 0) && (this->pair_distance_tail_quantile < 0.5),
           "Argument --pair-distance-tail-quantile must be in the range (0,0.5)")
    this->pair_distance_seed_z = from_string<double>(options["pair-distance-seed-z"]);
    ASSERT(this->pair_distance_seed_z >= 0, "Argument --pair-distance-seed-z must be >= 0")
    this->pair_distance_maximum_span = from_string<int32_t>(options["pair-distance-maximum-span"]);
    ASSERT(this->pair_distance_maximum_span >= 0, "Argument --pair-distance-maximum-span must be >= 0")
    this->pair_distance_minimum_pairs = from_string<int32_t>(options["pair-distance-minimum-pairs"]);
    ASSERT(this->pair_distance_minimum_pairs >= 0, "Argument --pair-distance-minimum-pairs must be >= 0")
    this->pair_distance_minimum_distinct = from_string<int32_t>(options["pair-distance-minimum-distinct"]);
    ASSERT(this->pair_distance_minimum_distinct >= 0, "Argument --pair-distance-minimum-distinct must be >= 0")
    this->pair_distance_minimum_shift = from_string<int32_t>(options["pair-distance-minimum-shift"]);
    ASSERT(this->pair_distance_minimum_shift >= 0, "Argument --pair-distance-minimum-shift must be >= 0")

    this->predict_soft_clipping = options.count("predict-soft-clipping");
    this->soft_clipping_minimum_bases = from_string<uint32_t>(options["soft-clipping-minimum-bases"]);
    this->soft_clipping_log10_e_value_cutoff = from_string<double>(options["soft-clipping-score-cutoff"]);
    ASSERT(this->soft_clipping_log10_e_value_cutoff >= 0, "Argument --soft-clipping-score-cutoff must be >= 0")

    // Only override the pre_option_initialize() values when the user actually supplied
    // the option -- see the note on the registrations above.
    if (options.count("soft-clipping-consensus-base-fraction"))
      this->soft_clipping_consensus_base_fraction = from_string<double>(options["soft-clipping-consensus-base-fraction"]);
    ASSERT((this->soft_clipping_consensus_base_fraction >= 0) && (this->soft_clipping_consensus_base_fraction <= 1),
           "Argument --soft-clipping-consensus-base-fraction must be in the range [0,1]")
    if (options.count("soft-clipping-dispersion-trim-frequency"))
      this->soft_clipping_dispersion_trim_frequency = from_string<double>(options["soft-clipping-dispersion-trim-frequency"]);
    ASSERT((this->soft_clipping_dispersion_trim_frequency >= 0) && (this->soft_clipping_dispersion_trim_frequency <= 1),
           "Argument --soft-clipping-dispersion-trim-frequency must be in the range [0,1]")
    if (options.count("soft-clipping-maximum-dispersion"))
      this->soft_clipping_maximum_dispersion = from_string<double>(options["soft-clipping-maximum-dispersion"]);
    ASSERT((this->soft_clipping_maximum_dispersion >= 0) && (this->soft_clipping_maximum_dispersion < 1),
           "Argument --soft-clipping-maximum-dispersion must be in the range [0,1)")
    if (options.count("soft-clipping-minimum-rate"))
      this->soft_clipping_minimum_rate = from_string<double>(options["soft-clipping-minimum-rate"]);
    ASSERT((this->soft_clipping_minimum_rate >= 0) && (this->soft_clipping_minimum_rate < 1),
           "Argument --soft-clipping-minimum-rate must be in the range [0,1)")

    this->soft_clipping_minimum_read_count = from_string<uint32_t>(options["soft-clipping-minimum-read-count"]);
    // NOTE: --soft-clipping-frequency-cutoff is applied further below, after the
    // prediction-mode defaults, because those overwrite it. See the block next to
    // --polymorphism-frequency-cutoff.
    if (options.count("soft-clipping-consensus-fraction-cutoff"))
      this->soft_clipping_consensus_fraction_cutoff = from_string<double>(options["soft-clipping-consensus-fraction-cutoff"]);
    ASSERT((this->soft_clipping_consensus_fraction_cutoff >= 0) && (this->soft_clipping_consensus_fraction_cutoff <= 1),
           "Argument --soft-clipping-consensus-fraction-cutoff must be in the range [0,1]")

    // Soft-clipping evidence needs partially-aligned reads to be kept so their clipped
    // ends are visible for tabulation, so loosen the match-fraction requirement when
    // enabled, unless the user explicitly set their own value.
    if (this->predict_soft_clipping && !options.count("require-match-fraction")) {
      this->require_match_fraction = 0.5;
    }

    //! Settings: bowtie2
    //  all have default that we only overwrite if present on command line
    //
    //  end-to-end is the default mapping mode: it swaps in a different, complete default for
    //  bowtie2-stage1/bowtie2-junction (using end-to-end scoring, including --ma 0 since bowtie2
    //  rejects a nonzero --ma whenever --score-min can go negative) rather than patching the
    //  --local-oriented pristine defaults at the point bowtie2 is actually invoked. The
    //  --local-mapping option disables this and keeps the --local defaults instead. An explicit
    //  --bowtie2-stage1/--bowtie2-junction always wins over either, checked first. bowtie2-stage2
    //  is never affected (always --local, --ma 1).
    this->end_to_end = !options.count("local-mapping");

    if (options.count("bowtie2-stage1")) {
      this->bowtie2_stage1 = options["bowtie2-stage1"];
    } else if (this->end_to_end) {
      string bowtie2_genome_alignment_reporting_parameters = (bowtie2_genome_maximum_alignments_to_consider_per_read > 0)
          ? " -k " + to_string(this->bowtie2_genome_maximum_alignments_to_consider_per_read) : " -a";
      this->bowtie2_stage1 = "--ma 0 --mp 3 --np 0 --rdg 2,3 --rfg 2,3 --ignore-quals --end-to-end -i S,1,0.25 --score-min L,0,-0.2" + bowtie2_genome_alignment_reporting_parameters;
    }

    if (options.count("bowtie2-stage2")) this->bowtie2_stage2 = options["bowtie2-stage2"];

    if (options.count("bowtie2-junction")) {
      this->bowtie2_junction = options["bowtie2-junction"];
    } else if (this->end_to_end) {
      string bowtie2_junction_alignment_reporting_parameters = (bowtie2_junction_maximum_alignments_to_consider_per_read > 0)
          ? " -k " + to_string(this->bowtie2_junction_maximum_alignments_to_consider_per_read) : " -a";
      this->bowtie2_junction = "--ma 0 --mp 3 --np 0 --rdg 2,3 --rfg 2,3 --ignore-quals --end-to-end -i S,1,0.25 --score-min L,1,-0.4" + bowtie2_junction_alignment_reporting_parameters;
    }

    //! Settings: Junction Prediction
    this->skip_new_junction_prediction = this->skip_new_junction_prediction || options.count("no-junction-prediction");
    this->minimum_candidate_junctions = from_string<int32_t>(options["junction-minimum-candidates"]);
    this->maximum_candidate_junctions = from_string<int32_t>(options["junction-maximum-candidates"]);
    this->maximum_candidate_junction_length_factor = from_string<double>(options["junction-candidate-length-factor"]);
    this->minimum_candidate_junction_pos_hash_score = from_string<double>(options["junction-minimum-candidate-pos-hash-score"]);
    this->concordant_pairs_to_make_unique = from_string<uint32_t>(options["concordant-pairs-to-make-unique"]);
    this->maximum_junction_sequence_passed_alignment_pairs_to_consider = from_string<uint64_t>(options["junction-alignment-pair-limit"]);
    this->junction_pos_hash_neg_log10_p_value_cutoff = from_string<double>(options["junction-score-cutoff"]);
    this->junction_allow_suboptimal_matches = options.count("junction-allow-suboptimal-matches");
    this->junction_weight_reads = !options.count("junction-no-read-weighting");
    
    //! Settings: Pipeline Control
    this->skip_read_alignment_and_missing_coverage_prediction = options.count("skip-RA-MC-prediction");
    this->skip_new_junction_prediction = this->skip_new_junction_prediction || options.count("skip-JC-prediction");
    this->skip_missing_coverage_prediction = options.count("skip-MC-prediction");
    this->skip_homologous_deletion_prediction = options.count("skip-homologous-DEL-prediction");
    
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
      // Recalibrated for the Clopper-Pearson tests. The frequency comparisons are now confidence
      // bounds rather than point estimates, which shifts the effective threshold: the lower bound
      // sits below the point estimate, so the floor comes down, and the upper bound sits above it,
      // so the consensus cutoff goes up. Keeping consensus = 1 - floor makes the two-sided
      // polymorphism test exactly symmetric.
      this->consensus_frequency_cutoff = 0.95; // zero is OFF - ensures any rejected poly with high freq move to consensus!
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
      // Polymorphism mode keeps its artifact guards on: it is explicitly hunting low-frequency
      // variants, where a one-strand indel in a homopolymer run is far more likely to be a
      // chemistry artifact than biology. These gate the POLYMORPHISM prediction only -- a call
      // that snaps to consensus above never reaches them.
      this->polymorphism_minimum_variant_coverage_each_strand = 2;
      this->polymorphism_minimum_total_coverage_each_strand = 0;
      
      this->mixed_base_prediction = false;
      // Homopolymer filters OFF pending a statistical replacement. These are hard length cutoffs
      // standing in for a model of why homopolymer runs are indel hotspots; a per-context error
      // model should express that as a rate, not as "reject any indel in a run of 3+". Still
      // available via --polymorphism-reject-indel-homopolymer-length / -surrounding-.
      this->polymorphism_reject_indel_homopolymer_length = 0; // zero is OFF!
      this->polymorphism_reject_surrounding_homopolymer_length = 0; // zero is OFF!
      this->polymorphism_fisher_strand_p_value_cutoff = 0.05;
      this->polymorphism_ks_quality_p_value_cutoff = 0;
      this->polymorphism_no_indels = false;
      this->polymorphism_precision_decimal = 0.000001;
      this->polymorphism_precision_places = 8;
      
      this->minimum_alignment_resolution_pos_hash_score = 3;
      // Was 6, matching consensus mode's 1 now that reads are weighted by how well they actually
      // discriminate. The cutoff was a positional stand-in for that: it measured how far a read
      // reached past the breakpoint, not whether those bases said anything. On 35 bp reads the two
      // correlate at only r = 0.6, so it admitted 9% of reads carrying essentially no evidence
      // while discarding 24% that carried two or more discriminating bases -- and it discarded
      // 31% of all spanning reads to do it. It also widens the read-counting windows by
      // 2*(value-1) (see minimum_side_match_correction), which on short reads distorts the
      // junction count far more than the side counts and biases the frequency itself.
      this->junction_minimum_side_match = 1;
      this->junction_pos_hash_neg_log10_p_value_cutoff = 0; // OFF
      
    }
    // This is strictly true if we are not in polymorphism mode...
    else /* if (this->mixed_base_prediction)*/ {

      this->minimum_mapping_quality = 0;
      
      this->mutation_log10_e_value_cutoff = 10;
      // Majority test: accept as consensus when the 95% lower confidence bound clears this, i.e.
      // we are confident the variant is more than half the reads. Not an upper-bound test.
      this->consensus_frequency_cutoff = 0.50;
      this->consensus_minimum_variant_coverage = 0;
      this->consensus_minimum_total_coverage = 0;
      this->consensus_minimum_variant_coverage_each_strand = 0;
      this->consensus_minimum_total_coverage_each_strand = 0;
      this->consensus_reject_indel_homopolymer_length = 0; // zero is OFF!
      this->consensus_reject_surrounding_homopolymer_length = 0; // zero is OFF!
      
      this->polymorphism_log10_e_value_cutoff = 10;
      this->polymorphism_frequency_cutoff = 0.10;
      this->polymorphism_minimum_variant_coverage = 0;
      this->polymorphism_minimum_total_coverage = 0;
      // These guard the POLYMORPHISM prediction only. Consensus mode asks one question first --
      // is the variant the majority allele -- and answers it from the score and the frequency
      // bound alone; a call that clears that never reaches these. They apply to what falls
      // through, where a ~15% single-strand indel is far more likely to be a chemistry artifact
      // than biology, so the variant must appear on both strands and must not fail Fisher's
      // strand-bias or the K-S base-quality test.
      this->polymorphism_minimum_variant_coverage_each_strand = 2;
      this->polymorphism_minimum_total_coverage_each_strand = 0;

      this->polymorphism_reject_indel_homopolymer_length = 0; // zero is OFF!
      this->polymorphism_reject_surrounding_homopolymer_length = 0; // zero is OFF!
      this->polymorphism_fisher_strand_p_value_cutoff = 0.05;
      // K-S quality bias also rejects strand-balanced, well-supported calls (e.g. a 48% variant
      // with an 80/59 strand split and Fisher p = 0.81), so it is not on by default.
      this->polymorphism_ks_quality_p_value_cutoff = 0;
      this->polymorphism_no_indels = false;
      this->polymorphism_precision_decimal = 0.000001;
      this->polymorphism_precision_places = 3;
      
      this->minimum_alignment_resolution_pos_hash_score = 3;
      this->junction_minimum_side_match = 1;
    }

    // en masse overrides come before specific overrides
    //
    // Must be after the prediction-mode blocks above, which reset these to their OFF values,
    // and before the specific overrides below, so that passing one of these options on its own
    // still wins over the preset.

    if (options.count("nanopore-fast")) {
      this->consensus_reject_indel_homopolymer_length = 4;
      this->polymorphism_reject_indel_homopolymer_length = 4;
      this->polymorphism_no_indels = true;
    }

    // override the default settings

    if (options.count("minimum-mapping-quality")) {
      this->minimum_mapping_quality = from_string<int32_t>(options["minimum-mapping-quality"]);
    }
    
    if (options.count("junction-indel-split-length") && (options["junction-indel-split-length"] != ""))
      this->junction_indel_split_length = from_string<int32_t>(options["junction-indel-split-length"]);

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

    // Every structural-evidence frequency gate FOLLOWS the polymorphism cutoff rather than carrying a
    // floor of its own, because they all measure the same quantity it governs: the variant allele
    // frequency at one locus -- the fraction of the molecules sampled there that carry the variant
    // (soft-clipped reads for SC, discordant pairs for DP, reads that lost their mate to novel
    // sequence for MP, distance-shifted pairs for PD). All are lower-bound tests on that frequency,
    // so a user who says "I care about variants down to X" gets structural evidence down to X as
    // well, and one number governs where the whole run's sensitivity sits.
    //
    // SC's and DP's floors used to be a separate 0.2 in consensus mode. That was a holdover: both were
    // tuned before either tested a Clopper-Pearson bound (SC against a point estimate, DP against an
    // exact test at a different scale), so neither number was calibrated for the test it now feeds.
    //
    // This must sit AFTER the option read just above, so an explicit --polymorphism-frequency-cutoff
    // is what gets tracked, and BEFORE the per-type option reads below, which still override it.
    this->soft_clipping_frequency_cutoff   = this->polymorphism_frequency_cutoff;
    this->discordant_pair_frequency_cutoff = this->polymorphism_frequency_cutoff;
    this->missing_pair_frequency_cutoff    = this->polymorphism_frequency_cutoff;
    this->pair_distance_frequency_cutoff   = this->polymorphism_frequency_cutoff;

    // Must come after the prediction-mode defaults above, which set this to track
    // polymorphism_frequency_cutoff and would otherwise clobber an explicit value.
    if (options.count("soft-clipping-frequency-cutoff"))
      this->soft_clipping_frequency_cutoff = from_string<double>(options["soft-clipping-frequency-cutoff"]);
    ASSERT((this->soft_clipping_frequency_cutoff >= 0) && (this->soft_clipping_frequency_cutoff <= 1),
           "Argument --soft-clipping-frequency-cutoff must be in the range [0,1]")

    // Same ordering requirement as soft-clipping: after the prediction-mode defaults.
    if (options.count("discordant-pair-frequency-cutoff"))
      this->discordant_pair_frequency_cutoff = from_string<double>(options["discordant-pair-frequency-cutoff"]);
    ASSERT((this->discordant_pair_frequency_cutoff >= 0) && (this->discordant_pair_frequency_cutoff <= 1),
           "Argument --discordant-pair-frequency-cutoff must be in the range [0,1]")

    // Same ordering requirement again: after the prediction-mode defaults.
    if (options.count("missing-pair-frequency-cutoff"))
      this->missing_pair_frequency_cutoff = from_string<double>(options["missing-pair-frequency-cutoff"]);
    ASSERT((this->missing_pair_frequency_cutoff >= 0) && (this->missing_pair_frequency_cutoff <= 1),
           "Argument --missing-pair-frequency-cutoff must be in the range [0,1]")

    // Same ordering requirement again: after the prediction-mode defaults.
    if (options.count("pair-distance-frequency-cutoff"))
      this->pair_distance_frequency_cutoff = from_string<double>(options["pair-distance-frequency-cutoff"]);
    ASSERT((this->pair_distance_frequency_cutoff >= 0) && (this->pair_distance_frequency_cutoff <= 1),
           "Argument --pair-distance-frequency-cutoff must be in the range [0,1]")


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
    if (options.count("polymorphism-bias-cutoff")) {
      double both = from_string<double>(options["polymorphism-bias-cutoff"]);
      this->polymorphism_fisher_strand_p_value_cutoff = both;
      this->polymorphism_ks_quality_p_value_cutoff = both;
    }
    if (options.count("polymorphism-strand-bias-cutoff"))
      this->polymorphism_fisher_strand_p_value_cutoff = from_string<double>(options["polymorphism-strand-bias-cutoff"]);
    if (options.count("polymorphism-quality-bias-cutoff"))
      this->polymorphism_ks_quality_p_value_cutoff = from_string<double>(options["polymorphism-quality-bias-cutoff"]);

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
    this->no_evidence_html = options.count("no-evidence-html");
    this->zip_html = !options.count("unzipped-html");

    // Backward compatibility: --zip-html is now the default. Accept it, but warn.
    if (options.count("zip-html")) {
      cerr << "WARNING: The --zip-html option is DEPRECATED. Bundling HTML evidence into a ZIP" << endl;
      cerr << "         archive is now enabled by default, so this flag has no effect. To turn" << endl;
      cerr << "         ZIP bundling OFF, use --unzipped-html instead." << endl;
      cerr << output_divider << endl;
    }
    this->no_javascript = options.count("no-javascript");
    this->output_unmapped_reads = options.count("output-unmapped-reads");
    
    if (options.count("header-genome-diff"))
      this->header_genome_diff_file_name = options["header-genome-diff"];

    if (options.count("mask-gd"))
      this->mask_genome_diff_file_name = options["mask-gd"];

    if (this->mask_genome_diff_file_name.empty()) {
      this->mask_mode = "NONE";
    } else {
      this->mask_mode = to_upper(options["mask-mode"]);
      if (this->mask_mode != "ALL" && this->mask_mode != "SMALL") {
        options.addUsage("");
        options.addUsage("Invalid --mask-mode value '" + options["mask-mode"] + "'. Valid choices are 'all' and 'small'.");
        options.printUsage();
        exit(-1);
      }
      // Fail fast: validate the mask GD file is readable before the pipeline starts
      cGenomeDiff mask_gd;
      mask_gd.read(this->mask_genome_diff_file_name);
    }

		this->post_option_initialize();
    
    //////// Check here for any conflicting options

    if (this->zip_html && this->no_javascript) {
      cerr << "WARNING: HTML evidence is not being bundled into a ZIP archive because --no-javascript was specified (ZIP bundling requires JavaScript)." << endl;
      this->zip_html = false;
    }

    if (this->zip_html && !file_exists((this->program_data_path + "/jszip.min.js").c_str())) {
      cerr << "WARNING: HTML evidence is not being bundled into a ZIP archive because a required data file (jszip.min.js) is missing from this breseq installation's data path: " << this->program_data_path << endl;
      this->zip_html = false;
    }

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
    cerr << color_green(output_divider) << endl;
    cerr << color_green(string(PACKAGE_STRING) + "  " + GITHUB_REVISION_STRING + "   " + PACKAGE_URL) << endl;
    fprintf(stderr, "\n");
    fprintf(stderr, "Active Developers: Barrick JE\n");
    fprintf(stderr, "Issues/Bugs: %s\n", PACKAGE_BUGREPORT);
    fprintf(stderr, "\n");
    fprintf(stderr, "License GPLv2+ (GPL-2.0-or-later) https://gnu.org/licenses/gpl.html\n");
    fprintf(stderr, "This is free software: you are free to change and redistribute it.\n");
    fprintf(stderr, "There is NO WARRANTY, to the extent permitted by law.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Copyright (c) 2008-2010 Michigan State University\n");
    fprintf(stderr, "Copyright (c) 2011-2025 The University of Texas at Austin\n");
    fprintf(stderr, "Copyright (c) 2025-     Michigan State University\n");
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
    cerr << color_green(output_divider) << endl;
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
    this->paired_mapping = true;
    this->predict_discordant_pairs = false;
    this->predict_missing_pairs = false;
    this->predict_pair_distance = false;

    //! Options that control which parts of the pipeline to execute
    this->skip_read_filtering = false;
    this->skip_new_junction_prediction = false;
		this->skip_read_alignment_and_missing_coverage_prediction = false;
		this->skip_missing_coverage_prediction = false;
    this->skip_homologous_deletion_prediction = false;
    this->no_evidence_html = false;
		this->predict_copy_number = false;
		this->do_periodicity = false;
    
    //! DEBUG options
    this->verbose = false;
		this->alignment_read_limit = 0;         
		this->resolve_alignment_read_limit = 0; 
    this->candidate_junction_read_limit = 0;
    this->output_unmapped_reads = false;
    this->keep_all_intermediates = false;

    //! Settings: Read Alignment and Candidate Junction Read Alignment
    this->require_match_length = 0;         
    this->require_match_fraction = 0.9; // CL value overrides
    this->maximum_read_mismatches = -1;

    this->bowtie2_junction_maximum_alignments_to_consider_per_read = 2000;
    this->bowtie2_genome_maximum_alignments_to_consider_per_read = 2000;

    // Different values for these:
    // We report ALL alignments for junctions
    string bowtie2_junction_alignment_reporting_parameters = (bowtie2_junction_maximum_alignments_to_consider_per_read > 0)
    ? " -k " + to_string(this->bowtie2_junction_maximum_alignments_to_consider_per_read) : " -a";

    // Need to report all matches for proper handling of
    string bowtie2_genome_alignment_reporting_parameters = (bowtie2_genome_maximum_alignments_to_consider_per_read > 0)
    ? " -k " + to_string(this->bowtie2_genome_maximum_alignments_to_consider_per_read) : " -a";

    // Note: this leaves off -L, since it is set based on read length

    // Set these defaults. Each is a complete, independent command-line string -- scoring
    // scheme bundled with alignment mode and mapping cutoffs -- rather than sharing one
    // scoring string across all three, since stage 2/junction always stay --local (--ma 1)
    // even when the default end-to-end mode switches stage 1/junction's own scoring during
    // command-line processing (see end_to_end above). These are the pristine --local baseline
    // defaults, kept as-is when --local-mapping is given.
    this->bowtie2_stage1 = "--ma 1 --mp 3 --np 0 --rdg 2,3 --rfg 2,3 --ignore-quals --local -i S,1,0.25 --score-min L,0,0.8" + bowtie2_genome_alignment_reporting_parameters; // "-L 22 -i S,1,0.25 --score-min L,0,0.9";
    this->end_to_end = false;

    this->bowtie2_stage2 = "--ma 1 --mp 3 --np 0 --rdg 2,3 --rfg 2,3 --ignore-quals --local -i S,1,0.25 --score-min L,6,0.2" + bowtie2_genome_alignment_reporting_parameters; // "-L 9 -i C,1,0 --score-min L,6,0.2";

    this->bowtie2_junction = "--ma 1 --mp 3 --np 0 --rdg 2,3 --rfg 2,3 --ignore-quals --local -i S,1,0.25 --score-min L,1,0.6" + bowtie2_junction_alignment_reporting_parameters; // "-L 9 -i C,1,0 --score-min L,6,0.2";

    
    this->num_processors = 1;
    
    //! Settings: Candidate Junction Prediction
		this->junction_indel_split_length = 5;
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
    this->concordant_pairs_to_make_unique = 3;
    this->penalize_negative_junction_overlap = true;

    //! Settings: Alignment Resolution
    this->add_split_junction_sides = true;
    this->junction_pos_hash_neg_log10_p_value_cutoff = 3;
    this->minimum_alignment_resolution_pos_hash_score = 3;
    this->minimum_pr_no_read_start_per_position = 0.1;
    this->junction_allow_suboptimal_matches = false;
    this->junction_weight_reads = true;

    //! Settings: Mutation Identification
    this->base_quality_cutoff = 3;
    this->quality_score_trim = 0;
    this->deletion_coverage_propagation_cutoff = 0;
    this->deletion_coverage_seed_cutoff = 0;
    this->call_mutations_overlapping_missing_coverage = false;
    this->discordant_pair_seed = 3;
    this->discordant_pair_skew_cutoff = 3.0;
    this->discordant_pair_minimum_pairs = 2;
    // Overwritten per prediction mode in the cmdline constructor to track polymorphism_frequency_cutoff.
    this->discordant_pair_frequency_cutoff = 0.1;
    // For a Poisson null with mean lambda, even k=0 only reaches -log10 P = lambda/ln(10), so the skew
    // test cannot clear a cutoff of 3 on merit below lambda ~ 6.9, and cannot reject k=1 below ~9.6.
    this->discordant_pair_minimum_crossing = 10.0;
    this->discordant_pair_background_e_value_cutoff = 0.05;

    this->missing_pair_seed = 3;
    this->missing_pair_seed_fraction = 0.25;
    this->missing_pair_minimum_reads = 3;
    this->missing_pair_minimum_distinct = 2;
    // Overwritten per prediction mode in the cmdline constructor to track polymorphism_frequency_cutoff.
    this->missing_pair_frequency_cutoff = 0.1;

    this->pair_distance_seed = 3;
    this->pair_distance_tail_quantile = 0.05;
    this->pair_distance_seed_z = 0.0;        // 0 = derive from the reference length
    this->pair_distance_maximum_span = 0;    // 0 = derive (2 * paired-mapping distance cutoff)
    this->pair_distance_minimum_pairs = 3;
    this->pair_distance_minimum_distinct = 2;
    this->pair_distance_minimum_shift = 0;   // 0 = derive from the distribution width
    // Overwritten per prediction mode in the cmdline constructor to track polymorphism_frequency_cutoff.
    this->pair_distance_frequency_cutoff = 0.1;

    this->predict_soft_clipping = false;
    this->soft_clipping_minimum_bases = 12;
    this->soft_clipping_log10_e_value_cutoff = 3.0;
    this->soft_clipping_consensus_base_fraction = 0.95;
    this->soft_clipping_dispersion_trim_frequency = 0.25;
    this->soft_clipping_maximum_dispersion = 0.005;
    this->soft_clipping_minimum_rate = 1e-5;
    this->soft_clipping_minimum_read_count = 3;
    // Overwritten per prediction mode in the cmdline constructor to track polymorphism_frequency_cutoff.
    this->soft_clipping_frequency_cutoff = 0.1;
    this->soft_clipping_consensus_fraction_cutoff = 0.5;

    this->polymorphism_prediction = false;
    this->mixed_base_prediction = true;
    
    this->consensus_minimum_variant_coverage = 0;
    this->consensus_minimum_total_coverage = 0;
    this->consensus_minimum_variant_coverage_each_strand = 0;
    this->consensus_minimum_total_coverage_each_strand = 0;
    
    this->polymorphism_log10_e_value_cutoff = this->mutation_log10_e_value_cutoff;
		this->polymorphism_fisher_strand_p_value_cutoff = 0;
		this->polymorphism_ks_quality_p_value_cutoff = 0;
		this->polymorphism_frequency_cutoff = 0.1;
    // Was never initialized here, which left it indeterminate for the no-command-line Settings used
    // by gdtools. The frequency tests branch on it, so give it the consensus-mode default.
    this->consensus_frequency_cutoff = 0.90;
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
    this->max_rejected_soft_clipping_evidence_to_show = 20;
    this->max_rejected_discordant_pair_evidence_to_show = 20;
    this->max_rejected_missing_pair_evidence_to_show = 20;
    this->max_rejected_pair_distance_evidence_to_show = 20;
		this->hide_circular_genome_junctions = true;
    
		this->lenski_format = false;
		this->no_evidence = false;
    this->no_javascript = false;
    this->zip_html = true;
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

		this->converted_fastq_file_name = this->sequence_conversion_path + "/#.converted.fastq.gz";
		this->unwanted_fasta_file_name = this->sequence_conversion_path + "/unwanted.fasta";
		this->reference_trim_file_name = this->sequence_conversion_path + "/@.trims";
		this->sequence_conversion_summary_file_name = this->sequence_conversion_path + "/summary.json";

    //! Paths: Read alignment
		this->reference_alignment_path = "02_reference_alignment";
		if (this->base_output_path.size() > 0) this->reference_alignment_path = this->base_output_path + "/" + this->reference_alignment_path;
		this->reference_alignment_done_file_name = this->reference_alignment_path + "/alignment.done";

		this->reference_hash_file_name = this->reference_alignment_path + "/reference";
		this->reference_sam_file_name = this->reference_alignment_path + "/#.reference.bam";

    this->stage1_reference_sam_file_name = this->reference_alignment_path + "/#.stage1.bam";
    this->stage1_unmapped_reads_fastq_file_name = this->reference_alignment_path + "/#.stage1.unmapped.fastq";
    // TODO: Enable when bowtie2 updated
    //this->stage1_unmapped_reads_fastq_file_name = this->reference_alignment_path + "/#.stage1.unmapped.fastq.gz";

    this->stage2_reference_sam_file_name = this->reference_alignment_path + "/#.stage2.mapped.bam";

    //! Paths: Junction Prediction
		this->candidate_junction_path = "03_candidate_junctions";
		if (this->base_output_path.size() > 0) this->candidate_junction_path = this->base_output_path + "/" + this->candidate_junction_path;

		this->preprocess_junction_done_file_name = this->candidate_junction_path + "/preprocess_junction_alignment.done";
		this->preprocess_junction_best_sam_file_name = this->candidate_junction_path + "/best.bam";
		this->preprocess_junction_split_sam_file_name = this->candidate_junction_path + "/#.split.bam";
		// Sidecar for --paired-mapping: for each split read whose mate is NOT split, records the mate's
		// alignment positions so candidate-junction identification can pin a redundant side to the copy
		// its mate maps concordantly next to. One per read file, matching the two split.bam files.
		this->preprocess_junction_split_pair_positions_file_name = this->candidate_junction_path + "/#.split_pair_positions.csv";

    this->paired_mapping_distance_done_file_name = this->candidate_junction_path + "/paired_mapping_distance.done";
    this->paired_mapping_distance_summary_file_name = this->candidate_junction_path + "/paired_mapping_distance.summary.json";
    // Transient intermediate: written during preprocessing, consumed by the distance-distribution
    // fit/plot, then deleted (tracked under paired_mapping_distance_done_file_name).
    this->paired_mapping_distance_distribution_file_name = this->candidate_junction_path + "/#.pair_stats.csv";
    // (paired_mapping_distance_histogram_file_name is assigned below, after data_path is set.)

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
    this->candidate_junction_trim_file_name = this->candidate_junction_path + "/@.trims";

    //! Paths: Junction Alignment
		this->candidate_junction_alignment_path = "04_candidate_junction_alignment";
		if (this->base_output_path.size() > 0) this->candidate_junction_alignment_path = this->base_output_path + "/" + this->candidate_junction_alignment_path;
		this->candidate_junction_alignment_done_file_name = this->candidate_junction_alignment_path + "/candidate_junction_alignment.done";

		this->candidate_junction_hash_file_name = this->candidate_junction_alignment_path + "/candidate_junction";
		this->candidate_junction_sam_file_name = this->candidate_junction_alignment_path + "/#.candidate_junction.bam";

    //! Paths: Alignment Resolution
		this->alignment_resolution_path = "05_alignment_correction";
		if (this->base_output_path.size() > 0) this->alignment_resolution_path = this->base_output_path + "/" + this->alignment_resolution_path;
		this->alignment_correction_done_file_name = this->alignment_resolution_path + "/alignment_resolution.done";

		this->resolved_reference_sam_file_name = this->alignment_resolution_path + "/reference.bam";
		this->resolved_junction_sam_file_name = this->alignment_resolution_path + "/junction.bam";
		this->alignment_resolution_summary_file_name = this->alignment_resolution_path + "/summary.json";
		this->resolved_paired_mapping_distance_summary_file_name = this->alignment_resolution_path + "/paired_mapping_distance.summary.json";
		this->jc_genome_diff_file_name = this->alignment_resolution_path + "/jc_evidence.gd";

    this->junction_debug_file_name = this->alignment_resolution_path + "/junction_debug.txt";
    this->junction_mapq_debug_file_name = this->alignment_resolution_path + "/junction_mapq_debug.tsv";

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
		this->soft_clipping_counts_file_name = this->error_calibration_path + "/soft_clipping_counts.tab";
		this->soft_clipping_summary_file_name = this->error_calibration_path + "/soft_clipping_summary.json";
		this->error_rates_base_qual_error_prob_file_name = this->error_calibration_path + "/base_qual_error_prob.#.tab";

		//! Paths: Mutation Identification
		this->mutation_identification_path = "08_mutation_identification";
		if (this->base_output_path.size() > 0) this->mutation_identification_path = this->base_output_path + "/" + this->mutation_identification_path;

    this->mutation_identification_done_file_name = this->mutation_identification_path + "/mutation_identification.done";
		this->mutation_identification_per_position_file_name = this->mutation_identification_path + "/per_position_file.tab";
		this->complete_coverage_text_file_name = this->mutation_identification_path + "/@.coverage.tab";
		this->ra_mc_genome_diff_file_name = this->mutation_identification_path + "/ra_mc_evidence.gd";
		this->dp_candidate_regions_file_name = this->mutation_identification_path + "/DP_candidate_regions.csv";
		this->discordant_pair_done_file_name = this->mutation_identification_path + "/discordant_pair.done";
		this->dp_genome_diff_file_name = this->mutation_identification_path + "/dp_evidence.gd";
		this->discordant_pair_summary_file_name = this->mutation_identification_path + "/discordant_pair_summary.json";
		this->mp_candidate_regions_file_name = this->mutation_identification_path + "/MP_candidate_regions.csv";
		this->missing_pair_done_file_name = this->mutation_identification_path + "/missing_pair.done";
		this->mp_genome_diff_file_name = this->mutation_identification_path + "/mp_evidence.gd";
		this->pd_candidate_regions_file_name = this->mutation_identification_path + "/PD_candidate_regions.csv";
		this->pair_distance_done_file_name = this->mutation_identification_path + "/pair_distance.done";
		this->pd_genome_diff_file_name = this->mutation_identification_path + "/pd_evidence.gd";
		this->pair_distance_summary_file_name = this->mutation_identification_path + "/pair_distance_summary.json";


    //! Paths: Copy Number Variation
		this->copy_number_variation_path = "09_copy_number_variation";
    if (this->base_output_path.size() > 0) this->copy_number_variation_path = this->base_output_path + "/" + this->copy_number_variation_path;
    this->copy_number_variation_done_file_name = this->copy_number_variation_path + "/copy_number_variation.done";
    this->copy_number_evidence_genome_diff_file_name = this->copy_number_variation_path + "/@.cn_evidence.gd";
    this->cnery_output_path = this->copy_number_variation_path + "/cnery_out";
    // CNery names every file it writes <basename of its -o argument><seq_id>_<what>, with no
    // separator anywhere, so its output prefix is the basename of the directory we hand it.
    {
      const string cnery_file_prefix = path_to_filename(this->cnery_output_path) + "@";
      this->cnery_cnv_csv_file_name = this->cnery_output_path + "/CNV_csv/" + cnery_file_prefix + "_CNV.csv";
      this->cnery_break_points_file_name = this->cnery_output_path + "/CNV_csv/" + cnery_file_prefix + "_break_pts.csv";
      this->cnery_otr_results_file_name = this->cnery_output_path + "/OTR_corr/" + cnery_file_prefix + "_otr_results.json";
    }

    //! Paths: Output
		this->output_path = "output";
		if (this->base_output_path.size() > 0) this->output_path = this->base_output_path + "/" + this->output_path;
		this->output_done_file_name = this->output_path + "/output.done";
    
		this->log_file_name = this->output_path + "/log.txt";
		this->index_html_file_name = this->output_path + "/index.html";
		this->summary_html_file_name = this->output_path + "/summary.html";
		this->marginal_html_file_name = this->output_path + "/marginal.html";

		this->local_evidence_path = "evidence";
		this->evidence_path = this->output_path + "/" + this->local_evidence_path;
    this->local_html_archive_file_name = "evidence.zip";
    this->html_archive_file_name = this->output_path + "/" + this->local_html_archive_file_name;
    this->local_html_evidence_file_name = "evidence.html";
    this->html_evidence_file_name = this->output_path + "/" + this->local_html_evidence_file_name;
		this->evidence_genome_diff_file_name = this->evidence_path + "/evidence.gd";
		this->final_genome_diff_file_name = this->output_path + "/output.gd";
    this->preannotated_genome_diff_file_name = this->evidence_path + "/preannotated.gd";
    this->output_mutation_prediction_done_file_name = this->evidence_path + "/mutation_prediction.done";
    this->output_mutation_annotation_done_file_name = this->evidence_path + "/mutation_annotation.done";

		this->local_coverage_plot_path = "evidence";
		this->coverage_plot_path = this->output_path + "/" + this->local_coverage_plot_path;
		this->overview_coverage_plot_file_name = this->coverage_plot_path + "/@.overview.svg";
		this->cn_overview_coverage_plot_file_name = this->coverage_plot_path + "/@.cn_overview.svg";

		this->output_calibration_path = this->evidence_path;
		this->unique_only_coverage_plot_file_name = this->output_calibration_path + "/@.unique_coverage.svg";
		this->concordant_pair_crossing_plot_file_name = this->output_calibration_path + "/concordant_pair_crossing.svg";
		this->concordant_pair_crossing_seq_plot_file_name = this->output_calibration_path + "/#.concordant_pair_crossing.svg";
		this->error_rates_plot_file_name = this->output_calibration_path + "/#.error_rates.svg";
		this->paired_mapping_distance_plot_file_name = this->output_calibration_path + "/#.pair_distance.svg";
		this->discordant_pairs_plot_file_name = this->output_calibration_path + "/discordant_pairs.svg";

		this->breseq_icon_graphic_from_file_name = this->program_data_path + "/breseq_icon.png";
		this->breseq_icon_graphic_to_file_name = this->output_path + "/breseq_icon.png";
    
    //! Paths: Data
		this->data_path = "data";
		if (this->base_output_path.size() > 0) this->data_path = this->base_output_path + "/" + this->data_path;

		// Per-seq_id concordant-pair crossing histogram tab (written in resolve, persists in data/).
		this->concordant_pair_crossing_distribution_file_name = this->data_path + "/#.concordant_pair_crossing.tab";

		// Per-read-file-set majority-orientation insert histogram (the empirical insert PMF for the DP Bayes
		// placement test), written by the fit in preprocessing and read back in Output; persists in data/.
		this->paired_mapping_distance_histogram_file_name = this->data_path + "/#.pair_distance_histogram.tab";

		this->reference_bam_file_name = this->data_path + "/reference.bam";
		this->reference_fasta_file_name = this->data_path + "/reference.fasta";
		this->reference_faidx_file_name = this->data_path + "/reference.fasta.fai";
		this->reference_gff3_file_name = this->data_path + "/reference.gff3";
		this->discordant_pairs_file_name = this->data_path + "/#.discordant_pairs.csv";
		this->unmapped_reads_fastq_file_name = this->data_path + "/unmapped_reads.fastq.gz";
    this->data_vcf_file_name = this->data_path + "/output.vcf";
    this->data_genome_diff_file_name = this->data_path + "/output.gd";
    this->data_annotated_genome_diff_file_name = this->data_path + "/annotated.gd";
    this->data_json_summary_file_name = this->data_path + "/summary.json";

    //


    //! Paths: Experimental
    this->long_pairs_file_name = this->output_path + "/long_pairs.tab";

    
    //! Multithreading
    Settings::pool.resize(this->num_processors);

    //! Display name for the HTML output. Last, so that it can use the paths
    //! and read file sets resolved above.
    this->resolve_sample_name();
  }

  void Settings::resolve_sample_name()
  {
    // 1. An explicit --name always wins.
    if (!this->print_custom_run_name.empty()) {
      this->sample_name = this->print_custom_run_name;
      return;
    }

    // 2. The title of the --header-genome-diff. Note that cGenomeDiff::read()
    //    seeds the title from the file's own base name, so this is populated
    //    even for a header GD with no explicit #=TITLE line. Guard on existence
    //    because read() asserts on a file it cannot open -- reporting that is
    //    the job of the check in do_breseq(), not of naming the sample.
    if (!this->header_genome_diff_file_name.empty()
        && file_exists(this->header_genome_diff_file_name.c_str())) {
      cGenomeDiff header_gd;
      header_gd.read(this->header_genome_diff_file_name, true);
      string title = header_gd.get_title();
      if (!title.empty()) {
        // set_title() stores spaces as underscores; undo that for display, the
        // same way print_custom_run_name is derived from custom_run_name.
        this->sample_name = substitute(title, "_", " ");
        return;
      }
    }

    // 3. The name of the output directory.
    string output_dir = this->base_output_path;
    while ((output_dir.size() > 1) && (output_dir[output_dir.size()-1] == '/'))
      output_dir.erase(output_dir.size()-1);
    string output_dir_name = path_to_filename(output_dir);
    if (!output_dir_name.empty() && (output_dir_name != ".") && (output_dir_name != "..")) {
      this->sample_name = output_dir_name;
      return;
    }

    // 4. Fall back to the read files. Uses the per-set names, which collapse a
    //    read pair to a single entry (sample_R1 + sample_R2 => sample_RX).
    vector<string> read_set_names;
    for (cReadFileSets::const_iterator it = this->read_file_sets.begin(); it != this->read_file_sets.end(); it++) {
      read_set_names.push_back(it->m_base_name);
    }
    // With no read files at all (e.g. the default-constructed Settings that
    // gdtools uses) leave the name empty so that nothing is shown, rather than
    // rendering a bare "Unnamed ()".
    if (!read_set_names.empty())
      this->sample_name = "Unnamed (" + join(read_set_names, ", ") + ")";
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
    this->installed["bowtie2"] = which("bowtie2");
    this->installed["bowtie2-build"] = which("bowtie2-build");
    this->installed["bowtie2_version"] = "";
    this->installed["bowtie2_version_string"] = "";
    
    if (this->installed["bowtie2"].size() != 0) {
      string version_string = SYSTEM_CAPTURE(this->installed["bowtie2"] + " --version", true);

      // The version is on the line that names the aligner executable. Do not assume a
      // specific executable name or a fixed offset -- newer bowtie2 releases dispatch to
      // a SIMD-specific binary, so the name varies. Examples of this line:
      //
      //   /usr/local/bin/bowtie2-align-s version 2.5.4
      //   /usr/local/bowtie2-2.5.5/bowtie2-align-s-v256 version 2.5.5   (GitHub issue #422)
      //   /usr/local/bin/bowtie2-align-s version 2.0.0-beta7
      //
      // Other lines may precede it and must be skipped, for example (GitHub issue #412):
      //
      //   [WARNING] Failed to launch x86-64-v3 version, staying with default
      //
      string new_version_string = "";
      vector<string> version_lines = split_on_any(version_string, "\n\r");
      for (vector<string>::iterator it=version_lines.begin(); it!=version_lines.end(); it++) {

        // Only the line naming the bowtie2 executable can tell us the version
        if (it->find("bowtie2") == string::npos) continue;

        // Take the word after the word "version", requiring it to begin with a digit
        vector<string> version_words = split_on_whitespace(*it);
        for (size_t i=0; i+1<version_words.size(); i++) {
          if (version_words[i] != "version") continue;
          if (version_words[i+1].find_first_of("0123456789") != 0) continue;
          new_version_string = version_words[i+1];
          break;
        }
        if (new_version_string.size() > 0) break;
      }

      if (new_version_string.size() > 0) {
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

    // detect gnuplot, used to render all diagnostic plots
    this->installed["gnuplot"] = which("gnuplot");
    this->installed["gnuplot_version_string"] = "";
    this->installed["gnuplot_version"] = "";

    if (this->installed["gnuplot"].size() > 0)
    {
      // version is reported as e.g. "gnuplot 5.4 patchlevel 1"
      string gnuplot_version = SYSTEM_CAPTURE("gnuplot --version", true);

      // default if output does not match our pattern
      this->installed["gnuplot_version"] = "0";

      size_t start_version_pos = gnuplot_version.find("gnuplot ");
      size_t patchlevel_pos = gnuplot_version.find("patchlevel ");
      if ((start_version_pos != string::npos) && (patchlevel_pos != string::npos))
      {
        start_version_pos += 8;
        size_t end_version_pos = gnuplot_version.find_first_of(" \t\r\n", start_version_pos);
        string major_minor = gnuplot_version.substr(start_version_pos, end_version_pos - start_version_pos);

        size_t patch_start_pos = patchlevel_pos + 11;
        size_t patch_end_pos = gnuplot_version.find_first_not_of("0123456789", patch_start_pos);
        if (patch_end_pos == string::npos) patch_end_pos = gnuplot_version.size();
        string patch = gnuplot_version.substr(patch_start_pos, patch_end_pos - patch_start_pos);

        vector<string> split_version_string = split(major_minor, ".");
        if (split_version_string.size() == 2)
        {
          string version_string = major_minor + "." + patch;
          this->installed["gnuplot_version_string"] = version_string;
          uint32_t gnuplot_numerical_version = from_string<uint32_t>(split_version_string[0]) * 1000000 + from_string<uint32_t>(split_version_string[1]) * 1000 + from_string<uint32_t>(patch);
          this->installed["gnuplot_version"] = to_string(gnuplot_numerical_version);
        }
      }
      else {
        this->installed["gnuplot_version_error_message"] = gnuplot_version;
      }
    }

    //cerr << "gnuplot: " << this->installed["gnuplot"] << " gnuplot version: " << this->installed["gnuplot_version_string"] << endl;

    // detect samtools, used for BAM sort/index/faidx and in the bowtie2 pipes
    this->installed["samtools"] = which("samtools");
    this->installed["samtools_version_string"] = "";
    this->installed["samtools_version"] = "";

    if (this->installed["samtools"].size() > 0)
    {
      // first line is reported as e.g. "samtools 1.21"
      string samtools_version = SYSTEM_CAPTURE("samtools --version", true);
      size_t start_version_pos = samtools_version.find("samtools ");
      if (start_version_pos != string::npos)
      {
        start_version_pos += 9;
        size_t end_version_pos = samtools_version.find_first_of(" \t\r\n", start_version_pos);
        string version_string = samtools_version.substr(start_version_pos, end_version_pos - start_version_pos);
        this->installed["samtools_version_string"] = version_string;

        // numeric form (major*1e6 + minor*1e3 + patch), consistent with gnuplot/bowtie2
        vector<string> sv = split(version_string, ".");
        while (sv.size() < 3) sv.push_back("0");
        uint32_t numeric = from_string<uint32_t>(sv[0]) * 1000000
                         + from_string<uint32_t>(sv[1]) * 1000
                         + from_string<uint32_t>(sv[2]);
        this->installed["samtools_version"] = to_string(numeric);
      }
      else {
        this->installed["samtools_version_error_message"] = samtools_version;
      }
    }

	}

	void Settings::check_installed()
	{
    // Developer's Note
    //
    // If you are running things through a debugger (like in XCode), your $PATH may need to be
    // set to include the paths where you have Bowtie2 and R installed within your IDE.
    
		bool good_to_go = true;

    cerr << endl << color_yellow("Checking for required executables") << endl;

    if (this->installed["bowtie2"].size() == 0)
    {
      good_to_go = false;
      cerr << color_red("---> ERROR Required executable \"bowtie2\" or \"bowtie2-build\" not found.") << endl;
      cerr << color_red("---> See http://bowtie-bio.sourceforge.net/bowtie2") << endl;
    }
    else if (this->installed.count("bowtie2_version_error_message")) {
      good_to_go = false;
      cerr << color_red("---> ERROR Could not determine version of installed executable \"bowtie2\" or \"bowtie2-build\".") << endl;
      cerr << color_red("---> For found executable installed at [" + this->installed["bowtie2"] + "]") << endl;
      cerr << color_red("---> Commands \"bowtie2 --version\" or \"bowtie2-build --version\" returned:") << endl;
      cerr << color_red(this->installed["bowtie2_version_error_message"]) << endl;
    }
    // version encoded in triplets of numbers
    else if (from_string<uint64_t>(this->installed["bowtie2_version"]) < 2000000007) {
      good_to_go = false;
      cerr << color_red("---> ERROR Required executable \"bowtie2\" version 2.0.0-beta7 or later not found.") << endl;
      cerr << color_red("---> Your version is " + this->installed["bowtie2_version_string"] + ".") << endl;
      cerr << color_red("---> See http://bowtie-bio.sourceforge.net/bowtie2") << endl;
    }
    else if ((from_string<uint64_t>(this->installed["bowtie2_version"]) == 2000003000)
             ||  (from_string<uint64_t>(this->installed["bowtie2_version"]) == 2000004000)) {
      good_to_go = false;
      cerr << color_red("---> ERROR bowtie2 versions 2.0.3 and 2.0.4 are known to have bugs in") << endl;
      cerr << color_red("---> ERROR SAM output that can cause breseq to crash. Please upgrade.") << endl;
      cerr << color_red("---> Your version is " + this->installed["bowtie2_version_string"] + ".") << endl;
      cerr << color_red("---> See http://bowtie-bio.sourceforge.net/bowtie2") << endl;
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
      cerr << color_red("---> ERROR bowtie2 version 2.3.1 is known to have a major bug in SAM") << endl;
      cerr << color_red("---> ERROR output that will cause breseq to crash. Please upgrade/downgrade.") << endl;
      cerr << color_red("---> Your version is " + this->installed["bowtie2_version_string"] + ".") << endl;
      cerr << color_red("---> See http://bowtie-bio.sourceforge.net/bowtie2") << endl;
    }
    else {
      cerr << color_green("---> bowtie2  :: version " + this->installed["bowtie2_version_string"] + " [" + this->installed["bowtie2"] + "]") << endl;
      // WARN if preferred version is not used
      if (from_string<uint64_t>(this->installed["bowtie2_version"]) != 2005005000) {
        cerr << "---> bowtie2  :: NOTE :: breseq output may vary slightly depending on your bowtie2 version," << endl;
        cerr << "---> bowtie2  :: NOTE :: and occasionally bowtie2 versions may have bugs that cause crashes." << endl;
        cerr << "---> bowtie2  :: NOTE :: bowtie2 version 2.5.5 is recommended with this breseq version." << endl;
      }
    }
		// gnuplot version 5.4.0 required: 'set datafile columnheaders', used by every
    // breseq-generated plot script, was not recognized by gnuplot until 5.4.0
    // (confirmed by testing 4.6.0/5.2.7, which reject it, against 5.4.1+, which accept it).
		if (this->installed["gnuplot"].size() == 0)
		{
			good_to_go = false;
			cerr << color_red("---> ERROR Required executable \"gnuplot\" not found.") << endl;
			cerr << color_red("---> See http://www.gnuplot.info") << endl;
		}
    else if (this->installed.count("gnuplot_version_error_message")) {
      good_to_go = false;
      cerr << color_red("---> ERROR Could not determine version of installed executable \"gnuplot\".") << endl;
      cerr << color_red("---> For found executable installed at [" + this->installed["gnuplot"] + "]") << endl;
      cerr << color_red("---> Command \"gnuplot --version\" returned:") << endl;
      cerr << color_red(this->installed["gnuplot_version_error_message"]) << endl;
    }
    else if (from_string<uint32_t>(this->installed["gnuplot_version"]) < 5004000)
		{
      good_to_go = false;
      cerr << color_red("---> ERROR Required executable \"gnuplot version 5.4.0 or later\" not found.") << endl;
      cerr << color_red("---> Your version is " + this->installed["gnuplot_version_string"] + ".") << endl;
      cerr << color_red("---> See http://www.gnuplot.info") << endl;
    }
    else {
      cerr << color_green("---> gnuplot  :: version " + this->installed["gnuplot_version_string"] + " [" + this->installed["gnuplot"] + "]") << endl;
    }
    
    if (this->installed["samtools"].size() == 0)
    {
      good_to_go = false;
      cerr << color_red("---> ERROR Required executable \"samtools\" not found.") << endl;
      cerr << color_red("---> See http://www.htslib.org") << endl;
    }
    else if (this->installed.count("samtools_version_error_message")) {
      good_to_go = false;
      cerr << color_red("---> ERROR Could not determine version of installed executable \"samtools\".") << endl;
      cerr << color_red("---> For found executable installed at [" + this->installed["samtools"] + "]") << endl;
      cerr << color_red("---> Command \"samtools --version\" returned:") << endl;
      cerr << color_red(this->installed["samtools_version_error_message"]) << endl;
    }
    else {
      cerr << color_green("---> samtools :: version " + this->installed["samtools_version_string"] + " [" + this->installed["samtools"] + "]") << endl;
    }

		if (!good_to_go) exit(0);
	}

	bool Settings::do_step(const string& done_key, const string& message)
	{
		string done_file_name = done_key;
		this->done_key_messages[done_key] = message;
    this->set_current_step_done_key(done_key);
		if (!file_exists(done_file_name.c_str()))
		{
      end_progress_line();
      cerr << endl << color_yellow("+++   NOW PROCESSING " + message) << endl;
      this->record_start_time(message);
      return true;
		}

		end_progress_line();
		cerr << endl << color_yellow("--- ALREADY COMPLETE " + message) << endl;

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

