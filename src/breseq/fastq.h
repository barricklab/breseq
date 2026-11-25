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

#ifndef _BRESEQ_FASTQ_H_
#define _BRESEQ_FASTQ_H_

#include "common.h"

#include "settings.h"
#include "summary.h"
#include "portable_random.h"

namespace breseq {

class cAnnotatedSequence;

	/*! Interface for loading and manipulating files in FASTQ format.
   */
  
  /*! estimate_fastq_file_bases

      Estimates the number of bases in a FASTQ file (gzipped or not) without reading all of it.
      Reads from the start of the file until FASTQ_ESTIMATE_COMPRESSED_BYTES of the file as stored
      on disk have been consumed, then scales the bases seen by the fraction of the file that was.
      A file smaller than that is read to its end, and its count is exact.

      Used by --limit-fold-coverage to turn the requested coverage into the fraction of reads to
      keep, so that the reads analyzed are drawn from the whole file and not from its start.
   */
  uint64_t estimate_fastq_file_bases(const string &file_name);

  /*! Chooses which reads to keep for --limit-fold-coverage.

      Reads (or read pairs) are kept at even intervals through the input: with a keep fraction of
      0.25, every fourth one. The subset is therefore spread through the whole input however it is
      ordered. Reads deposited from an aligned BAM arrive sorted by reference position, and taking
      them from the start of the file instead covers only the start of the genome.

      Nothing here is random, so the same input and the same limit always give the same reads. The
      fraction is held as a 32-bit fixed-point step added to an accumulator, and a read is kept each
      time the accumulator passes one. That is integer arithmetic only, so the choice is also the
      same on every platform, which the tests that compare a committed expected.gd depend on.
   */
  class cReadSubsampler {
    static const uint64_t ONE = static_cast<uint64_t>(1) << 32;
    uint64_t m_step;          // keep fraction, in units of 1/2^32; ONE keeps every read
    uint64_t m_accumulator;   // starts at one half, so the kept reads are centered in their intervals
  public:
    cReadSubsampler(const double keep_fraction)
    : m_step( (keep_fraction >= 1.0) ? ONE
              : (keep_fraction <= 0.0) ? 0
              : static_cast<uint64_t>(keep_fraction * 4294967296.0) )
    , m_accumulator(ONE / 2)
    { }

    bool keep() {
      m_accumulator += m_step;
      if (m_accumulator < ONE) return false;
      m_accumulator -= ONE;
      return true;
    }
  };

  /*! normalize_fastq
   
      Main function for this analysis. Prints summary information
      about fastq and converts file to SANGER format if necessary.
   */  
  AnalyzeFastqSummary normalize_fastq(
                                        const string &file_name,
                                        const string &convert_file_name,
                                        const uint32_t file_index,
                                        const int32_t trim_end_on_base_quality,
                                        const bool filter_reads,
                                        const double read_keep_fraction,
                                        const uint32_t read_length_min,
                                        const double max_same_base_fraction,
                                        const double max_N_fraction,
                                        const uint32_t _long_read_trigger_length,
                                        const uint32_t _long_read_split_length,
                                        const bool _long_read_distribute_remainder,
                                        // Quality score given to every base of a file whose bases
                                        // all have one quality score (0 = keep them as they are).
                                        // See kDefaultConstantBaseQualityReplacement below.
                                        const uint32_t constant_base_quality_replacement,
                                        // Cut the reads of a mate-pair library at its junction
                                        // adapter (see MATE-PAIR JUNCTION ADAPTERS below). The
                                        // library is detected and reported either way.
                                        const bool trim_mate_pair_junction_adapters,
                                        const uint32_t num_threads = 1,
                                        // EXPERIMENTAL synthetic pairs from long reads. When
                                        // long_read_pair_distance is nonzero, a split piece is
                                        // paired with the piece that many bases downstream on the
                                        // same long read. Mate 1 goes to lp1_convert_file_name with
                                        // read-name prefix file_index+1 and mate 2 to
                                        // lp2_convert_file_name with prefix file_index+2. Both keep
                                        // the long read's strand, so the pairs are FF with mate 1
                                        // upstream -- the arrangement the pieces really have. Pieces
                                        // without a partner stay in convert_file_name as before.
                                        const uint32_t long_read_pair_distance = 0,
                                        const string &lp1_convert_file_name = "",
                                        const string &lp2_convert_file_name = "",
                                        AnalyzeFastqSummary* lp1_summary = NULL,
                                        AnalyzeFastqSummary* lp2_summary = NULL
                                        );

  pair<AnalyzeFastqSummary, AnalyzeFastqSummary> normalize_fastq_paired(
                                        const string &r1_file_name,
                                        const string &r1_convert_file_name,
                                        const string &r2_file_name,
                                        const string &r2_convert_file_name,
                                        const uint32_t r1_file_index,
                                        const uint32_t r2_file_index,
                                        const int32_t trim_end_on_base_quality,
                                        const bool filter_reads,
                                        const double read_keep_fraction,
                                        const uint32_t read_length_min,
                                        const double max_same_base_fraction,
                                        const double max_N_fraction,
                                        const uint32_t long_read_trigger_length,
                                        const uint32_t long_read_split_length,
                                        const bool long_read_distribute_remainder,
                                        const uint32_t constant_base_quality_replacement,
                                        const bool trim_mate_pair_junction_adapters,
                                        const uint32_t num_threads = 1
                                        );

  //! Read files whose base quality scores carry no information
  //
  // Some read files give every base the same quality score. PacBio Sequel CLR reads downloaded from
  // the SRA, for example, are all '!' (Q0). One value cannot tell good bases from bad ones, and Q0 is
  // at or below --base-quality-cutoff, so as they stand every base would be ignored and the reads
  // would yield no base-level evidence. When a file has only one quality score, breseq instead gives
  // all of its bases one score above the cutoff, --constant-base-quality-replacement (default 20).
  // The exact value matters little: error calibration (stage 07) measures the error rate of each
  // quality score from the alignments, so the one score these bases share gets an error rate
  // measured from the data like any other.
  //
  // Detected and replaced where qualities enter breseq: normalize_fastq/normalize_fastq_paired for
  // FASTQ input, and, for --aligned-sam input, sam_file_constant_base_quality() during sequence
  // conversion with the replacement made as the SAM records are loaded (load_sam_only_alignments).
  const uint32_t kDefaultConstantBaseQualityReplacement = 20;
  // Highest score a SANGER quality character can hold ('~')
  const uint32_t kMaxSangerBaseQuality = 93;
  // htslib stores a SAM QUAL of '*' (qualities absent) as this value for every base
  const uint32_t kMissingBamBaseQuality = 0xff;

  //! Prints the warning shown for each read file whose qualities are replaced.
  //  original_quality is on the SANGER scale, or kMissingBamBaseQuality for a SAM QUAL of '*'.
  void warn_constant_base_quality(const string& file_name, const uint32_t original_quality, const uint32_t replacement_quality);

  // READ TYPES
  //
  // Every read file has a read type: the sequencing platform and, for long reads, the kind of
  // basecalling. It is either given on the command line as a prefix on the read file argument
  // (ont-hac:reads.fastq.gz) or detected from the first reads of the file (cReadTypeDetector).
  // The names are the ones brefito uses for the "type" column of its data.csv, synonyms included,
  // so a brefito type can be handed to breseq as it is.
  //
  // The type is RECORDED ONLY at present (stage 1 log, summary.json, summary.html, @RG PL). No
  // calling behavior depends on it yet; it is there so that per-type artifact handling has
  // something to key on. Consumers should ask cReadFile::m_read_type.
  //
  // POLICY: no setting is to be chosen from a nanopore accuracy (fast versus hac) whose source is
  // "quality" alone. That guess is kept because it is usually right and worth showing, but it
  // fails by the dataset (see kReadTypeNanoporeHACMeanQuality below). Anything that changes
  // behavior for fast-basecalled reads must have a model name from the read headers, a type
  // given on the command line, or the indel error rate measured after mapping. Check
  // cReadFile::m_read_type_source before acting on m_read_type.
  extern const char* const kReadTypeIllumina;      // "illumina"
  extern const char* const kReadTypeNanopore;      // "nanopore"       basecalling accuracy not known
  extern const char* const kReadTypeNanoporeFast;  // "nanopore-fast"
  extern const char* const kReadTypeNanoporeHAC;   // "nanopore-hac"   high or super accuracy: sup is not
                                                   //                  a type of its own, "nanopore-sup"
                                                   //                  is accepted as a synonym for this
  extern const char* const kReadTypePacBio;        // "pacbio"         CLR or HiFi not known
  extern const char* const kReadTypePacBioCLR;     // "pacbio-clr"
  extern const char* const kReadTypePacBioHiFi;    // "pacbio-hifi"
  extern const char* const kReadTypeUnknown;       // "unknown"        only for a file with no reads

  // How a file's type was decided
  extern const char* const kReadTypeSourceSpecified; // "specified"  prefix on the command line
  extern const char* const kReadTypeSourceHeader;    // "header"     read names/comments identify it
  extern const char* const kReadTypeSourceQuality;   // "quality"    guessed from lengths and qualities

  // Thresholds of the guess from lengths and qualities (cReadTypeDetector::detect)
  const uint32_t kReadTypeSampleReads = 10000;           // reads examined per file
  const uint32_t kReadTypeLongReadLength = 700;          // 99th percentile read length of a long-read file...
  const uint32_t kReadTypeLongestShortRead = 1000;       // ...or any read this long: no short-read platform makes one
  const uint32_t kReadTypeHiFiQuality = 60;              // a base quality only consensus reads reach...
  const double   kReadTypeHiFiQualityFraction = 0.1;     // ...and the fraction of bases that must have it
  const uint32_t kReadTypeBinnedQualityValues = 10;      // at most this many distinct scores: binned
  const uint32_t kReadTypeCLRMaximumQuality = 20;        // PacBio RS II subreads have no base this good...
  const double   kReadTypeCLRMaximumQualityFraction = 0.001; // ...(allowing this fraction), nanopore reads do
  // Nanopore accuracy from the mean quality (the quality of the mean error probability).
  //
  // Calibrated on samples of public runs whose basecalling mode is on record: ENA study ERP173262
  // (bacterial isolates, R9.4.1 and R10.4.1, each called fast/hac/sup with Guppy and with Dorado;
  // 16 fast, 16 hac and 16 sup samples), the three basecallers' calls of one R9.4.1 library in
  // SRP416531, four other runs described as fast-basecalled, and dorado v5 hac/sup calls of three
  // R10.4.1 flowcells (ONT open data giab_2025.01).
  //
  // Within ERP173262 fast is Q9.3-12.0 and hac Q11.9-13.7 (sup Q13.1-16.7): the classes touch, and
  // 12.2 misclassifies one hac sample of the 48. No other statistic tried did better (arithmetic
  // mean, median, mode, fraction of bases at Q20/Q30/Q40, median per-read quality).
  //
  // ACROSS datasets the mean quality reflects the library as much as the basecaller, and the
  // guess fails by the dataset: every call of the SRP416531 library, sup included, is Q5-7 and so
  // "fast"; one run described as Dorado fast is Q16 and so "hac". Over all 69 labelled samples
  // 21 of 22 fast and 40 of 47 hac/sup are called correctly. Treat a type whose source is
  // "quality" accordingly: a header that names the model, or a type given on the command line,
  // is always better. (Base quality inside the reads' own homopolymer runs was tried too: fast
  // calls score those bases relatively higher, but the ranges overlap as much.)
  //
  // What does separate the modes is the indel error rate AFTER mapping. For 14 E. coli samples of
  // ERP173262 mapped with minimap2, indels per 100 aligned bases inside reference homopolymers of
  // 4+ bases were 6.2-7.0 for fast against 1.7-4.8 for hac/sup, and deletions overall 1.8-2.8
  // against 0.6-1.5. A rule that needs the fast/hac distinction should measure that (stage 07
  // tabulates errors per read file) rather than rely on this guess.
  const double   kReadTypeNanoporeHACMeanQuality = 12.2; // mean quality below this: fast

  //! Canonical read type for a name or synonym (case-insensitive); "" if it is not one.
  string canonical_read_type(const string& name);
  //! Canonical read types, for help and error messages.
  vector<string> read_type_names();
  //! SAM @RG PL value for a read type ("ILLUMINA", "ONT", "PACBIO"); "" for unknown.
  string read_type_platform(const string& read_type);
  bool is_long_read_type(const string& read_type);
  bool is_nanopore_read_type(const string& read_type);
  //! Text for the stage 1 log and summary.html: "detected from read names", ...
  string read_type_source_description(const string& read_type_source);

  class cReadTypeDetector;

  // MATE-PAIR JUNCTION ADAPTERS
  //
  // A mate-pair library circularizes each fragment of several kb and sequences across the
  // junction, so a read is genome, then the junction adapter, then genome from the OTHER END of
  // the fragment, kilobases away. Two protocols leave a known adapter there:
  //
  //   Illumina Nextera mate pair   CTGTCTCTTATACACATCT|AGATGTGTATAAGAGACAG   (two Tn5 mosaic ends
  //                                back to back; its own reverse complement)
  //   Cre-lox (454 and Illumina    TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG   (a loxP site)
  //   long mate pair)
  //
  // Such libraries are deposited as ordinary paired runs, and left alone they wreck the junction
  // steps. Every read through the adapter is a chimera, breseq built candidate junctions whose
  // unique sequence was the adapter, and each adapter read then seeded on hundreds of them under
  // bowtie2 -k 2000: 37 minutes per file against 15 seconds for an ordinary library of the same
  // genome (Nextera, DRR315860 against DRR315859), and 52 predicted junctions that were every one
  // the Cre-lox adapter (SRR2084429, at 24-fold coverage).
  //
  // breseq protects itself in two places:
  //
  //  1. Sequence conversion counts the reads of the first sample that contain an exact adapter.
  //     A file with at least kMatePairJunctionAdapterFraction of them is a mate-pair library
  //     (26-51% of reads in five such runs; an ordinary library sequenced on the same flowcell
  //     as one had 0.1%, from index hopping, which is why the gate is not zero). Its reads are
  //     cut at the adapter, keeping the part before it: that piece keeps the mate-pair relation
  //     to the other read of the pair, where the part after it belongs to the other end of the
  //     fragment. breseq warns and reports it, since a tool made for mate-pair libraries
  //     recovers more of the data. Files below the gate are not touched at all, so nothing
  //     changes for any other library. --no-mate-pair-adapter-trimming keeps the reads as they
  //     are; the library is still reported.
  //
  //  2. Candidate junctions whose unique read sequence contains adapter are never made
  //     (CandidateJunctions::alignment_pair_to_candidate_junction), for the reads the trimming
  //     misses -- adapters with indels or several errors, about 0.6% of reads -- and for any
  //     library.
  //
  // A mate-pair library with no adapter at its junctions (the original Illumina protocol) cannot
  // be handled this way and is not detected.
  struct cMatePairJunctionAdapter {
    string name;       // for messages: "Nextera mate-pair junction adapter"
    string sequence;   // as given above; reads may contain it or its reverse complement
  };
  //! The known junction adapters.
  const vector<cMatePairJunctionAdapter>& mate_pair_junction_adapters();

  const double   kMatePairJunctionAdapterFraction = 0.02;    // of sampled reads with an exact adapter
  const uint32_t kJunctionAdapterMaxMismatches = 2;          // allowed in either half of the adapter
  const uint32_t kJunctionAdapterMinEndMatch = 10;           // bases of adapter a read may end in
  const uint32_t kJunctionAdapterKmer = 12;                  // adapter k-mer that disqualifies a junction

  //! Index (in mate_pair_junction_adapters()) of an adapter the sequence contains exactly, on
  //  either strand; -1 if none.
  int32_t exact_mate_pair_junction_adapter(const string& sequence);
  //! Where to cut a read to remove a junction adapter and everything after it; string::npos if
  //  it has none. Finds either half of an adapter, on either strand, with up to
  //  kJunctionAdapterMaxMismatches mismatches, and a read that ends partway into one (at least
  //  kJunctionAdapterMinEndMatch bases of it).
  size_t mate_pair_junction_adapter_position(const string& sequence);
  //! True if the sequence contains kJunctionAdapterKmer consecutive bases of a junction adapter.
  bool contains_mate_pair_junction_adapter_kmer(const string& sequence);

  //! Counts the junction adapters in the first reads of a file, to decide whether it is a
  //  mate-pair library.
  class cMatePairJunctionAdapterCounter {
  public:
    cMatePairJunctionAdapterCounter() : m_num_reads(0), m_counts(mate_pair_junction_adapters().size(), 0) {}
    void add(const string& sequence)
    {
      m_num_reads++;
      int32_t index = exact_mate_pair_junction_adapter(sequence);
      if (index >= 0) m_counts[index]++;
    }
    //! Fraction of the reads with the most common adapter, and that adapter's name ("" if none)
    double fraction(string& adapter_name) const
    {
      adapter_name = "";
      size_t best = 0;
      for (size_t i = 1; i < m_counts.size(); i++) if (m_counts[i] > m_counts[best]) best = i;
      if ((m_num_reads == 0) || m_counts.empty() || (m_counts[best] == 0)) return 0.0;
      adapter_name = mate_pair_junction_adapters()[best].name;
      return static_cast<double>(m_counts[best]) / static_cast<double>(m_num_reads);
    }
  private:
    uint64_t m_num_reads;
    vector<uint64_t> m_counts;
  };

  //! Feeds the first reads of a FASTQ file to a detector (breseq SUMMARIZE-FASTQ; the pipeline
  //  does this inside normalize_fastq, from the reads it has already buffered).
  void sample_fastq_file_for_read_type(const string& file_name, cReadTypeDetector& detector);

  //! Settles a file's read type and prints the stage 1 log line for it.
  //  read_type/read_type_source come in as detected and go out as decided: a type given on the
  //  command line (specified_read_type, "" if none) wins, with a warning when the reads say otherwise.
  void resolve_read_type(const string& file_name, const string& specified_read_type,
                         string& read_type, string& read_type_source);

  //! Splits a read file argument of the form <read type>:<path>.
  //  The text before the first ':' is a prefix only if it is a read type name, so a path that
  //  merely contains a colon is left alone. Returns the path; read_type is "" without a prefix.
  string split_read_type_prefix(const string& read_file_argument, string& read_type);
  //! Same, for callers that only want the path (AnyOption's positional path filter).
  string strip_read_type_prefix(const string& read_file_argument);

  /*! Guesses the read type of a file from its first reads.

      Add each sampled read's name, the rest of its header line, its base quality characters as
      they are in the file, and its length BEFORE any splitting; then ask for the type.

      Read headers are trusted over everything else when most reads have one that identifies the
      platform: a dorado/guppy basecaller model name, a PacBio movie/ZMW name, an Illumina
      instrument:run:flowcell:lane:tile:x:y name. SRA downloads replace the name and only sometimes
      keep the original after it, so every whitespace-separated field of the header is examined.

      Otherwise the type is guessed from read lengths and base qualities, and reported with the
      source "quality". The nanopore fast/hac/sup split by mean quality in particular is only a
      best guess: the distributions overlap across chemistries and basecaller versions.
   */
  class cReadTypeDetector {
  public:
    cReadTypeDetector();

    //! quality_chr_offset is 33 for SANGER-encoded characters, 64 for the old Illumina formats.
    //  It may be given later, to detect(), when it is not yet known as reads are added.
    void add(const string& name, const string& comment, const string& qualities, const size_t length);

    void detect(const uint8_t quality_chr_offset, string& read_type, string& read_type_source) const;

    // Statistics of the sample, for reporting and for calibrating the rules (SUMMARIZE-FASTQ)
    uint64_t num_reads() const { return m_num_reads; }
    uint32_t read_length_percentile(const double fraction) const;
    //! Quality of the mean error probability over all sampled bases
    double mean_quality(const uint8_t quality_chr_offset) const;
    //! Fraction of sampled bases with at least this quality
    double fraction_bases_at_least_quality(const uint8_t quality_chr_offset, const uint32_t quality) const;
    uint32_t num_distinct_qualities() const;

  private:
    enum HeaderEvidence { kNone = 0, kIllumina, kNanopore, kNanoporeFast, kNanoporeHAC,
                          kPacBioCLR, kPacBioHiFi, kNumHeaderEvidence };
    static HeaderEvidence classify_header(const string& name, const string& comment);

    uint64_t m_num_reads;
    uint64_t m_header_votes[kNumHeaderEvidence];
    uint64_t m_quality_chr_counts[256];
    mutable vector<uint32_t> m_read_lengths;   // sorted on demand
  };

  // Utility function for converting FASTQ files between formats
  void convert_fastq(const string &from_file_name, const string &to_file_name, const string &from_format, const string &to_format, bool reverse_complement = false);
  
  /*! Sequence class.
   */
   
  struct cFastqSequence {
    public:
      string   m_name;        //@NAME
      string   m_name_comment;//rest of the @NAME line after the first space (only read, never written)
      string   m_sequence;    //sequence
      string   m_name_plus;   //+NAME
      string   m_qualities;   //quality score characters
      bool     m_numerical_qualities; // quality scores read were read from numerical format
      uint32_t m_base_counts[base_list_including_N_size]; // number of each base, including N. Used for filtering.
    
    
    cFastqSequence() 
    : m_numerical_qualities(false)
    {
      for (uint8_t b=0; b<base_list_including_N_size; b++) {
        m_base_counts[b] = 0;
      }
    }
    
    bool identical(cFastqSequence& seq);
    
    size_t length() { return m_sequence.length(); }
   };

  void fastq_sequence_trim_end_on_base_quality(cFastqSequence& seq, const uint32_t base_quality);


  //! Per-cycle base quality model for simulated reads.
  //
  //  Deliberately parametric rather than an empirical table: no data files, no new dependencies, and
  //  integer-only so the output is byte-identical on every platform (see portable_random.h).
  //
  //  The mean quality declines QUADRATICALLY with the cycle number, which is the shape real Illumina
  //  data has -- essentially flat over the first half of the read, falling off fastest over the last
  //  few cycles. This replaces a hard-coded cumulative table taken from 2008 GAII runs that had two
  //  problems: it was position-independent (the same distribution at cycle 1 and cycle 100, which is
  //  the least realistic thing a quality model can be), and its marginal implied a ~3.3% per-base
  //  substitution rate, about 3 mismatches in every 100 bp read.
  //
  //  Constant quality needs no special case: set phred_start == phred_end and phred_stdev = 0.
  struct cSimQualityModel {
    uint32_t phred_start;   //!< mean Phred at the first cycle
    uint32_t phred_end;     //!< mean Phred at the last cycle
    uint32_t phred_stdev;   //!< per-base scatter around the cycle mean
    uint32_t phred_min;     //!< floor
    uint32_t phred_max;     //!< ceiling

    cSimQualityModel()
      : phred_start(38), phred_end(28), phred_stdev(4), phred_min(2), phred_max(41) {}

    //! Draw one quality score as a Sanger-encoded (offset 33) character. cycle is 0-based.
    char draw(cPortableRandom& rng, uint32_t cycle, uint32_t read_size) const;
  };

  //! Indel error rates for simulated reads, per million bases. The substitution rate is NOT here:
  //  it derives from each base's own drawn quality (see cSimFastqSequence::is_error_base), which is
  //  what makes the quality string meaningful to breseq's per-quality-bin error calibration.
  struct cSimErrorModel {
    uint32_t insertion_rate_per_million;
    uint32_t deletion_rate_per_million;

    cSimErrorModel() : insertion_rate_per_million(10), deletion_rate_per_million(10) {}  // 1e-5 each
  };

  class cSimFastqSequence {
    public:
      static int32_t SEED_VALUE;

      static map<char, string> random_snp_base_options;

      static char random_insertion_base_options[];

      static char get_random_error_base(cPortableRandom& rng, const char not_this_base);
      static char get_random_insertion_base(cPortableRandom& rng);

      //! P(substitution) = 10^(-phred/10), evaluated from an integer table rather than pow().
      static bool is_error_base(cPortableRandom& rng, uint32_t phred);

      //! Build one read whose footprint in the reference is [left_1, left_1 + read_size - 1].
      //
      //  strand == +1 returns that sequence; strand == -1 returns its reverse complement, so the
      //  read's 5' end sits at left_1 + read_size - 1. Stating the footprint directly (rather than a
      //  strand-dependent buffer origin, as the previous interface did) removes a genuinely
      //  confusing asymmetry: with the old contract a minus-strand read's footprint was actually
      //  [start_1 + read_size, start_1 + 2*read_size - 1], not [start_1, start_1 + read_size - 1].
      //
      //  Simulated deletions consume extra reference bases, so a tail beyond the footprint is
      //  fetched; near the end of a LINEAR sequence that tail is short, and the builder stops
      //  generating deletions rather than reading past it.
      static cFastqSequence simulate_read(const cAnnotatedSequence& ref_sequence,
                                          int32_t left_1,
                                          uint32_t read_size,
                                          int8_t strand,
                                          const cSimQualityModel& quality_model,
                                          const cSimErrorModel& error_model,
                                          cPortableRandom& rng,
                                          const string& read_name,
                                          char& min_quality_char,
                                          bool verbose = false);

      static void simulate_single_ends(const cAnnotatedSequence& sequence,
                                       uint32_t n_reads,
                                       uint32_t read_size,
                                       uint32_t seed,
                                       const cSimQualityModel& quality_model,
                                       const cSimErrorModel& error_model,
                                       string file_name,
                                       bool verbose = false);

      static void simulate_paired_ends(const cAnnotatedSequence& sequence,
                                       uint32_t n_pairs,
                                       uint32_t read_size,
                                       uint32_t fragment_mean,
                                       uint32_t fragment_stdev,
                                       uint32_t seed,
                                       const cSimQualityModel& quality_model,
                                       const cSimErrorModel& error_model,
                                       string pair_1_file_name,
                                       string pair_2_file_name,
                                       bool verbose = false);

      static void simulate_tiled(const cAnnotatedSequence& sequence,
                                 uint32_t read_size,
                                 uint32_t coverage,
                                 const cSimQualityModel& quality_model,
                                 string file_name,
                                 bool verbose);

      //! breseq infers a FASTQ's quality encoding from the LOWEST quality character in the whole
      //  file (see normalize_fastq): a minimum of 59 or more is read as SOLEXA and every score is
      //  then shifted by 31. A simulated library with no low-quality tail therefore gets silently
      //  reinterpreted. Callers pass the minimum character they actually wrote; this warns if it is
      //  high enough to be misread.
      static void warn_if_quality_format_ambiguous(char min_quality_char);
  };

  
  /*! Quality score conversion class.
   */
  
  struct cFastqQualityConverter : public vector<uint8_t> {
  public:    
    cFastqQualityConverter(const string &_from_quality_format, const string &_to_quality_format);
    ~cFastqQualityConverter() {};
    
    string from_quality_format;
    string to_quality_format;
    string from_quality_type;
    string to_quality_type;
    int32_t from_chr_offset;
    int32_t to_chr_offset;
    
    void convert_sequence(cFastqSequence &seq);
    
    static string predict_fastq_file_format(const string& file_name, uint64_t& num_original_reads, uint64_t& num_original_bases, uint32_t& read_length_min, uint32_t& read_length_max, uint8_t& min_quality_score, uint8_t& max_quality_score);

  };
   

	/*! File class.
	 */ 
  
  class cFastqFile : public flexgzfstream {
    
  protected:
    uint32_t  m_current_line;
    string    m_file_name;
    
  public:
    // keep track of duplicate read names one after another and append r# to later ones
    bool m_check_for_repeated_read_names;
    string m_last_read_name;
    uint32_t m_repeated_read_name_count;
    
  
    cFastqFile();
    cFastqFile(const string &file_name, ios_base::openmode mode, unsigned int num_threads = 1);
    ~cFastqFile() { };
      
    bool read_sequence(cFastqSequence &sequence, cFastqQualityConverter& fqc);
    void write_sequence(const cFastqSequence &sequence);
  };

  /*! General sequence helper function.
	 */ 
  
  extern char reverse_complement_lookup_table[256];
  
  inline string reverse_complement(const string& seq)
	{
		string retval(seq.length(), ' ');
		for (uint32_t i = 0; i < seq.size(); i++)      
      retval[i] = reverse_complement_lookup_table[static_cast<uint8_t>(seq[seq.size() - 1 - i])];
    return retval;
	}
  
  inline char reverse_complement(const char& seq)
	{
		char retval(' ');
    retval = reverse_complement_lookup_table[static_cast<uint8_t>(seq)];
    return retval;
	}
  
  inline cFastqSequence reverse_complement(const cFastqSequence& _seq) {
    cFastqSequence seq = _seq;
    seq.m_sequence = reverse_complement(seq.m_sequence);
    seq.m_qualities = reverse_string(seq.m_qualities);
    seq.m_name += "_RC";
    return seq;
  }
	
} // breseq namespace

#endif
