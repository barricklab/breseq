package Bio::DB::Sam;
# $Id: Sam.pm 23215 2010-05-14 21:10:10Z lstein $

our $VERSION = '1.19';

=head1 NAME

Bio::DB::Sam -- Read SAM/BAM database files

=head1 SYNOPSIS

 use Bio::DB::Sam;

 # high level API
 my $sam = Bio::DB::Sam->new(-bam  =>"data/ex1.bam",
                             -fasta=>"data/ex1.fa",
			     );

 my @targets    = $sam->seq_ids;
 my @alignments = $sam->get_features_by_location(-seq_id => 'seq2',
                                                 -start  => 500,
                                                 -end    => 800);
 for my $a (@alignments) {

    # where does the alignment start in the reference sequence
    my $seqid  = $a->seq_id;
    my $start  = $a->start;
    my $end    = $a->end;
    my $strand = $a->strand;
    my $cigar  = $a->cigar_str;
    my $paired = $a->get_tag_values('PAIRED');

    # where does the alignment start in the query sequence
    my $query_start = $a->query->start;     
    my $query_end   = $a->query->end;

    my $ref_dna   = $a->dna;        # reference sequence bases
    my $query_dna = $a->query->dna; # query sequence bases

    my @scores    = $a->qscore;     # per-base quality scores
    my $match_qual= $a->qual;       # quality of the match
 }

 my @pairs = $sam->get_features_by_location(-type   => 'read_pair',
                                            -seq_id => 'seq2',
                                            -start  => 500,
                                            -end    => 800);

 for my $pair (@pairs) {
    my $length                    = $pair->length;   # insert length
    my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
    my $f_start = $first_mate->start;
    my $s_start = $second_mate->start;
 }

 # low level API
 my $bam          = Bio::DB::Bam->open('/path/to/bamfile');
 my $header       = $bam->header;
 my $target_count = $header->n_targets;
 my $target_names = $header->target_name;
 while (my $align = $bam->read1) {
    my $seqid     = $target_names->[$align->tid];
    my $start     = $align->pos+1;
    my $end       = $align->calend;
    my $cigar     = $align->cigar_str;
 }

 my $index = Bio::DB::Bam->index_open('/path/to/bamfile');
 my $index = Bio::DB::Bam->index_open_in_safewd('/path/to/bamfile');

 my $callback = sub {
     my $alignment = shift;
     my $start       = $alignment->start;
     my $end         = $alignment->end;
     my $seqid       = $target_names->[$alignment->tid];
     print $alignment->qname," aligns to $seqid:$start..$end\n";
 }
 my $header = $index->header;
 $index->fetch($bam,$header->parse_region('seq2'),$callback);

=head1 DESCRIPTION

This module provides a Perl interface to the libbam library for
indexed and unindexed SAM/BAM sequence alignment databases. It
provides support for retrieving information on individual alignments,
read pairs, and alignment coverage information across large
regions. It also provides callback functionality for calling SNPs and
performing other base-by-base functions. Most operations are
compatible with the BioPerl Bio::SeqFeatureI interface, allowing BAM
files to be used as a backend to the GBrowse genome browser
application (gmod.sourceforge.net).

=head2 The high-level API

The high-level API provides a BioPerl-compatible interface to indexed
BAM files. The BAM database is treated as a collection of
Bio::SeqFeatureI features, and can be searched for features by name,
location, type and combinations of feature tags such as whether the
alignment is part of a mate-pair.

When opening a BAM database using the high-level API, you provide the
pathnames of two files: the FASTA file that contains the reference
genome sequence, and the BAM file that contains the query sequences
and their alignments. If either of the two files needs to be indexed,
the indexing will happen automatically. You can then query the
database for alignment features by combinations of name, position,
type, and feature tag.

The high-level API provides access to up to four feature "types":

 * "match": The "raw" unpaired alignment between a read and the
   reference sequence.

 * "read_pair": Paired alignments; a single composite
   feature that contains two subfeatures for the alignments of each 
   of the mates in a mate pair.

 * "coverage": A feature that spans a region of interest that contains
   numeric information on the coverage of reads across the region.

 * "region": A way of retrieving information about the reference
   sequence. Searching for features of type "region" will return a
   list of chromosomes or contigs in the reference sequence, rather
   than read alignments.

 * "chromosome": A synonym for "region".

B<Features> can be en masse in a single call, retrieved in a
memory-efficient streaming basis using an iterator, or interrogated
using a filehandle that return a series of TAM-format lines.

B<SAM alignment flags> can be retrieved using BioPerl's feature "tag"
mechanism. For example, to interrogate the FIRST_MATE flag, one
fetches the "FIRST_MATE" tag:

  warn "aye aye captain!" if $alignment->get_tag_values('FIRST_MATE');

The Bio::SeqFeatureI interface has been extended to retrieve all flags
as a compact human-readable string, and to return the CIGAR alignment
in a variety of formats.  

B<Split alignments>, such as reads that cover introns, are dealt with
in one of two ways. The default is to leave split alignments alone:
they can be detected by one or more "N" operations in the CIGAR
string. Optionally, you can choose to have the API split these
alignments across two or more subfeatures; the CIGAR strings of these
split alignments will be adjusted accordingly.

B<Interface to the pileup routines> The API provides you with access
to the samtools "pileup" API. This gives you the ability to write a
callback that will be invoked on every column of the alignment for the
purpose of calculating coverage, quality score metrics, or SNP
calling.

B<Access to the reference sequence> When you create the Bio::DB::Sam
object, you can pass the path to a FASTA file containing the reference
sequence. Alternatively, you may pass an object that knows how to
retrieve DNA sequences across a range via the seq() of fetch_seq()
methods, as described under new().

If the SAM/BAM file has MD tags, then these tags will be used to
reconstruct the reference sequence when necessary, in which case you
can completely omit the -fasta argument. Note that not all SAM/BAM
files have MD tags, and those that do may not use them correctly due
to the newness of this part of the SAM spec. You may wish to populate
these tags using samtools' "calmd" command.

If the -fasta argument is omitted and no MD tags are present, then the
reference sequence will be returned as 'N'.

The B<main object classes> that you will be dealing with in the
high-level API are as follows:

 * Bio::DB::Sam               -- A collection of alignments and reference sequences.
 * Bio::DB::Bam::Alignment    -- The alignment between a query and the reference.
 * Bio::DB::Bam::Query        -- An object corresponding to the query sequence in
                                  which both (+) and (-) strand alignments are
                                  shown in the reference (+) strand.
 * Bio::DB::Bam::Target       -- An interface to the query sequence in which
                                   (-) strand alignments are shown in reverse
                                   complement

You may encounter other classes as well. These include:

 * Bio::DB::Sam::Segment       -- This corresponds to a region on the reference
                                  sequence.
 * Bio::DB::Sam::Constants     -- This defines CIGAR symbol constants and flags.
 * Bio::DB::Bam::AlignWrapper  -- An alignment helper object that adds split
                                  alignment functionality. See Bio::DB::Bam::Alignment
                                  for the documentation on using it.
 * Bio::DB::Bam::ReadIterator  -- An iterator that mediates the one-feature-at-a-time 
                                  retrieval mechanism.
 * Bio::DB::Bam::FetchIterator -- Another iterator for feature-at-a-time retrieval.

=head2 The low-level API

The low-level API closely mirrors that of the libbam library. It
provides the ability to open TAM and BAM files, read and write to
them, build indexes, and perform searches across them. There is less
overhead to using the API because there is very little Perl memory
management, but the functions are less convenient to use. Some
operations, such as writing BAM files, are only available through the
low-level API.

The classes you will be interacting with in the low-level API are as
follows:

 * Bio::DB::Tam            -- Methods that read and write TAM (text SAM) files.
 * Bio::DB::Bam            -- Methods that read and write BAM (binary SAM) files.
 * Bio::DB::Bam::Header    -- Methods for manipulating the BAM file header.
 * Bio::DB::Bam::Index     -- Methods for retrieving data from indexed BAM files.
 * Bio::DB::Bam::Alignment -- Methods for manipulating alignment data.
 * Bio::DB::Bam::Pileup    -- Methods for manipulating the pileup data structure.
 * Bio::DB::Sam::Fai       -- Methods for creating and reading from indexed Fasta
                              files.
=head1 METHODS

We cover the high-level API first. The high-level API code can be
found in the files Bio/DB/Sam.pm, Bio/DB/Sam/*.pm, and
Bio/DB/Bam/*.pm.

=head2 Bio::DB::Sam Constructor and basic accessors

=over 4

=item $sam = Bio::DB::Sam->new(%options)

The Bio::DB::Sam object combines a Fasta file of the reference
sequences with a BAM file to allow for convenient retrieval of
human-readable sequence IDs and reference sequences. The new()
constructor accepts a -name=>value style list of options as
follows:

  Option         Description
  ------         -------------

  -bam           Path to the BAM file that contains the
                   alignments (required). When using samtools 0.1.6
                   or higher, an http: or ftp: URL is accepted.

  -fasta         Path to the Fasta file that contains
                   the reference sequences (optional). Alternatively,
                   you may pass any object that supports a seq()
                   or fetch_seq() method and takes the three arguments
                   ($seq_id,$start,$end).

  -expand_flags  A boolean value. If true then the standard
                   alignment flags will be broken out as 
                   individual tags such as 'M_UNMAPPED' (default
                   false).

  -split_splices A boolean value. If true, then alignments that
                  are split across splices will be broken out
                  into a single alignment containing two sub-
                  alignments (default false).

  -split          The same as -split_splices.

  -autoindex      Create a BAM index file if one does not exist
                   or the current one has a modification date
                   earlier than the BAM file.

An example of a typical new() constructor invocation is:
 
  $sam = Bio::DB::Sam->new(-fasta => '/home/projects/genomes/hu17.fa',
                           -bam   => '/home/projects/alignments/ej88.bam',
                           -expand_flags  => 1,
                           -split_splices => 1);

If the B<-fasta> argument is present, then you will be able to use the
interface to fetch the reference sequence's bases. Otherwise, calls
that return the reference sequence will return sequences consisting
entirely of "N".

B<-expand_flags> option, if true, has the effect of turning each of
the standard SAM flags into a separately retrievable B<tag> in the
Bio::SeqFeatureI interface. Otherwise, the standard flags will be
concatenated in easily parseable form as a tag named "FLAGS". See
get_all_tags() and get_tag_values() for more information.

Any two-letter extension flags, such as H0 or H1, will always appear
as separate tags regardless of the setting.

B<-split_splices> has the effect of breaking up alignments that
contain an "N" operation into subparts for more convenient
manipulation. For example, if you have both paired reads and spliced
alignments in the BAM file, the following code shows the subpart
relationships:

  $pair        = $sam->get_feature_by_name('E113:01:01:23');
  @mates       = $pair->get_SeqFeatures;
  @mate1_parts = $mates[0]->get_SeqFeatures;
  @mate2_parts = $mates[1]->get_SeqFeatures;

Because there is some overhead to splitting up the spliced alignments,
this option is false by default.

B<Remote access> to BAM files located on an HTTP or FTP server is
possible when using the Samtools library version 0.1.6 or
higher. Simply replace the path to the BAM file with the appropriate
URL. Note that incorrect URLs may lead to a core dump.

It is not currently possible to refer to a remote FASTA file. These
will have to be downloaded locally and indexed before using.

=item $flag = $sam->expand_flags([$new_value])

Get or set the expand_flags option. This can be done after object
creation and will have an immediate effect on all alignments fetched
from the BAM file.

=item $flag = $sam->split_splices([$new_value])

Get or set the split_splices option. This can be done after object
creation and will affect all alignments fetched from the BAM file
B<subsequently.>

=item $header = $sam->header

Return the Bio::DB::Bam::Header object associated with the BAM
file. You can manipulate the header using the low-level API.

=item $bam    = $sam->bam

Returns the low-level Bio::DB::Bam object associated with the opened
file.

=item $fai    = $sam->fai

Returns the Bio::DB::Sam::Fai object associated with the Fasta
file. You can then manipuate this object with the low-level API.

B<The index will be built automatically for you if it does not already
exist.> If index building is necessarily, the process will need write
privileges to the same directory in which the Fasta file resides.> If
the process does not have write permission, then the call will fail.
Unfortunately, the BAM library does not do great error recovery for
this condition, and you may experience a core dump. This is not
trappable via an eval {}.

=item $bai    = $sam->bam_index

Return the Bio::DB::Bam::Index object associated with the BAM file. 

B<The BAM file index will be built automatically for you if it does
not already exist.> In addition, if the BAM file is not already sorted
by chromosome and coordinate, it will be sorted automatically, an
operation that consumes significant time and disk space. The current
process must have write permission to the directory in which the BAM
file resides in order for this to work.> In case of a permissions
problem, the Perl library will catch the error and die. You can trap
it with an eval {}.

=item $sam->clone

Bio::DB::SAM objects are not stable across fork() operations. If you
fork, you must call clone() either in the parent or the child process
before attempting to call any methods.

=back

=head2 Getting information about reference sequences

The Bio::DB::Sam object provides the following methods for getting
information about the reference sequence(s) contained in the
associated Fasta file.

=over 4

=item @seq_ids = $sam->seq_ids

Returns an unsorted list of the IDs of the reference sequences (known
elsewhere in this document as seq_ids). This is the same as the
identifier following the ">" sign in the Fasta file (e.g. "chr1").

=item $num_targets = $sam->n_targets

Return the number of reference sequences.

=item $length = $sam->length('seqid')

Returns the length of the reference sequence named "seqid".

=item $seq_id = $sam->target_name($tid)

Translates a numeric target ID (TID) returned by the low-level API
into a seq_id used by the high-level API.

=item $length = $sam->target_len($tid)

Translates a numeric target ID (TID) from the low-level API to a
sequence length.

=item $dna    = $sam->seq($seqid,$start,$end)

Returns the DNA across the region from start to end on reference
seqid. Note that this is a string, not a Bio::PrimarySeq object. If
no -fasta path was passed when the sam object was created, then you
will receive a series of N nucleotides of the requested length.

=back

=head2 Creating and querying segments

Bio::DB::Sam::Segment objects refer regions on the reference
sequence. They can be used to retrieve the sequence of the reference,
as well as alignments that overlap with the region.

=over 4

=item $segment = $sam->segment($seqid,$start,$end);

=item $segment = $sam->segment(-seq_id=>'chr1',-start=>5000,-end=>6000);

Segments are created using the Bio:DB::Sam->segment() method. It can
be called using one to three positional arguments corresponding to the
seq_id of the reference sequence, and optionally the start and end
positions of a subregion on the sequence. If the start and/or end are
undefined, they will be replaced with the beginning and end of the
sequence respectively.

Alternatively, you may call segment() with named -seq_id, -start and
-end arguments.

All coordinates are 1-based.

=item $seqid = $segment->seq_id

Return the segment's sequence ID.

=item $start = $segment->start

Return the segment's start position.

=item $end  = $segment->end

Return the segment's end position.

=item $strand = $segment->strand

Return the strand of the segment (always 0).

=item $length = $segment->length

Return the length of the segment.

=item $dna    = $segment->dna

Return the DNA string for the reference sequence under this segment.

=item $seq    = $segment->seq

Return a Bio::PrimarySeq object corresponding to the sequence of the
reference under this segment. You can get the actual DNA string in
this redundant-looking way:

 $dna = $segment->seq->seq

The advantage of working with a Bio::PrimarySeq object is that you can
perform operations on it, including taking its reverse complement and
subsequences.

=item @alignments = $segment->features(%args)

Return alignments that overlap the segment in the associated BAM
file. The optional %args list allows you to filter features by name,
tag or other attributes. See the documentation of the
Bio::DB::Sam->features() method for the full list of options. Here are
some typical examples:

 # get all the overlapping alignments
 @all_alignments = $segment->features;  

 # get an iterator across the alignments
 my $iterator     = $segment->features(-iterator=>1);
 while (my $align = $iterator->next_seq) { do something }

 # get a TAM filehandle across the alignments
 my $fh           = $segment->features(-fh=>1);
 while (<$fh>) { print }

 # get only the alignments with unmapped mates
 my @unmapped    = $segment->features(-flags=>{M_UNMAPPED=>1});

 # get coverage across this region
 my ($coverage)       = $segment->features('coverage');
 my @data_points      = $coverage->coverage;

 # grep through features using a coderef
 my @reverse_alignments = $segment->features(
                           -filter => sub {
                                  my $a = shift;
                                  return $a->strand < 0;
                               });

=item $tag = $segment->primary_tag

=item $tag = $segment->source_tag

Return the strings "region" and "sam/bam" respectively. These methods
allow the segment to be passed to BioPerl methods that expect
Bio::SeqFeatureI objects.

=item $segment->name, $segment->display_name, $segment->get_SeqFeatures, $segment->get_tag_values

These methods are provided for Bio::SeqFeatureI compatibility and
don't do anything of interest.

=back

=head2 Retrieving alignments, mate pairs and coverage information

The features() method is an all-purpose tool for retrieving alignment
information from the SAM/BAM database. In addition, the methods
get_features_by_name(), get_features_by_location() and others provide
convenient shortcuts to features().

These methods either return a list of features, an iterator across a
list of features, or a filehandle opened on a pseudo-TAM file.

=over 4

=item @features   = $sam->features(%options)

=item $iterator   = $sam->features(-iterator=>1,%more_options)

=item $filehandle = $sam->features(-fh=>1,%more_options)

=item @features   = $sam->features('type1','type2'...)

This is the all-purpose interface for fetching alignments and other
types of features from the database. Arguments are a -name=>value
option list selected from the following list of options:

  Option         Description
  ------         -------------

  -type          Filter on features of a given type. You may provide
  		 either a scalar typename, or a reference to an 
                 array of desired feature types. Valid types are
                 "match", "read_pair", "coverage" and "chromosome."
		 See below for a full explanation of feature types.

  -name          Filter on reads with the designated name. Note that
                 this can be a slow operation unless accompanied by
                 the feature location as well.

  -seq_id        Filter on features that align to seq_id between start
  -start         and end. -start and -end must be used in conjunction
  -end           with -seq_id. If -start and/or -end are absent, they
                 will default to 1 and the end of the reference
		 sequence, respectively.

  -flags         Filter features that match a list of one or more
                 flags. See below for the format.

  -attributes    The same as -flags, for compatibility with other
  -tags          APIs.
 
  -filter        Filter on features with a coderef. The coderef will
                 receive a single argument consisting of the feature
                 and should return true to keep the feature, or false
                 to discard it.

  -iterator      Instead of returning a list of features, return an
                 iterator across the results. To retrieve the results,
		 call the iterator's next_seq() method repeatedly
                 until it returns undef to indicate that no more
		 matching features remain.

  -fh            Instead of returning a list of features, return a
                 filehandle. Read from the filehandle to retrieve 
                 each of the results in TAM format, one alignment
                 per line read. This only works for features of type
                 "match."

The high-level API introduces the concept of a B<feature "type"> in order
to provide several convenience functions. You specify types by using
the optional B<-type> argument. The following types are currently
supported:

B<match>. The "match" type corresponds to the unprocessed SAM
alignment. It will retrieve single reads, either mapped or
unmapped. Each match feature's primary_tag() method will return the
string "match." The features returned by this call are of type
Bio::DB::Bam::AlignWrapper.

B<read_pair>. The "paired_end" type causes the sam interface to find
and merge together mate pairs. Fetching this type of feature will
yield a series of Bio::SeqFeatureI objects, each as long as the total
distance on the reference sequence spanned by the mate pairs. The
top-level feature is of type Bio::SeqFeature::Lite; it contains two
Bio::DB::Bam::AlignWrapper subparts. 

Call get_SeqFeatures() to get the two individual reads. Example:

 my @pairs    = $sam->features(-type=>'read_pair');
 my $p        = $pairs[0];
 my $i_length = $p->length;
 my @ends     = $p->get_SeqFeatures;
 my $left     = $ends[0]->start;
 my $right    = $ends[1]->end;
 
B<coverage>. The "coverage" type causes the sam interface to calculate
coverage across the designated region. It only works properly if
accompanied by the desired location of the coverage graph; -seq_id is
a mandatory argument for coverage calculation, and -start and -end are
optional. The call will return a single Bio::SeqFeatureI object whose
primary_tag() is "coverage." To recover the coverage data, call the
object's coverage() method to obtain an array (list context) or
arrayref (scalar context) of coverage counts across the region of
interest:

 my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>'seq1');
 my @data       = $coverage->coverage;
 my $total;
 for (@data) { $total += $_ }
 my $average_coverage = $total/@data;

By default the coverage graph will be at the base pair level. So for a
region 5000 bp wide, coverage() will return an array or arrayref with
exactly 5000 elements. However, you also have the option of
calculating the coverage across larger bins. Simply append the number
of intervals you are interested to the "coverage" typename. For
example, fetching "coverage:500" will return a feature whose
coverage() method will return the coverage across 500 intervals.

B<chromosome> or B<region>. The "chromosome" or "region" type are
interchangeable. They ask the sam interface to construct
Bio::DB::Sam::Segment representing the reference sequences. These two
calls give similar results:

 my $segment = $sam->segment('seq2',1=>500);
 my ($seg)   = $sam->features(-type=>'chromosome',
		              -seq_id=>'seq2',-start=>1,-end=>500);

Due to an unresolved bug, you cannot fetch chromosome features in the
same call with matches and other feature types call. Specifically,
this works as expected:

 my @chromosomes = $sam->features (-type=>'chromosome');

But this doesn't (as of 18 June 2009):

 my @chromosomes_and_matches = $sam->features(-type=>['match','chromosome']);

If no -type argument is provided, then features() defaults to finding
features of type "match."

You may call features() with a plain list of strings (positional
arguments, not -type=>value arguments). This will be interpreted as a
list of feature types to return:

 my ($coverage) = $sam->features('coverage')

For a description of the methods available in the features returned
from this call, please see L<Bio::SeqfeatureI> and
L<Bio::DB::Bam::Alignment>.

You can B<filter> "match" and "read_pair" features by name, location
and/or flags. The name and flag filters are not very efficient. Unless
they are combined with a location filter, they will initiate an
exhaustive search of the BAM database.

Name filters are case-insensitive, and allow you to use shell-style
"*" and "?"  wildcards. Flag filters created with the B<-flag>,
B<-attribute> or B<-tag> options have the following syntax:

 -flag => { FLAG_NAME_1 => ['list','of','possible','values'],
            FLAG_NAME_2 => ['list','of','possible','values'],
            ...
          }

The value of B<-flag> is a hash reference in which the keys are flag
names and the values are array references containing lists of
acceptable values. The list of values are OR'd with each other, and
the flag names are AND'd with each other.

The B<-filter> option provides a completely generic filtering
interface. Provide a reference to a subroutine. It will be called
once for each potential feature. Return true to keep the feature, or
false to discard it. Here is an example of how to find all matches
whose alignment quality scores are greater than 80.

 @features = $sam->features(-filter=>sub {shift->qual > 80} );

By default, features() returns a list of all matching features. You
may instead request an iterator across the results list by passing
-iterator=>1. This will give you an object that has a single method,
next_seq():

  my $high_qual  = $sam->features(-filter  => sub {shift->qual > 80},
                                  -iterator=> 1 );
  while (my $feature = $high_qual->next_seq) {
    # do something with the alignment
  }

Similarly, by passing a true value to the argument B<-fh>, you can
obtain a filehandle to a virtual TAM file. This only works with the
"match" feature type:

  my $high_qual  = $sam->features(-filter  => sub {shift->qual > 80},
                                  -fh      => 1 );
  while (my $tam_line = <$high_qual>) {
    chomp($tam_line);
    # do something with it
  }

=item @features   = $sam->get_features_by_name($name)

Convenience method. The same as calling $sam->features(-name=>$name);

=item $feature    = $sam->get_feature_by_name($name)

Convenience method. The same as ($sam->features(-name=>$name))[0].

=item @features   = $sam->get_features_by_location($seqid,$start,$end)

Convenience method. The same as calling
$sam->features(-seq_id=>$seqid,-start=>$start,-end=>$end).

=item @features   = $sam->get_features_by_flag(%flags)

Convenience method. The same as calling
$sam->features(-flags=>\%flags). This method is also called
get_features_by_attribute() and get_features_by_tag(). Example:

 @features = $sam->get_features_by_flag(H0=>1)

=item $feature    = $sam->get_feature_by_id($id)

The high-level API assigns each feature a unique ID composed of its
read name, position and strand and returns it when you call the
feature's primary_id() method. Given that ID, this method returns the
feature.

=item $iterator   = $sam->get_seq_stream(%options)

Convenience method. This is the same as calling
$sam->features(%options,-iterator=>1).

=item $fh         = $sam->get_seq_fh(%options)

Convenience method. This is the same as calling
$sam->features(%options,-fh=>1).

=item $fh         = $sam->tam_fh

Convenience method. It is the same as calling $sam->features(-fh=>1).

=item @types      = $sam->types

This method returns the list of feature types (e.g. "read_pair")
returned by the current version of the interface.

=back

=head2 The generic fetch() and pileup() methods

Lastly, the high-level API supports two methods for rapidly traversing
indexed BAM databases.

=over 4

=item $sam->fetch($region,$callback)

This method, which is named after the native bam_fetch() function in
the C interface, traverses the indicated region and invokes a callback
code reference on each match. Specify a region using the BAM syntax
"seqid:start-end", or either of the alternative syntaxes
"seqid:start..end" and "seqid:start,end". If start and end are absent,
then the entire reference sequence is traversed. If end is absent,
then the end of the reference sequence is assumed.

The callback will be called repeatedly with a
Bio::DB::Bam::AlignWrapper on the argument list.

Example:

  $sam->fetch('seq1:600-700',
              sub {
                my $a = shift;
                print $a->display_name,' ',$a->cigar_str,"\n";
              });

Note that the fetch() operation works on reads that B<overlap> the
indicated region. Therefore the callback may be called for reads that
align to the reference at positions that start before or end after the
indicated region.

=item $sam->pileup($region,$callback [,$keep_level])

This method, which is named after the native bam_lpileupfile()
function in the C interfaces, traverses the indicated region and
generates a "pileup" of all the mapped reads that cover it. The
user-provided callback function is then invoked on each position of
the alignment along with a data structure that provides access to the
individual aligned reads.

As with fetch(), the region is specified as a string in the format
"seqid:start-end", "seqid:start..end" or "seqid:start,end".

The callback is a coderef that will be invoked with three arguments:
the seq_id of the reference sequence, the current position on the
reference (in 1-based coordinates!), and a reference to an array of
Bio::DB::Bam::Pileup objects. Here is the typical call signature:

  sub {
       my ($seqid,$pos,$pileup) = @_;
       # do something
  }

For example, if you call pileup on the region "seq1:501-600", then the
callback will be invoked for all reads that overlap the indicated
region. The first invocation of the callback will typically have a
$pos argument somewhat to the left of the desired region and the last
call will be somewhat to the right. You may wish to ignore positions
that are outside of the requested region. Also be aware that the
reference sequence position uses 1-based coordinates, which is
different from the low-level interface, which use 0-based coordinates.

The optional $keep_level argument, if true, asks the BAM library to
keep track of the level of the read in the multiple alignment, an
operation that generates some overhead. This is mostly useful for text
alignment viewers, and so is off by default.

The size of the $pileup array reference indicates the read coverage
at that position. Here is a simple average coverage calculator:

 my $depth      = 0;
 my $positions  = 0;
 my $callback = sub {
         my ($seqid,$pos,$pileup) = @_;
         next unless $pos >= 501 && $pos <= 600;
         $positions++;
         $depth += @$pileup;
 }
 $sam->pileup('seq1:501-600',$callback);
 print "coverage = ",$depth/$positions;

Each Bio::DB::Bam::Pileup object describes the position of a read in
the alignment. Briefly, Bio::DB::Bam::Pileup has the following
methods:

 $pileup->alignment  The alignment at this level (a
                     Bio::DB::Bam::AlignWrapper object).
 
 $pileup->qpos   The position of the read base at the pileup site,
                 in 0-based coordinates.

 $pileup->pos    The position of the read base at the pileup site,
                 in 1-based coordinates;

 $pileup->level  The level of the read in the multiple alignment
                 view. Note that this field is only valid when
                 $keep_level is true.

 $pileup->indel  Length of the indel at this position: 0 for no indel, positive
                 for an insertion (relative to the reference), negative for a
                 deletion (relative to the reference.)

 $pileup->is_del True if the base on the padded read is a deletion.

 $pileup->is_head Undocumented field in the bam.h header file.

 $pileup->is_tail Undocumented field in the bam.h header file.

See L</Examples> for a very simple SNP caller.

=item $sam->fast_pileup($region,$callback [,$keep_level])

This is identical to pileup() except that the pileup object returns
low-level Bio::DB::Bam::Alignment objects rather than the higher-level
Bio::DB::Bam::AlignWrapper objects. This makes it roughly 50% faster,
but you lose the align objects' seq_id() and get_tag_values()
methods. As a compensation, the callback receives an additional
argument corresponding to the Bio::DB::Sam object. You can use this to
create AlignWrapper objects on an as needed basis:

 my $callback = sub {
    my($seqid,$pos,$pileup,$sam) = @_;
    for my $p (@$pileup) {
       my $alignment = $p->alignment;
       my $wrapper   = Bio::DB::Bam::AlignWrapper->new($alignment,$sam);
       my $has_mate  = $wrapper->get_tag_values('PAIRED');
    }
  };

=back

The next sections correspond to the low-level API, which let you
create and manipulate Perl objects that correspond directly to data
structures in the C interface. A major difference between the high and
low level APIs is that in the high-level API, the reference sequence
is identified using a human-readable seq_id. However, in the low-level
API, the reference is identified using a numeric target ID
("tid"). The target ID is established during the creation of the BAM
file and is a small 0-based integer index. The Bio::DB::Bam::Header
object provides methods for converting from seq_ids to tids.

=head2 Indexed Fasta Files

These methods relate to the BAM library's indexed Fasta (".fai")
files.

=over 4

=item $fai = Bio::DB::Sam::Fai->load('/path/to/file.fa')

Load an indexed Fasta file and return the object corresponding to
it. If the index does not exist, it will be created
automatically. Note that you pass the path to the Fasta file, not the
index.

For consistency with Bio::DB::Bam->open() this method is also called
open().

=item $dna_string = $fai->fetch("seqid:start-end")

Given a sequence ID contained in the Fasta file and optionally a
subrange in the form "start-end", finds the indicated subsequence and
returns it as a string.

=back

=head2 TAM Files

These methods provide interfaces to the "TAM" text version of SAM
files; they often have a .sam extension.

=over 4

=item $tam = Bio::DB::Tam->open('/path/to/file.sam')

Given the path to a SAM file, opens it for reading. The file can be
compressed with gzip if desired.

=item $header = $tam->header_read2('/path/to/file.fa.fai')

Create and return a Bio::DB::Bam::Header object from the information
contained within the indexed Fasta file of the reference
sequences. Note that you have to pass the path to the .fai file, and
not the .fa file. The header object contains information on the
reference sequence names and lengths.

=item $bytes = $tam->read1($header,$alignment)

Given a Bio::DB::Bam::Header object, such as the one created by
header_read2(), and a Bio::DB::Bam::Alignment object created by
Bio::DB::Bam::Alignment->new(), reads one line of alignment information
into the alignment object from the TAM file and returns a status
code. The result code will be the number of bytes read.

=back

=head2 BAM Files

These methods provide interfaces to the "BAM" binary version of
SAM. They usually have a .bam extension.

=over 4

=item $bam = Bio::DB::Bam->open('/path/to/file.bam' [,$mode])

Open up the BAM file at the indicated path. Mode, if present, must be
one of the file stream open flags ("r", "w", "a", "r+", etc.). If
absent, mode defaults to "r".

Note that Bio::DB::Bam objects are not stable across fork()
operations. If you fork, and intend to use the object in both parent
and child, you must reopen the Bio::DB::Bam in either the child or the
parent (but not both) before attempting to call any of the object's
methods.

The path may be an http: or ftp: URL, in which case a copy of the
index file will be downloaded to the current working directory (see
below) and all accesses will be performed on the remote BAM file.

Example:

   $bam = Bio::DB::Bam->open('http://some.site.com/nextgen/chr1_bowtie.bam');

=item $header = $bam->header()

Given an open BAM file, return a Bio::DB::Bam::Header object
containing information about the reference sequence(s).

=item $status_code = $bam->header_write($header)

Given a Bio::DB::Bam::Header object and a BAM file opened in write
mode, write the header to the file. If the write fails the process
will be terminated at the C layer. The result code is (currently)
always zero.

=item $integer = $bam->tell()

Return the current position of the BAM file read/write pointer.

=item $bam->seek($integer)

Set the current position of the BAM file read/write pointer.

=item $alignment = $bam->read1()

Read one alignment from the BAM file and return it as a
Bio::DB::Bam::Alignment object.

=item $bytes = $bam->write1($alignment)

Given a BAM file that has been opened in write mode and a
Bio::DB::Bam::Alignment object, write the alignment to the BAM file
and return the number of bytes successfully written.

=item Bio::DB::Bam->sort_core($by_qname,$path,$prefix,$max_mem)

Attempt to sort a BAM file by chromosomal location or name and create a
new sorted BAM file. Arguments are as follows:

 Argument      Description
 --------      -----------

 $by_qname     If true, sort by read name rather than chromosomal
               location.

 $path         Path to the BAM file

 $prefix       Prefix to use for the new sorted file. For example,
               passing "foo" will result in a BAM file named 
	       "foo.bam".

 $max_mem      Maximum core memory to use for the sort. If the sort
               requires more than this amount of memory, intermediate
               sort files will be written to disk. The default, if not
               provided is 500M.

=back

=head2 BAM index methods

The Bio::DB::Bam::Index object provides access to BAM index (.bai)
files.

=over 4

=item $status_code = Bio::DB::Bam->index_build('/path/to/file.bam')

Given the path to a .bam file, this function attempts to build a
".bai" index. The process in which the .bam file exists must be
writable by the current process and there must be sufficient disk
space for the operation or the process will be terminated in the C
library layer. The result code is currently always zero, but in the
future may return a negative value to indicate failure.

=item $index = Bio::DB::Bam->index('/path/to/file.bam',$reindex)

Attempt to open the index for the indicated BAM file. If $reindex is
true, and the index either does not exist or is out of date with
respect to the BAM file (by checking modification dates), then attempt
to rebuild the index. Will throw an exception if the index does not
exist or if attempting to rebuild the index was unsuccessful.

=item $index = Bio::DB::Bam->index_open('/path/to/file.bam')

Attempt to open the index file for a BAM file, returning a
Bio::DB::Bam::Index object. The filename path to use is the .bam file,
not the .bai file.

=item $index = Bio::DB::Bam->index_open_in_safewd('/path/to/file.bam' [,$mode])

When opening a remote BAM file, you may not wish for the index to be
downloaded to the current working directory. This version of index_open
copies the index into the directory indicated by the TMPDIR
environment variable or the system-defined /tmp directory if not
present. You may change the environment variable just before the call
to change its behavior.

=item $code = $index->fetch($bam,$tid,$start,$end,$callback [,$callback_data])

This is the low-level equivalent of the $sam->fetch() function
described for the high-level API. Given a open BAM file object, the
numeric ID of the reference sequence, start and end ranges on the
reference, and a coderef, this function will traverse the region and
repeatedly invoke the coderef with each Bio::DB::Bam::Alignment
object that overlaps the region.

Arguments:

 Argument      Description
 --------      -----------

 $bam          The Bio::DB::Bam object that corresponds to the
               index object.

 $tid          The target ID of the reference sequence. This can
               be obtained by calling $header->parse_region() with
               an appropriate opened Bio::DB::Bam::Header object.

 $start        The start and end positions of the desired range on
               the reference sequence given by $tid, in 0-based 
 $end          coordinates. Like the $tid, these can be obtained from
               $header->parse_region().

 $callback     A coderef that will be called for each read overlapping
               the designated region.

 $callback_data  Any arbitrary Perl data that you wish to pass to the
               $callback (optional).

The coderef's call signature should look like this:

  my $callback = sub {
                    my ($alignment,$data) = @_;
                    ...
                 }

The first argument is a Bio::DB::Bam::Alignment object. The second is
the callback data (if any) passed to fetch().

Fetch() returns an integer code, but its meaning is not described in
the SAM/BAM C library documentation.

=item $index->pileup($bam,$tid,$start,$end,$callback [,$callback_data])

This is the low-level version of the pileup() method, which allows you
to invoke a coderef for every position in a BAM alignment. Arguments
are:

 Argument      Description
 --------      -----------

 $bam          The Bio::DB::Bam object that corresponds to the
               index object.

 $tid          The target ID of the reference sequence. This can
               be obtained by calling $header->parse_region() with
               an appropriate opened Bio::DB::Bam::Header object.

 $start        The start and end positions of the desired range on
               the reference sequence given by $tid, in 0-based 
 $end          coordinates. Like the $tid, these can be obtained from
               $header->parse_region().

 $callback     A coderef that will be called for each position of the
               alignment across the designated region.

 $callback_data  Any arbitrary Perl data that you wish to pass to the
               $callback (optional).

The callback will be invoked with four arguments corresponding to the
numeric sequence ID of the reference sequence, the B<zero-based>
position on the alignment, an arrayref of Bio::DB::Bam::Pileup
objects, and the callback data, if any. A typical call signature will
be this:

 $callback = sub {
       my ($tid,$pos,$pileups,$callback_data) = @_;
       for my $pileup (@$pileups) {
          # do something
       };

Note that the position argument is zero-based rather than 1-based, as
it is in the high-level API.

The Bio::DB::Bam::Pileup object was described earlier in the
description of the high-level pileup() method.

=item $coverage = $index->coverage($bam,$tid,$start,$end [,$bins])

Calculate coverage for the region on the target sequence given by $tid
between positions $start and $end (zero-based coordinates). This
method will return an array reference equal to the size of the region
(by default). Each element of the array will be an integer indicating
the number of reads aligning over that position. If you provide an
option binsize in $bins, the array will be $bins elements in length,
and each element will contain the average coverage over that region as
a floating point number.

=back

=head2 BAM header methods

The Bio::DB::Bam::Header object contains information regarding the
reference sequence(s) used to construct the corresponding TAM or BAM
file. It is most frequently used to translate between numeric target
IDs and human-readable seq_ids. Headers can be created either from
reading from a .fai file with the Bio::DB::Tam->header_read2() method,
or by reading from a BAM file using Bio::DB::Bam->header(). You can
also create header objects from scratch, although there is not much
that you can do with such objects at this point.

=over 4

=item $header = Bio::DB::Bam::Header->new()

Return a new, empty, header object.

=item $n_targets = $header->n_targets

Return the number of reference sequences in the database.

=item $name_arrayref = $header->target_name

Return a reference to an array of reference sequence names,
corresponding to the high-level API's seq_ids.

To convert from a target ID to a seq_id, simply index into this array:

 $seq_id = $header->target_name->[$tid];

=item $length_arrayref = $header->target_len

Return a reference to an array of reference sequence lengths. To get
the length of the sequence corresponding to $tid, just index into the
array returned by target_len():

 $length = $header->target_len->[$tid];

=item $text = $header->text
=item $header->text("new value")

Read the text portion of the BAM header. The text can be replaced by
providing the replacement string as an argument. Note that you should
follow the header conventions when replacing the header text. No
parsing or other error-checking is performed.

=item ($tid,$start,$end) = $header->parse_region("seq_id:start-end")

Given a string in the format "seqid:start-end" (using a human-readable
seq_id and 1-based start and end coordinates), parse the string and
return the target ID and start and end positions in 0-based
coordinates. If the range is omitted, then the start and end
coordinates of the entire sequence is returned. If only the end
position is omitted, then the end of the sequence is assumed.

=item $header->view1($alignment)

This method will accept a Bio::DB::Bam::Alignment object, convert it
to a line of TAM output, and write the output to STDOUT. In the
low-level API there is currently no way to send the output to a
different filehandle or capture it as a string.

=back

=head2 Bio::DB::Bam::Pileup methods

An array of Bio::DB::Bam::Pileup object is passed to the pileup()
callback for each position of a multi-read alignment. Each pileup
object contains information about the alignment of a single read at a
single position.

=over 4

=item $alignment = $pileup->alignment

Return the Bio::DB::Bam::Alignment object at this level. This provides
you with access to the aligning read.

=item $alignment = $pileup->b

An alias for alignment(), provided for compatibility with the C API.

=item $pos = $pileup->qpos

The position of the aligning base in the read in zero-based
coordinates.

=item $pos = $pileup->pos

The position of the aligning base in 1-based coordinates.

=item $level = $pileup->level

The "level" of the read in the BAM-generated text display of the
alignment.

=item $indel = $pileup->indel

Length of the indel at this position: 0 for no indel, positive for an
insertion (relative to the reference), negative for a deletion
(relative to the reference sequence.)

=item $flag = $pileup->is_del

True if the base on the padded read is a deletion.

=item $flag = $pileup->is_head

=item $flag = $pileup->is_del

These fields are undocumented in the BAM documentation, but are
exported to the Perl API just in case.

=back

=head2 The alignment objects

Please see L<Bio::DB::Bam::Alignment> for documentation of the
Bio::DB::Bam::Alignment and Bio::DB::Bam::AlignWrapper objects.

=cut

use strict;
use warnings;

use Carp 'croak';
use Bio::SeqFeature::Lite;
use Bio::PrimarySeq;

use base 'DynaLoader';
bootstrap Bio::DB::Sam;

use Bio::DB::Bam::Alignment;
use Bio::DB::Sam::Segment;
use Bio::DB::Bam::AlignWrapper;
use Bio::DB::Bam::PileupWrapper;
use Bio::DB::Bam::FetchIterator;
use Bio::DB::Bam::ReadIterator;

sub new {
    my $class         = shift;
    my %args          = $_[0] =~ /^-/ ? @_ : (-bam=>shift);
    my $bam_path      = $args{-bam}   or croak "-bam argument required";
    my $fa_path       = $args{-fasta};
    my $expand_flags  = $args{-expand_flags};
    my $split_splices = $args{-split} || $args{-split_splices};
    my $autoindex     = $args{-autoindex};

    # file existence checks
    unless ($class->is_remote($bam_path)) {
	-e $bam_path or croak "$bam_path does not exist";
	-r _  or croak "is not readable";
    }

    my $bam = Bio::DB::Bam->open($bam_path)      or croak "$bam_path open: $!";

    my $fai = $class->new_dna_accessor($fa_path) if $fa_path;

    my $self =  bless {
	fai           => $fai,
	bam           => $bam,
	bam_path      => $bam_path,
	fa_path       => $fa_path,
	expand_flags  => $expand_flags,
	split_splices => $split_splices,
	autoindex     => $autoindex,
    },ref $class || $class;
    $self->header;  # catch it

    return $self;
}

sub bam { shift->{bam} }

sub is_remote {
    my $self = shift;
    my $path = shift;
    return $path =~ /^(http|ftp):/;
}

sub clone {
    my $self = shift;
    $self->{bam} = Bio::DB::Bam->open($self->{bam_path})     if $self->{bam_path};
    $self->{fai} = $self->new_dna_accessor($self->{fa_path}) if $self->{fa_path};
}

sub header {
    my $self = shift;
    return $self->{header} ||= $self->{bam}->header;
}

sub fai { shift->{fai} }

sub new_dna_accessor {
    my $self     = shift;
    my $accessor  = shift;

    return unless $accessor;

    if (-e $accessor) {  # a file, assume it is a fasta file
	-r _  or croak "$accessor is not readable";
	my $a = Bio::DB::Sam::Fai->open($accessor)  or croak "$accessor open: $!"
	    or croak "Can't open FASTA file $accessor: $!";
	return $a;
    }

    if (ref $accessor && $self->can_do_seq($accessor)) {
	return $accessor;  # already built
    }

    return;
}

sub can_do_seq {
    my $self = shift;
    my $obj  = shift;
    return 
	UNIVERSAL::can($obj,'seq') ||
	UNIVERSAL::can($obj,'fetch_sequence');
}


sub seq {
    my $self = shift;
    my ($seqid,$start,$end) = @_;
    my $fai = $self->fai or return 'N' x ($end-$start+1);
    return $fai->can('seq')            ? $fai->seq($seqid,$start,$end) 
	  :$fai->can('fetch_sequence') ? $fai->fetch_sequence($seqid,$start,$end)
	  :'N' x ($end-$start+1);
}

sub expand_flags {
    my $self = shift;
    my $d    = $self->{expand_flags};
    $self->{expand_flags} = shift if @_;
    $d;
}

sub split_splices {
    my $self = shift;
    my $d    = $self->{split_splices};
    $self->{split_splices} = shift if @_;
    $d;
}

sub autoindex {
    my $self = shift;
    my $d    = $self->{autoindex};
    $self->{autoindex} = shift if @_;
    $d;
}

sub reset_read {
    my $self = shift;
    $self->{bam}->header;
}

sub n_targets {
    shift->header->n_targets;
}

sub target_name {
    my $self = shift;
    my $tid  = shift;
    $self->{target_name} ||= $self->header->target_name;
    return $self->{target_name}->[$tid];
}

sub target_len {
    my $self = shift;
    my $tid  = shift;
    $self->{target_len} ||= $self->header->target_len;
    return $self->{target_len}->[$tid];
}

sub seq_ids {
    my $self    = shift;
    return @{$self->header->target_name};
}

sub _cache_targets {
    my $self = shift;
    return $self->{targets} if exists $self->{targets};
    my @targets = map {lc $_} @{$self->header->target_name};
    my @lengths =             @{$self->header->target_len};
    my %targets;
    @targets{@targets}      = @lengths;  # just you try to figure out what this is doing!
    return $self->{targets} = \%targets;
}


sub length {
    my $self        = shift;
    my $target_name = shift;
    return $self->_cache_targets->{lc $target_name};
}

sub _fetch {
    my $self     = shift;
    my $region   = shift;
    my $callback = shift;

    my $header              = $self->{bam}->header;
    $region                 =~ s/\.\.|,/-/;

    my ($seqid,$start,$end) = $header->parse_region($region);

    return unless defined $seqid;
    my $index  = $self->bam_index;
    $index->fetch($self->{bam},$seqid,$start,$end,$callback,$self);
}

sub fetch {
    my $self     = shift;
    my $region   = shift;
    my $callback = shift;
    
    my $code     = sub {
	my ($align,$self) = @_;
	$callback->(Bio::DB::Bam::AlignWrapper->new($align,$self));
    };
    $self->_fetch($region,$code);
}

sub pileup {
    my $self   = shift;
    my ($region,$callback,$keep_level) = @_;

    my $header   = $self->header;
    $region      =~ s/\.\.|,/-/;
    my ($seqid,$start,$end) = $header->parse_region($region);
    return unless defined $seqid;

    my $refnames = $self->header->target_name;

    my $code = sub {
	my ($tid,$pos,$pileup) = @_;
	my $seqid = $refnames->[$tid];
	my @p = map {
	      Bio::DB::Bam::PileupWrapper->new($_,$self)
	      } @$pileup;
	$callback->($seqid,$pos+1,\@p);
    };

    my $index  = $self->bam_index;
    if ($keep_level) {
	$index->lpileup($self->{bam},$seqid,$start,$end,$code);
    } else {
	$index->pileup($self->{bam},$seqid,$start,$end,$code);
    }
}

sub fast_pileup {
    my $self   = shift;
    my ($region,$callback,$keep_level) = @_;

    my $header   = $self->header;
    $region      =~ s/\.\.|,/-/;
    my ($seqid,$start,$end) = $header->parse_region($region);
    return unless defined $seqid;

    my $refnames = $self->header->target_name;

    my $code = sub {
	my ($tid,$pos,$pileup) = @_;
	my $seqid = $refnames->[$tid];
	$callback->($seqid,$pos+1,$pileup,$self);
    };

    my $index  = $self->bam_index;
    if ($keep_level) {
	$index->lpileup($self->{bam},$seqid,$start,$end,$code);
    } else {
	$index->pileup($self->{bam},$seqid,$start,$end,$code);
    }
}

# segment returns a segment across the reference
# it will not work on a arbitrary aligned feature
sub segment {
    my $self = shift;
    my ($seqid,$start,$end) = @_;

    if ($_[0] =~ /^-/) {
	my %args = @_;
	$seqid = $args{-seq_id} || $args{-name};
	$start = $args{-start};
	$end   = $args{-stop}    || $args{-end};
    } else {
	($seqid,$start,$end) = @_;
    }

    my $targets = $self->_cache_targets;
    return unless exists $targets->{lc $seqid};

    $start = 1                     unless defined $start;
    $end   = $targets->{lc $seqid} unless defined $end;
    $start = 1 if $start < 1;
    $end   = $targets->{lc $seqid} if $end > $targets->{lc $seqid};

    return Bio::DB::Sam::Segment->new($self,$seqid,$start,$end);
}

sub get_features_by_location {
    my $self = shift;
    my %args;

    if ($_[0] =~ /^-/) { # named args
	%args = @_;
    } else {             # positional args
	$args{-seq_id} = shift;
	$args{-start}  = shift;
	$args{-end}    = shift;
    }
    $self->features(%args);
}

sub get_features_by_attribute {
  my $self       = shift;
  my %attributes = ref($_[0]) ? %{$_[0]} : @_;
  $self->features(-attributes=>\%attributes);
}

sub get_features_by_tag {
    shift->get_features_by_attribute(@_);
}

sub get_features_by_flag {
    shift->get_features_by_attribute(@_);
}

sub get_feature_by_name {
    my $self = shift;
    my %args;
    if ($_[0] =~ /^-/) {
	%args = @_;
    } else {
	$args{-name} = shift;
    }
    $self->features(%args);
}

sub get_features_by_name { shift->get_feature_by_name(@_) }

sub get_feature_by_id {
    my $self = shift;
    my $id   = shift;
    my ($name,$tid,$start,$end,$strand) = map {s/%3B/;/ig;$_} split ';',$id;
    return unless $name && defined $tid;
    my $seqid = $self->target_name($tid);
    my @features = $self->features(-name=>$name,
				   -seq_id=>$seqid,
				   -start=>$start,
				   -end=>$end,
				   -strand=>$strand);
    return unless @features;
    return $features[0];
}


sub get_seq_stream {
    my $self = shift;
    $self->features(@_,-iterator=>1);
}

sub get_seq_fh {
    my $self = shift;
    $self->features(@_,-fh=>1);
}

sub types {
    return qw(match read_pair coverage region chromosome);
}

sub features {
    my $self = shift;

    my %args;
    if (defined $_[0] && $_[0] !~ /^-/) {
	$args{-type} = \@_;
    } else {
	%args = @_;
    }

    my $seqid     = $args{-seq_id} || $args{-seqid};
    my $start     = $args{-start};
    my $end       = $args{-end}  || $args{-stop};
    my $types     = $args{-type} || $args{-types} || [];
    my $attributes = $args{-attributes} || $args{-tags} || $args{-flags};
    my $iterator  = $args{-iterator};
    my $fh        = $args{-fh};
    my $filter    = $args{-filter};
    my $max       = $args{-max_features};

    $types        = [$types] unless ref $types;
    $types        = [$args{-class}] if !@$types && defined $args{-class};
    my $use_index = defined $seqid;

    # we do some special casing to retrieve target (reference) sequences
    # if they are requested
     if (defined($args{-name})
 	&& (!@$types || $types->[0]=~/region|chromosome/) 
	 && !defined $seqid) {
 	my @results = $self->_segment_search(lc $args{-name});
 	return @results if @results;
     } elsif (@$types && $types->[0] =~ /region|chromosome/) {
 	return map {$self->segment($_)} $self->seq_ids;
     }

    my %seenit;
    my @types = grep {!$seenit{$_}++} ref $types ? @$types : $types;
    @types    = 'match' unless @types;

    # the filter is intended to be inserted into a closure
    # it will return undef from the closure unless the filter
    # criteria are satisfied
    if (!$filter) {
	$filter = '';
	$filter   .= $self->_filter_by_name(lc $args{-name})
	    if defined $args{-name};
	$filter   .= $self->_filter_by_attribute($attributes)
	    if defined $attributes;
    }

    # Special cases for unmunged data
        if (@types == 1 && $types[0] =~ /^match/) {

	# if iterator is requested, and no indexing is possible,
	# then we directly iterate through the database using read1()
	if ($iterator && !$use_index) {
	    $self->reset_read;
	    my $code = eval "sub {my \$a=shift;$filter;1}";
	    die $@ if $@;
	    return Bio::DB::Bam::ReadIterator->new($self,$self->{bam},$code);
	}

	# TAM filehandle retrieval is requested
	elsif ($fh) {
	    return $self->_features_fh($seqid,$start,$end,$filter);
	}

    }

    # otherwise we're going to do a little magic
    my ($features,@result);

    for my $t (@types) {

	if ($t =~ /^(match|read_pair)/) {

	    # fetch the features if type is 'match' or 'read_pair'
	    $features = $self->_filter_features($seqid,$start,$end,$filter,undef,$max);

	    # for "match" just return the alignments
	    if ($t =~ /^(match)/) {
		push @result,@$features;
	    } 

	    # otherwise aggregate mate pairs into two-level features
	    elsif ($t =~ /^read_pair/) {
		$self->_build_mates($features,\@result);
	    }
	    next;
	}

	# create a coverage graph if type is 'coverage'
	# specify coverage:N, to create a map of N bins
	# units are coverage per bp
	# resulting array will be stored in the "coverage" attribute
	if ($t =~ /^coverage:?(\d*)/) {
	    my $bins = $1;
	    push @result,$self->_coverage($seqid,$start,$end,$bins,$filter);
	}
	
    }

    return $iterator ? Bio::DB::Bam::FetchIterator->new(\@result,$self->last_feature_count)
	             : @result;
}

sub _filter_features {
    my $self = shift;
    my ($seqid,$start,$end,$filter,$do_tam_fh,$max_features) = @_;

    my @result;
    my $action = $do_tam_fh ? '\$self->header->view1($a)'
                            : $self->_push_features($max_features);

    my $user_code;
    if (ref ($filter) eq 'CODE') {
	$user_code = $filter;
	$filter = '';
    }

    my $callback = defined($seqid) ? <<INDEXED : <<NONINDEXED;
sub {
    my \$a = shift;
    $filter
    return unless defined \$a->start;
    $action;
}
INDEXED
sub {
    my \$a    = shift;
    $filter
    $action;
}
NONINDEXED
    ;

    my $code = eval $callback;
    die $@ if $@;

    if ($user_code) {
	my $new_callback = sub {
	    my $a = shift;
	    $code->($a) if $user_code->($a);
	};
	$self->_features($seqid,$start,$end,$new_callback);
    } else {
	$self->_features($seqid,$start,$end,$code);
    }

    return \@result;
}

sub _push_features {
    my $self = shift;
    my $max  = shift;

    # simple case -- no max specified. Will push onto an array called
    # @result.

    return 'push @result,Bio::DB::Bam::AlignWrapper->new($a,$self)'
	unless $max;

    $self->{_result_count} = 0;

    # otherwise we implement a simple subsampling
    my $code=<<END;
    my \$count = ++\$self->{_result_count};
    if (\@result < $max) {
	push \@result,Bio::DB::Bam::AlignWrapper->new(\$a,\$self);
    } else {
	\$result[rand \@result] = Bio::DB::Bam::AlignWrapper->new(\$a,\$self) 
	    if rand() < $max/\$count;
    }
END
    return $code;
}

sub last_feature_count { shift->{_result_count}||0 }

sub _features {
    my $self = shift;
    my ($seqid,$start,$end,$callback) = @_;

    if (defined $seqid) {
 	my $region = $seqid;
 	if (defined $start) { 
 	    $region   .= ":$start";
 	    $region   .= "-$end"   if defined $end;
 	}
 	$self->_fetch($region,$callback);
    } 

    else {
	$self->reset_read;
	while (my $b = $self->{bam}->read1) {
	    $callback->($b);
 	}
    }
}

# build mate pairs
sub _build_mates {
    my $self = shift;
    my ($src,$dest) = @_;

    my %read_pairs;
    for my $a (@$src) {
	my $name = $a->display_name;
	unless ($read_pairs{$name}) {
	    my $isize = $a->isize;
	    my $start = $isize >= 0 ? $a->start : $a->end+$isize+1;
	    my $end   = $isize <= 0 ? $a->end   : $a->start+$isize-1;
	    $read_pairs{$name} = 
		Bio::SeqFeature::Lite->new(
		    -display_name => $name,
		    -seq_id       => $a->seq_id,
		    -start => $start,
		    -end   => $end,
		    -type  => 'read_pair',
		    -class => 'read_pair',
		);
	}
	$read_pairs{$name}->add_SeqFeature($a);
    }
    push @$dest,values %read_pairs;
}

sub _coverage {
    my $self = shift;
    my ($seqid,$start,$end,$bins,$filter) = @_;

    # Currently filter is ignored. In reality, we should
    # turn filter into a callback and invoke it on each 
    # position in the pileup.
    croak "cannot calculate coverage unless a -seq_id is provided"
	unless defined $seqid;

    my $region = $seqid;
    if (defined $start) { 
	$region   .= ":$start";
	$region   .= "-$end"   if defined $end;
    }

    my $header     = $self->{bam}->header;
    my ($id,$s,$e) = $header->parse_region($region);
    return unless defined $id;

    # parse_region may return a very high value if no end specified
    $end   = $e >= 1<<29 ? $header->target_len->[$id] : $e;
    $start = $s+1;
    $bins ||= $end-$start+1;

    my $index      = $self->bam_index;
    my $coverage   = $index->coverage($self->{bam},
				      $id,$s,$e,
				      $bins);

    return Bio::SeqFeature::Coverage->new(
	-display_name => "$seqid coverage",
	-seq_id       => $seqid,
	-start        => $start,
	-end          => $end,
	-strand       => 0,
	-type         => "coverage:$bins",
	-class        => "coverage:$bins",
	-attributes   => { coverage => [$coverage] }
    );
}

sub _segment_search {
    my $self = shift;
    my $name = shift;

    my $targets = $self->_cache_targets;
    return $self->segment($name) if $targets->{$name};

    if (my $regexp = $self->_glob_match($name)) {
	my @results = grep {/^$regexp$/i} keys %$targets;
	return map {$self->segment($_)} @results;
    }

    return;
}

sub bam_index {
    my $self = shift;
    return $self->{bai} ||= Bio::DB::Bam->index($self->{bam_path},$self->autoindex);
}

sub _features_fh {
    my $self  = shift;
    my ($seqid,$start,$end,$filter) = @_;

    my $result = open my $fh,"-|";
    if (!$result) {  # in child
	$self->_filter_features($seqid,$start,$end,$filter,'do_fh'); # will print TAM to stdout
	exit 0;
    }
    return $fh;
    
}

sub tam_fh {
    my $self   = shift;
    return $self->features(-fh=>1);
}

# return a fragment of code that will be placed in the eval "" filter
# to eliminate alignments that don't match by name
sub _filter_by_name {
    my $self = shift;
    my $name = shift;

    my $frag = "my \$name=\$a->qname; defined \$name or return; ";

    if (my $regexp = $self->_glob_match($name)) {
	$frag .= "return unless \$name =~ /^$regexp\$/i;\n";
    } else {
	$frag .= "return unless lc \$name eq '$name';\n";
    }
}

# return a fragment of code that will be placed in the eval "" filter
# to eliminate alignments that don't match by attribute
sub _filter_by_attribute {
    my $self       = shift;
    my $attributes = shift;
    my $result;
    for my $tag (keys %$attributes) {
	$result .= "my \$value = lc \$a->get_tag_values('$tag');\n";
	$result .= "return unless defined \$value;\n";
	my @comps = ref $attributes->{$tag} eq 'ARRAY' 
	    ? @{$attributes->{$tag}} 
	    : $attributes->{$tag};
	my @matches;
	for my $c (@comps) {
	    if ($c =~ /^[+-]?[\deE.]+$/) { # numeric-looking argument
		push @matches,"CORE::length \$value && \$value == $c";
	    }
	    elsif (my $regexp = $self->_glob_match($c)) {
		push @matches,"\$value =~ /^$regexp\$/i";
	    }
	    else {
		push @matches,"\$value eq lc '$c'";
	    }
	}
	$result .= "return unless " . join (' OR ',@matches) . ";\n";
    }
    return $result;
}

# turn a glob expression into a regexp
sub _glob_match {
    my $self = shift;
    my $term = shift;
    return unless $term =~ /(?:^|[^\\])[*?]/;
    $term =~ s/(^|[^\\])([+\[\]^{}\$|\(\).])/$1\\$2/g;
    $term =~ s/(^|[^\\])\*/$1.*/g;
    $term =~ s/(^|[^\\])\?/$1./g;
    return $term;
}

package Bio::DB::Sam::Fai;

sub open { shift->load(@_) }

sub seq {
    my $self = shift;
    my ($seqid,$start,$end) = @_;
    my $region = $seqid;
    $region   .= ":$start" if defined $start;
    $region   .= "-$end"   if defined $end;
    return $self->fetch($region)
}

package Bio::SeqFeature::Coverage;

use base 'Bio::SeqFeature::Lite';

sub coverage {
    my $self       = shift;
    my ($coverage) = $self->get_tag_values('coverage');
    return wantarray ? @$coverage : $coverage;
}

sub source {
    my $self = shift;
    my $type = $self->type;
    my ($base,$width) = split ':',$type;
    return $width;
}

sub method {
    my $self = shift;
    my $type = $self->type;
    my ($base,$width) = split ':',$type;
    return $base;
}

sub gff3_string {
    my $self = shift;
    my $gff3 = $self->SUPER::gff3_string;
    my $coverage = $self->escape(join(',',$self->coverage));
    $gff3 =~ s/coverage=[^;]+/coverage=$coverage/g;
    return $gff3;
}

package Bio::DB::Bam;

use File::Spec;
use Cwd;
use Carp 'croak';

sub index {
    my $self = shift;
    my $path = shift;
    my $autoindex = shift;

    return $self->index_open_in_safewd($path) if Bio::DB::Sam->is_remote($path);

    if ($autoindex) {
	$self->reindex($path) unless
	    -e "${path}.bai" && mtime($path) <= mtime("${path}.bai");
    }

    croak "No index file for $path; try opening file with -autoindex" unless -e "${path}.bai";
    return $self->index_open($path);
}

sub reindex {
    my $self = shift;
    my $path = shift;

    # if bam file is not sorted, then index_build will exit.
    # we spawn a shell to intercept this eventuality
    print STDERR "[bam_index_build] creating index for $path\n" if -t STDOUT;

    my $result = open my $fh,"-|";
    die "Couldn't fork $!" unless defined $result;

    if ($result == 0) { # in child
	# dup stderr to stdout so that we can intercept messages from library
	open STDERR,">&STDOUT";  
	$self->index_build($path);
	exit 0;
    }

    my $mesg = <$fh>;
    $mesg  ||= '';
    close $fh;
    if ($mesg =~ /not sorted/i) {
	print STDERR "[bam_index_build] sorting by coordinate...\n" if -t STDOUT;
	$self->sort_core(0,$path,"$path.sorted");
	rename "$path.sorted.bam",$path;
	$self->index_build($path);
    } elsif ($mesg) {
	die $mesg;
    }
}

# same as index_open(), but changes current wd to TMPDIR to accomodate
# the C library when it tries to download the index file from remote
# locations.
sub index_open_in_safewd {
    my $self = shift;
    my $dir    = getcwd;
    my $tmpdir = File::Spec->tmpdir;
    chdir($tmpdir);
    my $result = $self->index_open(@_);
    chdir $dir;
    $result;
}

sub mtime {
    my $path = shift;
    (stat($path))[9];
}


1;
__END__


=head1 EXAMPLES

For illustrative purposes only, here is an extremely stupid SNP caller
that tallies up bases that are q>20 and calls a SNP if there are at
least 4 non-N/non-indel bases at the position and at least 25% of them
are a non-reference base.

 my @SNPs;  # this will be list of SNPs
 my $snp_caller = sub {
	my ($seqid,$pos,$p) = @_;
	my $refbase = $sam->segment($seqid,$pos,$pos)->dna;
        my ($total,$different);
	for my $pileup (@$p) {
	    my $b     = $pileup->alignment;
            next if $pileup->indel;  # don't deal with these ;-)

            my $qbase  = substr($b->qseq,$pileup->qpos,1);
            next if $qbase =~ /[nN]/;

            my $qscore = $b->qscore->[$pileup->qpos];
            next unless $qscore > 25;

            $total++;
            $different++ if $refbase ne $qbase;
	}
        if ($total >= 4 && $different/$total >= 0.25) {
           push @SNPs,"$seqid:$pos";
        }
    };

 $sam->pileup('seq1',$snp_caller);
 print "Found SNPs: @SNPs\n";

=head1 GBrowse Compatibility

The Bio::DB::Sam interface can be used as a backend to GBrowse
(gmod.sourceforge.net/gbrowse). GBrowse can calculate and display
coverage graphs across large regions, alignment cartoons across
intermediate size regions, and detailed base-pair level alignments
across small regions.

Here is a typical configuration for a BAM database that contains
information from a shotgun genomic sequencing project. Some notes:

 * It is important to set "search options = none" in order to avoid
   GBrowse trying to scan through the BAM database to match read
   names. This is a time-consuming operation.

 * The callback to "bgcolor" renders pairs whose mates are unmapped in
   red.

 * The callback to "balloon hover" causes a balloon to pop up with the
   read name when the user hovers over each paired read. Otherwise the
   default behavior would be to provide information about the pair as
   a whole.

 * When the user zooms out to 1001 bp or greaterp, the track switches
   to a coverage graph.

 [bamtest:database]
 db_adaptor    = Bio::DB::Sam
 db_args       = -bam   /var/www/gbrowse2/databases/bamtest/ex1.bam
 search options= default

 [Pair]
 feature       = read_pair
 glyph         = segments
 database      = bamtest
 draw_target   = 1
 show_mismatch = 1
 bgcolor      = sub {
	     	 my $f = shift;
		 return $f->get_tag_values('M_UNMAPPED') ? 'red' : 'green';
	       }
 fgcolor       = green
 height        = 3
 label         = sub {shift->display_name}
 label density = 50
 bump          = fast
 connector     = dashed
 balloon hover = sub {
	      	    my $f     = shift;
		    return '' unless $f->type eq 'match';
		    return 'Read: '.$f->display_name.' : '.$f->flag_str;
                }
 key          = Read Pairs

 [Pair:1000]
 feature      = coverage:1001
 glyph        = wiggle_xyplot
 height       = 50
 min_score    = 0
 autoscale    = local

To show alignment data correctly when the user is zoomed in, you
should also provide a pointer to the FASTA file containing the
reference genome. In this case, modify the db_args line to read:

 db_args       = -bam   /var/www/gbrowse2/databases/bamtest/ex1.bam
                 -fasta /var/www/gbrowse2/databases/bamtest/ex1.fa

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::Bam::Alignment>, L<Bio::DB::Bam::Constants>

=head1 AUTHOR

Lincoln Stein E<lt>lincoln.stein@oicr.on.caE<gt>.
E<lt>lincoln.stein@bmail.comE<gt>

Copyright (c) 2009 Ontario Institute for Cancer Research.

This package and its accompanying libraries is free software; you can
redistribute it and/or modify it under the terms of the GPL (either
version 1, or at your option, any later version) or the Artistic
License 2.0.  Refer to LICENSE for the full license text. In addition,
please see DISCLAIMER.txt for disclaimers of warranty.

=cut

