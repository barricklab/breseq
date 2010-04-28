package Bio::DB::Bam::Alignment;

# $Id: Alignment.pm 22988 2010-04-02 17:17:34Z lstein $

=head1 NAME

Bio::DB::Bam::Alignment -- The SAM/BAM alignment object

=head1 SYNOPSIS

 use Bio::DB::Sam;

 my $sam = Bio::DB::Sam->new(-fasta=>"data/ex1.fa",
			     -bam  =>"data/ex1.bam");

 my @alignments = $sam->get_features_by_location(-seq_id => 'seq2',
                                                 -start  => 500,
                                                 -end    => 800);
 for my $a (@alignments) {
    my $seqid  = $a->seq_id;
    my $start  = $a->start;
    my $end    = $a->end;
    my $strand = $a->strand;
    my $ref_dna= $a->dna;

    my $query_start  = $a->query->start;
    my $query_end    = $a->query->end;
    my $query_strand = $a->query->strand;
    my $query_dna    = $a->query->dna;
   
    my $cigar     = $a->cigar_str;
    my @scores    = $a->qscore;     # per-base quality scores
    my $match_qual= $a->qual;       # quality of the match

    my $paired = $a->get_tag_values('PAIRED');
 }

=head1 DESCRIPTION

The Bio::DB::Bam::Alignment and Bio::DB::Bam::AlignWrapper classes
together represent an alignment between a sequence read (the "query")
and a reference sequence (the "target"). Bio::DB::Bam::Alignment
adheres strictly to the C-level BAM library's definition of a bam1_t*
and is used in the Bio::DB::Sam low-level API The latter adds
convenience methods that make it similar to a BioPerl Bio::SeqFeatureI
object. This manual page describes both.

=head1 High-level Bio::DB::Bam::Alignment methods

These methods are provided by Bio::DB::Bam::Alignment, and are
intended to be compatible with the Bio::SeqFeatureI interfaces. Note
that these objects are B<not> compatible with Bio::Align::AlignI, as
the BAM API is fundamentally incompatible with the BioPerl API for
alignments (the first deals with the alignment of a single read
against the reference sequence, while the second deals with a multiple
alignment).

Note that the high-level API return Bio::DB::Bam::AlignWrapper objects
B<except> in the case of the callback to the fast_pileup() method. In
this case only, the object returned by calling $pileup->b() is a
Bio::DB::Bam::Alignment object for performance reasons.

=over 4

=item $seq_id = $align->seq_id

Return the seq_id of the reference (target) sequence. This method is only
available in the Bio::DB::Bam::AlignWrapper extension.

=item $start = $align->start

Return the start of the alignment in 1-based reference sequence
coordinates.

=item $end = $align->end

Return the end of the alignment in 1-based reference sequence
coordinates.

=item $len = $align->length

Return the length of the alignment on the reference sequence.

=item $strand = $align->strand

Return the strand of the alignment as -1 for reversed, +1 for
forward. 

NOTE: In versions 1.00-1.06, this method always returned +1. As of
version 1.07, this behavior is fixed.

=item $mstrand = $align->mstrand

If the read has a mate pair, return the strand of the mate in the
format -1 or +1. 

=item $ref_dna        = $align->dna

Returns the B<reference> sequence's DNA across the aligned region. If
an MD tag is present in the alignment, it will be used preferentially
to reconstruct the reference sequence. Otherwise the reference DNA
access object passed to Bio::DB::Sam->new() will be used.

=item $ref_dna        = $align->seq

The B<reference> sequence's DNA as a Bio::PrimarySeqI object (useful
for passing to BioPerl functions and for calculating subsequences and
reverse complements).

=item $query = $align->query

This method returns a Bio::DB::Alignment::Query object that can be
used to retrieve information about the query sequence. The next few
entries show how to use this object.

=item $read_name = $align->query->name

The name of the read.

=item $q_start   = $align->query->start

This returns the start position of the query (read) sequence in
1-based coordinates. It acts via a transient Bio::DB::Bam::Query
object that is provided for Bio::Graphics compatibility (see
L<Bio::Graphics>).

=item $q_end     = $align->query->end

This returns the end position of the query sequence in 1-based
coordinates.

=item $q_len     = $align->query->length

Return the length of the alignment on the read.

=item $scores = $align->query->score

Return an array reference containing the unpacked quality scores for
each base of the query sequence. The length of this array reference
will be equal to the length of the read.

=item $read_dna = $align->query->dna

The read's DNA string.

=item $read_seq = $align->query->seq

The read's DNA as a Bio::PrimarySeqI object.

=item $target  = $align->target;

The target() method is similar to query(), except that it follows
Bio::AlignIO conventions for how to represent minus strand
alignments. The object returned has start(), end(), qscore(), dna()
and seq() methods, but for minus strand alignments the sequence will
be represented as it appears on the reverse strand, rather than on the
forward strand. This has the advantage of giving you the read as it
came off the machine, before being reverse complemented for use in the
SAM file.

=item $query   = $align->hit

The hit() method is identical to target() and returns information
about the read. It is present for compatibility with some of the
Bio::Graphics glyphs, which use hit() to represent the non-reference
sequence in aligned sequences.

=item $primary_id = $align->primary_id

This method synthesizes a unique ID for the alignment which can be
passed to $sam->get_feature_by_id() to retrieve the alignment at a
later date.

=item @tags = $align->get_all_tags

Return all tag names known to this alignment. This includes SAM flags
such as M_UNMAPPED, as well as auxiliary flags such as H0. The
behavior of this method depends on the value of -expand_flags when the
SAM object was created. If false (the default), then the standard SAM
flags will be concatenated together into a single string and stored in
a tag named 'FLAGS'. The format of this tag value is the list of one
or more flag constants separated by the "|" character, as in:
"PAIRED|MAP_PAIR|REVERSED|SECOND_MATE". If -expand_flags was true,
then each flag becomes its own named tag, such as "MAP_PAIR".

=item @values = $align->get_tag_values($tag)

Given a tag name, such as 'PAIRED' or 'H0', return its
value(s). -expand_flags must be true in order to use the standard SAM
flag constants as tags. Otherwise, they can be fetched by asking for
the "FLAGS" tag, or by using the low-level methods described below.

=item $is_true = $align->has_tag($tag)

Return true if the alignment has the indicated tag.

=item $string = $align->cigar_str

Return the CIGAR string for this alignment in conventional human
readable format (e.g. "M34D1M1").

=item $arrayref = $align->cigar_array

Return a reference to an array representing the CIGAR string. This is
an array of arrays, in which each subarray consists of a CIGAR
operation and a count. Example:

 [ ['M',34], ['D',1], ['M1',1] ]

=item ($ref,$matches,$query) = $align->padded_alignment

Return three strings that show the alignment between the reference
sequence (the target) and the query. It will look like this:

 $ref     AGTGCCTTTGTTCA-----ACCCCCTTGCAACAACC
 $matches ||||||||||||||     |||||||||||||||||
 $query   AGTGCCTTTGTTCACATAGACCCCCTTGCAACAACC 

=item $tag = $align->primary_tag

This is provided for Bio::SeqFeatureI compatibility. Return the string
"match".

=item $tag = $align->source_tag

This is provided for Bio::SeqFeatureI compatibility. Return the string
"sam/bam".

=item @parts = $align->get_SeqFeatures

Return subfeatures of this alignment. If you have fetched a
"read_pair" feature, this will be the two mate pair objects (both of
type Bio::DB::Bam::AlignWrapper). If you have -split_splices set to
true in the Bio::DB::Sam database, calling get_SeqFeatures() will
return the components of split alignments. See
L<Bio::DB::Sam/Bio::DB::Sam Constructor and basic accessors> for an
example of how to use this.

=back

=head1 Low-level Bio::DB::Bam::Alignment methods

These methods are available to objects of type Bio::DB::Bam::Alignment
as well as Bio::DB::Bam::AlignWrapper and closely mirror the native C
API.

=over 4

=item $align = Bio::DB::Bam::Alignment->new

Create a new, empty alignment object. This is usually only needed when
iterating through a TAM file using Bio::DB::Tam->read1().

=item $tid = $align->tid( [$new_tid] )

Return the target ID of the alignment. Optionally you may change the
tid by providing it as an argument (currently this is the only field
that you can change; the functionality was implemented as a proof of
principle).

=item $read_name = $align->qname

Returns the name of the read.

=item $pos = $align->pos

0-based leftmost coordinate of the aligned sequence on the reference
sequence.

=item $end = $align->calend

The 0-based rightmost coordinate of the aligned sequence on the
reference sequence after taking alignment gaps into account.

=item $len = $align->cigar2qlen

The length of the query sequence calculated from the CIGAR string.

=item $quality = $align->qual

The quality score for the alignment as a whole.

=item $flag = $align->flag

The bitwise flag field (see the SAM documentation).

=item $n_cigar = $align->n_cigar

Number of CIGAR operations in this alignment.

=item $length = $align->l_qseq

The length of the query sequence (the read).

=item $dna = $align->qseq

The actual DNA sequence of the query. As in the SAM file, reads that
are aligned to the minus strand of the reference are returned in
reverse complemented form.

=item $score_str = $align->_qscore

A packed binary string containing the quality scores for each base of
the read. It will be the same length as the DNA. You may unpack it
using unpack('C*',$score_str), or use the high-level qscore() method.

=item $score_arry = $align->qscore

=item @score_arry = $align->qscore

In a scalar context return an array reference containing the unpacked
quality scores for each base of the query sequence. In a list context
return a list of the scores. This array is in the same orientation as
the reference sequence.

=item $length = $align->isize

The calculated insert size for mapped paired reads.

=item $length = $align->l_aux

The length of the align "auxiliary" data.

=item $value = $align->aux_get("tag")

Given an auxiliary tag, such as "H0", return its value.

=item @keys  = $align->aux_keys

Return the list of auxiliary tags known to this alignment.

=item $data = $align->data

Return a packed string containing the alignment data (sequence,
quality scores and cigar string).

=item $length = $align->data_len

Return the current length of the alignment data.

=item $length = $align->m_data

Return the maximum length of the alignment data.

=item $is_paired = $align->paired

Return true if the aligned read is part of a mate/read pair
(regardless of whether the mate mapped).

=item $is_proper = $align->proper_pair

Return true if the aligned read is part of a mate/read pair and both
partners mapped to the reference sequence.

=item $is_unmapped = $align->unmapped

Return true if the read failed to align.

=item $mate_is_unmapped = $align->munmapped

Return true if the read's mate failed to align.

=item $reversed = $align->reversed

Return true if the aligned read was reverse complemented prior to
aligning.

=item $mate_reversed = $align->mreversed

Return true if the aligned read's mate was reverse complemented prior
to aligning.

=item $mstart  = $align->mate_start

For paired reads, return the start of the mate's alignment in
reference sequence coordinates.


=item $mend  = $align->mate_end

For paired reads, return the end position of the mate's alignment. in
reference sequence coordinates.

-item $len   = $align->mate_len

For mate-pairs, retrieve the length of the mate's alignment on the
reference sequence. 

=item $isize = $align->isize

For mate-pairs, return the computed insert size.

=item $arrayref = $align->cigar

This returns the CIGAR data in its native BAM format. You will receive
an arrayref in which each operation and count are packed together into
an 8-bit structure. To decode each element you must use the following
operations:

 use Bio::DB::Sam::Constants;
 my $c   = $align->cigar;
 my $op  = $c->[0] & BAM_CIGAR_MASK;
 my $len = $c->[0] >> BAM_CIGAR_SHIFT;

=back

=cut

use strict;
use warnings;
use Bio::DB::Bam::Query;
use Bio::DB::Bam::Target;
use Bio::DB::Sam::Constants;

sub each_tag_value { shift->get_tag_values(@_)  }

sub get_tag_values {
    my $self = shift;
    my $tag  = shift;
    defined $tag or return;
    if (my $mask = RFLAGS->{uc $tag}) {  # special tag
	# to avoid warnings when making numeric comps
	return ($self->flag & $mask) == 0 ? 0 : 1; 
    } elsif ($tag eq 'FLAGS') {
	$self->flag_str;
    } else {
	$self->aux_get($tag);
    }
}

sub has_tag {
    my $self = shift;
    my $tag  = shift;
    defined $tag or return;
    if (my $mask = RFLAGS->{uc $tag}) {  # special tag
	return 1;
    } elsif ($tag eq 'FLAGS') {
	return 1;
    } else {
	my %keys = map {$_=>1} $self->aux_keys;
	return exists $keys{uc $tag};
    }
}

sub get_all_tags {
    my $self      = shift;
    my @aux_tags  = $self->aux_keys;
    my @flag_tags = keys %{RFLAGS()};
    return (@aux_tags,@flag_tags);
}


sub start {
    my $self = shift;
    return if $self->unmapped;
    return $self->pos+1;
}

sub end {
    my $self = shift;
    return if $self->unmapped;
    return $self->calend;
}

sub stop { shift->end }

# in SAM format, alignment is always to the forward strand
sub strand {
    my $self     = shift;
    return $self->reversed ? -1 : 1;
#     return 1;
}

sub abs_strand { shift->strand }

sub mstrand {
    my $self     = shift;
    return $self->mreversed ? -1 : 1;
}

sub display_name {
    return shift->qname;
}

sub qscore {
    my $self   = shift;
    my $scores = $self->_qscore;
    my @scores  = unpack('C*',$scores);
    return wantarray ? @scores : \@scores;
}

sub primary_id {
    my $self = shift;
    return join ';',
    map {s/;/%3B/g; $_}
    ($self->display_name,
     $self->tid,
     $self->start,
     $self->end,
     $self->strand);
}
sub cigar_str {
    my $self   = shift;
    my $cigar  = $self->cigar;
    my $result = '';
    for my $c (@$cigar) {
	my $op     = $c & BAM_CIGAR_MASK;
	my $l      = $c >> BAM_CIGAR_SHIFT();
	my $symbol = CIGAR_SYMBOLS()->[$op];
	$result .= "${symbol}${l}";
    }
    return $result;
}

sub cigar_array {
    my $self   = shift;
    my $cigar  = $self->cigar;
    my @result;
    for my $c (@$cigar) {
	my $op     = $c & BAM_CIGAR_MASK();
	my $l      = $c >> BAM_CIGAR_SHIFT();
	my $symbol = CIGAR_SYMBOLS()->[$op];
	push @result,[$symbol,$l];
    }
    return \@result;

}

sub flag_str {
    my $self  = shift;
    my $flag  = $self->flag;
    my $flags = FLAGS;
    return join '|',map {$flags->{$_}}
			 grep {$flag & $_}
			 sort {$a<=>$b}
			 keys %{$flags};
}

sub length {
    my $self = shift;
    my $end   = $self->end   || 0;
    my $start = $self->start || 0;
    return $end-$start+1;
}

sub mate_start {
    shift->mpos+1;
}

sub mate_len {
    my $self    = shift;
    my $ins_len = $self->isize or return;
    my $len     = $self->length;

    my $adjust = 0;
    my @cigar   = $self->cigar_str =~ /(\w)(\d+)/g;
    while (@cigar) {
	my ($op,$len) = splice(@cigar,0,2);
	$adjust += $len if $op eq 'I';
	$adjust -= $len if $op eq 'D';
    }

    return $adjust + $ins_len + ($self->start-$self->mate_start) if $ins_len > 0;
    return $adjust + $self->mate_start-($self->start+$ins_len)   if $ins_len < 0;
}

sub mate_end {
    my $self = shift;
    return unless $self->mate_len;
    return $self->mate_start+$self->mate_len-1;
}

sub query {
    my $self = shift;
    return Bio::DB::Bam::Query->new($self);
}

sub get_SeqFeatures { return; }

# Target is the same as Query, but with meaning of start() and end() reversed
# for compatibility with Bio::DB::GFF and its ilk. Please use Query if you can!
sub target {
    my $self = shift;
    return Bio::DB::Bam::Target->new($self);
}

sub primary_tag { 'match'   }
sub source_tag  { 'sam/bam' }

sub hit { shift->target(@_); }


1;

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::Sam>, L<Bio::DB::Bam::Constants>

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

