package Bio::DB::Bam::Query;

# $Id: Query.pm 22675 2010-02-08 21:57:02Z lstein $

=head1 NAME

Bio::DB::Bam::Query -- Object representing the query portion of a BAM/SAM alignment

=head1 SYNOPSIS

Given an alignment retrieved from a Bio::DB::Sam database,

 my $query = $alignment->query;

 my $name   = $query->display_name;
 my $start  = $query->start;
 my $end    = $query->end;
 my $dna    = $query->dna;    # dna string
 my $seq    = $query->seq;    # Bio::PrimarySeq object
 my @scores = $query->qscore; # quality score

=head1 DESCRIPTION

This is a simple Bio::SeqFeatureI object that represents the query
part of a SAM alignment.

=head2 Methods

=over 4

=cut

use strict;
use Bio::DB::Sam;
use Bio::DB::Sam::Constants qw(CIGAR_SYMBOLS BAM_CREF_SKIP BAM_CSOFT_CLIP BAM_CHARD_CLIP);

use constant CIGAR_SKIP      => {CIGAR_SYMBOLS->[BAM_CREF_SKIP]  => 1,
 				 CIGAR_SYMBOLS->[BAM_CSOFT_CLIP] => 1,
 				 CIGAR_SYMBOLS->[BAM_CHARD_CLIP] => 1};


sub new {
    my $self      = shift;
    my $alignment = shift;
    bless \$alignment,ref $self || $self;
}

=item $seqid = $query->seq_id

The name of the read.

=cut

sub seq_id {
    my $self = shift;
    $$self->qname;
}

=item $name = $query->name

The read name (same as seq_id in this case).

=cut

sub name {
    my $self = shift;
    $$self->qname;
}

=item $name = $query->display_name

The read display_name (same as seq_id in this case).

=cut

sub display_name {shift->name}

=item $tag = $query->primary_tag

The string "match".

=cut

sub primary_tag { ${shift()}->primary_tag }


=item $tag = $query->source_tag

The string "sam/bam".

=cut

sub source_tag  { ${shift()}->source_tag  }

=item $start = $query->start

The start of the match in read coordinates.

=cut

sub start {
    my $self = shift;
    return $self->low;
}

=item $end = $query->end

The end of the match in read coordinates;

=cut

sub end {
    my $self = shift;
    return $self->high;
}

sub low {
    my $self       = shift;
    my $cigar_arry = $$self->cigar_array;
    my $start      = 1;
    for my $c (@$cigar_arry) {
	last unless CIGAR_SKIP->{$c->[0]};
	$start += $c->[1];
    }
    $start;
}

sub high {
    my $self      = shift;
    my $len       = $$self->cigar2qlen;
    my $cigar_arry = $$self->cigar_array;

    # alignment stops at first non-clip CIGAR position
    my $i = $len - 1;
    for my $c (reverse @$cigar_arry) {
	last unless CIGAR_SKIP->{$c->[0]};
	$len -= $c->[1];
    }
    return $len;
}

=item $len = $query->length

The length of the read.

=cut

sub length {
    my $self = shift;
    $self->high-$self->low+1;
#    $$self->cigar2qlen;
}

=item $seq = $query->seq

A Bio::PrimarySeq representing the read sequence in REFERENCE
orientation.

=cut

sub seq { 
    my $self = shift;
    my $dna  = $self->dna;
    return Bio::PrimarySeq->new(-seq => $dna,
				-id  => $$self->qname);
}

=item $scores = $query->qscore

The read quality scores. In a list context, a list of integers equal
in length to the read sequence length. In a scalar context, an array
ref. The qscores are in REFERENCE sequence orientation.

=cut

sub qscore {
    my $self = shift;
    my @qscore = $$self->qscore;
    return wantarray ? @qscore : \@qscore;
}

=item $dna = $query->dna

The DNA string in reference sequence orientation.

=cut

sub dna {
    my $self = shift;
    return $$self->qseq;
}

=item $strand = $query->strand

If the query was reversed to align it, -1. Otherwise +1.

=cut

sub strand { 
    my $self = shift;
    return $$self->reversed ? -1 : 1;
}

=item $seq = $query->subseq($start,$end)

Return a Bio::PrimarySeq object representing the requested subsequence
on the read.

=cut

sub subseq {
    my $self = shift;
    my ($start,$end) = @_;
    $start = 1 if $start < 1;
    $end   = $self->high if $end > $self->high;
    ($end,$start) = ($start,$end) if $start > $end;
    return Bio::PrimarySeq->new(-seq=>substr($self->dna,
					     $start-1,
					     $end-$start+1)
				);
}


1;

=back

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::Sam>, L<Bio::DB::Bam::Alignment>, L<Bio::DB::Bam::Constants>

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
