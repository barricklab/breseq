package Bio::DB::Bam::Pileup;

# $Id: Pileup.pm 22410 2009-12-15 16:23:24Z lstein $

# documentation only

1;

=head1 NAME

Bio::DB::Bam::Pileup -- Object passed to pileup() callback

=head1 SYNOPSIS

See L<Bio::DB::Sam/The generic fetch() and pileup() methods> for how
this object is passed to pileup callbacks.

=head1 DESCRIPTION

A Bio::DB::Bam::Pileup object (or a Bio::DB::Bam::PileupWrapper
object) is passed to the callback passed to the Bio::DB::Sam->pileup()
method for each column in a sequence alignment. The only difference
between the two is that the latter returns the more convenient
Bio::DB::Bam::AlignWrapper objects in response to the alignment()
method, at the cost of some performance loss.

=head2 Methods

=over 4

=item $alignment = $pileup->alignment

Return the Bio::DB::Bam::Alignment or Bio::DB::Bam::AlignWrapper
object representing the aligned read.

=item $alignment = $pileup->b

This method is an alias for alignment(). It is available for
compatibility with the C API.

=item $qpos = $pileup->qpos

Return the position of this aligned column in read coordinates, using
zero-based coordinates.

=item $pos  = $pileup->pos

Return the position of this aligned column in read coordinates, using
1-based coordinates.

=item $indel = $pileup->indel

If this column is an indel, return a positive integer for an insertion
relative to the reference, a negative integer for a deletion relative
to the reference, or 0 for no indel at this column.

=item $is_del = $pileup->is_del

True if the base on the padded read is a deletion.

=item $level  = $pileup->level

If pileup() or fast_pileup() was invoked with the "keep_level" flag,
then this method will return a positive integer indicating the level
of the read in a printed multiple alignment.

=item $pileup->is_head

=item $pileup->is_tail

These fields are defined in bam.h but their interpretation is obscure.

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
