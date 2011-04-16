package Bio::DB::Sam::Constants;

# $Id: Constants.pm 22410 2009-12-15 16:23:24Z lstein $

=head1 NAME

Bio::DB::Sam::Constants -- Constants for use with SAM/BAM

=head1 SYNOPSIS

 use Bio::DB::Sam::Constants;
 my $pad_flag = BAM_CPAD;

=head1 DESCRIPTION

This module exports several constants for use with the SAM/BAM
module. See the SAM documentation for their interpretation.

=over 4

=item Cigar operations

  BAM_CIGAR_SHIFT
  BAM_CIGAR_MASK
  BAM_CMATCH
  BAM_CINS
  BAM_CDEL
  BAM_CREF_SKIP
  BAM_CSOFT_CLIP
  BAM_CHARD_CLIP
  BAM_CPAD

=item FLAGS

A hashref that maps flag values to human-readable names. For example:

 FLAGS->{0x0008} == 'M_UNMAPPED'

=item RFLAGS

The reverse of FLAGS:

 FLAGS->{M_UNMAPPED} == 0x0008

=back

=head1 SEE ALSO

L<Bio::Perl>, L<Bio::DB::Sam>, L<Bio::DB::Bam::Alignment>

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


use strict;

require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(CIGAR_SYMBOLS BAM_CIGAR_SHIFT BAM_CIGAR_MASK
                 BAM_CMATCH BAM_CINS BAM_CDEL BAM_CREF_SKIP
                 BAM_CSOFT_CLIP BAM_CHARD_CLIP BAM_CPAD FLAGS RFLAGS);
our @EXPORT_OK = @EXPORT;

use constant CIGAR_SYMBOLS   => [qw(M I D N S H P)];
use constant BAM_CIGAR_SHIFT => 4;
use constant BAM_CIGAR_MASK  => (1 << BAM_CIGAR_SHIFT) - 1;
use constant BAM_CMATCH      => 0;
use constant BAM_CINS        => 1;
use constant BAM_CDEL        => 2;
use constant BAM_CREF_SKIP   => 3;
use constant BAM_CSOFT_CLIP  => 4;
use constant BAM_CHARD_CLIP  => 5;
use constant BAM_CPAD        => 6;

use constant FLAGS => {
    0x0001 => 'PAIRED',
    0x0002 => 'MAP_PAIR',
    0x0004 => 'UNMAPPED',
    0x0008 => 'M_UNMAPPED',
    0x0010 => 'REVERSED',
    0x0020 => 'M_REVERSED',
    0x0040 => 'FIRST_MATE',
    0x0080 => 'SECOND_MATE',
    0x0100 => 'NOT_PRIMARY',
    0x0200 => 'QC_FAILED',
    0x0400 => 'DUPLICATE'
};
use constant RFLAGS => {reverse %{FLAGS()}};

1;
