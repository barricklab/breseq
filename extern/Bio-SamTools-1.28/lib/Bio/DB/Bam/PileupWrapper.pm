package Bio::DB::Bam::PileupWrapper;
#$Id: PileupWrapper.pm 22410 2009-12-15 16:23:24Z lstein $

=head1 NAME

Bio::DB::Bam::PileupWrapper -- Add high-level methods to Bio::DB::Bam::Pileup

=head1 SYNOPSIS

See L<Bio::DB::Sam/The generic fetch() and pileup() methods> for usage of the pileup() method.

=head1 DESCRIPTION

See L<Bio::DB::Bam::Pileup> for documentation of this object's
methods. This class is used by the high-level API to return
Bio::DB::Bam::AlignWrapper objects from the call to alignment() rather
than Bio::DB::Bam::Alignment.

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

use strict;
use Bio::DB::Bam::AlignWrapper;

our $AUTOLOAD;
use Carp 'croak';

sub new {
    my $package = shift;
    my ($align,$sam) = @_;
    return bless {sam    => $sam,
		  pileup => $align},ref $package || $package;

}

sub AUTOLOAD {
  my($pack,$func_name) = $AUTOLOAD=~/(.+)::([^:]+)$/;
  return if $func_name eq 'DESTROY';

  no strict 'refs';
  $_[0] or die "autoload called for non-object symbol $func_name";
  croak qq(Can't locate object method "$func_name" via package "$pack")
      unless $_[0]->{pileup}->can($func_name);

  *{"${pack}::${func_name}"} = sub { shift->{pileup}->$func_name(@_) };

  shift->$func_name(@_);
}

sub can {
    my $self = shift;
    return 1 if $self->SUPER::can(@_);
    return $self->{pileup}->can(@_);
}

sub alignment {
    my $self = shift;
    return Bio::DB::Bam::AlignWrapper->new($self->{pileup}->b,$self->{sam});
}

1;

