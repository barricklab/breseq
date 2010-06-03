package Bio::DB::Bam::ReadIterator;

use strict;

sub new {
    my $self = shift;
    my ($sam,$bam,$filter) = @_;
    return bless {sam   => $sam,
		  bam   => $bam,
		  filter=> $filter},ref $self || $self;
}
sub next_seq {
    my $self = shift;
    while (my $b = $self->{bam}->read1) {
	return Bio::DB::Bam::AlignWrapper->new($b,$self->{sam})
	    if $self->{filter}->($b);
    }
    return;
}

1;
