package Bio::DB::Bam::FetchIterator;
use strict;

sub new {
    my $self = shift;
    my $list = shift;
    my $total= shift;
    $total ||= @$list;
    return bless {list=>$list,
		  total=>$total},ref $self || $self;
}

sub next_seq {
    my $self = shift;
    return shift @{$self->{list}};
}

sub total {
    shift->{total};
}

1;
