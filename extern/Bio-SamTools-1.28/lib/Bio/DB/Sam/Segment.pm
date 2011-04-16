package Bio::DB::Sam::Segment;
use strict;

sub new {
    my $class                   = shift;
    my ($db,$seqid,$start,$end) = @_;
    return bless {
	db     => $db,
	seqid  => $seqid,
	start  => $start,
	end    => $end},ref $class || $class;
}

sub db       { shift->{db}    };

sub features {
    my $self = shift;
    my $db   = $self->db;
    my @args;

    if (@_ && $_[0] !~ /^-/) { # type list
	@args = (-types=>\@_);
    } else {
	@args = @_;
    }
    return $db->features(-seq_id => $self->seq_id,
			 -start  => $self->start,
			 -end    => $self->end,
			 @args);
}

#required by GBrowse1
sub get_feature_stream {
    my $self = shift;
    my @args = @_;

    my @features = $self->features(@args);

    return Bio::DB::Sam::Segment::Iterator->new(\@features);
}

# required by api
sub seq_id   { shift->{seqid} };
# required by GBrowse1
*ref = *abs_ref = *sourceseq = \&seq_id;
# required by api
sub start    { shift->{start} };
# required by api
sub end      { shift->{end}   };
# required by api
sub strand   { 0              };
# required by api
sub length   { 
    my $self = shift;
    return $self->end - $self->start + 1;
}
# required by api
sub seq      {
    my $self   = shift;
    return Bio::PrimarySeq->new(-seq => $self->dna,
				-id  => $self->seq_id);
}

sub primary_id {
    shift->seq_id;
}
# required by api
sub primary_tag {
    my $self = shift;
    return 'region';
}
sub dna {
    my $self = shift;
    my $db     = $self->db;
    $db->seq($self->seq_id,$self->start,$self->end);
}
# required by api
sub source_tag { return 'sam/bam' }
# required by api
sub name    { shift->seq_id }
# required by api
sub display_name    { shift->seq_id }
# required by api
sub factory { shift->db  }
# required by api
sub get_SeqFeatures { return;   }
# required by api
sub method { shift->primary_tag }
# required by GBrowse1
sub type { shift->primary_tag }
# required by api
sub get_tag_values {  return; }
# required by api
sub score { return;  }
# required by api
sub class { 'sequence'  }



package Bio::DB::Sam::Segment::Iterator;

sub new {
  my $package  = shift;
  my $features = shift;
  return bless $features,$package;
}

sub next_seq {
  my $self = shift;
  return unless @$self;
    my $next_feature = shift @$self;
  return $next_feature;
}

1;
