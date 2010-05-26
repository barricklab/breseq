#-*-Perl-*-

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use ExtUtils::MakeMaker;
use File::Temp qw(tempfile);
use FindBin '$Bin';
use constant TEST_COUNT => 101;

use lib "$Bin/../lib","$Bin/../blib/lib","$Bin/../blib/arch";

BEGIN {
  # to handle systems with no installed Test module
  # we include the t dir (where a copy of Test.pm is located)
  # as a fallback
  eval { require Test; };
  if( $@ ) {
    use lib 't';
  }
  use Test;
  plan test => TEST_COUNT;
}

use Bio::DB::Sam;

# low level tests (defined in lib/Bio/DB/Sam.xs) 
{
    my $bamfile = "$Bin/data/ex1.bam";
    my $bam     = Bio::DB::Bam->open($bamfile);
    ok($bam);

    my $header  = $bam->header;
    my $targets = $header->n_targets;
    ok($targets,2);

    my $target_names = $header->target_name;
    ok($target_names);
    ok(scalar @$target_names,2);
    ok($target_names->[0],'seq1');
    
    my $target_lens = $header->target_len;
    ok($target_lens);
    ok(scalar @$target_lens,2);
    ok($target_lens->[0],1575);
    
    my $text = $header->text;
    ok(length $text > 0);
    
    my $fai  = Bio::DB::Sam::Fai->open("$Bin/data/ex1.fa");
    my $seq  = $fai->fetch('seq2:51-1000');
    ok(length $seq,950);

    my $count;
    while (my $b = $bam->read1) {
	$count++;
    }
    ok($count,3307);
    
    my @result = $header->parse_region('seq2:51-1000');
    ok($result[0],1);
    ok($result[1],50);
    @result    = $header->parse_region('seq_invalid:51-1000');
    ok(scalar @result,0);
    
    my $index = Bio::DB::Bam->index($bamfile,1);
    ok($index);

    my @a;
    my $print_region = sub {
	my($alignment,$data) = @_;
	push @a,$alignment;
	return;
    };

    $index->fetch($bam,$header->parse_region('seq2'),$print_region,"foobar");
    ok(scalar @a > 1);

    my %matches;
    my $fetch_back = sub {
	my ($tid,$pos,$p) = @_;
	my $p1 = $pos+1;
	my $r = $fai->fetch($header->target_name->[$tid].":$p1-$p1");
	for my $pileup (@$p) {
	    my $b    = $pileup->b;
	    my $qpos = $pileup->qpos;
	    my $base = $pileup->indel == 0 ? substr($b->qseq,$qpos,1)
                      :$pileup->indel >  0 ? '*'
                      : '-';
	    $matches{matched}++ if $r eq $base;
	    $matches{total}++;
	}
    };
    
    $index->pileup($bam,$header->parse_region('seq2:1-100'),$fetch_back);
    ok($matches{matched}/$matches{total} > 0.99);

    # try to get coverage
    my $coverage = $index->coverage($bam,
				   $header->parse_region('seq2'),
				    100);
    ok(scalar @$coverage,100);
    my @c = sort {$a<=>$b} @$coverage;
    ok($c[0]  >= 0);
    ok($c[-1] < 1000);

    # try reading from a TAM (text sam) file
    my $sam = Bio::DB::Tam->open("$Bin/data/ex1.sam.gz");
    ok($sam);
    my $align = Bio::DB::Bam::Alignment->new();
    ok($align);

    # quench annoying stderr message from library here
    open my $saverr,">&STDERR";
    open STDERR,">/dev/null";
    my $head  = Bio::DB::Tam->header_read2("$Bin/data/ex1.fa.fai");
    open STDERR,">&",$saverr;
    ok($head);

    my $result = $sam->read1($head,$align);
    ok($result>0);
    ok($align->qseq,'CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCG');
    ok($align->start,1);
    ok($sam->read1($head,$align)>0);
    ok($align->start,3);
    ok($header->target_name->[$align->tid],'seq1');

    # test ability to write a BAM file
    my (undef,$filename) = tempfile('BAM_TEST_XXXXXXX',UNLINK=>1);
    $sam = Bio::DB::Tam->open("$Bin/data/ex1.sam.gz");
    $bam = Bio::DB::Bam->open($filename,'w');
    ok($bam);
    ok($bam->header_write($head),0);
    $count = 0;
    while ($sam->read1($head,$align) > 0) {
	$count++;
	$bam->write1($align);
    }
    ok($count,3307);
    undef $bam;

    $bam     = Bio::DB::Bam->open($filename);
    ok($bam);

    $header  = $bam->header;
    $targets = $header->n_targets;
    ok($targets,2);

    $target_names = $header->target_name;
    ok($target_names);
    ok(scalar @$target_names,2);
    ok($target_names->[0],'seq1');
    
    $target_lens = $header->target_len;
    ok($target_lens);
    ok(scalar @$target_lens,2);
    ok($target_lens->[0],1575);

    # try removing and regenerating index
    unlink "$Bin/data/ex1.bam.bai";
    ok(Bio::DB::Bam->index($bamfile,1));
    ok(-e "$Bin/data/ex1.bam.bai");

    Bio::DB::Bam->sort_core(1,"$Bin/data/ex1.bam","$Bin/data/ex1.sorted");
    ok(-e "$Bin/data/ex1.sorted.bam");
    ok(Bio::DB::Bam->index("$Bin/data/ex1.sorted.bam",1));
    ok(-e "$Bin/data/ex1.sorted.bam.bai");
    unlink ("$Bin/data/ex1.sorted.bam","$Bin/data/ex1.sorted.bam.bai");
}

# high level tests (defined in lib/Bio/DB/Sam.pm)
{
    my $sam = Bio::DB::Sam->new(-fasta=>"$Bin/data/ex1.fa",
			        -bam  =>"$Bin/data/ex1.bam",
				-expand_flags => 1,
				-autoindex => 1,
	);
    ok($sam);
    ok($sam->n_targets,2);

    ok($sam->length('seq1'),1575);
    ok(join $sam->seq_ids,'seq1 seq2');

    my $seg = $sam->segment('seq1');
    ok($seg);
    ok($seg->length,1575);
    my $seq = $seg->seq;
    ok($seq->isa('Bio::PrimarySeq'));
    ok(length $seq->seq,1575);

    my $fh = $sam->tam_fh;
    ok($fh);
    my $samline = <$fh>;
    ok($samline =~ /^B7_591:4:96:693:509/);
    $fh->close;

    my $dummy = eval {Bio::DB::Sam->new(-fasta=>"invalid_path.txt",
					-bam  =>"invalid_path.txt")};
    ok($dummy,undef);
    ok($@ =~ /does not exist/);
    
    my @alignments = 
	$sam->get_features_by_location(
	    -seq_id => 'seq2',
	    -start  => 500,
	    -end    => 800
	);
    ok(scalar @alignments,442);
    ok($alignments[0]->seq_id,'seq2');

    ok(scalar @{$alignments[0]->qscore},length $alignments[0]->dna);

    my @keys = $alignments[0]->get_all_tags;
    ok(scalar @keys,17);
    ok($alignments[0]->get_tag_values('MF'),18);

    my %att  = $alignments[0]->attributes;
    ok(scalar(keys %att),17);
    ok($alignments[0]->cigar_str,'M35');

    $sam->expand_flags(0);
    @keys = $alignments[0]->get_all_tags;
    ok(scalar @keys,7);

    my $query = $alignments[0]->query;
    ok($query);
    ok($query->start,1);
    ok($query->end,35);
    ok($query->length,35);
    ok($query->dna,$alignments[0]->dna);
    ok($alignments[0]->strand,-1);
    ok($query->strand,-1);

    my $target = $alignments[0]->target;
    ok($target);
    ok($target->start,35);
    ok($target->end,1);
    ok($target->length,35);
    ok($target->dna,reversec($alignments[0]->dna));

    ok($alignments[0]->get_tag_values('FLAGS'),$alignments[0]->flag_str);

    my @pads = $alignments[0]->padded_alignment;
    ok(@pads,3);
    ok($pads[0],$pads[2]);
    ok($pads[1]=~tr/|/|/,length($pads[0]));

    my @f = $sam->features(-name=>'EAS114_45:2:1:1140:1206');
    ok(scalar @f,2);

    @f    = $sam->features(-filter=> sub {
	                       my $a = shift;
			       return 1 if $a->display_name eq 'EAS114_45:2:1:1140:1206';
			   });
    ok(scalar @f,2);

    @f = $sam->features(-seq_id  => 'seq2',
			-filter  => sub {
			    my $a = shift;
			    return 1 if $a->display_name =~ /^EAS114/;
			});
    ok(scalar @f,306);
    @f = $sam->features(-filter  => sub {
			    my $a = shift;
			    return 1 if $a->display_name =~ /^EAS114/;
			});
    ok(scalar @f,534);

    @f   = $sam->features(-name=>'EAS114_45:2:1:1140:1206',
			  -tags=>{FIRST_MATE=>1});
    ok(scalar @f,1);

    # try iteration
    my $i = $sam->get_seq_stream(
	-seq_id => 'seq2',
	-start  => 500,
	-end    => 800
	);
    ok($i);
    my $count = 0;
    while ($i->next_seq) { $count++ }
    ok($count,442);

    # try tam fh retrieval
    $fh = $sam->get_seq_fh(
	-seq_id => 'seq2',
	-start  => 500,
	-end    => 800,
	);
    $count = 0;
    $count++ while <$fh>;
    ok($count,442);
    $fh->close;

    $i = $sam->get_seq_stream(); # all features!
    ok($i);
    $count = 0;
    while ($i->next_seq) { $count++ }
    ok($count,3307);

    $i = $sam->get_seq_stream(-max_features=>200,-seq_id=>'seq1');
    ok ($i);
    $count = 0;
    while ($i->next_seq) { $count++ }
    ok($count,200);
    ok($sam->last_feature_count,1482);

    # try the read_pair aggregation
    my @pairs = $sam->features(-type=>'read_pair',
			       -seq_id => 'seq2');
    ok(scalar @pairs,939);

    # try coverage
    my @coverage = $sam->features(-type   => 'coverage',
				  -seq_id => 'seq2');
    ok(scalar @coverage,1);
    my ($c) = $coverage[0]->get_tag_values('coverage');
    ok($c);
    ok($c->[0],3);
    ok($c->[1],4);
    ok($coverage[0]->type,"coverage:1584");

    # test high level API version of pileup
    my %matches;
    my $fetch_back = sub {
	my ($seqid,$pos,$p) = @_;
	my $r = $sam->segment($seqid,$pos,$pos)->dna;
	for my $pileup (@$p) {
	    my $a    = $pileup->alignment;
	    my $qpos = $pileup->qpos;
	    my $dna  = $a->query->dna;
	    my $base = $pileup->indel == 0 ? substr($dna,$qpos,1)
                      :$pileup->indel >  0 ? '*'
                      : '-';
	    $matches{matched}++ if $r eq $base;
	    $matches{total}++;
	}
    };
    
    $sam->pileup('seq2:1-100',$fetch_back);
    ok($matches{matched}/$matches{total} > 0.99);

    exit 0;

# this is not a unit test, but a piece of demo to show cigar string
# processing
    for my $a (@alignments) {
	warn $a->display_name,' :: ',$a->flag_str,' :: ',
	$a->start,'..',$a->end,' ',' mpos=',$a->mate_start,' ins=',$a->isize,
	' qlen=',$a->cigar2qlen,
	' [',$a->strand,'/',$a->mstrand,']',
	' ',
	$a->cigar_str,
	' ',
	$a->mate_start,'..',$a->mate_end,
	"\n";
    }
}

sub reversec {
    my $dna = shift;
    $dna    =~ tr/gatcGATC/ctagCTAG/;
    return scalar reverse $dna;
}
