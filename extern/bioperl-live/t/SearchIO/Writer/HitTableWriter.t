# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_HitTableWriter.t 14995 2008-11-16 06:20:00Z cjfields $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 8);
    
    use_ok('Bio::SearchIO');
    use_ok('Bio::SearchIO::Writer::HitTableWriter');
}

my ($searchio, $result, $hit, $hsp);

$searchio = Bio::SearchIO->new('-format' => 'blast',
    '-file'   => test_input_file('HUMBETGLOA.tblastx'));

$result = $searchio->next_result;

isa_ok($result,'Bio::Search::Result::ResultI');
$hit = $result->next_hit;
is($hit->accession, 'AE000479');
is($hit->bits, 33.6);
$hsp = $hit->next_hsp;
is($hit->hsp->bits,$hsp->bits);
isa_ok($hsp->get_aln,'Bio::Align::AlignI');

my $writer = Bio::SearchIO::Writer::HitTableWriter->new( 
    -columns => [qw(query_name
                    query_length
                    hit_name
                    hit_length
                    bits
                    score
                    frac_identical_query
                    expect
                    )]  );

my $outfile = test_output_file();
my $out = Bio::SearchIO->new(-writer => $writer,
             -file   => ">$outfile");
$out->write_result($result, 1);
ok(-s $outfile);

# tests checking file output?
