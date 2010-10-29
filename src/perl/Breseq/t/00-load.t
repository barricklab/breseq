#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Bio::Breseq' ) || print "Bail out!
";
}

diag( "Testing Bio::Breseq $Bio::Breseq::VERSION, Perl $], $^X" );
