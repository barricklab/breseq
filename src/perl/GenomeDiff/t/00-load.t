#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'CAMP' ) || print "Bail out!
";
}

diag( "Testing CAMP $CAMP::VERSION, Perl $], $^X" );
