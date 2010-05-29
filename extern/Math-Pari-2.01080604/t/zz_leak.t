#!/usr/bin/perl -w

use strict;
my $test = 0;
my $c;
$| = 1;
print "1..",&last,"\n";

sub test {
  $test++; 
  my $mess = shift() ? '' : 'not ';
  my $what = shift;
  print "${mess}ok $test # $what\n";
  return !$mess;
}

use Math::Pari;

test(1, "before the test");		# 1
my $x;

for my $by ([1, "easy case: 1"],
	    [0, "harder: 0"]) {
  my $t = PARI($by->[0]);
  for (1..100) {
    for my $i ((1) x 1e4) {		# Was leaking with 0, no leaks with 1
      $x = 256*$t;
    }
    #print '.' unless $_ % 100;		# Give chance to free temporaries
    $c++;
  }
  test(1, "after: $by->[1]");	# 2, 3: after the test
}

sub last {3}
