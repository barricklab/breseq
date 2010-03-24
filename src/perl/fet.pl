#!/usr/bin/perl -w

use FindBin;
use lib $FindBin::Bin;
use Statistics::FishersExactTest;
my $fisher_strand_p_value = Statistics::FishersExactTest::fishers_exact(100,200, 5, 17,0);
#$fisher_strand_p_value = sprintf "%.1e", $fisher_strand_p_value; #round immediately
print "$fisher_strand_p_value\n";
