# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..9\n"; }
END {print "not ok 1\n" unless $loaded;}
use Statistics::Distributions;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

if (Statistics::Distributions::chisqrdistr (2,.05) == 5.9915) { print "ok 2\n"; }
else { print "not ok 2\n"; }
if (Statistics::Distributions::udistr (.05) == 1.6449) { print "ok 3\n"; }
else { print "not ok 3\n"; }
if (Statistics::Distributions::tdistr (1,.005) == 63.657) { print "ok 4\n"; }
else { print "not ok 4\n"; }
if (Statistics::Distributions::fdistr (1,3,.01) == 34.116) { print "ok 5\n"; }
else { print "not ok 5\n"; }
if (Statistics::Distributions::uprob (-.85) == .80234) { print "ok 6\n"; }
else { print "not ok 6\n"; }
if (Statistics::Distributions::chisqrprob (3,6.25) == .10006) { print "ok 7\n"; }
else { print "not ok 7\n"; }
if (Statistics::Distributions::tprob (3,6.251) == .0041301) { print "ok 8\n"; }
else { print "not ok 8\n"; }
if (Statistics::Distributions::fprob (6,6,.625) == 0.70879) { print "ok 9\n"; }
else { print "not ok 9\n"; }



