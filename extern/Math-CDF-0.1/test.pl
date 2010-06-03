# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..22\n"; }
END {print "not ok 1\n" unless $loaded;}
use Math::CDF;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

print abs(&Math::CDF::pnorm(0) - 0.5) <= 0.000001 ? "ok 2\n" : "not ok 2\n";
print abs(&Math::CDF::qnorm(0.5)) <= 0.000001 ? "ok 3\n" : "not ok 3\n";

print abs(&Math::CDF::pt(0, 5) - 0.5) <= 0.000001 ? "ok 4\n" : "not ok 4\n";
print !&Math::CDF::pt(0, -1) ? "ok 5\n" : "not ok 5\n";
print abs(&Math::CDF::qt(0.5, 5)) <= 0.000001 ? "ok 6\n" : "not ok 6\n";
print !&Math::CDF::qt(0.5, -1) ? "ok 7\n" : "not ok 7\n";
print !&Math::CDF::qt(2, 5) ? "ok 8\n" : "not ok 8\n";

print abs(&Math::CDF::pbinom(5, 11, 0.5) - 0.5) <= 0.000001 ? "ok 9\n" : "not ok 9\n";
print !&Math::CDF::pbinom(12, 11, 0.5) ? "ok 10\n" : "not ok 10\n";
print !&Math::CDF::pbinom(5, 11, -0.5) ? "ok 11\n" : "not ok 11\n";

print abs(&Math::CDF::pnbinom(5, 10, 0.5) - 0.150878) <= 0.000001 ? "ok 12\n" : "not ok 12\n";

$p = &Math::CDF::pbeta(0.25, 2, 3);
print abs($p - 0.26171875) <= 0.000001 ? "ok 13\n" : "not ok 13\n";
print abs(0.25 - &Math::CDF::qbeta($p, 2, 3)) <= 0.000001 ? "ok 14\n" : "not ok 14\n";

$p = &Math::CDF::pchisq(5, 5);
print abs($p - 0.5841198) <= 0.000001 ? "ok 15\n" : "not ok 15\n";
print abs(5 - &Math::CDF::qchisq($p,5)) <= 0.000001 ? "ok 16\n" : "not ok 16\n";

$p = &Math::CDF::pf(5, 2, 3);
print abs($p - 0.889142) <= 0.000001 ? "ok 17\n" : "not ok 17\n";
print abs(5 - &Math::CDF::qf($p, 2, 3)) <= 0.000001 ? "ok 18\n" : "not ok 18\n";

$p = &Math::CDF::pgamma(1, 2, 3);
print abs($p - 0.8008517) <= 0.000001 ? "ok 19\n" : "not ok 19\n";
print abs(1 - &Math::CDF::qgamma($p, 2, 3)) <= 0.000001 ? "ok 20\n" : "not ok 20\n";

$p = &Math::CDF::ppois(5, 5);
print abs($p - 0.6159606) <= 0.000001 ? "ok 21\n" : "not ok 21\n";
print abs(5 - &Math::CDF::qpois($p, 5)) <= 0.000001 ? "ok 22\n" : "not ok 22\n";
