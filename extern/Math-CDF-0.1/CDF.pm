package Math::CDF;

use strict;
use Carp;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $AUTOLOAD);

require Exporter;
require DynaLoader;
require AutoLoader;

@ISA = qw(Exporter DynaLoader);

@EXPORT_OK = qw(pnorm qnorm pt qt pbeta qbeta pchisq qchisq
		pf qf pgamma qgamma ppois qpois
		pbinom pnbinom
);

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

$VERSION = '0.1';

sub AUTOLOAD {
    # This AUTOLOAD is used to 'autoload' constants from the constant()
    # XS function.  If a constant is not found then control is passed
    # to the AUTOLOAD in AutoLoader.

    my $constname;
    ($constname = $AUTOLOAD) =~ s/.*:://;
    croak "& not defined" if $constname eq 'constant';
    my $val = constant($constname, @_ ? $_[0] : 0);
    if ($! != 0) {
	if ($! =~ /Invalid/) {
	    $AutoLoader::AUTOLOAD = $AUTOLOAD;
	    goto &AutoLoader::AUTOLOAD;
	}
	else {
		croak "Your vendor has not defined CDF macro $constname";
	}
    }
    no strict 'refs';
    *$AUTOLOAD = sub () { $val };
    goto &$AUTOLOAD;
}

bootstrap Math::CDF $VERSION;

# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is the stub of documentation for your module. You better edit it!

=head1 NAME

Math::CDF - Generate probabilities and quantiles from several 
statistical probability functions

=head1 SYNOPSIS

use Math::CDF;

$prob = C<&Math::CDF::pnorm(1.96)>;

C<if( not defined($z = &Math::CDF::qnorm(0.975)) )> {
    die "C<qnorm()> failed"; }

or

use Math::CDF C<qw(:all)>;

$prob = C<pnorm(1.96)>;

=head1 DESCRIPTION

This module provides a perl interface to the DCDFLIB. See the section on
L<DCDFLIB> for more information.

Functions are available for 7 continuous distributions (Beta,
Chi-square, F, Gamma, Normal, Poisson and T-distribution) and for 
two discrete distributions (Binomial and Negative Binomial). Optional
non-centrality parameters are available for the Chi-square, F and 
T-distributions. Cumulative probabilities are available for all 9 
distributions and quantile functions are available for the 7 continuous 
distributions.

All cumulative probability function names begin with the character "p".
They give the probability of being less than or equal to the given value
S<[ C<P(X <= x)> ]>

All quantile function names begin with the character q. They give a value
of x such that S<C<P(X <= x)> = p> where the value of p is provided to the function.

Non-centrality parameters are always the last function argument when
available. You do not need to supply the non-centrality parameter in which case
it will be assumed to be 0.

All functions will return an undefined value if the function fails 
(probably due to parameters being out of allowed range) but will not otherwise
generate an error message. The user should check for valid output from the
Math::CDF functions with the C<defined()> function as demonstrated in the 
L<SYNOPSIS> section.

=head1 FUNCTION DESCRIPTIONS

In all, 16 functions are available via Math::CDF:

    pbeta(), qbeta()          [Beta Distribution]
    pchisq(), qchisq()        [Chi-square Distribution]
    pf(), qf()                [F Distribution]
    pgamma(), qgamma()        [Gamma Distribution]
    pnorm(), qnorm()          [Standard Normal Dist]
    ppois(), qpois()          [Poisson Distribution]
    pt(), qt()                [T-distribution]
    pbinom()                  [Binomial Distribution]
    pnbinom()                 [Negative Binomial Distribution]

=over 4

=item C<pbeta($x, $a, $b), qbeta($p, $a, $b)>

Generates cumulative probabilities and quantiles from the 
Beta distribution.
$x should be in the range of [0,1] and $p should be in [0,1].
$a and $b are parameters of the beta distribution. 
Both $a and $b must be in the range [0,Inf).

=item C<pchisq($x, $df, $ncp), qbeta($p, $df, $ncp)>

Generates cumulative probabilities and quantiles from the 
Chi-square distribution.
$x should be in the range of [0,Inf) and $p should be in [0,1].
$df is the degrees of freedom of the distribution and must be in the 
range (0,Inf).
$ncp is the optional non-centrality parameter and must be in the range [0,Inf).
The non-centrality parameter need not be specified and the calls 
C<pchisq(5, 5, 0.0)> and C<pchisq(5, 5)> will return identical values.

=item C<pf($x, $dfn, $dfd, $ncp), qf($p, $dfn, $dfd, $ncp)>

Generates cumulative probabilities and quantiles from the 
F distribution.
$x should be in the range of [0,Inf) and $p should be in [0,1].
$dfn and $dfd are the numerator and denominator degrees of freedom, 
respectively. Both must be in the range (0,Inf).
$ncp is the optional non-centrality parameter and must be in the range [0,Inf).
The non-centrality parameter need not be specified and the calls 
[CS]<pf(5, 2, 3, 0.0)> and [CS]<pf(5, 2, 3)> will return identical values.

=item C<pgamma($x, $shape, $scale), qgamma($p, $shape, $scale)>

Generates cumulative probabilities and quantiles from the 
Gamma distribution. The gamma density is proportional to 
[CS]<$x**($shape - 1) * EXP(- $scale * $x)>
$x should be in the range of [0,Inf) and $p should be in [0,1].
$shape and $scale are parameters of the Gamma distribution and
both must be in the range (0,Inf).

=item C<pnorm($x), qnorm($p)>

Generates cumulative probabilities and quantiles from the 
standard Normal distribution.
$x should be in the range of (-Inf,Inf) and $p should be in [0,1].

=item C<ppois($x, $lambda), qpois($p, $lambda)>

Generates cumulative probabilities and quantiles from the 
Poisson distribution.
$x should be in the range of [0,Inf) and $p should be in [0,1].
$lambda is the parameter of the Poisson distribution and
must be in the range [0,Inf).

=item C<pt($x, $df, $ncp), qt($p, $df, $ncp)>

Generates cumulative probabilities and quantiles from the 
T distribution.
$x should be in the range of (-Inf,Inf) and $p should be in [0,1].
$df is the degree of freedom of the distribution and 
must be in the range (0,Inf).
$ncp is the optional non-centrality parameter and must be in the 
range (-Inf,Inf).
The non-centrality parameter need not be specified and the calls 
[CS]<pt(0, 3, 0.0)> and [CS]<pt(0, 3)> will return identical values.

=item C<pbinom($x, $n, $p)>

Generates cumulative probabilities from the 
Binomial distribution.
This is the probability of having $x or fewer successes in $n trials
when each trial has a $p probability of success.
$x should be in the range of [0,$n], $n should be in the range (0,Inf)
and $p should be in [0,1].

=item C<pnbinom($x, $n, $p)>

Generates cumulative probabilities from the 
Negative Binomial distribution.
The is the probability of having $x or fewer failures before the 
$n'th success when each trial has a $p probability of success.
$x should be in the range of [0,Inf), $n should be in the range (0,Inf)
and $p should be in [0,1].

=back

=head1 DCDFLIB

DCDFLIB is a library of C routines for cumulative distribution functions, 
inverses, and other parameters written by Barry W. Brown, James Lovato and 
Kathy Russell of the Department of Biomathematics at The University of 
Texas M.D. Anderson Cancer Center in Houston Texas. 

Version 1.1 of DCDFLIB is included with this distribution and can be 
downloaded via ftp from odin.mda.uth.tmc.edu as /pub/src/dcdflib.c-1.1.tar.gz.
The library is also available in Fortran source.

Documentation for DCDFLIB is found in the cdflib/doc/ directory of the
source distribution of Math::CDF

DCDFLIB has been put in the public domain. However, some of the Algorithms
are from the ACM publication and is under their copyright. Generally this
means that those algorithms can be distributed and used freely for 
non-commercial purposes. See the file cdflib/doc/README in the source
distribution of Math::CDF for more information.

=head1 AUTHOR

=over 4

=item *

Math::CDF was put together by Ed Callahan, callahan@envstat.com

Environmental Statistics
PO Box 563
Fountain City, WI  54629

=item *

DCDFLIB was written by Barry W. Brown (bwb@odin.mdacc.tmc.edu), 
James Lovato and Kathy Russell

Department of Biomathematics, Box 237
The University of Texas, M.D. Anderson Cancer Center
1515 Holcombe Boulevard
Houston, TX  77030

=back

=head1 SEE ALSO

L<Math::Random>, L<Statistics::ChiSquare>, L<Statistics::OLS>

=cut
