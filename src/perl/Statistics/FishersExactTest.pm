=head1 NAME

Statistics::FET - Fisher's Exact Test statistic (2x2 case)

=head1 SYNOPSIS

  use Statistics::FET 'fishers_exact';

=head1 DESCRIPTION

This module exports only one function, C<fishers_exact>, which
computes the one- or two-sided Fisher's Exact Test statistic for the 2
x 2 case.  In the following documentation I will be referring to the
following family of 2 x 2 contingency tables


                        *   |   *   | r1 = a+b
                    --------+-------+----------
                        *   |   *   | r2 = c+d
                    --------+-------+----------
                       c1   |   c2  |    N
                     = a+c  | = b+d | = a+b+c+d

The *'s mark the cells, N is the total number of points represented by
the table, and r1, r2, c1, c2 are the marginals.  As suggested by the
equalities, the letters a, b, c, d refer to the various cells (reading
them left-to-right, top-to-bottom).  Depending on context, the letter
c (for example) will refer *either* to the lower left cell of the
table *or* to a specific value in that cell.  Same for a, b, and d.

In what follows, H(x) (or more precisely, H(x; r1, r2, c1)) refers
to the hypergeometric expression

                           r1!*r2!*c1!*c2!
                -------------------------------------
                (r1+r2)!*x!*(r1-x)!*(c1-x)!*(r2-c1+x)!

(I omit c2 from the parametrization of H because c2 = r1 + r2 - c1.)

=head1 FUNCTION

=over 4

=item fishers_exact( $a, $b, $c, $d, $two_sided )

The paramater C<$two_sided> is optional.  If missing or false
C<fishers_exact> computes the one-sided FET p-value.  More
specifically, it computes the sum of H(x; a+b, c+d, a+c) for x = a to
x = min(a+b, a+c) - "the right side".  (If you want "the left side", i.e. the sum of H(x; a+b, c+d, a+c) for x
= max(0, a-d) to x = a, then compute C<fishers_exact( $b, $a, $d, $c )
+> or
C<fishers_exact( $c, $d, $a, $b )> (these two are equivalent).)

If C<$two_sided> is true, the returned p-value will be for the two-sid
+ed
FET.

=back

=cut

## let the code begin...

use strict;
use warnings;

package Statistics::FishersExactTest;
use Exporter 'import';

use Math::Pari ();

my $Tolerance = 1;
$Tolerance /= 2 while 1 + $Tolerance/2 > 1;

@Statistics::FET::EXPORT_OK = 'fishers_exact';
sub fishers_exact {
  my ( $a, $b, $c, $d, $ts ) = @_;
 TEST:
  ## simplified test equation
  my $test = $a*$d - $b*$c;
  ## introduced switching around of input pairs
  if ($test < 0 && $ts){
    ($a, $b, $c, $d, ) = ($b, $a, $d, $c,);
    goto TEST;
  }

  if ($test < 0 and $ts){
    die "Cannot compute two tailed FET for input values given \( ".(join ", ", @_)."\)\n";
  }
  # below here, $test < 0 implies !$ts;

  my $p_val;
  if ( $test < 0 ) {
    if ( $d < $a ) {
      $p_val = _fishers_exact( $d, $c, $b, $a, $ts, 1 );
    }
    else {
      $p_val = _fishers_exact( $a, $b, $c, $d, $ts, 1 );
    }
  }
  else {
    if ( $b < $c ) {
      $p_val = _fishers_exact( $b, $a, $d, $c, $ts, 0 );
    }
    else {
      $p_val = _fishers_exact( $c, $d, $a, $b, $ts, 0 );
    }
  }
  return $p_val;
}

sub _fishers_exact {
  my ( $a, $b, $c, $d, $ts, $complement ) = @_;
  die "Bad args\n" if $ts && $complement;

  my ( $aa, $bb, $cc, $dd ) = ( $a, $b, $c, $d );
  my $first = my $delta = _single_term( $aa, $bb, $cc, $dd );
  my $p_val = 0;

  {
    $p_val += $delta;
    last if $aa < 1;
    $delta *= ( ( $aa-- * $dd-- )/( ++$bb * ++$cc ) );
    redo;
  }

  if ( $ts ) {
    my $m = $b < $c ? $b : $c;
    ($aa, $bb, $cc, $dd ) = ( $a + $m, $b - $m, $c - $m, $d + $m );
    $delta = _single_term( $aa, $bb, $cc, $dd );
    my $bound = -$Tolerance;
    while ( $bound <= ( $first - $delta )/$first && $aa > $a ) {

      $p_val += $delta;
      $delta *= ( ( $aa-- * $dd-- )/( ++$bb * ++$cc ) );
    }
  }
  elsif ( $complement ) {
    $p_val = 1 - $p_val + $first;
  }
  return $p_val;
}

sub _single_term {
  my ( $a, $b, $c, $d ) = @_;
  my ( $r1, $r2 ) = ($a + $b, $c + $d);
  my ( $c1, $c2 ) = ($a + $c, $b + $d);
  my $N = $r1 + $r2;

  return  exp( _ln_fact( $r1 ) + _ln_fact( $r2 ) +
               _ln_fact( $c1 ) + _ln_fact( $c2 ) -
               _ln_fact( $N ) -
               ( _ln_fact( $a ) + _ln_fact( $b ) +
                 _ln_fact( $c ) + _ln_fact( $d ) ) );
}

{
  my $two_pi;
  my $pi_over_3;
  my $half;

  BEGIN {
    $two_pi    = Math::Pari::PARI( 2 * atan2 0, -1 );
    $pi_over_3 = Math::Pari::PARI( atan2( 0, -1 )/3 );
    $half      = Math::Pari::PARI( 0.5 );
  }

  sub _ln_fact {
    my $n = Math::Pari::PARI( shift() );
    die "Bad args to _ln_fact: $n" if $n < 0;
    my $ln_fact;
    eval {
      $ln_fact = log Math::Pari::factorial( $n );
    };
    if ( $@ ) {
      die $@ unless $@ =~ /\QPARI:   ***   exponent overflow/;
      # Gosper's approximation; from
      # http://mathworld.wolfram.com/StirlingsApproximation.html
      $ln_fact = $half * log( $two_pi*$n + $pi_over_3 )
                 +  $n * log( $n )
                 -  $n;
    }
    return $ln_fact;
  }
}
1;

__END__
