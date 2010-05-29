#!/usr/bin/perl -w
use strict;

$| = 1;
print "1..1\nok 1\n";	# report success

# `centerlift' is usually the first failing symbol detected...
exit if eval 'use Math::Pari; use Math::Pari "centerlift"; 1';

# If failed, report build parameters

sub report_build_parameters ($) {
  my ($makefile, $in) = (shift, '');
  warn "# reporting $makefile header:\n# ==========================\n";
  open M, "< $makefile" or die "Can't open $makefile";
  $in = <M> while defined $in and $in !~ /MakeMaker \s+ Parameters/xi;
  $in = <M>;		# Reload till first non-empty line after "Param" one:
  $in = <M> while defined $in and $in !~ /\S/;
  warn $in and $in = <M> while defined $in and $in =~ /^#/;
  close M;
  warn "# ==========================\n";
}

my ($base_d, $in) = (-f "t/000_load-problem.t" ? '.' : '..', '');
report_build_parameters("$base_d/Makefile");
report_build_parameters("$base_d/libPARI/Makefile");
