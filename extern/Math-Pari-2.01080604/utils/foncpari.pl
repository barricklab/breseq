unshift @ARGV, 'Pari.xs' if @ARGV < 2 and -r 'Pari.xs';
unshift @ARGV, 'libPARI/anal.c' if @ARGV < 2 and -r 'libPARI/anal.c';
(@ARGV == 2) || &usage;

@known = split /,\s*/, 'label, while, goto, until, read, pprint, print, texprint, pprint1, print1, O, if, o';
@known{@known} = (1) x @known;

open(XSUB,$ARGV[1]) || die "Cannot open $ARGV[1]: $!";
while (<XSUB>) {
  $supported{$1}++ if /^\s*CASE_INTERFACE\s*\(\s*(\d+)\s*\)/;
}
close(XSUB) || die "Cannot close $ARGV[1]: $!";

open(ANAL,$ARGV[0]) || die "Cannot open $ARGV[0]: $!";
while (<ANAL>) {
  if (/^entree\s+fonctions\[/ ... /^\s*\}\s*;\s*$/) {
    next unless $i++;		# Skip first line
    last if /^\s*\}\s*;\s*$/;
    &warnl() unless /
		      ^ \s* \{ \s* " 
		      (
			[^""]+	# 1 Name
		      )
		      " \s* , \s* 
		      (
			\d+	# 2 Interface
		      )
		      \s* , \s* 
		      (
			[^,]+	# 3 C function pointer
		      )
		      \s* , \s* 
		      (
			\d+	# 4 Group
		      )
		      \s* , \s* 
		      (
			\d+	# 5 
		      )
		      ( \s* , \s* ((" [^"]* ") | NULL) \s* , \s* NULL )? # New fields
		      \s* \} \s* ,? \s* $
		    /x; # ";
    ($pari, $interface, $gp, $group, $code) = ($1, $2, $3, $4, $8);
    if ($gp eq "0") {
      if ($known{$pari}) {
	$builtin_known{$pari}++;
      } else {
	$builtin{$pari}++;
      }
    } else {
      $interface{$pari} = $interface;
      $code{$interface} ||= ($code || '');
      push @{$group{$group}}, $pari unless exists $supported{$interface};
      $interfaces{$interface}++;
    }
    # print "'$pari' <= '$gp' via $interface\n";
    # &warnl() unless $gp =~ /\b$pari\b/;
  }
}
close(ANAL) || die "Cannot close $ARGV[0]: $!";
print "Builtins, unsupported as functions (but available in Perl):\n\t", join(", ", keys %builtin_known), "\n\n"
  if %builtin_known;

print "Builtins, completely unsupported:\n\t", join(", ", keys %builtin), "\n\n"
  if %builtin;

for (keys %interfaces) {
  $unsupported{$_}++ unless $supported{$_};
}

@unsupported = sort {$interfaces{$a} <=> $interfaces{$b}} keys %unsupported;

print "\tTotal number of unsupported interfaces: ",scalar @unsupported,":\n";
for $i (sort {$a <=> $b} @unsupported) {
  print "Interface $i=$code{$i} used in $interfaces{$i} function(s): ",
     join(", ", @f=grep($interface{$_}==$i, keys %interface)), ".\n";
  if ($code{$i}) {
    $write = write_interface($i,$code{$i});
    $suggest{$i} = $write if defined $write;
  }
  $total += $interfaces{$i};
  push(@ff,@f);
}

print "\n\tTotal number of unsupported functions: $total:\n";
  #join(", ", sort @ff), "\n";

for $g (sort {$a <=> $b} keys %group) {
  print "group $g:\t", join(', ', sort @{$group{$g}}), "\n";
}

if (%suggest) {
  print "Suggested code for interfaces:\n\n";
  
  for $i (sort keys %suggest) {
    print $suggest{$i};
  }
  for $i (sort keys %suggest) {
    print <<EOI;
	   CASE_INTERFACE($i);
EOI
  }
  print "\n";
  for $i (sort keys %suggest) {
    print <<EOI;
		   case $i:
EOI
  }
}

sub usage {die "Usage: $0 [path/to/anal.c] [path/to/Pari.xs]\n"}

sub warnl {warn "Unrecognized line:\n$_"}

sub write_interface {
  my ($num, $interface) = @_;
  my ($int) = $interface =~ /^"(.*)"$/ or return;
  my @types;
  my @c_arg_names;
  my $ret_type = 'GEN';
  
  while (length $int) {
    if ($int =~ s/^s//) {
      push @types, 'char *';
      push @c_arg_names, 'arg' . scalar @types;
    } elsif ($int =~ s/^l//) {
      $ret_type = 'long';
    } elsif ($int =~ s/^L//) {
      push @types, 'long';
      push @c_arg_names, 'arg' . scalar @types;
    } elsif ($int =~ s/^V=//) {
      push @types, 'PariVar';
      push @c_arg_names, 'arg' . scalar @types;
    } elsif ($int =~ s/^I//) {
      push @types, 'PariExpr';
      push @c_arg_names, 'arg' . scalar @types;
    } elsif ($int =~ s/^G//) {
      push @types, 'GEN';
      push @c_arg_names, 'arg' . scalar @types;
    } elsif ($int =~ s/^p//) {
      push @c_arg_names, 'prec';
    } else {
      print "tail `$int' of interface$num unsupported\n";
      return;
    }
  }
  my @args = map {"arg$_"} 1 .. @types;
  my $args = join ', ', @args;
  my $c_args = join ', ', @c_arg_names;
  my $i = 0;
  $argdecl = join '', map {$i++; "    $_ arg$i\n"} @types;
  
  my $out = <<EOA;

$ret_type
interface$num($args)
long	oldavma=avma;
EOA
  $out .= $argdecl;
  $out .= <<EOA;
 CODE:
  {
    dFUNCTION($ret_type);

    if (!FUNCTION) {
      croak("XSUB call through interface did not provide *function");
    }

    RETVAL=FUNCTION($c_args);
  }
 OUTPUT:
   RETVAL
EOA
  if ($ret_type ne 'GEN') {
    $out .= <<EOA;
 CLEANUP:
   avma=oldavma;
EOA
  }
  $out .= "\n";
  $out;
}
