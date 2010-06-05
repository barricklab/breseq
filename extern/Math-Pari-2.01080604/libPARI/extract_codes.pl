#!/usr/bin/perl -wn
BEGIN {
  @ARGV = ("$ARGV[0]/src/language/init.c", "$ARGV[0]/src/gp/highlvl.c")
    if @ARGV == 1;
  @ARGV == 2 or die <<EOD;
Usage: $0  \$PARI_SRC_DIR
       $0  \$PARI_SRC_DIR/src/language/init.c \$PARI_SRC_DIR/src/gp/highlvl.c
EOD
  open EXP, 'expected_codes' or die "open expected_codes: $!";
  while (<EXP>) {
    die <<EOD unless ($ord, $descr) = /^(\d+)\t(.*)/;
cannot parse line of expected_codes: $_
EOD
    $expected{$ord} = $descr;
  }
  close EXP or die "close expected_codes: $!";
}
next unless /^\s*{\s*\"/;
chomp;
warn("Unrecognized line: `$_'\n"), next
  unless /^\s*\{\s*"(\w+)"\s*,\s*(\d+)\s*,[^,]*,\s*\d+\s*,(?:\s*\d+\s*,)?\s*("((?:\\.|[^\s"])*)"|NULL)\s*(,|\})/;
next unless defined $4;
#print;
(undef, $code, $descr) = ($1, $2, $4);

$descr{$code} = [] unless exists $descr{$code};
push @{$descr{$code}}, $descr unless grep $descr eq $_, @{$descr{$code}};

END {
  for $k (sort {$a <=> $b} keys %descr) {
    next if $k == 99;
    if (@{$descr{$k}} > 1) {
      warn "Multiple descriptors for code $k: @{$descr{$k}}\n"
	if $k;
    } elsif (@{$descr{$k}} == 0){
      warn "Empty descriptors for code $k."
    } elsif (exists $expected{$k} and $expected{$k} ne $descr{$k}[0]){
      warn <<EOW
Unexpected descriptor for code $k: `$descr{$k}[0]', expecting `$expected{$k}'
EOW
    } else {
      print "$k\t$descr{$k}[0]\n";
    }
  }
}
