#!/usr/bin/perl

use strict;
use FindBin '$Bin';
use lib "$Bin/../blib/lib","$Bin/../blib/arch";
use Bio::DB::Sam;

my $bamfile = shift or die "Usage $0 <bamfile>\n";

my $bam = Bio::DB::Sam->new(-bam=>$bamfile) 
    or die "Couldn't open $bamfile: $!";
$bam->coverage2BedGraph();

