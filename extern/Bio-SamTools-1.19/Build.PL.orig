#!/usr/bin/perl

use strict;
use Module::Build;

my $HeaderFile = "bam.h";
my $LibFile    = "libbam.a";
my $ReadLine;

my ($sam_include,$sam_lib) = find_sam(); # may exit with error here

my $build = Module::Build->new(
    module_name        => 'Bio::SamTools',
    dist_version_from  => 'lib/Bio/DB/Sam.pm',
    dist_author        => 'Lincoln Stein <lincoln.stein@gmail.com>',
    dist_abstract      => 'Perl interface to SamTools library for DNA sequencing',
    license            => 'perl',
    include_dirs       => [$sam_include],
    extra_linker_flags => ["-L$sam_lib",'-lbam','-lz'],

    extra_compiler_flags=>[

        # must match DFLAGS in Samtools Makefile
	'-D_IOLIB=2','-D_FILE_OFFSET_BITS=64', 

	# turn off warnings originating in Perl's Newx* calls
	'-Wformat=0'
    ],

    build_requires => {
      'ExtUtils::CBuilder' => 0,
    },
    requires     => {
	'perl'                => '5.008',
	'Bio::Root::Version'  => '1.005009001',
    },
#    create_makefile_pl => 'passthrough',
    );

$build->create_build_script;

exit 0;

sub find_sam {
    my ($sam_include,$sam_lib);

    if (my $samtools = $ENV{SAMTOOLS}) {
	$sam_include = $samtools
	    if -e "$samtools/$HeaderFile";
	$sam_lib     = $samtools
	    if -e "$samtools/$LibFile";
    }

    my @search_path = qw(/ /usr /usr/share /usr/local);

    unless ($sam_include) {
	for my $p (@search_path) {
	    $sam_include ||= "$p/include" if
		-e "$p/include/$HeaderFile";
	}
    }

    unless ($sam_lib) {
	for my $p (@search_path) {
	    $sam_lib ||= "$p/lib" if
		-e "$p/lib/$LibFile";
	}
    }


    unless ($sam_include && $sam_lib) {
	print STDOUT "This module requires samtools 0.1.4 or higher (samtools.sourceforge.net).\n";
	my $prompt = "Please enter the location of the bam.h and compiled libbam.a files: ";
	my $found;
	while (!$found) {
	    my $path = prompt($prompt);
	    print STDOUT "\n";
	    last unless $path;
	    $sam_include = $path
		if -e "$path/$HeaderFile";
	    $sam_lib = $path
		if -e "$path/$LibFile";
	    $found = $sam_include && $sam_lib;
	    unless ($found) {
		print STDOUT "That didn't seem to be right.\n";
		$prompt = "Try again, or hit <enter> to cancel: ";
	    }
	}
    }

    unless ($sam_include && $sam_lib) {
	die <<END;
Can\'t find $LibFile and/or $HeaderFile!

If you haven\'t done so already, please compile samtools 0.1.4 or
higher from http://samtools.sourceforge.net and set the SAMTOOLS
environment variable to point to a samtools distribution directory
containing the compiled $LibFile and $HeaderFile files.
END
}

    print STDOUT "Found $sam_include/$HeaderFile and $sam_lib/$LibFile.\n";
    return ($sam_include,$sam_lib);
}


sub prompt {
    my $msg = shift;

    unless (defined $ReadLine) {
	eval "require Term::ReadLine";
	$ReadLine   = Term::ReadLine->can('new') || 0;
	$ReadLine &&= Term::ReadLine->new(\*STDOUT);
	eval {readline::rl_set('TcshCompleteMode','On')};
    }

    
    unless ($ReadLine) {
	print STDOUT $msg;
	chomp (my $in = <>);
	return $in;
    }

    my $in = $ReadLine->readline($msg);
    chomp $in;
    $in=~ s/\s+//;
    $ReadLine->addhistory($in) if $in =~ /\S/;
    return $in;
}
