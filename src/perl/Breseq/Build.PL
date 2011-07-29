=head1 AUTHOR

Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>

=head1 LICENSE

Copyright (C) 2010 Michigan State University.

This is free software; you can redistribute it and/or modify it under the terms the 
GNU General Public License as published by the Free Software Foundation; either 
version 1, or (at your option) any later version.

=cut

use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Breseq',
    license             => 'gpl',
    dist_author         => q{Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>},
    dist_abstract		=> 'breseq program for genome re-sequencing data analysis',
    dist_version_from   => 'lib/Breseq.pm',
    build_requires => {
        'Test::More' => 0,
    },
    add_to_cleanup      => [ 'Breseq-*', '_build' ],
    create_makefile_pl => 'traditional',
);

## additional file types to copy within lib directory
$builder->add_build_element('png');
$builder->add_build_element('r');

$builder->create_build_script();
