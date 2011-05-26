use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'GenomeDiff',
    license             => 'gpl',
    dist_author         => q{Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>},
    dist_version_from   => 'lib/GenomeDiff.pm',
    build_requires => {
        'Test::More' => 0,
    },
    add_to_cleanup      => [ 'GenomeDiff-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
