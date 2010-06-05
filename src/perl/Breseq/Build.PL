use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Breseq',
    license             => 'perl',
    dist_author         => q{Jeffrey E. Barrick <jeffrey@barricklab.org>},
    dist_version_from   => 'lib/Breseq.pm',
    build_requires => {
        'Test::More' => 0,
    },
    add_to_cleanup      => [ 'Breseq-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
