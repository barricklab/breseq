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

## additional file types to copy within lib directory
$builder->add_build_element('png');
$builder->add_build_element('r');

$builder->create_build_script();
