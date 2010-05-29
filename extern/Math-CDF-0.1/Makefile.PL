use ExtUtils::MakeMaker;
WriteMakefile(
    'NAME'	   => 'Math::CDF',
    'VERSION_FROM' => 'CDF.pm',
    'LIBS'         => ['-lm'],
    'MYEXTLIB'     => 'cdflib/libcdflib$(LIB_EXT)',
    'dist'         => { COMPRESS=>"gzip -9", SUFFIX=>"gz" },
);

sub MY::postamble {
'
$(MYEXTLIB): cdflib/Makefile
	cd cdflib && $(MAKE) $(PASTHRU)
';
}
