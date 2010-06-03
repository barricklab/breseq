use ExtUtils::MakeMaker;
WriteMakefile(
	      NAME => 'Math::CDF::cdflib',
	      SKIP => [qw(all static static_lib dynamic 
                          dynamic_lib test_dynamic
			  test)],
	      clean => {'FILES' => 'libcdflib$(LIB_EXT)'},
	      );

sub MY::top_targets {
'
test ::

all :: static

static :: libcdflib$(LIB_EXT)

libcdflib$(LIB_EXT): $(O_FILES)
	$(AR) cr libcdflib$(LIB_EXT) $(O_FILES)
	$(RANLIB) libcdflib$(LIB_EXT)
';
}
