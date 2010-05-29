file=Makefile
echo "Extracting $file"
rm -f $file

dbgobjdir=`echo "$objdir" | sed -e 's/\.[dbgprf]*$//'`
dbgobjdir=$dbgobjdir.dbg

dosversion=`echo $version|sed -e 's/\.//g'`
dosversion="_$dosversion$patch"
__status__=$status
case "$status" in
  development*) __status__='development';;
esac

cat > $file << EOT
# This file was created by Configure. All changes made will be lost
# next time Configure is run.
#
SHELL = $make_sh
VERS = pari-$version.$patch $__status__
SVN=svn+ssh://kb@pari.math.u-bordeaux.fr/home/svn/pari
TAG=release-$VersionMajor-$VersionMinor-$patch

dft target::
	@echo "Possible targets are:"
	@echo "	all             	 Compilation + Documentation"
	@echo "	gp			 Compilation"
	@echo "	bench                    Compilation + quick test"
	@echo "	dobench                	 Test only"
	@echo "	doc			 Documentation"
	@echo "	install			 Installation"
	@echo "	clean, cleantest         Clean up"
	@echo "For pari maintainers..."
	@echo "	dbg	                Compile gp binary suitable for debugging"
	@echo "	prf	                Compile gp binary suitable for profiling"
	@echo "	distrib,alpha,beta,release  Create a new distribution"
	@echo "	ctags			Generate VI/VIM tags file in ./src"
	@echo "	etags			Generate Emacs  tags file in ./src"
	@echo "	tune	                Generate tuning utility"
	@echo "Available extra tests: test-all $top_test_extra"

all::
	@-cd doc && \$(MAKE) doc
	@\$(MAKE) gp

gp bench test-kernel test-all $top_test_extra $top_dotest_extra install cleanobj cleantest install-bin install-doc install-data install-lib-sta install-bin-sta dobench tune::
	@dir=\`config/objdir\`; echo "Making \$@ in \$\$dir";\\
	 if test ! -d \$\$dir; then echo "Please run Configure first!"; exit 1; fi;\\
	cd \$\$dir && \$(MAKE) \$@

dbg gp.dbg::
	@dir=\`config/objdir\`.dbg; echo "Making gp in \$\$dir";\\
	 if test ! -d \$\$dir; then echo "Please run Configure -g first!"; exit 1; fi;\\
	cd \$\$dir && \$(MAKE) gp
	
prf gp.prf::
	@dir=\`config/objdir\`.prf; echo "Making gp in \$\$dir";\\
	 if test ! -d \$\$dir; then echo "Please run Configure -pg first!"; exit 1; fi;\\
	cd \$\$dir && \$(MAKE) gp
	
doc docps docpdf gpman cleandoc::
	cd doc && \$(MAKE) \$@

clean:: cleandoc cleanobj

distrib:
	$config_dir/settar \$(VERS) $__status__
tag:
	@echo svn copy \$(SVN)/branches/release-2-3-patches \$(SVN)/tags/\$(TAG)

alpha: tag
	$config_dir/settar \$(VERS) alpha
beta: tag
	$config_dir/settar \$(VERS) beta
release: tag
	$config_dir/settar \$(VERS) released

dosdistrib: gpman
	-zip -kr GPB${dosversion}.ZIP GP.EXE README README.DOS ../EMX/README.DOC ../RSX/README.TXT
	-zip -k GPD${dosversion}.ZIP doc/*.tex doc/gphelp doc/gp.man examples misc/gprc* misc/gpalias

etags:
	cd src; sh make_emacs_tags

ctags:
	src/make_vi_tags "$TOP/src"

EOT
