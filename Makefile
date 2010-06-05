## Build and installation script

## Handle boost builds

ROOTDIR=$(PWD)

BRESEQDIR=$(ROOTDIR)/src/perl/Breseq
BIOSAMTOOLS=$(ROOTDIR)/extern/Bio-SamTools-1.19
SAMTOOLSDIR=$(ROOTDIR)/extern/samtools-0.1.7a
STATSDISTS=$(ROOTDIR)/extern/Statistics-Distributions-1.02
PARIDIR=$(ROOTDIR)/extern/Math-Pari-2.01080604
STAGEDIR=$(ROOTDIR)/stage
	
all :: make

make :

	bjam
	bjam install-libbam-for-perl

	## Breseq
	cd $(BRESEQDIR) ; \
	perl Build.PL --install_base=$(STAGEDIR) ; \
	./Build ; 

	## Bio::SamTools
	cd $(BIOSAMTOOLS) ; \
	export SAMTOOLS=$(ROOTDIR)/extern/samtools-0.1.7a ; \
	perl Build.PL --install_base=$(STAGEDIR) ; \
	./Build ; 
	
	## Statistics::Distributions
	cd $(STATSDISTS) ; \
	perl Makefile.PL INSTALL_BASE=$(STAGEDIR) ; \
	make ; 
	
	## Math::PARI
	cd $(PARIDIR) ; \
	perl Makefile.PL INSTALL_BASE=$(STAGEDIR) ; \
	make ; 
	
clean-breseq :
	cd $(BRESEQDIR) ; \
	./Build clean
	
make-breseq:
	cd $(BRESEQDIR) ; \
	perl Build.PL --install_base=$(STAGEDIR) ; \
	./Build ; 

install-breseq:
	cd $(BRESEQDIR) ; \
	./Build install
	
all-breseq :: clean-breseq make-breseq install-breseq
	
clean :
	bjam clean
	bjam clean install

	cd $(BRESEQDIR) ; \
	./Build clean
	
	cd $(BIOSAMTOOLS) ; \
	./Build clean
	
	cd $(STATSDISTS) ; \
	make clean
	
	cd $(PARIDIR) ; \
	make clean
	
	rm -rf $(STAGEDIR)
	
	tests/test.sh clean tests/	
	
install :
	bjam install

	cd $(BRESEQDIR) ; \
	./Build install

	cd $(BIOSAMTOOLS) ; \
	./Build install

	cd $(STATSDISTS) ; \
	make install	
	
	cd $(PARIDIR) ; \
	make install
	
test:
	tests/test.sh test tests
	
clean-test:
	tests/test.sh clean tests