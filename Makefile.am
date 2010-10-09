ACLOCAL_AMFLAGS = -I aux_build/m4

SUBDIRS = extern/samtools-0.1.7a src/new-c

ROOTDIR=$(PWD)
BRESEQDIR=$(ROOTDIR)/src/perl/Breseq
BIOSAMTOOLS=$(ROOTDIR)/extern/Bio-SamTools-1.19
SAMTOOLSDIR=$(ROOTDIR)/extern/samtools-0.1.7a
BIOPERL=$(ROOTDIR)/extern/bioperl-live

STAGEDIR=$(ROOTDIR)/stage

## Main

all-local :

	## Breseq
	cd $(BRESEQDIR) ; \
	perl Build.PL --install_base=$(STAGEDIR) ; \
	./Build ; 
	
	## Bioperl
	cd $(BIOPERL) ; \
	perl Build.PL --accept --install_base=$(STAGEDIR) ; \
	./Build ; 
	
	## Bio::SamTools
	cd $(BIOSAMTOOLS) ; \
	export SAMTOOLS=$(ROOTDIR)/extern/samtools-0.1.7a ; \
	perl Build.PL --install_base=$(STAGEDIR) ; \
	./Build ; 


#utility to avoid everything else
breseq :

	cd $(BRESEQDIR) ; \
	perl Build.PL --install_base=$(STAGEDIR) ; \
	./Build

	cd $(BRESEQDIR) ; \
	./Build install

install-exec-local :

	cd $(BRESEQDIR) ; \
	./Build install

	cd $(BIOPERL) ; \
	./Build install
	
	cd $(BIOSAMTOOLS) ; \
	./Build install
	
clean-local :

	cd $(BRESEQDIR) ; \
	./Build clean

	cd $(BIOPERL) ; \
	./Build clean
	
	cd $(BIOSAMTOOLS) ; \
	./Build clean
	
	rm -rf $(STAGEDIR)
	
	tests/test.sh clean tests/	

## breseq only	
make-breseq:
	cd $(BRESEQDIR) ; \
	perl Build.PL --install_base=$(STAGEDIR) ; \
	./Build ; 

install-breseq:
	cd $(BRESEQDIR) ; \
	./Build install
	
clean-breseq :
	cd $(BRESEQDIR) ; \
	./Build clean
	
remake-breseq :: clean-breseq make-breseq install-breseq	

## tests

test:
	tests/lambda_mult_ref_read/testcmd.sh  test
	tests/lambda_mixed_pop/testcmd.sh  test
	tests/REL606_fragment_manipulated/testcmd.sh  test

test-all:
	tests/test.sh test tests

clean-tests:
	tests/test.sh clean tests
	
