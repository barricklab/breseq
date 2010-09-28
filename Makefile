## Build and installation script

ROOTDIR=$(PWD)
BRESEQDIR=$(ROOTDIR)/src/perl/Breseq
BIOSAMTOOLS=$(ROOTDIR)/extern/Bio-SamTools-1.19
SAMTOOLSDIR=$(ROOTDIR)/extern/samtools-0.1.7a
BIOPERL=$(ROOTDIR)/extern/bioperl-live

STAGEDIR=$(ROOTDIR)/stage

LOCALUSERCONFIG=$(ROOTDIR)/user-config.jam
BJAMFLAGS=$(shell if [ -e $(LOCALUSERCONFIG) ]; then echo --user-config=$(LOCALUSERCONFIG); fi )
#echo "$(LOCALUSERCONFIG)"
#test -e $(LOCALUSERCONFIG) && BJAMFLAGS=--user-config=$(LOCALUSERCONFIG)


## Main
all : make

make :
	$(info $(BJAMFLAGS))

	bjam $(BJAMFLAGS)
	bjam $(BJAMFLAGS) install-libbam-for-perl

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


install :
	bjam $(BJAMFLAGS) install

	cd $(BRESEQDIR) ; \
	./Build install

	cd $(BIOPERL) ; \
	./Build install
	
	cd $(BIOSAMTOOLS) ; \
	./Build install
	
clean :
	bjam $(BJAMFLAGS) clean
	bjam $(BJAMFLAGS) clean install

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
	tests/test.sh test tests
	
clean-test:
	tests/test.sh clean tests
	
