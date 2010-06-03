all: error_count test

error_count:
	make -C src/c	
	cp src/c/error_count/error_count++ bin/

test:
	./bin/test.sh test tests/
