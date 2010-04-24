all: test

test:
	./bin/test.sh test tests

clean : 
	rm -r ${TEST_OUTPUT}