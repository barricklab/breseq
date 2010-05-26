all: samtools error_count

samtools:
	make -C samtools-0.1.7a

error_count: error_count/src/*.cpp
	g++ error_count/src/*.cpp \
		/opt/local/lib/libboost_program_options-mt.a \
		/opt/local/lib/libz.a \
		samtools-0.1.7a/libbam.a \
		-I error_count/include \
		-I samtools-0.1.7a \
		-I /opt/local/include \
		-o error_count/error_count++ \
