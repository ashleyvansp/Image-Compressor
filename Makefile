EXTRA_CXXFLAGS=
CXXFLAGS=-O3 -Wall -std=c++17 $(EXTRA_CXXFLAGS)

all: uvg_compress uvg_decompress

uvg_compress: uvg_compress.o
	g++ -o $@ $^

uvg_decompress: uvg_decompress.o
	g++ -o $@ $^

clean:
	rm -f uvg_compress uvg_decompress *.o
