include ../../SETUPFILE

# Flags for the compiles
CFLAGS += -Wno-sign-compare $(HTSLIB_CFLAGS)
LIBS= $(HTSLIB_LIBS) -lhts -pthread -lz

all: qaCompute

qaCompute: qaCompute.o
	$(CXX) qaCompute.o -o qaCompute $(LIBS)

qaCompute.o: qaCompute.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf *o qaCompute

