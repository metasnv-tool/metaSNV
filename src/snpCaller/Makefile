include ../../SETUPFILE

# -std=c++11 is required for older versions of G++
CFLAGS := $(CFLAGS) -std=c++0x -O2 -Wall -Wno-sign-compare

all: snpCall

snpCall: snpCall.o
	$(CXX) snpCall.o -o snpCall

snpCall.o: call_vC.cpp
	$(CXX) -c $(BOOST_CFLAGS) $(CFLAGS) $< -o $@

clean:
	rm -rf *o snpCall
