CXX ?= g++
CFLAGS = -std=c++11 -I /usr/local/include -I /usr/local/share/halide/tutorial/ -L /usr/local/lib/ -lHalide

%.standalone: %.cpp
	$(CXX) $(CFLAGS) -DSTANDALONE $< -g -o $@

%: %.cpp
	$(CXX) $(CFLAGS) $< -g -o $@
