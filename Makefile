CXX ?= g++
CFLAGS = -std=c++11 -I /usr/local/include -I /usr/local/share/halide/tutorial/ -L /usr/local/lib/ -lHalide

#need to fix this, not working
#-g is for debugging, check whether it should be in both statements
#%.standalone1: %.o
#	$(CXX) $(CFLAGS) -DSTANDALONE $< -g -o $@
#
#%.o: %.cpp %.h
#	$(CXX) $(CFLAGS) -DSTANDALONE $< -g -c $@

#original before adding .h file
#careful using it with .h file, can't tell when .h has been changed
%.standalone: %.cpp
	$(CXX) $(CFLAGS) -DSTANDALONE $< -g -o $@


#need to fix this so it works with LSST stack
%: %.cpp
	$(CXX) $(CFLAGS) $< -g -o $@
