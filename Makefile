#OPTIMIZATION=-ggdb -Og
OPTIMIZATION=-Ofast -funroll-loops

CPPFLAGS=$(OPTIMIZATION) -Wall -std=c++0x -I./includes/ -I/usr/include/hdf5/serial/ -fopenmp
LDFLAGS=-g -Wall -L/usr/lib/x86_64-linux-gnu/hdf5/serial
LDLIBS=-lm -lgsl -lgslcblas -lhdf5 -lpthread -fopenmp

SRCS=src/main.cpp src/abundances.cpp src/physics.cpp src/nuclear.cpp src/eos.cpp
SRCS_TESTS=src/tests.cpp src/abundances.cpp src/physics.cpp src/nuclear.cpp src/eos.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
OBJS_TESTS=$(subst .cpp,.o,$(SRCS_TESTS))
DEPS=$(wildcard includes/*.h)

all: $(OBJS) $(OBJS_TESTS) $(DEPS)
	$(CXX) $(LDFLAGS) -o run $(OBJS) $(LDLIBS) 
	$(CXX) $(LDFLAGS) -o tests $(OBJS_TESTS) $(LDLIBS) 

main.o: src/main.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/main.cpp

tests.o: src/tests.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/tests.cpp

eos.o: src/eos.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/eos.cpp

abundances.o: src/abundances.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/abundances.cpp

physics.o: src/physics.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/physics.cpp

nuclear.o: src/nuclear.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/nuclear.cpp

clean:
	rm -rf src/*.o
