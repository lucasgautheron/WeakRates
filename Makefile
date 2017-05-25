#OPTIMIZATION=-ggdb -Og
OPTIMIZATION=-Ofast
FC=gfortran

CPPFLAGS=$(OPTIMIZATION) `root-config --cflags` -Wall -std=c++0x -I./includes/ -I/usr/include/hdf5/serial/ -fopenmp -D DEBUG
LDFLAGS=-g -Wall -Wextra -L/usr/lib/x86_64-linux-gnu/hdf5/serial
LDLIBS=-lm -lgsl -lgslcblas -lhdf5 -lpthread -fopenmp -lgfortran

CPPFLAGS_TESTS=`root-config --cflags`
LDLIBS_TESTS=`root-config --glibs` -lgslcblas

SRCS=src/main.cpp src/abundances.cpp src/physics.cpp src/nuclear.cpp src/eos.cpp
SRCS_TESTS=src/tests.cpp src/abundances.cpp src/physics.cpp src/nuclear.cpp src/eos.cpp
SRCS_F=src/dz.for

OBJS=$(subst .cpp,.o,$(SRCS))
OBJS_TESTS=$(subst .cpp,.o,$(SRCS_TESTS))
OBJS_F=$(subst .for,.o,$(SRCS_F))
DEPS=$(wildcard includes/*.h)

all: $(OBJS) $(OBJS_F) $(DEPS)
	$(CXX) $(LDFLAGS) -o run $(OBJS_F) $(OBJS) $(LDLIBS)

tests: $(OBJS_F) $(OBJS_TESTS) $(DEPS)
	$(CXX) $(LDFLAGS) -o tests $(OBJS_F) $(OBJS_TESTS) $(LDLIBS) $(LDLIBS_TESTS) 

src/dz.o: src/dz.for
	$(FC) -c src/dz.for -o src/dz.o

main.o: src/main.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/main.cpp

tests.o: src/tests.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) $(CPPFLAGS_TESTS) -c src/tests.cpp

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
	rm run
	rm tests
