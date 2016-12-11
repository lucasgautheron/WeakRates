OPTIMIZATION=-ggdb -Og

CPPFLAGS=$(OPTIMIZATION) -Wall -std=c++0x -I./includes/ -I/usr/include/hdf5/serial/ -fopenmp
LDFLAGS=-g -Wall -L/usr/lib/x86_64-linux-gnu/hdf5/serial
LDLIBS=-lm -lgsl -lgslcblas -lhdf5 -lpthread -fopenmp

SRCS=src/main.cpp src/abundances.cpp src/physics.cpp src/nuclear.cpp
OBJS=$(subst .cpp,.o,$(SRCS))
DEPS=$(wildcard includes/*.h)

all: $(OBJS) $(DEPS)
	$(CXX) $(LDFLAGS) -o run $(OBJS) $(LDLIBS) 

main.o: src/main.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/main.cpp

abundances.o: src/abundances.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/abundances.cpp

physics.o: src/physics.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/physics.cpp

nuclear.o: src/nuclear.cpp $(DEPS)
	$(CXX) $(CPPFLAGS) -c src/nuclear.cpp

clean:
	rm -rf src/*.o
