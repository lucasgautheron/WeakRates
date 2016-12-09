CPPFLAGS=-ggdb -Og -std=c++0x -Wall -I./includes/ -I/usr/include/hdf5/serial/
LDFLAGS=-g -Wall -L/usr/lib/x86_64-linux-gnu/hdf5/serial
LDLIBS=-lm -lgsl -lgslcblas -lhdf5

SRCS=src/main.cpp src/abundances.cpp src/physics.cpp src/nuclear.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: $(OBJS)
	g++ $(LDFLAGS) -o run $(OBJS) $(LDLIBS) 

main.o: src/main.cpp includes/common.h
	g++ $(CPPFLAGS) -c src/main.cpp

abundances.o: src/abundances.cpp includes/common.h
	g++ $(CPPFLAGS) -c src/abundances.cpp

physics.o: src/physics.cpp includes/common.h
	g++ $(CPPFLAGS) -c src/physics.cpp

nuclear.o: src/nuclear.cpp includes/common.h
	g++ $(CPPFLAGS) -c src/nuclear.cpp

clean:
	rm -rf src/*.o
