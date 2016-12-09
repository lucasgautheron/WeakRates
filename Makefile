CPPFLAGS=-g -std=c++0x -Wall -I./includes/ -I/usr/include/hdf5/serial/
LDFLAGS=-g -Wall -L/usr/lib/x86_64-linux-gnu/hdf5/serial
LDLIBS=-lm -lgsl -lgslcblas -lhdf5

SRCS=src/main.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: $(OBJS)
	g++ $(LDFLAGS) -o run src/main.o $(LDLIBS) 

main.o: src/main.cpp includes/common.h
	g++ $(CPPFLAGS) -c src/main.cpp
