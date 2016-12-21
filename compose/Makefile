# Makefile for compiling Compose
# ==============================
#
#
#
# Options:
# --------
#
#
# BOUNDCHECK=(0|1)
#
# Note: Switch on bound checking for Fortran part.
#
# HDF5=(0|1)
# 
# Enable writing in HDF5-format for the data tables



NAME = compose
EXEC = $(NAME)

# default flag settings (BOUNDCHECK=0, HDF5 = 0)

BOUNDCHECK = 0
HDF5 = 1
ifeq ($(BOUNDCHECK),1)
   FC_FLAGS_BOUNDCHECK = -fbounds-check
endif	

ifeq ($(HDF5),1)
   HDF5_LIB = -lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial
   HDF5_C = -Dhdf5	
else   
   HDF5_LIB = 
   HDF5_C = 
endif
FC = gfortran
FC_FLAGS = -c -g  -ffree-form $(FC_FLAGS_BOUNDCHECK) 
CC = gcc
CC_FLAGS = -c -g -I/usr/include/hdf5/serial/

LINK = gfortran
LD_LIB = $(HDF5_LIB)

ifeq ($(HDF5),1)
  SRC_F = composemodules.f90 \
          hdf5compose.f90 \
          compose.f90 
  SRC_C = hdf5writecompose.c \
	  hdf5readcompose.c
else
  SRC_F = composemodules.f90 \
          compose.f90
  SRC_C =  
endif


OBJ_C := $(SRC_C:.c=.o)
OBJ_F := $(SRC_F:.f90=.o)

$(NAME):  $(OBJ_F) $(OBJ_C) 
	rm -f $(EXEC);
	@echo building compose;
	$(LINK) -o $(EXEC) $(OBJ_F) $(OBJ_C) $(LD_LIB)

compose.o : compose.f90 
	cpp $(HDF5_C) compose.f90 >compose.for;
	$(FC) $(FC_FLAGS) compose.for
composemodules.o : composemodules.f90 
	$(FC) $(FC_FLAGS) composemodules.f90
hdf5compose.o : hdf5compose.f90 
	$(FC) $(FC_FLAGS) hdf5compose.f90
hdf5writecompose.o : hdf5writecompose.c
	$(CC) $(CC_FLAGS) hdf5writecompose.c
hdf5readcompose.o : hdf5readcompose.c
	$(CC) $(CC_FLAGS) hdf5readcompose.c


clean:
	rm -f *.o
