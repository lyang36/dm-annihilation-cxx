FITSIOLIB=/Users/lyang/Tools/cfitsio
HEALPIXHOME=/Users/lyang/Tools/Healpix_3.00
#FITSIOLIB=/home/lyang/Tools/cfitsio
#HEALPIXHOME=/home/lyang/Tools/Healpix_3.00


HEALPIXSRC = $(HEALPIXHOME)/src/cxx/basic_gcc/include
HEALPIXLIB = $(HEALPIXHOME)/src/cxx/basic_gcc/lib
HEALPIXCXX = $(HEALPIXHOME)/src/cxx/cxxsupport
HEALPIXC  = $(HEALPIXHOME)/cxx/Healpix_cxx

FLUXDIR = ./flux
COREDIR = ./halocore

C=gcc
CC=g++
OBJS = angprofile.o

LFLAGS=$(HEALPIXLIB)/libhealpix_cxx.a\
	   $(HEALPIXLIB)/libcxxsupport.a\
	   $(FITSIOLIB)/libcfitsio.a

IFLAG=-I $(HEALPIXSRC) -I$(FITSIOLIB)

all: comp
	$(CC) -o angprofile ${OBJS} ${LFLAGS}

comp:
	$(CC) -c angprofile.cpp $(IFLAG)



