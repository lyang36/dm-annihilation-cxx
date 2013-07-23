#FITSIOLIB=/Users/lyang/Tools/cfitsio
#HEALPIXHOME=/Users/lyang/Tools/Healpix_3.00
FITSIOLIB=/home/lyang/Tools/cfitsio
HEALPIXHOME=/home/lyang/Tools/Healpix_3.00


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

IFLAG=-I $(HEALPIXSRC) -I$(FITSIOLIB) -I ../ -I ../halocore/

all: comp
	$(CC) -o angprofile angprofile.o ../datareader.o ../tipsy_io.o ../halocore.o  ${LFLAGS}
	$(CC) -o mapselect mapselect.o ../datareader.o ../tipsy_io.o ../halocore.o  ${LFLAGS}
	$(CC) -o rprofile radialprofile.o ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}

comp:
	$(CC) -c angprofile.cpp radialprofile.cpp mapselect.cpp $(IFLAG)



