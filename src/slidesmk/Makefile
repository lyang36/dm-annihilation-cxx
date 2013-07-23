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
	$(CC) -o slidemake slidemake.o point.o \
		  ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}
	$(CC) -o imagemake imagemake.o point.o \
		  ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}

comp: 
	$(CC) -c slidemake.cpp ../point.cpp $(IFLAG)
	$(CC) -c imagemake.cpp $(IFLAG)


clean:
	rm -f *.o
	rm -f slidemake

