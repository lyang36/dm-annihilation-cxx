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
	$(CC) -o rprof_adi rprof_adi.o ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}
	$(CC) -o rpf_test radialprofile_test.o ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}

comp:
	$(CC) -c  radialprofile_test.cpp angprofile.cpp rprof_adi.cpp radialprofile.cpp mapselect.cpp $(IFLAG)

clean:
	rm -f rprof_adi
	rm -f mapselect
	rm -f rprofile
	rm -f rprof_adi

