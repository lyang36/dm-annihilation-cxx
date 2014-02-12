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

IFLAG=-I $(HEALPIXSRC) -I$(FITSIOLIB) -I ../ -I ../halocore/

all: comp
	$(CC) -o angprofile angprofile.o ../datareader.o ../tipsy_io.o ../halocore.o  ${LFLAGS}
	$(CC) -o mapselect mapselect.o ../datareader.o ../tipsy_io.o ../halocore.o  ${LFLAGS}
	$(CC) -o rprofile radialprofile.o ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}
	$(CC) -o rprof_adi rprof_adi.o ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}
	$(CC) -o rpf_shell rpf_shell.o ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}
	$(CC) -o rpf_shell_adi rpf_shell_adi.o ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}
	$(CC) -o bulkvelprof bulkvelprof.o ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}
	$(CC) -o sampleparts sampleparts.o ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}
	$(CC) -o distribution distribution.o ../datareader.o ../tipsy_io.o ../halocore.o ${LFLAGS}

comp:
	$(CC) -c rpf_shell_adi.cpp bulkvelprof.cpp rpf_shell.cpp angprofile.cpp\
		   sampleparts.cpp rprof_adi.cpp radialprofile.cpp mapselect.cpp \
		   distribution.cpp $(IFLAG)

clean:
	rm -f *.o
	rm -f sampleparts
	rm -f rpf_shell_adi
	rm -f rpf_shell
	rm -f rprof_adi
	rm -f mapselect
	rm -f rprofile
	rm -f rprof_adi
	rm -f angprofile
	rm -f bulkvelprof
