#must setup these 2 lib to compile
FITSIOLIB=/Users/lyang/Tools/cfitsio
HEALPIXHOME=/Users/lyang/Tools/Healpix_3.00
#FITSIOLIB=/home/lyang/Tools/cfitsio
#HEALPIXHOME=/home/lyang/Tools/Healpix_3.00


HEALPIXSRC = $(HEALPIXHOME)/src/cxx/basic_gcc/include
HEALPIXLIB = $(HEALPIXHOME)/src/cxx/basic_gcc/lib
HEALPIXCXX = $(HEALPIXHOME)/src/cxx/cxxsupport
HEALPIXC  = $(HEALPIXHOME)/cxx/Healpix_cxx

C=gcc
CC=g++
OBJS=selectdm.o tipsy_io.o 

LFLAGS=$(HEALPIXLIB)/libhealpix_cxx.a\
	   $(HEALPIXLIB)/libcxxsupport.a\
	   $(FITSIOLIB)/libcfitsio.a

IFLAG=-I $(HEALPIXSRC) -I$(FITSIOLIB)



all: comp
	$(CC) -o selectdm $(OBJS)
	$(CC) -o id2mark id2mark.o
	$(CC) -o ascii2bin ascii2bin.o
	$(CC) -o bin2ascii bin2ascii.o
	$(CC) -o showdark tipsy_io.o showdark.o

comp:
	$(CC) -c selectdm.cpp showdark.cpp \
			id2mark.cpp ../src/tipsy_io.c \
			ascii2bin.cpp bin2ascii.cpp	\
			$(IFLAG)

clean:
	rm -fr *.o
	rm -fr selectdm
	rm -fr ascii2bin
	rm -fr id2mark
	rm -fr bin2ascii
	rm -fr showdark
