#must setup these 2 lib to compile
#FITSIOLIB=/Users/lyang/Tools/cfitsio
#HEALPIXHOME=/Users/lyang/Tools/Healpix_3.00
FITSIOLIB=/home/lyang/Tools/cfitsio
HEALPIXHOME=/home/lyang/Tools/Healpix_3.00


HEALPIXSRC = $(HEALPIXHOME)/src/cxx/basic_gcc/include
HEALPIXLIB = $(HEALPIXHOME)/src/cxx/basic_gcc/lib
HEALPIXCXX = $(HEALPIXHOME)/src/cxx/cxxsupport
HEALPIXC  = $(HEALPIXHOME)/cxx/Healpix_cxx

C=gcc
CC=g++
OBJS=main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
	 genmap.o anglefuncs.o

LFLAGS=$(HEALPIXLIB)/libhealpix_cxx.a\
	   $(HEALPIXLIB)/libcxxsupport.a\
	   $(FITSIOLIB)/libcfitsio.a

IFLAG=-I $(HEALPIXSRC) -I$(FITSIOLIB)



all: ann-nosop

comp:
	$(CC) -c main.cpp  datareader.cpp  tipsy_io.c \
		   	parameters.cpp mapgenerator.cpp genmap.cpp\
			anglefuncs.cpp $(IFLAG)
ann-nosop: comp
	$(CC) -c flux_ann_nosop.cpp
	$(CC) -o annmap_ann_nosop flux_ann_nosop.o $(OBJS) $(LFLAGS)
	@cd dataconverter; make
	@cd test; make


cuda:
	nvcc -c -arch=sm_13 kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
			genmap_cuda.o anglefuncs.o kernel.o $(LFLAGS)

cuda-debug:
	nvcc -c -arch=sm_13 -g -G kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
		genmap_cuda.o anglefuncs.o -g -G kernel.o $(LFLAGS)	

clean:
	rm -fr *.o
	rm -fr annmap_ann_nosop
	rm -fr annmap_cuda
	@cd dataconverter; make clean
	@cd test; make clean
