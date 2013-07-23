#must setup these 2 lib to compile
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
OBJS=main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
	 genmap.o anglefuncs.o

LFLAGS=$(HEALPIXLIB)/libhealpix_cxx.a\
	   $(HEALPIXLIB)/libcxxsupport.a\
	   $(FITSIOLIB)/libcfitsio.a

IFLAG=-I $(HEALPIXSRC) -I$(FITSIOLIB)



all: ann-nosop ann-sopv ann-sopv2 ann-core-nosop decay others

others:
	@cd dataconverter; make
	@cd test; make

comp:
	$(CC) -c main.cpp  datareader.cpp  tipsy_io.c \
		   	parameters.cpp mapgenerator.cpp genmap.cpp\
			anglefuncs.cpp $(IFLAG)

ann-nosop: comp
	$(CC) -c $(FLUXDIR)/flux_ann_nosop.cpp
	$(CC) -o annmap_nosop flux_ann_nosop.o $(OBJS) $(LFLAGS)

ann-core-nosop: comp
	$(CC) -c $(FLUXDIR)/flux_ann_core_nosop.cpp $(COREDIR)/halocore.cpp
	$(CC) -o annmap_core_nosop flux_ann_core_nosop.o halocore.o $(OBJS) $(LFLAGS)

ann-sopv: comp
	$(CC) -c $(FLUXDIR)/flux_ann_sopv.cpp
	$(CC) -o annmap_sopv flux_ann_sopv.o $(OBJS) $(LFLAGS)
	
ann-sopv2: comp
	$(CC) -c $(FLUXDIR)/flux_ann_sopv2.cpp
	$(CC) -o annmap_sopv2 flux_ann_sopv2.o $(OBJS) $(LFLAGS)
	
decay: comp
	$(CC) -c $(FLUXDIR)/flux_decay.cpp
	$(CC) -o decaymap_nocorr flux_decay.o $(OBJS) $(LFLAGS)




cuda-all: cuda-ann-nosop cuda-ann-sopv cuda-ann-sopv2 cuda-decay cuda-ann-core-nosop

cuda-ann-nosop:
	nvcc -c -arch=sm_13 kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_nosop main.o flux_ann_nosop.o datareader.o\
		   	tipsy_io.o parameters.o mapgenerator.o\
			genmap_cuda.o anglefuncs.o kernel.o $(LFLAGS)

cuda-ann-core-nosop:
	nvcc -c -arch=sm_13 kernel.cu genmap_cuda_core.cpp $(COREDIR)/cudahalocore.cu $(IFLAG)
	nvcc -o annmap_cuda_core_nosop main.o \
			halocore.o flux_ann_nosop.o cudahalocore.o\
			datareader.o tipsy_io.o parameters.o \
			mapgenerator.o genmap_cuda_core.o \
			anglefuncs.o kernel.o $(LFLAGS)	

cuda-ann-sopv:
	nvcc -c -arch=sm_13 kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_sopv main.o flux_ann_sopv.o datareader.o\
		   	tipsy_io.o parameters.o mapgenerator.o\
			genmap_cuda.o anglefuncs.o kernel.o $(LFLAGS)

cuda-ann-sopv2:
	nvcc -c -arch=sm_13 kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_sopv2 main.o flux_ann_sopv2.o datareader.o\
		   	tipsy_io.o parameters.o mapgenerator.o\
			genmap_cuda.o anglefuncs.o kernel.o $(LFLAGS)

cuda-decay:
	nvcc -c -arch=sm_13 kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o decaymap_cuda_nocorr main.o flux_decay.o datareader.o\
		   	tipsy_io.o parameters.o mapgenerator.o\
			genmap_cuda.o anglefuncs.o kernel.o $(LFLAGS)


cuda-debug:
	nvcc -c -arch=sm_13 -g -G kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_debug main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
		genmap_cuda.o anglefuncs.o -g -G kernel.o flux_ann_nosop.o $(LFLAGS)	
	nvcc -c -arch=sm_13 -g -G kernel.cu -g -G genmap_cuda_core.cpp \
			-g -G $(COREDIR)/cudahalocore.cu $(IFLAG)
	nvcc -o annmap_cuda_core_nosop_debug main.o \
			halocore.o flux_ann_nosop.o cudahalocore.o\
			datareader.o tipsy_io.o parameters.o \
			mapgenerator.o genmap_cuda_core.o \
			anglefuncs.o -g -G kernel.o $(LFLAGS)

clean:
	rm -fr *.o
	rm -fr annmap_nosop
	rm -fr annmap_sopv
	rm -fr annmap_sopv2
	rm -fr annmap_cuda_nosop
	rm -fr annmap_cuda_sopv
	rm -fr annmap_cuda_sopv2
	rm -fr annmap_cuda_core_nosop
	rm -fr annmap_core_nosop
	rm -fr decaymap_nocorr
	rm -fr decaymap_cuda_nocorr
	rm -fr annmap_cuda_debug
	rm -fr annmap_cuda_core_nosop_debug
	@cd dataconverter; make clean
	@cd test; make clean
