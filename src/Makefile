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

NVCCARGS=-arch=compute_13

C=gcc
CC=g++
OBJS=main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
	 genmap.o anglefuncs.o

LFLAGS=$(HEALPIXLIB)/libhealpix_cxx.a\
	   $(HEALPIXLIB)/libcxxsupport.a\
	   $(FITSIOLIB)/libcfitsio.a

IFLAG=-I $(HEALPIXSRC) -I$(FITSIOLIB)



all: ann-nosom ann-somv ann-somv2 ann-pwave ann-core-nosom decay others

others:
	@cd dataconverter; make
	@cd test; make

comp:
	$(CC) -c main.cpp  datareader.cpp  tipsy_io.c \
		   	parameters.cpp mapgenerator.cpp genmap.cpp\
			anglefuncs.cpp $(IFLAG)

ann-nosom: comp
	$(CC) -c $(FLUXDIR)/flux_ann_nosom.cpp
	$(CC) -o annmap_nosom flux_ann_nosom.o $(OBJS) $(LFLAGS)

ann-core-nosom: comp
	$(CC) -c $(FLUXDIR)/flux_ann_core_nosom.cpp $(COREDIR)/halocore.cpp
	$(CC) -o annmap_core_nosom flux_ann_core_nosom.o halocore.o $(OBJS) $(LFLAGS)

ann-somv: comp
	$(CC) -c $(FLUXDIR)/flux_ann_somv.cpp
	$(CC) -o annmap_somv flux_ann_somv.o $(OBJS) $(LFLAGS)
	
ann-somv2: comp
	$(CC) -c $(FLUXDIR)/flux_ann_somv2.cpp
	$(CC) -o annmap_somv2 flux_ann_somv2.o $(OBJS) $(LFLAGS)

ann-pwave: comp
	$(CC) -c $(FLUXDIR)/flux_ann_pwave.cpp
	$(CC) -o annmap_pwave flux_ann_pwave.o $(OBJS) $(LFLAGS)
	
decay: comp
	$(CC) -c $(FLUXDIR)/flux_decay.cpp
	$(CC) -o decaymap_nocorr flux_decay.o $(OBJS) $(LFLAGS)




cuda-all: cuda-ann-nosom cuda-ann-somv cuda-ann-somv2 cuda-ann-pwave cuda-decay cuda-ann-core-nosom

cuda-ann-nosom:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_nosom main.o flux_ann_nosom.o datareader.o\
		   	tipsy_io.o parameters.o mapgenerator.o\
			genmap_cuda.o anglefuncs.o kernel.o $(LFLAGS)

cuda-ann-core-nosom:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda_core.cpp $(COREDIR)/cudahalocore.cu $(IFLAG)
	nvcc -o annmap_cuda_core_nosom main.o \
			halocore.o flux_ann_nosom.o cudahalocore.o\
			datareader.o tipsy_io.o parameters.o \
			mapgenerator.o genmap_cuda_core.o \
			anglefuncs.o kernel.o $(LFLAGS)	

cuda-ann-somv:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_somv main.o flux_ann_somv.o datareader.o\
		   	tipsy_io.o parameters.o mapgenerator.o\
			genmap_cuda.o anglefuncs.o kernel.o $(LFLAGS)

cuda-ann-somv2:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_somv2 main.o flux_ann_somv2.o datareader.o\
		   	tipsy_io.o parameters.o mapgenerator.o\
			genmap_cuda.o anglefuncs.o kernel.o $(LFLAGS)

cuda-ann-pwave:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_pwave main.o flux_ann_pwave.o datareader.o\
		   	tipsy_io.o parameters.o mapgenerator.o\
			genmap_cuda.o anglefuncs.o kernel.o $(LFLAGS)


cuda-decay:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o decaymap_cuda_nocorr main.o flux_decay.o datareader.o\
		   	tipsy_io.o parameters.o mapgenerator.o\
			genmap_cuda.o anglefuncs.o kernel.o $(LFLAGS)


cuda-debug:
	nvcc -c -arch=sm_13 -g -G kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_debug main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
		genmap_cuda.o anglefuncs.o -g -G kernel.o flux_ann_nosom.o $(LFLAGS)	
	nvcc -c -arch=sm_13 -g -G kernel.cu -g -G genmap_cuda_core.cpp \
			-g -G $(COREDIR)/cudahalocore.cu $(IFLAG)
	nvcc -o annmap_cuda_core_nosom_debug main.o \
			halocore.o flux_ann_nosom.o cudahalocore.o\
			datareader.o tipsy_io.o parameters.o \
			mapgenerator.o genmap_cuda_core.o \
			anglefuncs.o -g -G kernel.o $(LFLAGS)

clean:
	rm -fr *.o
	rm -fr annmap_nosom
	rm -fr annmap_somv
	rm -fr annmap_somv2
	rm -fr annmap_pwave
	rm -fr annmap_cuda_nosom
	rm -fr annmap_cuda_somv
	rm -fr annmap_cuda_somv2
	rm -fr annmap_cuda_core_nosom
	rm -fr annmap_core_nosom
	rm -fr annmap_cuda_pwave
	rm -fr decaymap_nocorr
	rm -fr decaymap_cuda_nocorr
	rm -fr annmap_cuda_debug
	rm -fr annmap_cuda_core_nosom_debug
	@cd dataconverter; make clean
	@cd test; make clean
