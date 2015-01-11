#must setup these 2 lib to compile
#FITSIOLIB=/Users/lyang/Tools/cfitsio
#HEALPIXHOME=/Users/lyang/Tools/Healpix_3.00
FITSIOLIB=/home/lyang/Tools/cfitsio
HEALPIXHOME=/home/lyang/Tools/Healpix_3.00
#FITSIOLIB=/scratch0/hhpc-hn7/lyang/Tools/cfitsio
#HEALPIXHOME=/scratch0/hhpc-hn7/lyang/Tools/Healpix_3.00

HEALPIXSRC = $(HEALPIXHOME)/src/cxx/basic_gcc/include
HEALPIXLIB = $(HEALPIXHOME)/src/cxx/basic_gcc/lib
HEALPIXCXX = $(HEALPIXHOME)/src/cxx/cxxsupport
HEALPIXC  = $(HEALPIXHOME)/cxx/Healpix_cxx

FLUXDIR = ./flux
COREDIR = ./halocore

LFLAGS=$(HEALPIXLIB)/libhealpix_cxx.a\
	   $(HEALPIXLIB)/libcxxsupport.a\
	   $(FITSIOLIB)/libcfitsio.a

IFLAG=-I $(HEALPIXSRC) -I$(FITSIOLIB)


#OGLLFLAG= -framework Carbon -framework OpenGL -framework GLUT /usr/lib/libGLEW.a
OGLLFLAG=-lGL -lglut -lGLEW
#OGLLFLAG= -lGL -lglut /home/lyang/Tools/glew/lib/libGLEW.a
OGLIFLAG=-I /usr/include/GL/ -I./ -I../ -I$(HEALPIXSRC) -I$(FITSIOLIB)
#OGLIFLAG= -I /usr/include/GL/ -I./ -I ~/Tools/glew/include/ -I../

NVCCARGS=-arch=sm_20

C=gcc -O3 -Wall
CC=g++ -O3 -Wall
OBJS=main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
	 genmap.o anglefuncs.o datatypes.o

OBJS_CUDA=main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
	 genmap_cuda.o anglefuncs.o datatypes.o

OBJS_CUDA_RENDER=main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
	  anglefuncs.o datatypes.o

OBJS_OGL_RENDER=main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
	 genmap_ogl.o anglefuncs.o datatypes.o opengl/buffers.o opengl/render.o\
	 opengl/shaderManager.o opengl/shaders.o opengl/shadersrc.o



all: ann-nosom ann-somv ann-somv2 ann-pwave ann-core-nosom ann-adi decay others

others:
	@cd dataconverter; make
	@cd test; make

comp:
	$(CC) -c main.cpp  datareader.cpp  tipsy_io.c \
		   	parameters.cpp mapgenerator.cpp genmap.cpp\
			anglefuncs.cpp datatypes.cpp $(IFLAG)

ann-nosom: comp
	$(CC) -c $(FLUXDIR)/flux_ann_nosom.cpp
	$(CC) -o annmap_nosom flux_ann_nosom.o $(OBJS) $(LFLAGS)

ann-adi: comp
	$(CC) -c $(FLUXDIR)/flux_ann_adiabatic.cpp
	$(CC) -o annmap_adi flux_ann_adiabatic.o $(OBJS) $(LFLAGS)

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
	@echo "Warning: PLEASE USE THE NEW METHOD!!"

cuda-ann-nosom:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_nosom flux_ann_nosom.o \
			$(OBJS_CUDA) kernel.o $(LFLAGS)

cuda-ann-core-nosom:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda_core.cpp \
			$(COREDIR)/cudahalocore.cu $(IFLAG)
	nvcc -o annmap_cuda_core_nosom  \
			halocore.o flux_ann_nosom.o cudahalocore.o\
			genmap_cuda_core.o kernel.o $(LFLAGS)\
			main.o datareader.o tipsy_io.o parameters.o mapgenerator.o\
	 	    anglefuncs.o datatypes.o	

cuda-ann-somv:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_somv flux_ann_somv.o\
		   	$(OBJS_CUDA) kernel.o $(LFLAGS)

cuda-ann-somv2:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_somv2 flux_ann_somv2.o\
			$(OBJS_CUDA) kernel.o  $(LFLAGS)

cuda-ann-pwave:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o annmap_cuda_pwave flux_ann_pwave.o\
			$(OBJS_CUDA) kernel.o  $(LFLAGS)


cuda-decay:
	nvcc -c $(NVCCARGS) kernel.cu genmap_cuda.cpp $(IFLAG)
	nvcc -o decaymap_cuda_nocorr flux_decay.o\
			$(OBJS_CUDA) kernel.o $(LFLAGS)


cuda-render: cuda-ann-adi-render cuda-ann-nosom-render cuda-ann-somv-render cuda-ann-somv2-render cuda-ann-pwave-render cuda-decay-render

cuda-ann-nosom-render:
	nvcc $(NVCCARGS) -c chealpixrender.cu genmap_cuda_render.cpp $(IFLAG)
	nvcc -o annmap_cuda_nosom flux_ann_nosom.o $(OBJS_CUDA_RENDER)\
			genmap_cuda_render.o chealpixrender.o $(LFLAGS)

cuda-ann-adi-render:
	nvcc $(NVCCARGS) -c chealpixrender.cu genmap_cuda_render.cpp $(IFLAG)
	nvcc -o annmap_cuda_adi flux_ann_adiabatic.o $(OBJS_CUDA_RENDER)\
			genmap_cuda_render.o chealpixrender.o $(LFLAGS)

cuda-ann-somv-render:
	nvcc $(NVCCARGS) -c chealpixrender.cu genmap_cuda_render.cpp $(IFLAG)
	nvcc -o annmap_cuda_somv  flux_ann_somv.o $(OBJS_CUDA_RENDER)\
			genmap_cuda_render.o chealpixrender.o $(LFLAGS)

cuda-ann-somv2-render:
	nvcc $(NVCCARGS) -c chealpixrender.cu genmap_cuda_render.cpp $(IFLAG)
	nvcc -o annmap_cuda_somv2 flux_ann_somv2.o $(OBJS_CUDA_RENDER)\
			genmap_cuda_render.o chealpixrender.o $(LFLAGS)

cuda-ann-pwave-render:
	nvcc $(NVCCARGS) -c chealpixrender.cu genmap_cuda_render.cpp $(IFLAG)
	nvcc -o annmap_cuda_pwave flux_ann_pwave.o $(OBJS_CUDA_RENDER)\
			genmap_cuda_render.o chealpixrender.o $(LFLAGS)

cuda-decay-render:
	nvcc $(NVCCARGS) -c chealpixrender.cu genmap_cuda_render.cpp $(IFLAG)
	nvcc -o decaymap_cuda_nocorr flux_decay.o $(OBJS_CUDA_RENDER)\
			genmap_cuda_render.o chealpixrender.o $(LFLAGS)




ogl-comp:
	@cd opengl; make
	$(CC) -c genmap_ogl.cpp $(IFLAG) $(OGLIFLAG)

ogl-render: ogl-ann-nosom ogl-ann-adi ogl-ann-core-nosom ogl-ann-somv ogl-ann-somv2\
	ogl-ann-pwave ogl-decay

ogl-ann-nosom: ogl-comp
	$(CC) -c $(FLUXDIR)/flux_ann_nosom.cpp
	$(CC) -o annmap_ogl_nosom flux_ann_nosom.o $(OBJS_OGL_RENDER) $(LFLAGS) $(OGLLFLAG)

ogl-ann-adi: ogl-comp
	$(CC) -c $(FLUXDIR)/flux_ann_adiabatic.cpp
	$(CC) -o annmap_ogl_adi flux_ann_adiabatic.o $(OBJS_OGL_RENDER) $(LFLAGS) $(OGLLFLAG)

ogl-ann-core-nosom: ogl-comp
	$(CC) -c $(FLUXDIR)/flux_ann_core_nosom.cpp $(COREDIR)/halocore.cpp
	$(CC) -o annmap_ogl_core_nosom flux_ann_core_nosom.o halocore.o $(OBJS_OGL_RENDER) $(LFLAGS) $(OGLLFLAG)

ogl-ann-somv: ogl-comp
	$(CC) -c $(FLUXDIR)/flux_ann_somv.cpp
	$(CC) -o annmap_ogl_somv flux_ann_somv.o $(OBJS_OGL_RENDER) $(LFLAGS) $(OGLLFLAG)
	
ogl-ann-somv2: ogl-comp
	$(CC) -c $(FLUXDIR)/flux_ann_somv2.cpp
	$(CC) -o annmap_ogl_somv2 flux_ann_somv2.o $(OBJS_OGL_RENDER) $(LFLAGS) $(OGLLFLAG)

ogl-ann-pwave: ogl-comp
	$(CC) -c $(FLUXDIR)/flux_ann_pwave.cpp
	$(CC) -o annmap_ogl_pwave flux_ann_pwave.o $(OBJS_OGL_RENDER) $(LFLAGS) $(OGLLFLAG)
	
ogl-decay: ogl-comp
	$(CC) -c $(FLUXDIR)/flux_decay.cpp
	$(CC) -o decaymap_ogl_nocorr flux_decay.o $(OBJS_OGL_RENDER) $(LFLAGS) $(OGLLFLAG)

	

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
	rm -fr annmap_cuda_adi
	rm -fr annmap_cuda_core_nosom
	rm -fr annmap_core_nosom
	rm -fr annmap_cuda_pwave
	rm -fr decaymap_nocorr
	rm -fr decaymap_cuda_nocorr
	rm -fr annmap_adi
	rm -fr annmap_cuda_debug
	rm -fr annmap_cuda_core_nosom_debug
	rm -fr annmap_ogl_adi annmap_ogl_nosom annmap_ogl_somv\
		annmap_ogl_core_nosom  annmap_ogl_pwave  annmap_ogl_somv2\
		decaymap_ogl_nocorr
	@cd dataconverter; make clean
	@cd opengl; make clean
	@cd test; make clean
