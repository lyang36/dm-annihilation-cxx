
CC=g++
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


#LFLAG= -framework Carbon -framework OpenGL -framework GLUT
#LFLAG=-lGL -lglut -lGLEW
LFLAG= -lGL -lglut /home/lyang/Tools/glew/lib/libGLEW.a
IFLAG=-I /usr/include/GL/ -I./ -I../ -I$(HEALPIXSRC) -I$(FITSIOLIB)
#IFLAG= -I /usr/include/GL/ -I./ -I ~/Tools/glew/include/ -I../



SRC=flux.vert flux.frag buffers.cpp render.cpp shaderManager.cpp shaders.cpp shadersrc.cpp

all: ${SRC}
	$(CC) -c buffers.cpp render.cpp shaderManager.cpp shaders.cpp shadersrc.cpp $(IFLAG)


shadersrc.cpp: flux.vert flux.frag
	@echo "Create shadersrc.cpp"
	@echo "#include <rpc/types.h>" > shadersrc.cpp
	@echo "#include <rpc/xdr.h>" >> shadersrc.cpp
	@echo "#include <string>" >> shadersrc.cpp
	@echo "#include \"shaders.h\"" >> shadersrc.cpp
	@echo "const std::string fluxShader::VERTEX_SHADER_SRC = " >> shadersrc.cpp
	@cat flux.vert | sed "s/^/\"/" | sed 's/.$$/\\\n\"/' >> shadersrc.cpp
	@echo ";" >> shadersrc.cpp
	@echo  >> shadersrc.cpp
	@echo "const std::string fluxShader::FRAGMENT_SHADER_SRC = " >> shadersrc.cpp
	@cat flux.frag | sed "s/^/\"/" | sed 's/.$$/\\\n\"/' >> shadersrc.cpp
	@echo ";" >> shadersrc.cpp

clean:
	rm -f shadersrc.cpp
	rm -f *.o
