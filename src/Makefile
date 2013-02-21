C=gcc
CC=g++
OBJS=main.o datareader.o tipsy_io.o parameters.o mapgenerator.o

all: comp
	$(CC) -o annmap $(OBJS) 
	@cd dataconverter; make
	@cd test; make

comp:
	$(CC) -c main.cpp
	$(CC) -c datareader.cpp
	$(CC) -c readfiles.cpp
	$(CC) -c tipsy_io.c
	$(CC) -c parameters.cpp
	$(CC) -c mapgenerator.cpp

clean:
	rm -fr *.o
	rm -fr annmap
	@cd dataconverter; make clean
	@cd test; make clean
