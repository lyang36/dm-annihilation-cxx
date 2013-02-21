C=gcc
CC=g++

all: comp
	#$(CC) -o annmap main.o
	@cd dataconverter; make
	@cd test; make

comp:
	$(CC) -c main.cpp
	$(CC) -c datareader.cpp
	$(CC) -c readfiles.cpp
	$(CC) -c tipsy_io.c

clean:
	rm -fr *.o
	rm -fr annmap
	@cd dataconverter; make clean
	@cd test; make clean
