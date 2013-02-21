OBJS=../datareader.o ../tipsy_io.o ../parameters.o

all:
	g++ -c testdatareader.cpp
	g++ -c testparams.cpp
	g++ -o testdatareader testdatareader.o $(OBJS)
	g++ -o testparams testparams.o $(OBJS)

clean:
	rm -fr *.o
	rm -fr testdatareader
	rm -fr testparams
