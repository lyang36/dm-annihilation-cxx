OBJS=../datareader.o ../tipsy_io.o

all:
	g++ -c testdatareader.cpp
	g++ -o testdatareader testdatareader.o $(OBJS)

clean:
	rm -fr *.o
	rm -fr testdatareader
