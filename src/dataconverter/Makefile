OBJS=../datareader.o ../tipsy_io.o

all:
	g++ -c dataconverter.cpp
	g++ -o dataconverter dataconverter.o $(OBJS)

clean:
	rm -fr *.o
	rm -fr dataconverter
