
all: cicpy vlshow

cicpy:
	 swig -python -c++ -o _vl2readwrapper_module.cc vl2read.i
	 python setup.py build_ext --inplace

vlshow:
	g++ -c ../datareader.cpp ../tipsy_io.c vl2showdata.cpp
	g++ -o VLShow datareader.o tipsy_io.o vl2showdata.o

clean:
	rm -fr *.o *.so
	rm -fr build
	rm -fr _vl2readwrapper_module.cc
	rm -fr vl2read.py
	rm -fr VLShow	
