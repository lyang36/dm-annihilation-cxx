
all: cicpy

cicpy:
	 swig -python -c++ -o _vl2readwrapper_module.cc vl2read.i
	 python setup.py build_ext --inplace

clean:
	rm -fr *.o *.so
	rm -fr build
	rm -fr _vl2readwrapper_module.cc
	rm -fr VL2read.py	
