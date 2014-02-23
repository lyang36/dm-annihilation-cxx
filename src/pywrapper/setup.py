from distutils.core import setup, Extension

extension_mod = Extension("_vl2read", ["_vl2readwrapper_module.cc", "../datareader.cpp", "./vl2read.cpp", "../tipsy_io.c"])
setup(name = "vl2read", ext_modules=[extension_mod])
