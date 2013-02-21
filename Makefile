all: annmap

fitsio:
	@cd libs/cfitsio; make

healpix:
	@cd libs/Healpix; make

annmap:
	@cd src; make
	@mkdir -p bin
	@cp src/annmap bin/annmap
	@cp src/VL2_info.ini bin/.
	@cp src/dataconverter/dataconverter bin/.

clean:
	@cd src; make clean
	rm -fr bin/annmap
	rm -fr bin/dataconverter

clean-all:
	cd libs/cfitsio; make clean
	cd libs/Healpix; make clean
