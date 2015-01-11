all: annmap copy

annmap:
	@cd src; make
	@mkdir -p bin

cuda: cuda-a copy

cuda-a:
	@cd src; make cuda-render

copy:
	cp src/annmap_* bin/.
	cp src/decaymap_* bin/.
	cp src/dataconverter/dataconverter bin/.




clean:
	@cd src; make clean
	rm -fr bin/annmap_*
	rm -fr bin/decaymap_*
	rm -fr bin/dataconverter
