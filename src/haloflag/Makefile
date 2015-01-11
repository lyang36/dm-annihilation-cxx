all:
	nvcc haloflag_largemem.cu -o haloflag_l
	nvcc haloflag_smallmem.cu -o haloflag_s
	nvcc gethaloflags.cu -o gethaloflags


clean:
	rm -fr haloflag_l
	rm -fr haloflag_s
	rm -fr gethaloflags
