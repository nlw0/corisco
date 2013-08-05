all: corisco_aux.so

corisco_aux.so: corisco_aux_setup.py corisco_aux.pyx corisco_aux_external.o
	python2.7 corisco_aux_setup.py build_ext --inplace --force

corisco_aux_external.o: corisco_aux_external.c
	gcc -O3 -fPIC -msse2 -mfpmath=sse corisco_aux_external.c -c -o corisco_aux_external.o

CORISCO-SNAPSHOT:
	tar jcf corisco-`date +'%Y%m%d'`.tar.bz2 Makefile README corisco.py corisco_aux_setup.py corisco_aux.pyx corisco_aux_external.c
