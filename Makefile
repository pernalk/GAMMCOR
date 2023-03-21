MKL_ROOT = /opt/intel/mkl/

FCC = ifort 
WARNINGS 		=	-warn nounused
OPTIMIZATION 	=	-xCOMMON-AVX512 -O3
#OPTIMIZATION 	=	-xHost -O3
PARALLELIZATION	=	-mkl=parallel -qopenmp
XCFUN 			=	-I xcfun/fortran
TREXIO         = $(shell pkg-config --cflags trexio)
FFLAGS = $(WARNINGS) $(OPTIMIZATION) $(PARALLELIZATION) $(XCFUN) $(TREXIO) -assume byterecl -heap-arrays -g

MKL_LIB 		=	-L$(MKL_ROOT)lib/intel64/ -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
XCFUN_LIB 	=	-L./xcfun/lib/ -lxcfun
TREXIO_LIB = $(shell pkg-config --libs trexio) -Wl, -Xlinker="-rpath=$(shell pkg-config --variable=libdir trexio)"
LIBS = $(MKL_LIB) $(XCFUN_LIB) $(TREXIO_LIB) -limf

include Makefile.common
