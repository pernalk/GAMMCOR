MKL_ROOT = /opt/intel/mkl/

FCC = ifort 
WARNINGS 		=	-warn nounused
OPTIMIZATION 	=	-xCOMMON-AVX512 -O3
PARALLELIZATION	=	-mkl=parallel -qopenmp
XCFUN 			=	-I xcfun/fortran
FFLAGS = $(WARNINGS) $(OPTIMIZATION) $(PARALLELIZATION) $(XCFUN) -assume byterecl -heap-arrays -g

MKL_LIB 		=	-L$(MKL_ROOT)lib/intel64/ -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
XCFUN_LIB 		=	-L./xcfun/lib/ -lxcfun
LIBS = $(MKL_LIB) $(XCFUN_LIB) -limf

include Makefile.common
