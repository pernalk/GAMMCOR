MKL_ROOT = /opt/intel/mkl/

FCC = ifort 
WARNINGS 		=	-warn nounused
OPTIMIZATION 	=	-xHost -O3
PARALLELIZATION	=	-mkl=parallel -qopenmp
XCFUN 			=	-I xcfun/fortran
CHOLESKY       =  -I ./gammcor-cholesky/include
FFLAGS = $(WARNINGS) $(OPTIMIZATION) $(PARALLELIZATION) $(XCFUN) $(CHOLESKY) -assume byterecl -heap-arrays -g

MKL_LIB 		=	-L$(MKL_ROOT)lib/intel64/ -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
XCFUN_LIB 	=	-L./xcfun/lib/ -lxcfun
CHOLESKY_LIB = ./gammcor-cholesky/lib/cholesky.a

LIBS = $(MKL_LIB) $(XCFUN_LIB) $(CHOLESKY_LIB) -limf

include Makefile.common
