MKL_ROOT = /opt/intel/mkl/

FCC = ifort 
FFLAGS = -assume byterecl -heap-arrays -mkl=parallel -qopenmp -I xcfun/fortran -xCOMMON-AVX512 -O3 -g

MKL_LIB = -L$(MKL_ROOT)lib/intel64/ -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
XCFUN_LIB = -L./xcfun/lib/ -lxcfun

LIBS = $(MKL_LIB) $(XCFUN_LIB) -limf

include Makefile.common