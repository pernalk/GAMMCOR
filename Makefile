MKL_ROOT = /opt/intel/mkl/

FCC = ifort
WARNINGS          =     -warn nounused
OPTIMIZATION      =     -xHost -O3
PARALLELIZATION   = -coarray=single -mkl=parallel -qopenmp
XCFUN                   =     -I xcfun/fortran
CHOLESKY       =  -I ./gammcor-integrals/include

COMMON_FLAGS = $(WARNINGS) $(PARALLELIZATION) $(XCFUN) $(CHOLESKY) -assume byterecl -heap-arrays -g
FFLAGS = $(COMMON_FLAGS) $(OPTIMIZATION) 
DEBUG_FLAGS = -O0 -traceback -check all

MKL_LIB           =     -L$(MKL_ROOT)lib/intel64/ -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
XCFUN_LIB   =     -L./xcfun/lib/ -lxcfun
CHOLESKY_LIB = ./gammcor-integrals/lib/cholesky.a

LIBS = $(MKL_LIB) $(XCFUN_LIB) $(CHOLESKY_LIB) -limf 

include Makefile.common

.PHONY: all debug clean

all: $(PROG)

debug: FFLAGS = $(COMMON_FLAGS) $(DEBUG_FLAGS)
debug: all
