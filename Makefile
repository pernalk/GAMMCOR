MKL_ROOT = /opt/intel/mkl/

FCC = ifort -assume byterecl
FFLAGS = -mkl=parallel -heap-arrays -O3 -g -I xcfun/fortran -xCOMMON-AVX512
LIBS = -L$(MKL_ROOT)lib/intel64/ -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core \
-L./xcfun/lib -lxcfun -lopenblas

#FFLAGS = -O3 -g -I xcfun/fortran -funroll-loops
#LIBS   = -L./xcfun/lib -lxcfun -lopenblas
#FCC = gfortran

include Makefile.common
