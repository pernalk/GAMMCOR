#MKL_ROOT = /opt/intel/mkl/

#FCC = ifort -assume byterecl
#FFLAGS = -mkl=parallel -heap-arrays -O3 -g -I xcfun/fortran -xCOMMON-AVX512
#LIBS = -L$(MKL_ROOT)lib/intel64/ -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core \
#-L./xcfun/lib -lxcfun -lopenblas

FFLAGS = -O3 -Ixcfun/fortran
LIBS   = -L./xcfun/lib -lxcfun -L/home/mitek/programy/ATLAS-3.11.39/build/lib -llapack -lcblas -lf77blas -latlas
FCC = gfortran

include Makefile.common