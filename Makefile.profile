MKL_ROOT = /home/hemik/miniconda3/

FCC = gfortran
FFLAGS = -O2 -g -fno-omit-frame-pointer -march=native -fopenmp -fallow-argument-mismatch -I xcfun/fortran

MKL_LIB = -L$(MKL_ROOT)lib/ -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
XCFUN_LIB = -L./xcfun/lib/ -lxcfun

LIBS = $(MKL_LIB) $(XCFUN_LIB)

include Makefile.common
