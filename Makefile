MKL_ROOT = /opt/intel/mkl/
O = OBJ/

FCC = ifx
ifeq ($(findstring gfortran,$(FCC)),gfortran)
  MODFLAG = -J$(O)
else
  MODFLAG = -module $(O)
endif
WARNINGS          =     -warn nounused
OPTIMIZATION      =     -xHost -O3
PARALLELIZATION   = -coarray=single -qmkl=parallel -qopenmp
XCFUN                   =     -I xcfun/fortran
CHOLESKY       =  -I ./gammcor-integrals/include

HDF5_ROOT      =  $(CURDIR)/hdf5
HDF5           =  -I $(HDF5_ROOT)/include

COMMON_FLAGS = $(WARNINGS) $(PARALLELIZATION) $(XCFUN) $(CHOLESKY) $(HDF5) -assume byterecl -heap-arrays -fpp
FFLAGS = $(COMMON_FLAGS) $(OPTIMIZATION) 
DEBUG_FLAGS = -O0 -traceback -check all -DDEBUG

# MKL_LIB is not needed when -qmkl=parallel is used
#MKL_LIB           =     -L$(MKL_ROOT)lib/intel64/ -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
#
XCFUN_LIB   =     -L./xcfun/lib/ -lxcfun
CHOLESKY_LIB = ./gammcor-integrals/lib/cholesky.a
HDF5_LIB = -L$(HDF5_ROOT)/lib -Wl,-rpath,$(HDF5_ROOT)/lib -lhdf5_hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lsz -ldl -lm

#LIBS = $(MKL_LIB) $(XCFUN_LIB) $(CHOLESKY_LIB) -limf 

LIBS = $(XCFUN_LIB) $(CHOLESKY_LIB) $(HDF5_LIB) -limf 

include Makefile.common

.PHONY: all debug clean

all: $(PROG)

debug: FFLAGS = $(COMMON_FLAGS) $(DEBUG_FLAGS)
debug: all
