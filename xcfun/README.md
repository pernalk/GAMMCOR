# NEW README

## Building
* Using Intel compilers (icc/ipcp & ifort):

    `make`

* Using GCC (gcc/g++ & gfortran):

    `make -f Makefile.gcc`


# ORIGINAL README

** Arbitrary order Exchange-Correlation functional library **

Copyright Ulf Ekstrom <uekstrom@gmail.com> and contributors 2009-2010. 
See http://admol.org/xcfun for more information.

The main interface is in include/xcfun.h
(or fortran/xc_fun_module.f90 for Fortran bindings).

Copying:

The library is licensed under the LGPL license version 3, see
COPYING.LESSER for more information.


Configuration: 

Check that XC_MAX_ORDER is defined to the highest order derivatives
you need (and not higher) in src/config.h.  Using a too large value
for XC_MAX_ORDER makes compilation slow and the generated code huge.


Building a debug/development version:

Edit Makefile to set CXX (C++ compiler) and flags and run
make, which will create the library file lib/libxcfun.so.


Building an optimized version:

Edit Makefile and add -DNDEBUG to the compiler flags. Add optimization
compiler options. Make sure your compiler performs inlining (-O3 with
gcc).
