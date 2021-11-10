# GAMMCOR
Authors: Kasia Pernal and Michal Hapka

## Installation
#### 1. Clone the repository
* Using ssh:
```
git clone git@gitlab.com:qchem/gammcor.git
```
* Using https:
```
git clone https://gitlab.com/qchem/gammcor.git
```

#### 2. Create an OBJ directory inside the repository folder:
```
cd <repository_name>
mkdir OBJ
```

#### CMake build: switch to 'cmake' branch!
#### 3a. Build with CMake (IFort example)
env FC=ifort CC=gcc CXX=icc cmake -S. -Bbuild

#### 3b. Build with CMake (GNU example)
env FC=gfortran CC=gcc CXX=gcc cmake -S. -Bbuild

#### 3. Build XCFun
* ##### Using Intel compilers (icc/ipcp & ifort):
```
cd xcfun
make
```
* ##### Using GCC (gcc/g++ & gfortran):
```
cd xcfun
make -f Makefile.gcc
```
#### 4. Build GammCor
* ##### Using ifort:

Set the path to MKL (MKL_ROOT) in `Makefile` and build GammCor using:
```
cd ..
make
```
* ##### Using gfortran:

Set the path to MKL (MKL_ROOT) in `Makefile.gcc` and build GammCor using:
```
cd ..
make -f Makefile.gcc
```

This will create GammCor executable:
```
gammcor
```

## Third party software
Third party software used in GammCor:
#### Intel Math Kernel Library
Intel-optimized linear algebra library with low-level routines that operate on vectors and matrices
* Website: [https://software.intel.com](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html)
#### XCFun
XCFun is a library of exchange-correlation functionals with arbitrary-order derivatives.
* Website: https://github.com/dftlibs/xcfun
* Licence: XCFun is licensed under version 2.0 of the Mozilla Public License (MPLv2.0).
* Reference:
_Ulf Ekström, Lucas Visscher, Radovan Bast, Andreas J. Thorvaldsen and Kenneth Ruud, 
Arbitrary-Order Density Functional Response Theory from Automatic Differentiation, 
Journal of Chemical Theory and Computation 6, 1971 (2010), DOI: 10.1021/ct100117s_.
