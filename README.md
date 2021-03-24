# GAMMCOR
Authors: Kasia Pernal and Michal Hapka

## Installation
#### 1. Clone the repository
```
git clone git@gitlab.com:michal.hapka/pr-dmft.git
```

#### 2. Change directory to GAMMCOR and create a directory OBJ:
```
cd GAMMCOR
mkdir OBJ
```

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
#### 4. Build GAMMCOR
* ##### Using ifort:

Set the path to MKL (MKL_ROOT) in `Makefile` and build GAMMCOR using:
```
cd ..
make
```
* ##### Using gfortran:

Set the path to MKL (MKL_ROOT) in `Makefile.gcc` and build GAMMCOR using:
```
cd ..
make -f Makefile.gcc
```


This will create GAMMCOR executable:
```
gammcor
```
#### 5. Run the verification suite (optional):
1. Obtain test input files:
```
git clone https://gitlab.version.fz-juelich.de/trex/GammCor
```
2. Copy TESTS to a main GAMMCOR directory with a compiled code
3. Run test suit :

```
gammcor_verify.py
```

## Third party software
Third party software used in GAMMCOR:
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
