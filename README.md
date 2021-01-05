# GAMMCOR
Authors: Kasia Pernal and Michal Hapka

## Instalation
1. Clone the repository
```
git clone https://github.com/pernalk/GAMMCOR.git 
```

2. Change directory to GAMMCOR and create a directory OBJ:
```
chdir GAMMCOR
mkdir OBJ
```

3. Set the path to mkl libraries in `Makefile` and compile using ifort with `make`:
```
make
```
This will create GAMMCOR executable:
```
gammcor
```
4. Run the verification suite:
   1. Obtain test input files
   ```
   git clone https://gitlab.version.fz-juelich.de/trex/GammCor
   ```
   2. Copy TESTS to a main GAMMCORR directory with a compiled code
   3. Run test suit

   ```
   gammcor_verify.py
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


