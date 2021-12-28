<style>
H1{color:Red !important;}
H2{color:DodgerBlue !important;}
H3{color:SkyBlue !important;}
code{color:Green !important;}
</style>

# GAMMCOR
Authors: Kasia Pernal and Michal Hapka

<br>

## Installation

<br>

### 1. Clone the repository
* Using ssh:
    ```
    git clone git@gitlab.com:qchem/gammcor.git
    ```
* Using https:
    ```
    git clone https://gitlab.com/qchem/gammcor.git
    ```

<br>

### 2. Build GammCor
There are two ways to build GammCor:
* using CMake (recommended)
* using ready-made Makefiles (legacy)

<br>

### 2a. Build with CMake
#### a) Go to the GammCor main directory
#### b) Set building parameters in `config` file (optional)
#### c) Prepare the build system (`-C [config file]` is optional):
```
env FC=<fortran compiler> CXX=<C++ compiler> CC=<C compiler> cmake -S. -Bbuild -C [config file]
```
* Intel example:
    ```
    env FC=ifort CXX=icpc CC=icc cmake -S. -Bbuild -C config_intel
    ```
* GNU example:
    ```
    env FC=gfortran CXX=g++ CC=gcc cmake -S. -Bbuild -C config_gnu
    ```
#### d) Start building process:
```
cmake --build build
```

#### d) Congratulations! `gammcor` executable is in the main directory

<br>

#### Note. If you want to rebuild the project from scratch, delete `build` directory and run:
```
cmake --build build --target clean-build
```

<br>

### 2b. Build with ready-made Makefiles
#### a) Go to the GammCor main directory and create `OBJ` folder:
```
mkdir OBJ
```

#### b) Build XCFun
Go to `xcfun` subdirectory and run:
* for Intel compilers (icc/ipcp & ifort):
    ```
    make
    ```
* for GNU compilers (gcc/g++ & gfortran).    
    ```
    make -f Makefile.gcc
    ```


#### c) Build GammCor:
* Using ifort:

    Go to GammCor main directory, set path to MKL `(MKL_ROOT)` in `Makefile` and run:
    ```
    make
    ```
* Using gfortran:

    Go to GammCor main directory, set path to MKL `(MKL_ROOT)` in `Makefile.gcc` and run:
    ```
    make -f Makefile.gcc
    ```

#### d) Congratulations! `gammcor` executable is in the main directory

<br>

---

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
