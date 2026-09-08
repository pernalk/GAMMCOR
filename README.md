# GAMMCOR
Authors: Kasia Pernal and Michal Hapka<br>
User Manual: [link](https://qchem.gitlab.io/gammcor-manual/)


## Installation
* Using CMake: [link](https://qchem.gitlab.io/gammcor-manual/pages/introduction/cmake.html)
* Using ready-made Makefiles: [link](https://qchem.gitlab.io/gammcor-manual/pages/introduction/makefiles.html)


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

#### HDF5
Used to read the PySCF interface file (`pyscf_data.h5`). The compiled libraries
live in `hdf5/` and are committed; rebuild them with `hdf5/build_hdf5.sh` only
if your compiler cannot read them. See `hdf5/README.md`.
* Website: https://www.hdfgroup.org/solutions/hdf5/
* Licence: BSD-style.

#### TREXIO
* Website: [TREXIO repository](https://github.com/TREX-CoE/trexio)
