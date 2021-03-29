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


## Testing

### Preparation
1. Obtain test input files from repository:
```
git clone git@172.20.50.110:qchem/TESTS
```
2. Open _gammcor_verify.bash_ in the main GAMMCOR directory and assign variable ```TEST_SCRIPT_PATH``` the path to the file ```run_test.py``` located in your TESTS directory

### Manual testing
Run test script located in the main GAMMCOR directory:
```
./gammcor_verify.bash
```

### Automatic testing
* ##### Enabling
1. Go to ```GAMMCOR/auto_test``` directory
2. Run ```./install-hooks.bash``` to install test scripts: ```pre-commit.bash```, ```pre-push.bash``` and ```post-merge.bash```
3. From now on, the tests will be performed before each ```commit``` and ```push``` commands as well as after the ```merge``` command

* ##### Skipping
To bypass ```pre-commit``` or ```pre-push``` tests, you can use ```--no-verify``` option, e.g.
```
git commit --no-verify
```

* ##### Disabling
To permanently disable the automatic testing feature, go to ```GAMMCOR/auto_test``` directory and run ```./uninstall-hooks.bash```

* ##### Adjusting the level of testing
To customize the list of tests performed for a given _git_ command, you can change level parameter passed to ```gammcor_verify.bash``` in ```pre-commit.bash```, ```pre-push.bash``` and ```post-merge.bash``` scripts. There are three levels: ```'short'```, ```'long'``` and ```'full'``` (default), which correspond to the levels defined for individual tests listed in ```TESTS/scripts/test_list.py```.
* ```'full'``` will perfom tests with ```'full'```, ```'long'``` and ```'short'``` levels
* ```'long'``` will perfom tests with ```'long'``` and ```'short'``` levels
* ```'short'``` will only perfom tests with ```'short'``` level


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
