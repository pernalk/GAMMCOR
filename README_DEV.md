# 🛠️ GammCor Developer Zone
### 📚 [Documentation](#-documentation) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 🧪 [Testing](#-testing) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 🧰 [Toolbox](#-toolbox) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 🔗 [Links](#-links)
---

## 📚 Documentation
* GammCor documentation is available at: https://qchem.gitlab.io/gammcor-doc
* If you would like to contribute, please visit: https://gitlab.com/qchem/gammcor-doc

## 🧪 Testing
### 🔧 Preparation
1. Define _gitos_ in your `~/.ssh/config` file:
```
Host gitos
  HostName 172.20.50.110
  User git
  Port 22
  IdentityFile ~/.ssh/<public_key_name>
```
* If you do not have a key, please contact the server administrator. 
2. Obtain test input files from repository:
```
git clone gitos:qchem/TESTS
```
3. Open _gammcor_verify.bash_ in the main GammCor directory and assign variable ```TEST_SCRIPT_PATH``` the path to the file ```run_test.py``` located in your TESTS directory

### ✋ Manual testing
Run test script located in the main GammCor directory:
```
./gammcor_verify.bash
```

### 🚦 Automatic testing
* ##### Enabling
1. Go to ```gammcor/auto_test``` directory
2. Run ```./install-hooks.bash``` to install test scripts: ```pre-commit.bash```, ```pre-push.bash``` and ```post-merge.bash```
3. From now on, the tests will be performed before each ```commit``` and ```push``` commands as well as after the ```merge``` command

* ##### Skipping
To bypass ```pre-commit``` or ```pre-push``` tests, you can use ```--no-verify``` option, e.g.
```
git commit --no-verify
```

* ##### Disabling
To permanently disable the automatic testing feature, go to ```gammcor/auto_test``` directory and run ```./uninstall-hooks.bash```

* ##### Adjusting the level of testing
To customize the list of tests performed for a given _git_ command, you can change level parameter passed to ```gammcor_verify.bash``` in ```pre-commit.bash```, ```pre-push.bash``` and ```post-merge.bash``` scripts. There are three levels: ```'short'```, ```'long'``` and ```'full'``` (default), which correspond to the levels defined for individual tests listed in ```TESTS/scripts/test_list.py```.
* ```'full'``` will perfom tests with ```'full'```, ```'long'``` and ```'short'``` levels
* ```'long'``` will perfom tests with ```'long'``` and ```'short'``` levels
* ```'short'``` will only perfom tests with ```'short'``` level

## 🧰 Toolbox
##### Profilers
* ☕ [QAWA Fortran-Code-Profiler](https://github.com/SokolAK/QAWA-Fortran-Code-Profiler)
* [MAQAO](http://www.maqao.org/documentation.html)
* [gprof](https://sourceware.org/binutils/docs/gprof/)
* [OProfile](https://oprofile.sourceforge.io/docs/)

## 🔗 Links

##### Intel
* [Intel Fortran Compiler Options](https://software.intel.com/content/www/us/en/develop/documentation/fortran-compiler-oneapi-dev-guide-and-reference/top/compiler-reference/compiler-options/alphabetical-list-of-compiler-options.html)
* [oneAPI Math Kernel Library - Fortran](https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top.html)
