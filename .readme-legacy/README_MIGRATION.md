# Migration from ifort to ifx Compiler

## Overview

Intel's ifort compiler is being deprecated and replaced by ifx, the new LLVM-based Fortran compiler. This migration updates the gammcor codebase to use ifx while maintaining full compatibility with all Fortran standards (F77 through F2018/2023).

## Changes Made

### Compiler Updates

1. **Main gammcor Makefile**
   - Changed `FCC = ifort` to `FCC = ifx`
   - Updated MKL flag from `-mkl` to `-qmkl` (new syntax for Intel oneAPI)

2. **xcfun Library**
   - `xcfun/Makefile`:
     - Changed `FC = ifort` to `FC = ifx`
     - Removed icc-specific warning flags: `-wd981 -wd279 -wd383 -wd1572 -wd177`

3. **System Path Updates** (`xcfun/Makefile.common`)
   - Updated hardcoded GCC 4.9 paths to current system paths after dragon cluster update:
   ```makefile
   # Old:
   PATHS:=-Iinclude -Isrc -Isrc/taylor -Isrc/functionals -Llib -I/usr/include/x86_64-linux-gnu/c++/4.9/
   
   # New:
   PATHS:=-Iinclude -Isrc -Isrc/taylor -Isrc/functionals -Llib -isystem /usr/include/c++/11 -isystem /usr/include/x86_64-linux-gnu/c++/11 -isystem /usr/include/x86_64-linux-gnu

## Code Fixes

#### projector.f - Fixed parameter naming conflict detected by ifx:

Changed function RtBis to take integer argument ISum=1 or 2 instead of taking function XSum=XSum1 or XSum2 

##### RtBis(XSum,...) to RtBis(ISum,...)

Added conditional logic to call appropriate XSum function:

if (ISum.eq.1) then
    fmid=XSum1(XNorm,x2,Occ,NB)
    f=XSum1(XNorm,x1,Occ,NB)
else if (ISum.eq.2) then
    fmid=XSum2(XNorm,x2,Occ,NB)
    f=XSum2(XNorm,x1,Occ,NB)
end if

## Required Actions After Pulling

1. #### Rebuild xcfun Library

  bash

  ```
  cd xcfun
  make clean
  make
  cd ..
  ```

  

2. #### Rebuild gammcor-integrals Library

  The gammcor-integrals library must be recompiled with ifx:

  ```
  cd gammcor-integrals
  cd CompilerFlags
  cp -r ifort-gammcor ifx-gammcor
  cd ifx-gammcor
  # modfiy ifort ->ifx   in linker and compiler
  # modify -mkl -> -qmkl in linker and compiler
  cd ../..
  ./Build.py --clean
  ./Build.py -np 4 ifx-gammcor
  ```

  

3. #### Rebuild Main Code

  ```
  make clean
  make
  ```

  
