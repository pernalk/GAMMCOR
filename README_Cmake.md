# Compilation



   
1. Create a build directory:

   ```shell
       mkdir build
       cd build
   ```
    
2. Configure the project with CMake. Check out the additional [configuration options](#cmake-configuration-options).

    ```shell
    cmake .. -DXCFUN_PATH=$XC_FUN_PATH -DCMAKE_Fortran_COMPILER=$COMPILER
    ```

   *  ifort exapmle for Dragon.
   
       ```shell
       cmake .. -DXCFUN_PATH=/home/pkowalski/xcfun_intel -DCMAKE_Fortran_COMPILER=ifort
       ```

3. Build the project using [CMake Build Tool Mode]. This is a portable variant of `make`.

    ```shell
    make 
    ```
   
 4. Mixing Fortran with C
    
    Usefull links:
 
    1. http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
    
    2. http://www.bndhep.net/Software/Mixing/Mixing_Manual.html
    
    3. https://www.cae.tntech.edu/help/programming/mixed_languages
    
    4. http://indico.ictp.it/event/a13229/session/2/contribution/11/material/0/0.pdf
    
    5. http://www.math.chalmers.se/Math/Grundutb/CTH/tma881/0607/Assignments/Mixing_C_Fortran.html
    
    6. https://docs.oracle.com/cd/E19422-01/819-3685/11_cfort.html
    
    7. https://pl.qaru.tech/questions/25127493/calling-a-c-function-from-a-c-function-called-by-a-fortran-code
    
    Other intresting
    
    1. https://sites.google.com/site/akohlmey/
    
    2. https://sites.google.com/site/akohlmey/lectures/introduction-to-hpc-fall-2010
    
    
    
    