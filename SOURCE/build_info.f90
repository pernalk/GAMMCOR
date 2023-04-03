module build_info

    implicit none
    character(len=*), parameter :: FC_COMPILER_ID = "GNU"
    character(len=*), parameter :: CXX_COMPILER_ID = "GNU"
    character(len=*), parameter :: C_COMPILER_ID = "GNU"
    character(len=*), parameter :: FC_COMPILER = "/usr/bin/gfortran"
    character(len=*), parameter :: CXX_COMPILER = "/usr/bin/g++"
    character(len=*), parameter :: C_COMPILER = "/usr/bin/gcc"
    character(len=*), parameter :: COMPILE_OPTIONS = "-funroll-loops -O3 -ffree-line-length-none -fallow-argument-mismatch -march=native"
    character(len=*), parameter :: BLAS_LIBRARIES = "/usr/lib/x86_64-linux-gnu/libopenblas.so"

    contains

    subroutine build_print_info()
        write(*,*) ""
        write(*,*) "BUILD INFO"
        write(*,'(a)') "********************************************************************************"
        write(*,*) "Fortran compiler:   ", FC_COMPILER_ID, " ", FC_COMPILER
        write(*,*) "C++ compiler:       ", CXX_COMPILER_ID, " ", CXX_COMPILER
        write(*,*) "C compiler:         ", C_COMPILER_ID, " ", C_COMPILER
        write(*,*) "Compile options:    ", COMPILE_OPTIONS
        write(*,*) "BLAS libraries:     ", BLAS_LIBRARIES
        write(*,*) ""
    end subroutine

end module build_info
