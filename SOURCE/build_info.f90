module build_info

    implicit none
    character(len=*), parameter :: FC_COMPILER_ID = "Intel"
    character(len=*), parameter :: CXX_COMPILER_ID = "Intel"
    character(len=*), parameter :: C_COMPILER_ID = "Intel"
    character(len=*), parameter :: FC_COMPILER = "/opt/intel/oneapi/compiler/2021.1.2/linux/bin/intel64/ifort"
    character(len=*), parameter :: CXX_COMPILER = "/opt/intel/oneapi/compiler/2021.1.2/linux/bin/intel64/icpc"
    character(len=*), parameter :: C_COMPILER = "/opt/intel/oneapi/compiler/2021.1.2/linux/bin/intel64/icc"
    character(len=*), parameter :: COMPILE_OPTIONS = "-assume byterecl -heap-arrays -O3 -g -xHost"
    character(len=*), parameter :: BLAS_LIBRARIES = "/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.so /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_thread.so /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.so /usr/lib/x86_64-linux-gnu/libiomp5.so -lpthread -lm -ldl"

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
