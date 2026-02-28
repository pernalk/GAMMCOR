module build_info

    implicit none
    character(len=*), parameter :: FC_COMPILER_ID = "@CMAKE_Fortran_COMPILER_ID@"
    character(len=*), parameter :: CXX_COMPILER_ID = "@CMAKE_CXX_COMPILER_ID@"
    character(len=*), parameter :: C_COMPILER_ID = "@CMAKE_C_COMPILER_ID@"
    character(len=*), parameter :: FC_COMPILER = "@CMAKE_Fortran_COMPILER@"
    character(len=*), parameter :: CXX_COMPILER = "@CMAKE_CXX_COMPILER@"
    character(len=*), parameter :: C_COMPILER = "@CMAKE_C_COMPILER@"
    character(len=*), parameter :: COMPILE_OPTIONS = "@COMPILE_OPTIONS_STRING@"
    character(len=*), parameter :: BLAS_LIBRARIES = "@BLAS_LIBRARIES_STRING@"

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
