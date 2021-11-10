# TODO
# * detect available SIMD instruction set (AVX512, AVX2, SSE, ...)
# [ ] prepare CMakeLists.txt for XCFun
#     XCFun crashes for j>1 with Intel

cmake_minimum_required(VERSION 3.14)
project(gammcor LANGUAGES Fortran CXX C)

set(XCFUN_PATH "${CMAKE_SOURCE_DIR}/xcfun/")

set(vendors Intel10_64ilp OpenBLAS IBMESSL)
foreach(vendor IN LISTS vendors)
   set(BLA_VENDOR ${vendor})
   find_package(BLAS QUIET)
   message(STATUS "Looking for ${vendor}")
   if(${BLAS_FOUND})
       break()
   endif()
endforeach()

find_package(OpenMP REQUIRED COMPONENTS Fortran)

# detect available SIMD instruction set
include(CheckFortranCompilerFlag)
check_Fortran_compiler_flag("-march=native" _march_native_works)
check_Fortran_compiler_flag("-xHost" _xhost_works)
set(_Fortran_FLAGS)
if(_march_native_works)
    message(STATUS "Using processor's vector instructions (-march=native compiler flag set)")
    set(_Fortran_FLAGS "-march=native")
elseif(_xhost_works)
    message(STATUS "Using processor's vector instructions (-xHost compiler flag set)")
    set(_Fortran_FLAGS "-xHost")
else()
    message(STATUS "No suitable compiler flag found for vectorization")
endif()
message("Flags" ${_Fortran_FLAGS})

# WIP
#set(march)
#if(NOT DEFINED march)
#    # set(march xCORE-AVX2)
#    set(_Fortran_FLAGS "-xHost")
#    message(STATUS "not defined XXX")
#else()
#    set(_Fortran_FLAGS "-${march}")
#    message(STATUS "defined XXX")
#endif()
#message(STATUS "_Fortran_FLAGS: ${_Fortran_FLAGS}")

# INTEL
# ===========================================================================================================
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)

    set(CMAKE_Fortran_FLAGS "-assume byterecl -heap-arrays -mkl=parallel -O3 -g ${_Fortran_FLAGS}")
    set(XCFUN_MAKEFILE "Makefile")
    link_libraries(${BLAS_LIBRARIES})
#    link_libraries(mkl_intel_ilp64)
endif()
# ===========================================================================================================

# GNU
# ===========================================================================================================
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS "-funroll-loops -O3 -g ${_Fortran_FLAGS}")
    set(XCFUN_MAKEFILE "Makefile.gcc")
    link_libraries(${BLAS_LIBRARIES})
#    link_libraries(mkl_gf_lp64)
endif()

# XCFUN
# ===========================================================================================================
include_directories(${XCFUN_PATH}fortran)
link_directories(${XCFUN_PATH}lib)
link_libraries(xcfun)
set(XCFUN_SOURCE
    ${XCFUN_PATH}fortran/xcfun_module.f90
    ${XCFUN_PATH}fortran/xcfun_autogen.f90
)
add_custom_command(OUTPUT ${XCFUN_PATH}fortran/xcfun_autogen.f90 
    COMMAND make clean
    COMMAND make -f ${XCFUN_MAKEFILE}
    WORKING_DIRECTORY ${XCFUN_PATH}
)
# ===========================================================================================================

message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "C compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")


set(GAMMCOR_SOURCE
    SOURCE/abfofo.f90
    SOURCE/abmats.f90
    SOURCE/accas.f
    SOURCE/ac_exact_2el.f
    SOURCE/acfd.f
    SOURCE/ac_fofo.f90
    SOURCE/ac_iter.f
    SOURCE/caspidft.f
    SOURCE/constants.h
    SOURCE/cpdmft.f
    SOURCE/dftacg_pw92c.f
    SOURCE/dftfun_ecerf.f
    SOURCE/dftfun_ecerfpbe.f
    SOURCE/dftfun_exerf.f
    SOURCE/dftfun_exerfpbe.f
    SOURCE/dftgrid.f
    SOURCE/diis.f90
    SOURCE/dmscf.f
    SOURCE/eigenno.f
    SOURCE/ekt.f
    SOURCE/erpa.f
    SOURCE/exact2el.f
    SOURCE/exappr.f90
    SOURCE/exdpino.f90
    SOURCE/exi.f90
    SOURCE/exmisc.f90
    SOURCE/gridmolpro.f90
    SOURCE/hf.f
    SOURCE/initia.f
    SOURCE/inputfill.f90
    SOURCE/interpa.f
    SOURCE/lsd_sr.f
    SOURCE/mainp.f
    SOURCE/matrix.f
    SOURCE/matvec.f
    SOURCE/misc.f
    SOURCE/newton.f
    SOURCE/nonadia.f
    SOURCE/optapsg.f
    SOURCE/optnorb.f
    SOURCE/optocc.f
    SOURCE/projector.f
    SOURCE/sapt_exch.f90
    SOURCE/sapt_main.f90
    SOURCE/sapt_pol.f90
    SOURCE/sapt_utils.f90
    SOURCE/sorter.f90
    SOURCE/srefex.f90
    SOURCE/srlrdynamic.f
    SOURCE/systemdef.f90
    SOURCE/timing.f90
    SOURCE/tran.f90
    SOURCE/types.f90
    SOURCE/vv10.f
    SOURCE/xcfun.f90
)

add_executable(${PROJECT_NAME}
    ${XCFUN_SOURCE}
    ${GAMMCOR_SOURCE}
)

target_link_libraries(gammcor PRIVATE OpenMP::OpenMP_Fortran)

#set_target_properties(${PROJECT_NAME}
#    PROPERTIES ADDITIONAL_CLEAN_FILES "${XCFUN_PATH}fortran/xcfun_autogen.f90"
#)

