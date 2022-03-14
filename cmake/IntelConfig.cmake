if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)

    if("${CUSTOM_COMPILE_OPTIONS}" STREQUAL "")
        add_compile_options(-assume byterecl)
        add_compile_options(-heap-arrays)
        add_compile_options(-mkl=parallel)
        add_compile_options(-O3)
    else()
        add_compile_options(${CUSTOM_COMPILE_OPTIONS})
    endif()

    if(${FORTRAN_DEBUG_FLAG})
        add_compile_options(-g)
    endif()

    if(NOT SIMD_SET)
        set(SIMD_SET xHost)
    endif()

    set(XCFUN_MAKEFILE Makefile)

endif()