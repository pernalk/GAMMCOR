if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)

    if("${CUSTOM_COMPILE_OPTIONS}" STREQUAL "")
        add_compile_options(-funroll-loops)
        add_compile_options(-O3)
    else()
        add_compile_options(${CUSTOM_COMPILE_OPTIONS})
    endif()

    if(${FORTRAN_DEBUG_FLAG})
        add_compile_options(-g)
    endif()

    if(NOT DEFINED SIMD_SET)
        set(SIMD_SET march=native)
    endif()

    set(XCFUN_MAKEFILE Makefile.gcc)
    
endif()