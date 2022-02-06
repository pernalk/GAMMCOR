set(XCFUN_PATH "${CMAKE_SOURCE_DIR}/xcfun/")

include_directories(${XCFUN_PATH}fortran)
link_directories(${XCFUN_PATH}lib)
link_libraries(xcfun)

set(XCFUN_SOURCE
    ${XCFUN_PATH}fortran/xcfun_module.f90
    ${XCFUN_PATH}fortran/xcfun_autogen.f90
)

add_custom_command(OUTPUT ${XCFUN_PATH}fortran/xcfun_autogen.f90 
    COMMAND $(MAKE) clean
    COMMAND $(MAKE) -f ${XCFUN_MAKEFILE}
    WORKING_DIRECTORY ${XCFUN_PATH}
)