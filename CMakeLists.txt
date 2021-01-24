cmake_minimum_required(VERSION 3.17)
project(GAMMCOR)

enable_language(Fortran)
set(CMAKE_CXX_STANDARD 14)

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)

        set(CMAKE_Fortran_FLAGS
                #" -mkl -heap-arrays -fast -standard-semantics"
                " -mkl -heap-arrays "
                )

endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)

        link_libraries(openblas)
        set(CMAKE_Fortran_FLAGS
                " -O3 -fno-align-commons -fcheck=all -fallow-argument-mismatch "
                )

endif()



include_directories(
        ${XCFUN_PATH}/fortran
)
link_directories(
        ${XCFUN_PATH}/lib )

link_libraries(
        xcfun )

SET(XCFUN_SOURCE
        ${XCFUN_PATH}/fortran/xcfun_module.f90
        ${XCFUN_PATH}/fortran/xcfun_autogen.f90
        )

include_directories(SOURCE)

        SET( SOURCE_FILES
        SOURCE/abfofo.f90
        SOURCE/abmats.f90
        SOURCE/ac_exact_2el.f
        SOURCE/accas.f
        SOURCE/acfd.f
        SOURCE/berf.inc
        SOURCE/big
        SOURCE/caspidft.f
        SOURCE/cdft
        SOURCE/commons.inc
        SOURCE/constants.h
        SOURCE/cpdmft.f
        SOURCE/dberfda.inc
        SOURCE/debug
        SOURCE/dftacg_pw92c.f
        SOURCE/dftfun_ecerf.f
        SOURCE/dftfun_ecerfpbe.f
        SOURCE/dftfun_exerf.f
        SOURCE/dftfun_exerfpbe.f
        SOURCE/dftgrid.f
        SOURCE/dmscf.f
        SOURCE/eigenno.f
        SOURCE/ekt.f
 #       SOURCE/embedd.f
        SOURCE/erpa.f
        SOURCE/exact2el.f
        SOURCE/exdpino.f90
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
        SOURCE/mxrec
        SOURCE/newton.f
        SOURCE/nonadia.f
        SOURCE/optapsg.f
        SOURCE/optnorb.f
        SOURCE/optocc.f
        SOURCE/projector.f
        SOURCE/sapt.f90
        SOURCE/sapt_main.f90
        SOURCE/sorter.f90
        SOURCE/srlrdynamic.f
        SOURCE/systemdef.f90
        SOURCE/tapes
        SOURCE/timing.f90
        SOURCE/tran.f90
        SOURCE/types.f90
        SOURCE/xcfun.f90
        SOURCE/vv10.f
                )

add_executable( ${PROJECT_NAME}
        ${XCFUN_SOURCE}
        ${SOURCE_FILES}
        )