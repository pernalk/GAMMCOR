cmake_minimum_required(VERSION 3.16)
project(gammcor LANGUAGES Fortran CXX C)

# Set CMake paths
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
)

# Services
include(MessageService)
include(SettingsService)
include(CustomTargets)
include(AddCompileOptionsService)

# BLAS configuration
include(BLASConfig)

# OpenMP configuration
include(OpenMPConfig)

# Compilers configuration
include(IntelConfig)
include(GNUConfig)
include(CommonFortranConfig)

# Libraries configuration
include(XCFunConfig)

# Source files list
include(SourceFilesList)

# Info services
include(GitInfoService)
include(BuildInfoService)

# Create executable
add_executable(${PROJECT_NAME}
    ${XCFUN_SOURCE}
    ${GAMMCOR_SOURCE}
)

# Copy executable to main directory
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND cp -f ${CMAKE_SOURCE_DIR}/build/gammcor ${CMAKE_SOURCE_DIR}/gammcor
)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_Fortran)
