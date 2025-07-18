message(NOTICE "")
message(NOTICE "🧵  OpenMP info")

find_package(OpenMP COMPONENTS CXX Fortran QUIET)

if(OpenMP_Fortran_FOUND)
  message(STATUS "OpenMP for Fortran found")
else()
  message(WARNING "OpenMP for Fortran not found")
endif()

message(NOTICE "")
