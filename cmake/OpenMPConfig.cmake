message(NOTICE "")
message(NOTICE "🧵  OpenMP info")

find_package(OpenMP COMPONENTS Fortran)
# target_link_libraries(gammcor PRIVATE OpenMP::OpenMP_Fortran)

message(NOTICE "")