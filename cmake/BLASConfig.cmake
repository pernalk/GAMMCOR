message(NOTICE "")
message(NOTICE "🧮  BLAS info")

set(BLA_VENDOR ${BLAS_VENDOR})
find_package(BLAS)

if(NOT ${BLAS_FOUND})
    message(FATAL_ERROR "BLAS vendor not found!
    Please specify BLAS vendor in config file or
    using -DBLAS_VENDOR=<blas_vendor_name>")
else()
    message("BLAS libraries:")
    foreach(lib ${BLAS_LIBRARIES})
        message("* ${lib}")
    endforeach()
endif()

link_libraries(${BLAS_LIBRARIES})
string (REPLACE ";" " " BLAS_LIBRARIES_STRING "${BLAS_LIBRARIES}")

message(NOTICE "")