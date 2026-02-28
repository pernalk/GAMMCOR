message(NOTICE "")
message(NOTICE "🛠️   Settings")

# Options
option(FORTRAN_DEBUG_FLAG   "FORTRAN_DEBUG_FLAG"    OFF)

# Show values
message("SIMD_SET: ..................... ${SIMD_SET}")
message("BLAS_VENDOR: .................. ${BLAS_VENDOR}")
message("FORTRAN_DEBUG_FLAG: ........... ${FORTRAN_DEBUG_FLAG}")
message("ADDITIONAL_COMPILE_OPTIONS: ... ${ADDITIONAL_COMPILE_OPTIONS}")
message("CUSTOM_COMPILE_OPTIONS: ....... ${CUSTOM_COMPILE_OPTIONS}")

message(NOTICE "")
