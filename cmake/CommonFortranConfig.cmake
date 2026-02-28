message(NOTICE "")
message(NOTICE "⚙️   Compiler options")

# string (REPLACE ";" " " COMPILE_OPTIONS_STRING "${COMPILE_OPTIONS}")
# message("Compile options: .............. ${COMPILE_OPTIONS_STRING}")

# set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} "-${SIMD_SET}")
# message("SIMD set: ..................... ${SIMD_SET}")

# message("Debug info: ................... ${FORTRAN_DEBUG_FLAG}")

# add_compile_options(${ADDITIONAL_COMPILE_OPTIONS})
# message("Additional compile options: ... ${ADDITIONAL_COMPILE_OPTIONS}")

add_compile_options("-${SIMD_SET}")
add_compile_options(${ADDITIONAL_COMPILE_OPTIONS})

string (REPLACE ";" " " COMPILE_OPTIONS_STRING "${COMPILE_OPTIONS}")
message("Compile options: .............. ${COMPILE_OPTIONS_STRING}")

message(NOTICE "")
