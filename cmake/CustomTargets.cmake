# cmake --build build --target clean-build-verify
# Clean, build and verify on the current node
add_custom_target(clean-build-verify
    COMMAND cmake --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND cmake --build ${CMAKE_BINARY_DIR} -j 8
    COMMAND ../gammcor_verify_cmake.bash ${CMAKE_BINARY_DIR}/gammcor
)

# cmake --build build --target build-verify
# Build and verify on the current node
add_custom_target(build-verify
    COMMAND cmake --build ${CMAKE_BINARY_DIR} -j 8
    COMMAND ../gammcor_verify_cmake.bash ${CMAKE_BINARY_DIR}/gammcor
)

# cmake --build build --target clean-build
# Clean and build on the current node
add_custom_target(clean-build
    COMMAND cmake --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND cmake --build ${CMAKE_BINARY_DIR} -j 8
)
