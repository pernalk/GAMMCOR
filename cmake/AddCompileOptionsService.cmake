function(add_compile_options)
    set(COMPILE_OPTIONS ${COMPILE_OPTIONS} ${ARGV} PARENT_SCOPE)
    _add_compile_options(${ARGV})
endfunction()