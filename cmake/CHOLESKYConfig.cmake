message(NOTICE "")
message(NOTICE "🗃️  Cholesky info")

include_directories(./gammcor-cholesky/include)

find_library(
    CHOLESKY_LIBRARY 
    NAMES cholesky.a
    PATHS ../gammcor-cholesky/lib/
    NO_DEFAULT_PATH
)

message("CHOLESKY location: .............. ${CHOLESKY_LIBRARY}")
message(NOTICE "")
