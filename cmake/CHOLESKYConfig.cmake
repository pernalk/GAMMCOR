message(NOTICE "")
message(NOTICE "🗃️  Cholesky info")

include_directories(./gammcor-integrals/include)

find_library(
    CHOLESKY_LIBRARY 
    NAMES cholesky.a
    PATHS ./gammcor-integrals/lib/
    NO_DEFAULT_PATH
)

message("CHOLESKY location: .............. ${CHOLESKY_LIBRARY}")
message(NOTICE "")
