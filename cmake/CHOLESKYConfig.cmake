message(NOTICE "")
message(NOTICE "🗃️  Cholesky info")

include_directories(/home/michalhapka/CholeskyOnTheFly/include)

find_library(
    CHOLESKY_LIBRARY 
    NAMES cholesky.a
    PATHS /home/michalhapka/CholeskyOnTheFly/lib/
    NO_DEFAULT_PATH
)

message("CHOLESKY location: .............. ${CHOLESKY_LIBRARY}")
message(NOTICE "")
