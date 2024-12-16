message(NOTICE "")
message(NOTICE "🗃️  TREXIO info")

find_library(
    TREXIO_LIBRARY 
    NAMES trexio
    PATHS /home/hapka/.local/lib
    NO_DEFAULT_PATH
)

message("TREXIO location: .............. ${TREXIO_LIBRARY}")
message(NOTICE "")
