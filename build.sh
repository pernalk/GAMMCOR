#!/bin/bash
# Build GAMMCOR
# Usage: ./build.sh [-c] [-j N]
#   -c   clean before build
#   -j N parallel jobs (default: nproc)

set -e
cd "$(dirname "$0")"
source env.sh

JOBS=$(nproc)
CLEAN=0

while getopts "cj:" opt; do
    case $opt in
        c) CLEAN=1 ;;
        j) JOBS=$OPTARG ;;
    esac
done

if [ $CLEAN -eq 1 ]; then
    echo "Cleaning..."
    make -f Makefile.gcc clean
fi

echo "Building with $JOBS jobs..."
make -f Makefile.gcc -j"$JOBS"
echo "Done. Binary: $(pwd)/gammcor"
